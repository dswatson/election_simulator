# Set working directory
setwd('~/Documents/Omega/Economist/election_simulator')

# Load libraries, register cores
library(data.table)
library(Rfast)
library(LDATS)
library(matrixStats)
library(ggplot2)
library(ggsci)
library(foreach)

# Set seed
set.seed(123)

# Hyperparameters
n_mc <- 1e4L
n_sim <- 1e4L
n_parties <- 4L
n_states <- 10L
national_pop <- 1e6L
party_names <- c('pdal', 'cc', 'dgm', 'ssp')
super_mu <- c(0.3825, 0.31, 0.2575, 0.05)   
names(super_mu) <- party_names
mu0 <- rep(0, 3)

# Import fundamentals
nz <- fread('nz_econ.csv')
nz <- nz[year <= 2024]
nz[, inflation_delta := inflation - shift(inflation, 1)]
nz[, unemployment_delta := unemployment - shift(unemployment, 1)]
nz[, in_party_econ_impact := -(2 * unemployment_delta) - (0.5 * inflation_delta)]

# Master data frame 
dataland_melt <- merge(melt(
  fread('dataland.csv')[, -c('population', 'region', 'python_pop_share', 'cobolite_pop_share', 
                             'javarian_pop_share', 'elasticity'), with = F], id.vars = 'province',
  variable.name = 'party', value.name = 'share'
), fread('dataland.csv')[, .(region, population, province, python_pop_share, cobolite_pop_share, 
                             javarian_pop_share, elasticity)], by = 'province')
dataland_melt[, party_regional_vote_share := weighted.mean(share, population), by = .(region, party)]

# Calculate electoral votes
ec_df <- data.table(province = dataland_melt[, unique(province)],
                    ec_votes = 1L, population = dataland_melt[, unique(population)])
vote_being_allocated <- nrow(ec_df) + 1
while (vote_being_allocated <= 100) {
  ec_df[, priority_number := population / sqrt(ec_votes * (ec_votes + 1))]
  ec_df[which.max(priority_number), ec_votes := ec_votes + 1]
  vote_being_allocated <- vote_being_allocated + 1
}
ec_df[, ec_votes := ec_votes + 2]

### REGION LEVEL ###

# Compute coefficient of variation (requires Monte Carlo)
Sigma <- matrix(c(0.1, -0.025, 0, 
               -0.025,    0.1, 0,
                    0,      0, 0.15), nrow = 3, byrow = TRUE)
mu <- log(super_mu) - log(super_mu[1])
mu <- mu[-1]
y <- rmvnorm(n_mc, mu = mu, sigma = Sigma)
y <- cbind(rep(0, n_mc), y)
x <- softmax(y)
coef_var <- colSds(x) / super_mu
rm(x, y)

# Sigma function
sigma_fn <- function(sigma_in, mu_in, days_out, n_mc) {
  if (days_out == 0) {
    target_var <- (coef_var * mu_in)^2
  } else {
    target_var <- (coef_var / 20 * mu_in)^2 # THIS IS A KNOB TO TUNE
  }
  logit_var <- sapply(2:length(mu_in), function(k) {
    alpha <- ((1 - mu_in[k]) / target_var[k] - 1 / mu_in[k]) * mu_in[k]^2
    beta <- alpha * (1 / mu_in[k] - 1)
    p <- suppressWarnings(rbeta(n_mc, alpha, beta))
    q <- qlogis(p)
    return(var(q))
  })
  logit_var[is.na(logit_var)] <- 0
  sigma_out <- diag(logit_var)
  for (i in 2:(n_parties - 1)) {
    for (j in 1:(i - 1)) {
      sigma_out[i, j] <- sigma_out[j, i] <- 
        sigma_out[i, i] * (sigma_in[i, j] / sigma_in[i, i])
    }
  }
  return(sigma_out)
}

# Simulation function
sim_fn <- function(yr, 
                   in_party = NULL, 
                   in_party_num_terms = NULL, 
                   days_out = 0,
                   data_in = NULL,
                   b = n_sim) {
  
  # Scale variance for random walks
  if (days_out > 0) {
    this_sigma <- Sigma / 400 # THIS IS A KNOB TO TUNE
    if (days_out > 1) {
      b <- data_in[, .N / 20]
    }
  } else {
    this_sigma <- Sigma
  }
  
  # Grab economic data
  if (days_out == 0) {
    in_party_penalty <- ifelse(in_party_num_terms > 2, -2.5, 0)
    total_in_party_adj <- (nz[year == yr]$in_party_econ_impact + in_party_penalty) / 100
  }
  
  # Define national noise first
  ntl_noise <- rmvnorm(b, mu = mu0, sigma = this_sigma * 0.45)
  
  # Define ethnic noise 
  eth_noise1 <- rmvnorm(b, mu = mu0, sigma = this_sigma * 0.25)
  eth_noise2 <- rmvnorm(b, mu = mu0, sigma = this_sigma * 0.25)
  eth_noise3 <- rmvnorm(b, mu = mu0, sigma = this_sigma * 0.25) 
  
  # Region function
  region_fn <- function(r) {
    
    # Regional noise
    reg_share <- dataland_melt[region == r, unique(party_regional_vote_share)]
    sigma_r <- sigma_fn(this_sigma, reg_share, days_out, n_mc)
    if (r != 'Synapse Territories') {
      mu0 <- c(0, 0)
      sigma_r <- sigma_r[1:2, 1:2]
      ntl_noise <- ntl_noise[, 1:2]
    }
    reg_noise <- rmvnorm(b, mu = mu0, sigma = sigma_r * 0.2)
    
    prov_fn <- function(p) {
      
      # Provincial noise
      tmp <- dataland_melt[province == p]
      prov_share <- tmp[, share]
      sigma_p <- sigma_fn(this_sigma, prov_share, days_out, n_mc)
      if (r != 'Synapse Territories') {
        sigma_p <- sigma_p[1:2, 1:2]
      }
      prov_noise <- rmvnorm(b, mu = mu0, sigma = sigma_p * 0.1)
      
      # Elasticity factor
      ntl_noise <- ntl_noise * tmp[, unique(elasticity)]
      
      # Ethnic noise
      eth_noise1 <- eth_noise1 * tmp$python_pop_share[1]
      eth_noise2 <- eth_noise2 * tmp$python_pop_share[2]
      eth_noise3 <- eth_noise3 * tmp$python_pop_share[3]
      eth_noise <- eth_noise1 + eth_noise2 + eth_noise3
      if (r != 'Synapse Territories') {
        eth_noise <- eth_noise[, 1:2]
      }
      
      # Total noise
      eps <- ntl_noise + reg_noise + prov_noise + eth_noise
      
      if (days_out == 0) {
        # Add economic impact
        out_party_share_old <- tmp[party != in_party, sum(share)]
        tmp[party == in_party, share := share + total_in_party_adj]
        out_party_share_new <- 1 - tmp[party == in_party, share]
        fctr <- out_party_share_new / out_party_share_old
        tmp[party != in_party, share := share * fctr]
        prov_share <- tmp[, share]
      } else if (days_out == 1) {
        # Just a single vector of 4 vote shares
        prov_share <- data_in[year == yr & province == p, share]
      } else {
        # A potentially large matrix of n_sim * 4 vote shares
        prov_share <- dcast(data_in[year == yr & province == p], 
                            sim_idx ~ party, value.var = 'share')[, sim_idx := NULL]
      }
      
      # Transform, export
      if (days_out < 2) {
        mu_p <- log(prov_share) - log(prov_share[1])
        mu_p <- mu_p[-1]
        if (r != 'Synapse Territories') {
          mu_p <- mu_p[-3]
        }
        y <- cbind(rep(0L, b), t(mu_p + t(eps)))
      } else {
        mu_p <- log(prov_share) - prov_share[, log(pdal)]
        mu_p[, pdal := NULL]
        if (r != 'Synapse Territories') {
          mu_p[, ssp := NULL]
        }
        y <- cbind(rep(0L, b), as.matrix(mu_p) + eps)
      }
      x <- softmax(y)
      if (r != 'Synapse Territories') {
        x <- cbind(x, rep(0L, b))
      }
      x <- as.data.table(x)
      setnames(x, party_names)
      x[, province := p][, year := yr]
      out_p <- merge(x, unique(dataland_melt[, .(province, region)]), by = 'province')
      return(out_p)
    }
    out_r <- foreach(prov = dataland_melt[region == r, unique(province)], .combine = rbind) %do% 
      prov_fn(prov) 
    return(out_r)
  }
  out <- foreach(reg = dataland_melt[, unique(region)], .combine = rbind) %do% 
    region_fn(reg)
  out[, days_out := days_out]
  return(out)
}

# Start the Markov chain -- final results
in_power <- 'pdal'
yr <- 1984
no_terms <- 0
res <- data.table(
  'year' = integer(),
  'incumbent' = character(),
  'province' = character(),
  'region' = character(),
  'party' = character(),
  'share' = numeric(),
  'winner_prov' = character(),
  'winner_natl' = character(),
  'days_out' = integer()
)
cols <- colnames(res)
while (yr < 2024) {
  x <- sim_fn(yr, in_power, no_terms, days_out = 0, data_in = NULL, b = 1)
  x[, year := yr]
  prov_res <- melt(x, id.vars = c('province', 'year', 'region', 'days_out'), 
                   variable.name = 'party', value.name = 'share')
  prov_res[, winner_prov := party[which.max(share)], by = province]
  natl_res <- merge(unique(prov_res[, province, winner_prov]), ec_df[, .(province, ec_votes)],
                    by = 'province')
  total <- natl_res[, sum(ec_votes), by = winner_prov]
  winner <- total[which.max(V1), as.character(winner_prov)]
  # Update, loop
  out <- prov_res
  out[, incumbent := in_power][, winner_natl := winner]
  setcolorder(out, cols)
  res <- rbind(res, out)
  no_terms <- ifelse(winner == in_power, no_terms + 1, 0)
  in_power <- winner
  yr <- yr + 2
}

# New Markov chain -- random walks
d <- 1
random_walk <- function(yr, n_sim) {
  res_rw <- data.table(
    'year' = integer(),
    'province' = character(),
    'region' = character(),
    'party' = character(),
    'share' = numeric(),
    'days_out' = integer(),
    'sim_idx' = integer()
  )
  cols <- colnames(res_rw)
  df_init <- res[year == yr]
  while (d <= 200) {
    x <- sim_fn(yr, days_out = d, data_in = df_init, b = n_sim)
    x[, year := yr][, sim_idx := rep(1:n_sim, times = 5)]
    x <- melt(x, id.vars = c('province', 'year', 'region', 'days_out', 'sim_idx'), 
              variable.name = 'party', value.name = 'share')
    setorder(x, 'year', 'province', 'sim_idx')
    # Update, loop
    df_init <- x
    res_rw <- rbind(res_rw, x)
    setcolorder(res_rw, cols)
    d <- d + 1
  }
  return(res_rw)
}
df_1984 <- random_walk(yr = 1984, n_sim = 2000)


# Spot check a single trajectory for a given province
tmp <- df_1984[sim_idx == 2000 & province == 'Cerebrica']
ggplot(tmp, aes(-days_out, share, color = party)) + 
  geom_line(linewidth = 0.75) + 
  scale_color_d3() +
  theme_bw()

# Compare beginning and end of cycle vote shares
tmp <- df_1984[province == 'Cerebrica' & party == 'pdal' & days_out %in% c(1, 200)]
setorder(tmp, sim_idx)
tmp[, delta := ifelse(days_out == 200, share - shift(share, 1), NA_real_)]
hist(tmp[!is.na(delta), delta], breaks = 40)
tmp[!is.na(delta), mean(delta)]
tmp[!is.na(delta), sd(delta)]





