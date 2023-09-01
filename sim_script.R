# Set working directory
setwd('~/Documents/Omega/Economist/election_simulator')

# Load libraries, register cores
library(data.table)
library(Rfast)
library(LDATS)
library(matrixStats)
library(ggplot2)
library(foreach)

# Set seed
set.seed(123)

# Hyperparameters
n_mc <- 1e6L
n_sim <- 1e4L
n_parties <- 4L
n_states <- 10L
national_pop <- 1e6L
party_names <- c('pdal', 'cc', 'dgm', 'ssp')
super_mu <- c(0.3825, 0.31, 0.2575, 0.05)   
names(super_mu) <- party_names

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

# Sigma function
sigma_fn <- function(r) {
  mu_tmp <- unique(dataland_melt[region == r, .(party, party_regional_vote_share)])$party_regional_vote_share
  target_var <- (coef_var * mu_tmp)^2
  logit_var <- sapply(2:length(mu_tmp), function(k) {
    alpha <- ((1 - mu_tmp[k]) / target_var[k] - 1 / mu_tmp[k]) * mu_tmp[k]^2
    beta <- alpha * (1 / mu_tmp[k] - 1)
    p <- suppressWarnings(rbeta(n_mc, alpha, beta))
    q <- qlogis(p)
    return(var(q))
  })
  logit_var[is.na(logit_var)] <- 0
  sigma_out <- diag(logit_var)
  for (i in 2:(n_parties - 1)) {
    for (j in 1:(i - 1)) {
      sigma_out[i, j] <- sigma_out[j, i] <- 
        sigma_out[i, i] * (Sigma[i, j] / Sigma[i, i])
    }
  }
  return(sigma_out)
}

# State numbers
sim_prov <- function(prov, 
                     yr = NULL, 
                     in_party = NULL, 
                     in_party_num_terms = NULL, 
                     b = n_sim) {
  
  # Grab hyperparameters
  tmp <- dataland_melt[province == prov]
  pop_state <- national_pop * tmp[, unique(population)]
  sigma_r <- sigma_fn(tmp[, unique(region)])
  
  # Grab economic data
  in_party_penalty <- ifelse(in_party_num_terms > 2, -2.5, 0)
  total_in_party_adj <- (nz[year == yr]$in_party_econ_impact + in_party_penalty) / 100
  
  # Make adjustment
  out_party_share_old <- tmp[party != in_party, sum(share)]
  tmp[party == in_party, share := share + total_in_party_adj]
  out_party_share_new <- 1 - tmp[party == in_party, share]
  fctr <- out_party_share_new / out_party_share_old
  tmp[party != in_party, share := share * fctr]
  
  # Compute hyperparameters
  mu_state <- tmp[, share]
  zero_flag <- mu_state[4] == 0
  if (isTRUE(zero_flag)) {
    mu_state <- mu_state[1:3]
    sigma_r <- sigma_r[1:2, 1:2]
  }
  mu_r <- log(mu_state) - log(mu_state[1])
  mu_r <- mu_r[-1]
  
  # Simulate, transform, export
  y <- rmvnorm(b, mu_r, sigma_r)
  y <- cbind(rep(0, b), y)
  x <- softmax(y)
  x <- as.data.table(x)
  if (isTRUE(zero_flag)) {
    x[, ssp := 0L]
  }
  setnames(x, party_names)
  x[, province := prov][, year := yr]
  x <- merge(x, unique(dataland_melt[, .(province, region)]), by = 'province')
  return(x)
}

# National or regional noise function
noise_fn <- function(dat, noise_sd, regional = FALSE) {
  # Elasticity
  if (isTRUE(regional)) {
    dat$elasticity <- 1
  } else {
    dat <- merge(dat, unique(dataland_melt[, .(province, elasticity)]), 
                 by = 'province')
  }
  # Hyperparameters
  input_mat <- as.matrix(dat[, .(pdal, cc, dgm, ssp)])
  n <- nrow(input_mat)
  # Noise out
  out <- matrix(rnorm(n * 4, sd = noise_sd * dat$elasticity), ncol = 4)
  return(out)
}

# Ethnic noise is harder
grd <- setDT(expand.grid(province = dataland_melt[, unique(province)],
                         party = party_names,
                         ethnicity = c('python', 'cobolite', 'javarian')))
ethnic_noise_fn <- function(dat, noise_sd = 0.125) {
  n <- nrow(dat) / 5
  ethnic_noise <- matrix(rnorm(n * 12, sd = noise_sd), ncol = 12)
  inner_loop <- function(prov, party) {
    tmp <- grd[province == prov]
    tmp <- merge(tmp, dataland_melt[, .(province, party, python_pop_share,
                                        cobolite_pop_share, javarian_pop_share)],
                 by = c('province', 'party'))
    tmp <- melt(tmp, id.vars = c('province', 'party', 'ethnicity'),
                variable.name = 'weight', value.name = 'pop_share')
    tmp[, weight := gsub('_pop_share', '', weight)]
    tmp <- tmp[ethnicity == weight][, weight := NULL]
    if (party == 'pdal') {
      party_noise <- ethnic_noise[, 1:3]
    } else if (party == 'cc') {
      party_noise <- ethnic_noise[, 4:6]
    } else if (party == 'dgm') {
      party_noise <- ethnic_noise[, 7:9]
    } else if (party == 'ssp') {
      party_noise <- ethnic_noise[, 10:12]
    }
    eps <- as.numeric(party_noise %*% tmp[, unique(pop_share)])
    return(eps)
  }
  eps_eth <- foreach(provs = dataland_melt[, unique(province)], .combine = rbind) %:%
    foreach(parties = party_names, .combine = cbind) %do%
    inner_loop(provs, parties)
  colnames(eps_eth) <- party_names
  return(eps_eth)
}


# Start the Markov chain
in_power <- 'pdal'
yr <- 1984
no_terms <- 0
res <- data.table(
  'year' = numeric(),
  'incumbent' = character(),
  'province' = character(),
  'party' = character(),
  'share' = numeric(),
  'winner_prov' = character(),
  'winner_natl' = character()
)
cols <- colnames(res)
while (yr < 2024) {
  prov_sim1 <- foreach(p = dataland_melt[, unique(province)], .combine = rbind) %do%
    sim_prov(p, yr = yr, in_party = in_power, in_party_num_terms = no_terms, b = 1)
  eps_nat <- noise_fn(prov_sim1, noise_sd = 0.125)
  eps_reg <- foreach(r = prov_sim1[, unique(region)], .combine = rbind) %do%
    noise_fn(prov_sim1[region == r], noise_sd = 1/16, TRUE)
  eps_eth <- ethnic_noise_fn(prov_sim1)
  eps <- eps_nat + eps_reg + eps_eth
  x <- prov_sim1[, .(pdal, cc, dgm, ssp)]
  x <- x * exp(eps)
  prov_res <- as.data.table(x / rowSums(x))
  prov_res[, province := prov_sim1$province][, year := yr]
  prov_res <- melt(prov_res, id.vars = c('province', 'year'), variable.name = 'party', 
            value.name = 'share')
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

