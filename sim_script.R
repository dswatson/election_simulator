# Set working directory
setwd('~/Documents/Omega/Economist/election_simulator')

# Load libraries, register cores
library(data.table)
library(Rfast)
library(LDATS)
library(matrixStats)
library(ggplot2)

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

# Master data frame (ENSURE SUMMATION TO UNITY)
dataland_melt <- merge(melt(
  fread('dataland.csv')[, -c('population', 'region'), with = F], id.vars = 'province',
  variable.name = 'party', value.name = 'share'
), fread('dataland.csv')[, .(region, population, province)], by = 'province')

dataland_melt[, region_share := weighted.mean(share, population), by = .(region, party)]


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
  mu_tmp <- unique(dataland_melt[region == r, .(party, region_share)])$region_share
  target_var <- (coef_var * mu_tmp)^2
  logit_var <- sapply(2:length(mu_tmp), function(k) {
    alpha <- ((1 - mu_tmp[k]) / target_var[k] - 1 / mu_tmp[k]) * mu_tmp[k]^2
    beta <- alpha * (1 / mu_tmp[k] - 1)
    p <- rbeta(n_mc, alpha, beta)
    q <- qlogis(p)
    return(var(q))
  })
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
state_fn <- function(prov, 
                     yr = NULL, 
                     in_party = NULL, 
                     in_party_num_terms = NULL, 
                     b = n_sim) {
  
  # Grab hyperparameters
  tmp <- dataland_melt[province == prov]
  pop_state <- national_pop * tmp[, unique(population)]
  sigma_r <- sigma_fn(tmp[, unique(region)])
  
  # Grab economic data
  in_party_penalty <- ifelse(in_party_num_terms > 1, -2.5, 0)
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
  return(x)
}

# Example 
a <- state_fn('Cerebrica', 
              yr = 1984, 
              in_party = 'pdal', 
              in_party_num_terms = 1)
x <- a[, 1:4]
colMeans(x)
cov(x)








