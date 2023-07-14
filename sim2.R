# Load libraries, register cores
library(data.table)
library(Rfast)
library(LDATS)
library(ggplot2)

# Set seed
set.seed(123)

# Hyperparameters
n_sim <- 1e4L
n_parties <- 4L
n_states <- 10L
national_pop <- 1e6L
party_names <- c('juggernaut', 'big_oppo', 'small_oppo', 'pipsqueak')
super_mu <- c(0.3825, 0.31, 0.2575, 0.05)   
names(super_mu) <- party_names

### STATE LEVEL ###

# Average vote shares for states A, B, and C
mu_a <- c(0.25, 0.25, 0.225, 0.275) 
mu_b <- c(0.375, 0.275, 0.275, .075)
mu_c <- c(0.5, 0.125, 0.3, .075)
# Populations per state (this should change for each, no?)
state_pop <- national_pop * 0.175/3
# Region-wide covariance matrix
sigma_r <- matrix(c(0.1,    -0.005,       0,   0,
                 -0.005,       0.1,  -0.025,   0, 
                      0,    -0.025,     0.1,   0,
                      0,         0,       0,   0.15), nrow = 4, byrow = TRUE)
# State numbers
state_fn <- function(state_name, mu_state, sigma_state, pop_state, n_sim) {
  mu <- log(mu_state * pop_state)
  y <- as.data.table(exp(rmvnorm(n_sim, mu = mu, sigma = sigma_state)))
  colnames(y) <- party_names
  y[, state := state_name]
  return(y)
}
x_a <- state_fn('A', mu_a, sigma_r, state_pop, n_sim)
x_b <- state_fn('B', mu_a, sigma_r, state_pop, n_sim)
x_c <- state_fn('C', mu_a, sigma_r, state_pop, n_sim)
x <- rbind(x_a, x_b, x_c)

# Regional numbers
region_fn <- function(state_numbers, sigma_region) {
  mu_region <- rep(0, n_parties)
  vote_shares <- state_numbers / rowSums(state_numbers)
  logit_vote_shares <- qlogis(vote_shares)
  noise <- rmvnorm(n_sim, mu = mu_region, sigma = sigma_region)
  noised_logit_shares <- logit_vote_shares + noise 
  noised_shares <- plogis(noised_logit_shares)
  out <- noised_shares * rowSums(state_numbers)
  return(out)
}
# Compare vote shares for small oppo pre- and post-
x <- as.matrix(x_a[, state := NULL])
x2 <- region_fn(x, sigma_r)
df <- data.table('share' = c(a[, 3], b[, 3]), 
                 'status' = rep(c('before', 'after'), each = 1e4))
ggplot(df, aes(x = share, fill = status)) + 
  geom_density(alpha = 0.25)


# OPTION 1: do this for each state/region until we've got national numbers
# OPTION 2: start with the national and work top-down


### NATIONAL LEVEL ###
mu_n <- log(super_mu * national_pop)
sigma_n <- matrix(c(0.01,   -0.01,       0,        0,
                   -0.01,    0.06,   -0.05,    -0.02, 
                       0,   -0.05,    0.08,        0,
                       0,   -0.02,       0,     0.15), nrow = 4, byrow = TRUE)
# y is a matrix of log vote totals for each party
y <- rmvnorm(n_sim, mu = mu_n, sigma = sigma_n) 
colnames(y) <- party_names
# x is a data.table of vote shares for each party
x <- as.data.table(softmax(y))
colnames(x) <- party_names

# Battery of checks
cor(x)
x[, idx := .I]
x2 <- melt(x, measure.vars = seq_len(n_parties), variable.name = 'party',
           value.name = 'share', variable.factor = FALSE)
x2[, winner := party[which.max(share)], by = idx]
table(x2$winner)
x2[, mean(share), by = party]
x2[, sd(share), by = party]
x[, idx := NULL]



#sd_logit_space <- qlogis(super_mu + x[, sd(share), by = party]$V1 / 2) - qlogis(super_mu)
#var_logit_space <- sd_logit_space^2
#s[2, 3] <- s[3, 2] <- sd_logit_space[2] * sd_logit_space[3] * -0.8















