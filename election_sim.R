# Load libraries, register cores
library(data.table)
library(Rfast)
library(LDATS)

# Set seed
set.seed(123)

# Hyperparameters
n_parties <- 4
n_years <- 1e4
n_states <- 10

### Simulate from a multivariate logistic normal distribution ###
# First, define our "true" latent probability vector, with values for 
# juggernaut, big_oppo, small_oppo, and pipsqueak (respectively)
theta <- c(0.45, 0.3, 0.2, 0.05) 

# We set the reference to zero in log space 
mu <- log(theta) - log(theta[1])
mu <- mu[-1]

# Covariance matrix (excluding juggernaut)
Sigma <- matrix(c(0.1, -0.075, 0, 
                 -0.075,  0.1, 0,
                    0,      0, 0.1), nrow = 3, byrow = TRUE)

# Simulate from multivariate normal, add 0-vector
y <- rmvnorm(n_years, mu = mu, sigma = Sigma)
y <- cbind(rep(0, n_years), y)

# Transform
x <- softmax(y)
x <- as.data.table(x)
colnames(x) <- c('juggernaut', 'big_oppo', 'small_oppo', 'pipsqueak')
cor(x)
x <- melt(x, measure.vars = 1:4, variable.name = 'party', value.name = 'share')
x[, mean(share), by = party]
x[, sd(share), by = party]




# OG version
# Long run average: multivariate normal distro with a buncha fiddling
Sigma <- matrix(c(   1, -0.2, -0.2, -0.3,
                  -0.2,    1, -0.7,    0,
                  -0.2, -0.7,    1,    0,
                  -0.3,    0,    0,    1), nrow = 4, byrow = TRUE)
df <- as.data.table(
  rmvnorm(n_years, mu = c(0, 0, 0, 0), sigma = Sigma)
)
colnames(df) <- c('juggernaut', 'big_oppo', 'small_oppo', 'pipsqueak')
df[, juggernaut := plogis((juggernaut * .05) - 0.3)]
df[, big_oppo := plogis((big_oppo * 0.4) - 1)]
df[, small_oppo := plogis((small_oppo * 0.4) - 1.5)]
df[, pipsqueak := plogis((pipsqueak * 0.6) - 3)]
df[, tmp_sum := juggernaut + big_oppo + small_oppo + pipsqueak]
df[, juggernaut := juggernaut / tmp_sum]
df[, big_oppo := big_oppo / tmp_sum]
df[, small_oppo := small_oppo / tmp_sum]
df[, pipsqueak := pipsqueak / tmp_sum]
df[, tmp_sum := NULL]
cor(df)
df <- melt(df, measure.vars = 1:4, variable.name = 'party', value.name = 'share')
df[, mean(share), by = party]
df[, sd(share), by = party]

# Statewide variation


# brms package version 
# Doesn't preserve the means we want tho...?
x <- as.data.table(rlogistic_normal(n_years, mu, Sigma))











