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
theta <- c(0.3825, 0.31, 0.2575, 0.05) 

# We set the reference to zero in log space 
mu <- log(theta) - log(theta[1])
mu <- mu[-1]

# Covariance matrix (excluding juggernaut)
Sigma <- matrix(c(0.1, -0.025, 0, 
                  -0.025,  0.1, 0,
                  0,      0, 0.15), nrow = 3, byrow = TRUE)

# Simulate from multivariate normal, add 0-vector
y <- rmvnorm(n_years, mu = mu, sigma = Sigma)
y <- cbind(rep(0, n_years), y)

# Transform
x <- softmax(y)
x <- as.data.table(x)
colnames(x) <- c('juggernaut', 'big_oppo', 'small_oppo', 'pipsqueak')
cor(x)
x[, winner := names(x)[rowMaxs(as.matrix(.SD))]]
table(x$winner)
x <- melt(x, measure.vars = 1:4, variable.name = 'party', value.name = 'share')
x[, mean(share), by = party]
x[, sd(share), by = party]

