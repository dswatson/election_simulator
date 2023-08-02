# this simulates an election with n voters, k parties, d districts, over t days.
# The central idea is that voters have two-dimensional preferences over social and economic issues. These are shifted by being different by state, and geography (lat-long position). They are also shifted over time, by a random walk in each dimension. Parties also have ideal points in 2-dimensional space, which are fixed. 

# 1. model parameters
# voter n
n <- 10000

# days to model (note: t = 1 is the day of the election)
t <- 100

# parties
k <- 5

# districts
d <- 9

# volatility (change over time in each dimension, a non-interpretable parameter)
volatility <- 5

# district sameness (how similar are voters within a district - higher = greater)
district_sameness <- 0.2

# geography sameness (how predictive is lat / long of votes - higher = greater)
geo_sameness <- 0.2

# 2. simulate parties
parties <- data.frame(econ_pref = rnorm(k),
                      social_pref = rnorm(k))

# 3. simulate voters
voters <- data.frame(econ_pref = rnorm(n),
                     social_pref = rnorm(n),
                     lat = rnorm(n),
                     lng = rnorm(n),
                     district = rep(NA, n))

# 4. generate districts
#  a. generate centroids of districs
districts <- data.frame(lat = rnorm(d),
                        lng = rnorm(d))

#  b. assign voters to districts based on euclidian distance to district centroids. could obviously have districts/voters be at fixed lat / lng positions too.
for(i in 1:n){
  voters$district[i] <- which(sqrt((districts$lat - voters$lat[i])^2 + (districts$lng + voters$lng[i])^2) == min(sqrt((districts$lat - voters$lat[i])^2 + (districts$lng + voters$lng[i])^2)))
}

# 5. model voter preferences and vote choice
# voters social and economic preference is related to their geographical location, and district
voters$econ_pref_full <- voters$econ_pref + geo_sameness*rnorm(1)*voters$lat + geo_sameness*rnorm(1)*voters$lng + district_sameness*rnorm(1)*voters$district

voters$social_pref_full <- voters$social_pref + geo_sameness*rnorm(1)*voters$lat + geo_sameness*rnorm(1)*voters$lng + district_sameness*rnorm(1)*voters$district

# they vote for the party closest to their preference, given by euclidian distance to two-dimensional ideal point
for(i in 1:n){
voters$vote[i] <- which(sqrt((parties$econ_pref - voters$econ_pref_full[i])^2 + (parties$econ_pref - voters$econ_pref_full[i])^2) == min(sqrt((parties$econ_pref - voters$econ_pref_full[i])^2 + (parties$econ_pref - voters$econ_pref_full[i])^2)))
}

# 6. model change over time (working backward)
# next, we simulate the preferences and vote choice over time. we assume a random walk in both economics (economic conditions) and social (issue salience), backwards from election day.
econ_development <- c(0, rep(NA, t-1))
social_development <- c(0, rep(NA, t-1))
for(i in 2:t){
  econ_development[i] <- econ_development[i-1] + volatility*rnorm(1)/t
  social_development[i] <- social_development[i-1] + volatility*rnorm(1)/t
}

voters$t <- 1
voters_day_by_day <- list()
for(j in 1:t){
  temp <- voters
  temp$t <- j
  temp$econ_pref_full <- temp$econ_pref_full + econ_development[j]
  temp$social_pref_full <- temp$social_pref_full + social_development[j]

  for(i in 1:n){
    temp$vote[i] <- which(sqrt((parties$econ_pref - temp$econ_pref_full[i])^2 + (parties$social_pref - temp$social_pref_full[i])^2) == min(sqrt((parties$social_pref - temp$social_pref_full[i])^2 + (parties$econ_pref - temp$econ_pref_full[i])^2)))
  }
  voters_day_by_day <- c(voters_day_by_day, list(temp))
}

# We now have an array of voters for every day of the election.
# To get 'optimal' poll results:
voter_preference_by_day <- data.frame()
for(i in 1:t){
  voter_preference_by_day <- rbind(voter_preference_by_day, table(voters_day_by_day[[i]]$vote))
}
colnames(voter_preference_by_day) <- paste0('party_', 1:ncol(voter_preference_by_day))
head(voter_preference_by_day)

# Plot these for a few parties:
library(ggplot2)
p <- ggplot(voter_preference_by_day, aes(x=1:t))+
  geom_line(aes(y=party_1, col='party_1'))+
  geom_line(aes(y=party_2, col='party_2'))+
  geom_line(aes(y=party_3, col='party_3'))+
  geom_line(aes(y=party_4, col='party_4'))+
  geom_line(aes(y=party_5, col='party_5'))+ylab('')
p

# To get a poll result, simply sample from the full list of voters (and their districts) on any given day. If you would like, you could add polling house effects.

# To get voter preferences over time:
voter_preference_by_day <- data.frame(voters_day_by_day[[1]])
for(i in 1:t){
  voter_preference_by_day[, paste0('day_', i)] <- voters_day_by_day[[t]]$vote
}

# Get ideal points by day:
ideal_points_by_day <- data.frame(voters_day_by_day[[1]])
for(i in 1:t){
  ideal_points_by_day[, paste0('econ_pref_full_', i)] <- voters_day_by_day[[t]]$econ_pref_full
  ideal_points_by_day[, paste0('social_pref_full_', i)] <- voters_day_by_day[[t]]$social_pref_full
}

ggplot(ideal_points_by_day, aes(x=econ_pref_full_1, y=social_pref_full_1, col=as.factor(vote)))+geom_point(alpha = 0.5)+facet_wrap(.~district)+ggtitle('9 states in one country, col=party vote')


