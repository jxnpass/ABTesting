library(tidyverse)
library(bayesAB)
source("theme.R")
rename <- dplyr::rename

# Loading in Data
set.seed(39)
data <- read_csv("payment_fraud.csv") %>% 
  rename(X = numItems, Y = label) %>% 
  select(X, Y)

# Change to binary for groups if needed
data$Y <- ifelse(data$Y == 1, T, F) # is fraud
data$X <- ifelse(data$X > 1, T, F) # number of purchases > 1

data %>% 
  group_by(X, Y) %>% 
  summarise(N = n()) %>% 
  mutate(Proportion = N/sum(N)) 

#############################################################################
#                         BAYESIAN A/B Testing                              # 
#             Beta-binomial conjugate posterior distribution.               #
#############################################################################

### Proportions test on variable "Drained_after_socializing" ###

### Step 1: Prior Distribution

## Historical Prior
# We have a small amount of data prior to collection that suggests
# X = TRUE is more likely...

# Group = FALSE (~5% are y = TRUE)
f_a_hp <- 05
f_b_hp <- 95

# GROUP = TRUE (~10% are y = TRUE)
t_a_hp <- 10
t_b_hp <- 90

curve(dbeta(x, f_a_hp, f_b_hp), from = 0, to = 1, 
      col = "blue", lwd = 2, 
      xlab = "Probability of Fraud",
      ylab = "Density",
      main = "Historical prior distributions of fraud",
      )

curve(dbeta(x, t_a_hp, t_b_hp), from = 0, to = 1, 
      col = "red", lwd = 2, add = TRUE)

# Legend
legend("topright", 
       legend = c("X = 1 Purchase", 
                  "X = 2+ Purchases"),
       col = c("blue", "red"), lwd = 2, bty = "n")


## Non-informative prior
# We believe pre-study that there is no difference among groups. 
# We set up priors like the following (believe proportion is just small)
f_a_np <- t_a_np <- 2
f_b_np <- t_b_np <- 8

curve(dbeta(x, f_a_np, f_b_np), from = 0, to = 1, 
      col = "blue", lwd = 2, 
      xlab = "Probability of Fraud",
      ylab = "Density",
      main = "Non-informative prior distributions of fraud",
)

curve(dbeta(x, t_a_np, t_b_np), from = 0, to = 1, 
      col = "red", lwd = 1, add = TRUE)

# Legend
legend("topright", 
       legend = c("X = 1 Purchase", 
                  "X = 2+ Purchases"),
       col = c("blue", "red"), lwd = 2, bty = "n")

### Step 2: Generate conjugate posterior distribution

### Option 1: Direct conjugate
### directly compute posterior distribution:

f_s <- sum(data$Y[data$X == F])
f_astar_np <- f_a_np + f_s # non-informative
f_astar_hp <- f_a_hp + f_s # historical

f_n <- length(data$Y[data$X == F])
f_bstar_np <- f_b_np + f_n - f_s # non-informative
f_bstar_hp <- f_b_hp + f_n - f_s # historical

t_s <-  sum(data$Y[data$X == T])
t_astar_np <- t_a_np + t_s # non-informative
t_astar_hp <- t_a_hp + t_s # historical

t_n <- length(data$Y[data$X == T])
t_bstar_np <- t_b_np + t_n - t_s # non-informative
t_bstar_hp <- t_b_hp + t_n - t_s #historical

### Historical conjugate posterior curves
curve(dbeta(x, f_astar_hp, f_bstar_hp), from = 0, to = 1, col = "red", lwd = 2,
      xlab = "Probability of Fraud",
      ylab = "Density",
      main = "Historical posterior distributions of Fraud")
curve(dbeta(x, t_astar_hp, t_bstar_hp), from = 0, to = 1, col = "blue", lwd = 2, add = T)
legend("topright", 
       legend = c("X = 1 Purchase", 
                  "X = 2+ Purchases"),
       col = c("blue", "red"), lwd = 2, bty = "n")

# Distribution of differences
gA <- rbeta(10000, t_astar_hp, t_bstar_hp)
gB <- rbeta(10000, f_astar_hp, f_bstar_hp)
hist(gA - gB, main = "Group A - Group B Proportions",
     xlab = "Proportion Differences")


### Non-informative conjugate posterior curves
# note, depending on historical strength and dataset, this may look very similar
# to the previous graphic

curve(dbeta(x, f_astar_np, f_bstar_np), from = 0, to = 1, col = "red", lwd = 2,      
      xlab = "Probability of Fraud",
      ylab = "Density",
      main = "Non-Informative posterior distributions of Fraud")
curve(dbeta(x, t_astar_np, t_bstar_np), from = 0, to = 1, col = "blue", lwd = 2, add = T)
legend("topright", 
       legend = c("X = 1 Purchase", 
                  "X = 2+ Purchases"),
       col = c("blue", "red"), lwd = 2, bty = "n")

### Option 2: Metropolis-Hastings or Gibbs Sampler, using non-informative priors
# because there is no closed-form posterior distribution for two groups
# we use sampling to determine what the difference might actually be...

n_iters <- 1000
x <- c(t_s, f_s)
a_prior <- c(t_a_hp, f_a_hp)
b_prior <- c(t_b_hp, f_b_hp)
n <- c(t_n, f_n)
J <- length(x)
theta_samples <- matrix(NA, nrow = n_iters, ncol = J)

# Gibbs sampling
for (i in 1:n_iters) {
  for (j in 1:J) {
    a_post <- a_prior[j] + x[j]
    b_post <- b_prior[j] + n[j] - x[j]
    theta_samples[i, j] <- rbeta(1, a_post, b_post)
  }
}

# Create histograms of draws
h1 <- hist(theta_samples[,1], breaks = 15, plot = FALSE)
h2 <- hist(theta_samples[,2], breaks = 15, plot = FALSE)
plot(h1, col = rgb(1, 0, 0, 0.5), freq = FALSE, 
     xlim = range(h1$breaks, h2$breaks), 
     ylim = c(0, max(h1$density, h2$density)))
plot(h2, col = rgb(0, 0, 1, 0.5), freq = FALSE, add = TRUE)
curve(dbeta(x, f_astar_np, f_bstar_np), from = 0, to = 1, lwd = 4, col = "blue", add = T)
curve(dbeta(x, t_astar_np, t_bstar_np), from = 0, to = 1, lwd = 4, col = "red", add = T)

hist(theta_samples[,1] - theta_samples[,2], breaks = 100)

### Option 3: Use library bayesAB in R 
library(bayesAB)
t_data <- as.numeric(data$Y[data$X == T])
f_data <- as.numeric(data$Y[data$X == F])

abtest <- bayesTest(t_data, f_data,
          priors = c('alpha' = f_a_np, 'beta' = f_b_np),
          distribution = 'bernoulli')

plot(abtest)
summary(abtest)
# your final message here: probability group = TRUE proportion > group = FALSE proportion
summary(abtest)$probability$Probability


### Step 3: Treatment testing and comparison between A and B 
# what happens here is that we take the random draws from the MCMC simulation, and determine
# what percentage of draws has it such that the TRUE GROUP proportion is greater than
# the FALSE GROUP proportion

# probability that GROUP 1 > GROUP 2
mean(theta_samples[,1] > theta_samples[,2])
my_test <- mean(gA > gB)
quantile(gA - gB, probs = c(.025, .975))
ci <- quantile(gA - gB, probs = c(.025, .975))
paste("[", ci[1], ", ", ci[2], "]", sep = "")




