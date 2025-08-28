library(tidyverse)
library(invgamma)
source("theme.R")

# Loading in Data
# For this script, preset X as your groups to compare and Y as your target
data <- read_csv("medical_costs.csv") %>% 
  filter(Smoker == "no") %>% 
  rename(Y = `Medical Cost`, X = BMI) %>% 
  select(X, Y)

# Change to binary for groups if needed
data$X <- ifelse(data$X > 25, T, F) # unhealthy according to national BMI index

#############################################################################
#                            BAYESIAN A/B Testing                           # 
#                  A two-sample normal-normal distribution for means        #
#############################################################################

### Step 1: Prior Distribution

# We claim historicity for means, but no prior knowledge on variance 
# of each group

# GROUP MEANS (Historical)
# Group = FALSE (Mean is ~$5k)
f_m_hp <- 5000
f_k_hp <- .0025

# Group = TRUE (Mean is ~$6k)
t_m_hp <- 6000
t_k_hp <- .0025

curve(dnorm(x, f_m_hp, 1/f_k_hp), from = 3000, to = 8000, 
      col = "blue", lwd = 2, 
      xlab = "Medical Cost Means",
      ylab = "Density",
      main = "Historical prior distributions of Medical Cost Means",
)

curve(dnorm(x, t_m_hp, 1/t_k_hp), from = 3000, to = 8000, 
      col = "red", lwd = 2, add = TRUE)

legend("topright", 
       legend = c("X = Low-Normal BMI", 
                  "X = High BMI"),
       col = c("blue", "red"), lwd = 2, bty = "n")

# GROUP VARIANCES (Non-informative)
# can do this based on a inverse-gaussian, inverse-chi^2, or another normal dist.
# I prefer the first one, which typically models better than a normal dist. and
# a better conjugate result than an inv-chi^2

# GROUP = FALSE
f_a_np <- 0.005 
f_b_np <- 0.005 

# GROUP = TRUE
t_a_np <- 0.01 
t_b_np <- 0.01 

# Inv-Gamma distribution of non-informative prior
# curve(dinvgamma(x, f_a_np, f_b_np), from = 0, to = 5000, 
#       col = "blue", lwd = 2, 
#       xlab = "Medical Cost Variance",
#       ylab = "Density",
#       main = "Non-informative prior distributions of Medical Cost Variance",
# )
# 
# curve(dinvgamma(x, t_a_np, t_b_np), from = 0, to = 5000, 
#       col = "red", lwd = 2, add = TRUE)
# 
# # Legend
# legend("topright", 
#        legend = c("X = Low-Normal BMI", 
#                   "X = High BMI"),
#        col = c("blue", "red"), lwd = 2, bty = "n")


### Step 2: Generate conjugate posterior distribution -----------
### Option 1: Direct conjugate
# Using priors:
#   mu ~ Normal(m, v)  
#   sigma2 ~ InvGamma(a0, b0) 
# Posterior derivation 
#   A = n + 1/v0
#   m* = (n * ybar + mu0 / v0) / A
#   v* = 1 / A
#   a* = a0 + n/2
#   b* = b0 + 0.5 * sum((y - ybar)^2) + 0.5 * n * (ybar - mu0)^2 / (n * v0 + 1)
# marginal posterior variance of mu (if you want a normal approx) =
#    Var(mu) = b* / ((a* - 1) * A)   (only valid if a* > 1)

# GROUP = FALSE
f_n <- length(data$Y[data$X == F])
f_mu <- mean(data$Y[data$X == F])
f_sq_err <- sum((data$Y[data$X == F] - f_mu)^2)

# compute historical posterior ~N(m*, v*)
f_kstar_hp <- f_k_hp + f_n
f_mstar_hp <- (f_k_hp * f_m_hp + f_n * f_mu) / f_kstar_hp

# compute non-informative posterior ~INV-G(a*, b*)
f_astar_hp <- f_a_hp + f_n / 2
f_bstar_hp <- f_b_hp + 0.5 * f_sq_err + (f_k_hp * f_n) / (2 * f_kstar_hp) * (f_mu - f_m_hp)^2
f_sigma2 <- f_bstar_hp / (f_astar_hp - 1)


# GROUP = TRUE
t_n <- length(data$Y[data$X == T])
t_mu <- mean(data$Y[data$X == T])
t_sq_err <- sum((data$Y[data$X == T] - t_mu)^2)

# compute historical posterior ~N(m*, v*)
t_kstar_hp <- t_k_hp + t_n
t_mstar_hp <- (t_k_hp * t_m_hp + t_n * t_mu) / t_kstar_hp

# compute non-informative posterior ~INV-G(a*, b*)
t_astar_hp <- t_a_hp + t_n / 2
t_bstar_hp <- t_b_hp + 0.5 * t_sq_err + (t_k_hp * t_n) / (2 * t_kstar_hp) * (t_mu - t_m_hp)^2
t_sigma2 <- t_bstar_hp / (t_astar_hp - 1)

### Conjugate posterior curves for normal (mean)
curve(dnorm(x, t_mstar_hp, sqrt(t_sigma2/t_kstar_hp)), from = 5500, to = 6250, col = "red", lwd = 2,
      xlab = "Medical Cost Means",
      ylab = "Density",
      main = "Historical posterior distributions of Medical Cost Means",
)

curve(dnorm(x, f_mstar_hp, sqrt(f_sigma2/f_kstar_hp)), from = 5500, to = 6250, 
      col = "blue", lwd = 2, add = TRUE)

legend("topright", 
       legend = c("X = Low-Normal BMI", 
                  "X = High BMI"),
       col = c("blue", "red"), lwd = 2, bty = "n")

### Conjugate posterior curves for inv-gamma (variance)
curve(dinvgamma(x, t_astar_np, t_bstar_np), from = 600000, to = 800000, col = "red", lwd = 2,
      xlab = "Medical Cost Variance",
      ylab = "Density",
      main = "Historical posterior distributions of Medical Cost Variance",
)

curve(dinvgamma(x, f_astar_np, f_bstar_np), from = 600000, to = 800000, 
      col = "blue", lwd = 2, add = TRUE)

legend("topright", 
       legend = c("X = Low-Normal BMI", 
                  "X = High BMI"),
       col = c("blue", "red"), lwd = 2, bty = "n")

### Option 2: Metropolis-Hastings or Gibbs Sampler, using non-informative priors
# because there is no closed-form posterior distribution for two groups
# we use sampling to determine what the difference might actually be...

n_iters <- 10000
burn_in <- 100

y_t <- data$Y[data$X == T]
y_f <- data$Y[data$X == F]

t_n <- length(y_t)
f_n <- length(y_f)

m_prior <- c(t_m_hp, f_m_hp)
v_prior <- c(t_v_hp, f_v_hp)
k_prior <- 1 / v_prior

a_prior <- c(t_a_np, f_a_np)
b_prior <- c(t_b_np, f_b_np)

n <- c(t_n, f_n)

J <- 2

mu_samples <- matrix(NA, nrow = n_iters, ncol = J)
sigma_samples <- matrix(NA, nrow = n_iters, ncol = J)

# initialize values for sigma and mu for each group
init_sigma <- c(var(y_t), var(y_f))
init_mu <- c(mean(y_t), mean(y_f))

for (i in 1:n_iters) {
  for (j in 1:J) {
    y_j <- if(j == 1) y_t else y_f
    n_j <- n[j]
    ybar_j <- mean(y_j)
    
    k0 <- k_prior[j]
    k_n <- k0 + n_j
    
    # get current sigma
    sigma_cur <- ifelse(i == 1, init_sigma[j], sigma_samples[i-1, j])
    
    # sample mu | sigma2, data
    m_post <- (k0 * m_prior[j] + n_j * ybar_j) / k_n
    sd_post <- sqrt(sigma_cur / k_n)
    mu_new <- rnorm(1, mean = m_post, sd = sd_post)
    
    # sample sigma2 | mu, data
    a_post <- a_prior[j] + n_j / 2
    ss_mu <- sum((y_j - mu_new)^2)
    b_post <- b_prior[j] + 0.5 * ss_mu + 0.5 * k0 * (mu_new - m_prior[j])^2
    
    mu_samples[i, j] <- mu_new
    sigma_samples[i, j] <- rinvgamma(1, shape = a_post, rate = b_post)
  }
}

# Remove burn-in
mu_samples <- mu_samples[(burn_in + 1):n_iters,]
sigma_samples <- sigma_samples[(burn_in + 1):n_iters,]

# Create histograms of draws (treatment = col 1 red, control = col 2 blue)
h1 <- hist(mu_samples[,1], plot = FALSE)
h2 <- hist(mu_samples[,2], plot = FALSE)
plot(h1, col = rgb(1, 0, 0, 0.5), freq = FALSE, 
     xlim = range(h1$breaks, h2$breaks), 
     ylim = c(0, max(h1$density, h2$density)))
plot(h2, col = rgb(0, 0, 1, 0.5), freq = FALSE, add = TRUE)
curve(dnorm(x, t_mstar_hp, sqrt(t_sigma2/t_kstar_hp)), col = "red", lwd = 2, add = T)
curve(dnorm(x, f_mstar_hp, sqrt(f_sigma2/f_kstar_hp)), col = "blue", lwd = 2, add = TRUE)

# Create histograms of draws (treatment = col 1 red, control = col 2 blue)
h1 <- hist(sigma_samples[,1], plot = FALSE)
h2 <- hist(sigma_samples[,2], plot = FALSE)
plot(h1, col = rgb(1, 0, 0, 0.5), freq = FALSE, 
     xlim = range(h1$breaks, h2$breaks), 
     ylim = c(0, max(h1$density, h2$density)))
plot(h2, col = rgb(0, 0, 1, 0.5), freq = FALSE, add = TRUE)
curve(dinvgamma(x, t_astar_np, t_bstar_np), col = "red", lwd = 2, add = T)
curve(dinvgamma(x, f_astar_np, f_bstar_np), col = "blue", lwd = 2, add = TRUE)
legend("topright", 
       legend = c("X = Low-Normal BMI", 
                  "X = High BMI"),
       col = c("blue", "red"), lwd = 2, bty = "n")

### Option 3: Use library bayesAB in R 
library(bayesAB)
t_data <- as.numeric(data$Y[data$X == T])
f_data <- as.numeric(data$Y[data$X == F])

abtest <- bayesTest(t_data, f_data,
                    priors = c('mu' = f_m_hp, 'lambda' = f_k_hp,
                               'alpha' = f_a_np, 'beta' = f_b_np),
                    distribution = 'normal')

plot(abtest)
summary(abtest)
# your final message here: probability group = TRUE mean > group = FALSE mean
summary(abtest)$probability$Mu
# your final message here: probability group = TRUE var > group = FALSE var
summary(abtest)$probability$Sig_Sq

### Step 3: Treatment testing and comparison between A and B 
# what happens here is that we take the random draws from the MCMC simulation, and determine
# what percentage of draws has it such that the TRUE GROUP means is greater than
# the FALSE GROUP means

# probability that GROUP 1 > GROUP 2
mean(mu_samples[,1] > mu_samples[,2])
mean(sigma_samples[,1] > sigma_samples[,2])



