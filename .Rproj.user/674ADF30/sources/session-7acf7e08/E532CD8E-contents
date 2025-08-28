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

# Group = FALSE (~5% are y = TRUE)
f_a_hp <- 05
f_b_hp <- 95

# GROUP = TRUE (~10% are y = TRUE)
t_a_hp <- 10
t_b_hp <- 90

# Data
f_s <- sum(data$Y[data$X == F])
f_n <- length(data$Y[data$X == F])
t_s <-  sum(data$Y[data$X == T])
t_n <- length(data$Y[data$X == T])


n_iters <- 1000
x <- c(t_s, f_s)             # successes
a_prior <- c(t_a_hp, f_a_hp) # Beta prior alpha
b_prior <- c(t_b_hp, f_b_hp) # Beta prior beta
n <- c(t_n, f_n)             # total trials
J <- length(x)

theta_samples <- matrix(NA, nrow = n_iters, ncol = J)
theta_current <- rep(0.5, J)   # initialize at 0.5 for each group

# MH proposal standard deviation (tune this!)
proposal_sd <- 0.2

# function to compute log-posterior density
log_posterior <- function(theta, x, n, a, b) {
  if (theta <= 0 || theta >= 1) return(-Inf) # invalid
  log_lik <- dbinom(x, size = n, prob = theta, log = TRUE)
  log_prior <- dbeta(theta, a, b, log = TRUE)
  return(log_lik + log_prior)
}

for (i in 1:n_iters) {
  for (j in 1:J) {
    
    # propose new theta via logit-normal random walk
    logit_theta <- qlogis(theta_current[j])                  # current on logit scale
    logit_theta_prop <- rnorm(1, mean = logit_theta, sd = proposal_sd)
    theta_prop <- plogis(logit_theta_prop)                   # back-transform
    
    # log-posterior for current and proposed
    log_post_current <- log_posterior(theta_current[j], x[j], n[j], a_prior[j], b_prior[j])
    log_post_prop <- log_posterior(theta_prop, x[j], n[j], a_prior[j], b_prior[j])
    
    # symmetric proposal, so acceptance ratio simplifies
    log_accept_ratio <- log_post_prop - log_post_current
    
    if (log(runif(1)) < log_accept_ratio) {
      theta_current[j] <- theta_prop  # accept
    }
    
    theta_samples[i, j] <- theta_current[j]
  }
}
