library(tidyverse)
library(tidymodels)
source("theme.R")

# Loading in Data
set.seed(39)
data <- read_csv("payment_fraud.csv") %>% 
  rename(X = numItems, Y = label) %>% 
  select(X, Y)

# Change to binary for groups if needed
data$Y <- ifelse(data$Y == 1, T, F) # is fraud
data$X <- ifelse(data$X > 1, T, F) # number of purchases > 1

#############################################################################
#                         FREQUENTIST A/B Testing                           # 
#                    A two-sample z-test for proportions.                   #
#############################################################################

### Proportions test on variable Y ###

### Step 1: Hypothesis Test
# We believe pre-study that there is no difference among groups.
exp_diff <- 0 - 0

# The alternative would be that at least one group has greater exhaustion
# than another. This is called a two-way test. We will specify this test with
# a significance level of 0.05.
## NOTE: THIS IS NOT THE SAME AS HYPOTHESIZING THAT EXTROVERTS ARE LESS DRAINED!
alpha = .05

### Step 2: Post-study EDA
eda <- data %>% 
  group_by(X, Y) %>% 
  summarise(N = n()) %>% 
  mutate(Proportion = N/sum(N)) 

ggplot(data = eda %>% filter(Y == T), mapping = aes(x = X, 
                                 y = Proportion)) +
  geom_col(alpha = .5, color = "black", position = 'dodge' ) +
  scale_y_continuous(labels = scales::percent, limits = c(0,.08)) +
  geom_text(mapping = aes(label = scales::percent(Proportion)), 
            size = 5, vjust = -.5,
            position = position_dodge(width = .9)) +
  # Adjust titles as you see fit 
  labs(title = "Is a multiple-item purchase on a credit card a signal for fraud?", 
       y = "Fraud Proportion", 
       x = "Was there more than one purchase?") +
  my_theme

### Step 3: Compute statistics

# We will need the mean, standard deviation, and sample size of each group

group_stats <- data %>% 
  group_by(X) %>% 
  summarise(Prop = mean(Y), 
            Sum = sum(Y),
            Sd = sd(Y),
            N = n())

group_stats

g1 <- group_stats[2,] # TRUE
g2 <- group_stats[1,] # FALSE

# Can also verify assumptions
g1$N * g1$Prop > 5
g1$N * (1-g1$Prop) > 5
g2$N * g2$Prop > 5
g2$N * (1-g2$Prop) > 5

# The main assumption is that for each group X and response, 
# there should be at least five respondents in each category. 
# If this is met, proceed with the test 

### Step 4a: Compute z-score and p-value by hand

# We start with computing the observed 
obs_diff <- abs(g1$Prop - g2$Prop)

# Compute unpooled SE
x1 <- g1$Sum
x2 <- g2$Sum
n1 <- g1$N
n2 <- g2$N
p1 <- g1$Prop
p2 <- g2$Prop
SE_u <- sqrt(((p1) * (1 - p1)) / n1 + ((p2) * (1 - p2)) / n2)

# ... or compute pooled sample proportion (p_hat) and standard error
p_hat <- (g1$Sum + g2$Sum) / (n1 + n2)
SE_p <- sqrt(p_hat * (1-p_hat) * (1/n1 + 1/n2))

# Compute critical value and p-value
z_val <- (obs_diff - exp_diff) / SE_p
p_val <- 2*pnorm(q = abs(z_val), lower.tail = F) # multiply by 2 for two-sided
paste("Z-score", round(z_val, 2))
paste("p-value", p_val)  
paste("stat sig", p_val < alpha) #  

# Compute confidence interval of difference
lwr = obs_diff - qnorm(1-alpha/2) * SE_u
upr = obs_diff + qnorm(1-alpha/2) * SE_u
paste("Observed Difference")
paste(obs_difference)
paste("Observed Difference Confidence Interval:")
paste("[", lwr, ", ", upr, "]", sep = "")

### Step 4b: Can also use the prop.test() function 

prop.test(x = c(x1, x2), n = c(n1, n2), 
          alternative = "two.sided",
          correct = F)
# note X-squared = z_val^2

### Step 4c: Can also use tidymodels 
# clean and works directly with your dataframe

data %>% 
  prop_test(Y ~ X, order = c("TRUE", "FALSE"), 
            correct = F, z = TRUE)

### Step 5: Visualize and compute using bootstrapping methods
# Note that A/B tests incorporates assumptions regarding the sampling distribution
# of each population. Under frequentist assumptions, this sample of data is assumed
# to be one of an infinite number of possible results. While we can only observe what
# our data shows, we can infer what the sampling distribution would look like using
# bootstrap samples for alternative (observed) and null results.

# Set up the bootstrap simulation
nsamps <- 1000
g1_idx <- (1:nrow(data))[data$X == T]
g2_idx <- (1:nrow(data))[data$X == F]

# Compute bootstrap samples
sample_sizes <- c(30, 50, 100, 500, 1000, 3000) # for showing how sample size affects sampling distribution
obs_dist_results <- data.frame(matrix(NA, nrow = nsamps, ncol = length(sample_sizes)))
colnames(obs_dist_results) <- sample_sizes
null_dist_results <- obs_dist_results
p_value_results <- obs_dist_results[1,]

col <- 0
for (nsize in sample_sizes) {
  col <- col + 1
  gd_means <- numeric(nsamps)
  gd_null <- numeric(nsamps)
  gd_p <- numeric(nsamps)
  for (sim in 1:nsamps) {
    # which obs are in the bs
    g1_bs_idx <- sample(g1_idx, nsize, replace = T)
    g2_bs_idx <- sample(g2_idx, nsize, replace = T)
    g1_bs <- data$Y[g1_bs_idx]
    g2_bs <- data$Y[g2_bs_idx]
    
    # compute null distribution
    g1_null <- rbernoulli(nsize, p = .01)
    g2_null <- rbernoulli(nsize, p = .01)
    
    # compute and save bs means
    gd_means[sim] <- mean(g1_bs) - mean(g2_bs)
    gd_null[sim] <- mean(g1_null) - mean(g2_null)
    gd_p[sim] <- prop.test(x = c(sum(g1_bs), sum(g2_bs)), 
                           n = c(length(g1_bs), length(g2_bs)), 
                           alternative = "two.sided")$p.value
    gd_p[sim] <- ifelse(is.na(gd_p[sim]), 1, gd_p[sim])
  }
  obs_dist_results[col] <- gd_means
  null_dist_results[col] <- gd_null
  p_value_results[col] <- mean(gd_p)
}

# Visualize bootstrapped sampling distribution
obs_dist_results$type <- "Observed"
null_dist_results$type <- "Null"

p_value_format <- p_value_results %>% 
  pivot_longer(cols = everything()) %>% 
  rename(N_num = name) %>% 
  mutate(N_num = as.numeric(N_num)) %>% 
  mutate(p_value = ifelse(abs(value) >= 1e4 | abs(value) < 1e-4,
                          formatC(value, format = "e", digits = 1),
                          formatC(value, format = "f", digits = 4))) %>% 
  select(N_num, p_value)


gd_df <- rbind(obs_dist_results, null_dist_results) %>% 
  pivot_longer(cols = -type) %>% 
  mutate(
    N_num = as.numeric(name),  # Extract numeric value
    name = factor(
      paste("N =", name),
      levels = paste("N =", sort(unique(as.numeric(name))))
    )
  ) %>%
  left_join(p_value_format, by = 'N_num') %>%
  # Keep original order, just append p-values
  mutate(
    name = factor(
      paste0(as.character(name), " (p = ", p_value, ")"),
      levels = paste0(levels(name), " (p = ", p_value[match(levels(name), as.character(name))], ")")
    )
  )

ggplot(data = gd_df) +
  geom_density(mapping = aes(x = value, fill = type), 
               color = "black", alpha = .25) +
  geom_histogram(mapping = aes(x = value, y = after_stat(density), fill = type), 
                 bins = 50, color = "black", alpha = .25, position = 'identity') +
  facet_wrap(~ name, nrow = 2, scales = 'free_y') + 
  labs(title = "Bootstrapped Sampling Distributions",
       subtitle = "Each value represents a difference in sample means of size N",
       x = "Sample Proportion", y= "Density", fill = "Type") +
  my_theme


### Step 5b: can also show how p-value is associated with the null distribution


# 1. Compute observed difference in proportions
point_estimate <- data %>%
  specify(Y ~ X, success = 'TRUE') %>%
  calculate(stat = "diff in props", order = c("TRUE", "FALSE"))

obs_stat <- point_estimate$stat

# 2. Generate null distribution under independence
null_dist <- data %>%
  specify(Y ~ X, success = 'TRUE') %>%
  hypothesize(null = "independence") %>%
  generate(reps = 500, type = "permute") %>%
  calculate(stat = "diff in props", order = c("TRUE", "FALSE"))

# 3. Visualize null and shade p-value
null_dist %>%
  visualize() +
  shade_p_value(obs_stat = obs_stat, direction = "two_sided")


