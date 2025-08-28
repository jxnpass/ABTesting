library(tidyverse)
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
#                         FREQUENTIST A/B Testing                           # 
#                      A two-sample t-test for means.                       #
#############################################################################

### Step 1: Hypothesis Test
# We believe pre-study that there is no difference among X groups by Y.
# In other words, the expected difference between X in Y is 0
exp_difference <- 0

# The alternative would be that at least one group has significantly different
# costs than another. This is called a two-way test. We will specify this test with
# a significance level of 0.05.
alpha = .05

### Step 2: Post-study EDA
bw <- 250
n <- nrow(data)

ggplot(data = data, mapping = aes(x = Y, y = after_stat(count), fill = X)) +
  geom_histogram(alpha = .35, binwidth = bw, position = 'identity') +
  geom_density(show.legend = F, alpha = 0, mapping = aes(color = X, y = after_stat(bw*n*density))) + 
  scale_x_continuous(labels = scales::dollar, breaks = seq(0,25000,2500)) +
  # adjust these titles as you see fit
  labs(title = "Is BMI correlated to higher medical expenses?", 
       y = "Frequency", x = "Medical Cost", fill = "Overweight") +
  my_theme

### Step 3: Compute statistics

# We will need the mean, standard deviation, and sample size of each group
group_stats <- data %>% 
  group_by(X) %>% 
  summarise(Mean = mean(Y), 
            Sd = sd(Y),
            N = n())

group_stats

g1 <- group_stats[1,] # FALSE
g2 <- group_stats[2,] # TRUE

# Since each group has a sample greater than N > 30, we can use a
# two-sample t-test for means. 

### Step 4a: Compute t-score and p-value by hand

# We start with computing the observed 
obs_difference <- abs(g1$Mean - g2$Mean)

# ... and if the variances are different, calculate the following
SE <- sqrt((g1$Sd^2 / g1$N) + (g2$Sd^2 / g2$N))

# Conversely, the pooled standard deviation (since sd between groups are close)
df <- g1$N + g2$N - 2
Sp <- sqrt( ((g1$N - 1) * g1$Sd^2 + (g2$N - 1) * g2$Sd^2) / df )
SE <- Sp * sqrt(1/g1$N + 1/g2$N)

# Compute critical value and p-value
t_val <- (obs_difference - exp_difference) / SE
p_val <- 2*pt(q = t_val, df = df, lower.tail = F) # multiply by 2 for two-sided
paste("T-score", round(t_val, 2))
paste("p-value", p_val) # definitely significant 

# Compute confidence interval of difference
lwr = obs_difference - qt(1-alpha/2, df = df) * SE
upr = obs_difference + qt(1-alpha/2, df = df) * SE
paste("Observed Difference:")
paste(obs_difference)
paste("Observed Difference Confidence Interval:")
paste("[", lwr, ", ", upr, "]", sep = "")

### Step 4b: Compute using t.test
g1_df <- data %>% 
  filter(X == F) 
g2_df <- data %>% 
  filter(X == T) 

t.test(x = g1_df$Y,
       y = g2_df$Y,
       alternative = "two.sided",
       mu = exp_difference, 
       var.equal = T,
       conf.level = 1-alpha)

### Step 5: Visualize and compute using bootstrapping methods
# Note that A/B tests incorporates assumptions regarding the sampling distribution
# of each population. Under frequentist assumptions, this sample of data is assumed
# to be one of an infinite number of possible results. While we can only observe what
# our data shows, we can infer what the sampling distribution would look like using
# bootstrap samples for alternative (observed) and null results.

# Set up the bootstrap simulation
nsamps <- 5000
g1_idx <- (1:nrow(data))[data$X == F]
g2_idx <- (1:nrow(data))[data$X == T]

# Compute bootstrap samples
sample_sizes <- c(3, 10, 30, 100, 500, nrow(data)) # for showing how sample size affects sampling distribution
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
    g1_null <- rnorm(nsize, mean = 0, sd = g1$Sd)
    g2_null <- rnorm(nsize, mean = 0, sd = g2$Sd)
    
    # compute and save bs means
    gd_means[sim] <- mean(g1_bs) - mean(g2_bs)
    gd_null[sim] <- mean(g1_null) - mean(g2_null)
    gd_p[sim] <- t.test(g1_bs, g2_bs, alternative = "two.sided",
                        mu = exp_difference, var.equal = T, conf.level = 1-alpha)$p.value
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
       x = "Sample Means", y= "Density", fill = "Type") +
  my_theme

  

