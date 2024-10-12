## Co_occurrence of inversions in lines of Mimulus guttatus
#@ Complementary script
#@ Authors: Luis J. Madrigal-Roca & John K. Kelly
#@ Date: 2024-05-27

## 1) Loading the required libraries ----

library(foreach)
library(doParallel)
library(ggplot2)
library(CooccurrenceAffinity)
library(jaccard)
library(tidyverse)

## 2) Defining the working directory and a seed ----

set.seed(1998) # For reproducibility
setwd("/home/l338m483/scratch/Cooccurrence_Inv/R_directory/Simulation")

## 1) Simulating the tables ----

n_ind <- 400

# Register a parallel backend
registerDoParallel(cores = detectCores())

# 1.1) Mendelian expected frequencies ----

expected <- c(1/16, 1/8, 1/16, 1/8, 1/4, 1/8, 1/16, 1/8, 1/16)
expected <- matrix(expected, nrow = 3, ncol = 3)

# 1.2) Simulating the observed counts and checking the null hypothesis ----

obs <- sample(x = 1:9, size = n_ind, replace = T, prob = expected) # 1:9 represents all the cells in the table

obs <- table(obs)

# 1.2.1) Creating a function for generating multiple observed tables ----

table_sim_null <-function(rep, n_cells, n_ind, expected)
{
  simulated_m <- foreach(i = 1:rep) %dopar% 
    {
      matrix(tabulate(sample(x = 1:n_cells, size = n_ind, replace = T, 
                             prob = expected)), nrow = 3)
    }
  
  simulated_m <- array(unlist(simulated_m), 
                       dim = c(nrow(simulated_m[[1]]), ncol(simulated_m[[1]]), 
                               length(simulated_m)))
  
  return(simulated_m)
}

null_model <- table_sim_null(100000, 9, n_ind, expected)

# 1.2.2) Checking the null hypothesis ----

submatrix_extractor <- function(matrix, row_to_exclude, col_to_exclude)
{
  ## Extracting the submatrix of interest
  submatrix <- matrix[-row_to_exclude, -col_to_exclude]
  
  return(submatrix)
}

two_vectors_from_submatrix <- function(submatrix)
{
  ## Getting the two vectors of interest (presence/absebce of a state)
  
  total <- sum(submatrix)
  
  intersection <- submatrix[2, 2]
  intersection <- rep(1, intersection)
  
  both_absents <- submatrix[1, 1]
  both_absents <- rep(0, both_absents)
  
  only1_present <- submatrix[1, 2]
  only2_present <- submatrix[2, 1]
  only1_absent <- total - (length(both_absents) + length(intersection) + 
                             only1_present)
  only2_absent <- total - (length(both_absents) + length(intersection) +
                             only2_present)
  
  vector1 <- c(both_absents, intersection, rep(1, only1_present), 
               rep(0, only1_absent))
  vector2 <- c(both_absents, intersection, rep(0, only2_absent), 
               rep(1, only2_present))
  
  return(list(vector1, vector2))
}

p_value_extractor <- function(models, expected)
{
  ## Models is an array, where the third dimension represent the model
  
  p_values_X2 <- foreach(i = 1:dim(models)[3], .combine = 'c') %dopar% 
    {
      chisq.test(models[,,i], p = expected)$p.value
    }
  
  ## Getting the four submatrix of interest
  
  ## Submatrix 1: Both heterozygous contrast (exluding the last row and column)
  
  sub1 <- foreach(i = 1:dim(models)[3]) %dopar% 
    {
      submatrix_extractor(models[,,i], 3, 3)
    }
  
  ## Submatrix 2: Both homozygous contrast (excluding the second row and column)
  
  sub2 <- foreach(i = 1:dim(models)[3]) %dopar% 
    {
      submatrix_extractor(models[,,i], 2, 2)
    }
  
  ## Submatrix 3 and 4: One homozygous and the other hetwrozygous
  
  sub3 <- foreach(i = 1:dim(models)[3]) %dopar% 
    {
      submatrix_extractor(models[,,i], 2, 3)
    }
  
  sub4 <- foreach(i = 1:dim(models)[3]) %dopar% 
    {
      submatrix_extractor(models[,,i], 3, 2)
    }
  
  ## Getting the two vectors per submatrix (each vector is the presence absence)
  
  vectors1 <- foreach(i = 1:dim(models)[3]) %dopar% 
    {
      two_vectors_from_submatrix(sub1[[i]])
    }
  
  vectors2 <- foreach(i = 1:dim(models)[3]) %dopar% 
    {
      two_vectors_from_submatrix(sub2[[i]])
    }
  
  vectors3 <- foreach(i = 1:dim(models)[3]) %dopar% 
    {
      two_vectors_from_submatrix(sub3[[i]])
    }
  
  vectors4 <- foreach(i = 1:dim(models)[3]) %dopar% 
    {
      two_vectors_from_submatrix(sub4[[i]])
    }
  
  ## Getting the p-values for the four contrasts (affinity)
  
  affinity1 <- foreach(i = 1:dim(models)[3], .combine = 'c') %dopar% 
    {
      as.numeric(affinity(matrix(unlist(vectors1[[i]]), nrow = 2, byrow = T), 
                          row.or.col = 'row', datatype = 'binary', 
                          pvalType="midP", sigdigit = 3)$all[, 'p_value'])
    }
  
  affinity2 <- foreach(i = 1:dim(models)[3], .combine = 'c') %dopar% 
    {
      as.numeric(affinity(matrix(unlist(vectors2[[i]]), nrow = 2, byrow = T), 
                          row.or.col = 'row', datatype = 'binary', 
                          pvalType="midP", sigdigit = 3)$all[, 'p_value'])
    }
  
  affinity3 <- foreach(i = 1:dim(models)[3], .combine = 'c') %dopar% 
    {
      as.numeric(affinity(matrix(unlist(vectors3[[i]]), nrow = 2, byrow = T), 
                          row.or.col = 'row', datatype = 'binary', 
                          pvalType="midP", sigdigit = 3)$all[, 'p_value'])
    }
  
  affinity4 <- foreach(i = 1:dim(models)[3], .combine = 'c') %dopar% 
    {
      as.numeric(affinity(matrix(unlist(vectors4[[i]]), nrow = 2, byrow = T), 
                          row.or.col = 'row', datatype = 'binary', 
                          pvalType="midP", sigdigit = 3)$all[, 'p_value'])
    }
  
  p_values_affinity <- list(sub1 = affinity1, 
                            sub2 = affinity2, 
                            sub3 = affinity3, 
                            sub4 = affinity4)
  
  ## Getting the p-values for the four contrasts (jaccard)
  
  jaccard1 <- foreach(i = 1:dim(models)[3], .combine = 'c') %dopar% 
    {
      jaccard.test(vectors1[[i]][[1]], vectors1[[i]][[2]], method = "bootstrap",
              B = 10000, verbose = F)[['pvalue']]
    }
  
  jaccard2 <- foreach(i = 1:dim(models)[3], .combine = 'c') %dopar% 
    {
      jaccard.test(vectors2[[i]][[1]], vectors2[[i]][[2]], method = "bootstrap",
              B = 10000, verbose = F)[['pvalue']]
    }
  
  jaccard3 <- foreach(i = 1:dim(models)[3], .combine = 'c') %dopar% 
    {
      jaccard.test(vectors3[[i]][[1]], vectors3[[i]][[2]], method = "bootstrap",
              B = 10000, verbose = F)[['pvalue']]
    }
  
  jaccard4 <- foreach(i = 1:dim(models)[3], .combine = 'c') %dopar% 
    {
      jaccard.test(vectors4[[i]][[1]], vectors4[[i]][[2]], method = "bootstrap",
              B = 10000, verbose = F)[['pvalue']]
    }
  
  p_values_jaccard <- list(sub1 = jaccard1, 
                           sub2 = jaccard2, 
                           sub3 = jaccard3, 
                           sub4 = jaccard4)
  
  return(list(X2 = p_values_X2, affinity = p_values_affinity,
              jaccard = p_values_jaccard))
}

p_values_null_model <- p_value_extractor(null_model, expected)

# 1.2.3) Plotting the p_values associated to the null model ----

p_df <- data.frame(
  p_values_X2 = p_values_null_model[[1]],
  p_values_affinity_sub1 = p_values_null_model[[2]][[1]],
  p_values_affinity_sub2 = p_values_null_model[[2]][[2]],
  p_values_affinity_sub3 = p_values_null_model[[2]][[3]],
  p_values_affinity_sub4 = p_values_null_model[[2]][[4]],
  p_values_jaccard_sub1 = p_values_null_model[[3]][[1]],
  p_values_jaccard_sub2 = p_values_null_model[[3]][[2]],
  p_values_jaccard_sub3 = p_values_null_model[[3]][[3]],
  p_values_jaccard_sub4 = p_values_null_model[[3]][[4]]
)

p_df <- p_df |>
  pivot_longer(cols = everything(), names_to = "Test", values_to = "P-values")

# Calculate the mean for each group

mean_df <- p_df %>%
  group_by(Test) %>%
  summarise(mean_p_value = mean(`P-values`))

## Histograms

pdf("p_values_null_model.pdf", 10, 8)

p <- ggplot(data = p_df, 
            aes(x = `P-values`, fill = Test)) +
  geom_histogram(aes(y = after_stat(count/100000)),
                 breaks = seq(0, 1, by = 0.05),  # Manually set breaks
                 alpha = .5) +
  geom_vline(data = mean_df, aes(xintercept = mean_p_value), 
             color = "blue", linetype = "dashed") +
  facet_wrap(~ Test) +
  theme_bw() +
  labs(x = "P-values",
       y = "Density") +
  theme(axis.text = element_text(size = 12, color = 'black'),
        axis.title = element_text(size = 14, color = 'black'),
        legend.text = element_text(size = 12, color = 'black'),
        legend.title = element_text(size = 14, color = 'black'),
        panel.grid = element_blank())

print(p)

dev.off()

pdf("p_values_null_model_log10.pdf", 10, 8)

p <- ggplot(data = p_df, 
            aes(x = `P-values`, fill = Test)) +
  geom_histogram(aes(y = after_stat(count/100000)),
                 bins = 20,
                 boundary = 0,
                 alpha = .5) +
  geom_vline(data = mean_df, aes(xintercept = mean_p_value), 
             color = "blue", linetype = "dashed") +
  facet_wrap(~ Test) +
  theme_bw() +
  labs(x = "P-values_log10",
       y = "Density") +
  scale_x_log10() +
  theme(axis.text = element_text(size = 12, color = 'black'),
        axis.title = element_text(size = 14, color = 'black'),
        legend.text = element_text(size = 12, color = 'black'),
        legend.title = element_text(size = 14, color = 'black'),
        panel.grid = element_blank())

print(p)

dev.off()

# CDF plot

# Create a new data frame for the CDF plot
cdf_df <- p_df %>%
  arrange(`P-values`) %>%
  group_by(Test) %>%
  mutate(cdf = row_number() / n())

# Create the CDF plot
pdf("p_values_null_model_cdf_plot.pdf", 10, 8)

p <- ggplot(cdf_df, aes(x = `P-values`, y = cdf, color = Test)) +
  geom_line() +
  facet_wrap(~ Test) +
  theme_bw() +
  labs(x = "P-values",
       y = "Cumulative Distribution",
       color = "Test") +
  theme(axis.text = element_text(size = 12, color = 'black'),
        axis.title = element_text(size = 14, color = 'black'),
        legend.text = element_text(size = 12, color = 'black'),
        legend.title = element_text(size = 14, color = 'black'),
        panel.grid = element_blank())

print(p)

dev.off()

# Checking that the p values are uniformly distributed ----

uniform_ref <- runif(100000)

ks_res <-list()

for (level in unique(p_df$Test))
{
  ks_res[[level]] <- ks.test(p_df$`P-values`[p_df$Test == level], uniform_ref)
}
  
save.image('Env')
print('Checkpoint 1')

## 2) Simulating the tables with a specific pattern ----

# 2.1) Survival per cell ----

survival <- c(rep(1, 8), 0.2) # 9th cell deflated

survival <- matrix(survival, nrow = 3, ncol = 3)

expected_mod <- function(expected, survival)
{
  return(expected * survival / sum(expected * survival))
}

# 2.2) Simulating the tables the 9th cell deflated ----

expected_9_deflated <- expected_mod(expected, survival)

model_9_deflated <- table_sim_null(100000, 9, n_ind, expected_9_deflated)

# 2.3) Checking the p-values for the 9th cell deflated model ----

p_values_9_deflated_model <- p_value_extractor(model_9_deflated, expected)

# Plotting the p_values associated to the 9th cell deflated model ----

p_df_9_deflated <- data.frame(
  p_values_X2 = p_values_9_deflated_model[[1]],
  p_values_affinity_sub1 = p_values_9_deflated_model[[2]][[1]],
  p_values_affinity_sub2 = p_values_9_deflated_model[[2]][[2]],
  p_values_affinity_sub3 = p_values_9_deflated_model[[2]][[3]],
  p_values_affinity_sub4 = p_values_9_deflated_model[[2]][[4]],
  p_values_jaccard_sub1 = p_values_9_deflated_model[[3]][[1]],
  p_values_jaccard_sub2 = p_values_9_deflated_model[[3]][[2]],
  p_values_jaccard_sub3 = p_values_9_deflated_model[[3]][[3]],
  p_values_jaccard_sub4 = p_values_9_deflated_model[[3]][[4]]
)

p_df_9_deflated <- p_df_9_deflated |>
  pivot_longer(cols = everything(), names_to = "Test", values_to = "P-values")

# Calculate the mean for each group

mean_df_9_deflated <- p_df_9_deflated %>%
  group_by(Test) %>%
  summarise(mean_p_value = mean(`P-values`)) |>
  mutate(mean_p_value_log10 = log10(mean_p_value))

## Histograms

pdf("p_values_9_deflated_model.pdf", 10, 8)

p <- ggplot(data = p_df_9_deflated, 
            aes(x = `P-values`, fill = Test)) +
  geom_histogram(aes(y = after_stat(count/100000)),
                 breaks = seq(0, 1, by = 0.05),  # Manually set breaks
                 alpha = .5) +
  geom_vline(data = mean_df_9_deflated, aes(xintercept = mean_p_value), 
             color = "blue", linetype = "dashed") +
  facet_wrap(~ Test) +
  theme_bw() +
  labs(x = "P-values",
       y = "Density") +
  scale_x_continuous(limits = c(0, 1)) +
  theme(axis.text = element_text(size = 12, color = 'black'),
        axis.title = element_text(size = 14, color = 'black'),
        legend.text = element_text(size = 12, color = 'black'),
        legend.title = element_text(size = 14, color = 'black'),
        panel.grid = element_blank())

print(p)

dev.off()

# Log transformed histogram

pdf("p_values_9_deflated_model_log10.pdf", 10, 8)

p <- ggplot(data = p_df_9_deflated, 
            aes(x = `P-values`, fill = Test)) +
  geom_histogram(aes(y = after_stat(count/100000)),
                 bins = 20,
                 boundary = 0,
                 alpha = .5) +
  geom_vline(data = mean_df_9_deflated, aes(xintercept = mean_p_value), 
             color = "blue", linetype = "dashed") +
  facet_wrap(~ Test) +
  theme_bw() +
  labs(x = "P-values_log10",
       y = "Density") +
  scale_x_log10() +
  theme(axis.text = element_text(size = 12, color = 'black'),
        axis.title = element_text(size = 14, color = 'black'),
        legend.text = element_text(size = 12, color = 'black'),
        legend.title = element_text(size = 14, color = 'black'),
        panel.grid = element_blank())

print(p)

dev.off()

# CDF plot

# Create a new data frame for the CDF plot
cdf_df_9_deflated <- p_df_9_deflated %>%
  arrange(`P-values`) %>%
  group_by(Test) %>%
  mutate(cdf = row_number() / n())

# Create the CDF plot
pdf("p_values_9_deflated_model_cdf_plot.pdf", 10, 8)

p <- ggplot(cdf_df_9_deflated, aes(x = `P-values`, y = cdf, color = Test)) +
  geom_line() +
  facet_wrap(~ Test) +
  theme_bw() +
  labs(x = "P-values",
       y = "Cumulative Distribution",
       color = "Test") +
  theme(axis.text = element_text(size = 12, color = 'black'),
        axis.title = element_text(size = 14, color = 'black'),
        legend.text = element_text(size = 12, color = 'black'),
        legend.title = element_text(size = 14, color = 'black'),
        panel.grid = element_blank())

print(p)

dev.off()

# Box-plot for the unstransformed p-values (Zoom-in, 0 - 0.01)

pdf("boxplot_9_deflated_model.pdf", 10, 8)

p <- ggplot(data = p_df_9_deflated,
            aes(x = `P-values`, y = Test, fill = Test)) +
  geom_violin(alpha = .3) +
  geom_boxplot() +
  theme_bw() +
  labs(x = "P-values",
       y = "Test") +
  scale_x_continuous(limits = c(0, 0.01)) +
  theme(axis.text = element_text(size = 12, color = 'black'),
        axis.title = element_text(size = 14, color = 'black'),
        legend.text = element_text(size = 12, color = 'black'),
        legend.title = element_text(size = 14, color = 'black'),
        panel.grid = element_blank())

print(p)

dev.off()

# 3) Simulating a deficiency in the 8th cell and a inflation in the 9th one ----

# Position is based on the native filling order of a matrix in R (column-wise)

# 3.1) Survival per cell ----

survival_2 <- c(rep(0.7, 7), 0.2 , 1) # 8th cell deflated and 9th inflated

survival_2 <- matrix(survival_2, nrow = 3, ncol = 3)

# 2.2) Simulating the tables the 9th cell deflated ----

expected_9and8 <- expected_mod(expected, survival_2)

model_9and8  <- table_sim_null(100000, 9, n_ind, expected_9and8)

# 2.3) Checking the p-values for the 9th cell deflated model ----

p_values_9and8 <- p_value_extractor(model_9and8, expected)

# Plotting the p_values associated to the 9th cell deflated model ----

p_df_9and8 <- data.frame(
  p_values_X2 = p_values_9and8[[1]],
  p_values_affinity_sub1 = p_values_9and8[[2]][[1]],
  p_values_affinity_sub2 = p_values_9and8[[2]][[2]],
  p_values_affinity_sub3 = p_values_9and8[[2]][[3]],
  p_values_affinity_sub4 = p_values_9and8[[2]][[4]],
  p_values_jaccard_sub1 = p_values_9and8[[3]][[1]],
  p_values_jaccard_sub2 = p_values_9and8[[3]][[2]],
  p_values_jaccard_sub3 = p_values_9and8[[3]][[3]],
  p_values_jaccard_sub4 = p_values_9and8[[3]][[4]]
)

p_df_9and8 <- p_df_9and8 |>
  pivot_longer(cols = everything(), names_to = "Test", values_to = "P-values")

# Calculate the mean for each group

mean_df_9and8 <- p_df_9and8 %>%
  group_by(Test) %>%
  summarise(mean_p_value = mean(`P-values`)) |>
  mutate(mean_p_value_log10 = log10(mean_p_value))

## Histograms

pdf("p_values_9and8_model.pdf", 10, 8)

p <- ggplot(data = p_df_9and8, 
            aes(x = `P-values`, fill = Test)) +
  geom_histogram(aes(y = after_stat(count/100000)),
                 breaks = seq(0, 1, by = 0.05),  # Manually set breaks
                 alpha = .5) +
  geom_vline(data = mean_df_9and8, aes(xintercept = mean_p_value), 
             color = "blue", linetype = "dashed") +
  facet_wrap(~ Test) +
  theme_bw() +
  labs(x = "P-values",
       y = "Density") +
  scale_x_continuous(limits = c(0, 1)) +
  theme(axis.text = element_text(size = 12, color = 'black'),
        axis.title = element_text(size = 14, color = 'black'),
        legend.text = element_text(size = 12, color = 'black'),
        legend.title = element_text(size = 14, color = 'black'),
        panel.grid = element_blank())

print(p)

dev.off()

# Log transformed histogram

pdf("p_values_9and8_model_log10.pdf", 10, 8)

p <- ggplot(data = p_df_9and8, 
            aes(x = `P-values`, fill = Test)) +
  geom_histogram(aes(y = after_stat(count/100000)),
                 bins = 20,
                 boundary = 0,
                 alpha = .5) +
  geom_vline(data = mean_df_9and8, aes(xintercept = mean_p_value), 
             color = "blue", linetype = "dashed") +
  facet_wrap(~ Test) +
  theme_bw() +
  labs(x = "P-values_log10",
       y = "Density") +
  scale_x_log10() +
  theme(axis.text = element_text(size = 12, color = 'black'),
        axis.title = element_text(size = 14, color = 'black'),
        legend.text = element_text(size = 12, color = 'black'),
        legend.title = element_text(size = 14, color = 'black'),
        panel.grid = element_blank())

print(p)

dev.off()

# CDF plot

# Create a new data frame for the CDF plot
cdf_df_9and8 <- p_df_9and8 %>%
  arrange(`P-values`) %>%
  group_by(Test) %>%
  mutate(cdf = row_number() / n())

# Create the CDF plot
pdf("p_values_9and8_model_cdf_plot.pdf", 10, 8)

p <- ggplot(cdf_df_9and8, aes(x = `P-values`, y = cdf, color = Test)) +
  geom_line() +
  facet_wrap(~ Test) +
  theme_bw() +
  labs(x = "P-values",
       y = "Cumulative Distribution",
       color = "Test") +
  theme(axis.text = element_text(size = 12, color = 'black'),
        axis.title = element_text(size = 14, color = 'black'),
        legend.text = element_text(size = 12, color = 'black'),
        legend.title = element_text(size = 14, color = 'black'),
        panel.grid = element_blank())

print(p)

dev.off()

save.image('Env')