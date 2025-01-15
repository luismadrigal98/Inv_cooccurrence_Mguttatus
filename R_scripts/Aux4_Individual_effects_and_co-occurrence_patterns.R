##' @title Relation between the individual effects of the inversions and the
##' pattern of co-occurrence of the inversions.
##' 
##' @description This analysis will explore whether an inversion showing a significant
##' deviation from Mendelian segregation tends to co-occur or repel other inversions
##' more frequently than those that do not show a significant deviation.
##' 
##' @note This analysis is exploratory in nature. The logarithm of the raw p-values from the
##' co-occurrence analysis will be used to measure the strength of association. Given the large
##' number of tests performed (around 6000) and the small sample sizes, applying a correction
##' for multiple testing, such as the False Discovery Rate (FDR), may result in no significant
##' findings. However, the log-transformed p-values can still provide a measure of the strength
##' of the association, interpreted as the evidence against the null hypothesis. These findings
##' should be validated in a confirmatory analysis to account for multiple testing and reduce the
##' risk of false positives.
##' ____________________________________________________________________________

## >>>>>>>>>>>>>>>>>>>>>>>>>>>>>> START <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<

## *****************************************************************************
## 1) Relationship analysis ----
## _____________________________________________________________________________

# Setting the directory where the graphs are going to be saved
graphs <- "./Results/Plots"
output_dir <- "./Results/"

if (!dir.exists(graphs)) 
{
  dir.create(graphs)
}

# Load the and prepare the data
goodness_of_fit <- read.csv(file.path("Results", 
                                      "goodness_of_fit_individual_effect.csv"))
omnibus_X2_results <- readRDS(file.path("Results", "results_x2_df.rds"))
metadata <- read.csv(file.path("Data/inv_and_gene_metadata.csv"))

# Combining the information of both for the record

combined_data <- create_combined_dataframe(goodness_of_fit, omnibus_X2_results, 
                                           metadata)

write.csv(combined_data, file = file.path(output_dir, 
                                          "combined_data_Ind_and_Co-occurrence_X2.csv"))

# Summarize co-occurrence patterns for each inversion ----
inversion_summary_X2_summary <- 
  mapply(significant_X2_counter, omnibus_X2_results, 
         names(omnibus_X2_results), SIMPLIFY = FALSE)

inversion_summary_X2_summary <- do.call(rbind, inversion_summary_X2_summary)

# Merge with SDV data
analysis_data <- goodness_of_fit %>%
  dplyr::select(Line, Inv, SDV) %>%
  inner_join(inversion_summary_X2_summary, by = c("Line", "Inv"))

# Correlation analysis
cor_matrix <- cor(analysis_data %>% dplyr::select(-Line, -Inv), 
                  method = "spearman")

# Plotting the correlation matrix

pdf(file.path(graphs, "correlation_matrix.pdf"), width = 10, height = 10)

corrplot(cor_matrix, method = "color", type = "upper", order = "hclust", 
         tl.col = "black", tl.srt = 45)

dev.off()

# GLM
# It could be appropriate Poisson or the Negative Binomial (nevertheles, there is
# over dispersion in the data, so, negative binomial is the best choice)

quasi_poisson_model <- glm(sig_count ~ SDV + Line, 
                           data = analysis_data, family = quasipoisson())

out <- summary(quasi_poisson_model)
print(out)

write.csv(out$coefficients, file = file.path(output_dir, 
                                             "glm_results_Ind_and_Co-occurrence.csv"))

# Visualize relationships
long_data <- analysis_data %>%
  pivot_longer(cols = c(mean_p_value, sig_count),
               names_to = "Metric", values_to = "Value")

pdf(file.path(graphs, "scatterplot_metrics.pdf"), width = 10, height = 10)
plot_2 <- ggplot(long_data, aes(x = Value, y = SDV)) +
  geom_point() +
  stat_smooth(method = "glm", method.args = list(family = quasipoisson()), se = FALSE) +
  facet_wrap(~ Metric, scales = "free_x") +
  theme_minimal() +
  labs(title = "Relationship between SDV and omnibus X2 test results for inversion dosage independence",
       x = "Metric Value", y = "Segregation Distortion Value (SDV)")

print(plot_2)
dev.off()

## <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<< END >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> ##