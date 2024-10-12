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
## 1) Setting the working directory ----
## _____________________________________________________________________________

# NOTE: This directry should point to the cloned repository in your local machine.
setwd("/home/l338m483/scratch/Cooccurrence_Inv/R_directory/")

## *****************************************************************************
## 2) Environment and dependencies setup ----
## _____________________________________________________________________________

env_setup <- function()
{
  required_dirs <- "Results"

  if (!dir.exists(required_dirs)) 
  {
    dir.create(required_dirs)
  }
  
  # Load the required libraries
  required_packages <- c(
    "dplyr", 
    "ggplot2", 
    "tidyr", 
    "corrplot",
    "vegan",
    "MASS"
  )
  
  invisible(sapply(required_packages, function(pkg) {
    if (!require(pkg, character.only = TRUE)) {
      install.packages(pkg)
      library(pkg, character.only = TRUE)
    }
  }))
  
  # Validate data files (result of main analysis 1 and 2)
  required_files <- c(
    file.path("Results", "goodness_of_fit_individual_effect.csv"),
    file.path("Results", "results_x2_df.rds")
  )
  
  missing_files <- required_files[!file.exists(required_files)]
  if (length(missing_files) > 0) {
    stop(
      "Required data files are missing:\n",
      "You can manually set a work directory using setwd() or use the 'here' package."
    )
  }
  
  # Sourcing the required auxiliary functions
  
  source(file.path("Auxiliary_functions", "significant_X2_counter.R"))
}

## *****************************************************************************
## 3) Main analysis ----
## _____________________________________________________________________________

main <- function()
{
  # Preparing the environment
  env_setup()
  
  # Setting the directory where the graphs are going to be saved
  graphs <- "/home/l338m483/scratch/Cooccurrence_Inv/R_directory/Plots/Ind_vs_Co-occurrence"
  output_dir <- "/home/l338m483/scratch/Cooccurrence_Inv/R_directory/Results"
  
  if (!dir.exists(graphs)) 
  {
    dir.create(graphs)
  }
  
  # Load the and prepare the data
  goodness_of_fit <- read.csv(file.path("Results", 
                                        "goodness_of_fit_individual_effect.csv"))
  omnibus_X2_results <- readRDS(file.path("Results", "results_x2_df.rds"))
  
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
                                               "lm_results_Ind_and_Co-occurrence.csv"))
  
  # Visualize relationships
  long_data <- analysis_data %>%
    pivot_longer(cols = c(mean_p_value, sig_count),
                 names_to = "Metric", values_to = "Value")
  
  pdf(file.path(graphs, "scatterplot_metrics.pdf"), width = 10, height = 10)
  
  ggplot(long_data, aes(x = Value, y = SDV)) +
    geom_point() +
    geom_smooth(method = "lm", se = FALSE) +
    facet_wrap(~ Metric, scales = "free_x") +
    theme_minimal() +
    labs(title = "Relationship between SDV and omnibus X2 test results for inversion dosage independence",
         x = "Metric Value", y = "Segregation Distortion Value (SDV)")
  
  dev.off()
}
  
main()

## <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<< END >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> ##