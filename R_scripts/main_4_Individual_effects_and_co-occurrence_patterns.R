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

# NOTE: This directory should point to the cloned repository in your local machine.
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
    file.path("Results", "results_x2_df.rds"),
    file.path("Results", "inv_and_gene_metadata.csv")
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
  metadata <- read.csv(file.path("inv_and_gene_metadata.csv"))
  
  # Combining the information of both for the record
  ## Function to process metadata and create chromosome lookup
  create_chr_lookup <- function(metadata) {
    # Split the Line column if it contains multiple lines
    metadata_split <- data.frame(
      INV_ID = metadata$INV_ID,
      Chr = metadata$Chr,
      stringsAsFactors = FALSE
    )
    return(setNames(metadata_split$Chr, paste0("Inv_", metadata_split$INV_ID)))
  }
  
  # Function to check if two inversions are on different chromosomes
  are_different_chr <- function(inv1, inv2, chr_lookup) {
    chr1 <- chr_lookup[inv1]
    chr2 <- chr_lookup[inv2]
    return(!is.na(chr1) && !is.na(chr2) && chr1 != chr2)
  }
  
  # Main function to combine data and filter by chromosome
  create_combined_dataframe <- function(goodness_of_fit, omnibus_X2_results, metadata) {
    # Create chromosome lookup
    chr_lookup <- create_chr_lookup(metadata)
    
    # Create empty list to store results
    combined_results <- list()
    row_counter <- 1
    
    # Process each line's data
    for(line_name in names(omnibus_X2_results)) {
      # Get line-specific data
      line_data <- omnibus_X2_results[[line_name]]
      line_gof <- goodness_of_fit[goodness_of_fit$Line == line_name, ]
      
      # Process each row in the X2 results
      for(i in 1:nrow(line_data)) {
        # Get current row
        row <- line_data[i, ]
        
        # Check if inversions are on different chromosomes
        if(!are_different_chr(row$INV_1, row$INV_2, chr_lookup)) {
          next  # Skip this pair if they're on the same chromosome
        }
        
        # Get goodness of fit data for both inversions
        inv1_gof <- line_gof[line_gof$Inv == row$INV_1, ]
        inv2_gof <- line_gof[line_gof$Inv == row$INV_2, ]
        
        # Create a list with all information
        result_list <- list(
          Line = line_name,
          Inv1 = row$INV_1,
          Inv2 = row$INV_2,
          Inv1_Chr = chr_lookup[row$INV_1],
          Inv2_Chr = chr_lookup[row$INV_2],
          X2 = row$X2,
          X2_p_value = row$p,
          # Goodness of fit metrics for Inv1
          Inv1_G = if(nrow(inv1_gof) > 0) inv1_gof$G else NA,
          Inv1_df = if(nrow(inv1_gof) > 0) inv1_gof$df else NA,
          Inv1_p_value = if(nrow(inv1_gof) > 0) inv1_gof$p_value else NA,
          Inv1_p_corrected = if(nrow(inv1_gof) > 0) inv1_gof$p_corrected else NA,
          Inv1_SDV = if(nrow(inv1_gof) > 0) inv1_gof$SDV else NA,
          # Goodness of fit metrics for Inv2
          Inv2_G = if(nrow(inv2_gof) > 0) inv2_gof$G else NA,
          Inv2_df = if(nrow(inv2_gof) > 0) inv2_gof$df else NA,
          Inv2_p_value = if(nrow(inv2_gof) > 0) inv2_gof$p_value else NA,
          Inv2_p_corrected = if(nrow(inv2_gof) > 0) inv2_gof$p_corrected else NA,
          Inv2_SDV = if(nrow(inv2_gof) > 0) inv2_gof$SDV else NA,
          # Deviation metrics
          dev_1r_1c = row$dev_1r_1c,
          dev_2r_1c = row$dev_2r_1c,
          dev_3r_1c = row$dev_3r_1c,
          dev_1r_2c = row$dev_1r_2c,
          dev_2r_2c = row$dev_2r_2c,
          dev_3r_2c = row$dev_3r_2c,
          dev_1r_3c = row$dev_1r_3c,
          dev_2r_3c = row$dev_2r_3c,
          dev_3r_3c = row$dev_3r_3c,
          # Standardized residuals
          sr_1r_1c = row$sr_1r_1c,
          sr_2r_1c = row$sr_2r_1c,
          sr_3r_1c = row$sr_3r_1c,
          sr_1r_2c = row$sr_1r_2c,
          sr_2r_2c = row$sr_2r_2c,
          sr_3r_2c = row$sr_3r_2c,
          sr_1r_3c = row$sr_1r_3c,
          sr_2r_3c = row$sr_2r_3c,
          sr_3r_3c = row$sr_3r_3c,
          # Relative contributions
          rel_cont_1r_1c = row$rel_cont_1r_1c,
          rel_cont_2r_1c = row$rel_cont_2r_1c,
          rel_cont_3r_1c = row$rel_cont_3r_1c,
          rel_cont_1r_2c = row$rel_cont_1r_2c,
          rel_cont_2r_2c = row$rel_cont_2r_2c,
          rel_cont_3r_2c = row$rel_cont_3r_2c,
          rel_cont_1r_3c = row$rel_cont_1r_3c,
          rel_cont_2r_3c = row$rel_cont_2r_3c,
          rel_cont_3r_3c = row$rel_cont_3r_3c
        )
        
        combined_results[[row_counter]] <- result_list
        row_counter <- row_counter + 1
      }
    }
    
    # Convert list to data frame
    if(length(combined_results) > 0) {
      combined_df <- do.call(rbind, lapply(combined_results, data.frame))
      return(combined_df)
    } else {
      return(data.frame())  # Return empty data frame if no pairs found
    }
  }
  
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
}
  
main()

## <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<< END >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> ##