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
    "vegan"
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
    file.path("Results", "final_result_co_occurrence.csv")
  )
  
  missing_files <- required_files[!file.exists(required_files)]
  if (length(missing_files) > 0) {
    stop(
      "Required data files are missing:\n",
      "You can manually set a work directory using setwd() or use the 'here' package."
    )
  }
}

