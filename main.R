##' @title Are you with me? Co-occurrence tests from community ecology can 
##' identify positive and negative epistasis between inversions in Mimulus 
##' guttatus
##' 
##' @description This script will source and execute the different analysis employed
##' in this study, and will generate the figures in the Results subdirectory.
##' 
##' @author Luis Javier Madrigal-Roca & John K. Kelly
##' 
##' @date 2024-11-03
##' 
##' @note It is important to notice that the script in the current con figuration
##' was executed in a high performance cluster. Most of the individual analysis
##' does not require a lot of resources, but in case of execution errors, be sure
##' to adjust or eliminate the parallelization in the code (it can allocate a lot
##' of resources during the simulation step)
##' 
##' >>>>>>>>>>>>>>>>>>>>>>>>>>> IMPORTANT <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
##' If an installation from a package using CRAN fails, try the github repository
##' using remotes. 
##' ____________________________________________________________________________

## *****************************************************************************
## 1) Setting up the environment, the seed, and the libraries ----
## _____________________________________________________________________________

## 1.1) Setting the working directory ----

setwd("./") ## Set to the root directory of the project. All source commands are
            ## relative to this directory. Outputs are also relative to this pointer.

## 1.2) Loading the required libraries
required_libraries <- c("dplyr", "ggplot2", "foreach", "doParallel",
                        "CooccurrenceAffinity", "jaccard", "tidyverse",
                        "parallel", "DescTools", "readxl", "reshape",
                        "gplots", "igraph", "tidygraph", "ggraph",
                        "reshape2", "Cairo", "vegan", "rcompanion",
                        "chisq.posthoc.test", "MASS", "repmod",
                        "gtools", "viridis", "car", "nortest",
                        "tidyr", "corrplot", "genetics", "copula",
                        "gridExtra", "mgcv", "qvalue")

## 1.3) Sourcing the setup function ----
source("src/set_environment.R")

set_environment(required_libraries, 
                personal_seed = 1998, 
                parallel_backend = TRUE)

## *****************************************************************************
## 2) Sourcing all auxiliary functions ----
## _____________________________________________________________________________

source("src/pull_aux_functions.R")
pull_aux_functions("./src")

## *****************************************************************************
## 3) Importing the data ----
## _____________________________________________________________________________

# Creating a list with one slot per genetic family or cross
Data <- list()

for(file in list.files('Data/CSVs', 
                       full.names = T))
{
  # Getting the name of the cross or genetic family
  name <- strsplit(x = file, split = '/')[[1]][3]
  name <- strsplit(x = name, split = '_')[[1]][4]
  name <- name <- sub("\\..*$", "", name)
  
  Data[[paste0("Data_", name)]] <- read.csv(file, header = T)
  
  rm(file)
  rm(name)
}

# Transforming the dosage information into presence / absence
Data_splitted <- lapply(X = Data, FUN = dosage_splitter)

# Create an empty list to store the matrices of presence/absence
Data_p_a <- list()

# Iterate over each object in the current R environment
for (object in names(Data_splitted)) 
{
  # Transform the data frame
  new_matrix <- as.matrix(t(Data_splitted[[object]][, -1]))
  colnames(new_matrix) <- Data_splitted[[object]]$Probes
  
  # Add the transformed matrix to the list
  Data_p_a[[paste0('L_', strsplit(object, "_")[[1]][2])]] <- new_matrix
  
  rm(new_matrix)
  rm(object)
}

# Create an empty list to store the matrices of dosage levels

Data_d_l <- list()

for (object in names(Data))
{
  # Transform the data frame
  new_matrix <- as.matrix(t(Data[[object]][, -1]))
  colnames(new_matrix) <- Data[[object]]$Probes
  
  # Add the transformed matrix to the list
  Data_d_l[[paste0('L_', strsplit(object, "_")[[1]][2])]] <- new_matrix
  
  rm(new_matrix)
  rm(object)
}

## Importing the information related to position and chromosome of the INV

metadata <- read.csv("Data/inv_and_gene_metadata.csv", 
                     header = T)

# 3.1) Checking the data ----

#' Additional checking step for compliance in the marginals.
#' If there is an inversion with a dosage level of 0, it will be removed.
#' Its effect will be captured in the individual effect on survival analysis

# Searching for problematic in the data

problematic_pairs <- lapply(Data_d_l, marginal_compliance_checking)

print(problematic_pairs)

#' There is a particularly interesting case in the data. Inv_49, in the line 
#' 444 does not exhibit the homozygous state for the inversion (dosage 2). It is
#' pertinent then to exclude that inversion from posterior analysis related to
#' co-occurrence to avoid problems with marginals 0.

## *****************************************************************************
## 4) Simulation analysis ----
## _____________________________________________________________________________

#' This step is designed to show the performance of the 2x2 and 3x3 table-based
#' frameworks. Essentially, it shows the advantages of using a more granular 
#' approach (2x2 tables generated from the omnibus 3x3 table) through the co-
#' occurrence indexes.

source("R_scripts/Aux1_Contingency_tables_simulation.R")

## *****************************************************************************
## 5) Segregation distortion analysis (SD) ----
## _____________________________________________________________________________

#' This step is designed to test for segregation distortion of the individual
#' inversions. SD can be understood as the individual effect of each inversion in
#' terms of survival. It is an analysis based on deviations from the Mendelian
#' expectations.

source("R_scripts/Aux2_Segregation_distortion_analysis.R")

## *****************************************************************************
## 6) Co-occurrence analysis ----
## _____________________________________________________________________________

#' This is the main step of the paper. It will explore the pairwise interactome
#' between the inversions identified by exploiting co-occurrence indexes
#' traditionally used in ecological studies. Also, network visualizations will
#' be employed to present the results in a more intuitive manner.

source("R_scripts/Aux3_Co_occurrence_per_line.R")

## *****************************************************************************
## 7) SDV and co-occurrence metrics' relation ----
## _____________________________________________________________________________



source("R_scripts/Aux4_Individual_effects_and_co-occurrence_patterns.R")

## *****************************************************************************
## 8) Complementary comparison between LD metrics and co-occurrence ones at allele level ----
## _____________________________________________________________________________

source("R_scripts/Aux5_LD_vs_Co_occurence_metrics.R")

## *****************************************************************************
## 9) Close the parallel backend and save the environment for further inspection ----
## _____________________________________________________________________________

stopCluster(cl)

save.image('Env') ## In case you want to inspect the results closely.