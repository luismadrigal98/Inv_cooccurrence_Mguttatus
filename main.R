##' @title Are you with me? Co-occurrence tests from community ecology can 
##' identify positive and negative epistasis between inversions in Mimulus 
##' guttatus
##' 
##' @description This script will source and execute the different analysis employed
##' in this study, and will generate the figures in the Results sub directory.
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
                        "gtools", "viridis", "car", "nortest")

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

Data <- list()

for(file in list.files('Data/CSVs', 
                       full.names = T))
{
  name <- strsplit(x = file, split = '/')[[1]][7]
  name <- strsplit(x = name, split = '_')[[1]][4]
  name <- name <- sub("\\..*$", "", name)
  
  Data[[paste0("Data_", name)]] <- read.csv(file, header = T)
  
  rm(file)
  rm(name)
}

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

## *****************************************************************************
## 4) Simulation analysis ----
## _____________________________________________________________________________

source("R_scripts/Aux1_Contingency_tables_simulation.R")

## *****************************************************************************
## 5) Segregation distortion analysis (SD) ----
## _____________________________________________________________________________

source("R_scripts/Aux2_Segregation_distortion_analysis.R")

save.image('Env') ## In case you want to inspect the results closely.