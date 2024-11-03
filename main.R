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
                        )

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
## 3) Simulation analysis ----
## _____________________________________________________________________________

source("R_scripts/Contingency_tables_simulation.R")

save.image('Env') ## In case you want to inspect the results closely.