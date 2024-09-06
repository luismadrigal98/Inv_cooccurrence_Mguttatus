##' Individual effects of inversion over survival
##' 
##' First step in the survival based analysis of inversion over M. guttatus
##' plants. Next step will be co-occurrence analysis inversion pairs to explore
##' patterns of repulsion and attraction.
##' 
##' @Author: Luis Javier Madriga-Roca & John K. Kelly
##' @Date: 09/03/2024
##' ____________________________________________________________________________
##' <<<<<<<<<<<<<<<<<<<<<<<<<<<< Start >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

## *****************************************************************************
## 1) Load the required libraries ----
## _____________________________________________________________________________

library(parallel)
library(dplyr)
library(ggplot2)
library(DescTools)

cl <- makeCluster(detectCores())

set.seed(1998)

## *****************************************************************************
## 2) Setting the work directory ----
## _____________________________________________________________________________

setwd("/home/l338m483/scratch/R_Directory")

## *****************************************************************************
## 3) Importing the data ----
## _____________________________________________________________________________

Data <- list()

for(file in list.files('/home/l338m483/scratch/Cooccurrence_Inv/CSVs', 
                       full.names = T))
{
  name <- strsplit(x = file, split = '/')[[1]][7]
  name <- strsplit(x = name, split = '_')[[1]][4]
  name <- name <- sub("\\..*$", "", name)
  
  Data[[paste0("Data_", name)]] <- read.csv(file, header = T)
  
  rm(file)
  rm(name)
}

dosage_splitter <- function(Data, probes_col = 1)
{
  Probes <- Data[, probes_col]
  
  new_df <- data.frame()
  
  for (i in names(Data[, -probes_col]))
  {
    old_column <- Data[, i]
    df <- setNames(data.frame(ifelse(test = old_column == 1, 1, 0),
                              ifelse(test = old_column == 2, 1, 0)),
                   c(paste0(i, "_1"), paste0(i, "_2")))
    
    if (length(new_df) != 0) 
    {
      new_df <- cbind(new_df, df)
    }
    
    else
    {
      new_df <- df
    }
  }
  
  return(cbind(Probes, new_df))
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

metadata <- read.csv("/home/l338m483/scratch/Cooccurrence_Inv/R_directory/inv_and_gene_metadata.csv", 
                     header = T)

## *****************************************************************************
## 4) Distortion segregation analysis ----
## _____________________________________________________________________________

G_calculator <- function(row, dosage_matrix, exp = c(1/4, 1/2, 1/4)) {
  INV <- rownames(dosage_matrix)[row]
  
  # Manually count the occurrences of 0, 1, and 2
  counts <- c(sum(dosage_matrix[row, ] == 0),
              sum(dosage_matrix[row, ] == 1),
              sum(dosage_matrix[row, ] == 2))
  
  # Print the counts for debugging
  print(counts)
  
  # Ensure that counts and exp have the same length
  if (length(counts) != length(exp)) {
    warning("Observed counts and expected probabilities do not match in length.")
    return(c(Inv = INV, G = NA, df = NA, p_value = NA))
  }
  
  test_results <- GTest(counts, p = exp, correct = 'williams')
  
  return(data.frame(Inv = INV, G = test_results$statistic, 
                    df = test_results$parameter,
           p_value = test_results$p.value))
}

# Export the G_calculator function and any necessary objects to the cluster
clusterExport(cl, varlist = c("G_calculator", "Data_d_l", "GTest"))

# Perform the parallel computation
goodness_of_fit <- parLapply(cl, Data_d_l, function(observed_matrix) {
  mapply(G_calculator,
         1:nrow(observed_matrix),
         MoreArgs = list(dosage_matrix = observed_matrix),
         SIMPLIFY = F)
})

#' Condensing the results of each list into a unique dataframe (one data frame)
#' per list.

goodness_of_fit_df <- lapply(goodness_of_fit, function(x) do.call(rbind, x))

# Adding the line information into each data.frame

for (i in 1:length(goodness_of_fit_df))
{
  goodness_of_fit_df[[i]][["Line"]] <- rep(names(goodness_of_fit_df)[i], 
                                      nrow(goodness_of_fit_df[[i]]))
  
  # Multiple testing correction
  
  goodness_of_fit_df[[i]][["p_corrected"]] <- 
    p.adjust(goodness_of_fit_df[[i]]$p_value, 
            method = "BH")
  rm(i)
}

# Merging all the data.frames into a unique one

goodness_of_fit_df <- do.call(rbind, goodness_of_fit_df)

## *****************************************************************************
## 5) Pattern of selection analysis ----
## _____________________________________________________________________________

##' For this section we are going to perform the sequential X^2 analysis 
##' proposed by Fu et al. (2020). For that, we are going to use the information
##' related to those inversions that showed a significant deviation from the
##' Mendelian expectations.

# ~ Extracting the inversions that showed a significant deviation ----

significant_inv <- goodness_of_fit_df[goodness_of_fit_df$p_corrected < 0.05, ]

observed_count_extractor <- function(dosage_list, line, Inv_ID)
{
  #' This function will take as an input a list of dosage matrices (each matrix)
  #' refers to a particular line. THen, using the line and Inv_ID, the dosage 
  #' level per plant is going to be extracted.
  #' 
  #' @Arguments: dosage_list: list of dosage matrices
  #'             line: line of the plant (sublist)
  #'             Inv_ID: Inversion ID (row ID in the dosage matrix)
  #' @Retunrs: A vector with the dosage levels of the plants
  
  obs <- dosage_list[[line]][Inv_ID, ]
  
  RR <- sum(obs == 0) # Homozygous for the reference allele
  RA <- sum(obs == 1) # Heterozygous
  AA <- sum(obs == 2) # Homozygous for the alternative allele (inversion)
  
  return(data.frame(Inv_ID = Inv_ID, Line = line,
                    RR = RR, RA = RA, AA = AA))
}

# ~ Extract the counts per category ----

# Export the G_calculator function and any necessary objects to the cluster
clusterExport(cl, varlist = c("observed_count_extractor", "Data_d_l", 
                              "significant_inv"))

observed_counts_deviants <- parLapply(cl, significant_inv$Line, 
                                      function(line) {
  mapply(observed_count_extractor,
         significant_inv$Inv,
         MoreArgs = list(dosage_list = Data_d_l),
         line = line,
         SIMPLIFY = F)
})

# Stop the cluster after use
stopCluster(cl)

##' <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<< End >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>