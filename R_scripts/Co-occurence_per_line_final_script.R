## Co_occurrence of inversions in lines of Mimulus guttatus
#@ Main Script
#@ Authors: Luis J. Madrigal-Roca & John K. Kelly
#@ Date: 2024-05-27

## *****************************************************************************
## 1) Loading the required libraries ----
## _____________________________________________________________________________

library(ggplot2)
library(readxl)
library(reshape)
library(gplots)
library(parallel)
library(igraph)
library(tidygraph)
library(ggraph)
library(reshape2)
library(Cairo)
library(vegan)
library(CooccurrenceAffinity)
library(jaccard)
library(dplyr)
library(tidyverse)
library(reshape2)
library(gprofiler2)
library(rcompanion)
library(chisq.posthoc.test)
library(MASS)
library(repmod)
library(gtools)
library(viridis)
library(car)
library(nortest)

set.seed(1998)

## *****************************************************************************
## 2) Setting the work directory ----
## _____________________________________________________________________________

setwd("/home/l338m483/scratch/Cooccurrence_Inv/R_directory")

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

## Importing the name of the genes contained in each inversion

## Additional checking step for compliance in the marginals.
## If there is an inversion with a dosage level of 0, it will be removed.
## Its effect will be captured in the individual effect on survival analysis

# Searching for problematic in the data

# Function to check for zero marginals in the contingency table
check_zero_marginals <- function(data, row_index, col_index) {
  row_name <- rownames(data)[row_index]
  col_name <- rownames(data)[col_index]

  observed <- matrix(NA, 3, 3)

  for (x in 0:2) {
    for (y in 0:2) {
      observed[x + 1, y + 1] <- sum(data[row_index, ] == x & data[col_index, ] == y)
    }
  }

  # Check for zero marginals
  if (any(rowSums(observed) == 0) || any(colSums(observed) == 0)) {
    return(list(row_name = row_name, col_name = col_name, observed = observed))
  } else {
    return(NULL)
  }
}

marginal_compliance_checking <- function(data) {
  problematic_pairs <- list()
  for (i in 1:(nrow(data) - 1)) {
    for (j in (i + 1):nrow(data)) {
      result <- check_zero_marginals(data, i, j)
      if (!is.null(result)) {
        problematic_pairs <- append(problematic_pairs, list(result))
      }
    }
  }

  return(problematic_pairs)
}

problematic_pairs <- lapply(Data_d_l, marginal_compliance_checking)

#' There is a particularly interesting case in the data. Inv_49, in the line 
#' 444 does not exhibit the homozygous state for the inversion (dosage 2). It is
#' pertinent then to exclude that inversion from posterior analysis to avoid 
#' problems with marginals 0.

## *****************************************************************************
## 3) Chi-square test for independence of inversions ----
## _____________________________________________________________________________

expec_freq <- matrix(data = c(1/16, 2/16, 1/16, 2/16, 4/16, 2/16,
                              1/16, 2/16, 1/16), nrow = 3, ncol = 3, 
                     dimnames = list(c('0', '1', '2'),
                                     c('0', '1', '2')))

# 3.2) X2 square test ----

x2_calculator <- function(i, j, dosage_matrix)
{
  #' This function calculates the X2 test for independence between two inversions
  #' in the dataset. The function will return a dataframe with the results of the
  #' test.
  #' 
  #' @param i The row index of the first inversion
  #' @param j The row index of the second inversion
  #' @param dosage_matrix The matrix containing the dosage levels of the inversions
  #' 
  #' @return A dataframe with the results of the X2 test
  #' 
  #' @note In some instances, the result test will not have some of the elements.
  #' This happen because for INV_49, there is one dosage level that is not present,
  #' and this result in a contingency analysis of a 2x3 table. This is not a problem
  #' for the X2 test, but it is for the post-hoc analysis. We will preserve the
  #' structure of the output, but for those missing elements, the values will be
  #' NA's.
  #' ___________________________________________________________________________
  
  INV1_name <- rownames(dosage_matrix)[i]
  INV2_name <- rownames(dosage_matrix)[j]
  
  observed <- table(dosage_matrix[i,], dosage_matrix[j,])
  
  test <- chisq.test(x = observed, simulate.p.value = T,
                     B = 10000)
  
  x2 <- test$statistic
  p_value <- test$p.value
  dev <- test$observed - test$expected
  standardized_residuals <- dev / sqrt(test$expected)
  relative_contribution <- (dev^2 / (test$expected * x2)) * 100
  
  return(data.frame(INV_1 = INV1_name,
                    INV_2 = INV2_name,
                    X2 = x2,
                    p = p_value,
                    dev_1r_1c = dev[1, 1],
                    dev_2r_1c = dev[2, 1],
                    dev_3r_1c = if(nrow(dev) == 3) dev[3, 1] else NA,
                    dev_1r_2c = dev[1, 2],
                    dev_2r_2c = dev[2, 2],
                    dev_3r_2c = if(nrow(dev) == 3) dev[3, 2] else NA,
                    dev_1r_3c = if(ncol(dev) == 3) dev[1, 3] else NA,
                    dev_2r_3c = if(ncol(dev) == 3) dev[2, 3] else NA,
                    dev_3r_3c = if(ncol(dev) == 3 & nrow(dev) == 3) dev[3, 3] else NA,
                    sr_1r_1c = standardized_residuals[1, 1],
                    sr_2r_1c = standardized_residuals[2, 1],
                    sr_3r_1c = if(nrow(standardized_residuals) == 3) standardized_residuals[3, 1] else NA,
                    sr_1r_2c = standardized_residuals[1, 2],
                    sr_2r_2c = standardized_residuals[2, 2],
                    sr_3r_2c = if(nrow(standardized_residuals) == 3) standardized_residuals[3, 2] else NA,
                    sr_1r_3c = if(ncol(standardized_residuals) == 3) standardized_residuals[1, 3] else NA,
                    sr_2r_3c = if(ncol(standardized_residuals) == 3) standardized_residuals[2, 3] else NA,
                    sr_3r_3c = if (ncol(standardized_residuals) == 3 & 
                                   nrow(standardized_residuals) == 3) standardized_residuals[3, 3] else NA,
                    rel_cont_1r_1c = relative_contribution[1, 1],
                    rel_cont_2r_1c = relative_contribution[2, 1],
                    rel_cont_3r_1c = if(nrow(relative_contribution) == 3) relative_contribution[3, 1] else NA,
                    rel_cont_1r_2c = relative_contribution[1, 2],
                    rel_cont_2r_2c = relative_contribution[2, 2],
                    rel_cont_3r_2c = if(nrow(relative_contribution) == 3) relative_contribution[3, 2] else NA,
                    rel_cont_1r_3c = if(ncol(relative_contribution) == 3) relative_contribution[1, 3] else NA,
                    rel_cont_2r_3c = if(ncol(relative_contribution) == 3) relative_contribution[2, 3] else NA,
                    rel_cont_3r_3c = if(ncol(relative_contribution) == 3 &
                                        nrow(relative_contribution) == 3) relative_contribution[3, 3] else NA))
}

results_x2 <- lapply(X = Data_d_l, FUN = function(observed_matrix)
{
  mapply(x2_calculator,
         rep(1:nrow(observed_matrix), 
             each = nrow(observed_matrix)), 
         rep(1:nrow(observed_matrix), 
             nrow(observed_matrix)),
         MoreArgs = list(dosage_matrix = observed_matrix),
         SIMPLIFY = F)
})

# Condensing the results into a dataframe per line

results_x2_df <- lapply(results_x2, function(x)
{
  do.call(rbind, x)
})

# Exporting these results to perform analysis in main 4 script

saveRDS(results_x2_df, file = "Results/results_x2_df.rds")

# Extracting the p_values as an square matrix ---

x2p_to_square <- function(x2_result, col_names = c("INV_1", "INV_2"), 
                          p_col = 'p') 
{
  x2_result <- as.data.frame(x2_result)
  
  unique_values <- unique(c(x2_result[, col_names[1]], 
                            x2_result[, col_names[2]]))
  
  p_values <- matrix(nrow = length(unique_values), 
                     ncol = length(unique_values))
  
  rownames(p_values) <- colnames(p_values) <- unique_values
  
  for (i in unique_values) {
    for (j in unique_values) {
      p_values[i, j] <- x2_result[x2_result[, col_names[1]] == i & 
                                    x2_result[, col_names[2]] == j , p_col]
    }
  }
  
  return(p_values)
}

x2_p_square_matrix <- lapply(results_x2_df, x2p_to_square)

# Filtering the results

contingency_filter <- function(results, metainfo)
{
  results <- results |>
    filter(INV_1 != INV_2) |>
    mutate(ID1 = sapply(strsplit(INV_1, "_"), `[`, 2),
           ID2 = sapply(strsplit(INV_2, "_"), `[`, 2)) |>
    mutate(chrom1 = metainfo[match(ID1, metainfo$INV_ID), "Chr"],
           chrom2 = metainfo[match(ID2, metainfo$INV_ID), "Chr"]) |>
    filter(chrom1 != chrom2)
  
  results$INV_combination <- apply(results[, c("INV_1", "INV_2")], 1, 
                                   function(x) paste(sort(x), collapse = ""))
  
  results <- results[!duplicated(results$INV_combination), ]
  
  return(results)
}

results_x2_filtered_df <- sapply(X = results_x2_df,
                                 FUN = contingency_filter,
                                 metadata,
                                 simplify = F)

# Correction of the p-values ----

p_corrector <- function(df, var = "p", method = "BH")
{
  df <- df |> mutate(p_corrected = p.adjust(df[, var], method = method))
  
  return(df)
}

results_x2_filtered_df <- sapply(X = results_x2_filtered_df, 
                                 FUN = p_corrector, var = "p", 
                        simplify = F)

# Getting the number of combinations without dosage
combos_without_d <- unique(do.call(
  'rbind', results_x2_filtered_df)$INV_combination)

## Additional analysis ----

# Expanding the dataframe to decompose each test into individual components
# This will preserve the omnibus p-value of the X2 test for the four contrasts
# derived from the test

contingency_expanded <- function(cont_res_df) 
{
  # Define a helper function
  
  process_cols <- function(df, cols, value_name, prefix) {
    df |>
      dplyr::select(all_of(c("INV_1", "INV_2", "X2", "p", "p_corrected", "ID1", 
                      "ID2", "chrom1", 
                      "chrom2", "INV_combination", cols))) |>
      pivot_longer(cols = all_of(cols), names_to = "cell_location", 
                   values_to = value_name) |>
      mutate(cell_location = gsub(prefix, "", cell_location))
  }
  
  # Define the sets of columns
  cols_list <- list(
    dev = list(cols = c("dev_2r_2c", "dev_2r_3c", "dev_3r_2c", "dev_3r_3c"), 
               value_name = "deviation_value", prefix = "dev_"),
    sr = list(cols = c("sr_2r_2c", "sr_2r_3c", "sr_3r_2c", "sr_3r_3c"), 
              value_name = "standardized_residual", prefix = "sr_"),
    rel_cont = list(cols = c("rel_cont_2r_2c", "rel_cont_2r_3c", 
                             "rel_cont_3r_2c", "rel_cont_3r_3c"), 
                    value_name = "relative_contribution", prefix = "rel_cont_")
  )
  
  # Apply the helper function to each set of columns
  df_list <- lapply(cols_list, function(x) process_cols(cont_res_df, x$cols, 
                                                        x$value_name, x$prefix))
  
  # Merge the data frames together
  df <- Reduce(function(x, y) merge(x, y, by = c("INV_1", "INV_2", "X2", "p",
                                                 "p_corrected",
                                                 "ID1", "ID2", "chrom1", 
                                                 "chrom2", "INV_combination", 
                                                 "cell_location")), df_list)
  
  return(df)
}

x2_results_filtered_expanded <- lapply(results_x2_filtered_df, 
                                       contingency_expanded)

for(i in seq_along(x2_results_filtered_expanded))
{
  x2_results_filtered_expanded[[i]] <- x2_results_filtered_expanded[[i]] |>
    mutate(cell_location = as.factor(cell_location))
  
  levels(x2_results_filtered_expanded[[i]]$cell_location) <- 
    c("_1-_1", "_1-_2", "_2-_1", "_2-_2")
  
  x2_results_filtered_expanded[[i]] <- x2_results_filtered_expanded[[i]] |>
    mutate(cell_location = as.character(cell_location)) |>
    mutate(INV_1 = paste0(INV_1, sapply(strsplit(cell_location, '-'), `[`, 1)),
           INV_2 = paste0(INV_2, sapply(strsplit(cell_location, '-'), `[`, 2)))
}

for (i in 1:length(names(x2_results_filtered_expanded))) {
  x2_results_filtered_expanded[[i]] <- x2_results_filtered_expanded[[i]] |>
    rowwise() |>
    mutate(INV_combination = paste0(sort(c(INV_1, INV_2)), collapse = "")) |>
    ungroup()
}

for (i in 1:length(names(x2_results_filtered_expanded))) {
  x2_results_filtered_expanded[[i]] <- x2_results_filtered_expanded[[i]] |>
    dplyr::rename(X2_global = X2, p_X2_global = p, p_X2_corrected = p_corrected)
  
  rm(i)
}

# 3.1) Post-hoc analysis for the X2 test ----

x2_posthoc_calculator <- function(i, j, dosage_matrix, expectation)
{
  result <- list()
  
  INV1_name <- rownames(dosage_matrix)[i]
  INV2_name <- rownames(dosage_matrix)[j]
  
  observed <- matrix(NA, 3, 3)
  
  for (x in 0:2) {
    for (y in 0:2) {
      observed[x + 1, y + 1] <- sum(dosage_matrix[i,] == x & 
                                      dosage_matrix[j,] == y)
    }
  }
  
  rownames(observed) <- c(paste0(INV1_name, "_0"), paste0(INV1_name, "_1"), 
                          paste0(INV1_name, "_2"))
  colnames(observed) <- c(paste0(INV2_name, "_0"), paste0(INV2_name, "_1"), 
                          paste0(INV2_name, "_2"))
  
  res <- chisq.posthoc.test(observed, method = 'none', p = expectation, 
                            simulate.p.value = T, 
                            B = 10000)
  
  df_melt <- melt(res, id.vars = c("Dimension", "Value"))
  
  # Split the data into two data frames
  df_residuals <- df_melt |> filter(Value == "Residuals") |> 
    dplyr::rename(Residual = value) |> dplyr::select(-Value)
  df_pvalues <- df_melt |> filter(Value == "p values") |> 
    dplyr::rename(p_value = value) |> dplyr::select(-Value)
  
  # Join the data frames
  df_final <- full_join(df_residuals, df_pvalues, 
                        by = c("Dimension", "variable"))
  
  # Rename the columns
  colnames(df_final) <- c("INV1", "INV2", "Residual", "p_value")
  
  return(df_final)
}

results_posthoc_x2 <- lapply(X = Data_d_l, FUN = function(observed_matrix)
{
  mapply(x2_posthoc_calculator,
         rep(1:nrow(observed_matrix), 
             each = nrow(observed_matrix)), 
         rep(1:nrow(observed_matrix), 
             nrow(observed_matrix)),
         MoreArgs = list(dosage_matrix = observed_matrix,
                         expectation = expec_freq),
         SIMPLIFY = F)
})

results_posthoc_x2_df <- lapply(results_posthoc_x2, function(x)
{
  do.call(rbind, x)
})

for (i in names(results_posthoc_x2_df))
{
  results_posthoc_x2_df[[i]] <- results_posthoc_x2_df[[i]] |>
    mutate(INV2 = as.character(INV2)) |>
    rowwise() |>
    mutate(INV_combination = paste0(sort(c(INV1, INV2)), collapse = '')) |>
    filter(! any(grepl("_0", INV1) | grepl("_0", INV2))) |>
    ungroup()
  
  rm(i)
}

## Extracting the p_values as an square matrix ----

x2_p_post_hoc_square_matrix <- lapply(results_posthoc_x2_df, x2p_to_square,
                                      c("INV1", "INV2"), 'p_value')

## Filtering out irrelevant contrasts ----

results_posthoc_x2_df_filtered <- sapply(X = results_posthoc_x2_df, 
                                         FUN = function(x)
{
  x |>
    mutate(ID1 = sapply(strsplit(INV1, "_"), `[`, 2),
           ID2 = sapply(strsplit(INV2, "_"), `[`, 2)) |>
    mutate(chrom1 = metadata[match(ID1, metadata$INV_ID), "Chr"],
           chrom2 = metadata[match(ID2, metadata$INV_ID), "Chr"]) |>
    filter(chrom1 != chrom2)
}, simplify = F)

for (i in names(results_posthoc_x2_df_filtered))
{
  results_posthoc_x2_df_filtered[[i]] <- 
    results_posthoc_x2_df_filtered[[i]][!duplicated(
      results_posthoc_x2_df_filtered[[i]]$INV_combination), ]
  
  rm(i)
}

for (i in names(results_posthoc_x2_df_filtered))
{
  results_posthoc_x2_df_filtered[[i]] <- results_posthoc_x2_df_filtered[[i]] |>
    mutate(p_adjusted = p.adjust(p_value, method = "BY"))
  
  rm(i)
}

## Condensing all the results (main and post hoc) into a single dataframe

for (i in 1: length(names(x2_results_filtered_expanded)))
{
  x2_results_filtered_expanded[[i]] <- merge(x2_results_filtered_expanded[[i]], 
                                             results_posthoc_x2_df_filtered[[i]], 
                                             by = "INV_combination", 
                                             suffixes = c("", "_posthoc"))
  
  rm(i)
}

for (i in 1: length(names(x2_results_filtered_expanded)))
{
  x2_results_filtered_expanded[[i]] <- x2_results_filtered_expanded[[i]] |>
    dplyr::select(!c("ID1", "ID2", "cell_location", 
              "standardized_residual", "INV1", "INV2", "ID1_posthoc",
              "ID2_posthoc", "chrom1_posthoc", "chrom2_posthoc")) |>
    dplyr::rename(p_X2_posthoc = p_value,
           p_X2_posthoc_corrected = p_adjusted)
}

## Getting the number of combinations with dosage
combos_with_d <- unique(do.call(
  'rbind', x2_results_filtered_expanded)$INV_combination)

## *****************************************************************************
## 4) Calculation of alternative indexes used in community studies ----
## _____________________________________________________________________________

# Filtering out rows with 0 marginals. If something is not there in any instance,
# we cannot test for association with another inversion.

Data_p_a <- lapply(Data_p_a, function(x) x[!apply(x, 1, 
                                                  function(y) all(y == 0)), ])

# Instance INV_49_2 in line 444 was eliminated from the analysis

# 4.1) Jaccard ----

jaccard_tanimoto_results <- sapply(X = Data_p_a, 
                                   FUN = jaccard.test.pairwise,
                                   method = 'bootstrap',
                                   compute.qvalue = F,
                                   simplify = F,
                                   B = 10000)

## Correcting the dim names of the outputs

for (i in 1:length(Data_p_a))
{
  denominations <- attributes(Data_p_a[[i]])[[2]][[1]]
  
  for (j in names(jaccard_tanimoto_results[[i]]))
  {
    dimnames(jaccard_tanimoto_results[[i]][[j]]) <- list(denominations,
                                                         denominations)
  }
  
  rm(i)
  rm(j)
}

# Extracting the p_values as a square matrix ----

jaccard_tanimoto_p_square_matrix <- list()

for (i in 1:length(jaccard_tanimoto_results)) {
  # Ensure that pvalues is a matrix
  mat <- as.matrix(jaccard_tanimoto_results[[i]]$pvalues)
  
  # Make the matrix symmetric
  mat[lower.tri(mat)] <- t(mat)[lower.tri(mat)]
  
  # Assign the symmetric matrix back to the list
  jaccard_tanimoto_p_square_matrix[[i]] <- mat
}

## Extracting the results as a df ---

jaccard_tanimoto_condenser <- function(results, metainfo)
{
  df <- data.frame(entity_1 = NA,
                   Dosage_1 = NA,
                   ID_1 = NA,
                   entity_2 = NA,
                   Dosage_2 = NA,
                   ID_2 = NA, 
                   chrom_1 = NA, 
                   chrom_2 = NA, 
                   jaccard = NA, 
                   expectation = NA, 
                   p_value = NA)
  
  ref <- which(! is.na(results[[2]]), arr.ind = T)
  
  for (i in 1:nrow(ref))
  {
    row_i <- ref[i, 1]
    col_i <- ref[i, 2]
    
    INV_1 <- rownames(results[[1]])[row_i]
    Dosage_1 <- strsplit(INV_1, "_")[[1]][3]
    ID_1 <- strsplit(INV_1, "_")[[1]][2]
    
    INV_2 <- colnames(results[[1]])[col_i]
    Dosage_2 <- strsplit(INV_2, "_")[[1]][3]
    ID_2 <- strsplit(INV_2, "_")[[1]][2]
    
    new_obs <- data.frame(entity_1 = INV_1,
                          Dosage_1 = Dosage_1,
                          ID_1 = ID_1,
                          entity_2 = INV_2,
                          Dosage_2 = Dosage_2,
                          ID_2 = ID_2,
                          chrom_1 = metainfo[metainfo$INV_ID == ID_1, 'Chr'], 
                          chrom_2 = metainfo[metainfo$INV_ID == ID_2, 'Chr'], 
                          jaccard = results[['statistics']][row_i, col_i], 
                          expectation = results[['expectation']][row_i, col_i], 
                          p_value = results[['pvalues']][row_i, col_i])
    
    df <- rbind2(df, new_obs)
  }
  
  df = df[-1, ]
  
  return(df)
}

jaccard_tanimoto_results_df <- sapply(X = jaccard_tanimoto_results,
                                      FUN = jaccard_tanimoto_condenser,
                                      metadata,
                                      simplify = F)

# Filtering out elements in the same chromosome

jaccard_tanimoto_results_df_filtered <- jaccard_tanimoto_results_df

for (i in 1:length(jaccard_tanimoto_results_df_filtered))
{
  jaccard_tanimoto_results_df_filtered[[i]] <- 
    jaccard_tanimoto_results_df_filtered[[i]] |>
    dplyr::filter(chrom_1 != chrom_2)
  
  rm(i)
}

# Correcting the p_values ----
jaccard_tanimoto_results_df_filtered <- sapply(X = jaccard_tanimoto_results_df_filtered, 
                                      FUN = p_corrector, var = "p_value", 
                                      method = "BY",
                                      simplify = F)

# 4.2) Affinity measure ----

results_affinity <- lapply(X = Data_p_a, FUN = function(x)
{
  x <- affinity(x, row.or.col = 'row', 
  datatype = 'binary', sigdigit = 3)
  x
})

for (i in 1:length(results_affinity))
{
  results_affinity[[i]]$all <- results_affinity[[i]]$all |>
    mutate(p_value = as.numeric(p_value))
  
  rm(i)
}

# Extracting the affinity values from the output as a matrix for plotting

alpha_values_extractor <- function(affinity_results) {
  
  df <- affinity_results$all
  inv_ID <- names(affinity_results$occur_mat)
  
  A_matrix <- matrix(nrow = length(inv_ID), ncol = length(inv_ID))
  
  dimnames(A_matrix) <- list(inv_ID, inv_ID)
  
  for (i in inv_ID) {
    for(j in inv_ID) {
      if(i == j) {
        A_matrix[i, j] <- 1 # Maximum value (given the capped behavior of affinity calculator)
      } else if (substr(i, 1, nchar(i) - 2) == substr(j, 1, nchar(j) - 2)) {
        A_matrix[i, j] <- 1 # A plant cannot be both homozygous and heterozygous at the same time
      } else {
        alpha_p <- df[df[, "entity_1"] == i & df[, "entity_2"] == j, "p_value"]
        if (length(alpha_p) == 0) {
          alpha_p <- df[df[, "entity_1"] == j & df[, "entity_2"] == i, "p_value"]
        }
        A_matrix[i, j] <- ifelse(length(alpha_p) == 0, NA, alpha_p)
      }
    }
  }
  
  return(A_matrix)
}

A_matrices <- lapply(results_affinity, alpha_values_extractor)

for (i in 1:length(names(A_matrices))) {
  old_name <- names(A_matrices)[i]
  new_name <- paste0("A_", old_name)
  names(A_matrices)[i] <- new_name
  rm(i)
  rm(new_name)
  rm(old_name)
}

affinity_filter <- function(affinity_res, metadata, filter_LG = T, 
                            filter_CS = T)
{
  # filter_LG: Removes entries that belongs to the same chromosome
  # filter_CS: Removes entries from the same inversion with mutually exclusive
  #states
  
  df <- affinity_res$all
  
  df <- df |> separate(col = entity_1, into = c("Generic_name_1", "INV_ID_1", 
                                                "State_1"), sep = '_', 
                       remove = F) |>
    separate(col = entity_2, into = c("Generic_name_2", "INV_ID_2", 
                                      "State_2"), sep = '_', 
             remove = F) |>
    dplyr::select(-Generic_name_1, -State_1, -Generic_name_2, -State_2)
  
  df <- df |> 
    left_join(metadata |> dplyr::select(INV_ID, Chr), by = c("INV_ID_1" = "INV_ID")) |> 
    rename(chrom1 = Chr)
  
  df <- df |> 
    left_join(metadata |> dplyr::select(INV_ID, Chr), by = c("INV_ID_2" = "INV_ID")) |> 
    rename(chrom2 = Chr)
  
  if (filter_LG == T)
  {
    df <- df |> filter(chrom1 != chrom2)
  }
  
  if (filter_CS == T)
  {
    df <- df |> filter(INV_ID_1 != INV_ID_2)
  }
  
  return(df)
}

results_affinity_filtered <- sapply(X = results_affinity,
                                    FUN = affinity_filter,
                                    metadata,
                                    simplify = F)

# Correcting the p_values ----

results_affinity_filtered <- sapply(X = results_affinity_filtered, 
                                    FUN = p_corrector, 
                                    var = "p_value", 
                                    method = "BY",
                                    simplify = F)

# 5) Extracting the significant values detected for the three scores ----

# Common field for all datasets: INV_combination

for (i in 1:length(names(results_affinity_filtered))) {
  results_affinity_filtered[[i]] <- results_affinity_filtered[[i]] |>
    rowwise() |>
    mutate(INV_combination = paste0(sort(c(entity_1, entity_2)), collapse = "")) |>
    ungroup()
}

for (i in 1:length(names(jaccard_tanimoto_results_df_filtered))) {
  jaccard_tanimoto_results_df_filtered[[i]] <- 
    jaccard_tanimoto_results_df_filtered[[i]] |>
    rowwise() |>
    mutate(INV_combination = paste0(sort(c(entity_1, entity_2)), 
                                    collapse = "")) |>
    ungroup()
  
  rm(i)
}

## Merging all dataframes into a single one per line----

merged_results <- list()

for (i in 1:9) {
  merged_results[[i]] <- merge(x2_results_filtered_expanded[[i]], 
                               results_affinity_filtered[[i]], 
                               by = "INV_combination", 
                               suffixes = c("_x2", "_affinity"))
  
  merged_results[[i]] <- merge(merged_results[[i]], 
                               jaccard_tanimoto_results_df_filtered[[i]], 
                               by = "INV_combination", 
                               suffixes = c("", "_jaccard_tanimoto"))
}

names(merged_results) <- names(x2_results_filtered_expanded)

## Combining all data frames in an unique one ----

# Add a new column to each data frame with the line information
for (i in 1:length(merged_results)) {
  merged_results[[i]]$line <- names(merged_results)[i]
}

# Combine all data frames into one
final_result <- do.call(rbind, merged_results)

# Depuration of the final data frame

final_result <- final_result |>
  dplyr::select(line, INV_1, INV_2, chrom_1, chrom_2, X2_global, p_X2_global, 
         p_X2_posthoc, p_X2_posthoc_corrected, relative_contribution, 
         Residual, alpha_mle, p_value, 
         p_corrected, jaccard_jaccard_tanimoto, p_value_jaccard_tanimoto,
         p_corrected_jaccard_tanimoto) |>
  rename(Line = line,
         Chr_1 = chrom_1, Chr_2 = chrom_2,
         X2_SR = Residual, X2_RC = relative_contribution,
         A_p = p_value, A_alpha = alpha_mle, 
         A_p_corrected = p_corrected,
         Jaccard = jaccard_jaccard_tanimoto, J_p = p_value_jaccard_tanimoto, 
         J_p_corrected = p_corrected_jaccard_tanimoto)

# 5.1) Relaxed results ----

final_result_relaxed_any <- final_result |> 
  filter(p_X2_posthoc < 0.05 | A_p < 0.05 | J_p < 0.05)

final_result_relaxed_all <- final_result |>
  filter(p_X2_global < 0.05 & p_X2_posthoc < 0.05 & A_p < 0.05 & J_p < 0.05)

# 5.2) Stringent results ----

final_result_stringent_any <- final_result |> 
  filter(p_X2_posthoc_corrected < 0.05 | A_p_corrected < 0.05 | 
           J_p_corrected < 0.05)

final_result_stringent_all <- final_result |>
  filter(p_X2_global < 0.05 & p_X2_posthoc_corrected < 0.05 & 
           A_p_corrected < 0.05 & J_p_corrected < 0.05)

## Saving the data frame sinto csv files

write.csv(final_result_relaxed_any, "Results/final_result_relaxed_any.csv", 
          row.names = FALSE)
write.csv(final_result_relaxed_all, "Results/final_result_relaxed_all.csv", 
          row.names = FALSE)
write.csv(final_result_stringent_any, "Results/final_result_stringent_any.csv", 
          row.names = FALSE)
write.csv(final_result_stringent_all, "Results/final_result_stringent_all.csv", 
          row.names = FALSE)

# 6) Visualization of the relevant patterns in terms of evidence ----

# 6.1) Evidence according to X2_global ----

for (i in names(x2_p_square_matrix)) 
{
  pdf(file = paste0("Plots/", i, "_omnibus_X2_support_heatmap.pdf"), width = 6, 
      height = 6)
  
  logical_matrix <- x2_p_square_matrix[[i]] < 0.05
  info <- matrix(as.numeric(logical_matrix), nrow = nrow(logical_matrix), 
                 ncol = ncol(logical_matrix))
  rownames(info) <- colnames(info) <- rownames(logical_matrix)
  
  diag(info) <- NA
  
  heatmap.2(x = info, dendrogram = 'none', breaks = 3, 
            symbreaks = F,
            col = c("white", "#64AC59"),
            key = T, density.info = 'histogram', key.title = "Support",
            main = paste0("Pairwise", "_", i), xlab = "Inversions",
            ylab = "Inversions", trace = "none", denscol = "black",
            densadj = 0.50)
  
  dev.off()
  
  rm(i)
  rm(info)
}

# 6.2) Evidence according to X2_posthoc, Jaccard, and Affinity ----

support_counter <- function(supp_X2_posthoc,
                            supp_Affinity,
                            supp_Jaccard)
{
  # Remove rows and columns with NA values
  row_na <- apply(supp_X2_posthoc, 1, function(x) any(!is.na(x)))
  col_na <- apply(supp_X2_posthoc, 2, function(x) any(!is.na(x)))
  supp_X2_posthoc <- supp_X2_posthoc[row_na, col_na]
  
  names <- list(rownames(supp_X2_posthoc), colnames(supp_X2_posthoc))
  
  supp_X2_posthoc <- matrix(as.numeric(supp_X2_posthoc < 0.05), 
                            nrow = nrow(supp_X2_posthoc), 
                            ncol = ncol(supp_X2_posthoc))
  supp_Affinity <- matrix(as.numeric(supp_Affinity < 0.05), 
                          nrow = nrow(supp_Affinity), 
                          ncol = ncol(supp_Affinity))
  supp_Jaccard <- matrix(as.numeric(supp_Jaccard < 0.05), 
                         nrow = nrow(supp_Jaccard), 
                         ncol = ncol(supp_Jaccard))
  
  concensus <- supp_X2_posthoc + supp_Affinity + supp_Jaccard
  dimnames(concensus) <- names
  
  return(concensus)
}

support_matrix <- mapply(support_counter, 
                         supp_X2_posthoc = x2_p_post_hoc_square_matrix,
                         supp_Affinity = A_matrices,
                         supp_Jaccard = jaccard_tanimoto_p_square_matrix,
                         SIMPLIFY = F)

for (i in names(support_matrix)) 
{
  pdf(file = paste0("Plots/", i, "_concensus_support_heatmap.pdf"), width = 7, 
      height = 7)
  
  count_matrix <- support_matrix[[i]]
  
  heatmap.2(x = count_matrix, dendrogram = 'none', breaks = 5, 
            symbreaks = F,
            col = c("white", "#EB4F48", "#F4D166", "#64AC59"),
            key = T, density.info = 'histogram', key.title = "Support",
            main = paste0("Pairwise", "_", i), xlab = "Inversions",
            ylab = "Inversions", trace = "none", denscol = "black",
            densadj = 0.50)
  
  dev.off()
  
  rm(i)
}

## 6.3) Detecting Linkage Groups ----

# Helper function for making NA the p_values associated to inversion pairs in
# different chromosomes (masking the p-values)

mask_chromosomes <- function(p_matrix, metadata)
{
  row_names <-  rownames(p_matrix)
  row_names <- as.character(sapply(strsplit(row_names, "_"), `[`, 2))
  
  col_names <- colnames(p_matrix)
  col_names <- as.character(sapply(strsplit(col_names, "_"), `[`, 2))
  
  for (i in 1:length(row_names))
  {
    for (j in 1:length(col_names))
    {
      if (metadata[metadata$INV_ID == row_names[i], "Chr"] != 
          metadata[metadata$INV_ID == col_names[j], "Chr"])
      {
        p_matrix[i, j] <- NA
      }
    }
  }
  
  return(p_matrix)
}

# Masking the p_values of contrast in different chromosomes

x2_p_square_matrix_masked <- lapply(x2_p_square_matrix, 
                                    mask_chromosomes, metadata)

jaccard_tanimoto_p_square_matrix_masked <- lapply(jaccard_tanimoto_p_square_matrix, 
                                                 mask_chromosomes, metadata)

A_matrices_masked <- lapply(A_matrices, mask_chromosomes, metadata)

x2_p_post_hoc_square_matrix_masked <- lapply(x2_p_post_hoc_square_matrix, 
                                      mask_chromosomes, metadata)

# 6.3.1) Heatmaps for the p_values of the different metrics ----

# The inversion in the heatmap will appear sorted by chromosomes

# 6.3.1.1) Linkage groups according to omnibus X2 test ----

for (i in names(x2_p_square_matrix_masked)) 
{
  pdf(file = paste0("Plots/", i, "_omnibus_X2_linkage_heatmap.pdf"), width = 6, 
      height = 6)
  
  logical_matrix <- x2_p_square_matrix_masked[[i]] < 0.05
  info <- matrix(as.numeric(logical_matrix), nrow = nrow(logical_matrix), 
                 ncol = ncol(logical_matrix))
  
  rownames(info) <- colnames(info) <- rownames(logical_matrix)
  
  sorted_rownames <- mixedsort(rownames(info))
  
  info <- info[sorted_rownames, sorted_rownames]
  
  diag(info) <- NA
  
  info[is.na(info)] <- F
  
  heatmap.2(x = info, dendrogram = 'none', breaks = 3, 
            Rowv = F,
            Colv = F,
            symbreaks = F,
            col = c("white", "blue"),
            key = F, density.info = 'none', key.title = "",
            main = paste0("Pairwise", "_", i), xlab = "Inversions",
            ylab = "Inversions", trace = "none",
            )
  
  dev.off()
  
  rm(i)
  rm(info)
}

# 6.3.1.1) Linkage groups according to the 2x2 based tests ----

support_matrix_LG <- mapply(support_counter, 
                         supp_X2_posthoc = x2_p_post_hoc_square_matrix_masked,
                         supp_Affinity = A_matrices_masked,
                         supp_Jaccard = jaccard_tanimoto_p_square_matrix_masked,
                         SIMPLIFY = F)

# Impossible contrast cleaner (same inversion in different dosage) ----

for (i in names(support_matrix_LG)) {
  for (j in 1:nrow(support_matrix_LG[[i]])) {
    for (k in 1:ncol(support_matrix_LG[[i]])) {
      # Splitting row names and comparing the specific part
      if (strsplit(rownames(support_matrix_LG[[i]])[j], split = "_")[[1]][2] ==
          strsplit(rownames(support_matrix_LG[[i]])[k], split = "_")[[1]][2]) {
        support_matrix_LG[[i]][j, k] <- NA
      }
    }
  }
  
  rm(i)
  rm(j)
  rm(k)
}

for (i in names(support_matrix_LG)) 
{
  pdf(file = paste0("Plots/", i, "_concensus_linkage_heatmap.pdf"), width = 7, 
      height = 7)
  
  count_matrix <- support_matrix_LG[[i]]
  
  sorted_rownames <- mixedsort(rownames(count_matrix))
  
  count_matrix <- count_matrix[sorted_rownames, sorted_rownames]
  
  heatmap.2(x = count_matrix, dendrogram = 'none', breaks = 5, 
            symbreaks = F,
            Rowv = F,
            Colv = F,
            col = c("white", "#EB4F48", "#F4D166", "#64AC59"),
            key = T, density.info = 'histogram', key.title = "Support",
            main = paste0("Pairwise", "_", i), xlab = "Inversions",
            ylab = "Inversions", trace = "none", denscol = "black",
            densadj = 0.50)
  
  dev.off()
  
  rm(i)
}

# 7) Relation between statistics employed in terms of linear regression ----
# This will be done taking the overall metric across all crosses
# Linear models:

# 7.1) X2_SR vs Jaccard ----

model1 <- lm(X2_SR ~ Jaccard, data = final_result)

# 7.2) X2_SR vs Affinity ----

model2 <- lm(X2_SR ~ A_alpha, data = final_result)

# 7.3) Jaccard vs Affinity ----

model3 <- lm(Jaccard ~ A_alpha, data = final_result)

# 7.5) Plotting the results ----

# 7.5.1) Model 1 ----

pdf(file = "Plots/Linear_models_X2_SR_Jaccard.pdf", width = 8, height = 6)

ggplot(data = final_result, aes(x = Jaccard, y = X2_SR)) +
  geom_point() +
  geom_smooth(method = "lm", se = FALSE, color = "red") +
  theme_bw() +
  annotate(geom = 'text', x = -0.0475, y = 3.8, 
           label = "X2_SR ~ 27.53 * Jaccard - 0.03, df = 4122, p < 2e-16, R^2 = 0.98") +
  theme(axis.text = element_text(size = 12, color = 'black'),
        axis.title = element_text(size = 14, color = 'black'),
        legend.text = element_text(size = 12, color = 'black'),
        legend.title = element_text(size = 14, color = 'black'),
        panel.grid = element_blank())

dev.off()

# 7.5.2) Model 2 ----

pdf(file = "Plots/Linear_models_X2_SR_Affinity.pdf", width = 8, height = 6)

ggplot(data = final_result, aes(x = A_alpha, y = X2_SR)) +
  geom_point() +
  geom_smooth(method = "lm", se = FALSE, color = "red") +
  theme_bw() +
  annotate(geom = 'text', x = -4.25, y = 7, 
           label = "X2_SR ~ 0.79 * Alpha + 0.03, df = 4122, p < 2e-16, R^2 = 0.35") +
  theme(axis.text = element_text(size = 12, color = 'black'),
        axis.title = element_text(size = 14, color = 'black'),
        legend.text = element_text(size = 12, color = 'black'),
        legend.title = element_text(size = 14, color = 'black'),
        panel.grid = element_blank())

dev.off()

# 7.5.3) Model 3 ----

pdf(file = "Plots/Linear_models_Jaccard_Affinity.pdf", width = 8, height = 6)

ggplot(data = final_result, aes(x = A_alpha, y = Jaccard)) +
  geom_point() +
  geom_smooth(method = "lm", se = FALSE, color = "red") +
  theme_bw() +
  annotate(geom = 'text', x = 0, y = 1, 
           label = "Jaccard ~ 0.03 * Alpha + 0.002, df = 4122, p < 2e-16, R^2 = 0.32") +
  theme(axis.text = element_text(size = 12, color = 'black'),
        axis.title = element_text(size = 14, color = 'black'),
        legend.text = element_text(size = 12, color = 'black'),
        legend.title = element_text(size = 14, color = 'black'),
        panel.grid = element_blank()) +
  ylim(c(-1, 1))

dev.off()

## After checking the premise of normality of residuals, and given the detected
# deviations, a robust linear model will be employed.

lillie.test(model1$residuals)
lillie.test(model2$residuals)
lillie.test(model3$residuals)

# 7.6) Robust linear models ----

# 7.6.1) X2_SR vs Jaccard ----
# Robust model
robust_model1 <- rlm(X2_SR ~ Jaccard, data = final_result)

pdf(file = "Plots/Robust_Linear_models_X2_SR_Jaccard.pdf", width = 8, 
    height = 6)

ggplot(data = final_result[!robust_model1$w < 0.15,], 
       aes(x = Jaccard, y = X2_SR)) +
  geom_point() +
  geom_smooth(method = "rlm", se = FALSE, color = "red") +
  theme_bw() +
  annotate(geom = 'text', x = -0.0475, y = 3.8, 
           label = "X2_SR ~ 28.09 * Jaccard - 0.02, df = 4122, p < 0.001") +
  theme(axis.text = element_text(size = 12, color = 'black'),
        axis.title = element_text(size = 14, color = 'black'),
        legend.text = element_text(size = 12, color = 'black'),
        legend.title = element_text(size = 14, color = 'black'),
        panel.grid = element_blank())

dev.off()

# 7.6.2) X2_SR vs Affinity ----
# Robust model
robust_model2 <- rlm(X2_SR ~ A_alpha, data = final_result)

pdf(file = "Plots/Robust_Linear_models_X2_SR_Affinity.pdf", width = 8, 
    height = 6)

ggplot(data = final_result[!robust_model2$w < 0.15,], 
       aes(x = A_alpha, y = X2_SR)) +
  geom_point() +
  geom_smooth(method = "rlm", se = FALSE, color = "red") +
  theme_bw() +
  annotate(geom = 'text', x = -0.0475, y = 3.8, 
           label = "X2_SR ~ 2.45 * Affinity + 0.01, df = 4122, p < 0.01") +
  theme(axis.text = element_text(size = 12, color = 'black'),
        axis.title = element_text(size = 14, color = 'black'),
        legend.text = element_text(size = 12, color = 'black'),
        legend.title = element_text(size = 14, color = 'black'),
        panel.grid = element_blank())

dev.off()

# 7.6.3) Jaccard vs Affinity ----
# Robust model
robust_model3 <- rlm(Jaccard ~ A_alpha, data = final_result)

pdf(file = "Plots/Robust_Linear_models_Jaccard_Affinity.pdf", width = 8, 
    height = 6)

ggplot(data = final_result[!robust_model3$w < 0.15,], 
       aes(x = A_alpha, y = Jaccard)) +
  geom_point() +
  geom_smooth(method = "rlm", se = FALSE, color = "red") +
  theme_bw() +
  annotate(geom = 'text', x = -0.0475, y = 0.8, 
           label = "Jaccard ~ 0.09 * Affinity + 0.001, df = 4122, p < 0.001") +
  theme(axis.text = element_text(size = 12, color = 'black'),
        axis.title = element_text(size = 14, color = 'black'),
        legend.text = element_text(size = 12, color = 'black'),
        legend.title = element_text(size = 14, color = 'black'),
        panel.grid = element_blank()) +
  ylim(c(-1, 1))

dev.off()

# p-values for the robust models

rob.pvals(robust_model1)
rob.pvals(robust_model2)
rob.pvals(robust_model3)

# 8) Network analysis ----

# Splitting the final result data frame into a list of data frames per line

final_result_split <- split(final_result, final_result$Line)

# 8.1) Network builder function ----

nodes <- paste0("Inv_", metadata$INV_ID)

nodes <- sapply(nodes, function(x) {
  c(paste0(x, "_1"), paste0(x, "_2"))
})

nodes <- as.vector(unlist(nodes))

sig_network_builder <- function(nodes, meta_nodes, edges, type = "both",
                                co_occ_color = "blue",
                                rep_color = "red",
                                non_sig_color = "grey")
{
  # Create the nodes
  nodes <- data.frame(name = nodes)
  nodes$chromosome <- meta_nodes[, 'Chr'][match(
    sapply(strsplit(nodes$name, "_"), `[`, 2), meta_nodes$INV_ID)]
  
  nodes$name <- sapply(nodes$name, function(x) {
    paste0(strsplit(x, "_")[[1]][2], "_", strsplit(x, "_")[[1]][3])
  })
  
  # Initialize edges_net as NULL
  edges_net <- NULL
  
  # Change the names of the inversions
  
  edges$INV_1 <- sapply(edges$INV_1, function(x) {
    paste0(strsplit(x, "_")[[1]][2], "_", strsplit(x, "_")[[1]][3])
  })
  
  edges$INV_2 <- sapply(edges$INV_2, function(x) {
    paste0(strsplit(x, "_")[[1]][2], "_", strsplit(x, "_")[[1]][3])
  })
  
  # Create the edges
  edges <- edges %>% 
    rowwise() %>% 
    mutate(color = ifelse(p_X2_global > 0.05 & J_p > 0.05, 
                          non_sig_color, 
                          ifelse(all(Jaccard > 0, J_p < 0.05), 
                                 co_occ_color, 
                                 ifelse(all(Jaccard < 0, J_p < 0.05), 
                                        rep_color,
                                        ifelse(all(Jaccard > 0, J_p > 0.05), 
                                               "yellow", 
                                               "#5E4FA2")))))
  
  edges_net <- data.frame(from = edges$INV_1, 
                          to = edges$INV_2, 
                          weight = edges$Jaccard, 
                          color = edges$color,
                          color.order = factor(edges$color, 
                                               levels = c(non_sig_color,
                                                          co_occ_color, 
                                                          rep_color
                                                          )))
  
  edges_net$weight <- abs(edges_net$weight) + 0.0001
  
  edges_net <- edges_net |> mutate(color2 = ifelse(color == co_occ_color, 
                                                   co_occ_color, NA),
                                   color3 = ifelse(color == rep_color, 
                                                   rep_color, NA),
                                   color4 = ifelse(color == "yellow", 
                                                   "yellow", NA),
                                   color5 = ifelse(color == "#5E4FA2", 
                                                   "#5E4FA2", NA),
                                   color1 = ifelse(color == non_sig_color, 
                                                   non_sig_color, NA))
  
  # Create a graph object
  network_global <- graph_from_data_frame(d = edges_net, 
                                                     vertices = nodes, 
                                                     directed = FALSE)
  
  # Add the color attribute to the edges
  E(network_global)$color <- edges_net$color
  
  # Simplify the graph to remove loops
  network_simplified_global <- igraph::simplify(network_global, 
                                                remove.loops = T, 
                                         edge.attr.comb = "first")
  
  # Calculate network metrics
  V(network_simplified_global)$community <- 
    cluster_fast_greedy(network_simplified_global)$membership
  V(network_simplified_global)$degree <- 
    degree(network_simplified_global, mode = "all")
  V(network_simplified_global)$closeness <- 
    closeness(network_simplified_global)
  V(network_simplified_global)$betweenness <- 
    betweenness(network_simplified_global)
  V(network_simplified_global)$eigenvector <- 
    eigen_centrality(network_simplified_global)$vector
  
  # Auxiliar for masking nodes that does not exist in a particular line
  
  chrom_number <- length(unique(V(network_simplified_global)$chromosome))
  
  palette <- hcl.colors(14, "Set3", rev = TRUE)
  palette <- setNames(palette, unique(V(network_simplified_global)$chromosome))
  
  V(network_simplified_global)$chrom_color1 <- 
    palette[V(network_simplified_global)$chromosome]
  V(network_simplified_global)$chrom_color1[
    V(network_simplified_global)$degree == 0] <- NA
  
  # Auxiliar for masking node labels that does not exist in a particular line
  V(network_simplified_global)$name[
    V(network_simplified_global)$degree == 0] <- ""
  
  # Convert to a tidygraph object
  tidy_network_global <- as_tbl_graph(network_simplified_global, 
                               directed = F)
  
  return(tidy_network_global)
}

sig_network_builder_lite <- function(nodes, meta_nodes, edges, type = "both",
                                co_occ_color = "blue",
                                rep_color = "red")
{
  # Create the nodes
  nodes <- data.frame(name = nodes)
  nodes$chromosome <- meta_nodes[, 'Chr'][match(
    sapply(strsplit(nodes$name, "_"), `[`, 2), meta_nodes$INV_ID)]
  
  nodes$name <- sapply(nodes$name, function(x) {
    paste0(strsplit(x, "_")[[1]][2], "_", strsplit(x, "_")[[1]][3])
  })
  
  # Initialize edges_net as NULL
  edges_net <- NULL
  
  # Change the names of the inversions
  
  edges$INV_1 <- sapply(edges$INV_1, function(x) {
    paste0(strsplit(x, "_")[[1]][2], "_", strsplit(x, "_")[[1]][3])
  })
  
  edges$INV_2 <- sapply(edges$INV_2, function(x) {
    paste0(strsplit(x, "_")[[1]][2], "_", strsplit(x, "_")[[1]][3])
  })
  
  # Create the edges
  edges <- edges %>% 
    rowwise() %>% 
    mutate(color = ifelse(Jaccard < 0, rep_color, co_occ_color))
  
  edges_net <- data.frame(from = edges$INV_1, 
                          to = edges$INV_2, 
                          weight = edges$Jaccard, 
                          color = edges$color)
  
  edges_net$weight <- abs(edges_net$weight) + 0.0001
  
  # Create a graph object
  network_global <- graph_from_data_frame(d = edges_net, 
                                          vertices = nodes, 
                                          directed = FALSE)
  
  # Add the color attribute to the edges
  E(network_global)$color <- edges_net$color
  
  # Simplify the graph to remove loops
  network_simplified_global <- igraph::simplify(network_global, 
                                                remove.loops = T, 
                                                edge.attr.comb = "first")
  
  # Calculate network metrics
  V(network_simplified_global)$degree <- 
    degree(network_simplified_global, mode = "all")
  V(network_simplified_global)$eigenvector <- 
    eigen_centrality(network_simplified_global)$vector
  
  # Auxiliar for masking nodes that does not exist in a particular line
  
  chrom_number <- length(unique(V(network_simplified_global)$chromosome))
  
  palette <- hcl.colors(14, "Set3", rev = TRUE)
  palette <- setNames(palette, unique(V(network_simplified_global)$chromosome))
  
  V(network_simplified_global)$chrom_color1 <- 
    palette[V(network_simplified_global)$chromosome]
  
  # Convert to a tidygraph object
  tidy_network_global <- as_tbl_graph(network_simplified_global, 
                                      directed = F)
  
  return(tidy_network_global)
}

# Apply the function to your data
networks <- mapply(sig_network_builder, 
                   final_result_split,
                   MoreArgs = list(nodes = nodes,
                                   meta_nodes = metadata),
                   SIMPLIFY = FALSE)

# 8.2) Plotting the networks ----

plot_network <- function(networks, 
                         edge_breaks = c(0, 0.01, 0.02, 0.05, 0.10, 0.14),
                         range = c(0, 5)) {
  for (i in 1:length(names(networks))) {
    pdf_file <- paste0("Plots/Networks/Network_", names(networks)[i], ".pdf")
    CairoPDF(file = pdf_file, width = 12, height = 10)
    
    p <- ggraph(networks[[i]], layout = "circle") +
      geom_edge_link(aes(edge_width = weight, 
                         edge_color = color1),
                     lineend = 'round',
                     linejoin = 'bevel') +
      geom_edge_link(aes(edge_width = weight, 
                         edge_color = color4),
                     lineend = 'round',
                     linejoin = 'bevel') +
      geom_edge_link(aes(edge_width = weight, 
                         edge_color = color5),
                     lineend = 'round',
                     linejoin = 'bevel') +
      geom_edge_link(aes(edge_width = weight, 
                         edge_color = color2),
                     lineend = 'round',
                     linejoin = 'bevel') +
      geom_edge_link(aes(edge_width = weight, 
                         edge_color = color3),
                     lineend = 'round',
                     linejoin = 'bevel') +
      geom_node_point(aes(color = chrom_color1, 
                          size = eigenvector),
                      na.rm = TRUE) +  # Fill the nodes with a color based on the chromosome attribute
      geom_node_text(aes(label = name), 
                     size = 5,
                     repel = T,
                     hjust = 0,
                     vjust = 0) +  # Use the sorted names for the labels
      scale_edge_width_continuous(range = range, 
                                  breaks = edge_breaks,
                                  labels = as.character(edge_breaks)) +
      scale_radius(range = c(1, 20), breaks = c(0.0, 0.15, 0.2, 0.55, 1.0),
                   labels = as.character(c(0.0, 0.15, 0.2, 0.55, 1.0))) +
      scale_edge_color_identity() +
      scale_color_identity() +  # Use a color palette for the ring
      theme_graph()
    
    print(p)  # Print the plot to the PDF file
    dev.off()
  }
}

plot_network_lite <- function(networks, 
                         edge_breaks = c(0, 0.01, 0.02, 0.05, 0.10, 0.14),
                         range = c(0, 10)) {
  for (i in 1:length(names(networks))) {
    pdf_file <- paste0("Plots/Networks/Network_", names(networks)[i], ".pdf")
    CairoPDF(file = pdf_file, width = 8, height = 6)
    
    p <- ggraph(networks[[i]], layout = "circle") +
      geom_edge_link(aes(edge_width = weight, 
                         edge_color = color),
                     lineend = 'round',
                     linejoin = 'bevel') +
      geom_node_point(aes(color = chrom_color1, 
                          size = eigenvector),
                      na.rm = TRUE) +  # Fill the nodes with a color based on the chromosome attribute
      geom_node_text(aes(label = name), 
                     size = 5,
                     repel = T,
                     hjust = 0,
                     vjust = 0) +  # Use the sorted names for the labels
      scale_edge_width_continuous(range = range, 
                                  breaks = edge_breaks,
                                  labels = as.character(edge_breaks)) +
      scale_radius(range = c(1, 15), breaks = c(0.0, 0.15, 0.2, 0.55, 1.0),
                   labels = as.character(c(0.0, 0.15, 0.2, 0.55, 1.0))) +
      scale_edge_color_identity() +
      scale_color_identity() +  # Use a color palette for the ring
      theme_graph()
    
    print(p)  # Print the plot to the PDF file
    dev.off()
  }
}

plot_network(networks)

# Legend for the colors of chromosomes:

# Create a data frame with the chromosomes and their corresponding colors
df <- data.frame(chromosome = unique(V(networks[[1]])$chromosome), 
                 color = hcl.colors(14, "Set3", rev = TRUE))

# Create a plot that only contains the legend
pdf("Plots/Networks/Chromosome_legend.pdf", width = 2, height = 8)

legend_plot <- ggplot(df, aes(x = 1, fill = chromosome)) +
  geom_tile(aes(y = chromosome)) +
  scale_fill_manual(values = df$color, 
                    guide = guide_legend(title = "Chromosome", 
                                         direction = "vertical", 
                                         title.position = "top", 
                                         label.position = "right")) +
  theme_void()

# Print the plot
print(legend_plot)

dev.off()

# 8.3) Study case of inversions 29, 32, and 40 ----

# Building a network with only 6 nodes (two dosages for the study cases)

# Function to count the occurrences of each inversion across the lines
count_inversions <- function(data_list) {
  # Initialize an empty list to store inversion counts
  inversion_counts <- list()
  
  # Loop through each matrix in the list
  for (name in names(data_list)) {
    # Get the row names (inversion names) for the current matrix
    inversions <- rownames(data_list[[name]])
    
    # Count the occurrences of each inversion
    for (inversion in inversions) {
      if (inversion %in% names(inversion_counts)) {
        inversion_counts[[inversion]] <- inversion_counts[[inversion]] + 1
      } else {
        inversion_counts[[inversion]] <- 1
      }
    }
  }
  
  # Convert the list to a data frame for easier manipulation
  inversion_counts_df <- data.frame(
    inversion = names(inversion_counts),
    count = unlist(inversion_counts)
  )
  
  return(inversion_counts_df)
}

# Function to select inversions that appear multiple times
select_frequent_inversions <- function(inversion_counts_df, min_count = 2) {
  # Filter inversions that appear at least min_count times
  frequent_inversions <- inversion_counts_df[
    inversion_counts_df$count >= min_count, ]
  
  return(frequent_inversions)
}

# Example usage
inversion_counts_df <- count_inversions(Data_d_l)
frequent_inversions <- select_frequent_inversions(inversion_counts_df, 
                                                  min_count = 7)

# Print the frequent inversions
print(frequent_inversions)

nodes_29_32_40 <- c("Inv_29_1", "Inv_29_2", "Inv_32_1", "Inv_32_2", 
                    "Inv_40_1", "Inv_40_2")

final_result_29_32_40 <- final_result |>
  filter(INV_1 %in% c("Inv_29_1", "Inv_29_2", "Inv_32_1", "Inv_32_2", 
                      "Inv_40_1", "Inv_40_2") & 
           INV_2 %in% c(c("Inv_29_1", "Inv_29_2", "Inv_32_1", "Inv_32_2", 
                          "Inv_40_1", "Inv_40_2")))

final_result_29_32_40 <- split(final_result_29_32_40, 
                               final_result_29_32_40$Line)

networks_29_32_40 <- mapply(sig_network_builder_lite, 
                          final_result_29_32_40,
                          MoreArgs = list(nodes = nodes_29_32_40,
                                          meta_nodes = metadata),
                          SIMPLIFY = FALSE)

names(networks_29_32_40) <- c("SUB_L_1034", "SUB_L_1192", "SUB_L_155",  
                              "SUB_L_444",  "SUB_L_502",
                              "SUB_L_541", "SUB_L_62", 
                              "SUB_L_664", "SUB_L_909")

plot_network_lite(networks_29_32_40, 
             edge_breaks = c(0.001, 0.008, 0.03, 0.05, 0.07, 0.10),
             range = c(1, 10))

# 9) Network analysis ----

# Function to compare two networks
compare_networks <- function(g1, g2) {
  # Check if the networks are isomorphic
  is_isomorphic <- igraph::isomorphic(g1, g2)
  
  # Compare network statistics
  avg_path_length_diff <- abs(igraph::mean_distance(g1, directed = FALSE) - igraph::mean_distance(g2, directed = FALSE))
  clustering_coeff_diff <- abs(igraph::transitivity(g1) - igraph::transitivity(g2))
  
  # Compare spectra
  spectrum_g1 <- eigen(igraph::as_adjacency_matrix(g1))
  spectrum_g2 <- eigen(igraph::as_adjacency_matrix(g2))
  spectrum_diff <- sum((spectrum_g1$values - spectrum_g2$values)^2)
  
  # Return a list with the comparison results
  return(list(is_isomorphic = is_isomorphic,
              avg_path_length_diff = avg_path_length_diff,
              clustering_coeff_diff = clustering_coeff_diff,
              spectrum_diff = spectrum_diff))
}

# Function to compare all pairs of networks
compare_all_pairs <- function(networks) {
  n <- length(networks)
  comparisons <- list()
  
  for (i in 1:(n-1)) {
    for (j in (i+1):n) {
      comparison <- compare_networks(networks[[i]], networks[[j]])
      comparisons[[paste0(names(networks)[i], "_vs_", names(networks)[j])]] <- comparison
    }
  }
  
  return(comparisons)
}

# Compare all pairs of networks
network_comparisons <- compare_all_pairs(networks)

# Condensing all the results into a single data frame

# Convert the list to a data frame
network_comparisons_df <- do.call(rbind, lapply(network_comparisons, function(x) {
  data.frame(
    is_isomorphic = x$is_isomorphic,
    avg_path_length_diff = x$avg_path_length_diff,
    clustering_coeff_diff = x$clustering_coeff_diff,
    spectrum_diff = x$spectrum_diff,
    stringsAsFactors = FALSE
  )
}))

# Add the comparison names as a new column
network_comparisons_df$comparison <- rownames(network_comparisons_df)

# Reset the row names
rownames(network_comparisons_df) <- NULL

# 9.1) Eigenvalue as a function of the number of genes in the inversion ----

# Extract the eigenvectors from each node in each line

eigenvectors <- lapply(networks, function(network) {
  V(network)$eigenvector
})

# Collapsing everything in a dataframe, keeping track of the line (names of the list)

eigenvectors_df <- do.call(rbind, lapply(names(eigenvectors), function(line) {
  data.frame(
    line = line,
    eigenvector = eigenvectors[[line]],
    stringsAsFactors = FALSE
  )
}))

# Incorporating the information of the nodes

eigenvectors_df$node <- rep(nodes, length(eigenvectors))

# Incorporating the information about the presence or absence of an inversion
# (using the degree information of the nodes from the networks)

degrees <- lapply(networks, function(network) {
  V(network)$degree
})

degrees_df <- do.call(rbind, lapply(names(degrees), function(line) {
  data.frame(
    line = line,
    degree = degrees[[line]],
    stringsAsFactors = FALSE
  )
}))

degrees_df$node <- rep(nodes, length(eigenvectors))

# Merging the eigenvectors and the degrees

eigenvectors_df <- merge(eigenvectors_df, degrees_df, by = c("line", "node"))

eigenvectors_df$Presence <- ifelse(eigenvectors_df$degree > 0, "Present", 
                                   "Absent")

# Incorporating the information about the number of genes in the inversion

eigenvectors_df$INV_ID <- sapply(strsplit(eigenvectors_df$node, "_"), `[`, 2)

eigenvectors_df <- eigenvectors_df |>
  mutate(INV_ID = as.numeric(INV_ID))

eigenvectors_df <- eigenvectors_df %>%
  left_join(genes_inv_metadata[, c("INV_ID", "INV_gene_number")], 
            by = "INV_ID")

# Filtering out the nodes with degree 0

eigenvectors_df <- eigenvectors_df |>
  filter(degree > 0)

# Taking the log of the number of genes in the inversion
eigenvectors_df$log_INV_gene_number <- log(eigenvectors_df$INV_gene_number)

eigenvectors_df <- eigenvectors_df |> 
  mutate(line = as.factor(line),
         node = as.factor(node))

# 9.2) Adjust a robust linear model ----

# Model selection based on AIC
# Fit a linear regression model
# Define the full model
full_model <- lm(eigenvector ~ ., data = eigenvectors_df[, c("eigenvector", 
                                                              "degree", 
                                                              "log_INV_gene_number", 
                                                              "INV_ID",
                                                              "INV_gene_number", 
                                                              "line")])

# Perform backward stepwise selection based on AIC
# Define the scope for stepwise selection

scope_list <- list(
  lower = ~ 1,  # The simplest model, only the intercept
  upper = ~ degree + log_INV_gene_number + INV_ID + INV_gene_number + line + 
    I(degree^2) + I(log_INV_gene_number^2) + I(INV_gene_number^2) + 
    degree:log_INV_gene_number + degree:INV_ID + degree:INV_gene_number + 
    degree:line + 
    log_INV_gene_number:INV_ID + log_INV_gene_number:INV_gene_number + 
    log_INV_gene_number:line + 
    INV_ID:INV_gene_number + INV_ID:line + INV_gene_number:line  # An example of a more complex model
)

# Perform stepwise selection based on AIC, correcting the scope usage
model_selected_aic <- step(full_model, direction = "both", 
                           scope = scope_list, steps = 100000)

final_model <- rlm(eigenvector ~ degree + INV_gene_number + line +
                     I(log_INV_gene_number^2), data = eigenvectors_df)

# Estimating the p-values

rob.pvals(final_model)

# Checking for multicollinearity

vif_result <- vif(final_model)

# Print the selected model
summary(final_model)

# Plotting the model (color coding degree)

eigenvectors_df$predicted_eigenvector <- predict(final_model, eigenvectors_df)

pdf("Plots/Genes_Eigenvector_Degree_Relationship.pdf", width = 8, height = 6)

ggplot(data = eigenvectors_df, aes(x = INV_gene_number, y = eigenvector, 
                                   color = degree)) +
  geom_point() + # Plot points
  geom_line(aes(y = predicted_eigenvector), color = "black") + # Add a regression line
  scale_color_gradient(low = "blue", high = "red") + # Color gradient for 'degree'
  labs(x = "Number of Genes (INV_gene_number)", y = "Eigenvector", 
       color = "Degree") +
  theme_bw() +
  theme(axis.text = element_text(size = 12, color = 'black'),
        axis.title = element_text(size = 14, color = 'black'),
        legend.text = element_text(size = 12, color = 'black'),
        legend.title = element_text(size = 14, color = 'black'),
        panel.grid = element_blank()) +
  xlim(0, 110)

dev.off()

# Plotting the model (color coding line)

pdf("Plots/Genes_Eigenvector_Line_Relationship.pdf", width = 8, height = 6)

ggplot(data = eigenvectors_df, aes(x = INV_gene_number, y = eigenvector)) +
  geom_point(aes(color = line)) + # Color points by 'line'
  geom_line(aes(y = predicted_eigenvector, 
                group = line, 
                color = line)) + # Separate regression lines for each 'line'
  labs(x = "Number of Genes (INV_gene_number)", y = "Eigenvector", color = "Line") +
  theme_bw() +
  theme(axis.text = element_text(size = 12, color = 'black'),
        axis.title = element_text(size = 14, color = 'black'),
        legend.text = element_text(size = 12, color = 'black'),
        legend.title = element_text(size = 14, color = 'black'),
        panel.grid = element_blank()) +
  xlim(0, 110)

dev.off()

##' <<<<<<<<<<<<<<<<<<<<<<<<<<<< End of the script >>>>>>>>>>>>>>>>>>>>>>>>>>>>>