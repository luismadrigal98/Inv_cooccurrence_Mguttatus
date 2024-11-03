##' @title Are you with me? Co-occurrence tests from community ecology can 
##' identify positive and negative epistasis between inversions in Mimulus 
##' guttatus
##' 
##' @description This script will explore the segregation distortion analysis
##' for the inversions considered in the study. Essentaily, it is an exploration
##' of the individual effect of each inversion in the survival of Mimulus 
##' guttatus.
##' 
##' @Author: Luis Javier Madriga-Roca & John K. Kelly
##' 
##' @Date: 09/03/2024
##' ____________________________________________________________________________
##' <<<<<<<<<<<<<<<<<<<<<<<<<<<< Start >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

## *****************************************************************************
## 1) Distortion segregation analysis ----
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
  #' refers to a particular line. Then, using the line and Inv_ID, the dosage 
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

observed_counts_deviants <- mapply(significant_inv$Line, significant_inv$Inv, 
                                   MoreArgs = list(dosage_list = Data_d_l),
                                   FUN = observed_count_extractor,
                                   SIMPLIFY = F)

# Merging individual data.frames into a unique one
observed_counts_deviants <- do.call(rbind, observed_counts_deviants)

# Creating the functions to calculate the X^2_sub1 and X^2_sub2

X2_sub1_calculator <- function(RR, RA, AA)
{
  #' Allele frequency homogeneity.
  #' 
  #' This function is designed to check the balanced frequency of both alleles
  #' (p = q, X^2_sub1) for the sample analyzed. Is one of the measurements that 
  #' are going to be used to know the pattern of selection (gametyc or zygotic).
  #' 
  #' @Arguments: RR: Homozygous for the reference allele
  #'             RA: Heterozygous
  #'             AA: Homozygous for the alternative allele (inversion)
  #'
  #' @Returns: X^2_sub1 value
  #' 
  #' @Note: df of the test:
  #' 2 categories (two alleles) - 1 parameter estimated (p) - 1 = 1 df
  #' ___________________________________________________________________________

  n <- RR + RA + AA # Total number of plants
  p <- (2 * RR + RA) / (2 * n) # Calculate the frequency of the reference allele
  q <- 1 - p # Calculate the frequency of the alternative allele
  
  X2_sub1 <- ((2 * n * p - n) ^ 2 + (2 * n * q - n) ^ 2) / n
  
  return(X2_sub1)
}

X2_sub2_calculator <- function(RR, RA, AA)
{
  #' F2 distribution.
  #' 
  #' This function is designed to check the distribution of different genotypic
  #' frequencies. Is one of the measurements that are going to be used to know 
  #' the pattern of selection (gametyc or zygotic).
  #' 
  #' @Arguments: RR: Homozygous for the reference allele
  #'             RA: Heterozygous
  #'             AA: Homozygous for the alternative allele (inversion)
  #'
  #' @Returns: X^2_sub1 value
  #' 
  #' @Note: df of the test:
  #' 3 categories (three genotypes) - 1 parameter estimated (p) - 1 = 1 df
  #' ___________________________________________________________________________

  n <- RR + RA + AA # Total number of plants
  p <- (2 * RR + RA) / (2 * n) # Calculate the frequency of the reference allele
  q <- 1 - p # Calculate the frequency of the alternative allele
  
  X2_sub2 <- ((AA - n * p ^ 2) ^ 2) / (n * p ^ 2)  + 
    ((RA - 2 * n * p * q) ^ 2) / (2 * n * p * q) +
    ((RR - n * q ^ 2) ^ 2) / (n * q ^ 2)
  
  return(X2_sub2)
}

# Including the measures in the data.frame

observed_counts_deviants$X2_sub1 <- mapply(observed_counts_deviants$RR, 
                                          observed_counts_deviants$RA, 
                                          observed_counts_deviants$AA, 
                                          FUN = X2_sub1_calculator)

observed_counts_deviants$X2_sub2 <- mapply(observed_counts_deviants$RR,
                                          observed_counts_deviants$RA,
                                          observed_counts_deviants$AA,
                                          FUN = X2_sub2_calculator)

# ~ Classifying the selection based on the X^2_sub1 and X^2_sub2 values ----

observed_counts_deviants$X2_sub1_p <- pchisq(observed_counts_deviants$X2_sub1, 
                                             df = 1, lower.tail = F)

observed_counts_deviants$X2_sub2_p <- pchisq(observed_counts_deviants$X2_sub2,
                                             df = 1, lower.tail = F)

observed_counts_deviants$Selection <- 
  ifelse(observed_counts_deviants$X2_sub1_p < 0.05 & 
           observed_counts_deviants$X2_sub2_p < 0.05, "Zygotic",
         ifelse(observed_counts_deviants$X2_sub1_p < 0.05 & 
                  observed_counts_deviants$X2_sub2_p > 0.05, "Gametic",
                "Zygotic"))

# ~ Calculating the segregation distortion value (SDV) ----

#' The segregation distortion value is the natural logarithm of the omnibus X2
#' p-value.

significant_inv <- significant_inv |>
  mutate(SDV = -log(p_value))

#' Adding the distortion segregation value to the dataframe with all inversions

goodness_of_fit_df <- goodness_of_fit_df |>
 mutate(SDV = -log(p_value))

# Exporting the results

write.csv(goodness_of_fit_df, 
          file = "/home/l338m483/scratch/Cooccurrence_Inv/R_directory/Results/goodness_of_fit_individual_effect.csv",
          row.names = F)

write.csv(significant_inv, 
          file = "/home/l338m483/scratch/Cooccurrence_Inv/R_directory/Results/significant_inv_individual_effect.csv",
          row.names = F)

# Stop the cluster after use
stopCluster(cl)

##' <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<< End >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>