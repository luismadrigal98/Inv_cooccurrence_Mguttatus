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