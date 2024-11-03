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