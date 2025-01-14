calculate_LD <- function(vector1, vector2)
{
  #' This function will calculate the statistics D and D' between two corresponding
  #' vectors of presence/absence between two alleles of two different loci. The
  #' presence is for the alleles of interest (e. g. the presence of the recessive
  #' allele in locus 1 and the presence of the recessive allele in locus 2).
  #' 
  #' @param vector1 Numeric vector. Presence/absence of the first allele.
  #' @param vector2 Numeric vector. Presence/absence of the second allele.
  #' 
  #' @return A named vector with the statistics D and D'.
  #' 
  #' @note
  #' It is necessary to transfor the genotype vector into an allele vector
  #' ___________________________________________________________________________
  
  N <- length(vector1)
  
  # Getting the frequency of the alleles (sum presence/absence vector and divide
  # by total)
  
  pA = sum(vector1) / N
  pB = sum(vector2) / N
  
  # Getting the frequency of the haplotype AB
  
  mask <- (vector1 == 1) & (vector2 == 1) # Getting those cases where both alleles
  # are present
  
  pAB <- sum(mask) / N
  
  D <- pAB - pA * pB
  
  # Calculating D_prime
  
  if (D < 0)
  {
    denominator <- max(-(pA*pB), -((1 - pB) * (1 - pA)))
  }
  
  else
  {
    denominator <- min(pA * (1 - pB), pB * (1 - pA))
  }
  
  D_prime <- D / denominator
  
  res <- c(D, D_prime)
  names(res) <- c('D', 'D_prime')
  
  return(res)
}