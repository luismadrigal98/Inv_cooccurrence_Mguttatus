simulate_frequency_data <- function(n_samples, freq1, freq2, association_strength) {
  #' Simulate genotype data for two inversions with specified frequencies and association strength
  #'
  #' @param n_samples Integer. Number of samples to generate.
  #' @param freq1 Numeric. Frequency of the first inversion (0 to 1).
  #' @param freq2 Numeric. Frequency of the second inversion (0 to 1).
  #' @param association_strength Numeric. Strength of association between inversions (-1 to 1).
  #'
  #' @return Data frame with two columns representing genotypes of the two inversions.
  #' @examples
  #' simulate_frequency_data(100, 0.3, 0.5, 0.7)
  
  # Validate inputs
  if (!is.numeric(n_samples) || n_samples <= 0 || n_samples != round(n_samples)) {
    stop("n_samples must be a positive integer.")
  }
  if (!is.numeric(freq1) || freq1 < 0 || freq1 > 1) {
    stop("freq1 must be a numeric value between 0 and 1.")
  }
  if (!is.numeric(freq2) || freq2 < 0 || freq2 > 1) {
    stop("freq2 must be a numeric value between 0 and 1.")
  }
  if (!is.numeric(association_strength) || association_strength < -1 || association_strength > 1) {
    stop("association_strength must be a numeric value between -1 and 1.")
  }
  
  # Base independent probabilities (inv is the alternative allele, q)
  prob1 <- c((1 - freq1)^2, 2 * freq1 * (1 - freq1), freq1^2)  # HWE proportions
  prob2 <- c((1 - freq2)^2, 2 * freq2 * (1 - freq2), freq2^2)
  
  # Generate independent uniform random variables
  u <- runif(n_samples)
  v <- runif(n_samples)
  
  # Introduce correlation using a copula
  cop <- normalCopula(param = association_strength, dim = 2, dispstr = "un")
  uv <- cCopula(cbind(u, v), copula = cop, inverse = TRUE)
  
  # Convert uniform random variables to genotypes
  geno1 <- cut(uv[, 1], breaks = c(0, cumsum(prob1)), labels = 0:2, include.lowest = TRUE)
  geno2 <- cut(uv[, 2], breaks = c(0, cumsum(prob2)), labels = 0:2, include.lowest = TRUE)
  geno1 <- as.numeric(as.character(geno1))
  geno2 <- as.numeric(as.character(geno2))
  
  data.frame(inv1 = geno1, inv2 = geno2)
}

geno_to_allele <- function(geno_vector) {
  #' This function will take a genotype vector with 0, 1, or 2 describing the
  #' number of copies of the alternative allele, and it will return an allele
  #' vector describing the presence / absence of the recessive / alternative 
  #' allele.
  #' 
  #' @param geno_vector Numeric vector. Genotype vector.
  #' 
  #' @return Numeric vector. Allele vector.
  #' 
  #' @example 
  #' geno_vector <- c(0, 1, 2, 1, 0)
  #' allele_vector <- geno_to_allele(geno_vector)
  #' print(allele_vector)
  #' ___________________________________________________________________________
  
  # Preallocate the allele vector to twice the length of the genotype vector
  allele_vector <- integer(length(geno_vector) * 2)
  
  # Fill in the allele vector based on the genotype vector
  for (i in seq_along(geno_vector)) {
    index <- (i - 1) * 2 + 1
    if (geno_vector[i] == 0) {
      allele_vector[index] <- 0
      allele_vector[index + 1] <- 0
    } else if (geno_vector[i] == 1) {
      allele_vector[index] <- 0
      allele_vector[index + 1] <- 1
    } else if (geno_vector[i] == 2) {
      allele_vector[index] <- 1
      allele_vector[index + 1] <- 1
    }
  }
  
  return(allele_vector)
}

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
  
  mask <- (vector1 == 1) == (vector2 == 1) # Getting those cases where both alleles
                                           # are together
  
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

calculate_jaccard_for_alleles <- function(vector1, vector2)
{
  #' This function will calculate the Centered Jaccard-Tanimoto index for a given
  #' pair of alleles.
  #' 
  #' @param vector1 Numeric vector. Presence/absence of the first allele.
  #' @param vector2 Numeric vector. Presence/absence of the second allele.
  #' 
  #' @return Numeric. Jaccard-Tanimoto index.
  #' ___________________________________________________________________________
  
  jaccard.test(vector1, vector2, method = 'mca')$statistics
}
