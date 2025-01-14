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