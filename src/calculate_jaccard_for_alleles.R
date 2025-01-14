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
  
  jaccard.test(vector1, vector2, method = 'bootstrap')$statistics
}