calculate_affinity_for_alleles <- function(data, which_dim = 'col')
{
  #' Calculate the affinity for each allele in the data.
  #' 
  #' @param data A data frame of allele data. Each column represent one of the
  #' alleles of interest
  #' 
  #' @param which_dim A string indicating whether the affinity should be calculated
  #' for the rows or columns of the data frame. Default is 'col'.
  #' 
  #' @return A numeric vector of affinities for each allele in the data.
  #' ___________________________________________________________________________
  
  CooccurrenceAffinity::affinity(data = data, 
                                 row.or.col = which_dim)$all$alpha_mle
}