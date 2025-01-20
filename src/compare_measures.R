compare_measures <- function(freq_range,
                             n_samples = 400,
                             n_reps = 100,
                             association_strength = 0.5) 
{
  #' Compare LD, Jaccard and Affinity measures
  #' 
  #' @param freq_range A vector of allele frequencies to simulate
  #' @param n_samples Number of samples to simulate
  #' @param n_reps Number of replicates to simulate
  #' @param association_strength The strength of association between alleles
  #' 
  #' @return A data frame with the results of the simulation
  #' ___________________________________________________________________________
  
  # Create parameter combinations
  param_grid <- expand.grid(freq1 = freq_range,
                            freq2 = freq_range,
                            rep = 1:n_reps)
  
  # Run in parallel
  results_list <- foreach(i = 1:nrow(param_grid), 
                          .packages = c('foreach', 'copula'),
                          .export = c('simulate_frequency_data',
                                      'geno_to_allele',
                                      'calculate_LD',
                                      'calculate_jaccard_for_alleles',
                                      'calculate_affinity_for_alleles',
                                      'normalCopula',
                                      'cCopula',
                                      'jaccard.test')) %dopar% 
    {
      data <- simulate_frequency_data(n_samples, param_grid$freq1[i], 
                                      param_grid$freq2[i], 
                                      association_strength)
      
      # Getting the allele vectors 
      allele1 <- geno_to_allele(data$inv1)
      allele2 <- geno_to_allele(data$inv2)
      allele_data <- data.frame(allele1, allele2)
      
      # Calculate all metrics
      ld_res <- calculate_LD(allele1, allele2)
      cJaccard <- calculate_jaccard_for_alleles(allele1, allele2)
      affinity <- calculate_affinity_for_alleles(allele_data)
      
      # Return results as vector
      c(
        freq1 = param_grid$freq1[i],
        freq2 = param_grid$freq2[i],
        rep = param_grid$rep[i],
        ld_res['D'],
        ld_res['D_prime'],
        cJaccard = cJaccard,
        Affinity = affinity
      )
    }
  
  # Convert list to dataframe
  results <- do.call(rbind, results_list)
  results <- as.data.frame(results)
  
  return(results)
}