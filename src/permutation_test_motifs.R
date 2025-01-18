permutation_test_motifs <- function(obs_res, 
                                    cross, motif_name, 
                                    n_permutations, parallel = T) 
{
  #' Test the significance of a motif in a network using permutation tests
  #' 
  #' @param obs_res A data frame with the observed results
  #' @param cross The name of the network
  #' @param motif_name The name of the motif
  #' @param n_permutations The number of permutations to perform
  #' @param parallel Logical indicating whether to run the permutation in parallel
  #' 
  #' @return A data frame with the obs_res of the permutation test
  #' ___________________________________________________________________________
  
  # Pre-fetch network and motif outside the loop
  network <- networks[[cross]]
  motif <- motifs[[motif_name]]
  observed <- obs_res$counts[obs_res$Cross == cross & 
                               obs_res$motif_names == motif_name]
  
  # Pre-allocate vector for results
  null_counts <- numeric(n_permutations)
  
  # Calculate rewire iterations once
  rewire_iters <- gsize(network) * 10
  
  if (parallel)
  {
    n_cores <- parallel::detectCores() - 1
      
    null_counts <- unlist(mclapply(1:n_permutations, function(i) {
      rewired <- rewire(network, keeping_degseq(niter = rewire_iters))
      count_subgraph_isomorphisms(pattern = motif, target = rewired)
    }, mc.cores = n_cores))
  }
  
  else
  {
    for(i in seq_len(n_permutations)) {
      rewired <- rewire(network, keeping_degseq(niter = rewire_iters))
      null_counts[i] <- count_subgraph_isomorphisms(pattern = motif, target = rewired)
    }  
  }
  
  # Vectorized calculations for statistics
  mean_null <- mean(null_counts)
  sd_null <- sd(null_counts)
  z_score <- (observed - mean_null) / sd_null
  p_value <- sum(null_counts >= observed) / n_permutations  # More precise calculation
  
  # Create result data frame in one step
  result <- data.frame(
    Cross = cross,
    motif_names = motif_name,
    observed = observed,
    mean_null = mean_null,
    sd_null = sd_null,
    z_score = z_score,
    p_value = p_value
  )
  
  return(result)
}