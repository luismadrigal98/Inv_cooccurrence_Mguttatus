permutation_test_motifs <- function(obs_res, cross, motif_name, n_permutations) 
{
  #' Test the significance of a motif in a network using permutation tests
  #' 
  #' @param obs_res A data frame with the observed results
  #' @param cross The name of the network
  #' @param motif_name The name of the motif
  #' @param n_permutations The number of permutations to perform
  #' 
  #' @return A data frame with the obs_res of the permutation test
  #' ___________________________________________________________________________
  
  # Get original network and motif
  network <- networks[[cross]]
  motif <- motifs[[motif_name]]
  
  # Get observed count
  observed <- obs_res$counts[obs_res$Cross == cross & 
                               obs_res$motif_names == motif_name]
  
  # Generate null distribution through permutations
  null_counts <- replicate(n_permutations, {
    # Create rewired network preserving degree distribution
    rewired <- rewire(network, keeping_degseq(niter = gsize(network) * 10))
    # Count motifs in permuted network
    subgraphs <- graph.get.subisomorph(rewired, motif)
    length(subgraphs)
  })
  
  # Calculate statistics
  z_score <- (observed - mean(null_counts)) / sd(null_counts)
  p_value <- mean(null_counts >= observed)  # One-sided test
  
  return(data.frame(
    Cross = cross,
    motif_names = motif_name,
    observed = observed,
    mean_null = mean(null_counts),
    sd_null = sd(null_counts),
    z_score = z_score,
    p_value = p_value
  ))
}