permutation_test_motifs <- function(obs_res, networks, motifs,
                                    cross, motif_name, 
                                    n_permutations, parallel = TRUE) 
{
  #' Test the significance of a motif in a network using permutation tests
  #' @param obs_res A data frame with the observed results
  #' @param networks A list of networks
  #' @param motifs A list of motifs
  #' @param cross The name of the network
  #' @param motif_name The name of the motif
  #' @param n_permutations The number of permutations to perform
  #' @param parallel Logical indicating whether to run the permutation in parallel
  
  # Pre-fetch network and motif
  network <- networks[[cross]]
  motif <- motifs[[motif_name]]
  observed <- obs_res$counts[obs_res$Cross == cross & 
                               obs_res$motif_names == motif_name]
  
  # Calculate degree sequence once
  degree_seq <- degree(network)
  
  if (parallel) {
    # Use mclapply for parallel processing
    null_counts <- unlist(parallel::mclapply(1:n_permutations, function(i) {
      rand_graph <- sample_degseq(degree_seq, method = "simple")
      rand_graph <- simplify(rand_graph)
      count_subgraph_isomorphisms(pattern = motif, target = rand_graph)
    }, mc.cores = parallel::detectCores() - 1))
  } else {
    null_counts <- numeric(n_permutations)
    for(i in seq_len(n_permutations)) {
      rand_graph <- sample_degseq(degree_seq, method = "simple")
      rand_graph <- simplify(rand_graph)
      null_counts[i] <- count_subgraph_isomorphisms(pattern = motif, 
                                                    target = rand_graph)
    }
  }
  
  # Calculate statistics
  mean_null <- mean(null_counts)
  sd_null <- sd(null_counts)
  sd_null <- ifelse(sd_null == 0, 1e-10, sd_null)
  z_score <- (observed - mean_null) / sd_null
  p_value <- sum(null_counts >= observed) / n_permutations
  
  data.frame(
    Cross = cross,
    motif_names = motif_name,
    observed = observed,
    mean_null = mean_null,
    sd_null = sd_null,
    z_score = z_score,
    p_value = p_value
  )
}