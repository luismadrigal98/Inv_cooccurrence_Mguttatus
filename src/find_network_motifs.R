find_network_motifs <- function(networks, motif_size = 3) 
{
  #' Count network motifs and calculate z-scores
  #' 
  #' This function counts the number of network motifs in a set of networks and
  #' calculates z-scores by comparing to random networks.
  #' 
  #' @param networks A list of igraph objects.
  #' 
  #' @param motif_size The size of the motifs to count.
  #' 
  #' @return A list of lists, each containing the motif counts and z-scores for
  #' each network.
  #' ___________________________________________________________________________
  
  lapply(networks, function(g) {
    # Get motif counts and z-scores
    motifs <- count_motifs(g, size = motif_size)
    
    # Compare to random networks
    null_motifs <- replicate(1000, {
      g_random <- rewire(g, keeping_degseq(niter = vcount(g) * 10))
      count_motifs(g_random, size = motif_size)
    })
    
    # Calculate z-scores
    z_scores <- (motifs - rowMeans(null_motifs)) / apply(null_motifs, 1, sd)
    
    list(
      motif_counts = motifs,
      z_scores = z_scores
    )
  })
}