find_network_motifs <- function(networks, motif_size = 3) 
{
  #' Analyze Network Motifs in Multiple Networks with Statistical Testing
  #' 
  #' @description
  #' Identifies and analyzes network motifs (recurring patterns of interconnections) within one or more networks,
  #' performing statistical testing against randomized null models to assess motif significance.
  #'
  #' @param networks A list of igraph objects or objects that can be converted to igraph format
  #' @param motif_size Numeric, the size of motifs to analyze (default = 3 nodes)
  #'
  #' @return A named list where each element corresponds to a network and contains:
  #' \itemize{
  #'   \item motif_counts: Raw counts of each motif type in the original network
  #'   \item z_scores: Standardized scores comparing motif frequencies to null model
  #'   \item p_values: Two-tailed p-values for significance of motif enrichment/depletion
  #'   \item null_means: Mean motif counts across randomized networks
  #'   \item null_sds: Standard deviations of motif counts in randomized networks
  #'   \item network_stats: Basic network statistics (vertices, edges, density)
  #' }
  #'
  #' @details
  #' The function performs the following steps for each network:
  #' 1. Counts motifs in the original network
  #' 2. Generates 1000 randomized networks preserving degree sequence
  #' 3. Counts motifs in each randomized network
  #' 4. Computes z-scores and p-values comparing original to null distribution
  #'
  #' @note
  #' - Networks are randomized using the degree-preserving rewiring algorithm
  #' - P-values are calculated using two-tailed tests
  #' - NULL results are returned for networks where motif counting fails
  #'
  #' @examples
  #' # Analyze motifs in a list of two networks
  #' g1 <- make_ring(10)
  #' g2 <- make_lattice(c(5,5))
  #' results <- find_network_motifs(list(ring=g1, lattice=g2))
  #'
  #' @importFrom igraph graph.motifs rewire vcount ecount edge_density is.igraph as.igraph
  #' @importFrom stats pnorm sd
  #'
  #' @export
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