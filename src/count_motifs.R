count_motifs <- function(network, motifs)
{
  #' This function will count the presence of each motif in an individual network
  #' 
  #' @param network An igraph object
  #' @param motifs A list of igraph objects representing the motifs of interest
  #' 
  #' @return A data frame with the motifs and their counts
  #' ___________________________________________________________________________
  
  motifs_names <- names(motifs)
  
  counts <- sapply(motifs, function(x)
    {
      count_subgraph_isomorphisms(pattern = x, target = network) 
    }, simplify = T)
  
  data.frame(motifs_names, counts)
}