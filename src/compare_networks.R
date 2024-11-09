compare_networks <- function(g1, g2) 
{
  #' Compare two networks
  #' 
  #' This function compares two networks in terms of isomorphism, average path 
  #' length, clustering coefficient, and spectrum.
  #' 
  #' @param g1 An igraph object representing the first network.
  #' @param g2 An igraph object representing the second network.
  #' 
  #' @return A list with the comparison results.
  #' ___________________________________________________________________________
  
  # Check if the networks are isomorphic
  is_isomorphic <- igraph::isomorphic(g1, g2)
  
  # Compare network statistics
  avg_path_length_diff <- abs(igraph::mean_distance(g1, directed = FALSE) - 
                                igraph::mean_distance(g2, directed = FALSE))
  clustering_coeff_diff <- abs(igraph::transitivity(g1) - 
                                 igraph::transitivity(g2))
  
  # Compare spectra
  spectrum_g1 <- eigen(igraph::as_adjacency_matrix(g1))
  spectrum_g2 <- eigen(igraph::as_adjacency_matrix(g2))
  spectrum_diff <- sum((spectrum_g1$values - spectrum_g2$values)^2)
  
  # Return a list with the comparison results
  return(list(is_isomorphic = is_isomorphic,
              avg_path_length_diff = avg_path_length_diff,
              clustering_coeff_diff = clustering_coeff_diff,
              spectrum_diff = spectrum_diff))
}