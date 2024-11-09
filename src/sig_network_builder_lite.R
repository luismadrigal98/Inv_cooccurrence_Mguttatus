sig_network_builder_lite <- function(nodes, meta_nodes, edges, type = "both",
                                     co_occ_color = "blue",
                                     rep_color = "red")
{
  #' Create a network from a list of nodes and edges
  #' 
  #' This function creates a network from a list of nodes and edges. The nodes
  #' are the inversions and the edges are the Jaccard index between them. This is
  #' a more light version dedicated to zoom-in certain contrasts.
  #' 
  #' @param nodes A character vector with the names of the nodes
  #' @param meta_nodes A data frame with the metadata of the nodes
  #' @param edges A data frame with the edges
  #' @param type A character vector with the type of the edges
  #' @param co_occ_color A character vector with the color of the co-occurrence edges
  #' @param rep_color A character vector with the color of the repulsion edges
  #' 
  #' @return A tidygraph object
  #' ___________________________________________________________________________
  
  # Create the nodes
  nodes <- data.frame(name = nodes)
  nodes$chromosome <- meta_nodes[, 'Chr'][match(
    sapply(strsplit(nodes$name, "_"), `[`, 2), meta_nodes$INV_ID)]
  
  nodes$name <- sapply(nodes$name, function(x) {
    paste0(strsplit(x, "_")[[1]][2], "_", strsplit(x, "_")[[1]][3])
  })
  
  # Initialize edges_net as NULL
  edges_net <- NULL
  
  # Change the names of the inversions
  
  edges$INV_1 <- sapply(edges$INV_1, function(x) {
    paste0(strsplit(x, "_")[[1]][2], "_", strsplit(x, "_")[[1]][3])
  })
  
  edges$INV_2 <- sapply(edges$INV_2, function(x) {
    paste0(strsplit(x, "_")[[1]][2], "_", strsplit(x, "_")[[1]][3])
  })
  
  # Create the edges
  edges <- edges %>% 
    rowwise() %>% 
    mutate(color = ifelse(Jaccard < 0, rep_color, co_occ_color))
  
  edges_net <- data.frame(from = edges$INV_1, 
                          to = edges$INV_2, 
                          weight = edges$Jaccard, 
                          color = edges$color)
  
  edges_net$weight <- abs(edges_net$weight) + 0.0001
  
  # Create a graph object
  network_global <- graph_from_data_frame(d = edges_net, 
                                          vertices = nodes, 
                                          directed = FALSE)
  
  # Add the color attribute to the edges
  E(network_global)$color <- edges_net$color
  
  # Simplify the graph to remove loops
  network_simplified_global <- igraph::simplify(network_global, 
                                                remove.loops = T, 
                                                edge.attr.comb = "first")
  
  # Calculate network metrics
  V(network_simplified_global)$degree <- 
    degree(network_simplified_global, mode = "all")
  V(network_simplified_global)$eigenvector <- 
    eigen_centrality(network_simplified_global)$vector
  
  # Auxiliar for masking nodes that does not exist in a particular line
  
  chrom_number <- length(unique(V(network_simplified_global)$chromosome))
  
  palette <- hcl.colors(14, "Set3", rev = TRUE)
  palette <- setNames(palette, unique(V(network_simplified_global)$chromosome))
  
  V(network_simplified_global)$chrom_color1 <- 
    palette[V(network_simplified_global)$chromosome]
  
  # Convert to a tidygraph object
  tidy_network_global <- as_tbl_graph(network_simplified_global, 
                                      directed = F)
  
  return(tidy_network_global)
}