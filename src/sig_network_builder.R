sig_network_builder <- function(nodes, meta_nodes, edges, type = "both",
                                co_occ_color = "blue",
                                rep_color = "red",
                                non_sig_color = "grey")
{
  #' Network builder function
  #' 
  #' This function builds a network from a list of nodes and edges.
  #' 
  #' @param nodes A character vector with the names of the nodes.
  #' @param meta_nodes A data frame with the metadata of the nodes.
  #' @param edges A data frame with the edges of the network.
  #' @param type A character vector with the type of network to be built.
  #' @param co_occ_color A character vector with the color of the edges that represent co-occurrence.
  #' @param rep_color A character vector with the color of the edges that represent repulsion.
  #' @param non_sig_color A character vector with the color of the edges that are not significant.
  #' 
  #' @return A tidygraph object with the network.
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
    mutate(color = ifelse(p_X2_global > 0.05 & J_p > 0.05, 
                          non_sig_color, 
                          ifelse(all(Jaccard > 0, J_p < 0.05), 
                                 co_occ_color, 
                                 ifelse(all(Jaccard < 0, J_p < 0.05), 
                                        rep_color,
                                        ifelse(all(Jaccard > 0, J_p > 0.05), 
                                               "yellow", 
                                               "#5E4FA2")))))
  
  edges_net <- data.frame(from = edges$INV_1, 
                          to = edges$INV_2, 
                          weight = edges$Jaccard, 
                          color = edges$color,
                          color.order = factor(edges$color, 
                                               levels = c(non_sig_color,
                                                          co_occ_color, 
                                                          rep_color
                                               )))
  
  edges_net$weight <- abs(edges_net$weight) + 0.0001
  
  edges_net <- edges_net |> mutate(color2 = ifelse(color == co_occ_color, 
                                                   co_occ_color, NA),
                                   color3 = ifelse(color == rep_color, 
                                                   rep_color, NA),
                                   color4 = ifelse(color == "yellow", 
                                                   "yellow", NA),
                                   color5 = ifelse(color == "#5E4FA2", 
                                                   "#5E4FA2", NA),
                                   color1 = ifelse(color == non_sig_color, 
                                                   non_sig_color, NA))
  
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
  V(network_simplified_global)$community <- 
    cluster_fast_greedy(network_simplified_global)$membership
  V(network_simplified_global)$degree <- 
    degree(network_simplified_global, mode = "all")
  V(network_simplified_global)$closeness <- 
    closeness(network_simplified_global)
  V(network_simplified_global)$betweenness <- 
    betweenness(network_simplified_global)
  V(network_simplified_global)$eigenvector <- 
    eigen_centrality(network_simplified_global)$vector
  
  # Auxiliar for masking nodes that does not exist in a particular line
  
  chrom_number <- length(unique(V(network_simplified_global)$chromosome))
  
  palette <- hcl.colors(14, "Set3", rev = TRUE)
  palette <- setNames(palette, unique(V(network_simplified_global)$chromosome))
  
  V(network_simplified_global)$chrom_color1 <- 
    palette[V(network_simplified_global)$chromosome]
  V(network_simplified_global)$chrom_color1[
    V(network_simplified_global)$degree == 0] <- NA
  
  # Auxiliar for masking node labels that does not exist in a particular line
  V(network_simplified_global)$name[
    V(network_simplified_global)$degree == 0] <- ""
  
  # Convert to a tidygraph object
  tidy_network_global <- as_tbl_graph(network_simplified_global, 
                                      directed = F)
  
  return(tidy_network_global)
}