plot_network_lite <- function(networks, 
                              edge_breaks = c(0, 0.01, 0.02, 0.05, 0.10, 0.14),
                              range = c(0, 10)) 
{
  #' Plot a network using the ggraph package
  #' 
  #' @param networks A list of igraph objects
  #' @param edge_breaks A vector of breaks for the edge width scale
  #' @param range A vector of the range for the edge width scale
  #' 
  #' @return A PDF file with the network plot
  #' ___________________________________________________________________________
  
  for (i in 1:length(names(networks))) {
    pdf_file <- paste0("Results/Plots/Networks/Network_", names(networks)[i], ".pdf")
    CairoPDF(file = pdf_file, width = 8, height = 6)
    
    p <- ggraph(networks[[i]], layout = "circle") +
      geom_edge_link(aes(edge_width = weight, 
                         edge_color = color),
                     lineend = 'round',
                     linejoin = 'bevel') +
      geom_node_point(aes(color = chrom_color1, 
                          size = eigenvector),
                      na.rm = TRUE) +  # Fill the nodes with a color based on the chromosome attribute
      geom_node_text(aes(label = name), 
                     size = 5,
                     repel = T,
                     hjust = 0,
                     vjust = 0) +  # Use the sorted names for the labels
      scale_edge_width_continuous(range = range, 
                                  breaks = edge_breaks,
                                  labels = as.character(edge_breaks)) +
      scale_radius(range = c(1, 15), breaks = c(0.0, 0.15, 0.2, 0.55, 1.0),
                   labels = as.character(c(0.0, 0.15, 0.2, 0.55, 1.0))) +
      scale_edge_color_identity() +
      scale_color_identity() +  # Use a color palette for the ring
      theme_graph()
    
    print(p)  # Print the plot to the PDF file
    dev.off()
  }
}