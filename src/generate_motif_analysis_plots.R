generate_motif_analysis_plots <- function(data) 
{
  #' Generate motif analysis plots
  #' 
  #' This function generates a series of plots to summarize the results of the
  #' motif analysis. The plots include a heatmap of the significance of the
  #' motifs, a barplot of the enrichment of the motifs, a heatmap of the z-scores
  #' of the motifs, and a barplot of the counts of the motifs.
  #' 
  #' @param data A data frame containing the results of the motif analysis
  #' 
  #' @return A list containing the plots
  #' ___________________________________________________________________________
  
  p1 <- plot_significance_heatmap(data)
  p2 <- plot_enrichment(data)
  p3 <- plot_zscores(data)
  p4 <- plot_counts_comparison(data)
  p <- gridExtra::grid.arrange(p1, p2, p3, p4, nrow = 2)
  
  # Save the plots
  ggsave("./Results/Plots/Motif_analysis_results.pdf", p, 
         width = 10, height = 8)
  
  # Return the plots as a list
  list(
    heatmap = p1,
    enrichment = p2,
    zscores = p3,
    counts = p4
  )
}