plot_significance_heatmap <- function(data) 
{
  #' Plot a heatmap of motif significance across networks
  #' 
  #' @param data A data frame with columns motif_names, Cross, p_adjusted
  #' @return A ggplot2 object
  #' ___________________________________________________________________________
  
  ggplot(data, aes(x = motif_names, y = Cross, fill = -log10(p_adjusted))) +
    geom_tile() +
    scale_fill_viridis(name = "-log10(adj.P)") +
    theme_bw() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1),
          axis.text = element_text(color = 'black'),
          panel.grid = element_blank()) +
    labs(title = "Motif Significance Across Networks",
         x = "Motif", y = "Network")
}