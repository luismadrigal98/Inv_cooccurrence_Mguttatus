plot_zscores <- function(data) 
{
  #' This function plots the z-scores of motifs across networks.
  #' 
  #' @param data A data frame with columns Cross, motif_names, z_score, 
  #' and p_adjusted.
  #' 
  #' @return A ggplot object.
  #' ___________________________________________________________________________
  
  ggplot(data, aes(x = reorder(paste(Cross, motif_names), z_score), 
                   y = z_score, fill = p_adjusted < 0.05)) +
    geom_bar(stat = "identity") +
    scale_fill_manual(values = c("grey70", "dodgerblue"),
                      name = "Significant") +
    theme_bw() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1),
          axis.text = element_text(color = "black"),
          panel.grid = element_blank()) +
    labs(title = "Motif Z-scores Across Networks",
         x = "Network-Motif Pair", y = "Z-score")
}