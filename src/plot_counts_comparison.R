plot_counts_comparison <- function(data) 
{
  #' This function plots the observed vs expected counts of motifs.
  #' 
  #' @param data A data frame with the following columns:
  #'  - motif_names: The names of the motifs.
  #'  - observed: The observed counts of the motifs.
  #'  - mean_null: The expected counts of the motifs.
  #'  - p_adjusted: The adjusted p-values of the motifs.
  #'  - z_score: The z-scores of the motifs.
  #'  
  #'  @returns ggplot2 object
  #' ___________________________________________________________________________
  
  data %>%
    ggplot(aes(x = mean_null, y = observed, group = motif_names)) +
    geom_point(aes(color = p_adjusted < 0.05, size = abs(z_score))) +
    geom_abline(linetype = "dashed", color = "grey50") +
    scale_color_manual(values = c("grey70", "dodgerblue"),
                       name = "Significant") +
    theme_bw() +
    theme(axis.text = element_text(color = "black"),
          panel.grid = element_blank()) +
    labs(title = "Observed vs Expected Motif Counts",
         x = "Expected Count", y = "Observed Count")
}