plot_enrichment <- function(data) 
{
  #' Plot motif enrichment by cross
  #' 
  #' @param data A data frame with columns motif_names, Cross, observed, 
  #' mean_null, p_adjusted 
  #' 
  #' @return A ggplot object
  #' ___________________________________________________________________________
  
  data %>%
    mutate(enrichment = observed / mean_null,
           significant = p_adjusted < 0.05,
           motif_cross = interaction(motif_names, Cross, sep = " - ")) %>%
    ggplot(aes(x = reorder(motif_cross, enrichment), 
               y = enrichment, fill = significant)) +
    geom_bar(stat = "identity", position = "dodge") +
    scale_fill_manual(values = c("grey70", "dodgerblue")) +
    theme_bw() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1),
          axis.text = element_text(color = "black"),
          panel.grid = element_blank()) +
    labs(title = "Motif Enrichment by Cross",
         x = "Motif - Cross", y = "Observed/Expected Ratio")
}