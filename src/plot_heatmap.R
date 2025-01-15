plot_heatmap <- function(data, measure) 
{
  #' Plot a heatmap of the mean of a measure by two inversion frequencies
  #' 
  #' @param data A data frame with columns freq1, freq2, and measure
  #' 
  #' @param measure The name of the column in data to plot
  #' 
  #' @return A ggplot2 object
  #' ___________________________________________________________________________
  
  data %>%
    group_by(freq1, freq2) %>%
    summarize(
      mean_value = mean(!!sym(measure)),
      sd_value = sd(!!sym(measure))
    ) %>%
    ggplot(aes(x=freq1, y=freq2, fill=mean_value)) +
    geom_tile() +
    scale_fill_gradient2(
      low="blue", 
      mid="white",
      high="red",
      midpoint=0
    ) +
    labs(
      title=paste("Mean", measure, "by Inversion Frequencies"),
      x="Frequency of Inversion 1",
      y="Frequency of Inversion 2"
    ) +
    theme_bw() +
    theme(panel.grid = element_blank(),
          axis.text = element_text(color = 'black'))
}