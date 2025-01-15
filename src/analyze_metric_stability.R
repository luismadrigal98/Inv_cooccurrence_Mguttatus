analyze_metric_stability <- function(results) 
{
  #' Analyze the stability of metrics across different frequency combinations
  #' 
  #' This function takes a data frame of metric values calculated across different
  #' frequency combinations and assesses the stability of these metrics. Stability
  #' is assessed in terms of the variation of medians and the complexity of the
  #' relationship between metrics and frequencies.
  #' 
  #' @param results A data frame containing metric values calculated across different
  #' frequency combinations. The data frame should contain columns for the two
  #' frequencies (freq1, freq2) and the metric values (D, D_prime, cJaccard, Affinity).
  #' 
  #' @return A list containing the following elements:
  #' 
  #' - median_variation: A data frame containing the standard deviation, mean, and
  #' relative variation of the medians for each metric.
  #' 
  #' - edf_scores: A numeric vector containing the effective degrees of freedom for
  #' each metric's relationship with the product of frequencies.
  #' 
  #' - freq_correlations: A matrix containing the correlation coefficients between
  #' the frequencies and the metrics.
  #' 
  #' - freq_response_plot: A ggplot object showing the relationship between metric
  #' values and the product of frequencies.
  #' 
  #' - freq_ranges: A data frame containing the normalized range of metric values
  #' for each frequency combination.
  #' ___________________________________________________________________________
  
  # Calculate normalized range for each frequency combination
  freq_ranges <- results %>%
    group_by(freq1, freq2) %>%
    summarize(
      D_median = median(D),
      D_range = diff(range(D)),
      D_prime_median = median(D_prime),
      D_prime_range = diff(range(D_prime)),
      cJaccard_median = median(cJaccard),
      cJaccard_range = diff(range(cJaccard)),
      Affinity_median = median(Affinity),
      Affinity_range = diff(range(Affinity)),
      .groups = 'drop'
    )
  
  # Calculate variation of medians for each metric
  median_stats <- data.frame(
    metric = c("D", "D_prime", "cJaccard", "Affinity"),
    sd = c(
      sd(freq_ranges$D_median),
      sd(freq_ranges$D_prime_median),
      sd(freq_ranges$cJaccard_median),
      sd(freq_ranges$Affinity_median)
    ),
    mean = c(
      mean(abs(freq_ranges$D_median)),
      mean(abs(freq_ranges$D_prime_median)),
      mean(abs(freq_ranges$cJaccard_median)),
      mean(abs(freq_ranges$Affinity_median))
    )
  ) %>%
    mutate(relative_variation = sd/mean)
  
  # Create visualization of metric values vs frequencies
  freq_response_plot <- results %>%
    pivot_longer(
      cols = c(D, D_prime, cJaccard, Affinity),
      names_to = "Metric",
      values_to = "Value"
    ) %>%
    ggplot(aes(x = freq1 * freq2, y = Value, color = Metric)) +
    geom_point(alpha = 0.1) +
    geom_smooth(method = "loess", se = TRUE) +
    facet_wrap(~Metric, scales = "free_y") +
    theme_minimal() +
    labs(
      title = "Metric Values vs Product of Frequencies",
      x = "Frequency Product (freq1 * freq2)",
      y = "Metric Value"
    )
  
  # Calculate mean values for each frequency combination
  freq_sensitivity <- results %>%
    group_by(freq1, freq2) %>%
    summarize(
      D_mean = mean(D),
      D_prime_mean = mean(D_prime),
      cJaccard_mean = mean(cJaccard),
      Affinity_mean = mean(Affinity),
      .groups = 'drop'
    ) %>%
    mutate(freq_prod = freq1 * freq2)
  
  # Fit GAMs to assess complexity of frequency relationship
  gam_fits <- list(
    D = gam(D_mean ~ s(freq_prod), data = freq_sensitivity),
    D_prime = gam(D_prime_mean ~ s(freq_prod), data = freq_sensitivity),
    cJaccard = gam(cJaccard_mean ~ s(freq_prod), data = freq_sensitivity),
    Affinity = gam(Affinity_mean ~ s(freq_prod), data = freq_sensitivity)
  )
  
  # Extract effective degrees of freedom as measure of relationship complexity
  edf_scores <- sapply(gam_fits, function(x) summary(x)$edf)
  
  # Create correlation matrix between metrics and frequencies
  freq_correlations <- cor(
    results %>% 
      select(freq1, freq2, D, D_prime, cJaccard, Affinity)
  )[1:2, 3:6]
  
  # Return all results
  list(
    median_variation = median_stats,
    edf_scores = edf_scores,
    freq_correlations = freq_correlations,
    freq_response_plot = freq_response_plot,
    freq_ranges = freq_ranges
  )
}