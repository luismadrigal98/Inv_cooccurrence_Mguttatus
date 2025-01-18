find_network_motifs <- function(networks, motifs, plot_individual_motifs = TRUE,
                                test_significance = TRUE, parallel = TRUE,
                                prevalence_thr = NULL,
                                n_permutations = 1000) 
{
  #' Find motifs in networks and test their significance
  #' 
  #' @param networks A list of igraph objects
  #' @param motifs A list of igraph objects
  #' @param plot_individual_motifs Logical indicating whether to plot motifs
  #' @param test_significance Logical indicating whether to test significance
  #' @param parallel Logical indicating whether to run in parallel
  #' @param prevalence_thr Numeric threshold for motif counts
  #' @param n_permutations Number of permutations for testing
  
  net_names <- names(networks)
  
  # Get observed counts
  observed_counts <- lapply(seq_along(networks), function(i) {
    counts <- count_motifs_in_net(networks[[i]], motifs)
    counts$Cross <- net_names[i]
    counts
  })
  
  results <- do.call('rbind', observed_counts)
  
  # Remove unused motifs
  unused_motifs <- results %>%
    dplyr::group_by(motif_names) %>%
    dplyr::summarise(missing = all(counts == 0))
  
  motifs <- motifs[!names(motifs) %in% unused_motifs$motif_names[unused_motifs$missing]]
  results <- results %>% dplyr::filter(motif_names %in% names(motifs))
  
  # Apply prevalence threshold if specified
  if (!is.null(prevalence_thr)) {
    results <- results %>% dplyr::filter(counts > prevalence_thr)
  }
  
  # Plot motifs if requested
  if (plot_individual_motifs) {
    plot_dir <- file.path(getwd(), "Results", "Plots", "Networks", "Motifs")
    dir.create(plot_dir, recursive = TRUE, showWarnings = FALSE)
    
    lapply(names(motifs), function(x) {
      pdf(file.path(plot_dir, paste0(x, "_Motif.pdf")), height = 4, width = 4)
      plot(motifs[[x]])
      dev.off()
    })
  }
  
  # Perform significance testing
  if (test_significance && nrow(results) > 0) {
    # Get unique combinations to test
    combinations <- unique(results[, c("Cross", "motif_names")])
    
    if (parallel) {
      # Use mclapply for parallel processing
      significance_results <- parallel::mclapply(
        1:nrow(combinations),
        function(i) {
          permutation_test_motifs(
            results,
            networks,
            motifs,
            combinations$Cross[i],
            combinations$motif_names[i],
            n_permutations,
            parallel = FALSE  # Avoid nested parallelization
          )
        },
        mc.cores = parallel::detectCores() - 1
      )
    } else {
      significance_results <- lapply(1:nrow(combinations), function(i) {
        permutation_test_motifs(
          results,
          networks,
          motifs,
          combinations$Cross[i],
          combinations$motif_names[i],
          n_permutations,
          parallel = FALSE
        )
      })
    }
    
    # Combine results and adjust p-values
    significance_results <- do.call(rbind, significance_results)
    significance_results$p_adjusted <- stats::p.adjust(
      significance_results$p_value,
      method = "BH"
    )
    
    # Merge with original results
    results <- merge(results, significance_results,
                     by = c("Cross", "motif_names"))
  }
  
  return(results)
}