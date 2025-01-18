find_network_motifs <- function(networks, motifs, plot_individual_motifs = TRUE,
                                test_significance = TRUE, parallel = TRUE,
                                prevalence_thr = NULL,
                                n_permutations = 1000) 
{
  #' Find motifs in networks and test their significance
  #' 
  #' @param networks A list of igraph objects
  #' @param motifs A list of igraph objects
  #' @param plot_individual_motifs Logical indicating whether to plot the 
  #' individual motifs
  #' @param test_significance Logical indicating whether to test the significance 
  #' the motifs
  #' @param parallel Logical indicating whether to run the permutation tests 
  #' in parallel
  #' @param prevalence_thr Numeric indicating the minimum prevalence of a motif 
  #' to be considered
  #' @param n_permutations The number of permutations to perform
  #' ___________________________________________________________________________
  
  net_names <- names(networks)
  
  observed_counts <- lapply(networks, count_motifs_in_net, motifs)
  
  for (i in 1:length(observed_counts)) {
    observed_counts[[i]] <- observed_counts[[i]] |>
      dplyr::mutate(Cross = net_names[[i]])
  }
  
  results <- do.call('rbind', observed_counts)
  
  # Removing motifs that are missing in all crosses
  unused_motifs <- results |>
    dplyr::group_by(motif_names) |>
    dplyr::summarise(missing = all(counts == 0))
  
  motifs <- motifs[!names(motifs) %in% unused_motifs$motif_names[
    unused_motifs$missing]]
  
  # Filtering observed counts to remove the missing motifs
  results <- results |>
    dplyr::filter(motif_names %in% names(motifs))
  
  if (!is.null(prevalence_thr))
  {
    results <- results |>
      dplyr::filter(counts > prevalence_thr)
  }
  
  # Plotting the individual motifs as reference
  if (plot_individual_motifs) {
    plot_dir <- file.path(getwd(), "Results", "Plots", "Networks", "Motifs")
    
    if (!dir.exists(plot_dir)) {
      dir.create(plot_dir, recursive = TRUE)
    }
    
    lapply(names(motifs), function(x) {
      pdf(file = file.path(plot_dir, paste0(x, "_Motif.pdf")), 
          height = 4, width = 4)
      plot(motifs[[x]])
      dev.off()
    })
  }
  
  # Hypothesis testing framework
  if (test_significance) {
    # Create combinations of networks and motifs to test
    combinations <- results[,c('motif_names', 'Cross')]
    
    if (parallel) {

      # Export required objects to cluster
      parallel::clusterExport(cl, 
                              c("networks", "motifs", "results", 
                                "count_motifs_in_net", "permutation_test_motifs",
                                "rewire", "keeping_degseq", "gsize",
                                "mclapply", "degree"), 
                              envir = environment())
      
      # Run permutation tests
      significance_results <- parallel::parLapply(cl, 1:nrow(combinations), 
                                                  function(i) {
                                                    permutation_test_motifs(
                                                      results,
                                                      networks,
                                                      motifs,
                                                      combinations$Cross[i], 
                                                      combinations$motif_names[i],
                                                      n_permutations,
                                                      parallel)
                                                  })
    } else {
      # Run tests sequentially
      significance_results <- lapply(1:nrow(combinations), function(i) {
        permutation_test_motifs(results,
                                networks,
                                motifs,
                                combinations$Cross[i], 
                                combinations$motif_names[i],
                                n_permutations,
                                parallel)
      })
    }
    
    # Combine results
    significance_results <- do.call(rbind, significance_results)
    
    # Add multiple testing correction
    significance_results$p_adjusted <- stats::p.adjust(
      significance_results$p_value, 
      method = "BH"
    )
    
    # Add significance results to output
    results <- merge(results, significance_results, 
                     by = c("Cross", "motif_names"))
  }
  
  return(results)
}