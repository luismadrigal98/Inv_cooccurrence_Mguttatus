find_network_motifs <- function(networks, motifs, plot_individual_motifs = TRUE,
                                test_significance = TRUE, parallel = TRUE,
                                n_permutations = 1000) 
{
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
  
  # Plotting the individual motifs as reference
  if (plot_individual_motifs) {
    plot_dir <- file.path(getwd(), "Results", "Plots", "Individual_motifs", 
                          "Networks", "Motifs")
    
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
    combinations <- expand.grid(
      Cross = results$Cross,
      motif_names = results$motif_names,
      stringsAsFactors = FALSE
    )
    
    if (parallel) {

      # Export required objects to cluster
      parallel::clusterExport(cl, 
                              c("networks", "motifs", "results", 
                                "count_motifs_in_net", "permutation_test_motifs",
                                "rewire", "keeping_degseq", "gsize",
                                "graph.get.subisomorph"), 
                              envir = environment())
      
      # Run permutation tests
      significance_results <- parallel::parLapply(cl, 1:nrow(combinations), 
                                                  function(i) {
                                                    permutation_test_motifs(
                                                      results,
                                                      combinations$Cross[i], 
                                                      combinations$motif_names[i],
                                                      n_permutations)
                                                  })
    } else {
      # Run tests sequentially
      significance_results <- lapply(1:nrow(combinations), function(i) {
        permutation_test_motifs(combinations$Cross[i], 
                                combinations$motif_names[i],
                                n_permutations)
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