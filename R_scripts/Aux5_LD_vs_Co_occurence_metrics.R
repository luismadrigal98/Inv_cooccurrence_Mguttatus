##' >>>>>>>>>>>>>>>>>>>>>>>> AUXILIARY SCRIPT 5 <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
##' 
##' @title Are you with me? Co-occurrence tests from community ecology can 
##' identify positive and negative epistasis between inversions in Mimulus 
##' guttatus.
##' 
##' @description This script is designed to compare the LD-related metrics D and
##' D' with the co-occurrance metrics of this paper (centered Jaccard-Tanimoto
##' and affinity)
##' 
##' @author Luis Javier Madrigal-Roca & John K. Kelly
##' 
##' @date 2025-01-14
##' ____________________________________________________________________________

## 1) Defining the frequencies to test in the simulation ----

freqs <- seq(0.1, 0.9, 0.1)

## 2) Getting the results for the different frequencies and replicates ----

results <- compare_measures(freqs, n_reps = 1000, association_strength = 0.5)

## 3) Plotting the mean value of the statistics as a function of the allele freq ----

p0 <- plot_heatmap(results, "D")
p1 <- plot_heatmap(results, "D_prime")
p2 <- plot_heatmap(results, "cJaccard")
p3 <- plot_heatmap(results, "Affinity")

## 3.1) Getting all the plots arranged in a 2x2 grid ----

p <- grid.arrange(p0, p1, p2, p3, ncol = 2)

ggsave(file = "Results/Plots/Mean_metrics_per_freq_combination.pdf", plot = p,
       width = 10, height = 8)

## 4) Analyzing the relationship of LD and co-occurrence metrics with allele freq ----

stability_analysis <- analyze_metric_stability(results)

write.table(x = stability_analysis$median_variation, 
            file = "Results/median_variation_stability_analysis.tsv",
            sep = '\t')

write.table(stability_analysis$edf_scores, 
            file = "Results/edf_scores_metrics.tsv",
            sep = '\t')

write.table(stability_analysis$freq_correlations, 
            file = "Results/edf_scores_metrics.tsv",
            sep = '\t')