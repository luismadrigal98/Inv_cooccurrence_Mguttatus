## Co_occurrence of inversions in lines of Mimulus guttatus
#@ Main Script
#@ Authors: Luis J. Madrigal-Roca & John K. Kelly
#@ Date: 2024-05-27

## *****************************************************************************
## 3) Chi-square test for independence of inversions ----
## _____________________________________________________________________________

expec_freq <- matrix(data = c(1/16, 2/16, 1/16, 2/16, 4/16, 2/16,
                              1/16, 2/16, 1/16), nrow = 3, ncol = 3, 
                     dimnames = list(c('0', '1', '2'),
                                     c('0', '1', '2')))

# 3.2) X2 square test ----

results_x2 <- lapply(X = Data_d_l, FUN = function(observed_matrix)
{
  mapply(x2_calculator,
         rep(1:nrow(observed_matrix), 
             each = nrow(observed_matrix)), 
         rep(1:nrow(observed_matrix), 
             nrow(observed_matrix)),
         MoreArgs = list(dosage_matrix = observed_matrix),
         SIMPLIFY = F)
})

# Condensing the results into a dataframe per line

results_x2_df <- lapply(results_x2, function(x)
{
  do.call(rbind, x)
})

# Exporting these results to perform analysis in Aux script 4

saveRDS(results_x2_df, file = "Results/results_x2_df.rds")

# Extracting the p_values as an square matrix ---

x2_p_square_matrix <- lapply(results_x2_df, x2p_to_square)

# Filtering the results

results_x2_filtered_df <- sapply(X = results_x2_df,
                                 FUN = contingency_filter,
                                 metadata,
                                 simplify = F)

# Correction of the p-values ----

results_x2_filtered_df <- sapply(X = results_x2_filtered_df, 
                                 FUN = p_corrector, var = "p", 
                        simplify = F)

# Getting the number of combinations without dosage
combos_without_d <- unique(do.call(
  'rbind', results_x2_filtered_df)$INV_combination)

## Additional analysis ----

# Expanding the dataframe to decompose each test into individual components
# This will preserve the omnibus p-value of the X2 test for the four contrasts
# derived from the test

x2_results_filtered_expanded <- lapply(results_x2_filtered_df, 
                                       contingency_expanded)

for(i in seq_along(x2_results_filtered_expanded))
{
  x2_results_filtered_expanded[[i]] <- x2_results_filtered_expanded[[i]] |>
    mutate(cell_location = as.factor(cell_location))
  
  levels(x2_results_filtered_expanded[[i]]$cell_location) <- 
    c("_1-_1", "_1-_2", "_2-_1", "_2-_2")
  
  x2_results_filtered_expanded[[i]] <- x2_results_filtered_expanded[[i]] |>
    mutate(cell_location = as.character(cell_location)) |>
    mutate(INV_1 = paste0(INV_1, sapply(strsplit(cell_location, '-'), `[`, 1)),
           INV_2 = paste0(INV_2, sapply(strsplit(cell_location, '-'), `[`, 2)))
}

for (i in 1:length(names(x2_results_filtered_expanded))) {
  x2_results_filtered_expanded[[i]] <- x2_results_filtered_expanded[[i]] |>
    rowwise() |>
    mutate(INV_combination = paste0(sort(c(INV_1, INV_2)), collapse = "")) |>
    ungroup()
}

for (i in 1:length(names(x2_results_filtered_expanded))) {
  x2_results_filtered_expanded[[i]] <- x2_results_filtered_expanded[[i]] |>
    dplyr::rename(X2_global = X2, p_X2_global = p, p_X2_corrected = p_corrected)
  
  rm(i)
}

# 3.1) Post-hoc analysis for the X2 test ----

results_posthoc_x2 <- lapply(X = Data_d_l, FUN = function(observed_matrix)
{
  mapply(x2_posthoc_calculator,
         rep(1:nrow(observed_matrix), 
             each = nrow(observed_matrix)), 
         rep(1:nrow(observed_matrix), 
             nrow(observed_matrix)),
         MoreArgs = list(dosage_matrix = observed_matrix,
                         expectation = expec_freq),
         SIMPLIFY = F)
})

results_posthoc_x2_df <- lapply(results_posthoc_x2, function(x)
{
  do.call(rbind, x)
})

for (i in names(results_posthoc_x2_df))
{
  results_posthoc_x2_df[[i]] <- results_posthoc_x2_df[[i]] |>
    mutate(INV2 = as.character(INV2)) |>
    rowwise() |>
    mutate(INV_combination = paste0(sort(c(INV1, INV2)), collapse = '')) |>
    filter(! any(grepl("_0", INV1) | grepl("_0", INV2))) |>
    ungroup()
  
  rm(i)
}

## Extracting the p_values as an square matrix ----

x2_p_post_hoc_square_matrix <- lapply(results_posthoc_x2_df, x2p_to_square,
                                      c("INV1", "INV2"), 'p_value')

## Filtering out irrelevant contrasts ----

results_posthoc_x2_df_filtered <- sapply(X = results_posthoc_x2_df, 
                                         FUN = function(x)
{
  x |>
    mutate(ID1 = sapply(strsplit(INV1, "_"), `[`, 2),
           ID2 = sapply(strsplit(INV2, "_"), `[`, 2)) |>
    mutate(chrom1 = metadata[match(ID1, metadata$INV_ID), "Chr"],
           chrom2 = metadata[match(ID2, metadata$INV_ID), "Chr"]) |>
    filter(chrom1 != chrom2)
}, simplify = F)

for (i in names(results_posthoc_x2_df_filtered))
{
  results_posthoc_x2_df_filtered[[i]] <- 
    results_posthoc_x2_df_filtered[[i]][!duplicated(
      results_posthoc_x2_df_filtered[[i]]$INV_combination), ]
  
  rm(i)
}

for (i in names(results_posthoc_x2_df_filtered))
{
  results_posthoc_x2_df_filtered[[i]] <- results_posthoc_x2_df_filtered[[i]] |>
    mutate(p_adjusted = p.adjust(p_value, method = "BY"))
  
  rm(i)
}

## Condensing all the results (main and post hoc) into a single dataframe

for (i in 1: length(names(x2_results_filtered_expanded)))
{
  x2_results_filtered_expanded[[i]] <- merge(x2_results_filtered_expanded[[i]], 
                                             results_posthoc_x2_df_filtered[[i]], 
                                             by = "INV_combination", 
                                             suffixes = c("", "_posthoc"))
  
  rm(i)
}

for (i in 1: length(names(x2_results_filtered_expanded)))
{
  x2_results_filtered_expanded[[i]] <- x2_results_filtered_expanded[[i]] |>
    dplyr::select(!c("ID1", "ID2", "cell_location", 
              "standardized_residual", "INV1", "INV2", "ID1_posthoc",
              "ID2_posthoc", "chrom1_posthoc", "chrom2_posthoc")) |>
    dplyr::rename(p_X2_posthoc = p_value,
           p_X2_posthoc_corrected = p_adjusted)
}

## Getting the number of combinations with dosage
combos_with_d <- unique(do.call(
  'rbind', x2_results_filtered_expanded)$INV_combination)

## *****************************************************************************
## 4) Calculation of alternative indexes used in community studies ----
## _____________________________________________________________________________

# Filtering out rows with 0 marginals. If something is not there in any instance,
# we cannot test for association with another inversion.

Data_p_a <- lapply(Data_p_a, function(x) x[!apply(x, 1, 
                                                  function(y) all(y == 0)), ])

# Instance INV_49_2 in line 444 was eliminated from the analysis

# 4.1) Jaccard ----

jaccard_tanimoto_results <- sapply(X = Data_p_a, 
                                   FUN = jaccard.test.pairwise,
                                   method = 'bootstrap',
                                   compute.qvalue = F,
                                   simplify = F,
                                   B = 10000)

## Correcting the dim names of the outputs

for (i in 1:length(Data_p_a))
{
  denominations <- attributes(Data_p_a[[i]])[[2]][[1]]
  
  for (j in names(jaccard_tanimoto_results[[i]]))
  {
    dimnames(jaccard_tanimoto_results[[i]][[j]]) <- list(denominations,
                                                         denominations)
  }
  
  rm(i)
  rm(j)
}

# Extracting the p_values as a square matrix ----

jaccard_tanimoto_p_square_matrix <- list()

for (i in 1:length(jaccard_tanimoto_results)) {
  # Ensure that pvalues is a matrix
  mat <- as.matrix(jaccard_tanimoto_results[[i]]$pvalues)
  
  # Make the matrix symmetric
  mat[lower.tri(mat)] <- t(mat)[lower.tri(mat)]
  
  # Assign the symmetric matrix back to the list
  jaccard_tanimoto_p_square_matrix[[i]] <- mat
}

## Extracting the results as a df ---

jaccard_tanimoto_results_df <- sapply(X = jaccard_tanimoto_results,
                                      FUN = jaccard_tanimoto_condenser,
                                      metadata,
                                      simplify = F)

# Filtering out elements in the same chromosome

jaccard_tanimoto_results_df_filtered <- jaccard_tanimoto_results_df

for (i in 1:length(jaccard_tanimoto_results_df_filtered))
{
  jaccard_tanimoto_results_df_filtered[[i]] <- 
    jaccard_tanimoto_results_df_filtered[[i]] |>
    dplyr::filter(chrom_1 != chrom_2)
  
  rm(i)
}

# Correcting the p_values ----
jaccard_tanimoto_results_df_filtered <- sapply(X = jaccard_tanimoto_results_df_filtered, 
                                      FUN = p_corrector, var = "p_value", 
                                      method = "BY",
                                      simplify = F)

# 4.2) Affinity measure ----

results_affinity <- lapply(X = Data_p_a, FUN = function(x)
{
  x <- affinity(x, row.or.col = 'row', 
  datatype = 'binary', sigdigit = 3)
  x
})

for (i in 1:length(results_affinity))
{
  results_affinity[[i]]$all <- results_affinity[[i]]$all |>
    mutate(p_value = as.numeric(p_value))
  
  rm(i)
}

# Extracting the affinity values from the output as a matrix for plotting

A_matrices <- lapply(results_affinity, alpha_values_extractor)

for (i in 1:length(names(A_matrices))) {
  old_name <- names(A_matrices)[i]
  new_name <- paste0("A_", old_name)
  names(A_matrices)[i] <- new_name
  rm(i)
  rm(new_name)
  rm(old_name)
}

results_affinity_filtered <- sapply(X = results_affinity,
                                    FUN = affinity_filter,
                                    metadata,
                                    simplify = F)

# Correcting the p_values ----

results_affinity_filtered <- sapply(X = results_affinity_filtered, 
                                    FUN = p_corrector, 
                                    var = "p_value", 
                                    method = "BY",
                                    simplify = F)

# 5) Extracting the significant values detected for the three scores ----

# Common field for all datasets: INV_combination

for (i in 1:length(names(results_affinity_filtered))) {
  results_affinity_filtered[[i]] <- results_affinity_filtered[[i]] |>
    rowwise() |>
    mutate(INV_combination = paste0(sort(c(entity_1, entity_2)), collapse = "")) |>
    ungroup()
}

for (i in 1:length(names(jaccard_tanimoto_results_df_filtered))) {
  jaccard_tanimoto_results_df_filtered[[i]] <- 
    jaccard_tanimoto_results_df_filtered[[i]] |>
    rowwise() |>
    mutate(INV_combination = paste0(sort(c(entity_1, entity_2)), 
                                    collapse = "")) |>
    ungroup()
  
  rm(i)
}

## Merging all dataframes into a single one per line----

merged_results <- list()

for (i in 1:9) {
  merged_results[[i]] <- merge(x2_results_filtered_expanded[[i]], 
                               results_affinity_filtered[[i]], 
                               by = "INV_combination", 
                               suffixes = c("_x2", "_affinity"))
  
  merged_results[[i]] <- merge(merged_results[[i]], 
                               jaccard_tanimoto_results_df_filtered[[i]], 
                               by = "INV_combination", 
                               suffixes = c("", "_jaccard_tanimoto"))
}

names(merged_results) <- names(x2_results_filtered_expanded)

## Combining all data frames in an unique one ----

# Add a new column to each data frame with the line information
for (i in 1:length(merged_results)) {
  merged_results[[i]]$line <- names(merged_results)[i]
}

# Combine all data frames into one
final_result <- do.call(rbind, merged_results)

# Depuration of the final data frame

final_result <- final_result |>
  dplyr::select(line, INV_1, INV_2, chrom_1, chrom_2, X2_global, p_X2_global, 
         p_X2_posthoc, p_X2_posthoc_corrected, relative_contribution, 
         Residual, alpha_mle, p_value, 
         p_corrected, jaccard_jaccard_tanimoto, p_value_jaccard_tanimoto,
         p_corrected_jaccard_tanimoto) |>
  rename(Line = line,
         Chr_1 = chrom_1, Chr_2 = chrom_2,
         X2_SR = Residual, X2_RC = relative_contribution,
         A_p = p_value, A_alpha = alpha_mle, 
         A_p_corrected = p_corrected,
         Jaccard = jaccard_jaccard_tanimoto, J_p = p_value_jaccard_tanimoto, 
         J_p_corrected = p_corrected_jaccard_tanimoto)

# 5.1) Relaxed results ----

final_result_relaxed_any <- final_result |> 
  filter(p_X2_posthoc < 0.05 | A_p < 0.05 | J_p < 0.05)

final_result_relaxed_all <- final_result |>
  filter(p_X2_global < 0.05 & p_X2_posthoc < 0.05 & A_p < 0.05 & J_p < 0.05)

# 5.2) Stringent results ----

final_result_stringent_any <- final_result |> 
  filter(p_X2_posthoc_corrected < 0.05 | A_p_corrected < 0.05 | 
           J_p_corrected < 0.05)

final_result_stringent_all <- final_result |>
  filter(p_X2_global < 0.05 & p_X2_posthoc_corrected < 0.05 & 
           A_p_corrected < 0.05 & J_p_corrected < 0.05)

## Saving the data frame sinto csv files

write.csv(final_result_relaxed_any, "Results/final_result_relaxed_any.csv", 
          row.names = FALSE)
write.csv(final_result_relaxed_all, "Results/final_result_relaxed_all.csv", 
          row.names = FALSE)
write.csv(final_result_stringent_any, "Results/final_result_stringent_any.csv", 
          row.names = FALSE)
write.csv(final_result_stringent_all, "Results/final_result_stringent_all.csv", 
          row.names = FALSE)

# 5.3) Summarizing overall patterns ----

deviants_summary <- deviant_signal_summarizer(final_result_relaxed_all)

write.table(deviants_summary, "Results/deviants_summary.txt", 
            row.names = FALSE, quote = FALSE, sep = "\t")

# 6) Visualization of the relevant patterns in terms of evidence ----

# 6.1) Evidence according to X2_global ----

for (i in names(x2_p_square_matrix)) 
{
  pdf(file = paste0("Results/Plots/", i, "_omnibus_X2_support_heatmap.pdf"), width = 6, 
      height = 6)
  
  logical_matrix <- x2_p_square_matrix[[i]] < 0.05
  info <- matrix(as.numeric(logical_matrix), nrow = nrow(logical_matrix), 
                 ncol = ncol(logical_matrix))
  rownames(info) <- colnames(info) <- rownames(logical_matrix)
  
  diag(info) <- NA
  
  heatmap.2(x = info, dendrogram = 'none', breaks = 3, 
            symbreaks = F,
            col = c("white", "#64AC59"),
            key = T, density.info = 'histogram', key.title = "Support",
            main = paste0("Pairwise", "_", i), xlab = "Inversions",
            ylab = "Inversions", trace = "none", denscol = "black",
            densadj = 0.50)
  
  dev.off()
  
  rm(i)
  rm(info)
}

# 6.2) Evidence according to X2_posthoc, Jaccard, and Affinity ----

support_matrix <- mapply(support_counter, 
                         supp_X2_posthoc = x2_p_post_hoc_square_matrix,
                         supp_Affinity = A_matrices,
                         supp_Jaccard = jaccard_tanimoto_p_square_matrix,
                         SIMPLIFY = F)

for (i in names(support_matrix)) 
{
  pdf(file = paste0("Results/Plots/", i, "_concensus_support_heatmap.pdf"), width = 7, 
      height = 7)
  
  count_matrix <- support_matrix[[i]]
  
  heatmap.2(x = count_matrix, dendrogram = 'none', breaks = 5, 
            symbreaks = F,
            col = c("white", "#EB4F48", "#F4D166", "#64AC59"),
            key = T, density.info = 'histogram', key.title = "Support",
            main = paste0("Pairwise", "_", i), xlab = "Inversions",
            ylab = "Inversions", trace = "none", denscol = "black",
            densadj = 0.50)
  
  dev.off()
  
  rm(i)
}

## 6.3) Detecting Linkage Groups ----

# Masking the p_values of contrast in different chromosomes

x2_p_square_matrix_masked <- lapply(x2_p_square_matrix, 
                                    mask_chromosomes, metadata)

jaccard_tanimoto_p_square_matrix_masked <- lapply(jaccard_tanimoto_p_square_matrix, 
                                                 mask_chromosomes, metadata)

A_matrices_masked <- lapply(A_matrices, mask_chromosomes, metadata)

x2_p_post_hoc_square_matrix_masked <- lapply(x2_p_post_hoc_square_matrix, 
                                      mask_chromosomes, metadata)

# 6.3.1) Heatmaps for the p_values of the different metrics ----

# The inversion in the heatmap will appear sorted by chromosomes

# 6.3.1.1) Linkage groups according to omnibus X2 test ----

for (i in names(x2_p_square_matrix_masked)) 
{
  pdf(file = paste0("Results/Plots/", i, "_omnibus_X2_linkage_heatmap.pdf"), width = 6, 
      height = 6)
  
  logical_matrix <- x2_p_square_matrix_masked[[i]] < 0.05
  info <- matrix(as.numeric(logical_matrix), nrow = nrow(logical_matrix), 
                 ncol = ncol(logical_matrix))
  
  rownames(info) <- colnames(info) <- rownames(logical_matrix)
  
  sorted_rownames <- mixedsort(rownames(info))
  
  info <- info[sorted_rownames, sorted_rownames]
  
  diag(info) <- NA
  
  info[is.na(info)] <- F
  
  heatmap.2(x = info, dendrogram = 'none', breaks = 3, 
            Rowv = F,
            Colv = F,
            symbreaks = F,
            col = c("white", "blue"),
            key = F, density.info = 'none', key.title = "",
            main = paste0("Pairwise", "_", i), xlab = "Inversions",
            ylab = "Inversions", trace = "none",
            )
  
  dev.off()
  
  rm(i)
  rm(info)
}

# 6.3.1.1) Linkage groups according to the 2x2 based tests ----

support_matrix_LG <- mapply(support_counter, 
                         supp_X2_posthoc = x2_p_post_hoc_square_matrix_masked,
                         supp_Affinity = A_matrices_masked,
                         supp_Jaccard = jaccard_tanimoto_p_square_matrix_masked,
                         SIMPLIFY = F)

# Impossible contrast cleaner (same inversion in different dosage) ----

for (i in names(support_matrix_LG)) {
  for (j in 1:nrow(support_matrix_LG[[i]])) {
    for (k in 1:ncol(support_matrix_LG[[i]])) {
      # Splitting row names and comparing the specific part
      if (strsplit(rownames(support_matrix_LG[[i]])[j], split = "_")[[1]][2] ==
          strsplit(rownames(support_matrix_LG[[i]])[k], split = "_")[[1]][2]) {
        support_matrix_LG[[i]][j, k] <- NA
      }
    }
  }
  
  rm(i)
  rm(j)
  rm(k)
}

for (i in names(support_matrix_LG)) 
{
  pdf(file = paste0("Results/Plots/", i, "_concensus_linkage_heatmap.pdf"), width = 7, 
      height = 7)
  
  count_matrix <- support_matrix_LG[[i]]
  
  sorted_rownames <- mixedsort(rownames(count_matrix))
  
  count_matrix <- count_matrix[sorted_rownames, sorted_rownames]
  
  heatmap.2(x = count_matrix, dendrogram = 'none', breaks = 5, 
            symbreaks = F,
            Rowv = F,
            Colv = F,
            col = c("white", "#EB4F48", "#F4D166", "#64AC59"),
            key = T, density.info = 'histogram', key.title = "Support",
            main = paste0("Pairwise", "_", i), xlab = "Inversions",
            ylab = "Inversions", trace = "none", denscol = "black",
            densadj = 0.50)
  
  dev.off()
  
  rm(i)
}

# 7) Relation between statistics employed in terms of linear regression ----
# This will be done taking the overall metric across all crosses
# Linear models:

# 7.1) X2_SR vs Jaccard ----

model1 <- lm(X2_SR ~ Jaccard, data = final_result)

# 7.2) X2_SR vs Affinity ----

model2 <- lm(X2_SR ~ A_alpha, data = final_result)

# 7.3) Jaccard vs Affinity ----

model3 <- lm(Jaccard ~ A_alpha, data = final_result)

# 7.5) Plotting the results ----

# 7.5.1) Model 1 ----

pdf(file = "Results/Plots/Linear_models_X2_SR_Jaccard.pdf", width = 8, height = 6)

ggplot(data = final_result, aes(x = Jaccard, y = X2_SR)) +
  geom_point() +
  geom_smooth(method = "lm", se = FALSE, color = "red") +
  theme_bw() +
  annotate(geom = 'text', x = -0.0475, y = 3.8, 
           label = "X2_SR ~ 27.53 * Jaccard - 0.03, df = 4122, p < 2e-16, R^2 = 0.98") +
  theme(axis.text = element_text(size = 12, color = 'black'),
        axis.title = element_text(size = 14, color = 'black'),
        legend.text = element_text(size = 12, color = 'black'),
        legend.title = element_text(size = 14, color = 'black'),
        panel.grid = element_blank())

dev.off()

# 7.5.2) Model 2 ----

pdf(file = "Results/Plots/Linear_models_X2_SR_Affinity.pdf", width = 8, height = 6)

ggplot(data = final_result, aes(x = A_alpha, y = X2_SR)) +
  geom_point() +
  geom_smooth(method = "lm", se = FALSE, color = "red") +
  theme_bw() +
  annotate(geom = 'text', x = -4.25, y = 7, 
           label = "X2_SR ~ 0.79 * Alpha + 0.03, df = 4122, p < 2e-16, R^2 = 0.35") +
  theme(axis.text = element_text(size = 12, color = 'black'),
        axis.title = element_text(size = 14, color = 'black'),
        legend.text = element_text(size = 12, color = 'black'),
        legend.title = element_text(size = 14, color = 'black'),
        panel.grid = element_blank())

dev.off()

# 7.5.3) Model 3 ----

pdf(file = "Results/Plots/Linear_models_Jaccard_Affinity.pdf", width = 8, height = 6)

ggplot(data = final_result, aes(x = A_alpha, y = Jaccard)) +
  geom_point() +
  geom_smooth(method = "lm", se = FALSE, color = "red") +
  theme_bw() +
  annotate(geom = 'text', x = 0, y = 1, 
           label = "Jaccard ~ 0.03 * Alpha + 0.002, df = 4122, p < 2e-16, R^2 = 0.32") +
  theme(axis.text = element_text(size = 12, color = 'black'),
        axis.title = element_text(size = 14, color = 'black'),
        legend.text = element_text(size = 12, color = 'black'),
        legend.title = element_text(size = 14, color = 'black'),
        panel.grid = element_blank()) +
  ylim(c(-1, 1))

dev.off()

## After checking the premise of normality of residuals, and given the detected
# deviations, a robust linear model will be employed.

lillie.test(model1$residuals)
lillie.test(model2$residuals)
lillie.test(model3$residuals)

# 7.6) Robust linear models ----

# 7.6.1) X2_SR vs Jaccard ----
# Robust model
robust_model1 <- rlm(X2_SR ~ Jaccard, data = final_result)

pdf(file = "Results/Plots/Robust_Linear_models_X2_SR_Jaccard.pdf", width = 8, 
    height = 6)

ggplot(data = final_result[!robust_model1$w < 0.15,], 
       aes(x = Jaccard, y = X2_SR)) +
  geom_point() +
  geom_smooth(method = "rlm", se = FALSE, color = "red") +
  theme_bw() +
  annotate(geom = 'text', x = -0.0475, y = 3.8, 
           label = "X2_SR ~ 28.09 * Jaccard - 0.02, df = 4122, p < 0.001") +
  theme(axis.text = element_text(size = 12, color = 'black'),
        axis.title = element_text(size = 14, color = 'black'),
        legend.text = element_text(size = 12, color = 'black'),
        legend.title = element_text(size = 14, color = 'black'),
        panel.grid = element_blank())

dev.off()

# 7.6.2) X2_SR vs Affinity ----
# Robust model
robust_model2 <- rlm(X2_SR ~ A_alpha, data = final_result)

pdf(file = "Results/Plots/Robust_Linear_models_X2_SR_Affinity.pdf", width = 8, 
    height = 6)

ggplot(data = final_result[!robust_model2$w < 0.15,], 
       aes(x = A_alpha, y = X2_SR)) +
  geom_point() +
  geom_smooth(method = "rlm", se = FALSE, color = "red") +
  theme_bw() +
  annotate(geom = 'text', x = -0.0475, y = 3.8, 
           label = "X2_SR ~ 2.45 * Affinity + 0.01, df = 4122, p < 0.01") +
  theme(axis.text = element_text(size = 12, color = 'black'),
        axis.title = element_text(size = 14, color = 'black'),
        legend.text = element_text(size = 12, color = 'black'),
        legend.title = element_text(size = 14, color = 'black'),
        panel.grid = element_blank())

dev.off()

# 7.6.3) Jaccard vs Affinity ----
# Robust model
robust_model3 <- rlm(Jaccard ~ A_alpha, data = final_result)

pdf(file = "Results/Plots/Robust_Linear_models_Jaccard_Affinity.pdf", width = 8, 
    height = 6)

ggplot(data = final_result[!robust_model3$w < 0.15,], 
       aes(x = A_alpha, y = Jaccard)) +
  geom_point() +
  geom_smooth(method = "rlm", se = FALSE, color = "red") +
  theme_bw() +
  annotate(geom = 'text', x = -0.0475, y = 0.8, 
           label = "Jaccard ~ 0.09 * Affinity + 0.001, df = 4122, p < 0.001") +
  theme(axis.text = element_text(size = 12, color = 'black'),
        axis.title = element_text(size = 14, color = 'black'),
        legend.text = element_text(size = 12, color = 'black'),
        legend.title = element_text(size = 14, color = 'black'),
        panel.grid = element_blank()) +
  ylim(c(-1, 1))

dev.off()

# p-values for the robust models

rob.pvals(robust_model1)
rob.pvals(robust_model2)
rob.pvals(robust_model3)

# 8) Network analysis ----

# Splitting the final result data frame into a list of data frames per line

final_result_split <- split(final_result, final_result$Line)

# 8.1) Network builder function ----

nodes <- paste0("Inv_", metadata$INV_ID)

nodes <- sapply(nodes, function(x) {
  c(paste0(x, "_1"), paste0(x, "_2"))
})

nodes <- as.vector(unlist(nodes))

# Apply the function to your data
networks <- mapply(sig_network_builder, 
                   final_result_split,
                   MoreArgs = list(nodes = nodes,
                                   meta_nodes = metadata),
                   SIMPLIFY = FALSE)

# 8.2) Plotting the networks ----

plot_network(networks)

# Legend for the colors of chromosomes:

# Create a data frame with the chromosomes and their corresponding colors
df <- data.frame(chromosome = unique(V(networks[[1]])$chromosome), 
                 color = hcl.colors(14, "Set3", rev = TRUE))

# Create a plot that only contains the legend
pdf("Results/Plots/Networks/Chromosome_legend.pdf", width = 2, height = 8)

legend_plot <- ggplot(df, aes(x = 1, fill = chromosome)) +
  geom_tile(aes(y = chromosome)) +
  scale_fill_manual(values = df$color, 
                    guide = guide_legend(title = "Chromosome", 
                                         direction = "vertical", 
                                         title.position = "top", 
                                         label.position = "right")) +
  theme_void()

# Print the plot
print(legend_plot)

dev.off()

# 8.3) Study case of inversions 29, 32, and 40 ----

# Building a network with only 6 nodes (two dosages for the study cases)

inversion_counts_df <- count_inversions(Data_d_l)
frequent_inversions <- select_frequent_inversions(inversion_counts_df, 
                                                  min_count = 7)

# Print the frequent inversions
print(frequent_inversions)

nodes_29_32_40 <- c("Inv_29_1", "Inv_29_2", "Inv_32_1", "Inv_32_2", 
                    "Inv_40_1", "Inv_40_2")

final_result_29_32_40 <- final_result |>
  filter(INV_1 %in% c("Inv_29_1", "Inv_29_2", "Inv_32_1", "Inv_32_2", 
                      "Inv_40_1", "Inv_40_2") & 
           INV_2 %in% c(c("Inv_29_1", "Inv_29_2", "Inv_32_1", "Inv_32_2", 
                          "Inv_40_1", "Inv_40_2")))

final_result_29_32_40 <- split(final_result_29_32_40, 
                               final_result_29_32_40$Line)

networks_29_32_40 <- mapply(sig_network_builder_lite, 
                          final_result_29_32_40,
                          MoreArgs = list(nodes = nodes_29_32_40,
                                          meta_nodes = metadata),
                          SIMPLIFY = FALSE)

names(networks_29_32_40) <- c("SUB_L_1034", "SUB_L_1192", "SUB_L_155",  
                              "SUB_L_444",  "SUB_L_502",
                              "SUB_L_541", "SUB_L_62", 
                              "SUB_L_664", "SUB_L_909")

plot_network_lite(networks_29_32_40, 
             edge_breaks = c(0.001, 0.008, 0.03, 0.05, 0.07, 0.10),
             range = c(1, 10))

# Heterogeneity analysis ----

# Analysis for non-dosage levels
mask1 <- c("Inv_29", "Inv_32", "Inv_40")
G_replicated_omnibus <- perform_heterogeneity_analysis(
  mask = mask1,
  data = Data_d_l,
  metadata = metadata,
  degrees_of_freedom = 4
)

# Analysis for dosage levels
mask2 <- c("Inv_29_1", "Inv_29_2", "Inv_32_1", "Inv_32_2", 
           "Inv_40_1", "Inv_40_2")
G_replicated <- perform_heterogeneity_analysis(
  mask = mask2,
  data = Data_p_a,
  metadata = metadata,
  degrees_of_freedom = 1
)

# 9) Network analysis ----

# Compare all pairs of networks
network_comparisons <- compare_all_pairs(networks)

# Condensing all the results into a single data frame

# Convert the list to a data frame
network_comparisons_df <- do.call(rbind, lapply(network_comparisons, function(x) {
  data.frame(
    is_isomorphic = x$is_isomorphic,
    avg_path_length_diff = x$avg_path_length_diff,
    clustering_coeff_diff = x$clustering_coeff_diff,
    spectrum_diff = x$spectrum_diff,
    stringsAsFactors = FALSE
  )
}))

# Add the comparison names as a new column
network_comparisons_df$comparison <- rownames(network_comparisons_df)

# Reset the row names
rownames(network_comparisons_df) <- NULL

# 9.1) Eigenvalue as a function of the number of genes in the inversion ----

# Extract the eigenvectors from each node in each line

eigenvectors <- lapply(networks, function(network) {
  V(network)$eigenvector
})

# Collapsing everything in a dataframe, keeping track of the line (names of the list)

eigenvectors_df <- do.call(rbind, lapply(names(eigenvectors), function(line) {
  data.frame(
    line = line,
    eigenvector = eigenvectors[[line]],
    stringsAsFactors = FALSE
  )
}))

# Incorporating the information of the nodes

eigenvectors_df$node <- rep(nodes, length(eigenvectors))

# Incorporating the information about the presence or absence of an inversion
# (using the degree information of the nodes from the networks)

degrees <- lapply(networks, function(network) {
  V(network)$degree
})

degrees_df <- do.call(rbind, lapply(names(degrees), function(line) {
  data.frame(
    line = line,
    degree = degrees[[line]],
    stringsAsFactors = FALSE
  )
}))

degrees_df$node <- rep(nodes, length(eigenvectors))

# Merging the eigenvectors and the degrees

eigenvectors_df <- merge(eigenvectors_df, degrees_df, by = c("line", "node"))

eigenvectors_df$Presence <- ifelse(eigenvectors_df$degree > 0, "Present", 
                                   "Absent")

# Incorporating the information about the number of genes in the inversion

eigenvectors_df$INV_ID <- sapply(strsplit(eigenvectors_df$node, "_"), `[`, 2)

eigenvectors_df <- eigenvectors_df %>%
  left_join(metadata[, c("INV_ID", "Genes")], 
            by = "INV_ID")

# Filtering out the nodes with degree 0

eigenvectors_df <- eigenvectors_df |>
  filter(degree > 0)

# Taking the log of the number of genes in the inversion
eigenvectors_df$log_INV_gene_number <- log(eigenvectors_df$Genes)

eigenvectors_df <- eigenvectors_df |> 
  mutate(line = as.factor(line),
         node = as.factor(node))

# 9.2) Adjust a robust linear model ----

# Model selection based on AIC
# Fit a linear regression model
# Define the full model
full_model <- lm(eigenvector ~ ., data = eigenvectors_df[, c("eigenvector", 
                                                              "degree", 
                                                              "log_INV_gene_number", 
                                                              "INV_ID",
                                                              "Genes", 
                                                              "line")])

# Perform backward stepwise selection based on AIC
# Define the scope for stepwise selection

scope_list <- list(
  lower = ~ 1,  # The simplest model, only the intercept
  upper = ~ degree + log_INV_gene_number + INV_ID + Genes + line + 
    I(degree^2) + I(log_INV_gene_number^2) + I(Genes^2) + 
    degree:log_INV_gene_number + degree:INV_ID + degree:Genes + 
    degree:line + 
    log_INV_gene_number:INV_ID + log_INV_gene_number:Genes + 
    log_INV_gene_number:line + 
    INV_ID:Genes + INV_ID:line + Genes:line  # An example of a more complex model
)

# Perform stepwise selection based on AIC, correcting the scope usage
model_selected_aic <- step(full_model, direction = "both", 
                           scope = scope_list, steps = 100000)

final_model <- rlm(eigenvector ~ degree + Genes + line +
                     I(log_INV_gene_number^2), data = eigenvectors_df)

# Estimating the p-values

rob.pvals(final_model)

# Checking for multicollinearity

vif_result <- vif(final_model)

# Print the selected model
summary(final_model)

# Plotting the model (color coding degree)

eigenvectors_df$predicted_eigenvector <- predict(final_model, eigenvectors_df)

pdf("Results/Plots/Genes_Eigenvector_Degree_Relationship.pdf", width = 8, height = 6)

ggplot(data = eigenvectors_df, aes(x = Genes, y = eigenvector, 
                                   color = degree)) +
  geom_point() + # Plot points
  geom_line(aes(y = predicted_eigenvector), color = "black") + # Add a regression line
  scale_color_gradient(low = "blue", high = "red") + # Color gradient for 'degree'
  labs(x = "Number of Genes (INV_gene_number)", y = "Eigenvector", 
       color = "Degree") +
  theme_bw() +
  theme(axis.text = element_text(size = 12, color = 'black'),
        axis.title = element_text(size = 14, color = 'black'),
        legend.text = element_text(size = 12, color = 'black'),
        legend.title = element_text(size = 14, color = 'black'),
        panel.grid = element_blank()) +
  xlim(0, 110)

dev.off()

# Plotting the model (color coding line)

pdf("Results/Plots/Genes_Eigenvector_Line_Relationship.pdf", width = 8, height = 6)

ggplot(data = eigenvectors_df, aes(x = Genes, y = eigenvector)) +
  geom_point(aes(color = line)) + # Color points by 'line'
  geom_line(aes(y = predicted_eigenvector, 
                group = line, 
                color = line)) + # Separate regression lines for each 'line'
  labs(x = "Number of Genes (INV_gene_number)", y = "Eigenvector", color = "Line") +
  theme_bw() +
  theme(axis.text = element_text(size = 12, color = 'black'),
        axis.title = element_text(size = 14, color = 'black'),
        legend.text = element_text(size = 12, color = 'black'),
        legend.title = element_text(size = 14, color = 'black'),
        panel.grid = element_blank()) +
  xlim(0, 110)

dev.off()

##' <<<<<<<<<<<<<<<<<<<<<<<<<<<< End of the script >>>>>>>>>>>>>>>>>>>>>>>>>>>>>