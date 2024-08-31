## Co_occurrence of inversions in lines of Mimulus guttatus
#@ Supplementary Script 2
#@ Author: Luis J. Madrigal-Roca
#@ Date: 2024-05-27
#@ Dependencies: Affinity results from the main script

# ******************************************************************************
# Exploring the extreme affinity values not significant ----
# ______________________________________________________________________________

# Extracting the relevant observations (only present in line 1034)

df_inquiry <- final_result[((final_result$A_alpha < -5 | 
                               final_result$A_alpha > 5) & 
                              final_result$A_p > 0.05),]

affinity_inquiry <-results_affinity[['L_1034']]$occur_mat

# Contingency table visualization

matrix_inv <- as.matrix(df_inquiry[,c('INV_1','INV_2')])
colnames(matrix_inv) <- c('INV_1','INV_2')

disparity_tables <- apply(matrix_inv, 1, function(x, co_occurrences) {
  table(co_occurrences[, x[1]], co_occurrences[, x[2]])
}, affinity_inquiry, simplify = F)

# ******************************************************************************
# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> END <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
# ______________________________________________________________________________