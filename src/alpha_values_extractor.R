alpha_values_extractor <- function(affinity_results) 
{
  #' Extracts the alpha values from the affinity results
  #' 
  #' @param affinity_results The affinity results object
  #' 
  #' @return A matrix of alpha values
  #' ___________________________________________________________________________
  
  df <- affinity_results$all
  inv_ID <- names(affinity_results$occur_mat)
  
  A_matrix <- matrix(nrow = length(inv_ID), ncol = length(inv_ID))
  
  dimnames(A_matrix) <- list(inv_ID, inv_ID)
  
  for (i in inv_ID) {
    for(j in inv_ID) {
      if(i == j) {
        A_matrix[i, j] <- 1 # Maximum value (given the capped behavior of affinity calculator)
      } else if (substr(i, 1, nchar(i) - 2) == substr(j, 1, nchar(j) - 2)) {
        A_matrix[i, j] <- 1 # A plant cannot be both homozygous and heterozygous at the same time
      } else {
        alpha_p <- df[df[, "entity_1"] == i & df[, "entity_2"] == j, "p_value"]
        if (length(alpha_p) == 0) {
          alpha_p <- df[df[, "entity_1"] == j & df[, "entity_2"] == i, "p_value"]
        }
        A_matrix[i, j] <- ifelse(length(alpha_p) == 0, NA, alpha_p)
      }
    }
  }
  
  return(A_matrix)
}