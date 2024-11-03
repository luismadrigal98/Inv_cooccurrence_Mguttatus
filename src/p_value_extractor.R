p_value_extractor <- function(models, expected) 
{
  #' Extract p-values from chi-square tests and submatrix tests
  #' 
  #' This function extracts p-values from chi-square tests and submatrix tests
  #' for a set of models.
  #' 
  #' @param models A 3D array with the contingency tables of the models.
  #' 
  #' @param expected A 2x2 matrix with the expected values for the chi-square test.
  #' 
  #' @return A list with the p-values for the chi-square tests, affinity tests, 
  #' and Jaccard tests.
  #' ___________________________________________________________________________
  
  # Calculate chi-square test p-values
  n_models <- dim(models)[3]
  p_values_X2 <- parallel::mclapply(1:n_models, function(i) {
    chisq.test(models[,,i], p = expected)$p.value
  })
  
  # Extract all submatrices at once
  submatrices <- list(
    sub1 = parallel::mclapply(1:n_models, 
                              function(i) submatrix_extractor(models[,,i], 3, 3)),
    sub2 = parallel::mclapply(1:n_models, 
                              function(i) submatrix_extractor(models[,,i], 2, 2)),
    sub3 = parallel::mclapply(1:n_models, 
                              function(i) submatrix_extractor(models[,,i], 2, 3)),
    sub4 = parallel::mclapply(1:n_models, 
                              function(i) submatrix_extractor(models[,,i], 3, 2))
  )
  
  # Calculate affinity p-values for all submatrices
  p_values_affinity <- lapply(submatrices, process_submatrix_tests, 
                              test_type = "affinity")
  
  # Calculate Jaccard p-values for all submatrices
  p_values_jaccard <- lapply(submatrices, process_submatrix_tests, 
                             test_type = "jaccard")
  
  # Return results
  list(
    X2 = unlist(p_values_X2),
    affinity = p_values_affinity,
    jaccard = p_values_jaccard
  )
}