support_counter <- function(supp_X2_posthoc,
                            supp_Affinity,
                            supp_Jaccard)
{
  #' @title support_counter
  #' @description This function will count the number of approaches that support
  #' the hypothesis of epistasis between inversions in Mimulus guttatus.
  #' 
  #' @param supp_X2_posthoc A matrix with the results of the Chi-square post-hoc
  #' test.
  #' @param supp_Affinity A matrix with the results of the Affinity test.
  #' @param supp_Jaccard A matrix with the results of the Jaccard test.
  #'
  #' @return A matrix with the number of approaches that support the hypothesis
  #' of epistasis between inversions in Mimulus guttatus.
  #' ___________________________________________________________________________
  
  # Remove rows and columns with NA values
  row_na <- apply(supp_X2_posthoc, 1, function(x) any(!is.na(x)))
  col_na <- apply(supp_X2_posthoc, 2, function(x) any(!is.na(x)))
  supp_X2_posthoc <- supp_X2_posthoc[row_na, col_na]
  
  names <- list(rownames(supp_X2_posthoc), colnames(supp_X2_posthoc))
  
  supp_X2_posthoc <- matrix(as.numeric(supp_X2_posthoc < 0.05), 
                            nrow = nrow(supp_X2_posthoc), 
                            ncol = ncol(supp_X2_posthoc))
  supp_Affinity <- matrix(as.numeric(supp_Affinity < 0.05), 
                          nrow = nrow(supp_Affinity), 
                          ncol = ncol(supp_Affinity))
  supp_Jaccard <- matrix(as.numeric(supp_Jaccard < 0.05), 
                         nrow = nrow(supp_Jaccard), 
                         ncol = ncol(supp_Jaccard))
  
  concensus <- supp_X2_posthoc + supp_Affinity + supp_Jaccard
  dimnames(concensus) <- names
  
  return(concensus)
}