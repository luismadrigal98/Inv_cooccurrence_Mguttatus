marginal_compliance_checking <- function(data) 
{
  #' Check if the marginal distribution of two variables is zero
  #' 
  #' @param data A data frame
  #' 
  #' @return A list of problematic pairs of variables
  #' ___________________________________________________________________________
  
  problematic_pairs <- list()
  for (i in 1:(nrow(data) - 1)) {
    for (j in (i + 1):nrow(data)) {
      result <- check_zero_marginals(data, i, j)
      if (!is.null(result)) {
        problematic_pairs <- append(problematic_pairs, list(result))
      }
    }
  }
  
  return(problematic_pairs)
}