x2p_to_square <- function(x2_result, col_names = c("INV_1", "INV_2"), 
                          p_col = 'p') 
{
  #' Convert the output of the chisq.test() function to a square matrix
  #' of p-values.
  #' 
  #' @param x2_result A data frame with the output of the chisq.test() function.
  #' @param col_names A vector of two strings with the names of the columns
  #'                 that contain the two variables of the chi-squared test.
  #'                 Default is c("INV_1", "INV_2").
  #'                 The order of the columns is important.
  #'                 The first column will be the row names of the matrix.
  #'                 The second column will be the column names of the matrix.
  #'                 The p-values will be the values of the matrix.
  #' @param p_col A string with the name of the column that contains the p-values.
  #'             Default is 'p'.
  #'             
  #' @return A square matrix with the p-values of the chi-squared test.
  #' ___________________________________________________________________________
  
  x2_result <- as.data.frame(x2_result)
  
  unique_values <- unique(c(x2_result[, col_names[1]], 
                            x2_result[, col_names[2]]))
  
  p_values <- matrix(nrow = length(unique_values), 
                     ncol = length(unique_values))
  
  rownames(p_values) <- colnames(p_values) <- unique_values
  
  for (i in unique_values) {
    for (j in unique_values) {
      p_values[i, j] <- x2_result[x2_result[, col_names[1]] == i & 
                                    x2_result[, col_names[2]] == j , p_col]
    }
  }
  
  return(p_values)
}