check_zero_marginals <- function(data, row_index, col_index) 
{
  #' Check for zero marginals in a 2x2 contingency table
  #' 
  #' @param data A matrix of data
  #' 
  #' @param row_index The row index of the data matrix
  #' 
  #' @param col_index The column index of the data matrix
  #' 
  #' @return A list containing the row name, column name, and observed values 
  #' of the 2x2 table if zero marginals are found, otherwise NULL.
  #' ___________________________________________________________________________
  
  row_name <- rownames(data)[row_index]
  col_name <- rownames(data)[col_index]
  
  observed <- matrix(NA, 3, 3)
  
  for (x in 0:2) {
    for (y in 0:2) {
      observed[x + 1, y + 1] <- sum(data[row_index, ] == x & data[col_index, ] == y)
    }
  }
  
  # Check for zero marginals
  if (any(rowSums(observed) == 0) || any(colSums(observed) == 0)) {
    return(list(row_name = row_name, col_name = col_name, observed = observed))
  } else {
    return(NULL)
  }
}