submatrix_extractor <- function(matrix, row_to_exclude, col_to_exclude)
{
  #' Extracts a submatrix from a matrix by excluding a row and a column
  #' 
  #' @param matrix A matrix from which to extract the submatrix
  #' 
  #' @param row_to_exclude The row to exclude from the matrix
  #' 
  #' @param col_to_exclude The column to exclude from the matrix
  #' 
  #' @return A submatrix of the input matrix with the specified row and column 
  #' excluded
  #' ___________________________________________________________________________
  
  ## Extracting the submatrix of interest
  submatrix <- matrix[-row_to_exclude, -col_to_exclude]
  
  return(submatrix)
}