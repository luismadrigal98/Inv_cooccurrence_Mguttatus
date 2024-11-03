submatrix_extractor <- function(matrix, row_to_exclude, col_to_exclude)
{
  ## Extracting the submatrix of interest
  submatrix <- matrix[-row_to_exclude, -col_to_exclude]
  
  return(submatrix)
}