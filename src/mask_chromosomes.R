mask_chromosomes <- function(p_matrix, metadata)
{
  #' Mask the p_matrix by setting the p-values to NA if the two inversions are 
  #' on different chromosomes.
  #' 
  #' @param p_matrix A matrix of p-values.
  #' @param metadata A data frame containing metadata for the inversions.
  #' 
  #' @return A matrix of p-values with the p-values for inversions on different
  #' chromosomes set to NA.
  #' ___________________________________________________________________________
  
  row_names <-  rownames(p_matrix)
  row_names <- as.character(sapply(strsplit(row_names, "_"), `[`, 2))
  
  col_names <- colnames(p_matrix)
  col_names <- as.character(sapply(strsplit(col_names, "_"), `[`, 2))
  
  for (i in 1:length(row_names))
  {
    for (j in 1:length(col_names))
    {
      if (metadata[metadata$INV_ID == row_names[i], "Chr"] != 
          metadata[metadata$INV_ID == col_names[j], "Chr"])
      {
        p_matrix[i, j] <- NA
      }
    }
  }
  
  return(p_matrix)
}