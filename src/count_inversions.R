count_inversions <- function(data_list) 
{
  #' Count the occurrences of each inversion in a list of matrices
  #' 
  #' This function takes a list of matrices as input, where each matrix 
  #' represents a set of inversions. It counts the occurrences of each
  #' inversion across all matrices and returns a data frame with the inversion
  #' names and their corresponding counts.
  #' 
  #' @param data_list A list of matrices, where each matrix represents a set of inversions
  #' @return A data frame with the inversion names and their corresponding counts
  #' ___________________________________________________________________________
  
  # Initialize an empty list to store inversion counts
  inversion_counts <- list()
  
  # Loop through each matrix in the list
  for (name in names(data_list)) {
    # Get the row names (inversion names) for the current matrix
    inversions <- rownames(data_list[[name]])
    
    # Count the occurrences of each inversion
    for (inversion in inversions) {
      if (inversion %in% names(inversion_counts)) {
        inversion_counts[[inversion]] <- inversion_counts[[inversion]] + 1
      } else {
        inversion_counts[[inversion]] <- 1
      }
    }
  }
  
  # Convert the list to a data frame for easier manipulation
  inversion_counts_df <- data.frame(
    inversion = names(inversion_counts),
    count = unlist(inversion_counts)
  )
  
  return(inversion_counts_df)
}