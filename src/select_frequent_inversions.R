select_frequent_inversions <- function(inversion_counts_df, min_count = 2) 
{
  #' Select frequent inversions
  #' 
  #' This function filters out inversions that appear less than min_count times
  #' 
  #' @param inversion_counts_df A data frame with columns "inversion" and "count"
  #' @param min_count The minimum number of times an inversion must appear to be 
  #' considered frequent
  #' 
  #' @return A data frame with columns "inversion" and "count" containing only 
  #' frequent inversions
  #' ___________________________________________________________________________
  
  # Filter inversions that appear at least min_count times
  frequent_inversions <- inversion_counts_df[
    inversion_counts_df$count >= min_count, ]
  
  return(frequent_inversions)
}