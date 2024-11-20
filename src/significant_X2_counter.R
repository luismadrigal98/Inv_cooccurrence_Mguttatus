significant_X2_counter <- function(df, line) 
{
  #' This function will take a data frame containing the p_values associated with
  #' an independence test for the inversions, obtained from the main script 2,
  #' and it will count how many times a given inversion is not independent from
  #' other inversion dosage levels.
  #' 
  #' @param df A data frame containing the inversion pair, defined by two
  #' columns (INV_1 and INV_2), and the statistic and p_value associated with the
  #' omnibus X2 test.
  #' @param line The genetic family to be analyzed. This is necessary to match
  #' the results obtained with the output of main script 1.
  #' 
  #' @return A data frame containing the inversion ID, the line, 
  #' the mean p_value, and the number of times the X2 test was significant.
  #' ___________________________________________________________________________
  
  # Look for the inversions present in the line
  inversions <- unique(df$INV_1)
  
  # Subset the df to get only the relevant entries
  results <- lapply(inversions, function(inv) {
    # Get the p_values for the inversions
    # Remove the current inversion from the secondary columns (INV_1 != INV_2)
    p_values <- df %>%
      filter(INV_1 == inv & INV_2 != inv) %>%
      pull(p)
    
    # Get the mean p_value
    mean_p_value <- mean(p_values, na.rm = TRUE)
    
    # Count the number of significant X2 tests
    sig_count <- sum(p_values < 0.05, na.rm = TRUE)
    
    # Return the results
    return(data.frame(Inv = inv, Line = line, mean_p_value = mean_p_value, 
                      sig_count = sig_count))
  })
  
  # Consolidate the results
  results <- do.call(rbind, results)
  
  return(results)
}