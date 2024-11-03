G_calculator <- function(row, dosage_matrix, exp = c(1/4, 1/2, 1/4)) 
{
  #' Calculate the G statistic for a single row of a dosage matrix
  #' 
  #' @param row A row index of the dosage matrix
  #' @param dosage_matrix A dosage matrix
  #' @param exp A vector of expected probabilities for 0, 1, and 2 (dosage of the
  #' inversion)
  #' 
  #' @return A data frame with the following columns:
  #'  - Inv: The name of the row
  #'  - G: The G statistic
  #'  - df: The degrees of freedom
  #'  - p_value: The p-value
  #'  __________________________________________________________________________
  
  INV <- rownames(dosage_matrix)[row]
  
  # Manually count the occurrences of 0, 1, and 2
  counts <- c(sum(dosage_matrix[row, ] == 0),
              sum(dosage_matrix[row, ] == 1),
              sum(dosage_matrix[row, ] == 2))
  
  # Print the counts for debugging
  print(counts)
  
  # Ensure that counts and exp have the same length
  if (length(counts) != length(exp)) {
    warning("Observed counts and expected probabilities do not match in length.")
    return(c(Inv = INV, G = NA, df = NA, p_value = NA))
  }
  
  test_results <- GTest(counts, p = exp, correct = 'williams')
  
  return(data.frame(Inv = INV, G = test_results$statistic, 
                    df = test_results$parameter,
                    p_value = test_results$p.value))
}