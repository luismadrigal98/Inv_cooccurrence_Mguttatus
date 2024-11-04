x2_posthoc_calculator <- function(i, j, dosage_matrix, expectation)
{
  #' This function will calculate the post-hoc chi-square test for two entities
  #' 
  #' @param i The index of the first entity
  #' @param j The index of the second entity
  #' @param dosage_matrix The matrix of dosage levels
  #' @param expectation The expected values for the chi-square test
  #' 
  #' @return A data frame with the residuals and p-values
  #' ___________________________________________________________________________
  
  result <- list()
  
  INV1_name <- rownames(dosage_matrix)[i]
  INV2_name <- rownames(dosage_matrix)[j]
  
  observed <- matrix(NA, 3, 3)
  
  for (x in 0:2) {
    for (y in 0:2) {
      observed[x + 1, y + 1] <- sum(dosage_matrix[i,] == x & 
                                      dosage_matrix[j,] == y)
    }
  }
  
  rownames(observed) <- c(paste0(INV1_name, "_0"), paste0(INV1_name, "_1"), 
                          paste0(INV1_name, "_2"))
  colnames(observed) <- c(paste0(INV2_name, "_0"), paste0(INV2_name, "_1"), 
                          paste0(INV2_name, "_2"))
  
  res <- chisq.posthoc.test(observed, method = 'none', p = expectation, 
                            simulate.p.value = T, 
                            B = 10000)
  
  df_melt <- melt(res, id.vars = c("Dimension", "Value"))
  
  # Split the data into two data frames
  df_residuals <- df_melt |> filter(Value == "Residuals") |> 
    dplyr::rename(Residual = value) |> dplyr::select(-Value)
  df_pvalues <- df_melt |> filter(Value == "p values") |> 
    dplyr::rename(p_value = value) |> dplyr::select(-Value)
  
  # Join the data frames
  df_final <- full_join(df_residuals, df_pvalues, 
                        by = c("Dimension", "variable"))
  
  # Rename the columns
  colnames(df_final) <- c("INV1", "INV2", "Residual", "p_value")
  
  return(df_final)
}