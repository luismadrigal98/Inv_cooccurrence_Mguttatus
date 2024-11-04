x2_calculator <- function(i, j, dosage_matrix)
{
  #' This function calculates the X2 test for independence between two inversions
  #' in the dataset. The function will return a dataframe with the results of the
  #' test.
  #' 
  #' @param i The row index of the first inversion
  #' @param j The row index of the second inversion
  #' @param dosage_matrix The matrix containing the dosage levels of the inversions
  #' 
  #' @return A dataframe with the results of the X2 test
  #' 
  #' @note In some instances, the result test will not have some of the elements.
  #' This happen because for INV_49, there is one dosage level that is not present,
  #' and this result in a contingency analysis of a 2x3 table. This is not a problem
  #' for the X2 test, but it is for the post-hoc analysis. We will preserve the
  #' structure of the output, but for those missing elements, the values will be
  #' NA's.
  #' ___________________________________________________________________________
  
  INV1_name <- rownames(dosage_matrix)[i]
  INV2_name <- rownames(dosage_matrix)[j]
  
  observed <- table(dosage_matrix[i,], dosage_matrix[j,])
  
  test <- chisq.test(x = observed, simulate.p.value = T,
                     B = 10000)
  
  x2 <- test$statistic
  p_value <- test$p.value
  dev <- test$observed - test$expected
  standardized_residuals <- dev / sqrt(test$expected)
  relative_contribution <- (dev^2 / (test$expected * x2)) * 100
  
  return(data.frame(INV_1 = INV1_name,
                    INV_2 = INV2_name,
                    X2 = x2,
                    p = p_value,
                    dev_1r_1c = dev[1, 1],
                    dev_2r_1c = dev[2, 1],
                    dev_3r_1c = if(nrow(dev) == 3) dev[3, 1] else NA,
                    dev_1r_2c = dev[1, 2],
                    dev_2r_2c = dev[2, 2],
                    dev_3r_2c = if(nrow(dev) == 3) dev[3, 2] else NA,
                    dev_1r_3c = if(ncol(dev) == 3) dev[1, 3] else NA,
                    dev_2r_3c = if(ncol(dev) == 3) dev[2, 3] else NA,
                    dev_3r_3c = if(ncol(dev) == 3 & nrow(dev) == 3) dev[3, 3] else NA,
                    sr_1r_1c = standardized_residuals[1, 1],
                    sr_2r_1c = standardized_residuals[2, 1],
                    sr_3r_1c = if(nrow(standardized_residuals) == 3) standardized_residuals[3, 1] else NA,
                    sr_1r_2c = standardized_residuals[1, 2],
                    sr_2r_2c = standardized_residuals[2, 2],
                    sr_3r_2c = if(nrow(standardized_residuals) == 3) standardized_residuals[3, 2] else NA,
                    sr_1r_3c = if(ncol(standardized_residuals) == 3) standardized_residuals[1, 3] else NA,
                    sr_2r_3c = if(ncol(standardized_residuals) == 3) standardized_residuals[2, 3] else NA,
                    sr_3r_3c = if (ncol(standardized_residuals) == 3 & 
                                   nrow(standardized_residuals) == 3) standardized_residuals[3, 3] else NA,
                    rel_cont_1r_1c = relative_contribution[1, 1],
                    rel_cont_2r_1c = relative_contribution[2, 1],
                    rel_cont_3r_1c = if(nrow(relative_contribution) == 3) relative_contribution[3, 1] else NA,
                    rel_cont_1r_2c = relative_contribution[1, 2],
                    rel_cont_2r_2c = relative_contribution[2, 2],
                    rel_cont_3r_2c = if(nrow(relative_contribution) == 3) relative_contribution[3, 2] else NA,
                    rel_cont_1r_3c = if(ncol(relative_contribution) == 3) relative_contribution[1, 3] else NA,
                    rel_cont_2r_3c = if(ncol(relative_contribution) == 3) relative_contribution[2, 3] else NA,
                    rel_cont_3r_3c = if(ncol(relative_contribution) == 3 &
                                        nrow(relative_contribution) == 3) relative_contribution[3, 3] else NA))
}