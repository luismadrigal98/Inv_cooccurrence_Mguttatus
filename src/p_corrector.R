p_corrector <- function(df, var = "p", FDR = 0.1,
                        pi0.method = "bootstrap", adj = 1.2)
{
  #' @title p_corrector
  #' @description This function corrects the p-values for multiple testing.
  #' @param df A data frame with the p-values to be corrected.
  #' @param var The name of the column containing the p-values.
  #' @param FDR The false discovery rate to be used in the correction.
  #' @param pi0.method The method to estimate the proportion of true null hypotheses.
  #' @param adj The adjustment factor to be used in the correction.
  #' @return A data frame with the corrected p-values.
  #' @examples
  #' p_corrector(df = df, var = "p", method = "BH")
  #' ___________________________________________________________________________
  
  df <- df |> mutate(q_value = qvalue(df[, var], fdr.level = FDR,
                                      pi0.method = pi0.method, 
                                      adj = adj)$qvalues)
  
  return(df)
}