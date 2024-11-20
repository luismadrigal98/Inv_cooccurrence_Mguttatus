p_corrector <- function(df, var = "p", method = "BH")
{
  #' @title p_corrector
  #' @description This function corrects the p-values for multiple testing.
  #' @param df A data frame with the p-values to be corrected.
  #' @param var The name of the column containing the p-values.
  #' @param method The method to be used for the correction.
  #' @return A data frame with the corrected p-values.
  #' @examples
  #' p_corrector(df = df, var = "p", method = "BH")
  #' ___________________________________________________________________________
  
  df <- df |> mutate(p_corrected = p.adjust(df[, var], method = method))
  
  return(df)
}