G_per_contrast_per_line <- function(df_contrast, Data_p_a,
                                    correction = "williams")
{
  #' This function will get the G value associated to an independence test for 
  #' each contrast and each line. The returned result is a df that will take as
  #' basis the df_contrast input. So, this function will add as many columns as
  #' lines are int he Data_p_a list.
  #' 
  #' @param df_contrast A data frame with the contrasts to be tested. It should
  #' have the columns "Inv_1" and "Inv_2".
  #' 
  #' @param Data_p_a A list with the matrices of presence/absence of the inversions.
  #' Each slot represent a line. So, the columns added to the data frame correspond to
  #' the G_value for each contrast in a line.
  #' 
  #' @param correction A character string with the name of the correction to be used
  #' in the G test. The default is "williams".
  #' 
  #' @return A data frame with the G values for each contrast in each line.
  #' ___________________________________________________________________________
  
  Lines <- names(Data_p_a)
  
  G_values_as_list <- lapply(Lines, function(line)
  {
    G_values <- apply(df_contrast, 1, function(contrast)
    {
      G_value <- GTest(Data_p_a[[line]][contrast[1],], 
                       Data_p_a[[line]][contrast[2],],
                       correct = correction)
      return(G_value$statistic)
    })
    
    return(G_values)
  })
  
  names(G_values_as_list) <- Lines
  
  for(line in Lines)
  {
    df_contrast <- dplyr::mutate(df_contrast, 
                                 !!paste0("G_", line) := G_values_as_list[[line]])
  }
  
  return(df_contrast)
}