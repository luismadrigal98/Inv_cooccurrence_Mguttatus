dosage_splitter <- function(Data, probes_col = 1)
{
  #' Function toread the dosage levels and transform them into presence/absence
  #' 
  #' This function is designed to read the dosage levels of the plants and
  #' transform them into presence/absence. This is a necessary step to analyze
  #' the segregation distortion of the inversion and the co-occurrence of the
  #' inversions in a pairwise fashion.
  #' 
  #' @param Data A data frame with the dosage levels of the plants.
  #' @param probes_col The column where the probes are located.
  #' 
  #' @return A data frame with the presence/absence of the probes.
  #' ___________________________________________________________________________
  
  Probes <- Data[, probes_col]
  
  new_df <- data.frame()
  
  for (i in names(Data[, -probes_col]))
  {
    old_column <- Data[, i]
    df <- setNames(data.frame(ifelse(test = old_column == 1, 1, 0),
                              ifelse(test = old_column == 2, 1, 0)),
                   c(paste0(i, "_1"), paste0(i, "_2")))
    
    if (length(new_df) != 0) 
    {
      new_df <- cbind(new_df, df)
    }
    
    else
    {
      new_df <- df
    }
  }
  
  return(cbind(Probes, new_df))
}