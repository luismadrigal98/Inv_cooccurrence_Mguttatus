contingency_expanded <- function(cont_res_df) 
{
  #' Expand the contingency table results data frame
  #' 
  #' This function takes the contingency table results data frame and expands it
  #' into a long format, with separate columns for each cell location and the
  #' corresponding values.
  #' 
  #' @param cont_res_df A data frame containing the contingency table results
  #' 
  #' @return A data frame in long format with separate columns for each cell
  #' location and the corresponding values
  #' ___________________________________________________________________________
  
  # Define the sets of columns
  cols_list <- list(
    dev = list(cols = c("dev_2r_2c", "dev_2r_3c", "dev_3r_2c", "dev_3r_3c"), 
               value_name = "deviation_value", prefix = "dev_"),
    sr = list(cols = c("sr_2r_2c", "sr_2r_3c", "sr_3r_2c", "sr_3r_3c"), 
              value_name = "standardized_residual", prefix = "sr_"),
    rel_cont = list(cols = c("rel_cont_2r_2c", "rel_cont_2r_3c", 
                             "rel_cont_3r_2c", "rel_cont_3r_3c"), 
                    value_name = "relative_contribution", prefix = "rel_cont_")
  )
  
  # Apply the helper function to each set of columns
  df_list <- lapply(cols_list, function(x) process_cols(cont_res_df, x$cols, 
                                                        x$value_name, x$prefix))
  
  # Merge the data frames together
  df <- Reduce(function(x, y) merge(x, y, by = c("INV_1", "INV_2", "X2", "p",
                                                 "p_corrected",
                                                 "ID1", "ID2", "chrom1", 
                                                 "chrom2", "INV_combination", 
                                                 "cell_location")), df_list)
  
  return(df)
}