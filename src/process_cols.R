process_cols <- function(df, cols, value_name, prefix) 
{
  #' Process columns for contingency_expanded.
  #' ___________________________________________________________________________
  
  df |>
    dplyr::select(all_of(c("INV_1", "INV_2", "X2", "p", "p_corrected", "ID1", 
                           "ID2", "chrom1", 
                           "chrom2", "INV_combination", cols))) |>
    pivot_longer(cols = all_of(cols), names_to = "cell_location", 
                 values_to = value_name) |>
    mutate(cell_location = gsub(prefix, "", cell_location))
}