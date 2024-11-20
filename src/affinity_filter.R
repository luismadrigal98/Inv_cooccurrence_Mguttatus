affinity_filter <- function(affinity_res, metadata, filter_LG = T, 
                            filter_CS = T)
{
  #' @title affinity_filter
  #' @description Filters the affinity results to remove entries that belong to
  #' the same chromosome or to the same inversion with mutually exclusive states
  #' @param affinity_res A list with the affinity results
  #' @param metadata A data frame with the metadata of the inversions
  #' @param filter_LG A logical value to filter the entries that belong to the
  #' same chromosome
  #' @param filter_CS A logical value to filter the entries that belong to the
  #' same inversion with mutually exclusive states
  #' @return A data frame with the filtered affinity results
  #' ___________________________________________________________________________
  
  df <- affinity_res$all
  
  df <- df |> separate(col = entity_1, into = c("Generic_name_1", "INV_ID_1", 
                                                "State_1"), sep = '_', 
                       remove = F) |>
    separate(col = entity_2, into = c("Generic_name_2", "INV_ID_2", 
                                      "State_2"), sep = '_', 
             remove = F) |>
    dplyr::select(-Generic_name_1, -State_1, -Generic_name_2, -State_2)
  
  df <- df |> 
    left_join(metadata |> dplyr::select(INV_ID, Chr), by = c("INV_ID_1" = "INV_ID")) |> 
    rename(chrom1 = Chr)
  
  df <- df |> 
    left_join(metadata |> dplyr::select(INV_ID, Chr), by = c("INV_ID_2" = "INV_ID")) |> 
    rename(chrom2 = Chr)
  
  if (filter_LG == T)
  {
    df <- df |> filter(chrom1 != chrom2)
  }
  
  if (filter_CS == T)
  {
    df <- df |> filter(INV_ID_1 != INV_ID_2)
  }
  
  return(df)
}