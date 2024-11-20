perform_heterogeneity_analysis <- function(mask, 
                                           data, 
                                           metadata, 
                                           degrees_of_freedom) 
{
  # Helper function to generate and filter pairs
  generate_valid_pairs <- function(mask_vector) {
    grid <- expand.grid(mask_vector, mask_vector)
    pairs <- grid[grid$Var1 != grid$Var2, ]
    pairs <- pairs[!duplicated(t(apply(pairs, 1, sort))), ]
    names(pairs) <- c("Inv_1", "Inv_2")
    
    pairs |>
      mutate(valid = mapply(is.valid, 
                            as.character(Inv_1), 
                            as.character(Inv_2),
                            MoreArgs = list(metadata = metadata,
                                            fields = c(1, 3)),
                            SIMPLIFY = TRUE)) |>
      dplyr::filter(valid) |>
      dplyr::select(-valid)
  }
  
  # Subset data
  data_subset <- lapply(data, function(x) x[mask,])
  
  # Generate valid pairs
  unique_pairs <- generate_valid_pairs(mask)
  
  # Calculate G values and perform replicated test
  G_results <- G_per_contrast_per_line(unique_pairs, data_subset) |>
    mutate(Inv_1 = as.character(Inv_1),
           Inv_2 = as.character(Inv_2)) |>
    replicated_G_test(degrees_of_fred_replicates = degrees_of_freedom,
                      Data = data_subset)
  
  return(G_results)
}