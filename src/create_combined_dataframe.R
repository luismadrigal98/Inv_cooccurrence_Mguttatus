create_combined_dataframe <- function(goodness_of_fit, 
                                      omnibus_X2_results, 
                                      metadata) 
{
  #' Combine goodness of fit and omnibus X2 results into a single data frame
  #' 
  #' This function combines the goodness of fit and omnibus X2 results into a 
  #' single data frame.
  #' 
  #' @param goodness_of_fit A data frame with goodness of fit from the SD analysis 
  #'                        for each inversion
  #' @param omnibus_X2_results A list of data frames with omnibus X2 results for 
  #'                           each line
  #' @param metadata A data frame with metadata for each inversion
  #' 
  #' @return A data frame with combined results
  #' ___________________________________________________________________________
  
  # Create chromosome lookup
  chr_lookup <- create_chr_lookup(metadata)
  
  # Create empty list to store results
  combined_results <- list()
  row_counter <- 1
  
  # Process each line's data
  for(line_name in names(omnibus_X2_results)) {
    # Get line-specific data
    line_data <- omnibus_X2_results[[line_name]]
    line_gof <- goodness_of_fit[goodness_of_fit$Line == line_name, ]
    
    # Process each row in the X2 results
    for(i in 1:nrow(line_data)) {
      # Get current row
      row <- line_data[i, ]
      
      # Check if inversions are on different chromosomes
      if(!are_different_chr(row$INV_1, row$INV_2, chr_lookup)) {
        next  # Skip this pair if they're on the same chromosome
      }
      
      # Get goodness of fit data for both inversions
      inv1_gof <- line_gof[line_gof$Inv == row$INV_1, ]
      inv2_gof <- line_gof[line_gof$Inv == row$INV_2, ]
      
      # Create a list with all information
      result_list <- list(
        Line = line_name,
        Inv1 = row$INV_1,
        Inv2 = row$INV_2,
        Inv1_Chr = chr_lookup[row$INV_1],
        Inv2_Chr = chr_lookup[row$INV_2],
        X2 = row$X2,
        X2_p_value = row$p,
        # Goodness of fit metrics for Inv1
        Inv1_G = if(nrow(inv1_gof) > 0) inv1_gof$G else NA,
        Inv1_df = if(nrow(inv1_gof) > 0) inv1_gof$df else NA,
        Inv1_p_value = if(nrow(inv1_gof) > 0) inv1_gof$p_value else NA,
        Inv1_p_corrected = if(nrow(inv1_gof) > 0) inv1_gof$p_corrected else NA,
        Inv1_SDV = if(nrow(inv1_gof) > 0) inv1_gof$SDV else NA,
        # Goodness of fit metrics for Inv2
        Inv2_G = if(nrow(inv2_gof) > 0) inv2_gof$G else NA,
        Inv2_df = if(nrow(inv2_gof) > 0) inv2_gof$df else NA,
        Inv2_p_value = if(nrow(inv2_gof) > 0) inv2_gof$p_value else NA,
        Inv2_p_corrected = if(nrow(inv2_gof) > 0) inv2_gof$p_corrected else NA,
        Inv2_SDV = if(nrow(inv2_gof) > 0) inv2_gof$SDV else NA,
        # Deviation metrics
        dev_1r_1c = row$dev_1r_1c,
        dev_2r_1c = row$dev_2r_1c,
        dev_3r_1c = row$dev_3r_1c,
        dev_1r_2c = row$dev_1r_2c,
        dev_2r_2c = row$dev_2r_2c,
        dev_3r_2c = row$dev_3r_2c,
        dev_1r_3c = row$dev_1r_3c,
        dev_2r_3c = row$dev_2r_3c,
        dev_3r_3c = row$dev_3r_3c,
        # Standardized residuals
        sr_1r_1c = row$sr_1r_1c,
        sr_2r_1c = row$sr_2r_1c,
        sr_3r_1c = row$sr_3r_1c,
        sr_1r_2c = row$sr_1r_2c,
        sr_2r_2c = row$sr_2r_2c,
        sr_3r_2c = row$sr_3r_2c,
        sr_1r_3c = row$sr_1r_3c,
        sr_2r_3c = row$sr_2r_3c,
        sr_3r_3c = row$sr_3r_3c,
        # Relative contributions
        rel_cont_1r_1c = row$rel_cont_1r_1c,
        rel_cont_2r_1c = row$rel_cont_2r_1c,
        rel_cont_3r_1c = row$rel_cont_3r_1c,
        rel_cont_1r_2c = row$rel_cont_1r_2c,
        rel_cont_2r_2c = row$rel_cont_2r_2c,
        rel_cont_3r_2c = row$rel_cont_3r_2c,
        rel_cont_1r_3c = row$rel_cont_1r_3c,
        rel_cont_2r_3c = row$rel_cont_2r_3c,
        rel_cont_3r_3c = row$rel_cont_3r_3c
      )
      
      combined_results[[row_counter]] <- result_list
      row_counter <- row_counter + 1
    }
  }
  
  # Convert list to data frame
  if(length(combined_results) > 0) {
    combined_df <- do.call(rbind, lapply(combined_results, data.frame))
    return(combined_df)
  } else {
    return(data.frame())  # Return empty data frame if no pairs found
  }
}