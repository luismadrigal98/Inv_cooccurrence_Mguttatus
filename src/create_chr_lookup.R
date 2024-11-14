create_chr_lookup <- function(metadata) 
{
  #' Create a lookup table for chromosome names
  #' 
  #' This function creates a lookup table for chromosome names based on the metadata
  #' 
  #' @param metadata A data frame containing the metadata
  #' 
  #' @return A named vector with the chromosome names
  #' ___________________________________________________________________________
  
  # Split the Line column if it contains multiple lines
  metadata_split <- data.frame(
    INV_ID = metadata$INV_ID,
    Chr = metadata$Chr,
    stringsAsFactors = FALSE
  )
  return(setNames(metadata_split$Chr, paste0("Inv_", metadata_split$INV_ID)))
}