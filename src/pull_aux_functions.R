pull_aux_functions <- function(src_dir)
{
  #' Pull auxiliary functions from a source directory
  #' 
  #' This function pulls all auxiliary functions from a source directory
  #' and sources them into the current environment.
  #' 
  #' @param src_dir The source directory containing the auxiliary functions.
  #' 
  #' @return invisible
  #' ___________________________________________________________________________
  
  # List all files in the source directory
  files <- list.files(src_dir, pattern = "\\.R$", full.names = TRUE)
  
  # Source all the files
  lapply(files, source)
  
  return(invisible())
}