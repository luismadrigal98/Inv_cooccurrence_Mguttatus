set_environment <- function(required_pckgs, 
                            personal_seed = as.numeric(Sys.time()),
                            parallel_backend = FALSE)
{
  #' This fucntion will set up the working environment for performing all the
  #' analysis. It will load all the required packages, set the seed for
  #' reproducibility and set the parallel backend if required.
  #' 
  #' @param required_pckgs A character vector of the required packages.
  #' 
  #' @param personal_seed An integer to set the seed for reproducibility. Default
  #' set to system time.
  #' 
  #' @param parallel_backend A logical value to set the parallel backend. Default
  #' set to FALSE.
  #' 
  #' @return invisible
  #' ___________________________________________________________________________
  
  ## Loading the required libraries ----
  
  message("Loading the required libraries")
  
  tryCatch(
  {
    for(pckg in required_pckgs)
    {
      if(!require(pckg, character.only = TRUE))
      {
        install.packages(pckg)
        library(pckg, character.only = TRUE)
      }
      else
      {
        library(pckg, character.only = TRUE)
      }
    }
    
    message("The required libraries have been loaded")
  }, error = function(e) {
    message("Some packages cannot be installed through CRAN, trying with remotes")
    print(e)
  })
  
  ## Setting the seed ----
  set.seed(personal_seed)
  
  ## Setting the parallel backend ----
  if(parallel_backend == TRUE)
  {
    require(doParallel)
    
    ## 1.2) Setting the parallel backend
    cl <<- makeCluster(detectCores() - 1)
    registerDoParallel(cl)
  }
  
  return(invisible())
}