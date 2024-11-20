compare_all_pairs <- function(networks) 
{
  #' Compare all pairs of networks in a list
  #' 
  #' @param networks A list of networks
  #' @return A list of comparisons
  #' ___________________________________________________________________________
  
  n <- length(networks)
  comparisons <- list()
  
  for (i in 1:(n-1)) {
    for (j in (i+1):n) {
      comparison <- compare_networks(networks[[i]], networks[[j]])
      comparisons[[paste0(names(networks)[i], "_vs_", 
                          names(networks)[j])]] <- comparison
    }
  }
  
  return(comparisons)
}
