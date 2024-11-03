table_sim_null <-function(rep, n_cells, n_ind, expected)
{
  #' Simulate a table of counts under the null hypothesis of independence
  #' 
  #' @param rep Number of replicates
  #' @param n_cells Number of cells in the table
  #' @param n_ind Number of individuals to distribute across cells
  #' @param expected Expected proportions of each cell
  #' 
  #' @return A 3D array of counts. Third dimension is the replicate tables.
  #' ___________________________________________________________________________
  
  simulated_m <- foreach(i = 1:rep) %dopar% 
    {
      matrix(tabulate(sample(x = 1:n_cells, size = n_ind, replace = T, 
                             prob = expected)), nrow = 3)
    }
  
  simulated_m <- array(unlist(simulated_m), 
                       dim = c(nrow(simulated_m[[1]]), ncol(simulated_m[[1]]), 
                               length(simulated_m)))
  
  return(simulated_m)
}