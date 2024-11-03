expected_mod <- function(expected, survival)
{
  #' Calculate the expected number of deaths for a given survival probability
  #' 
  #' @param expected A vector of expected number of deaths
  #' @param survival A vector of survival probabilities
  #' 
  #' @return A vector of expected number of deaths for a given survival probability
  #' ___________________________________________________________________________
  
  return(expected * survival / sum(expected * survival))
}