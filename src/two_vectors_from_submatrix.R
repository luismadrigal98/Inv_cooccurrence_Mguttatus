two_vectors_from_submatrix <- function(submatrix)
{
  ## Getting the two vectors of interest (presence/absebce of a state)
  
  total <- sum(submatrix)
  
  intersection <- submatrix[2, 2]
  intersection <- rep(1, intersection)
  
  both_absents <- submatrix[1, 1]
  both_absents <- rep(0, both_absents)
  
  only1_present <- submatrix[1, 2]
  only2_present <- submatrix[2, 1]
  only1_absent <- total - (length(both_absents) + length(intersection) + 
                             only1_present)
  only2_absent <- total - (length(both_absents) + length(intersection) +
                             only2_present)
  
  vector1 <- c(both_absents, intersection, rep(1, only1_present), 
               rep(0, only1_absent))
  vector2 <- c(both_absents, intersection, rep(0, only2_absent), 
               rep(1, only2_present))
  
  return(list(vector1, vector2))
}