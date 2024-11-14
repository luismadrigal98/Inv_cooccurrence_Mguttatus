are_different_chr <- function(inv1, inv2, chr_lookup) 
{
  #' Check if two inversions are on different chromosomes
  #' 
  #' @param inv1 Inversion 1
  #' @param inv2 Inversion 2
  #' @param chr_lookup A lookup table for chromosome names
  #' 
  #' @return TRUE if inv1 and inv2 are on different chromosomes, FALSE otherwise
  #' ___________________________________________________________________________
  
  chr1 <- chr_lookup[inv1]
  chr2 <- chr_lookup[inv2]
  return(!is.na(chr1) && !is.na(chr2) && chr1 != chr2)
}