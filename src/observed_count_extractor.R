observed_count_extractor <- function(dosage_list, line, Inv_ID)
{
  #' This function will take as an input a list of dosage matrices (each matrix)
  #' refers to a particular line. Then, using the line and Inv_ID, the dosage 
  #' level per plant is going to be extracted.
  #' 
  #' @Arguments: dosage_list: list of dosage matrices
  #'             line: line of the plant (sublist)
  #'             Inv_ID: Inversion ID (row ID in the dosage matrix)
  #' @Retunrs: A vector with the dosage levels of the plants
  #' ___________________________________________________________________________
  
  obs <- dosage_list[[line]][Inv_ID, ]
  
  RR <- sum(obs == 0) # Homozygous for the reference allele
  RA <- sum(obs == 1) # Heterozygous
  AA <- sum(obs == 2) # Homozygous for the alternative allele (inversion)
  
  return(data.frame(Inv_ID = Inv_ID, Line = line,
                    RR = RR, RA = RA, AA = AA))
}