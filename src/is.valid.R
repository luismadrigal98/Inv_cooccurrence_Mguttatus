is.valid <- function(inv1, inv2, metadata, fields)
{
  #' Auxiliary function to determine if a given contrast between inversions is
  #' valid. I this study, valid conctrasts are those between inversions in different
  #' chromosomes. Also, constrasts between different dosages of the same inversion
  #' are prohibited (cathed by the previous filter too) - One individual cannot 
  #' be heterozygous and homozygous at the same time for an inversion, for example.
  #' 
  #'  @param inv1 A character string with the name of the first inversion. If
  #'  with_dosage is TRUE, the name should be in the format "Inv_X_Y", where X is
  #'  the inversion number and Y is the dosage level.
  #'  
  #'  @param inv2 A character string with the name of the second inversion. If
  #'  with_dosage is TRUE, the name should be in the format "Inv_X_Y", where X is
  #'  the inversion number and Y is the dosage level.
  #'  
  #'  @param metadata A data frame with the information about the inversions and
  #'  genes. It should have the columns "INV_ID" and "Chr". You can refer to these
  #'  columns as indexes or strings.
  #'  
  #'  @param fields A character vector with the names of the columns in metadata
  #'  that should be used to identify the inversions. The first element should be
  #'  the name of the column with the inversion ID and the second element should be
  #'  the name of the column with the chromosome. Indexes are also supported.
  #'  
  #'  @return A logical value indicating if the contrast is valid.
  #'  __________________________________________________________________________
  
  inv1 <- strsplit(inv1, "_")[[1]][2]
  inv2 <- strsplit(inv2, "_")[[1]][2]
  
  chr1 <- ifelse(all(is.numeric(fields)), 
                 metadata[metadata[, fields[1]] == inv1, fields[2]], 
                 metadata[metadata[[fields[1]]] == inv1, fields[2]])
  
  chr2 <- ifelse(all(is.numeric(fields)),
                 metadata[metadata[, fields[1]] == inv2, fields[2]],
                 metadata[metadata[[fields[1]]] == inv2, fields[2]])
  
  if(chr1 != chr2)
  {
    return(TRUE)
  }
  else
  {
    return(FALSE)
  }
}