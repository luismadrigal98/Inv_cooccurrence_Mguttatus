jaccard_tanimoto_condenser <- function(results, metainfo)
{
  #' Condenses the results of the jaccard-tanimoto analysis
  #' 
  #' @param results A list containing the results of the jaccard-tanimoto analysis
  #' @param metainfo A data frame containing the meta information of the inversions
  #' 
  #' @return A data frame containing the condensed results
  #' ___________________________________________________________________________
  
  df <- data.frame(entity_1 = NA,
                   Dosage_1 = NA,
                   ID_1 = NA,
                   entity_2 = NA,
                   Dosage_2 = NA,
                   ID_2 = NA, 
                   chrom_1 = NA, 
                   chrom_2 = NA, 
                   jaccard = NA, 
                   expectation = NA, 
                   p_value = NA)
  
  ref <- which(! is.na(results[[2]]), arr.ind = T)
  
  for (i in 1:nrow(ref))
  {
    row_i <- ref[i, 1]
    col_i <- ref[i, 2]
    
    INV_1 <- rownames(results[[1]])[row_i]
    Dosage_1 <- strsplit(INV_1, "_")[[1]][3]
    ID_1 <- strsplit(INV_1, "_")[[1]][2]
    
    INV_2 <- colnames(results[[1]])[col_i]
    Dosage_2 <- strsplit(INV_2, "_")[[1]][3]
    ID_2 <- strsplit(INV_2, "_")[[1]][2]
    
    new_obs <- data.frame(entity_1 = INV_1,
                          Dosage_1 = Dosage_1,
                          ID_1 = ID_1,
                          entity_2 = INV_2,
                          Dosage_2 = Dosage_2,
                          ID_2 = ID_2,
                          chrom_1 = metainfo[metainfo$INV_ID == ID_1, 'Chr'], 
                          chrom_2 = metainfo[metainfo$INV_ID == ID_2, 'Chr'], 
                          jaccard = results[['statistics']][row_i, col_i], 
                          expectation = results[['expectation']][row_i, col_i], 
                          p_value = results[['pvalues']][row_i, col_i])
    
    df <- rbind2(df, new_obs)
  }
  
  df = df[-1, ]
  
  return(df)
}