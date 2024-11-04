contingency_filter <- function(results, metainfo)
{
  #' Filter the results of the inversion analysis to remove inversions that are 
  #' on the same chromosome or are the same inversion with a different dosage
  #' level.
  #' 
  #' @param results A data frame containing the results of the inversion analysis.
  #' 
  #' @param metainfo A data frame containing the metadata of the inversions.
  #' 
  #' @return A data frame with the filtered results.
  #' ___________________________________________________________________________
  
  results <- results |>
    filter(INV_1 != INV_2) |>
    mutate(ID1 = sapply(strsplit(INV_1, "_"), `[`, 2),
           ID2 = sapply(strsplit(INV_2, "_"), `[`, 2)) |>
    mutate(chrom1 = metainfo[match(ID1, metainfo$INV_ID), "Chr"],
           chrom2 = metainfo[match(ID2, metainfo$INV_ID), "Chr"]) |>
    filter(chrom1 != chrom2)
  
  results$INV_combination <- apply(results[, c("INV_1", "INV_2")], 1, 
                                   function(x) paste(sort(x), collapse = ""))
  
  results <- results[!duplicated(results$INV_combination), ]
  
  return(results)
}