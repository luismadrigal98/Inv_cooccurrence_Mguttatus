p_value_extractor <- function(models, expected)
{
  ## Models is an array, where the third dimension represent the model
  
  p_values_X2 <- foreach(i = 1:dim(models)[3], .combine = 'c') %dopar% 
    {
      chisq.test(models[,,i], p = expected)$p.value
    }
  
  ## Getting the four submatrix of interest
  
  ## Submatrix 1: Both heterozygous contrast (exluding the last row and column)
  
  sub1 <- foreach(i = 1:dim(models)[3]) %dopar% 
    {
      submatrix_extractor(models[,,i], 3, 3)
    }
  
  ## Submatrix 2: Both homozygous contrast (excluding the second row and column)
  
  sub2 <- foreach(i = 1:dim(models)[3]) %dopar% 
    {
      submatrix_extractor(models[,,i], 2, 2)
    }
  
  ## Submatrix 3 and 4: One homozygous and the other hetwrozygous
  
  sub3 <- foreach(i = 1:dim(models)[3]) %dopar% 
    {
      submatrix_extractor(models[,,i], 2, 3)
    }
  
  sub4 <- foreach(i = 1:dim(models)[3]) %dopar% 
    {
      submatrix_extractor(models[,,i], 3, 2)
    }
  
  ## Getting the two vectors per submatrix (each vector is the presence absence)
  
  vectors1 <- foreach(i = 1:dim(models)[3]) %dopar% 
    {
      two_vectors_from_submatrix(sub1[[i]])
    }
  
  vectors2 <- foreach(i = 1:dim(models)[3]) %dopar% 
    {
      two_vectors_from_submatrix(sub2[[i]])
    }
  
  vectors3 <- foreach(i = 1:dim(models)[3]) %dopar% 
    {
      two_vectors_from_submatrix(sub3[[i]])
    }
  
  vectors4 <- foreach(i = 1:dim(models)[3]) %dopar% 
    {
      two_vectors_from_submatrix(sub4[[i]])
    }
  
  ## Getting the p-values for the four contrasts (affinity)
  
  affinity1 <- foreach(i = 1:dim(models)[3], .combine = 'c') %dopar% 
    {
      as.numeric(affinity(matrix(unlist(vectors1[[i]]), nrow = 2, byrow = T), 
                          row.or.col = 'row', datatype = 'binary', 
                          pvalType="midP", sigdigit = 3)$all[, 'p_value'])
    }
  
  affinity2 <- foreach(i = 1:dim(models)[3], .combine = 'c') %dopar% 
    {
      as.numeric(affinity(matrix(unlist(vectors2[[i]]), nrow = 2, byrow = T), 
                          row.or.col = 'row', datatype = 'binary', 
                          pvalType="midP", sigdigit = 3)$all[, 'p_value'])
    }
  
  affinity3 <- foreach(i = 1:dim(models)[3], .combine = 'c') %dopar% 
    {
      as.numeric(affinity(matrix(unlist(vectors3[[i]]), nrow = 2, byrow = T), 
                          row.or.col = 'row', datatype = 'binary', 
                          pvalType="midP", sigdigit = 3)$all[, 'p_value'])
    }
  
  affinity4 <- foreach(i = 1:dim(models)[3], .combine = 'c') %dopar% 
    {
      as.numeric(affinity(matrix(unlist(vectors4[[i]]), nrow = 2, byrow = T), 
                          row.or.col = 'row', datatype = 'binary', 
                          pvalType="midP", sigdigit = 3)$all[, 'p_value'])
    }
  
  p_values_affinity <- list(sub1 = affinity1, 
                            sub2 = affinity2, 
                            sub3 = affinity3, 
                            sub4 = affinity4)
  
  ## Getting the p-values for the four contrasts (jaccard)
  
  jaccard1 <- foreach(i = 1:dim(models)[3], .combine = 'c') %dopar% 
    {
      jaccard.test(vectors1[[i]][[1]], vectors1[[i]][[2]], method = "bootstrap",
                   B = 10000, verbose = F)[['pvalue']]
    }
  
  jaccard2 <- foreach(i = 1:dim(models)[3], .combine = 'c') %dopar% 
    {
      jaccard.test(vectors2[[i]][[1]], vectors2[[i]][[2]], method = "bootstrap",
                   B = 10000, verbose = F)[['pvalue']]
    }
  
  jaccard3 <- foreach(i = 1:dim(models)[3], .combine = 'c') %dopar% 
    {
      jaccard.test(vectors3[[i]][[1]], vectors3[[i]][[2]], method = "bootstrap",
                   B = 10000, verbose = F)[['pvalue']]
    }
  
  jaccard4 <- foreach(i = 1:dim(models)[3], .combine = 'c') %dopar% 
    {
      jaccard.test(vectors4[[i]][[1]], vectors4[[i]][[2]], method = "bootstrap",
                   B = 10000, verbose = F)[['pvalue']]
    }
  
  p_values_jaccard <- list(sub1 = jaccard1, 
                           sub2 = jaccard2, 
                           sub3 = jaccard3, 
                           sub4 = jaccard4)
  
  return(list(X2 = p_values_X2, affinity = p_values_affinity,
              jaccard = p_values_jaccard))
}