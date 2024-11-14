process_submatrix_tests <- function(submatrices, 
                                    test_type = c("affinity", "jaccard")) 
{
  #' Process submatrices for affinity or Jaccard tests
  #' 
  #' @param submatrices A list of submatrices
  #' @param test_type The type of test to perform
  #' 
  #' @return A numeric vector of p-values
  #' ___________________________________________________________________________
  
  test_type <- match.arg(test_type)
  
  sapply(submatrices, function(mat) {
    vectors <- two_vectors_from_submatrix(mat)
    if (test_type == "affinity") {
      as.numeric(
        affinity(matrix(unlist(vectors), nrow = 2, byrow = TRUE),
                 row.or.col = 'row', datatype = 'binary',
                 pvalType = "midP", sigdigit = 3)$all[, 'p_value']
      )
    } else {
      jaccard.test(vectors[[1]], vectors[[2]],
                   method = "bootstrap", B = 10000, verbose = FALSE)$pvalue
    }
  })
}