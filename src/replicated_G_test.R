replicated_G_test <- function(df_replicates, degrees_of_fred_replicates,
                              Data_p_a, correction = "williams")
{
  #' This function will perform a replicated independence G-test for the
  #' contingency table of the replicates. The function will return the
  #' p-value of the test and general information associated to the total, pooled,
  #' and heterogeneity tests.
  #' 
  #' @param df_replicates A data frame with the replicates for each contrast 
  #' inversion (each row is a contrast). First two columns defined the two 
  #' inversions in the contrast. The rest of the columns integrate the G statistic
  #' for each line.
  #' 
  #' @param degrees_of_fred_replicates degrees of freedom for each individual test.
  #' 
  #' @param Data_p_a A list with the matrices of presence/absence for each inversion.
  #' 
  #' @param correction A character string with the name of the correction to be
  #' applied. The options are "williams", "yates", and "none". Default is "williams".
  #' 
  #' @return Expanded dataframe with total, pooled, and het G values, the df for
  #' those, and the p_values associated.
  #' ___________________________________________________________________________
  
  cases <- nrow(df_replicates)
  tables <- sum(grepl("G_", colnames(df_replicates)))
  
  # Obtaining the G_total and the df_total
  df_replicates <- df_replicates %>%
    rowwise() %>%
    mutate(G_total = sum(c_across(starts_with("G_"))),
           df_total = degrees_of_fred_replicates * tables) %>%
    ungroup()
  
  # Getting the p_value for G_total
  df_replicates <- df_replicates %>%
    mutate(p_value_total = 1 - pchisq(q = G_total, df = df_total))
  
  # Obtaining the G_pooled and the df_pooled
  G_pool <- vector(length = nrow(df_replicates), mode = "numeric")
  df_pool <- vector(length = nrow(df_replicates), mode = "numeric")
  
  for (i in 1:nrow(df_replicates))
  {
    inv1 <- as.character(df_replicates[i, 1])
    inv2 <- as.character(df_replicates[i, 2])
    
    # Building the pooled table
    pooled_vector_1 <- vector(length = 0, mode = "numeric")
    pooled_vector_2 <- vector(length = 0, mode = "numeric")
    
    for (line in names(Data_p_a))
    {
      pooled_vector_1 <- c(pooled_vector_1, Data_p_a[[line]][inv1, ])
      pooled_vector_2 <- c(pooled_vector_2, Data_p_a[[line]][inv2, ])
    }
    
    pooled_table <- table(pooled_vector_1, pooled_vector_2)
    
    df_pool[i] <- (nrow(pooled_table) - 1) * (ncol(pooled_table) - 1)
    
    # Performing the G-test
    G_pool[i] <- GTest(pooled_table, correct = correction)$statistic
  }
  
  # Adding the G_pooled and the df_pooled to the data frame
  df_replicates <- df_replicates %>%
    mutate(G_pooled = G_pool,
           df_pooled = df_pool)
  
  # Getting the p_value for G_pooled
  df_replicates <- df_replicates %>%
    mutate(p_value_pooled = 1 - pchisq(q = G_pooled, df = df_pooled))
  
  # Getting the G_het, df_het, and p_value_het
  df_replicates <- df_replicates %>%
    mutate(G_het = G_total - G_pooled,
           df_het = df_total - df_pooled,
           p_value_het = 1 - pchisq(q = G_het, df = df_het))
  
  return(df_replicates)
}