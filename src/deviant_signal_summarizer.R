deviant_signal_summarizer <- function(result_df, 
                                      jaccard_threshold = 0, 
                                      return_details = FALSE) {
  #' Summarize Deviant Signals in Inversion Interaction Analysis
  #'
  #' @description
  #' This function analyzes and summarizes deviant signals detected in inversion 
  #' interaction data. It categorizes interactions based on dosage combinations 
  #' (Het-Het, Hom-Hom, Hom-Het) and the nature of their epistatic effects 
  #' (positive or negative).
  #'
  #' @param result_df A data frame containing the analysis results. Must have columns:
  #'        INV_1, INV_2, and Jaccard
  #' @param jaccard_threshold Numeric value to determine positive/negative deviation 
  #'        (default: 0)
  #' @param return_details Logical; if TRUE, returns detailed information about each
  #'        interaction (default: FALSE)
  #'
  #' @return If return_details = FALSE, returns a summary data frame with counts.
  #'         If return_details = TRUE, returns a list containing both summary and 
  #'         detailed classification of each interaction.
  #'
  #' @examples
  #' # Basic usage
  #' summary <- deviant_signal_summarizer(result_df)
  #' 
  #' # With custom Jaccard threshold and detailed output
  #' detailed_results <- deviant_signal_summarizer(result_df, 
  #'                                              jaccard_threshold = 0.1,
  #'                                              return_details = TRUE)
  
  # Input validation
  if (!is.data.frame(result_df)) {
    stop("Input must be a data frame")
  }
  
  required_cols <- c("INV_1", "INV_2", "Jaccard")
  missing_cols <- setdiff(required_cols, names(result_df))
  if (length(missing_cols) > 0) {
    stop("Missing required columns: ", paste(missing_cols, collapse = ", "))
  }
  
  if (!is.numeric(jaccard_threshold)) {
    stop("jaccard_threshold must be numeric")
  }
  
  # Create result template
  result <- data.frame(
    submatrix = rep(1:4, 2),
    description = rep(c('Het_Het', 'Hom_Hom', 'Hom_Het_1', 'Hom_Het_2'), 2),
    deviant = rep(c('Negative', 'Positive'), each = 4),
    count = rep(0, 8)
  )
  
  # Initialize detailed results if requested
  if (return_details) {
    detailed_results <- data.frame(
      INV_1 = result_df$INV_1,
      INV_2 = result_df$INV_2,
      Jaccard = result_df$Jaccard,
      Category = character(nrow(result_df)),
      Effect = character(nrow(result_df)),
      stringsAsFactors = FALSE
    )
  }
  
  # Function to safely extract dosage
  extract_dosage <- function(inv_id) {
    parts <- tryCatch({
      strsplit(inv_id, '_')[[1]]
    }, error = function(e) {
      warning("Invalid inversion ID format: ", inv_id)
      return(NA)
    })
    if (length(parts) >= 3) {
      return(parts[3])
    } else {
      warning("Invalid inversion ID format: ", inv_id)
      return(NA)
    }
  }
  
  # Extract dosage information
  Inv1_dosage <- sapply(result_df$INV_1, extract_dosage)
  Inv2_dosage <- sapply(result_df$INV_2, extract_dosage)
  
  # Process each interaction
  for (i in seq_along(Inv1_dosage)) {
    if (is.na(Inv1_dosage[i]) || is.na(Inv2_dosage[i])) {
      next
    }
    
    # Determine category and effect
    category <- if (Inv1_dosage[i] == '1' && Inv2_dosage[i] == '1') {
      "Het_Het"
    } else if (Inv1_dosage[i] == '2' && Inv2_dosage[i] == '2') {
      "Hom_Hom"
    } else if (Inv1_dosage[i] == '1' && Inv2_dosage[i] == '2') {
      "Hom_Het_1"
    } else if (Inv1_dosage[i] == '2' && Inv2_dosage[i] == '1') {
      "Hom_Het_2"
    } else {
      warning("Unexpected dosage combination: ", 
              Inv1_dosage[i], "-", Inv2_dosage[i])
      next
    }
    
    effect <- if (result_df$Jaccard[i] < jaccard_threshold) "Negative" else "Positive"
    
    # Update counts
    result[result$description == category & 
             result$deviant == effect, "count"] <- 
      result[result$description == category & 
               result$deviant == effect, "count"] + 1
    
    # Update detailed results if requested
    if (return_details) {
      detailed_results$Category[i] <- category
      detailed_results$Effect[i] <- effect
    }
  }
  
  # Add summary statistics
  result$percentage <- round(result$count / sum(result$count) * 100, 2)
  
  # Return results
  if (return_details) {
    return(list(
      summary = result,
      details = detailed_results,
      stats = list(
        total_interactions = sum(result$count),
        positive_ratio = sum(result$count[result$deviant == "Positive"]) / 
          sum(result$count),
        negative_ratio = sum(result$count[result$deviant == "Negative"]) / 
          sum(result$count)
      )
    ))
  } else {
    return(result)
  }
}