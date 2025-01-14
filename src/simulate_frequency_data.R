simulate_frequency_data <- function(n_samples, freq1, freq2, association_strength) 
{
  #' Simulate genotype data for two inversions with specified frequencies and 
  #' association strength
  #'
  #' @param n_samples Integer. Number of samples to generate.
  #' @param freq1 Numeric. Frequency of the first inversion (0 to 1).
  #' @param freq2 Numeric. Frequency of the second inversion (0 to 1).
  #' @param association_strength Numeric. Strength of association between inversions (-1 to 1).
  #'
  #' @return Data frame with two columns representing genotypes of the two inversions.
  #' @examples
  #' simulate_frequency_data(100, 0.3, 0.5, 0.7)
  
  # Validate inputs
  if (!is.numeric(n_samples) || n_samples <= 0 || n_samples != round(n_samples)) {
    stop("n_samples must be a positive integer.")
  }
  if (!is.numeric(freq1) || freq1 < 0 || freq1 > 1) {
    stop("freq1 must be a numeric value between 0 and 1.")
  }
  if (!is.numeric(freq2) || freq2 < 0 || freq2 > 1) {
    stop("freq2 must be a numeric value between 0 and 1.")
  }
  if (!is.numeric(association_strength) || association_strength < -1 || association_strength > 1) {
    stop("association_strength must be a numeric value between -1 and 1.")
  }
  
  # Base independent probabilities (inv is the alternative allele, q)
  prob1 <- c((1 - freq1)^2, 2 * freq1 * (1 - freq1), freq1^2)  # HWE proportions
  prob2 <- c((1 - freq2)^2, 2 * freq2 * (1 - freq2), freq2^2)
  
  # Generate independent uniform random variables
  u <- runif(n_samples)
  v <- runif(n_samples)
  
  # Introduce correlation using a copula
  cop <- normalCopula(param = association_strength, dim = 2, dispstr = "un")
  uv <- cCopula(cbind(u, v), copula = cop, inverse = TRUE)
  
  # Convert uniform random variables to genotypes
  geno1 <- cut(uv[, 1], breaks = c(0, cumsum(prob1)), labels = 0:2, 
               include.lowest = TRUE)
  geno2 <- cut(uv[, 2], breaks = c(0, cumsum(prob2)), labels = 0:2, 
               include.lowest = TRUE)
  geno1 <- as.numeric(as.character(geno1))
  geno2 <- as.numeric(as.character(geno2))
  
  data.frame(inv1 = geno1, inv2 = geno2)
}