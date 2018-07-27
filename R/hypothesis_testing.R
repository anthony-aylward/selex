#===============================================================================
# hypothesis_testing.R
#===============================================================================

# Hypothesis testing on multinomial models of SELEX count data



# Imports ======================================================================

#' @import parallel




# Functions ====================================================================

#' @title Randomize SELEX count data
#'
#' @description
#' Produces a randomized version of read counts from a SELEXexperiment
#'
#' @param counts A matrix or data frame with alleles for column names
#' @return A matrix of randomized counts
#' @export
#' @seealso \code{\link{sample_coefficients}}
randomize_counts <- function(counts) {
  t(
    sapply(
      rowSums(counts),
      function(n_reads) rmultinom(1, n_reads, colSums(counts))
    )
  )
}

#' Generate a random sample of SELEX multinomial regression coefficients
#'
#' Randomizes the input counts and generates null coefficients
#'
#' @param counts A matrix or data frame with alleles for column names
#' @param weights A numeric vector of regression weights for the 5 cycles
#' @param n An integer giving number of coefficient samples to generate
#' @param cores An integer giving the number of cpu cores to use
#' @return A matrix of coefficients as drawn from the null distribution
#' @export
#' @seealso \code{\link{estimate_standard_errors}}
sample_coefficients <- function(
  counts,
  weights = default_weights,
  n = 1,
  cores = detectCores()
) {
  n_alleles = ncol(counts)
  matrix(
    unlist(
      mclapply(
        lapply(1:n, function(x) randomize_counts(counts)),
        function(c) {
          (
            selex_multinom(c, weights = weights)
            [["coefficients"]]
            [n_alleles:(2*(n_alleles-1))]
          )
        },
        mc.cores = cores
      )
    ),
    nrow = n_alleles - 1
  )
}

#' Estimate standard errors for multinomial regression coefficients
#'
#' Generates a sample of randomized SELEX data and computes "null" coefficients
#'
#' @param counts A matrix or data frame with alleles for column names
#' @param weights A numeric vector of regression weights for the 5 cycles
#' @param n An integer giving number of coefficient samples to generate
#' @param cores An integer giving the number of cpu cores to use
#' @return A numeric vector containing the estimated SE's
#' @seealso \code{\link{estimate_z_scores}}
#' @export
selex_standard_errors <- function(
  counts,
  weights = default_weights,
  n = 100,
  cores = detectCores()
) {
  sample <- sample_coefficients(counts, weights = weights, n = n, cores = cores)
  sapply(1:(ncol(counts)-1), function(row) sd(sample[row,]))
}

#' Estimate Z-scores for SELEX multinomial regression coefficients
#'
#' Compares multinomial logistic regression parameters to estimated SE's
#'
#' @param fit A matrix or data frame with alleles for column names
#' @param estimated_se A numeric vector giving standard errors to use
#' @param n An integer giving number of coefficient samples to generate
#' @param cores An integer giving the number of cpu cores to use
#' @return A numeric vector containing the estimated Z-scores
#' @seealso \code{\link{estimate_pvals}}
#' @export
selex_z_scores <- function(
  fit,
  estimated_se = NULL,
  n = 100,
  cores = detectCores()
) {
  if (is.null(estimated_se)) {
    estimated_se <- selex_standard_errors(
      fit[["counts"]],
      weights = fit[["input.weights"]],
      n = n,
      cores = cores
    )
  }
  n_alleles = ncol(fit[["counts"]])
  setNames(
    fit[["coefficients"]][n_alleles:(2*(n_alleles-1))] / estimated_se,
    fit[["lab"]][2:n_alleles]
  )
}

#' A two-tailed Z-test
#'
#' @param z A Z-score
#' @return A p-value
#' @seealso \code{\link{estimate_pvals}}
two_tailed_z_test <- function(z) {
  (1 - pnorm(abs(z))) * 2
}

#' Estimate p-values for SELEX multinomial regression coefficients
#'
#' Compares multinomial logistic regression parameters to estimated SE's
#'
#' @param fit A matrix or data frame with alleles for column names
#' @param estimated_se A numeric vector giving standard errors to use
#' @param n An integer giving number of coefficient samples to generate
#' @param cores An integer giving the number of cpu cores to use
#' @return A numeric vector containing the estimated p-values
#' @export
selex_pvals <- function(
  fit,
  estimated_se = NULL,
  n = 100,
  cores = detectCores()
) {
  p.adjust(
    two_tailed_z_test(selex_z_scores(fit, estimated_se, n, cores = cores)),
    method = "bonferroni"
  )
}
