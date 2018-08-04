#===============================================================================
# hypothesis_testing.R
#===============================================================================

# Hypothesis testing on multinomial models of SELEX count data



# Imports ======================================================================

#' @import parallel
#' @import R.utils




# Functions ====================================================================

#' @title Randomize SELEX count data
#'
#' @description
#' Produces a randomized version of read counts from a SELEXexperiment
#'
#' @param counts A matrix or data frame with alleles for column names
#' @return A matrix of randomized counts
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
#' @param timeout Numeric. The time limit for fitting a model during null
#'   sampling, in seconds. (workaround for divergent edge cases)
#' @return A matrix of coefficients as drawn from the null distribution
#' @seealso \code{\link{selex_standard_errors}}
sample_coefficients <- function(
  counts,
  weights = default_weights,
  n = 1,
  cores = detectCores(),
  timeout = 1
) {
  n_alleles = ncol(counts)
  matrix(
    unlist(
      Filter(
        Negate(is.null),
        mclapply(
          lapply(1:n, function(x) randomize_counts(counts)),
          function(c) {
            (
              evalWithTimeout(
                selex_multinom(c, weights = weights),
                timeout = timeout,
                onTimeout = "warning"
              )
              [["coefficients"]]
              [n_alleles:(2*(n_alleles-1))]
            )
          },
          mc.cores = cores
        )
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
#' @param timeout Numeric. The time limit for fitting a model during null
#'   sampling, in seconds. (workaround for divergent edge cases)
#' @return A numeric vector containing the estimated SE's
#' @export
#' @seealso \code{\link{selex_z_scores}}
selex_standard_errors <- function(
  counts,
  weights = default_weights,
  n = 100,
  cores = detectCores(),
  timeout = 1
) {
  sample <- sample_coefficients(
    counts,
    weights = weights,
    n = n,
    cores = cores,
    timeout = timeout
  )
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
#' @param timeout Numeric. The time limit for fitting a model during null
#'   sampling, in seconds. (workaround for divergent edge cases)
#' @return A numeric vector containing the estimated Z-scores
#' @export
#' @seealso \code{\link{selex_pvals}}
selex_z_scores <- function(
  fit,
  estimated_se = NULL,
  n = 100,
  cores = detectCores(),
  timeout = 1
) {
  if (is.null(estimated_se)) {
    estimated_se <- selex_standard_errors(
      fit[["counts"]],
      weights = fit[["input.weights"]],
      n = n,
      cores = cores,
      timeout = timeout
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
#' @export
#' @seealso \code{\link{selex_pvals}}
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
#' @param timeout Numeric. The time limit for fitting a model during null
#'   sampling, in seconds. (workaround for divergent edge cases)
#' @export
#' @return A numeric vector containing the estimated p-values
selex_pvals <- function(
  fit,
  estimated_se = NULL,
  n = 100,
  cores = detectCores(),
  timeout = 1
) {
  p.adjust(
    two_tailed_z_test(
      selex_z_scores(fit, estimated_se, n, cores = cores, timeout = timeout)
    ),
    method = "bonferroni"
  )
}
