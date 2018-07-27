#===============================================================================
# multinomial_regression.R
#===============================================================================

# A function for building a multinomial regression model on read counts from
# SELEX results




# Imports ======================================================================

#' @import nnet




# Functions ====================================================================

#' @title Preprocess count data
#'
#' @description Prepares count data for multinomial logistic regression
#'
#' @details
#' Converts the input data to matrix format if necessary and ensures that 
#' reference allele counts occupy the first column.
#'
#' @param counts A matrix or data frame with alleles for column names
#' @param ref_allele A character string representing the reference allele
#' @return Matrix containing SELEX count data
#' @seealso \code{\link{selex_multinom}}
preprocess_counts <- function(counts, ref_allele = NULL) {
  counts <- data.matrix(counts)
  if (!is.null(ref_allele)) {
    alleles <- colnames(counts)
    alt_alleles <- alleles[alleles != ref_allele]
    counts <- counts[,c(ref_allele, alt_alleles)]
  }
  counts
}

#' Multinomial logistic regression on count data from SELEX experiments
#'
#' Built around the \code{multinom} function from \code{nnet}
#'
#' @param counts A matrix or data frame with alleles for column names
#' @param ref_allele A character string representing the reference allele
#' @param weights A numeric vector of regression weights for the 5 cycles
#' @param sink A logical, when TRUE messages from \code{multinom} will be hidden
#' @return A multinom fit object
#' @export
selex_multinom <- function(
  counts,
  ref_allele = NULL,
  weights = default_weights,
  sink = TRUE
) {
  counts <- preprocess_counts(counts, ref_allele = ref_allele)
  cycle <- c(0, 1, 2, 3, 4)
  if (sink) {sink("/dev/null")}
  fit <- multinom(counts ~ cycle, weights = weights)
  if (sink) {sink()}
  fit[["counts"]] <- counts
  fit[["ref.allele"]] <- colnames(counts)[[1]]
  fit[["input.weights"]] <- weights
  fit[["coefficients"]] <- summary(fit)[["coefficients"]]
  fit[["standard.errors"]] <- summary(fit)[["standard.errors"]]
  fit[["entropy"]] <- -1 * rowSums(
    fit[["fitted.values"]] * log2(fit[["fitted.values"]])
  )
  fit[["info.content"]] <- 2 - fit[["entropy"]]
  fit
}
