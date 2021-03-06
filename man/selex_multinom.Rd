% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/multinomial_regression.R
\name{selex_multinom}
\alias{selex_multinom}
\title{Multinomial logistic regression on count data from SELEX experiments}
\usage{
selex_multinom(counts, ref_allele = NULL, weights = default_weights,
  sink = TRUE)
}
\arguments{
\item{counts}{A matrix or data frame with alleles for column names}

\item{ref_allele}{A character string representing the reference allele}

\item{weights}{A numeric vector of regression weights for the 5 cycles}

\item{sink}{A logical, when TRUE messages from \code{multinom} will be hidden}
}
\value{
A multinom fit object
}
\description{
Built around the \code{multinom} function from \code{nnet}
}
