#===============================================================================
# contrast.R
#===============================================================================

#' Extract a contrast from a SELEX multinomial regression
#'
#' The contrast for a provided ref/alt allele pair
#'
#' @param fit The multinom fit object from selex_multinom
#' @param ref_allele Character string indicating the reference allele
#' @param alt_allele Character string indicating the alternate allele
#' @return Intercept and coefficient for the contrast
contrast <- function(fit, ref_allele, alt_allele) {
  if (ref_allele == fit[["ref.allele"]]) {
    fit[["coefficients"]][alt_allele,]
  } else if (alt_allele == fit[["ref.allele"]]) {
    -1 * fit[["coefficients"]][ref_allele,]
  } else {
    fit[["coefficients"]][alt_allele,] - fit[["coefficients"]][ref_allele,]
  }
}