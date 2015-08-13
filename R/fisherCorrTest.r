#' Fisher's test for comparing Pearson's correlation coefficients 
#' 
#' This function computes a p-value of the difference between
#' two Pearson's correlation coefficients
#' 
#' @param cor1 Correlation coefficient
#' @param cor2 Correlation coefficient
#' @param N1 Sample size of cor1
#' @param N2 Sample size of cor2
#'
#' @keywords Humanzee
#' 
#' @export
#' 
#' @examples
#' fisherCorrTest()

fisherCorrTest <- function(cor1, cor2, N1, N2) {
  Z1 <- 0.5*log( (1+cor1)/(1-cor1) )
  Z2 <- 0.5*log( (1+cor2)/(1-cor2) )
  SE <- sqrt( 1/(N1-3) + 1/(N2-3) )
  z <- (Z1 - Z2)/SE 
  pval <- 2*(1-pnorm(z))
  list(Z1 = Z1, Z2 = Z2, SE = SE, 
       z = z, pval = pval)
}
