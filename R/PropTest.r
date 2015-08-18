#' Two-proportion z-test
#' 
#' This function compares two indepenent sample proportions.
#' The null hypothesis is that the two proportions are equal. 
#' 
#' @param prop1 sample proportion
#' @param prop2 sample proportion
#' @param N1 Sample size of prop1
#' @param N2 Sample size of prop2
#'
#' @keywords Humanzee
#' 
#' @export
#' 
#' @examples
#' PropCorrTest()

PropTest <- function(prop1, prop2, N1, N2) {
  prop_all <- (prop1*N1 + prop2*N2)/(N1 + N2)
  z <- (prop1 - prop2)/ sqrt( prop_all*(1-prop_all)*(1/N1 + 1/N2) )
  pval <- 2*(1-pnorm(z))
  list(z = z, pval = pval)
}
