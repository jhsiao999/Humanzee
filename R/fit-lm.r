#' Differential expression analysis using likelihood
#' ratio test statistic
#' 
#' This function tests for differential gene expression
#' using a nested linear model framework.
#' 
#' @param exp
#' @param species
#'
#' @keywords Humanzee
#' 
#' @export
#' 
#' @examples
#' fit_lm()

fit_lm <- function(exp,species) {

  null <- lm(exp ~ 1)
  model <- lm(exp ~ species)
  coef <- model$coef[2]
  se <- summary(model)$coef[2,2]
  t <- summary(model)$coef[2,3]
  pval <- lrtest(model,null)$Pr[2]
  
  return(list(coef = coef, se = se, 
              t = t, pval = pval) ) 
}
