#' Compute voom weights using customized log2-cpm values. 
#'
#' @param log2cpm  log2cpm values. 
#' @param lib.size Library size. 
#' @param design Experimental design of the data. Required to be an R 
#'                design.matrix object 
#' @export
#'
#' @examples
#' 
#' voomWeightsCustomized()
#'
voomWeightsCustomized <- function(log2cpm, lib.size, design) {
  if (is.null(design)) {
    design <- model.matrix(~ 1 )
  }
  fit <- lmFit(log2cpm, design)
  xx <- fit$Amean + mean(log2(lib.size + 1)) - log2(1e+06)
  yy <- sqrt(fit$sigma)
  l <- lowess(xx, yy, f = .5)
  f <- approxfun(l, rule = 2)
  fitted.values <- fit$coef %*% t(fit$design)
  fitted.cpm <- 2^fitted.values
  fitted.count <- 1e-06 * t(t(fitted.cpm) * (lib.size + 1))
  fitted.logcount <- log2(fitted.count)
  w <- 1/f(fitted.logcount)^4
  dim(w) <- dim(fitted.logcount)
  return(w)
}
