#' Compute voom weights using customized log2 counts
#'
#' @param log2counts counts on log2 scale 
#' @param lib.size Library size. 
#' @param design Experimental design of the data. Required to be an R 
#'                design.matrix object 
#' @param is.cpm if the data is CPM normalized.
#' 
#' @export
#'
#' @examples
#' voomWeightsCustomized()
#'
voomWeightsCustomized <- function(log2counts, lib.size = NULL, design, is.cpm = FALSE) {
  if (is.null(design)) {
    design <- model.matrix(~ 1 )
  }

  if (is.cpm == TRUE) {
    fit <- lmFit(log2counts, design)
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
  }

  if (is.cpm == FALSE) {
    fit <- lmFit(log2counts, design)
    xx <- fit$Amean 
    yy <- sqrt(fit$sigma)
    l <- lowess(xx, yy, f = .5)
    f <- approxfun(l, rule = 2)
    fitted.values <- fit$coef %*% t(fit$design)
    fitted.logcount <- log2(fitted.values)
    w <- 1/f(fitted.logcount)^4
    dim(w) <- dim(fitted.logcount)
  }

  return(w)
}
