#' @export

voom_NA <- function (counts, design = NULL, lib.size = NULL, normalize.method = "none", 
          plot = FALSE, span = 0.5, ...) 
{
  out <- list()
  if (is(counts, "DGEList")) {
    out$genes <- counts$genes
    out$targets <- counts$samples
    if (is.null(design) && diff(range(as.numeric(counts$sample$group))) > 
        0) 
      design <- model.matrix(~group, data = counts$samples)
    if (is.null(lib.size)) 
      lib.size <- with(counts$samples, lib.size * norm.factors)
    counts <- counts$counts
  }
  else {
    isExpressionSet <- suppressPackageStartupMessages(is(counts, 
                                                         "ExpressionSet"))
    if (isExpressionSet) {
      if (length(Biobase::fData(counts))) 
        out$genes <- Biobase::fData(counts)
      if (length(Biobase::pData(counts))) 
        out$targets <- Biobase::pData(counts)
      counts <- Biobase::exprs(counts)
    }
    else {
      counts <- as.matrix(counts)
    }
  }
  if (is.null(design)) {
    design <- matrix(1, ncol(counts), 1)
    rownames(design) <- colnames(counts)
    colnames(design) <- "GrandMean"
  }
  if (is.null(lib.size)) 
    lib.size <- colSums(counts)
  y <- t(log2(t(counts + 0.5)/(lib.size + 1) * 1e+06))
  y <- normalizeBetweenArrays(y, method = normalize.method)
  fit <- lmFit(y, design, ...)
  if (is.null(fit$Amean)) 
    fit$Amean <- rowMeans(y, na.rm = TRUE)
  sx <- fit$Amean + mean(log2(lib.size + 1)) - log2(1e+06)
  sy <- sqrt(fit$sigma)
  allzero <- rowSums(counts, na.rm = TRUE) == 0
  if (any(allzero)) {
    sx <- sx[!allzero]
    sy <- sy[!allzero]
  }
  l <- lowess(sx, sy, f = span)
  if (plot) {
    plot(sx, sy, xlab = "log2( count size + 0.5 )", ylab = "Sqrt( standard deviation )", 
         pch = 16, cex = 0.25)
    title("voom: Mean-variance trend")
    lines(l, col = "red")
  }
  f <- approxfun(l, rule = 2)
  if (fit$rank < ncol(design)) {
    j <- fit$pivot[1:fit$rank]
    fitted.values <- fit$coef[, j, drop = FALSE] %*% t(fit$design[, 
                                                                  j, drop = FALSE])
  }
  else {
    fitted.values <- fit$coef %*% t(fit$design)
  }
  fitted.cpm <- 2^fitted.values
  fitted.count <- 1e-06 * t(t(fitted.cpm) * (lib.size + 1))
  fitted.logcount <- log2(fitted.count)
  w <- 1/f(fitted.logcount)^4
  dim(w) <- dim(fitted.logcount)
  out$E <- y
  out$weights <- w
  out$design <- design
  if (is.null(out$targets)) 
    out$targets <- data.frame(lib.size = lib.size)
  else out$targets$lib.size <- lib.size
  new("EList", out)
}