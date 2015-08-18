#' Remove unwanted variation using mixed model
#'
#' This function is a modified version of limma gls.series
#' function. We changed to gls.series function to include
#' model residuals that are computed after removing variation
#' due to the random effect component. Otherwise, no changes was
#' made to the gls.series algorithm.
#'
#' @param M Data matrix of genes by samples.
#' @param design Design matrix of the model fixed effects.
#' @param block A vector that uniquely identifies levels of each random
#'          effect.
#' @param correlation Consensus correlation coefficient describing the
#'          average correlation between biological replicates. This is
#'          usually computed via the duplicateCorrelation() function in
#'          the limma package.
#'
#' @keywords batch
#'
#' @export
#'
#' @examples
#' ruv_mixed_model()
#'

ruv_mixed_model <-
  function (M, design = NULL, ndups = 1, spacing = 1,
            block = NULL, correlation = NULL, weights = NULL) {
  M <- as.matrix(M)
  narrays <- ncol(M)
  if (is.null(design))
    design <- matrix(1, narrays, 1)
  design <- as.matrix(design)
  if (nrow(design) != narrays)
    stop("Number of rows of design matrix does not match number of arrays")
  if (is.null(correlation))
    correlation <- duplicateCorrelation(M, design = design,
                    ndups = ndups, spacing = spacing, block = block,
                    weights = weights, ...)$consensus.correlation
  if (!is.null(weights)) {
    weights[is.na(weights)] <- 0
    weights <- asMatrixWeights(weights, dim(M))
    M[weights < 1e-15] <- NA
    weights[weights < 1e-15] <- NA
  }
  nbeta <- ncol(design)
  coef.names <- colnames(design)
  if (is.null(block)) {
    if (ndups < 2) {
      warning("No duplicates: correlation between duplicates set to zero")
      ndups <- 1
      correlation <- 0
    }
    if (is.null(spacing))
      spacing <- 1
    cormatrix <- diag(rep(correlation, len = narrays), nrow = narrays,
                      ncol = narrays) %x% array(1, c(ndups, ndups))
    M <- unwrapdups(M, ndups = ndups, spacing = spacing)
    if (!is.null(weights))
      weights <- unwrapdups(weights, ndups = ndups, spacing = spacing)
    design <- design %x% rep(1, ndups)
    colnames(design) <- coef.names
  }
  else {

    if (ndups > 1) {
      stop("Cannot specify ndups>2 and non-null block argument")
    }
    else {
      ndups <- spacing <- 1
    }
    # Prepare correlation matrix of the randome effects
    block <- as.vector(block)
    if (length(block) != narrays)
      stop("Length of block does not match number of arrays")
    ub <- unique(block)
    nblocks <- length(ub)
    Z <- matrix(block, narrays, nblocks) == matrix(ub, narrays,
                                                   nblocks, byrow = TRUE)
    cormatrix <- Z %*% (correlation * t(Z))
  }
  diag(cormatrix) <- 1
  ngenes <- nrow(M)
  stdev.unscaled <- matrix(NA, ngenes, nbeta, dimnames = list(rownames(M),
                                                              coef.names))
  NoProbeWts <- !any(is.na(M)) && (is.null(weights) || !is.null(attr(weights,
                  "arrayweights")))
  if (NoProbeWts) {
    V <- cormatrix
    if (!is.null(weights)) {
      wrs <- 1/sqrt(weights[1, ])
      V <- wrs * t(wrs * t(V))
    }
    cholV <- chol(V)
    y <- backsolve(cholV, t(M), transpose = TRUE)
    dimnames(y) <- rev(dimnames(M))
    X <- backsolve(cholV, design, transpose = TRUE)
    dimnames(X) <- dimnames(design)
    fit <- lm.fit(X, y)
    if (fit$df.residual > 0) {
      if (is.matrix(fit$effects))
        fit$sigma <- sqrt(colMeans(fit$effects[-(1:fit$rank),
                                               , drop = FALSE]^2))
      else fit$sigma <- sqrt(mean(fit$effects[-(1:fit$rank)]^2))
    }
    else fit$sigma <- rep(NA, ngenes)
    fit$fitted.values <- fit$residuals <- fit$effects <- NULL
    fit$coefficients <- t(fit$coefficients)
    fit$cov.coefficients <- chol2inv(fit$qr$qr, size = fit$qr$rank)
    est <- fit$qr$pivot[1:fit$qr$rank]
    dimnames(fit$cov.coefficients) <- list(coef.names[est],
                                           coef.names[est])
    stdev.unscaled[, est] <- matrix(sqrt(diag(fit$cov.coefficients)),
                                    ngenes, fit$qr$rank, byrow = TRUE)
    fit$stdev.unscaled <- stdev.unscaled
    fit$df.residual <- rep.int(fit$df.residual, ngenes)
    dimnames(fit$stdev.unscaled) <- dimnames(fit$stdev.unscaled) <- dimnames(fit$coefficients)
    fit$pivot <- fit$qr$pivot
    fit$ndups <- ndups
    fit$spacing <- spacing
    fit$block <- block
    fit$correlation <- correlation
  #  return(fit)
  }
  beta <- stdev.unscaled
  sigma <- rep(NA, ngenes)
  df.residual <- rep(0, ngenes)
  residuals <- vector("list", ngenes)
  for (i in 1:ngenes) {
    y <- drop(M[i, ])
    o <- is.finite(y)
    y <- y[o]
    n <- length(y)
    if (n > 0) {
      X <- design[o, , drop = FALSE]
      V <- cormatrix[o, o]
      if (!is.null(weights)) {
        wrs <- 1/sqrt(drop(weights[i, o]))
        V <- wrs * t(wrs * t(V))
      }
      cholV <- chol(V)
      y <- backsolve(cholV, y, transpose = TRUE)
      if (all(X == 0)) {
        df.residual[i] <- n
        sigma[i] <- sqrt(array(1/n, c(1, n)) %*% y^2)
      }
      else {
        X <- backsolve(cholV, X, transpose = TRUE)
        out <- lm.fit(X, y)
        est <- !is.na(out$coefficients)
        beta[i, ] <- out$coefficients
        stdev.unscaled[i, est] <- sqrt(diag(chol2inv(out$qr$qr,
                                                     size = out$rank)))
        df.residual[i] <- out$df.residual
        if (df.residual[i] > 0)
          sigma[i] <- sqrt(array(1/out$df.residual,
                        c(1,n)) %*% out$residuals^2)
        residuals[[i]] <- out$residuals
      }
    }
  }
  cholV <- chol(cormatrix)
  QR <- qr(backsolve(cholV, design, transpose = TRUE))
  cov.coef <- chol2inv(QR$qr, size = QR$rank)
  est <- QR$pivot[1:QR$rank]
  dimnames(cov.coef) <- list(coef.names[est], coef.names[est])

  return(
    list(coefficients = beta, stdev.unscaled = stdev.unscaled,
       sigma = sigma, df.residual = df.residual, ndups = ndups,
       spacing = spacing, block = block, correlation = correlation,
       cov.coefficients = cov.coef, pivot = QR$pivot,
       residuals = do.call(rbind, residuals) ) )
}
