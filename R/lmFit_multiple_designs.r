lmFit <- function (object, design = NULL, ndups = 1, spacing = 1, block = NULL, 
          correlation, weights = NULL, method = "ls", ...) 
{
    y <- getEAWP(object)
    if (is.null(design)) 
        design <- y$design
    if (is.null(design)) 
        design <- matrix(1, ncol(y$exprs), 1)
    else {
        design <- as.matrix(design)
        if (mode(design) != "numeric") 
            stop("design must be a numeric matrix")
        if (nrow(design) != ncol(y$exprs)) 
            stop("row dimension of design doesn't match column dimension of data object")
    }
    ne <- nonEstimable(design)
    if (!is.null(ne)) 
        cat("Coefficients not estimable:", paste(ne, collapse = " "), 
            "\n")
    if (missing(ndups) && !is.null(y$printer$ndups)) 
        ndups <- y$printer$ndups
    if (missing(spacing) && !is.null(y$printer$spacing)) 
        spacing <- y$printer$spacing
    if (missing(weights) && !is.null(y$weights)) 
        weights <- y$weights
    method <- match.arg(method, c("ls", "robust"))
    if (ndups > 1) {
        if (!is.null(y$probes)) 
            y$probes <- uniquegenelist(y$probes, ndups = ndups, 
                                       spacing = spacing)
        if (!is.null(y$Amean)) 
            y$Amean <- rowMeans(unwrapdups(as.matrix(y$Amean), 
                                           ndups = ndups, spacing = spacing), na.rm = TRUE)
    }
    if (method == "robust") 
        fit <- mrlm(y$exprs, design = design, ndups = ndups, 
                    spacing = spacing, weights = weights, ...)
    else if (ndups < 2 && is.null(block)) 
        fit <- lm.series(y$exprs, design = design, ndups = ndups, 
                         spacing = spacing, weights = weights)
    else {
        if (missing(correlation)) 
            stop("the correlation must be set, see duplicateCorrelation")
        fit <- gls.series(y$exprs, design = design, ndups = ndups, 
                          spacing = spacing, block = block, correlation = correlation, 
                          weights = weights, ...)
    }
    if (NCOL(fit$coef) > 1) {
        i <- is.na(fit$coef)
        i <- apply(i[, 1] == i[, -1, drop = FALSE], 1, all)
        n <- sum(!i)
        if (n > 0) 
            warning("Partial NA coefficients for ", n, " probe(s)", 
                    call. = FALSE)
    }
    fit$genes <- y$probes
    fit$Amean <- y$Amean
    fit$method <- method
    fit$design <- design
    new("MArrayLM", fit)
}