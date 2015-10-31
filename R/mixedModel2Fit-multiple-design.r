#' Fitting mixed models of one random effect
#' 
#' Adapted from statmod package's mixedModel2Fit function. This function
#' includes the additional feature of fitting a different design matrix
#' for each gene (row) of the input data matrix.
#'
#' 
#' @param design Design matrix.
#' @param block Phentype for the blocking factor (vector).
#' @param yy Gene (feature) by sample matrix.
#'
#' @keywords Humanzee
#' 
#' @export
#' @examples
#' mixedModel2Fit_multiple_design()
#' data_matrix <- pheno_data
#' design <- model.matrix(~ 1+ as.factor(individual) + protein, 
#'                        data = pheno_data)
#' block <- as.factor(pheno_data$biorep)
#' yy <- pheno_data$phos

mixedModel2Fit_multiple_design <- 
    function(design, block, yy) {
        o <- is.finite(yy)
        A <- factor(block[o])
        nobs <- sum(o)
        nblocks <- length(levels(A))
        nbeta <- NCOL(design)
        nafun <- function(e) NA

#        if (nobs > (nbeta + 2) && nblocks > 1 && nblocks < nobs - 1) {
            yy <- yy[o]
            X <- design[o, , drop = FALSE]
            Z <- model.matrix(~0 + A)
            s <- tryCatch(statmod::mixedModel2Fit(yy, X, 
                          Z, only.varcomp = TRUE, maxit = 20)$varcomp, 
                          error = nafun)
            if (!is.na(s[1])) 
                rho <- s[2]/sum(s)
        return(rho)
#        }
    }




