#' Interaction model for comparing inter-species divergence
#' 
#' This function is a wrapper for an interaction model testing
#' inter-species differences in divergence between molecular 
#' phenotypes under a nested linear model framework. We use
#' likelihood ratio statistic to determine statistical significance
#' of the difference between inter-species divergence. 
#' 
#' @param eSet data set contains the two molecular phenotypes
#'
#' @keywords Humanzee
#' 
#' @export
#' @examples
#' interact2way()
#' eSet <- eSetRRP.RP.Q.log2[ ,eSetRRP.RP.Q.log2$seqData!="protein" & 
#'                             eSetRRP.RP.Q.log2$species!="rhesus"]

interact2way <- function(eSet, coef = FALSE) {

  fNames <- featureNames(eSet)
  
  require(nlme)
  
  gls.res <- lapply(1:length(fNames), function(i) {
  
    mat1 <- eSet[ featureNames(eSet)==fNames[i], ]
    mat1aov.temp <- data.frame( cov=c(t(exprs(mat1))), seqData=mat1$seqData,
                              species = mat1$species, celline=mat1$celline )
    mat1aov.temp$species <- as.factor(mat1aov.temp$species)
    mat1aov.temp$seqData <- as.factor(mat1aov.temp$seqData)
    
    # Fit two-way main effect model
    fit_null_try <- tryCatch(
                    fit_null <- gls( cov ~ species + seqData,
                                       weights = varIdent(form=~1|seqData), 
                                       data=mat1aov.temp,
                                       na.action = na.omit ),
                    condition = function(c) c   
                    )
    
    if( inherits( fit_null_try, "condition") ) {

        res <- data.frame(LRatio = NA, LR_pval = NA)
  
    } else {
        # Fit interaction model  
        fit_interact_try <- tryCatch(
                            fit_interact <- gls( cov ~ species * seqData,
                                                  weights = varIdent(form = ~1| seqData),
                                                  data = mat1aov.temp, 
                                                  na.action = na.omit ),
                            condition = function(c) c   
                            )
        if( inherits( fit_interact_try, "condition") ) {
            res <- data.frame(LRatio = NA, LR_pval = NA)
        } else {
        # Likelihood ratio test
            aov <- anova(fit_null_try, fit_interact_try)
            res <- data.frame( LRatio = aov[2, "L.Ratio"],
                               LR_pval = aov[2, "p-value"])
            
            if(coef==TRUE) {
              coef_table <- summary(fit_interact_try)$tTable
              coef_res <- c(coef_table[2, ], coef_table[3, ], coef_table[4, ])
              names(coef_res) <- c(paste(rownames(coef_table)[2], colnames(coef_table), sep = "_"),
                                   paste(rownames(coef_table)[3], colnames(coef_table), sep = "_"),
                                   paste(rownames(coef_table)[4], colnames(coef_table), sep = "_") )
              res <- cbind(res, coef_res)
            }
          } 
    
    return(res) 
    }
    
  })
  LR_res <- do.call(rbind, gls.res)  
  
  int.qval <- get_qval(LR_res$LR_pval)

  return(data.frame(ENSGID=fNames, 
                    LRatio = LR_res$LRatio,
                    int.pval = LR_res$LR_pval,
                    int.qval = int.qval))
}
