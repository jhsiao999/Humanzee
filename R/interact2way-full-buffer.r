#' Original version of interaction model for comparing inter-speciees divergence
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
#' interact2way_full_buffer()

interact2way_full_buffer <- function(eSet) {
# eSet = eSetRRP.log2.Qmed[,eSetRRP.log2.Qmed$species!="rhesus" & eSetRRP.log2.Qmed$seqData!="ribo"]
#  source(file.path(codedir,"DEtesting.r"))
  fNames = featureNames(eSet)
  
  require(nlme)

  gls.res = lapply(1:length(fNames), function(i) {
    mat1 = eSet[featureNames(eSet)==fNames[i],]
    mat1aov.temp = data.frame(cov=c(t(exprs(mat1))),seqData=mat1$seqData,
                species=mat1$species,celline=mat1$celline,
                dummy.var=as.integer(mat1$seqData==head(rev(unique(mat1$seqData)),1)))
    mat1aov.temp$dummy.var = as.factor(mat1aov.temp$dummy.var)
    mat1aov.temp$species = as.factor(mat1aov.temp$species)
    
    tmp0 = tryCatch(fit0 <- gls(cov~species+dummy.var,
                          weights=varIdent(form=~1|dummy.var),data=mat1aov.temp,
                          na.action=na.omit),
                    error = function(c) list("error", conditionMessage(c)) )

    if (length(tmp0)==2) {
      res <- data.frame(LR.pval=NA,coef.pval=NA,dummy1=NA,dummy0=NA)
    } else {
      fit1 = tryCatch(foo <- gls(cov~species*dummy.var,
                 weights=varIdent(form=~1|dummy.var),data=mat1aov.temp,na.action=na.omit),
                 error = function(c) list("error", conditionMessage(c)))
      if (length(fit1)==2) {
          res = data.frame(LR.pval=NA,coef.pval=NA,dummy1=NA,dummy0=NA)
        } else {
          aov <- anova(tmp0,fit1)
          res <- data.frame(LR.pval = aov[2,"p-value"],
                            coef.pval = summary(fit1)$tTable[4,4])
          fit <- gls(cov~species*dummy.var-1,
                    weights=varIdent(form=~1|dummy.var),
                    data=mat1aov.temp,na.action=na.omit)
          res = data.frame(res)          
        }
    }
    
    return(res)
  })
  pval = do.call(rbind,gls.res)  

  int.qval = get_qval(pval$LR.pval)
  
  return(data.frame(ENSGID=fNames,int.pval=pval$LR.pval,int.qval=int.qval))
}



# # extract results on the regression coefficient
# interact2way.2 <- function(eSet) {
# #eSet=eSetNormNone_HumanChimp
#   suppressMessages(require(nlme))
#   fNames = featureNames(eSet)
#   nn = dim(exprs(eSet))[2]
#   gls.res = lapply(1:length(fNames), function(i) {
# #i=1
#     mat1 = eSet[featureNames(eSet)==fNames[i],]
#     mat1aov.temp = data.frame(cov=c(t(exprs(mat1))),seqData=mat1$seqData,
#                               species=mat1$species,celline=mat1$celline,
#                               dummy.var=as.integer(mat1$seqData==head(rev(unique(mat1$seqData)),1)))
#     
#       fit1 = gls(cov~species*dummy.var,
#                  weights=varIdent(form=~1|dummy.var),data=mat1aov.temp)
#         
#       res = data.frame(ENSGID=fNames[i],int.coef=fit1$coef[4],int.se=sqrt(fit1$varBeta[4,4]),
#                        int.t=fit1$coef[4]/sqrt(fit1$varBeta[4,4]),int.df=nn-4,
#                        int.pval=2*(1-pt(abs(fit1$coef[4]/sqrt(fit1$varBeta[4,4])),nn-4)) ,
#                        rna.coef=fit1$coef[2],rna.se=sqrt(fit1$varBeta[2,2]),
#                        rna.t=fit1$coef[2]/sqrt(fit1$varBeta[2,2]),rna.df=nn-4,
#                        rna.pval=2*(1-pt(abs(fit1$coef[2]/sqrt(fit1$varBeta[2,2])),nn-4)),
#                        ribo.coef=fit1$coef[2]+fit1$coef[4],
#                        ribo.se=sqrt(fit1$varBeta[2,2]+fit1$varBeta[4,4]+2*fit1$varBeta[2,4]),
#                        ribo.t=(fit1$coef[2]+fit1$coef[4])/sqrt(fit1$varBeta[2,2]+fit1$varBeta[4,4]+2*fit1$varBeta[2,4]),
#                        ribo.df=nn-4,
#                        ribo.pval=2*(1-pt((fit1$coef[2]+fit1$coef[4])/sqrt(fit1$varBeta[2,2]+fit1$varBeta[4,4]+2*fit1$varBeta[2,4]),nn-4))
#                        )
#     return(res)
#   })
#   gls.res = do.call(rbind,gls.res)  
# 
#   suppressMessages(require(fdrtool))  
#   gls.res$int.qval = fdrtool(gls.res$int.pval,statistic="pvalue",verbose=FALSE,plot=FALSE)$qval
#   gls.res$ribo.qval = fdrtool(gls.res$ribo.pval,statistic="pvalue",verbose=FALSE,plot=FALSE)$qval
#   gls.res$rna.qval = fdrtool(gls.res$rna.pval,statistic="pvalue",verbose=FALSE,plot=FALSE)$qval
# 
#   return(gls.res)
# }
# 
# 
# 


#=============================
# # coefficient intepretation
# ribo = 1; rna= 0 (dummy.var)
# chimp = 1; human = 2
# 
# cov = b1 (int) + b2 (human) + b3 (ribo) + b4 (ribo,human)
# 
# rna, chimp = b1
# int 
# mat1aov.temp$cov[1] - fit1$coef[1]
# fit1$resid[1]
# 
# ribo, chimp = b1 + b3
# int + 1*dummy.var 
# mat1aov.temp$cov[12] - fit1$coef[1] - fit1$coef[3]
# fit1$resid[12]
# 
# rna, human = b1 + b2
# int + speciesHuman 
# mat1aov.temp$cov[7] - fit1$coef[1] - fit1$coef[2]
# fit1$resid[7]
# 
# ribo,human = b1 + b2 + b3 + b4
# int + speciesHuman + dummy.var + speciesHuman:dummy.var
# mat1aov.temp[18] - fit1$coef[1] - fit1$coef[2] - fit1$coef[3] - fit1$coef[4]
#
# 
# rna, chimp = b1
# ribo, chimp = b1 + b3
# rna, human = b1 + b2
# ribo,human = b1 + b2 + b3 + b4
#
# contrasts 
# recall the data is on log scale
# human ribo/rna = b3 + b4
# chimp ribo/rna = b3
# human ratio / chimp ratio = b4
#
# rna human = b1+b2
# rna chimp = b1
# rna human/chimp = b2
# 
# ribo human = b1+b2+b3+b4
# ribo chimp = b1+b3
# ribo human/chimp = b2+b4

