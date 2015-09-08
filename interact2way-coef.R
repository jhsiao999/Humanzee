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
#' 
#' @examples
#' interact2way()

eSet <- eSetRRP.RP.Q.log2[,eSetRRP.RP.Q.log2$seqData!="protein" & eSetRRP.RP.Q.log2$species!="rhesus"]

interact2way <- function(eSet) {

  fNames <- featureNames(eSet)
  
  require(nlme)
  
  gls.res <- lapply(1:length(fNames), function(i) {
  
    gls.res <- lapply(1:10, function(i) {
      
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
                               LR_pval = aov[2, "p-value"] )
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

# interact2way_full_buffer <- function(eSet) {
# # eSet = eSetRRP.log2.Qmed[,eSetRRP.log2.Qmed$species!="rhesus" & eSetRRP.log2.Qmed$seqData!="ribo"]
# #  source(file.path(codedir,"DEtesting.r"))
#   fNames = featureNames(eSet)
#   
#   require(nlme)
# #  require(parallel)
#   require(contrast)
# #   gls.res = lapply(1:10, function(i) {
#   gls.res = lapply(1:length(fNames), function(i) {
#     mat1 = eSet[featureNames(eSet)==fNames[i],]
#     mat1aov.temp = data.frame(cov=c(t(exprs(mat1))),seqData=mat1$seqData,
#                 species=mat1$species,celline=mat1$celline,
#                 dummy.var=as.integer(mat1$seqData==head(rev(unique(mat1$seqData)),1)))
#     mat1aov.temp$dummy.var = as.factor(mat1aov.temp$dummy.var)
#     mat1aov.temp$species = as.factor(mat1aov.temp$species)
#     
#     tmp0 = tryCatch(fit0 <- gls(cov~species+dummy.var,
#                           weights=varIdent(form=~1|dummy.var),data=mat1aov.temp,
#                           na.action=na.omit),
#                     error = function(c) list("error", conditionMessage(c)) )
# 
#     if (length(tmp0)==2) {
#       res <- data.frame(LR.pval=NA,coef.pval=NA,dummy1=NA,dummy0=NA)
#     } else {
#       fit1 = tryCatch(foo <- gls(cov~species*dummy.var,
#                  weights=varIdent(form=~1|dummy.var),data=mat1aov.temp,na.action=na.omit),
#                  error = function(c) list("error", conditionMessage(c)))
#       if (length(fit1)==2) {
#           res = data.frame(LR.pval=NA,coef.pval=NA,dummy1=NA,dummy0=NA)
#         } else {
#           aov <- anova(tmp0,fit1)
#           res <- data.frame(LR.pval = aov[2,"p-value"],
#                             coef.pval = summary(fit1)$tTable[4,4])
#           fit <- gls(cov~species*dummy.var-1,
#                     weights=varIdent(form=~1|dummy.var),
#                     data=mat1aov.temp,na.action=na.omit)
#           contrast.res1 <- contrast(fit,a=list(dummy.var="1",species="human"),
#                                  b=list(dummy.var="1",species="chimp"))
#           contrast.res0 <- contrast(fit,a=list(dummy.var="0",species="human"),
#                                  b=list(dummy.var="0",species="chimp"))
#           # dummy.var = 0 for protein and dummy.var = 1 for ribo
#           # b1,b2,b3,b4 below correspond to the coef. in the model matrix of fit
#           # pro_chimp (b3-b4) + ribo_human (b2-b4)- pro_human (b3-(b3-b4)) - ribo_chimp (b1-(b3-b4))= 0
#           # equivalent to: -b1+b2+2b3-4b4=0
# #           contrast.buff=anova(fit,L=c(-1,1,2,-4))
# #           buff.dir = mean(mat1aov.temp$cov[6:10],na.rm=TRUE)-mean(mat1aov.temp$cov[1:5],
# #                   na.rm=TRUE)-(mean(mat1aov.temp$cov[16:20],na.rm=TRUE)-mean(mat1aov.temp$cov[11:15],na.rm=TRUE))
# #           if (buff.dir > 0) {  buff.pval= contrast.buff$p/2 }
# #           if (buff.dir < 0) {  buff.pval= 1-contrast.buff$p/2 }
# #           res = data.frame(res,dummy1=contrast.res1$Pvalue,
# #                            dummy0=contrast.res0$Pvalue,
# #                            buff=buff.pval)          
#           res = data.frame(res,dummy1=contrast.res1$Pvalue,
#                            dummy0=contrast.res0$Pvalue,
#                            int.SE=summary(fit1)$tTable[4,2])          
#         }
#     }
#     
#     return(res)
#   })
#   pval = do.call(rbind,gls.res)  
# 
#   int.qval = get.qval(pval$LR.pval)
#   dummy0.qval = get.qval(pval$dummy0)
#   dummy1.qval = get.qval(pval$dummy1)
# #   buff.qval = get.qval(pval$buff)
#   
# #   return(data.frame(ENSGID=fNames,int.pval=pval$LR.pval,int.qval=int.qval,
# #                     dummy0.qval=dummy0.qval,dummy1.qval=dummy1.qval,
# #                     buff.pval=pval$buff,buff.qval=buff.qval))
#   return(data.frame(ENSGID=fNames,int.pval=pval$LR.pval,int.qval=int.qval,
#                   dummy0.qval=dummy0.qval,dummy1.qval=dummy1.qval,
#                   dummy0.pval=pval$dummy0,dummy1.pval=pval$dummy1,
#                   int.SE=pval$int.SE))
# }
# 


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