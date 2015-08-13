# # Differential expression analysis of measurements
# # Example
# load(file.path(rdadir,"eSetRRP.rda"))
# massRes = testDE(eSetRRPall,seqData="protein",species=c("human","chimp"))$massRes
testDE <- function(eSet,seqData,species) {
# eSet=eSetRRP.Q.log2[which(riboRes.Q$qval<.01),];seqData="rna";species=c("human","chimp")
# eSet=eSetRRP.Q.log2;seqData="rna";species=c("human","chimp")

#  require(BiocParallel)
  require(lmtest)
  require(qvalue)

  submat = eSet[,eSet$seqData==seqData & eSet$species%in% species]
  species = as.vector(submat$species)
  emat = t(exprs(submat))
  
  res = lapply(1:NCOL(emat),function(i) {
# i=1
    tmp=tryCatch(fit <- fit.lm(emat[,i],species),
                    error=function(c) list("error", conditionMessage(c)))
    if (length(tmp)==2) return(rep(NA,4))
    if (length(tmp)==4) return(data.frame(coef=tmp$coef,se=tmp$se,t=tmp$t,pval=tmp$pval))
    })
  res = data.frame(ENSGID=featureData(submat)$ENSGID,do.call(rbind,res))
  
  res$qval = get.qval(res$pval)

  return(list(res=res))
}

fit.lm <- function(exp,species) {
# exp=emat[,3]
  null=lm(exp~1)
  model=lm(exp~species)
  coef = model$coef[2]
  se = summary(model)$coef[2,2]
  t = summary(model)$coef[2,3]
  pval = lrtest(model,null)$Pr[2]
  return(list(coef=coef,se=se,t=t,pval=pval))
}

get.qval <- function(pvalues) {
  require(qvalue)
  if(sum(is.na(pvalues))==0) {
    fdr = qvalue(pvalues)
    qval = fdr$qval
  } else {
    pval.temp = pvalues[which(!is.na(pvalues))]
    fdr = qvalue(pval.temp)
    fdr$pvalues = pvalues
    qval.temp = rep(NA,length(pvalues))
    qval.temp[which(!is.na(pvalues))] = fdr$qvalues
    qval=qval.temp
  }
  return(qval)
}