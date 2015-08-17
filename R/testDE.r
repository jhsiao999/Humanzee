#' Inter-species divergence likelihood ratio test
#' 
#' This function tests for divergence between two species
#' in one molecular phenotype at a time. 
#' 
#' @param eSet
#' @param seqData 
#' @param species
#'
#' @keywords Humanzee
#' 
#' @export
#' 
#' @examples
#' load(file.path(rdadir,"eSetRRP.rda"))
#' massRes = testDE(eSetRRPall,seqData="protein",species=c("human","chimp"))$massRes
testDE <- function(eSet,seqData,species) {

  require(lmtest)

  submat = eSet[,eSet$seqData==seqData & eSet$species%in% species]
  species = as.vector(submat$species)
  emat = t(exprs(submat))
  
  res = lapply(1:NCOL(emat),function(i) {
    tmp=tryCatch(fit <- fit_lm(emat[,i],species),
                    error=function(c) list("error", conditionMessage(c)))
    if (length(tmp)==2) return(rep(NA,4))
    if (length(tmp)==4) return(data.frame(coef=tmp$coef,se=tmp$se,t=tmp$t,pval=tmp$pval))
    })
  res = data.frame(ENSGID=featureData(submat)$ENSGID,do.call(rbind,res))
  
  res$qval = get_qval(res$pval)

  return(list(res=res))
}


