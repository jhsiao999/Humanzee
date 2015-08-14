#' Multiple hypothesis testings correction
#' 
#' This function uses qvalue package and computes false
#' discovery rate for each gene. 
#' 
#' @param pvalues
#'
#' @keywords Humanzee
#' 
#' @export
#' 
#' @examples
#' get_qval()

get_qval <- function(pvalues) {

  require(qvalue)

  if( sum(is.na(pvalues))==0 ) {

  fdr <- qvalue(pvalues)
  qval <- fdr$qval

  } else {

  pval.temp <- pvalues[which(!is.na(pvalues))]
  fdr <- qvalue(pval.temp)
  fdr$pvalues <- pvalues
  qval.temp <- rep(NA,length(pvalues))
  qval.temp[which(!is.na(pvalues))] <- fdr$qvalues
  qval <- qval.temp

  }

  return(qval)
}