#' Convert Prepare data for visualizing GO terms in a heatmap
#'
#' @param go_list a list that contains summary statistics of the GO terms in GOstats
#'                package. This is used with the GOtest function.
#'
#' @return heatmap_allgos a matrix that contains GO terms, associated annotations,
#'                        and the p-values of enrichment analysis.
#'
#' @export
#'
#' @examples
#' golist <- list(go_riboGTrna = summary(go_riboGTrna$CC, pvalue = 1),
#'                 go_rnaGTribo = summary(go_rnaGTribo$CC, pvalue = 1),
#'                 go_riboGTpro = summary(go_riboGTpro$CC, pvalue = 1),
#'                 go_proGTribo = summary(go_proGTribo$CC, pvalue = 1) )
#' heatmap_cc <- goHeatmapData (golist)
goHeatmapData <- function(golist){
  nlist <- length(golist)
  goIDlist <- lapply( 1:nlist, function(i) {
    golist[[i]][ , 1]
  })
  names(goIDlist) <- names(golist)

  goTermlist <- lapply( 1:nlist, function(i) {
    golist[[i]]$Term
  })

  goset <- unique( unlist(goIDlist) )
  gosetTerm <- unique( unlist(goTermlist) )

  heatmap_allgos <- data.frame(ID = goset)
  heatmap_allgos$Term <- gosetTerm

  pvals <- lapply(1:nlist, function(i) {
    foo <- rep(1-10^(-6), length(gosetTerm))
    foo[ which(goset %in% golist[[i]][ , 1]) ] <- golist[[i]]$Pvalue[ which( golist[[i]][ , 1] %in% goset )]
    return(foo)
  })
  pvals <- do.call(cbind, pvals)
  colnames(pvals) <- names(golist)
  heatmap_allgos <- cbind(heatmap_allgos, pvals)
  colnames(heatmap_allgos)[1] <- colnames(golist[[1]])[1]

  ind <- which( is.na(heatmap_allgos), arr.ind = TRUE)
  heatmap_allgos[ ind ] <- 1
  rownames(heatmap_allgos) <- heatmap_allgos[ , 1]

  iisig <- pvals > .01
  iisig_gos <- which(rowSums(iisig)!=nlist)

  heatmap_allgos <- heatmap_allgos[iisig_gos, ]

  return(heatmap_allgos)
}


