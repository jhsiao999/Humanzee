#' Wapper function for running interaction model on null data sts
#'
#' @param eset_full
#' @param datatypes
#' @param permute_labels
#' 
#' @export
#'
#' @examples 
#' permute_interact()

permute_interact <- function(eset_full, datatypes, permute_labels) {
  
  eset_sub <-   eset_full[ ,eset_full$seqData != exclude_datatypes
                           & eset_full$species != "rhesus"]
  
  n_permute <- length(permute_labels)
  order_datatypes <- match(datatypes, c("rna", "ribo", "protein"))
  exclude_datatypes <- c("rna", "ribo", "protein")[setdiff(c(1:3), order_datatypes)]
  
  null_interact <- lapply(1:n_permute, function(each_null) {
    emat1 <- exprs(eset_sub)[ permute_labels[[each_null]], 
                              eset_sub$seqData == datatypes[1] ]
    emat2 <- exprs(eset_sub)[ permute_labels[[each_null]][,2], 
                              eset_sub$seqData == datatypes[2] ]
    emat_per_null <- cbind(emat1, emat2)
    eset_per_null <- ExpressionSet(assayData = as.matrix(emat_per_null),
                                   phenoData = phenoData(eset_sub),
                                   experimentData = experimentData(eset_sub))
    featureData(eset_sub) = featureData(eset_sub)
    return(interact2way(eset_sub) )
  })
  null_interact
}
