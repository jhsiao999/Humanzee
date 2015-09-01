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
  
  n_permute <- length(permute_labels)
  order_datatypes <- match(datatypes, c("rna", "ribo", "protein"))
  exclude_datatypes <- c("rna", "ribo", "protein")[setdiff(c(1:3), order_datatypes)]

  eset_sub <-   eset_full[ ,eset_full$seqData != exclude_datatypes
                           & eset_full$species != "rhesus"]
  
  null_interact <- lapply(1:n_permute, function(each_null) {
    emat1 <- exprs(eset_sub)[ permute_labels[[each_null]][,order_datatypes[1]], 
                              eset_sub$seqData == datatypes[1] ]
    emat2 <- exprs(eset_sub)[ permute_labels[[each_null]][,order_datatypes[2]], 
                              eset_sub$seqData == datatypes[2] ]
    emat_per_null <- cbind(emat1, emat2)
    eset_per_null <- ExpressionSet(assayData = as.matrix(emat_per_null))
    phenoData(eset_per_null) <- phenoData(eset_sub)
    featureData(eset_per_null) <- featureData(eset_sub)
    return(interact2way(eset_per_null) )
  })
  null_interact
}
