#' Adjust coefficients of variation
#'
#' @param group_cv CVs per batch computed use compute_cv().
#' @param log2counts log2 count matrix of gene by cell.
#' 
#' @export
#' @examples
#' normalize_cv()
#'
normalize_cv <- function(group_cv, log2counts, anno) {
  library(zoo)
  # Compute a data-wide coefficient of variation on counts.
  data_cv <- apply(2^log2counts, 1, sd)/apply(2^log2counts, 1, mean)

  # Order genes by mean expression levels
  order_gene <- order(apply(2^log2counts, 1, mean))

  # Rolling medians of log10 squared CV by mean expression levels
  # Avoid warning message introduced by NA in the rollapply results,
  # which are introduced by coercion
  roll_medians <- suppressWarnings(
    rollapply( log10(data_cv^2)[order_gene],
               width = 50, by = 25,
               FUN = median, fill = list("extend", "extend", "NA") )
  )
  ii_na <- which( is.na(roll_medians) )
  roll_medians[ii_na] <- median( log10(data_cv^2)[order_gene][ii_na] )
  names(roll_medians) <- rownames(molecules_ENSG)[order_gene]

  # Order rolling medians according to the count matrix
  reorder_gene <- match(rownames(log2counts), names(roll_medians) )
  roll_medians <- roll_medians[ reorder_gene ]
  stopifnot( all.equal(names(roll_medians), rownames(log2counts) ) )

  group_cv_adj <- lapply(1:length(group_cv), function(ii_batch) {
    # Adjusted coefficient of variation on log10 scale
    log10cv2_adj <- log10(group_cv[[ii_batch]]$cv^2) - roll_medians
    # combine the adjusted cv with the unadjusted cv
    data.frame(group_cv[[ii_batch]],
               log10cv2_adj = log10cv2_adj )
  })
  names(group_cv_adj) <- names(group_cv)
  group_cv_adj
}
