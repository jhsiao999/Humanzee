#' Compute coefficient of variation 
#' 
#' Compute CV across cells belong to each level of the grouping variable. The 
#' function was developed for Tung et al. (2016) comparing cell-to-cell 
#' heterogeneity between individuals.
#''
#'
#' @param log2counts log2 count matrix of gene by cell.
#' @param grouping_vector Compute per gene CV for cells belonging to each level of
#'        the grouping vector.
#'
#' @export
#' @examples
#' compute_cv()
#'
compute_cv <- function(log2counts, grouping_vector) {

  group_cv <- lapply( unique(grouping_vector), function(per_group) {
    # Convert log2cpm to counts
    counts_per_group <- 2^log2counts[ , grouping_vector == per_group ]
    mean_per_gene <- apply(counts_per_group, 1, mean, na.rm = TRUE)
    sd_per_gene <- apply(counts_per_group, 1, sd, na.rm = TRUE)
    cv_per_gene <- data.frame(mean = mean_per_gene,
                              sd = sd_per_gene,
                              cv = sd_per_gene/mean_per_gene,
                              group = rep(per_group, dim(counts_per_group)[1]) )
    rownames(cv_per_gene) <- rownames(counts_per_group)

    return(cv_per_gene)
  })
  names(group_cv) <- unique(grouping_vector)
  group_cv
}
