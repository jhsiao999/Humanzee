#' Bootstrap cell labels for computing confidence intervals
#'
#' Bootstrap observations (e.g., cells) at each level of the grouping
#' variable (e.g., individual) to compute confidence interval
#' for a statistic that has only 1 measurement per individual. This
#' function was developed for Tung et al. (2006) in comparing
#' adjusted coefficients of variation computed from three subsets of
#' cells, each of which is quantified for a LCL.
#'
#'
#' @param log2counts log2counts matrix of gene by cells
#' @param grouping_vector the grouping vector corresponds to variable of interest
#' @param num_bootstrap number of bootstrap samples
#'
#' @export
#' @examples
#' bootstrap_cells()
#'
bootstrap_cells <- function(log2counts, grouping_vector,
                            number_bootstrap) {
#   log2counts <- molecules_ENSG
#   grouping_vector <- anno_filter$individual
  number_cells <- data.frame(table(anno_filter$individual))

  bootstrap_log2counts <- lapply(1:number_bootstrap, function(ii_boot) {
    # creata a sequence of random numbers
    per_group <- lapply(1:3, function(ii_boot) {
      num_cells <- number_cells$Freq[which(number_cells$Var1 == unique(grouping_vector)[ii_boot]) ]
      ind_log2counts <- log2counts[ , grouping_vector == unique(grouping_vector)[ii_boot]]
      num_cells <- ncol(ind_log2counts)
      bootstrap_data <- ind_log2counts[ , sample(1:num_cells, replace = TRUE)]
      bootstrap_data
    })
    per_group <- do.call(cbind, per_group)
    per_group
  })
  return(bootstrap_log2counts)
}




