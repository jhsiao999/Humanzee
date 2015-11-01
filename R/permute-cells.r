#' Permute cell labels for computing empirical p-values
#'
#' Permute individual labels to compute empirical p-value
#' for similarity metrics comparing the observed
#' individual level coefficients of variation (CV). This
#' function was developed for Tung et al. (2006) in comparing
#' adjusted coefficients of variation computed from three subsets of
#' cells, each of which is quantified for a LCL.
#'
#' @param log2counts log2counts matrix of gene by cells
#' @param grouping_vector the grouping vector corresponds to variable of interest
#' @param number_permute number of permuted samples
#'
#' @export
#' @examples
#' permuate_cells()
permute_cells <- function(log2counts, grouping_vector, number_permute) {
#     log2counts <- molecules_ENSG
#     grouping_vector <- anno_qc_filter$individual
#     num_permute <- 10
  number_cells <- dim(log2counts)[2]

  permuted_log2counts <- lapply(1:number_permute, function(ii_permute) {
    # creata a random sequence of labels
    perm_labels <- sample(number_cells, replace = FALSE)

    # reorder columns (cells) in the data matrix
    # according to perm_labels
    perm_data <- log2counts[ , perm_labels]

    # relabel the columns with the original individual labels
    # now the cells for individual A in the permuted data
    # can come from indivdual A, B, or C
    colnames(perm_data) <- grouping_vector
    perm_data
  })
  return(permuted_log2counts)
}
