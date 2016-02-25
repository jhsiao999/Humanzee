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
#' @family single-cell
#' @export
#' @examples
#' # tmp <- permute_cells(log2counts = molecules_final_df[1:5,],
#' #                      grouping_vector = anno_filter$individual,
#' #                      number_permute = 2,
#' #                      subset_matrix = molecules_expressed_df[1:5, ])
permute_cells <- function(log2counts,
                          grouping_vector,
                          number_permute,
                          subset_matrix = NULL) {
  if (!is.null(subset_matrix)) make_subset <- TRUE

  if (make_subset) {
    if (!all.equal(dim(log2counts), dim(subset_matrix))) {
      stop("dimension of count matrix does not match dimension of subset matrix")
    }
  }

  # number of cells
  number_cells <- dim(log2counts)[2]

  permuted_log2counts <- lapply(1:number_permute, function(ii_permute) {
    # creata a random sequence of labels
    perm_labels <- sample(number_cells, replace = FALSE)

    # reorder columns (cells) in the data matrix
    # according to perm_labels
    perm_data <- log2counts[ , perm_labels]

    if(make_subset) {
      perm_subset <- subset_matrix[ , perm_labels]
      perm_filtered <- perm_data*perm_subset
      perm_final <- perm_filtered
    } else {
      perm_final <- perm_data
    }

    # relabel the columns with the original individual labels
    # now the cells for individual A in the permuted data
    # can come from indivdual A, B, or C
    colnames(perm_final) <- grouping_vector
    perm_final
  })
  return(permuted_log2counts)
}
