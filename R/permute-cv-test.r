#' Permute cell labels for computing empirical p-values
#'
#' Permute individual labels to compute empirical p-value
#' for similarity metrics comparing the observed
#' individual level coefficients of variation (CV). This
#' function was developed for Tung et al. (2006) in comparing
#' adjusted coefficients of variation computed from three subsets of
#' cells, each of which is quantified for a LCL.
#'
#' @param log2counts log2counts matrix of gene by cells.
#' @param grouping_vector the grouping vector corresponds to variable of interest.
#' @param anno Annotation matrix (cell-by-phenotype).
#' @param number_permute number of permuted samples.
#'
#' @export
#' @examples
#' permute_cv_test()
#'
permute_cv_test <- function(log2counts, 
                            grouping_vector, 
                            anno, n
                            umber_permute,
                            subset_matrix = NULL,
                            output_rda = FALSE,
                            do_parallel = FALSE,
                            number_cores = NULL) {
  require(matrixStats)

  if (!is.null(subset_matrix)) make_subset <- TRUE

  if (make_subset) {
    if (!all.equal(dim(log2counts), dim(subset_matrix))) {
      stop("dimension of count matrix does not match dimension of subset matrix")
    }
  }

  permuted_data <- Humanzee::permute_cells(log2counts,
                                            grouping_vector,
                                            number_permute,
                                            subset_matrix)

  # Compute adjusted CV for the permuted data
  # no parallelization
  if(do_parallel == FALSE) {
  permuted_cv_adj <- lapply(permuted_data, function(per_data) {
    perm_cv <- Humanzee::compute_cv(log2counts = per_data,
                                    grouping_vector = grouping_vector)

    perm_cv_adj <- Humanzee::normalize_cv(group_cv = perm_cv,
                                          log2counts = per_data,
                                          anno = anno)

    df_perm <- cbind(perm_cv_adj[[1]]$log10cv2_adj,
                     perm_cv_adj[[2]]$log10cv2_adj,
                     perm_cv_adj[[3]]$log10cv2_adj)

    return(df_perm)
  })
  }
  if(do_parallel == TRUE) {
    require(doParallel)
    registerDoParallel(cores = number_cores)

    permuted_cv_adj <- foreach(ind_data = 1:number_permute) %dopar% {

      per_data <- permuted_data[[ind_data]]
      perm_cv <- Humanzee::compute_cv(log2counts = per_data,
                                      grouping_vector = grouping_vector)

      perm_cv_adj <- Humanzee::normalize_cv(group_cv = perm_cv,
                                           log2counts = per_data,
                                           anno = anno)

      df_perm <- cbind(perm_cv_adj[[1]]$log10cv2_adj,
                       perm_cv_adj[[2]]$log10cv2_adj,
                       perm_cv_adj[[3]]$log10cv2_adj)

      return(df_perm)
    }
  }
rm(permuted_data)

permuted_distance <- lapply(permuted_cv_adj, function(per_data) {
  mad <- rowMedians(abs( per_data - rowMedians(as.matrix(per_data)) ))
  list(mad = mad)
})
rm(permuted_cv_adj)

if (output_rda == TRUE) {
save(permuted_distance,
     file = "permuted-distance.rda")
}
if (output_rda == FALSE) {
  return(list(permuted_data = permuted_data,
              permuted_cv_adj = permuted_cv_adj,
              permuted_distance = permuted_distance$mad))
}

}

