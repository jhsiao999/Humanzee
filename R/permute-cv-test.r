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
permute_cv_test <- function(log2counts, grouping_vector, anno, number_permute,
                              output_rda = FALSE,
                              do_parallel = FALSE,
                              number_cores = NULL) {
  require(Humanzee)
  require(matrixStats)

  permuted_data <- Humanzee::permute_cells(log2counts,
                                            grouping_vector,
                                            number_permute)

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

    # Standardize CVs for each indivdual's CV mean and variance across genes
    df_norm <- sweep(df_perm, MARGIN = 2, STATS = colMeans(as.matrix(df_perm)), FUN = "-")
    df_norm <- sweep(df_perm, MARGIN = 2, STATS = sqrt(colVars(as.matrix(df_perm))), FUN = "/")
    return(df_norm)
  })
  }
  # no parallelization
  if(do_parallel == TRUE) {
    require(doParallel)
    registerDoParallel(cores = number_cores)
    
    permuted_cv_adj <- foreach(ind_data = 1:number_permute) %dopar% {
      perm_cv <- Humanzee::compute_cv(log2counts = per_data,
                                      grouping_vector = grouping_vector)

      perm_cv_adj <- Humanzee::normalize_cv(group_cv = perm_cv,
                                           log2counts = per_data,
                                           anno = anno)

      df_perm <- cbind(perm_cv_adj[[1]]$log10cv2_adj,
                       perm_cv_adj[[2]]$log10cv2_adj,
                       perm_cv_adj[[3]]$log10cv2_adj)

      # Standardize CVs for each indivdual's CV mean and variance across genes
      df_norm <- sweep(df_perm, MARGIN = 2, STATS = colMeans(as.matrix(df_perm)), FUN = "-")
      df_norm <- sweep(df_perm, MARGIN = 2, STATS = sqrt(colVars(as.matrix(df_perm))), FUN = "/")
      return(df_norm)
    }
  }


permuted_distance <- lapply(permuted_cv_adj, function(per_data) {
  squared_dev <- rowSums( ( per_data - rowMedians(as.matrix(per_data)) )^2 )
  abs_dev <- rowSums(abs( per_data - rowMedians(as.matrix(per_data)) ))
  list(squared_dev = squared_dev,
       abs_dev = abs_dev)
})

if (output_rda == TRUE) {
save(permuted_data, permuted_cv_adj,
     file = "permuted-data.rda")
save(permuted_distance,
     file = "permuted-distance.rda")
}
if (output_rda == FALSE) {
  return(list(permuted_data = permuted_data,
              permuted_cv_adj = permuted_cv_adj,
              permuted_distance = permuted_distance))
}

}

