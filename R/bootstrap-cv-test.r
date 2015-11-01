#' Bootstrapping to estimate confidence interval for
#' individual coefficients of variation
#'
#' For each gene, resample within each individual and compute a
#' confidence interval of each individual coefficient of variation
#' based on the boostrapped samples. The method was developed for
#' Tung et al. (2016) for comparing cell-to-cell heterogeneity between
#' individuals
#'
#' @param group_cv CVs per batch computed use compute_cv().
#' @param log2counts log2 count matrix of gene by cell.
#' @param number_bootstrap Number of boostrapped samples.
#'
#' @export
#' @examples
#' bootstrap_cv_test()
#'
bootstrap_cv_test <- function(log2counts, grouping_vector, anno, number_bootstrap,
                              output_rda = FALSE,
                              do_parallel = FALSE,
                              number_cores = NULL) {
  require(Humanzee)
  require(matrixStats)

  bootstrap_data <- bootstrap_cells(log2counts,
                                    grouping_vector,
                                    number_bootstrap)

  # Compute adjusted CV for the bootstrapped data
  # no parallelization
  if(do_parallel == FALSE) {
  bootstrap_cv_adj <- lapply(bootstrap_data, function(per_data) {
    boot_cv <- Humanzee::compute_cv(log2counts = per_data,
                          grouping_vector = grouping_vector)

    boot_cv_adj <- Humanzee::normalize_cv(group_cv = boot_cv,
                               log2counts = per_data,
                               anno = anno)

    df_boot <- cbind(boot_cv_adj[[1]]$log10cv2_adj,
                     boot_cv_adj[[2]]$log10cv2_adj,
                     boot_cv_adj[[3]]$log10cv2_adj)

    # Standardize CVs for each indivdual's CV mean and variance across genes
    df_norm <- sweep(df_boot, MARGIN = 2, STATS = colMeans(as.matrix(df_boot)), FUN = "-")
    df_norm <- sweep(df_boot, MARGIN = 2, STATS = sqrt(colVars(as.matrix(df_boot))), FUN = "/")
    return(df_norm)
  })
  }
  # no parallelization
  if(do_parallel == TRUE) {
    require(doParallel)
    registerDoParallel(cores = number_cores)
    bootstrap_cv_adj <- foreach(ind_data = 1:number_bootstrap) %dopar% {
      per_data <- bootstrap_data[[ind_data]]
      boot_cv <- compute_cv(log2counts = per_data,
                            grouping_vector = grouping_vector)

      boot_cv_adj <- normalize_cv(group_cv = boot_cv,
                                  log2counts = per_data,
                                  anno = anno)

      df_boot <- cbind(boot_cv_adj[[1]]$log10cv2_adj,
                       boot_cv_adj[[2]]$log10cv2_adj,
                       boot_cv_adj[[3]]$log10cv2_adj)

      # Standardize CVs for each indivdual's CV mean and variance across genes
      df_norm <- sweep(df_boot, MARGIN = 2, STATS = colMeans(as.matrix(df_boot)), FUN = "-")
      df_norm <- sweep(df_boot, MARGIN = 2, STATS = sqrt(colVars(as.matrix(df_boot))), FUN = "/")
      return(df_norm)
    }
  }

bootstrap_distance <- lapply(bootstrap_cv_adj, function(per_data) {
  squared_dev <- rowSums( ( per_data - rowMedians(as.matrix(per_data)) )^2 )
  abs_dev <- rowSums(abs( per_data - rowMedians(as.matrix(per_data)) ))
  list(squared_dev = squared_dev,
       abs_dev = abs_dev)
})

#rm(bootstrap_data)
#rm(bootstrap_cv_adj)

if (output_rda == TRUE) {
save(bootstrap_data, bootstrap_cv_adj,
     file = "bootstrap-data.rda")
save(bootstrap_distance,
     file = "bootstrap-distance.rda")
}
if (output_rda == FALSE) {
  return(list(bootstrap_data = bootstrap_data,
              bootstrap_cv_adj = bootstrap_cv_adj,
              bootstrap_distance = bootstrap_distance))
}

}

