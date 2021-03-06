% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/permute-cells.r
\name{permute_cells}
\alias{permute_cells}
\title{Permute cell labels for computing empirical p-values}
\usage{
permute_cells(log2counts, grouping_vector, number_permute,
  subset_matrix = NULL)
}
\arguments{
\item{log2counts}{log2counts matrix of gene by cells}

\item{grouping_vector}{the grouping vector corresponds to variable of interest}

\item{number_permute}{number of permuted samples}
}
\description{
Permute individual labels to compute empirical p-value
for similarity metrics comparing the observed
individual level coefficients of variation (CV). This
function was developed for Tung et al. (2006) in comparing
adjusted coefficients of variation computed from three subsets of
cells, each of which is quantified for a LCL.
}
\examples{
# tmp <- permute_cells(log2counts = molecules_final_df[1:5,],
#                      grouping_vector = anno_filter$individual,
#                      number_permute = 2,
#                      subset_matrix = molecules_expressed_df[1:5, ])
}
\seealso{
Other single.cell: \code{\link{cell_phase_assign}},
  \code{\link{plot_density_overlay}}
}

