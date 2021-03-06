% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/voomWeightsCustomized.r
\name{voomWeightsCustomized}
\alias{voomWeightsCustomized}
\title{Compute voom weights using customized log2 counts}
\usage{
voomWeightsCustomized(log2counts, lib.size = NULL, design, is.cpm = FALSE)
}
\arguments{
\item{log2counts}{counts on log2 scale}

\item{lib.size}{Library size.}

\item{design}{Experimental design of the data. Required to be an R 
design.matrix object}

\item{is.cpm}{if the data is CPM normalized.}
}
\description{
Compute voom weights using customized log2 counts
}
\examples{
voomWeightsCustomized()

}

