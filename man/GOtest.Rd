% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/GOtest.R
\name{GOtest}
\alias{GOtest}
\title{Convert ENSG IDs to ENTREZ IDs
adapted from http://davetang.org/muse/2010/11/10/gene-ontology-enrichment-analysis/}
\usage{
GOtest(my_ensembl_gene_universe, my_ensembl_gene_test, pval_cutoff = 0.01,
  ontology = c("BP", "CC", "MF"), conditional.method = F)
}
\arguments{
\item{my_ensembl_gene_universe}{Gene universe}

\item{my_ensembl_gene_test}{Genes of interest}

\item{pval_cutoff}{cutoff when determining significant enrichment GO terms}

\item{ontology}{A list of GO ontologies that will be tested in the enrichment 
anlaysis. Default: ontology=c("BP","CC","MF").}

\item{pval_cutoff}{Cutoff when determining significant enrichment GO terms. 
Default: pval_cutoff = .01.}
}
\value{
hg_over Lists of enrichment test results, each one is a GOHyperGResult 
       instance of the GOstats package.
}
\description{
Convert ENSG IDs to ENTREZ IDs
adapted from http://davetang.org/muse/2010/11/10/gene-ontology-enrichment-analysis/
}
\examples{
GOtest()
}

