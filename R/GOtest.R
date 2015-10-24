#' Convert ENSG IDs to ENTREZ IDs
# adapted from http://davetang.org/muse/2010/11/10/gene-ontology-enrichment-analysis/
#'   
#' @param my_ensembl_gene_universe: gene universe
#' @param my_ensembl_gene_test: genes of interest
#' @param pval_cutoff: cutoff when determining significant enrichment GO terms
#' @param ontology: a list of GO ontologies that will be tested in the enrichment anlaysis. 
#'        Default: ontology=c("BP","CC","MF").
#' @param pval_cutoff: cutoff when determining significant enrichment GO terms. 
#'        Default: pval_cutoff = .01.
#' @return hg_over: lists of enrichment test results, each one is a GOHyperGResult 
#'        instance of the GOstats package.
#' 
#' @export
#' 
#' @examples
#' GOtest()

GOtest <- function(my_ensembl_gene_universe,my_ensembl_gene_test,
                   pval_cutoff = .01,ontology=c("BP","CC","MF"),
                   conditional.method = F) {
    require(GO.db)
    require(GOstats)
    require(biomaRt)
    ensembl <- useMart("ensembl",dataset="hsapiens_gene_ensembl")
    # Convert Ensembl IDs to Entrez IDs
    ensembl_to_entrez <- getBM(attributes=c('ensembl_gene_id', 'entrezgene'), 
                                filters = 'ensembl_gene_id', 
                                values = my_ensembl_gene_universe, mart = ensembl)
    ii_entrez_gene <- which(ensembl_to_entrez$ensembl_gene_id %in% my_ensembl_gene_test)
    my_entrez_gene <- ensembl_to_entrez[ii_entrez_gene, ]
    
    ensembl_to_entrez <- unique(ensembl_to_entrez)
    my_entrez_gene <- unique(my_entrez_gene)
    
    # #get some more info on the entrez_gene
    # my_attribute <- c('entrezgene',
    #                   'hgnc_symbol',
    #                   'chromosome_name',
    #                   'start_position',
    #                   'end_position',
    #                   'strand')
    # my_entrez_gene_info  <- getBM(attributes=my_attribute,
    #                         filters = c('entrezgene', 'chromosome_name'),
    #                         values = list(entrezgene=my_entrez_gene$entrezgene, 
    #                                       chromosome_name=my_chr),
    #                         mart = ensembl)
    
    # GENE to GO BP test for over-representation
    hg_over <- lapply(seq_along(ontology), function(i) {
      params <- new("GOHyperGParams",
                    geneIds=my_entrez_gene$entrezgene,
                    universeGeneIds = ensembl_to_entrez$entrezgene,
                    ontology=ontology[i],
                    pvalueCutoff=pval_cutoff,
                    conditional = conditional.method,
                    testDirection="over",
                    annotation="org.Hs.eg.db")
      hg_over_foo <- hyperGTest(params)
      return(hg_over_foo)
    })
    names(hg_over) <- ontology
    list(GO = hg_over, geneIds = ensembl_to_entrez )
}


