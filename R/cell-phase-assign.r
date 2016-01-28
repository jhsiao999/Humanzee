#' Assign cell-phase using an ad hoc approach described 
#' in Macosco et al. (2015) (http://www.ncbi.nlm.nih.gov/pubmed/26000488)
#'
#' This approach was develped for single-cell RNA sequencing data (scRNA-seq)
#' collected using Drop-seq protocol.
#' 
#' @param cell_cycle_gene a list that is comprised of cell-cycle genes (Ensemble IDs)
#' @param expression matrix gene-by-cell matrix. Currently, the method can be applied
#'                          to both count data and log-transformed counts.
#' 
#' @return score_final_normed Phase-specific score matrix. Scores are standardized 
#'                            to z-scores across cells and across phases. 
#' 
#' @family single-cell
#' 
#' @export
#' @examples
#' cell_phase_assign()
#'
cell_phase_assign <- function(cell_cycle_genes, expression_matrix) {
    
    ## Output a gene-by-phase score matrix
    cell_phase_score <- sapply(cell_cycle_genes, function(xx){
        
        ## Extract phase-specific genes
        molecules_phase <- expression_matrix[rownames(expression_matrix) %in% unlist(xx) ,]
        
        ## Computes average expression across cells
        ## Then append the mean expression vector to the last row
        combined_matrix <- rbind(molecules_phase, 
                                 average = apply(molecules_phase,2,mean) )
        
        ## Compute correlation matrix between the genes
        cor_matrix <- cor(t(combined_matrix))
        
        ## Extract the vector that contains correlation
        ## Between the mean gene expression vector and
        ## Every single genes included in the gene set
        cor_vector <- cor_matrix[ ,dim(cor_matrix)[1]]
        
        ## Select genes which mean gene expression level correlates
        ## with every other gene at r >= 0.3 
        molecules_phase_restricted <- molecules_phase[rownames(molecules_phase) %in% names(cor_vector[cor_vector >= 0.3]),]
        
        ## Output the phase specific scores of each cell 
        ## which is the mean of phase-specific gene expression levels 
        apply(molecules_phase_restricted, 2, mean)
    })
    
    ## Two-step normalization (by row and then by column)
    ## For each cell, compute score mean and standard deviation 
    ## across the phases
    row_mean <- apply(cell_phase_score, 1, mean)
    row_sd   <- apply(cell_phase_score, 1, sd)
    score_row_normed <- do.call(rbind, 
                                lapply(1:dim(cell_phase_score)[1], function(i) {
                                    (cell_phase_score[i,] - row_mean[i])/row_sd[i]
                                })  )
    
    ## For every phase, compute score mean and standard deviation
    ## across cells
    col_mean <- apply(score_row_normed, 2, mean)
    col_sd   <- apply(score_row_normed, 2, sd)
    score_final_normed <- do.call(cbind, 
                                  lapply(1:dim(score_row_normed)[2], function(i) {
                                      (score_row_normed[, i] - col_mean[i])/col_sd[i]
                                  })
    )
    return(score_final_normed)
}