#' Overlay density curves using base package
#'
#' @param expression_matrix gene-by-sample matrix with rownames consisting of 
#'                          Ensemble gene IDs.
#' @param annotation meta-data of the samples.
#' @param which_gene Ensemble gene ID of the gene to plot.
#' @param labels secondary title right below gene symbol
#' @param gene_symbols a gene-by-c(Ensemble ID, gene symbol) matrix. Gene symbols can be
#'                     replaced with other gene information of interest.
#' 
#' @family single-cell
#' @export
#'
#' @examples
#' density_overlay()

plot_density <- function(molecules_ENSG, annotation,
                         individuals, batches = NULL,
                         which_gene, labels, 
                         xlims = NULL, ylims = NULL, gene_symbols) {
  if_present <- which(rownames(molecules_ENSG) == which_gene)
  if(length(if_present) == 0) {
    stop("Gene not present in the data")
  }
  
  library(scales)
  library(broman)
  crayon <- brocolors("crayon")
  
  if (is.null(batches)) {
    individuals <- unique(annotation$individual)
    colors <- c("Sunset Orange", "Tropical Rain Forest", "Denim")
    dens <- lapply(1:3, function(per_individual) {
      which_individual <- annotation$individual == individuals[per_individual]
      density(unlist( molecules_ENSG[ rownames(molecules_ENSG) == which_gene, 
                                      which_individual] ) )
    })

    if (is.null(xlims)) xlims <- range(sapply(dens, function(obj) obj$x))
    if (is.null(ylims)) ylims <- range(sapply(dens, function(obj) obj$y))
    
    plot(dens[[1]], 
         xlab = "log2 gene expression", main = "",
         ylab = "Density", axes = F, lwd = 0, xlim = xlims, ylim = ylims)
    for (i in 1:length(individuals)) {
      polygon(dens[[i]], 
              col = alpha(crayon[colors[i]], .4), 
              border = "grey40")
    }
    axis(1); axis(2)
    mtext(text = labels, side = 3)
    title(main = with(gene_symbols, 
                      external_gene_name[which(ensembl_gene_id == which_gene)]) )
  }
  
  if (!is.null(batches)) {
    
    individuals <- unique(annotation$individual)
    dens <- lapply(1:length(individuals), function(per_individual) {
      which_individual <- annotation$individual == individuals[per_individual]
      annotation_sub <- annotation[which_individual, ]
      molecules_sub <- molecules_ENSG[ , which_individual]
      replicates <- unique(annotation_sub$replicate)
      dens_batch <- lapply(1:length(replicates), function(per_replicate) {
        which_replicate <- annotation_sub$replicate == replicates[per_replicate]
        density(unlist( molecules_sub[ rownames(molecules_ENSG) == which_gene, 
                                       which_replicate] ) )
      })
    })
    
    if (is.null(xlims)) {
      xlims <- range( c( sapply(dens, function(obj_individual) {
        c( sapply(obj_individual, function(obj) {
          range(obj$x)
        }) )
      }) ) )
    }
    if (is.null(ylims)) {
      ylims <- range( c( sapply(dens, function(obj_individual) {
        c( sapply(obj_individual, function(obj) {
          range(obj$y)
        }) )
      }) ) )
    }
    
    colors <- c("Sunset Orange", "Tropical Rain Forest", "Denim")
    for (i in 1:length(dens)) {
      
      if (i == 1) col <- crayon["Sunset Orange"]
      if (i == 2) col <- crayon["Tropical Rain Forest"]
      if (i == 3) col <- crayon["Denim"]
 
      plot(dens[[i]][[1]], 
           xlab = "log2 gene expression", main = "",
           ylab = "Density", axes = F, lwd = 0, xlim = xlims, ylim = ylims)
      for (j in 1:length(dens[[i]])) {
        polygon(dens[[i]][[j]], 
                col = alpha(col, .4), 
                border = "grey40")
      }
    }
      # first individual
#       plot(dens[[1]][[1]], 
#            xlab = "log2 gene expression", main = "",
#            ylab = "Density", axes = F, lwd = 0, xlim = xlims, ylim = ylims)
#       for (j in 1:length(dens[[1]]) ) {
#         polygon(dens[[1]][[j]], col = alpha(crayon[colors[1]], .4), 
#                 border = "grey40") }
#       # second individuadl
#       plot(dens[[2]][[1]], 
#            xlab = "log2 gene expression", main = "",
#            ylab = "Density", axes = F, lwd = 0, xlim = xlims, ylim = ylims)
#       for (j in 1:length(dens[[2]]) ) {
#         polygon(dens[[2]][[j]], col = alpha(crayon[colors[2]], .4), 
#                 border = "grey40") }
#       # third individual
#       plot(dens[[3]][[1]], 
#            xlab = "log2 gene expression", main = "",
#            ylab = "Density", axes = F, lwd = 0, xlim = xlims, ylim = ylims)
#       for (j in 1:length(dens[[3]]) ) {
#         polygon(dens[[3]][[j]], col = alpha(crayon[colors[3]], .4), 
#                 border = "grey40") }
#       
      
        
    axis(1); axis(2)
    mtext(text = labels, side = 3)
    title(main = with(gene_symbols, 
                      external_gene_name[which(ensembl_gene_id == which_gene)]) )
  }
}


