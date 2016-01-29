#' Overlay density curves using base package
#'
#' @param expression_matrix gene-by-sample matrix with rownames consisting of 
#'                          Ensemble gene IDs.
#' @param annotation meta-data of the samples.
#' @param which_gene Ensemble gene ID of the gene to plot.
#' @param xlab x-axis label
#' @param gene_symbols a gene-by-c(Ensemble ID, gene symbol) matrix. Gene symbols can be
#'                     replaced with other gene information of interest.
#' @param labels secondary title right below gene symbol
#' 
#' @family single-cell
#' @export
#'
#' @examples
#' density_overlay()

density_overlay <- function(molecules_ENSG,
                         annotation,
                         which_gene, 
                         xlab, labels = NULL, gene_symbols = NULL) {
    library(scales)
    library(broman)
    crayon <- brocolors("crayon")
    dens <- 
        lapply(1:3, function(per_individual) {
            which_individual <- annotation$individual == unique(annotation$individual)[per_individual]
            density(unlist( molecules_ENSG[ rownames(molecules_ENSG) == which_gene, which_individual] ) )
        })
    xlims <- range(sapply(dens, function(obj) obj$x))
    ylims <- range(sapply(dens, function(obj) obj$y))
    plot(dens[[1]], 
         xlab = xlab, main = "",
         ylab = "Density", axes = F, lwd = 0, xlim = xlims, ylim = ylims)
    polygon(dens[[1]], col = alpha(crayon["Sunset Orange"], .4), border = "grey40")
    polygon(dens[[2]], col = alpha(crayon["Tropical Rain Forest"], .6), border = "grey40")
    polygon(dens[[3]], col = alpha(crayon["Denim"], .3), border = "grey40")
    axis(1); axis(2)
#    mtext(text = labels, side = 3)
#     title(main = with(gene_symbols, 
#                       external_gene_name[which(ensembl_gene_id == which_gene)]) )
}
