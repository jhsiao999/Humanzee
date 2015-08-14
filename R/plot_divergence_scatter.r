#' Scatter plot of per gene inter-species divergence estimates
#' 
#' This function is a wrapper for a scatter plot displaying
#' per gene inter-species divergence estimates
#' with each gene marked for two divergence categories
#' by color and shape
#' 
#' @param xy
#' @param int_points1
#'
#' @keywords Humanzee
#' 
#' @export
#' 
#' @examples
#' plot_divergence_scatter()

plot_divergence_scatter <- function(xy, 
                            ind_points1 = NULL, ind_points2 = NULL, 
                            xlab, ylab,
                            col_bkpoints = "black",
                            cex_bkpoints = .4,
                            lwd_bkpoints = .5,
                            pch_bkpoints = 1,
                            col_points1 = brocolors("crayons")["Red Orange"],
                            col_points2 = brocolors("crayons")["Red Orange"],
                            col_diag = "grey", 
                            cex_points1 = .13, cex_points2 = .13,
                            pch_points1 = 20, pch_points2 = 20,
                            col_hline = brocolors("crayons")["Pine Green"], 
                            col_vline = brocolors("crayons")["Pine Green"], 
                            lwd_hline = .5,
                            lwd_vline = .5,
                            cex.axis = .3, cex.lab = .3,
                            xlim = c(-6, 10), ylim = c(-6, 10) ) {
    require(scales)
    par( mar = c(2.2, 2, 1, .5), mgp = c(.7, .1, 0) )
    plot( xy, pch = "", 
          axes = F, ylim = ylim, xlim = xlim, xlab = "", ylab = "" )

    points(xy, col = col_bkpoints, pch = pch_bkpoints, cex = cex_bkpoints, lwd = .5 )

    if ( !is.null(ind_points2) ) {
      points(xy[ind_points2, ], pch = pch_points2, col = col_points2,
             cex = cex_points2)
    }

    points(xy[ind_points1, ], pch = pch_points1, col = col_points1,
           cex = cex_points1)
    
    if (xlab == "log2 (human/chimp) \n Ribosome occupancy") {
      title( xlab = xlab, cex.lab = cex.lab, line = 1.3)
    } else {
      title( xlab = xlab, cex.lab = cex.lab)
    }

    title( ylab = ylab, cex.lab = cex.lab)

    axis(1, tck = -.06, cex.axis = cex.axis, lwd = .6, at = c(-4, 0, 4, 8)
         , label = c(-4, 0, 4, 8) )
    axis(2, tck = -.03, cex.axis = cex.axis, lwd = .6, at = c(-4, 0, 4, 8)
         , label = c(-4, 0, 4, 8) )
    
    abline(h = 0, lty = 1, col = col_hline, lwd = lwd_hline)
    abline(v = 0, lty = 1, col = col_vline, lwd = lwd_vline)
    abline(a = 0, b = 1, lty = 1, col = col_diag, lwd = .5)
}


