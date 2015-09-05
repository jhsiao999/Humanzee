# draw_legend = function(color, breaks, legend, ...){
#     height = min(unit(1, "npc"), unit(100, "bigpts"))
    
#     legend_pos = (legend - min(breaks)) / (max(breaks) - min(breaks))
#     legend_pos = height * legend_pos + (unit(1, "npc") - height)
    
#     breaks = (breaks - min(breaks)) / (max(breaks) - min(breaks))
#     breaks = height * breaks + (unit(1, "npc") - height)
    
#     h = breaks[-1] - breaks[-length(breaks)]
    
#     rect = rectGrob(x = 0, y = breaks[-length(breaks)], width = unit(10, "bigpts"), height = h, hjust = 0, vjust = 0, gp = gpar(fill = color, col = "#FFFFFF00"))
#     text = textGrob(names(legend), x = unit(14, "bigpts"), y = legend_pos, hjust = 0, gp = gpar(...))
    
#     res = grobTree(rect, text)
    
#     return(res)
# }


# generate_breaks = function(x, n, center = F){
#     if(center){
#         m = max(abs(c(min(x, na.rm = T), max(x, na.rm = T))))
#         res = seq(-m, m, length.out = n + 1)
#     }
#     else{
#         res = seq(min(x, na.rm = T), max(x, na.rm = T), length.out = n + 1)
#     }
    
#     return(res)
# }

# ##############
# plotDivergenceDensity <- function( dens, main, plot_legend = FALSE, 
#                                    xlims_subset = NULL, ylims_subset = NULL) {
  
#   if (is.null(xlims_subset)) {
#     xlims = range(sapply(dens,function(xx) xx$x))
#   } else {
#     xlims = xlims_subset
#   }

#   if (is.null(ylims_subset)) {
#     ylims = range(sapply(dens,function(xx) xx$y))
#   } else {
#     ylims = ylims_subset
#   }
  
#   plot( dens[[1]], col = "grey20", ylim = ylims, xlim = xlims,
#         main = "", xlab="",lwd = 2, axes = F )
#   title( xlab = "Divergence FC", cex.lab = .6)
#   title(main=main, cex.main = .7)
#   axis(1, tck = -.04, cex.axis = .4)
#   axis(2, tck = -.03, cex.axis = .4)
#   polygon(dens[[1]],col=alpha("grey20",.6), lwd = 2 )
#   lines( dens[[2]], col = "grey70", lty = 2, lwd = 2 )
#   polygon( dens[[2]],col=alpha("grey70", .6), lty = 1, lwd = 2 )
#   lines( dens[[3]], col = "black", lty = 2, lwd = 2 )

#   if (plot_legend==TRUE) {
#     legend("topright", col=c("red","blue","grey20"), lty=1,
#          legend = c("RNA divergence", "Ribo divergence", "Protein divergence"))
#    }
  
# }



# plotDivergenceDensity_v2 <- function( dens, main, plot_legend = FALSE, 
#                                    xlims_subset = NULL, ylims_subset = NULL,
#                                    main.adj = 0, main.line = 0) {
  
#   if (is.null(xlims_subset)) {
#     xlims = range(sapply(dens,function(xx) xx$x))
#   } else {
#     xlims = xlims_subset
#   }

#   if (is.null(ylims_subset)) {
#     ylims = range(sapply(dens,function(xx) xx$y))
#   } else {
#     ylims = ylims_subset
#   }
  
#   require(broman)
#   crayon <- brocolors("crayons")
#   plot( dens[[1]], col = crayon["Scarlet"], ylim = ylims, xlim = xlims,
#         main = "", xlab="",lwd = 2, axes = F, ylab="" )
#   title( xlab = "Divergence effect size", cex.lab = .6)
#   title( ylab = "Density", cex.lab = .6)
#   title(main=main, cex.main = .7, adj = main.adj, line = main.line)
#   axis(1, tck = -.04, cex.axis = .6, at = c(0, 1.5, 3), label = c(0, 1.5, 3))
#   axis(2, tck = -.03, cex.axis = .6, at = c(0, 1, 2, 3), label = c(0, 1, 2, 3))
#   polygon( dens[[2]],col = alpha(crayon["Laser Lemon"], .6), lty = 1, lwd = 1,
#         border = "dodgerblue" )
#   lines( dens[[3]], col = crayon["Pine Green"], lty = 1, lwd = 2 )
# #  polygon(dens[[2]],col= alpha(crayon["Sea Green"],.6), lwd = 1 )
# #  lines( dens[[3]], col = crayon["Lemon Yellow"], lty = 2, lwd = 1 )

#   if (plot_legend==TRUE) {
#     legend("topright", col=c("red","blue","grey20"), lty=1,
#          legend = c("RNA divergence", "Ribo divergence", "Protein divergence"))
#    }
  
# }




# ##############
# require(broman)
# plotDivergenceRatioDensity <- function(dens, 
#                                        bk_col, fn_lcol = crayon["Green"],
#                                        col.main ) {
# crayon <- brocolors("crayon")
#   par( mfrow = c(1, 1), mar = c(2, 1, 1, 2), mgp = c(.9, .2, 0) )
#   plot(dens[[1]], lwd = 1, col = "white", xlab ="", ylab="",
#        main = "", axes = F,sub="",ylim=c(0,.8),
#        xlim = rev(c(-6,8)) )
#   title(xlab="Divergence ratio", cex.lab = .5)
#   axis(1, at= rev(seq(-6,8,by=2) ), labels = rev( seq(-6,8,by=2)), 
#        tck = -.03, cex.axis = .3, lwd = .6, las = 3)
#   axis(4, at=seq(0, 1, 0.2), labels=seq(0,1,0.2), 
#        tck = -.03, cex.axis = .3, lwd = .6)
#   mtext("Proportion", side = 4, cex = .5, line = .5)
#   polygon(dens[[1]], col = bk_col, lwd = 1, border = NA )
#   lines(dens[[2]], col = fn_lcol, lty = 1, lwd = 1.2, xlim=c(-8,8))
# #   legend("topleft",pch=21,col="black",pt.bg=c(col.main,"white"),
# #   legend=c("log2 Ribo/RNA divergence ratio","log2 Protein/RNA divergence ratio"),
# #   box.col="white")
# }


# #############
# plotLegend.sigpro <- function(col.riborna,bg.riborna,
#                               col.rnapro,bg.rnapro,col.both,numsigs) {
#   plot(0,pch="",xlab="",ylab="",axes=F,xlim=c(0,10),ylim=c(0,10))
#   points(1,5,pch=21,bg=bg.rnapro[1],col=col.rnapro[1],cex=2)
#   points(1,4.3,pch=21,bg=bg.rnapro[2],col=col.rnapro[2],cex=2)
#   points(1,3.6,pch=16,col=col.both[1],cex=2)
#   points(1,2.9,pch=16,col=col.both[2],cex=2)
#   text(2,5,label=paste("mRNA > Protein divergence",numsigs[1],"genes"),adj=0,cex=1)
#   text(2,4.3,label=paste("Protein > mRNA divergence",numsigs[2],"genes"),adj=0,cex=1)
#   text(2,3.6,label=paste("mRNA > Ribo and mRNA > Protein",numsigs[3],"genes"),adj=0,cex=1)
#   text(2,2.9,label=paste("Ribo > mRNA & Protein > mRNA",numsigs[4],"genes"),adj=0,cex=1)
# }

# plotLegend.sigribo <- function(col.riborna,bg.riborna,
#                                col.rnapro,bg.rnapro,col.both,numsigs) {
#   plot(0,pch="",xlab="",ylab="",axes=F,xlim=c(0,10),ylim=c(0,10))
#   points(1,5,pch=21,bg=bg.riborna[1],col=col.riborna[1],cex=2)
#   points(1,4.3,pch=21,bg=bg.riborna[2],col=col.riborna[2],cex=2)
#   points(1,3.6,pch=16,col=col.both[1],cex=2)
#   points(1,2.9,pch=16,col=col.both[2],cex=2)
#   text(2,5,label=paste("mRNA > Ribo divergence",numsigs[1],"genes"),adj=0,cex=1)
#   text(2,4.3,label=paste("Ribo > mRNA divergence",numsigs[2],"genes"),adj=0,cex=1)
#   text(2,3.6,label=paste("mRNA > Ribo and mRNA > Protein",numsigs[3],"genes"),adj=0,cex=1)
#   text(2,2.9,label=paste("Ribo > mRNA & Protein > mRNA",numsigs[4],"genes"),adj=0,cex=1)
# }


# ##############
# require(broman)
# plotDivergenceScatter <- function(xy, 
#                             ind_points1 = NULL, ind_points2 = NULL, 
#                             xlab, ylab,
#                             col_bkpoints = "black",
#                             cex_bkpoints = .4,
#                             lwd_bkpoints = .5,
#                             pch_bkpoints = 1,
#                             col_points1 = brocolors("crayons")["Red Orange"],
#                             col_points2 = brocolors("crayons")["Red Orange"],
#                             col_diag = "grey", 
#                             cex_points1 = .13, cex_points2 = .13,
#                             pch_points1 = 20, pch_points2 = 20,
#                             col_hline = brocolors("crayons")["Pine Green"], 
#                             col_vline = brocolors("crayons")["Pine Green"], 
#                             lwd_hline = .5,
#                             lwd_vline = .5,
#                             cex.axis = .3, cex.lab = .3,
#                             xlim = c(-6, 10), ylim = c(-6, 10) ) {
#     par( mar = c(2.2, 2, 1, .5), mgp = c(.7, .1, 0) )
#     plot( xy, pch = "", 
#           axes = F, ylim = ylim, xlim = xlim, xlab = "", ylab = "" )

#     points(xy, col = col_bkpoints, pch = pch_bkpoints, cex = cex_bkpoints, lwd = .5 )

#     if ( !is.null(ind_points2) ) {
#       points(xy[ind_points2, ], pch = pch_points2, col = col_points2,
#              cex = cex_points2)
#     }

#     points(xy[ind_points1, ], pch = pch_points1, col = col_points1,
#            cex = cex_points1)
    
#     if (xlab == "log2 (human/chimp) \n Ribosome occupancy") {
#       title( xlab = xlab, cex.lab = cex.lab, line = 1.3)
#     } else {
#       title( xlab = xlab, cex.lab = cex.lab)
#     }

#     title( ylab = ylab, cex.lab = cex.lab)

#     axis(1, tck = -.06, cex.axis = cex.axis, lwd = .6, at = c(-4, 0, 4, 8)
#          , label = c(-4, 0, 4, 8) )
#     axis(2, tck = -.03, cex.axis = cex.axis, lwd = .6, at = c(-4, 0, 4, 8)
#          , label = c(-4, 0, 4, 8) )
    
#     abline(h = 0, lty = 1, col = col_hline, lwd = lwd_hline)
#     abline(v = 0, lty = 1, col = col_vline, lwd = lwd_vline)
#     abline(a = 0, b = 1, lty = 1, col = col_diag, lwd = .5)
# }





# ##############
# require(broman)
# plotDivergenceScatter_goRiboRNA <- function(xy, 
#                                   ind_points1 = NULL, ind_points2 = NULL, 
#                                   ind_points3 = NULL,
#                                   xlab, ylab,
#                                   col_bkpoints = "black",
#                                   cex_bkpoints = .4,
#                                   lwd_bkpoints = .5,
#                                   pch_bkpoints = 1,
#                                   col_points1 = brocolors("crayons")["Red Orange"],
#                                   col_points2 = brocolors("crayons")["Red Orange"],
#                                   col_points3 = brocolors("crayons")["Black"],
#                                   col_diag = "grey", 
#                                   cex_points1 = .13, cex_points2 = .13,
#                                   cex_points3 = .13,
#                                   pch_points1 = 20, pch_points2 = 20,
#                                   pch_points3 = 20, 
#                                   col_hline = brocolors("crayons")["Pine Green"], 
#                                   col_vline = brocolors("crayons")["Pine Green"], 
#                                   lwd_hline = .5,
#                                   lwd_vline = .5,
#                                   cex.axis = .3, cex.lab = .3,
#                                   xlim = c(-6, 10), ylim = c(-6, 10) ) {
#   par( mar = c(2.2, 2, 1, .5), mgp = c(.7, .1, 0) )
#   plot( xy, pch = "", 
#         axes = F, ylim = ylim, xlim = xlim, xlab = "", ylab = "" )
  
#   points(xy, col = col_bkpoints, pch = pch_bkpoints, cex = cex_bkpoints, lwd = .5 )
  
#   if ( !is.null(ind_points2) ) {
#     points(xy[ind_points2, ], pch = pch_points2, col = col_points2,
#            cex = cex_points2)
#   }
  
#   points(xy[ind_points1, ], pch = pch_points1, col = col_points1,
#          cex = cex_points1)
  
#   points(xy[ind_points3, ], pch = pch_points3, col = col_points3,
#          cex = cex_points3)
  
#   if (xlab == "log2 (human/chimp) \n Ribosome occupancy") {
#     title( xlab = xlab, cex.lab = cex.lab, line = 1.3)
#   } else {
#     title( xlab = xlab, cex.lab = cex.lab)
#   }
  
#   title( ylab = ylab, cex.lab = cex.lab)
  
#   axis(1, tck = -.06, cex.axis = cex.axis, lwd = .6, at = c(-4, 0, 4, 8)
#        , label = c(-4, 0, 4, 8) )
#   axis(2, tck = -.03, cex.axis = cex.axis, lwd = .6, at = c(-4, 0, 4, 8)
#        , label = c(-4, 0, 4, 8) )
  
#   abline(h = 0, lty = 1, col = col_hline, lwd = lwd_hline)
#   abline(v = 0, lty = 1, col = col_vline, lwd = lwd_vline)
#   abline(a = 0, b = 1, lty = 1, col = col_diag, lwd = .5)
# }


# # plotDivergenceScatter.riborna.sigpro <- 
# #   function(xy.riborna,xy.rnapro,res.riborna,res.rnapro,pch,
# #            bg.riborna,col.riborna,
# #            bg.rnapro,col.rnapro,col.both,cex) {
# # #    dcols=densCols(xy.riborna,colramp=colorRampPalette(brewer.pal(9, "Greys")[-(1:2)]))
# #     plot(xy.riborna,pch="",
# #          xlab="mRNA log2 RPKM (Human/Chimp)", ylab="Ribo log2 RPKM (Human/Chimp)",
# #          axes=F,xlim=c(-7,7),ylim=c(-6,10))
# #     axis(1,at=seq(-6,6,2),labels=seq(-6,6,2))
# #     axis(2,at=seq(-4,8,2),labels=seq(-4,8,2))
# #     points(xy.riborna[res.rnapro$int.qval<.01 & (abs(xy.rnapro$rna)>abs(xy.rnapro$pro)),],
# #            pch=pch[1],bg=bg.rnapro[1],col=col.rnapro[1],cex=cex[1])
# #     points(xy.riborna[res.rnapro$int.qval<.01 & (abs(xy.rnapro$rna)<abs(xy.rnapro$pro)),],
# #            pch=pch[2],bg=bg.rnapro[2],col=col.rnapro[2],cex=cex[1])
# #     ii1 = res.riborna$int.qval<.01 & (abs(xy.riborna$rna)>abs(xy.riborna$ribo)) & res.rnapro$int.qval<.01 & (abs(xy.rnapro$rna) > abs(xy.rnapro$pro))
# #     ii2 = res.riborna$int.qval<.01 & (abs(xy.riborna$rna)<abs(xy.riborna$ribo)) & res.rnapro$int.qval<.01 & (abs(xy.rnapro$rna) < abs(xy.rnapro$pro))
# #     points(xy.riborna[ii1,],
# #            pch=pch[3],col=col.both[1],cex=cex[2])
# #     points(xy.riborna[ii2,],
# #            pch=pch[3],col=col.both[2],cex=cex[2])
# #     abline(h=0,v=0,lty=2,col="grey")
# # }
# # 
# # 
# # plotDivergenceScatter.rnapro.sigpro <- 
# #   function(xy.riborna,xy.rnapro,res.riborna,res.rnapro,pch,
# #            bg.riborna,col.riborna,
# #            bg.rnapro,col.rnapro,col.both,cex) {
# #     #    dcols=densCols(xy.rnapro,colramp=colorRampPalette(brewer.pal(9, "Greys")[-(1:2)]))
# #     plot(xy.rnapro,pch="",
# #          xlab="mRNA log2 RPKM (Human/Chimp)", ylab="Protein log2 FC (Human/Chimp)",
# #          axes=F,xlim=c(-7,7),ylim=c(-6,10))
# #     axis(1,at=seq(-6,6,2),labels=seq(-6,6,2))
# #     axis(2,at=seq(-4,8,2),labels=seq(-4,8,2))
# #     points(xy.rnapro[res.rnapro$int.qval<.01 & (abs(xy.rnapro$rna)>abs(xy.rnapro$pro)),],
# #            pch=pch[1],bg=bg.rnapro[1],col=col.rnapro[1],cex=cex[1])
# #     points(xy.rnapro[res.rnapro$int.qval<.01 & (abs(xy.rnapro$rna)<abs(xy.rnapro$pro)),],
# #            pch=pch[2],bg=bg.rnapro[2],col=col.rnapro[2],cex=cex[1])
# #     ii1 = res.riborna$int.qval<.01 & (abs(xy.riborna$rna)>abs(xy.riborna$ribo)) & res.rnapro$int.qval<.01 & (abs(xy.rnapro$rna) > abs(xy.rnapro$pro))
# #     ii2 = res.riborna$int.qval<.01 & (abs(xy.riborna$rna)<abs(xy.riborna$ribo)) & res.rnapro$int.qval<.01 & (abs(xy.rnapro$rna) < abs(xy.rnapro$pro))
# #     points(xy.rnapro[ii1,],
# #            pch=pch[3],col=col.both[1],cex=cex[2])
# #     points(xy.rnapro[ii2,],
# #            pch=pch[3],col=col.both[2],cex=cex[2])
# #     abline(h=0,v=0,lty=2,col="grey")
# #   }
# # 
# # 



# ###############
# # plotDivergenceScatter.riborna.sigribo <- 
# #   function(xy.riborna,xy.rnapro,res.riborna,res.rnapro,pch,
# #            bg.riborna,col.riborna,
# #            bg.rnapro,col.rnapro,col.both,cex) {
# # #    dcols=densCols(xy.riborna,colramp=colorRampPalette(brewer.pal(9, "Greys")[-(1:2)]))
# #     plot(xy.riborna,pch="",
# #          xlab="mRNA log2 RPKM (Human/Chimp)", ylab="Ribo log2 RPKM (Human/Chimp)",
# #          axes=F,xlim=c(-7,7),ylim=c(-6,10))
# #     axis(1,at=seq(-6,6,2),labels=seq(-6,6,2))
# #     axis(2,at=seq(-4,8,2),labels=seq(-4,8,2))
# #     points(xy.riborna[res.riborna$int.qval<.01 & (abs(xy.riborna$rna)>abs(xy.riborna$ribo)),],
# #            pch=pch[1],bg=bg.riborna[1],col=col.riborna[1],cex=cex[1])
# #     points(xy.riborna[res.riborna$int.qval<.01 & (abs(xy.riborna$rna)<abs(xy.riborna$ribo)),],
# #            pch=pch[2],bg=bg.riborna[2],col=col.riborna[2],cex=cex[1])
# #     ii1 = res.riborna$int.qval<.01 & (abs(xy.riborna$rna)>abs(xy.riborna$ribo)) & res.rnapro$int.qval<.01 & (abs(xy.rnapro$rna) > abs(xy.rnapro$pro))
# #     ii2 = res.riborna$int.qval<.01 & (abs(xy.riborna$rna)<abs(xy.riborna$ribo)) & res.rnapro$int.qval<.01 & (abs(xy.rnapro$rna) < abs(xy.rnapro$pro))
# #     points(xy.riborna[ii1,],
# #            pch=pch[3],col=col.both[1],cex=cex[2])
# #     points(xy.riborna[ii2,],
# #            pch=pch[3],col=col.both[2],cex=cex[2])
# #     abline(h=0,v=0,lty=2,col="grey")
# #   }
# # 
# # 
# # plotDivergenceScatter.rnapro.sigribo <- 
# #   function(xy.riborna,xy.rnapro,res.riborna,res.rnapro,pch,
# #            bg.riborna,col.riborna,
# #            bg.rnapro,col.rnapro,col.both,cex) {
# # #    dcols=densCols(xy.rnapro,colramp=colorRampPalette(brewer.pal(9, "Greys")[-(1:2)]))
# #     plot(xy.rnapro,pch="",
# #          xlab="mRNA log2 RPKM (Human/Chimp)", ylab="Protein log2 FC (Human/Chimp)",
# #          axes=F,xlim=c(-7,7),ylim=c(-6,10))
# #     axis(1,at=seq(-6,6,2),labels=seq(-6,6,2))
# #     axis(2,at=seq(-4,8,2),labels=seq(-4,8,2))
# #     points(xy.rnapro[res.riborna$int.qval<.01 & (abs(xy.riborna$rna)>abs(xy.riborna$ribo)),],
# #            pch=pch[1],bg=bg.riborna[1],col=col.riborna[1],cex=cex[1])
# #     points(xy.rnapro[res.riborna$int.qval<.01 & (abs(xy.riborna$rna)<abs(xy.riborna$ribo)),],
# #            pch=pch[2],bg=bg.riborna[2],col=col.riborna[2],cex=cex[1])
# #     ii1 = res.riborna$int.qval<.01 & (abs(xy.riborna$rna)>abs(xy.riborna$ribo)) & res.rnapro$int.qval<.01 & (abs(xy.rnapro$rna) > abs(xy.rnapro$pro))
# #     ii2 = res.riborna$int.qval<.01 & (abs(xy.riborna$rna)<abs(xy.riborna$ribo)) & res.rnapro$int.qval<.01 & (abs(xy.rnapro$rna) < abs(xy.rnapro$pro))
# #     points(xy.rnapro[ii1,],
# #            pch=pch[3],col=col.both[1],cex=cex[2])
# #     points(xy.rnapro[ii2,],
# #            pch=pch[3],col=col.both[2],cex=cex[2])
# #     abline(h=0,v=0,lty=2,col="grey")





# goHeatmapData <- function(golist){
# golist_test <- list(golist[[3]], golist[[4]])

# golist <- golist_test
# golist <- golist_bp

#   nlist <- length(golist)
#   goIDlist <- lapply( 1:nlist, function(i) {
#     golist[[i]][ , 1]
#   })
#   names(goIDlist) <- names(golist)

#   goTermlist <- lapply( 1:nlist, function(i) {
#     golist[[i]]$Term
#     })

#   goset <- unique( unlist(goIDlist) )
#   gosetTerm <- unique( unlist(goTermlist) )

#   heatmap_allgos <- data.frame(ID = goset)
#   heatmap_allgos$Term <- gosetTerm 
  
#   pvals <- merge(golist[[1]][, c(1,2,7)], golist[[2]][, c(1,2,7)], 
#             all = TRUE, by = "GOBPID")
#   pvals
#   # pvals <- lapply(1:nlist, function(i) {
#   #   foo <- rep(1-10^(-6), length(gosetTerm))       
#   #   foo[ which(goset %in% golist[[i]][ , 1]) ] <- golist[[i]]$Pvalue[ which( golist[[i]][ , 1] %in% goset )]
#   #   return(foo)
#   #   })
#   # pvals <- do.call(cbind, pvals)

# pvals_bp <- pvals
# pvals_test <- pvals

#   colnames(pvals) <- names(golist)
#   heatmap_allgos <- cbind(heatmap_allgos, pvals)
#   colnames(heatmap_allgos)[1] <- colnames(golist[[1]])[1]

#   ind <- which( is.na(heatmap_allgos), arr.ind = TRUE)
#   heatmap_allgos[ ind ] <- 1
#   rownames(heatmap_allgos) <- heatmap_allgos[ , 1]

#   iisig <- pvals > .01
#   iisig_gos <- which(rowSums(iisig)!=nlist)

#   heatmap_allgos <- heatmap_allgos[iisig_gos, ]

#   return(heatmap_allgos)  
#   }


# #' @example hmap <- plotHeatmap(heatmap_allgos, 
# #'                      labCol = c("Ribo > mRNA", "RNA > Ribo",
# #'                                "Ribo > Protein", "Protein > Ribo"))
# plotHeatmap <- function(heatmapData, mar =  c(0, 0, 0, 0),
#                         oma = c(0, 0, 0, 0), keysize = 1, cexRow = 1.3, cexCol = 1.3,
#                         labCol = NULL, key = TRUE, margin = c(5,5),
#                         dendrogram = "row", breaks = NULL) {  
# #heatmapData <- heatmap_allgos_bp
# # mar =  c(8, 1, 0, 20)
# # oma = c(2, 5, 0, 40) 
#     require(broman)
#     crayon <- brocolors("crayon")
#     heatcols <-  colorRampPalette(c("white", crayon["Banana Mania"],
#                                crayon["Laser Lemon"],
#                                crayon["Burnt Orange"], crayon["Orange Red"],
#                                crayon["Plum"]))
#     par( mfrow = c(1,1), mar = mar, oma = oma )
#     require(gplots)
#     mat <- heatmapData[ , -c(1:2)]
    
#     rownames(mat) <- heatmapData$Term
#     hmap <- heatmap.2(as.matrix(-log10(mat)), main = "",
#               labCol = labCol, cexCol = cexCol, 
#               cexRow = cexRow, 
#               breaks = breaks, 
#               sepcol = "grey20",
#               sepwidth=c(0.01, 0.01), lwd = .2,
#               colsep = 1:ncol(mat), Colv = FALSE, key = key,
#               rowsep = 1:nrow(mat), dendrogram = dendrogram, margin = margin,
#               srtCol = 0, adjCol = c(.5,NA), offsetCol = 0, offsetRow = .1,
#               keysize = keysize, lwid = c(1.5, 7), lhei = c(.5, 5),
#               col = heatcols, trace ="none", key.title = "-log p-value" )
#    return(hmap)
# }





# plotBoxes2 <- function(datavec1, datavec2, main, cols,
#                       labels, ylim = NULL, cex.main = .4, 
#                       main.line = .5) {
#   if (!is.null(ylim)) {
#   boxplot(datavec1, datavec2, 
#           outline=F, xlab="",lwd = .4,
#           axes = F, col= cols, border = "black", notch = TRUE, ylim = ylim,
#           varwidth = TRUE) 
#   } else {
#   boxplot(datavec1, datavec2, 
#           outline=F, xlab="",lwd = .4,
#           axes = F, col= cols, border = "black", notch = TRUE, varwidth = TRUE) 
#   }
#   title(main = main, cex.main = cex.main, line = main.line )
# #   axis(1, tck = -.03, cex.axis = .6, lwd = .6,
# #        at = 1:2, labels = labels, col= "white")
#   axis(2, tck = -.03, cex.axis = .4, lwd = .4)
# }


# #   }