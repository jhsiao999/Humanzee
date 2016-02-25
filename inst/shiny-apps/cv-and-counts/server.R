# Define input and output of each plot
library(shiny)

shinyServer(function(input, output) {
  output$contents <- renderTable({
    
    # input$file1 will be NULL initially. After the user selects
    # and uploads a file, it will be a data frame with 'name',
    # 'size', 'type', and 'datapath' columns. The 'datapath'
    # column will contain the local filenames where the data can
    # be found.

    inFile <- input$file1

    if (is.null(inFile))
      return(NULL)
    
    read.csv(inFile$datapath, header=input$header, sep=input$sep, 
         quote=input$quote)
  })
})

# shinyServer(function(input, output) {
#   selectPoint <- reactiveValues(idx = NULL)

#   observe({
#     xy <- c(input$plotCV_click$x, input$plotCV_click$y)
#     if (!is.null(xy)) {
#       pt <- nearPoints(df = df_cvplot, coordinfo = input$plotCV_click,
#                        xvar = "mean", yvar = "cv",
#                        maxpoints = 1)
#       selectPoint$idx <- match(rownames(pt), rownames(df_cvplot))
#     }
#   })

#   output$plotCV <- renderPlot({
#     par( mar = c(5, 5, 3, 2), cex.lab=2, cex.main = 2, cex.axis = 1.5)
#     plot(df_cvplot, xlab = "log2 mean gene count",
#          ylab = "Coefficient of variation", pch = 16, cex = .8,
#          axes = F)
#     title(main = per_person)
#     axis(1); axis(2)

#     idx <- selectPoint$idx
#     if(!is.null(idx)) {
#       points( df_cvplot$mean[idx], df_cvplot$cv[idx], col="dodgerblue", cex=3, lwd=3 )
#     }
#   })

#   # counts plot
#   output$plotGene <- renderPlot({
#     idx <- selectPoint$idx
#     par(mar = c(5,5,3,4), cex.lab = 2, cex.main = 1.8, cex.axis=1.5)
#     # Density plot for the selected gene

#     if(is.null(idx)) {
#     plot(0, pch = "", axes = F, ann = F, xlim = c(1,5), ylim = c(1,5))
#     text(3, 3, "Select a point!", cex = 3)
#     }

#     if (!is.null(idx)) {
#     dens <- density(unlist(df_geneplot[ idx, ]) )
#     gene <- rownames(df_geneplot)[idx]
#     plot(dens, ylab = "Density", xlab = "Count",
#          main = paste(gene, "\n count density across cells"), axes = F)
#     axis(1); axis(2)
#     }
#   })

# })
