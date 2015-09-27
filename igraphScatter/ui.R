# Set up the layout of the shiny page
library(shiny)
shinyUI(fluidPage(

  # Application title
  titlePanel("Coefficient of variation and count data"),

  # Page layout
  fluidRow(
    column(width = 5,
           # CV-versus-mean plot
           plotOutput("plotCV",
                      click = clickOpts(id = "plotCV_click", clip = TRUE),
                      hover = hoverOpts(id = "plotCV_hover", nullOutside = FALSE),
                      width = 500, height = 500)
    ),
    column(width = 4, offset = 1,
          # the counts plot for a single gene
          plotOutput("plotGene", width=500, height=500)
    )
  )

))


