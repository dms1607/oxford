## server.R

library(shiny)

library(searchable)
library(devtools)
library(rentrez)
library(europepmc)
library(shiny)
library(shinydashboard)
library(DT)
library(ggplot2)
library(stringr)


shinyServer(function(input, output) {

  output$distPlot <- renderPlot({

    # generate bins based on input$bins from ui.R
    x    <- faithful[, 2]
    bins <- seq(min(x), max(x), length.out = input$bins + 1)

    # draw the histogram with the specified number of bins
    hist(x, breaks = bins, col = 'darkgray', border = 'white')

  })

})

