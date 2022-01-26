# load the required packages
library(shiny)
require(shinydashboard)

library(org.Mm.eg.db)
library(searchable)

library(org.Mm.eg.db)
library(searchable)

library(devtools)
##install_github("ropensci/rentrez")
library(rentrez)

EntrezMap <- as.list(org.Mm.egALIAS2EG)
EntrezMap <- EntrezMap[!is.na(EntrezMap)]
EntrezMap_sl <- searchable(EntrezMap)
EntrezMap_sl <- ignore.case(EntrezMap_sl)

header <- dashboardHeader(title = "NextExperiment")  

sidebar <- dashboardSidebar(

  sidebarMenu(

	# input field
    	textInput("user_text", label = "Gene name", placeholder = "STAT3 for example"),

    	# submit button
    	actionButton("submit", label = "Submit")
    )
 )

body <- dashboardBody(
	
	    	textOutput("text")
		)

ui <- dashboardPage(title = 'This is my Page title', header, sidebar, body, skin='black')

server <- function(input, output) {

  # reactive expression
  text_reactive <- eventReactive( input$submit, {
    
	symbol <- input$user_text
	map <- EntrezMap_sl[symbol]
	map <- map[[1]]

	genes <- entrez_summary(db="gene", id=c(map))
	geneFunc <- genes$summary

	map
  })

  # text output
  output$text <- renderText({
    text_reactive()
  })
}

  
shinyApp(ui, server)

