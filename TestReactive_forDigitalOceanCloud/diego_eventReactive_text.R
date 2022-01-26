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

library(ggplot2)

load ("/Users/diego/Documents/DELFOS/pipeline.RData")


EntrezMap <- as.list(org.Mm.egALIAS2EG)
EntrezMap <- EntrezMap[!is.na(EntrezMap)]
EntrezMap_sl <- searchable(EntrezMap)
EntrezMap_sl <- ignore.case(EntrezMap_sl)

randomPhrase <- c("Alice in Wonderland")

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
	
 fluidPage (
	    	textOutput("text"),

		textOutput("text2"),

		box(title = "EVOLUTION IN PAPERS", status="warning", height=300, width="3", solidHeader=T, plotOutput("plot",height=225))
		)
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

 text_reactive2 <- eventReactive( input$submit, {

	ggplot(hits,aes(x=year,y=query_hits, group=1)) + geom_line(color='black') + xlab("") + ylab("") 
  })

  # text output
  output$text <- renderText({
    text_reactive()
  })

output$text2 <- renderText({
    text_reactive()
  })

output$plot <- renderPlot({ text_reactive2() })

}

  
shinyApp(ui, server)

