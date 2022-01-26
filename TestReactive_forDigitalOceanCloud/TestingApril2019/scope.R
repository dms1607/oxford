# load the required packages
library(shiny)
require(shinydashboard)
library(ggplot2)

library(org.Mm.eg.db)
library(searchable)
library(devtools)
library(rentrez)

EntrezMap <- as.list(org.Mm.egALIAS2EG)
EntrezMap <- EntrezMap[!is.na(EntrezMap)]
EntrezMap_sl <- searchable(EntrezMap)
EntrezMap_sl <- ignore.case(EntrezMap_sl)

n=0
map <- c("test")

##########################

##load ("/Users/diego/Documents/DELFOS/pipeline.RData")

header <- dashboardHeader(title = "DELPHi 0.7")  

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

                box(title = "probando", status="warning", height="300", width="3", solidHeader=T,textOutput("textFunction"),style = "height:250px; overflow-y: scroll;overflow-x: scroll;")
                ),

 	fluidPage (

                box(title = "probando", status="warning", height="300", width="3", solidHeader=T,textOutput("textFunction2"),style = "height:250px; overflow-y: scroll;overflow-x: scroll;")
		)
	)


ui <- dashboardPage(title='DELPHi 0.7', header, sidebar, body, skin='black')

server <- function(input, output, session) {


	text_reactive <- eventReactive( input$submit, {
		
		symbol <- input$user_text
		map <<- EntrezMap_sl[symbol]
		map <<- map[[1]]

		genes <- entrez_summary(db="gene",id=c(map))
		geneFunc <- genes$summary
		geneFunc
		})		

	output$textFunction <- renderText({ text_reactive() })


 	text_reactive2 <- eventReactive( input$submit, {

		#symbol <- input$user_text
                #map <- EntrezMap_sl[symbol]
                #map <- map[[1]]


                genes <- entrez_summary(db="gene",id=c(map))
                geneFunc <- genes$summary
                geneFunc
                })

        output$textFunction2 <- renderText({ text_reactive2() })

    }

  
shinyApp(ui, server)

