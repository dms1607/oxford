#################
## DESCRIPTION ##
#################



###############
## LIBRARIES ##
###############

library(easypackages)
libraries("org.Mm.eg.db","searchable","devtools","rentrez","UniProt.ws","europepmc","KEGGREST","InterMineR","shiny","shinydashboard","DT","ggplot2")

######################
## GLOBAL VARIABLES ##
######################
mouseUp <- UniProt.ws(10090)
humanUp <- UniProt.ws(9606)

EntrezMap <- as.list(org.Mm.egALIAS2EG)
EntrezMap <- EntrezMap[!is.na(EntrezMap)]
EntrezMap_sl <- searchable(EntrezMap)
EntrezMap_sl <- ignore.case(EntrezMap_sl)

n=0

symbol <- c("test")
map <- c("test")

entrezGene <- c("test")
entryName <- c("test")
uniprotKB <- c("test")
kegg <- c("test")
mgi <- c("test")

im <- initInterMine(mine=listMines()["MouseMine"])
template=getTemplates(im)

##########################

header <- dashboardHeader(title = "DELPHi 0.7")  

sidebar <- dashboardSidebar(

  sidebarMenu(

	# input field
    	textInput("user_text", label = "Gene name", placeholder = "STAT3 for example"),

    	# submit button
    	actionButton("submit", label = "Submit & Wait")
    )
 )

body <- dashboardBody(
	
	fluidPage (


                box(title = "FUNCTION", status="warning", height="300", width="3", solidHeader=T,textOutput("geneFunction"),style = "height:250px; overflow-y: scroll;overflow-x: scroll;"),

                column(width = 3, height="300",

                	infoBox("SIGNALING PATHWAYS", n,color="yellow",width="NULL",icon = icon("angle-right"),fill = TRUE),
                	infoBox("TISSUES EXPRESSED IN", n,color="yellow",width="NULL",icon = icon("angle-right"),fill = TRUE),
                	infoBox("DEVELOPMENTAL EXPRESSION", n,color="yellow",width="NULL",icon = icon("angle-right"),fill = TRUE)                
                	),

                column(width = 3, height="300",

                	infoBox("GENETIC PHENOTYPES", n,color="yellow",width="NULL",icon = icon("angle-right"),fill = TRUE),
                	infoBox("MOST REPORTED INTERACTION", n,color="yellow",width="NULL",icon = icon("angle-right"),fill = TRUE),
                	infoBox("MOST REPORTED DISEASE", n,color="yellow",width="NULL",icon = icon("angle-right"),fill = TRUE)
                	),

                box(title = "EVOLUTION IN PAPERS", status="warning", height=300, width="3", solidHeader=T, plotOutput("plot",height=225))
                ),

 	fluidPage (

		box( title = "pathwayTitle", status="info", height = "400",width = "4",solidHeader = T, 
                        DT::dataTableOutput("pathways"),style = "height:350px; overflow-y: scroll;overflow-x: scroll;"),

                box( title = "ppiTitle", status="info", height = "400",width = "4",solidHeader = T,
                        DT::dataTableOutput("ppi"),style = "height:350px; overflow-y: scroll;overflow-x: scroll;"),

                box( title = "diseaseTitle", status="info", height = "400",width = "4",solidHeader = T,
                        DT::dataTableOutput("diseases"),style = "height:350px; overflow-y: scroll;overflow-x: scroll;")
		),

	fluidPage (
                box( title = "tissueTitle", status="info", height = "400",width = "6",solidHeader = T,
                        DT::dataTableOutput("tissues"),style = "height:350px; overflow-y: scroll;overflow-x: scroll;"),
               
		box( title = "phenotypesTitle", status="info", height = "400",width = "6",solidHeader = T,
                        DT::dataTableOutput("phenotypes"),style = "height:350px; overflow-y: scroll;overflow-x: scroll;")
                ),

        fluidPage (
                box( title = "100 MOST CITED EVER", status="info", height = "400",width = "6",solidHeader = T, 
                        DT::dataTableOutput("top100papers"),style = "height:350px;overflow-y: scroll;overflow-x: scroll;"),

                box( title = "20 MOST CITED OF 2018", status="info", height = "400",width = "6",solidHeader = T, 
                        DT::dataTableOutput("top20papers"),style = "height:350px;overflow-y: scroll;overflow-x: scroll;")
                )
	)


ui <- dashboardPage(title='DELPHi 0.7', header, sidebar, body, skin='black')

