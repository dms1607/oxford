  
## https://stackoverflow.com/questions/47505893/adding-a-vertical-and-horizontal-scoll-bar-to-the-dt-table-in-r-shiny

rm(list = ls())

load("/Users/diego/Documents/OXFORD2/Project/Rcode/pipeline.RData")

## app.R ##

####### Header titles ###############


functionTitle <- paste(c("FUNCTION ::",symbol), collapse=" ")

pathwayNo <- nrow(kegg.df)
pathwayTitle <- paste(c("SIGNALING PATHWAYS","(",pathwayNo,"pathways )"), collapse=" ")

ppiNo <- nrow(ppi.df)
ppiTitle <- paste(c("PROTEIN INTERACTIONS","(",ppiNo,"proteins )"), collapse=" ")
mostReportedPPI <- ppi.df[1,2]

diseaseNo <- nrow(diseases.df)
diseaseTitle <- paste(c("DISEASES","(",diseaseNo,"diseases )"), collapse=" ")
mostReportedDisease <- diseases.df[1,1]

devStageNo <- length(unique(tissues.df[,1]))
tissueNo <- length(unique(tissues.df[,2]))

assayNo <- length(unique(tissues.df[,3]))
tissueTitle <- paste(c("TISSUE EXPRESSION","(",devStageNo,"developmental stages,",tissueNo,"tissues,",assayNo,"assays )"), collapse=" ")

phenotypeNo <- nrow(resPhenotype.df)
phenotypesTitle <- paste(c("GENETIC PHENOTYPES","(",phenotypeNo,"phenotypes )"), collapse=" ")


######################################

library(shiny)
library(shinydashboard)
library(DT)
library(ggplot2)
library(PubMedWordcloud)

createLink_pathways <- function(keggId) {

        url <- paste('https://www.genome.jp/kegg-bin/show_pathway?',keggId,'+20848',sep='')

        line <- paste('<a href="',url,'"',' target=\"_blank\" class=\"btn btn-primary\">',keggId,'</a>',sep='')
        print(line)
        }

createLink_ppi <- function(url,papers) {

        line <- paste('<a href="',url,'"',' target=\"_blank\" class=\"btn btn-primary\">',papers,'</a>',sep='')
        print(line)
 	}

createLink_DOID <- function(doid) {

        url <- paste('http://disease-ontology.org/?id=',doid,sep='')
  
        line <- paste('<a href="',url,'"',' target=\"_blank\" class=\"btn btn-primary\">',doid,'</a>',sep='')
        print(line)
        }

createLink_diseasePapers <- function(diseaseName,papers) {

	url <- paste('https://europepmc.org/search?query=',symbol,'%20AND%20',diseaseName,sep='')

        line <- paste('<a href="',url,'"',' target=\"_blank\" class=\"btn btn-primary\">',papers,'</a>',sep='')
        print(line)
        }

createLink_tissueAssay <- function(assay) {

        url <- paste('http://www.informatics.jax.org/assay/',assay,sep='')

        line <- paste('<a href="',url,'"',' target=\"_blank\" class=\"btn btn-primary\">',assay,'</a>',sep='')
        print(line)
        }
	

createLink_phenotypeId <- function(phenot) {

	url <- paste('http://www.informatics.jax.org/allele/genoview/',phenot,sep='')

        line <- paste('<a href="',url,'"',' target=\"_blank\" class=\"btn btn-primary\">',phenot,'</a>',sep='')
        print(line)
        }

createLink_Papers <- function(pmid) {

        url <- paste('https://europepmc.org/abstract/MED/',pmid,sep='')

        line <- paste('<a href="',url,'"',' target=\"_blank\" class=\"btn btn-primary\">',pmid,'</a>',sep='')
        print(line)
        }

body <- dashboardBody(

	fluidPage (

                box(title = functionTitle, status="warning", height="300", width="3", solidHeader=T,geneFunc,
		style = "height:250px; overflow-y: scroll;overflow-x: scroll;"),

                column(width = 3, height="300",

                infoBox("SIGNALING PATHWAYS", pathwayNo,color="yellow",width="NULL",icon = icon("angle-right"),fill = TRUE),
                infoBox("TISSUES EXPRESSED IN", tissueNo,color="yellow",width="NULL",icon = icon("angle-right"),fill = TRUE),
                infoBox("DEVELOPMENTAL STAGES - EXPRESSION", devStageNo,color="yellow",width="NULL",icon = icon("angle-right"),fill = TRUE)                
		),

		column(width = 3, height="300",

                infoBox("GENETIC PHENOTYPES", phenotypeNo,color="yellow",width="NULL",icon = icon("angle-right"),fill = TRUE),
                infoBox("MOST REPORTED INTERACTION", mostReportedPPI,color="yellow",width="NULL",icon = icon("angle-right"),fill = TRUE),
		infoBox("MOST REPORTED DISEASE", mostReportedDisease,color="yellow",width="NULL",icon = icon("angle-right"),fill = TRUE)
		),

                box(title = "EVOLUTION IN PAPERS", status="warning", height=300, width="3", solidHeader=T, plotOutput("plot",height=225))
                ),

    		
	fluidPage (
                box( title = pathwayTitle, status="primary", height = "400",width = "4",solidHeader = T, 
                        DT::dataTableOutput("pathways"),style = "height:350px; overflow-y: scroll;overflow-x: scroll;"),

                box( title = ppiTitle, status="primary", height = "400",width = "4",solidHeader = T,
                        DT::dataTableOutput("ppi"),style = "height:350px; overflow-y: scroll;overflow-x: scroll;"),

	 	box( title = diseaseTitle, status="primary", height = "400",width = "4",solidHeader = T,
                        DT::dataTableOutput("diseases"),style = "height:350px; overflow-y: scroll;overflow-x: scroll;")
                ),

	fluidPage (
                box( title = tissueTitle, status="primary", height = "400",width = "6",solidHeader = T,
                        DT::dataTableOutput("tissues"),style = "height:350px; overflow-y: scroll;overflow-x: scroll;"),
		box( title = phenotypesTitle, status="primary", height = "400",width = "6",solidHeader = T,
                        DT::dataTableOutput("phenotypes"),style = "height:350px; overflow-y: scroll;overflow-x: scroll;")
                ),

	fluidPage (
		box( title = "100 MOST CITED PAPERS", status="primary", height = "400",width = "6",solidHeader = T, 
                        DT::dataTableOutput("top100papers"),style = "height:350px;overflow-y: scroll;overflow-x: scroll;"),

		box( title = "20 MOST CITED PAPERS OF 2021", status="primary", height = "400",width = "6",solidHeader = T, 
                        DT::dataTableOutput("top20papers"),style = "height:350px;overflow-y: scroll;overflow-x: scroll;")
		)
	      )


sidebar <- dashboardSidebar(

  sidebarMenu(
    menuItem("Dashboard", tabName = "dashboard", icon = icon("dashboard")),

    textInput("searchText", label = "Enter gene name", value = "STAT3")
    )
 )


ui <- dashboardPage(
  skin="blue",  
  dashboardHeader(title = "DELPHi 0.7"),
  sidebar,
  body
  )

shinyApp (

  ui=ui,

  server=function(input, output) { 

        output$plot <- renderPlot({ ggplot(hits,aes(x=year,y=query_hits, group=1)) + geom_line(color='black') + xlab("") + ylab("") });

       	output$pathways <- renderDataTable({ 

			kegg.df$"View pathway" <- createLink_pathways(kegg.df$keggId)
                        kegg.df$keggId <- NULL
			kegg.df$URL <- NULL
                        return(kegg.df) }, escape = FALSE);
  
	output$ppi <- renderDataTable({ 
			ppi.df$"Papers citing both" <- createLink_ppi(ppi.df$URL,ppi.df[,3])
			ppi.df$URL <- NULL
			return(ppi.df) }, escape = FALSE);

	output$diseases <- renderDataTable({
                        
			diseases.df$"Disease info" <- createLink_DOID(diseases.df[,2])
			diseases.df$"Papers" <- createLink_diseasePapers(diseases.df[,1],diseases.df[,3])
			diseases.df[,2] <- NULL

                        return(diseases.df) }, escape = FALSE);

        output$tissues <- renderDataTable({ 

			tissues.df$"Assay ID" <- createLink_tissueAssay(tissues.df$"Assay ID")
                        return(tissues.df) }, escape = FALSE);

	output$phenotypes <- renderDataTable({

                        resPhenotype.df$"Phenotype ID" <- createLink_phenotypeId(resPhenotype.df$'Phenotype ID')
			resPhenotype.df$"PubMed" <- createLink_Papers(resPhenotype.df$PubMed)
                        return(resPhenotype.df) }, escape = FALSE);

        output$top100papers <- renderDataTable({

                        top100.df$"PubMed" <- createLink_Papers(top100.df$PMID)
                        top100.df$PMID <- NULL

                        return(top100.df) }, escape = FALSE);

        output$top20papers <- renderDataTable({

                        ##datatable(top20.df, options = list(paging = FALSE)) });

                        top20.df$"PubMed" <- createLink_Papers(top20.df$PMID)
                        top20.df$PMID <- NULL
                        top20.df$MeSH <- NULL
                        return(top20.df) }, escape = FALSE);
        }

     )

