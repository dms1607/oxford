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

##load ("/Users/diego/Documents/DELFOS/pipeline.RData")

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

server <- function(input, output) {


	geneFunction_reactive <- eventReactive( input$submit, {
		
		symbol <<- input$user_text
		map <- EntrezMap_sl[symbol]
		map <<- map[[1]]

		genes <- entrez_summary(db="gene",id=c(map))
		geneFunc <- genes$summary
		geneFunc
		})		

	output$geneFunction <- renderText({ geneFunction_reactive() })

 	papers_reactive2 <- eventReactive( input$submit, {

		load ("/Users/diego/Documents/DELFOS/pipeline.RData")

		ggplot(hits,aes(x=year,y=query_hits, group=1)) + geom_line(color='black') + xlab("") + ylab("") 
  		})
	output$plot <- renderPlot({ papers_reactive2() })


### SIGNALING PATHWAYS

	pathways_reactive <- eventReactive ( input$submit, {

		keggret <- keggGet(kegg)
		pathways <- keggret[[1]]$PATHWAY

		kegg.df <- setNames(data.frame(matrix(ncol=2,nrow=length(pathways))),c("keggId","Pathway"))

		##diseases.df <- setNames(data.frame(matrix(ncol = 3, nrow = length(diseases))), c("Disease name","Disease ontology ID","Papers"))

		for (i in 1:length(pathways)) {
        
        		temp <- pathways[i]

        		kegg.df[i,1] <- names(temp)
        		kegg.df[i,2] <- temp[[1]]
        		}
		kegg.df 
		})

	output$pathways <- renderDataTable({ pathways_reactive() })


### PROTEIN INTERACTIONS

	ppi_reactive <- eventReactive ( input$submit, {

		ret <- select(mouseUp, keys =map, columns = c("ENTRY-NAME","UNIPROTKB","KEGG","MGI","INTERACTOR","ORGANISM","ORGANISM-ID"), keytype = "ENTREZ_GENE")

		entrezGene <<- ret$ENTREZ_GENE[1]
		entryName <<- ret$"ENTRY-NAME"[1]
		uniprotKB <<- ret$UNIPROTKB[1]
		kegg <<- ret$KEGG[1]
		mgi <<- ret$MGI[1]

		ppi <- ret$INTERACTOR[1]
		ppi <- gsub(" ", "", ppi, fixed = TRUE)

		ppi.v <- strsplit(ppi,";")
		ppi.df <- setNames(data.frame(matrix(ncol = 5, nrow = length(ppi.v[[1]]))), c("Query", "UniprotKB","Interaction","Papers citing both","URL")) 
		ppi.df$UniprotKB <- ppi.v[[1]]

		for (i in 1:nrow(ppi.df)) { 

        		ppi.df[i,1] <- entryName
        
        		UniprotId <- ppi.df[i,2]
        		if (grepl("Itself", UniprotId)) { UniprotId <- uniprotKB }

        		temp <- try (temp <- select(mouseUp, keys = UniprotId, columns = c("ENTRY-NAME"), keytype = "UNIPROTKB"))
       		 	stderr <- grep ("Error",temp)

        		if (length(stderr) > 0) {
                		temp <- try (temp <- select(humanUp, keys = UniprotId, columns = c("ENTRY-NAME"), keytype = "UNIPROTKB"))
                		}

        		ppi.df[i,3] <- temp$"ENTRY-NAME"[1]
        
        		secondProt <- gsub("_.*","",ppi.df$Interaction[i])
   
        		queryPPI <- paste(symbol,'AND',secondProt,sep=' ')

        		if (grepl("Itself", ppi.df$Uniprot[i])) { queryPPI <- paste(symbol,'AND','homodimer',sep=' ') }

        		ppi.df[i,4] <- epmc_hits(queryPPI)

 		       	###
        		### AQUI TENGO QUE PONER LO DE HOMODIMER EN LA LINEA PARA QUE ME LO CALCULE
        		###

		        URLline <- c("url") ##paste('https://europepmc.org/search?query=',symbol,'%20AND%20',secondProt,sep='')

        		ppi.df[i,5] <- URLline
        		}

		ppi.df <- ppi.df[order(ppi.df[,4],decreasing=TRUE),]
		ppi.df$UniprotKB <- NULL
		row.names(ppi.df) <- NULL

		ppi.df
		})
	
        output$ppi <- renderDataTable({ ppi_reactive() })


### DISEASES

	diseases_reactive <- eventReactive (input$submit, {

		queryDisease = getTemplateQuery(im = im, name = "MFeature_HDisease")

		queryDisease$where[[4]]$value <- symbol
		resDisease <- runQuery(im, queryDisease)

		diseases.df <- setNames(data.frame(matrix(ncol = 3, nrow = length(resDisease))), c("Disease name","Disease ontology ID","Papers"))

		for (i in 1:nrow(diseases.df)) {

        		diseases.df[i,1] <- resDisease[i,4]
        		diseases.df[i,2] <- resDisease[i,3]

        		diseaseString <- paste (symbol,"AND",diseases.df[i,1],sep=" ")
        		diseases.df[i,3] <- epmc_hits(diseaseString)
        		}

		diseases.df <- diseases.df[order(diseases.df$Papers, decreasing=TRUE),]
		diseases.df <- diseases.df[!duplicated(diseases.df), ]
		
		#row.names(diseases.df) <- NULL

		diseases.df
		
		})

	output$diseases <- renderDataTable({ diseases_reactive() })

### TISSUE EXPRESSION

	tissues_reactive <- eventReactive (input$submit, {

		queryGeneExp = getTemplateQuery(im = im,name = "Gene_Expression")

		## modify query for your gene

		queryGeneExp$where[[2]]$value <- symbol

		## Running the query

		resGeneExp <- runQuery(im, queryGeneExp)

		resGeneExp.df <- resGeneExp[,c("GXDExpression.age","GXDExpression.structure.name","GXDExpression.assayId")]
		colnames(resGeneExp.df) <- c("Dev stage", "Tissue", "Assay ID")
		tissues.df <- resGeneExp.df[!duplicated(resGeneExp.df), ]

		tissues.df
                })

        output$tissues <- renderDataTable({ tissues_reactive() })


### GENETIC PHENOTYPES

        phenotypes_reactive <- eventReactive (input$submit, {

		queryPhenotype = getTemplateQuery(im = im, name = "_Genotype_Phenotype")

		## modify query for your gene

		queryPhenotype$where[[4]]$value <- symbol

		## Running the query

		resPhenotype <- runQuery(im, queryPhenotype)

		resPhenotype.df <- resPhenotype[,c("OntologyAnnotation.subject.primaryIdentifier","OntologyAnnotation.subject.symbol","OntologyAnnotation.ontologyTerm.name","OntologyAnnotation.evidence.publications.pubMedId")]

		colnames(resPhenotype.df) <- c("Phenotype ID","Alleles","Phenotype","PubMed")

		resPhenotype.df <- unique(resPhenotype.df[c("Phenotype", "Alleles", "Phenotype ID","PubMed")])

		row.names(resPhenotype.df) <- NULL

		resPhenotype.df

                })

        output$phenotypes <- renderDataTable({ phenotypes_reactive() })


### PAPERS

	mostCited100_reactive <- eventReactive( input$submit, {

		PMCcited <- epmc_search(symbol,limit=100,synonym=FALSE,sort='cited')
		PMCcited_ann.df <- as.data.frame(PMCcited[,c('pubYear','pmid','authorString','title','journalTitle','citedByCount')])

		top100.df <- setNames(data.frame(matrix(ncol=6, nrow= nrow(PMCcited_ann.df))), c("Year","PMID","First author","Title","Journal","Citations"))


		for (i in 1:nrow(top100.df)) {


        		top100.df[i,]$Year <- PMCcited_ann.df[i,]$pubYear

        		top100.df[i,]$PMID <- PMCcited_ann.df[i,]$pmid

        		temp <- strsplit(PMCcited_ann.df[i,]$authorString,",")
        		top100.df[i,]$FirstAuthor <- temp[[1]][1] 
 
        		top100.df[i,]$Year <- PMCcited_ann.df[i,]$pubYear
 
        		temp <- strsplit(PMCcited_ann.df[i,]$authorString,",")
        		top100.df[i,]$FirstAuthor <- temp[[1]][1]

        		top100.df[i,]$Title <- PMCcited_ann.df[i,]$title
 
        		top100.df[i,]$Journal <- PMCcited_ann.df[i,]$journalTitle
 
        		top100.df[i,]$Citations <- PMCcited_ann.df[i,]$citedByCount
        		}
		top100.df

		})

	output$top100papers <- renderDataTable({ mostCited100_reactive() })
     

        mostCited20_reactive <- eventReactive( input$submit, {

		query2018 <- c("Stat3 AND FIRST_PDATE:2018")

		top20cited <- epmc_search(query = query2018,limit=20, synonym = FALSE, sort='cited')
		top20cited_ann.df <- as.data.frame(top20cited[,c('pmid','authorString','title','journalTitle','citedByCount')])

		top20.df <- setNames(data.frame(matrix(ncol = 6, nrow = nrow(top20cited_ann.df))), c("PMID","FirstAuthor","Title","Journal","Citations","MeSH"))

		for (i in 1:nrow(top20.df)) {

        		top20.df[i,]$PMID <- top20cited_ann.df[i,]$pmid

        		MeSH.list <- epmc_details(ext_id=top20.df$PMID[i])
        		MeSH_topic.v <- MeSH.list$mesh_topic$descriptorName
        
        		topicstring <- c("")

        		for (j in 1:length(MeSH_topic.v)) {

                		word <- gsub(" ","",MeSH_topic.v[j])            
                		word <- gsub(",","-",word)

                		topicstring <- paste(topicstring,word,sep=", ")
                		}

        	top20.df[i,]$MeSH <- substring(topicstring, 2) 

        	temp <- strsplit(top20cited_ann.df[i,]$authorString,",")
        	top20.df[i,]$FirstAuthor <- temp[[1]][1]

        	top20.df[i,]$Title <- top20cited_ann.df[i,]$title

        	top20.df[i,]$Journal <- top20cited_ann.df[i,]$journalTitle

        	top20.df[i,]$Citations <- top20cited_ann.df[i,]$citedByCount
        	}

 		top20cited_ann.df

                })

        output$top20papers <- renderDataTable({ mostCited20_reactive() })

	}

  
shinyApp(ui, server)

