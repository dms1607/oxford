#################
## DESCRIPTION ##
#################

## DELPHi v. 0.5
## Date: 3 June 2019

###############
## LIBRARIES ##
###############

library(easypackages)
libraries("org.Mm.eg.db","searchable","devtools","rentrez","UniProt.ws","europepmc","KEGGREST","InterMineR","shiny","shinydashboard","DT","ggplot2","stringr")

###################################################
## GLOBAL VARIABLES ###############################
###################################################

#######################
## PERMANENT OBJECTS ## 
#######################
mouseUp <- UniProt.ws(10090)
humanUp <- UniProt.ws(9606)

EntrezMap <- as.list(org.Mm.egALIAS2EG)
EntrezMap <- EntrezMap[!is.na(EntrezMap)]
EntrezMap_sl <- searchable(EntrezMap)
EntrezMap_sl <- ignore.case(EntrezMap_sl)

im <- initInterMine(mine=listMines()["MouseMine"])
template=getTemplates(im)

##########################
## FUNCTION DEFINITIONS ##
##########################

######################################
## 1. Mapping symbol to Entrez gene ##
######################################

EntrezMap <- function(x) {

        map <- EntrezMap_sl[x]
        map <<- map[[1]]        
	}

##################################################
## 2. UniProt.ws mapping to all other databases ##
##################################################

UniProtMapping <- function(x) {
	
	temp <- select(mouseUp, keys=x, columns = c("ENTRY-NAME","SCORE","UNIPROTKB","LENGTH","KEGG","MGI","INTERACTOR","ORGANISM","ORGANISM-ID"), keytype = "ENTREZ_GENE")

	for (i in 1:nrow(temp)){

        	if (temp$SCORE[i]=="5 out of 5") { temp$SCORE[i] <- c(5) }
		else if (temp$SCORE[i]=="4 out of 5") { temp$SCORE[i] <- c(4) }
                else if (temp$SCORE[i]=="3 out of 5") { temp$SCORE[i] <- c(3) }
                else if (temp$SCORE[i]=="2 out of 5") { temp$SCORE[i] <- c(2) }
                else if (temp$SCORE[i]=="1 out of 5") { temp$SCORE[i] <- c(1) }
        	}

	temp.df <- temp [order(temp$SCORE,decreasing=TRUE),] 

	return (temp.df)
        }

#######################################################
## 3. Get RefSeq function annotation using Entrez id ##
#######################################################

GeneFunction <- function(x) {

 
        genes <- entrez_summary(db="gene", id=c(x))
        geneFunc <- genes$summary

	if(geneFunc == "") { geneFunc <- genes$otherdesignations }
	if(geneFunc == "") { geneFunc <- genes$description }
        if(geneFunc == "") { geneFunc <- genes$name }

	return(geneFunc)
        }


#####################################################
## 4. Map KEGG identifier to annotate all pathways ##
#####################################################

KEGGmapping <- function(x) {

        keggret <- keggGet(x)
        pathways <- keggret[[1]]$PATHWAY

        kegg.df <- setNames(data.frame(matrix(ncol=2,nrow=length(pathways))),c("keggId","Pathway"))

        for (i in 1:length(pathways)) {

                temp <- pathways[i]

                kegg.df[i,1] <- names(temp)
                kegg.df[i,2] <- temp[[1]]
                }
	return(kegg.df)        
        }


############################################################## 
## 5. Protein interactions and papers citing both proteins ###
##############################################################

PPImapping <- function(x,y,z) {

        ppi <- x
	symbol <- y
	queryid <- z

        ppi <- gsub(" ", "", ppi, fixed = TRUE)

 	ppi.df <- matrix()

	if (is.na(ppi)) { ppi.df[1,1]  <- c("No interacting proteins found in the UniProt annotation.") } else {

		ppi.v <- strsplit(ppi,";")

		for (i in 1:length(ppi.v[[1]])) {
	                
			if (ppi.v[[1]][i]=="Itself") { ppi.v[[1]][i] <- c("Itself")}  else  {

				ppi.v[[1]][i] <- str_extract(ppi.v[[1]][i], "([A-Z0-9]{6})")
                		}
			}

		ppi.v[[1]] <- unique(ppi.v[[1]])

        	ppi.df <- setNames(data.frame(matrix(ncol = 5, nrow = length(ppi.v[[1]]))), c("Query", "UniprotKB","Interaction","Papers citing both","URL")) 
        	ppi.df$UniprotKB <- ppi.v[[1]]

       		for (i in 1:nrow(ppi.df)) { 

                		ppi.df[i,1] <- symbol
                		UniprotId <- ppi.df[i,2]
                		if (grepl("Itself", UniprotId)) { UniprotId <- queryid }

                		temp <- try (temp <- select(mouseUp, keys = UniprotId, columns = c("ENTRY-NAME"), keytype = "UNIPROTKB"))
                		stderr <- grep ("Error",temp)

                		if (length(stderr) > 0) {
                        		temp <- try (temp <- select(humanUp, keys = UniprotId, columns = c("ENTRY-NAME"), keytype = "UNIPROTKB"))
                        		}
				stderr <- grep ("Error",temp)

				if (length(stderr) > 0) {          
					ppi.df <- ppi.df[-c(i), ] 
					} 

				if (length(stderr) == 0) {

					ppi.df[i,3] <- temp$"ENTRY-NAME"[1]
        
			                secondProt <- gsub("_.*","",ppi.df$Interaction[i])
        	        		queryPPI <- paste(symbol,'AND',secondProt,sep=' ')

                			if (grepl("Itself", ppi.df$UniprotKB[i])) { queryPPI <- paste(symbol,'AND','homodimer',sep=' ') }

                			ppi.df[i,4] <- epmc_hits(queryPPI)

		                	###
                			### Search for HOMODIMER in case the interacting partner is the same as the query
                			###

					if (ppi.df$UniprotKB[i] == "Itself") { secondProt <- c("HOMODIMER") }

                			URLline <- paste('https://europepmc.org/search?query=',symbol,'%20AND%20',secondProt,sep='')

		                	ppi.df[i,5] <- URLline
                			}
			}	
		 	
		ppi.df <- ppi.df[order(ppi.df[,4],decreasing=TRUE),]    
        	ppi.df$UniprotKB <- NULL
        	row.names(ppi.df) <- NULL
		}

	linesDelete <- c()
	status=0
		
	for (i in 1:nrow(ppi.df)) { 
		if (is.na(ppi.df[i,1])) { 
			linesDelete <- c(linesDelete,i) 
			status=1 
			}
		       }
	if (status==1) {

		linesDelete <- rev(linesDelete)
		for (i in 1:length(linesDelete)) { ppi.df <- ppi.df[-c(linesDelete[i]), ] } 
		}
	
        return(ppi.df)
        }

#########################################################################
## InterMineR for Anatomical Expression, Genetic Phenotypes, Diseases  ##
#########################################################################

##############################
## 6. InterMineR - Diseases ##
##############################

InterMineR_dis <- function(x) {

        queryDisease = getTemplateQuery(
                im = im,
                name = "MFeature_HDisease"
                )

        queryDisease$where[[4]]$value <- x
        resDisease <- runQuery(im, queryDisease)

	diseases.df <- matrix()

	if (is.null(resDisease)) { diseases.df[1,1]  <- c("No human diseases have been found that are modeled by mutations in this gene.") } else {

		diseases.df <- setNames(data.frame(matrix(ncol = 3, nrow = nrow(resDisease))), c("Disease name","Disease ontology ID","Papers"))

        	for (i in 1:nrow(diseases.df)) {
                	diseases.df[i,1] <- resDisease[i,4]
                	diseases.df[i,2] <- resDisease[i,3]

                	diseaseString <- paste (x,"AND","mouse","AND",diseases.df[i,1],sep=" ")
                	diseases.df[i,3] <- epmc_hits(diseaseString)
                	}


        	diseases.df <- diseases.df[order(diseases.df$Papers, decreasing=TRUE),]
        	diseases.df <- diseases.df[!duplicated(diseases.df), ]
        	row.names(diseases.df) <- NULL
		}

        return(diseases.df)
        }

################################
## 7. InterMineR - Expression ##
################################

InterMineR_exp <- function(x) {

        queryGeneExp = getTemplateQuery(
                im = im,
                name = "Gene_Expression"
                )

        ## modify query for your gene

        queryGeneExp$where[[2]]$value <- x

        ## Running the query

        resGeneExp <- runQuery(im, queryGeneExp)

	tissues.df <- matrix()

	if (is.null(resGeneExp)) { tissues.df[1,1]  <- c("No differential tissue expression has been found, probably because this gene is constituutively expressed.") } else {

        	resGeneExp.df <- resGeneExp[,c("GXDExpression.age","GXDExpression.structure.name","GXDExpression.assayId")]
        	colnames(resGeneExp.df) <- c("Dev stage", "Tissue", "Assay ID")
        	tissues.df <- resGeneExp.df[!duplicated(resGeneExp.df), ]
		}
        return(tissues.df)
        }

########################################
## 8. InterMineR - Genetic phenotypes ##
########################################

InterMineR_pheno <- function(x) {

        queryPhenotype = getTemplateQuery(
                im = im,
                name = "_Genotype_Phenotype"
                )

        ## modify query for your gene

        queryPhenotype$where[[4]]$value <- x

        ## Running the query

        resPhenotype <- runQuery(im, queryPhenotype)

        resPhenotype.df <- resPhenotype[,c("OntologyAnnotation.subject.primaryIdentifier","OntologyAnnotation.subject.symbol","OntologyAnnotation.ontologyTerm.name","OntologyAnnotation.evidence.publications.pubMedId")]
        colnames(resPhenotype.df) <- c("Phenotype ID","Alleles","Phenotype","PubMed")
        resPhenotype.df <- unique(resPhenotype.df[c("Phenotype", "Alleles", "Phenotype ID","PubMed")])
        row.names(resPhenotype.df) <- NULL

        return(resPhenotype.df)
        }

##################
## Publications ##
##################

###########################
## 9. Top 100 most cited ##
###########################

Top100papers <- function(x) {

        PMCcited <- epmc_search(x,limit=100, synonym=FALSE, sort='cited')
        PMCcited_ann.df <- as.data.frame(PMCcited[,c('pubYear','pmid','authorString','title','journalTitle','citedByCount')])

        top100.df <- setNames(data.frame(matrix(ncol = 6, nrow = nrow(PMCcited_ann.df))), c("Year","PMID","FirstAuthor","Title","Journal","Citations"))

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
 
        return(top100.df)
        }

##########################
## 10. Top 20 last year ##
##########################

Top20papers <- function(x) {

        query <- paste(x,"AND FIRST_PDATE:2018", sep=" ")

        top20cited <- epmc_search(query = query,limit=20, synonym = FALSE, sort='cited')
        top20cited_ann.df <- as.data.frame(top20cited[,c('pmid','authorString','title','journalTitle','citedByCount')])

        top20.df <- setNames(data.frame(matrix(ncol = 5, nrow = nrow(top20cited_ann.df))), c("PMID","FirstAuthor","Title","Journal","Citations"))

        for (i in 1:nrow(top20.df)) {

                top20.df[i,]$PMID <- top20cited_ann.df[i,]$pmid

                temp <- strsplit(top20cited_ann.df[i,]$authorString,",")
                top20.df[i,]$FirstAuthor <- temp[[1]][1]
                top20.df[i,]$Title <- top20cited_ann.df[i,]$title
                top20.df[i,]$Journal <- top20cited_ann.df[i,]$journalTitle
                top20.df[i,]$Citations <- top20cited_ann.df[i,]$citedByCount
                }
        
        return (top20.df)
        }


############################
## SHINY DASHBOARD - BODY ##
############################

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

                        infoBox("SIGNALING PATHWAYS",textOutput("A1"),color="yellow",width="NULL",icon = icon("angle-right"),fill = TRUE),
                        infoBox("PROTEIN INTERACTIONS",textOutput("A2"),color="yellow",width="NULL",icon = icon("angle-right"),fill = TRUE),
                        infoBox("DISEASES", textOutput("A3"),color="yellow",width="NULL",icon = icon("angle-right"),fill = TRUE)
                        ),
                column(width = 3, height="300",

                        infoBox("GENETIC PHENOTYPES",textOutput("A4"),color="yellow",width="NULL",icon = icon("angle-right"),fill = TRUE),
                        infoBox("MOST REPORTED INTERACTION",textOutput("A5"),color="yellow",width="NULL",icon = icon("angle-right"),fill = TRUE),
                        infoBox("MOST STUDIED DISEASE",textOutput("A6"),color="yellow",width="NULL",icon = icon("angle-right"),fill = TRUE)
                        ),

                box(title = "EVOLUTION IN PAPERS", status="warning", height=300, width="3", solidHeader=T, plotOutput("plot",height=225))
                ),

 	fluidPage (

		box( title = "SIGNALLING PATHWAYS", status="info", height = "400",width = "4",solidHeader = T, 
                        DT::dataTableOutput("pathways"),style = "height:350px; overflow-y: scroll;overflow-x: scroll;"),

                box( title = "PROTEIN INTERACTIONS", status="info", height = "400",width = "4",solidHeader = T,
                        DT::dataTableOutput("ppi"),style = "height:350px; overflow-y: scroll;overflow-x: scroll;"),

                box( title = "HUMAN DISEASES MODELLED", status="info", height = "400",width = "4",solidHeader = T,
                        DT::dataTableOutput("diseases"),style = "height:350px; overflow-y: scroll;overflow-x: scroll;")
		),

	fluidPage (
                box( title = "TISSUE EXPRESSION", status="info", height = "400",width = "6",solidHeader = T,
                        DT::dataTableOutput("tissues"),style = "height:350px; overflow-y: scroll;overflow-x: scroll;"),
               
		box( title = "GENETIC PHENOTYPES", status="info", height = "400",width = "6",solidHeader = T,
                        DT::dataTableOutput("phenotypes"),style = "height:350px; overflow-y: scroll;overflow-x: scroll;")
                ),

        fluidPage (
                box( title = "100 MOST CITED", status="info", height = "700",width = "6",solidHeader = T, 
                        DT::dataTableOutput("top100papers"),style = "height:650px;overflow-y: scroll;overflow-x: scroll;"),

                box( title = "20 MOST CITED OF 2018", status="info", height = "700",width = "6",solidHeader = T, 
                        DT::dataTableOutput("top20papers"),style = "height:650px;overflow-y: scroll;overflow-x: scroll;")
                )
	)


ui <- dashboardPage(title='DELPHi 0.7', header, sidebar, body, skin='black')


#######################
## Create hyperlinks ##
#######################

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

createLink_diseasePapers <- function(diseaseName,papers,symbol) {

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

createLink_TOPpapers <- function(pmid,journal) {

        url <- paste('https://europepmc.org/abstract/MED/',pmid,sep='')

        line <- paste('<a href="',url,'"',' target=\"_blank\" class=\"btn btn-primary\">',journal,'</a>',sep='')
        print(line)
        }

##############################    
## SHINY DASHBOARD - SERVER ##
##############################

server <- function(input, output, session) {

	SYMBOL <- c()
	MAP <- c()	
	UNIPROTRET.df <- matrix()
	UNIPROTMAP.df <- matrix()
	SIGNALLINGPATHWAYSNO <- c()
	KEGG.df <- matrix()
	PROTEININTERACTIONSNO <- c()
	TOPINTERACTION <- c()
	PPI.df <- matrix()
	DISEASESNO <- c()
	TOPDISEASE <- c()
	DISEASES.df <- matrix()
	TISSUES.df <- matrix()
	GENETICPHENOTYPESNO <- c()
	PHENOTYPES.df <- matrix()
	TOP100.df <- matrix()
	TOP20.df <- matrix()
	HITS.df <- matrix()
   
	## 1. Gene function & UniProt mapping

 	geneFunction_reactive <- eventReactive( input$submit, {

                SYMBOL <<- input$user_text
                MAP <<- EntrezMap(SYMBOL)

		UNIPROTRET.df <<- UniProtMapping(MAP)

  		retMatrix <- matrix(nrow=5,ncol=1)
        	retMatrix[1,1] <- UNIPROTRET.df$ENTREZ_GENE[1]
        	retMatrix[2,1] <- UNIPROTRET.df$"ENTRY-NAME"[1]
        	retMatrix[3,1] <- UNIPROTRET.df$UNIPROTKB[1]
        	retMatrix[4,1] <- UNIPROTRET.df$KEGG[1]
        	retMatrix[5,1] <- UNIPROTRET.df$MGI[1]

		UNIPROTMAP.df <<- retMatrix

                geneFunc <- GeneFunction(MAP)
                geneFunc
                })

        output$geneFunction <- renderText({ geneFunction_reactive() })

	## 2. 

	pathways_reactive <- eventReactive( input$submit, {

		kegg.df <- KEGGmapping(UNIPROTMAP.df[4])
		kegg.df$"View pathway" <- createLink_pathways(kegg.df$keggId)
                kegg.df$keggId <- NULL
                kegg.df$URL <- NULL

        	SIGNALLINGPATHWAYSNO <<- nrow(kegg.df)
		KEGG.df <<- kegg.df
                KEGG.df
		})

        output$pathways <- renderDataTable({ 

		pathways_reactive()
		return(KEGG.df) }, escape = FALSE)

   	## 3.

        ppi_reactive <- eventReactive( input$submit, {

                ppi.df <- PPImapping(UNIPROTRET.df$INTERACTOR[1],SYMBOL,UNIPROTMAP.df[3])

                if (ppi.df[1,1]=="No interacting proteins found in the UniProt annotation.") {

			PROTEININTERACTIONSNO <<- c(0)
			TOPINTERACTION <<- c("None") 
			colnames(ppi.df)<- c("Interaction")		
			PPI.df <<- ppi.df } else {

				ppi.df$URL <- createLink_ppi (ppi.df$URL,ppi.df$"Papers citing both")
				ppi.df$"Papers citing both" <- NULL
				colnames(ppi.df)[3] <- "Papers"

        			PROTEININTERACTIONSNO <<- nrow(ppi.df)
        			TOPINTERACTION <<- ppi.df$Interaction[1]
		                PPI.df <<- ppi.df
				}                
		PPI.df
		})

        output$ppi <- renderDataTable({ 

		ppi_reactive() 
		return(PPI.df) }, escape=FALSE)

  	## 4.

        diseases_reactive <- eventReactive( input$submit, {

                diseases.df <- InterMineR_dis(SYMBOL)

		if (diseases.df[1,1]=="No human diseases have been found that are modeled by mutations in this gene.") {

			DISEASESNO <<- c(0)
                	TOPDISEASE <<- c("None") 
			DISEASES.df <<- diseases.df } else {

       				DISEASESNO <<- nrow(diseases.df)
        			TOPDISEASE <<- diseases.df[1,1]

				diseases.df$"Disease ontology ID" <- createLink_DOID(diseases.df[,2])
				diseases.df$"Papers" <- createLink_diseasePapers(diseases.df[,1],diseases.df[,3],SYMBOL)
          			DISEASES.df <<- diseases.df
				}
		DISEASES.df
		})

        output$diseases <- renderDataTable({ 

		diseases_reactive() 
		return (DISEASES.df) }, escape=FALSE)

	## 5.

        tissueExpression_reactive <- eventReactive( input$submit, {

                tissues.df <- InterMineR_exp(SYMBOL)

		if (tissues.df[1,1]=="No differential tissue expression has been found, probably because this gene is constituutively expressed.") {

                        TISSUES.df <<- tissues.df } else {

				tissues.df$"Assay ID" <- createLink_tissueAssay(tissues.df$"Assay ID")
                		TISSUES.df <<- tissues.df
				}

    		TISSUES.df
		})

        output$tissues <- renderDataTable({ 

		tissueExpression_reactive() 
		return (TISSUES.df) }, escape=FALSE)

 	## 6.

        phenotypes_reactive <- eventReactive( input$submit, {

                resPhenotype.df <- InterMineR_pheno(SYMBOL)

		resPhenotype.df$"Phenotype ID" <- createLink_phenotypeId(resPhenotype.df$"Phenotype ID")
                resPhenotype.df$PubMed <- createLink_Papers(resPhenotype.df$PubMed)
		GENETICPHENOTYPESNO <<- nrow(resPhenotype.df)
                PHENOTYPES.df <<- resPhenotype.df
                PHENOTYPES.df
		})

        output$phenotypes <- renderDataTable({ 

		phenotypes_reactive()
		return (PHENOTYPES.df) },escape=FALSE)


	## 7.

	top100papers_reactive <- eventReactive( input$submit, {

		top100.df <- Top100papers(SYMBOL)

 		top100.df$"PubMed" <- createLink_TOPpapers(top100.df$PMID,top100.df$Journal)
		colnames(top100.df)[colnames(top100.df) == 'PubMed'] <- 'Paper'
                top100.df$PMID <- NULL
                top100.df$Journal <- NULL
		TOP100.df <<- top100.df
		TOP100.df
		})
	
	output$top100papers <- renderDataTable({ 

		top100papers_reactive()
		return (TOP100.df) },escape=FALSE)

	## 8.

        top20papers_reactive <- eventReactive( input$submit, {

       		top20.df <- Top20papers(SYMBOL)

		top20.df$"PubMed" <- createLink_TOPpapers(top20.df$PMID,top20.df$Journal)
                colnames(top20.df)[colnames(top20.df) == 'PubMed'] <- 'Paper'
               
               	top20.df$PMID <- NULL
                top20.df$Journal <- NULL

               	TOP20.df <<- top20.df
                TOP20.df
		})

        output$top20papers <- renderDataTable({ 

		top20papers_reactive()
		return (TOP20.df) },escape=FALSE)

	## 9. Plot output
	
	plot_reactive <- eventReactive( input$submit, {

		hits.df <- epmc_hits_trend(SYMBOL, period = 2008:2018, synonym = FALSE)
		HITS.df <<- hits.df
		HITS.df
		})

        output$plot <- renderPlot({

                plot_reactive()
                ggplot(HITS.df,aes(x=year,y=query_hits, group=1)) + geom_line(color='black') + xlab("") + ylab("") 
	        });

  	## 10. Boxes

        A1_reactive <- eventReactive( input$submit, {
                n <- SIGNALLINGPATHWAYSNO
                n
                })

        output$A1 <- renderText({ A1_reactive() })

        A2_reactive <- eventReactive( input$submit, {
                n <- PROTEININTERACTIONSNO
                n
                })
        output$A2 <- renderText({ A2_reactive() })

        A3_reactive <- eventReactive( input$submit, {
                n <- DISEASESNO
                n
                })

        output$A3 <- renderText({ A3_reactive() })

        A4_reactive <- eventReactive( input$submit, {
                n <- GENETICPHENOTYPESNO
                n
                })

        output$A4 <- renderText({ A4_reactive() })

        A5_reactive <- eventReactive( input$submit, {
                n <- TOPINTERACTION
                n
                })
        output$A5 <- renderText({ A5_reactive() })

        A6_reactive <- eventReactive( input$submit, {
                n <- TOPDISEASE
                n
                })

        output$A6 <- renderText({ A6_reactive() })

        }

shinyApp(ui, server)

