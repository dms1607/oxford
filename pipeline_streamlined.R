## This is the R pipeline for the initial analysis for DELPHI

if (!require("BiocManager", quietly = TRUE))
    install.packages("BiocManager")

#BiocManager::install("org.Mm.eg.db")
#BiocManager::install("annotate")


############################
## 0-The gene of interest ##
############################

symbol <- c("Stat3")

#####################################
## 1-Mapping Symbol to Entrez gene ##
#####################################

library(org.Mm.eg.db)
map <- mapIds(org.Mm.eg.db,symbol,'ENTREZID','SYMBOL')
map <- map[[1]]


##########################################################
## 2-get function via RefSeq annotation using Entrez id ##
##########################################################

library(devtools)
#install_github("ropensci/rentrez")
library(rentrez)

genes <- entrez_summary(db="gene", id=c(map))
geneFunc <- genes$summary

#####################################################
## 3-get UniProt.ws mapping to all other databases ##
#####################################################

#BiocManager::install("UniProt.ws")
#install.packages("europepmc",dependencies=TRUE)

library(UniProt.ws)
library(europepmc)
mouseUp <- UniProt.ws(10090)
humanUp <- UniProt.ws(9606)

ret <- select(mouseUp, keys =map, columns = c("ENTRY-NAME","UNIPROTKB","KEGG","MGI","INTERACTOR","ORGANISM","ORGANISM-ID"), keytype = "ENTREZ_GENE")

entrezGene <- ret$ENTREZ_GENE[1]
entryName <- ret$"ENTRY-NAME"[1]
uniprotKB <- ret$UNIPROTKB[1]
kegg <- ret$KEGG[1]
mgi <- ret$MGI[1]

########## PPIs AND PAPERS CITING BOTH INTERACTING PROTEINS ##########

ppi <- ret$INTERACTOR[1]
ppi <- gsub(" ", "", ppi, fixed = TRUE)

ppi.v <- strsplit(ppi,";")
ppi.df <- setNames(data.frame(matrix(ncol = 5, nrow = length(ppi.v[[1]]))), c("Query", "UniprotKB","Interaction","Papers citing both","URL")) 
ppi.df$UniprotKB <- ppi.v[[1]]

createLink <- function(entryname, url) {
	line <- paste(url,'target=\"_blank\"','class=\"btn btn-primary\">',entryname,'</a>',sep=' ')
  	sprintf(line)
 	}

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

	URLline <- paste('https://europepmc.org/search?query=',symbol,'%20AND%20',secondProt,sep='')

	ppi.df[i,5] <- URLline
	}

ppi.df <- ppi.df[order(ppi.df[,4],decreasing=TRUE),]
ppi.df$UniprotKB <- NULL
row.names(ppi.df) <- NULL

######################################################################
## 4-map KEGG identifier using KEGGREST(R) to retrieve all pathways ##
######################################################################

#BiocManager::install("KEGGREST")
library(KEGGREST)

keggret <- keggGet(kegg)
pathways <- keggret[[1]]$PATHWAY

kegg.df <- setNames(data.frame(matrix(ncol=2,nrow=length(pathways))),c("keggId","Pathway"))

for (i in 1:length(pathways)) {
	
	temp <- pathways[i]

        kegg.df[i,1] <- names(temp)
        kegg.df[i,2] <- temp[[1]]
	}


###########################################################################
## 5-InterMineR for Anatomical Expression, Genetic Phenotypes, Diseases  ##
###########################################################################

#BiocManager::install("InterMineR")
library(InterMineR)

im <- initInterMine(mine=listMines()["MouseMine"])
template = getTemplates(im)

##################
### Expression ###
##################

#template[grep("expression", template$name, ignore.case=TRUE),]

queryGeneExp = getTemplateQuery(
  im = im,
  name = "Gene_Expression"
 )

## modify query for your gene

queryGeneExp$where[[2]]$value <- symbol

## Running the query

resGeneExp <- runQuery(im, queryGeneExp)

resGeneExp.df <- resGeneExp[,c("GXDExpression.age","GXDExpression.structure.name","GXDExpression.assayId")]
colnames(resGeneExp.df) <- c("Dev stage", "Tissue", "Assay ID")
tissues.df <- resGeneExp.df[!duplicated(resGeneExp.df), ]

##########################
### Genetic phenotypes ###
##########################

#template[grep("phenotype", template$name, ignore.case=TRUE),]

queryPhenotype = getTemplateQuery(
  im = im,
  name = "_Genotype_Phenotype"
 )

## modify query for your gene

queryPhenotype$where[[4]]$value <- symbol

## Running the query

resPhenotype <- runQuery(im, queryPhenotype)

resPhenotype.df <- resPhenotype[,c("OntologyAnnotation.subject.primaryIdentifier","OntologyAnnotation.subject.symbol","OntologyAnnotation.ontologyTerm.name","OntologyAnnotation.evidence.publications.pubMedId")]
colnames(resPhenotype.df) <- c("Phenotype ID","Alleles","Phenotype","PubMed")
resPhenotype.df <- unique(resPhenotype.df[c("Phenotype", "Alleles", "Phenotype ID","PubMed")])
row.names(resPhenotype.df) <- NULL

################
### Diseases ###
################

#template[grep("disease", template$name, ignore.case=TRUE),]

queryDisease = getTemplateQuery(
  im = im,
  name = "MFeature_HDisease"
 )

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
row.names(diseases.df) <- NULL

######################
## 7-Hottest papers ##
######################

library(europepmc)

##### Get top 100 papers in the field ###############

PMCcited <- epmc_search(symbol,limit=100, synonym=FALSE, sort='cited')
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
 
## Top 20 this year

query2021 <- paste(symbol,"AND FIRST_PDATE:2021")

top20cited <- epmc_search(query = query2021,limit=20, synonym = FALSE, sort='cited')
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

##### Get publications 2000-2022 to build a line plot *****

library(ggplot2)
hits <- epmc_hits_trend(symbol, period = 2000:2022, synonym = FALSE)
hits$year <- factor (hits$year)

png(filename="/Users/diego/Documents/OXFORD2/Project/Rcode/hits.png")
ggplot(hits,aes(x=year,y=query_hits, group=1)) + geom_line(color='blue') + xlab("") + ylab("")
dev.off()

################# Top paper of 2018 ####################

topPaper2021 <- top20.df[1,]$Title

