library(ggplot2) # v2.2.1
library(clusterProfiler)# v3.2.11
library(GO.db) # v2.9
library(org.Hs.eg.db) # v3.4.0
library(topGO) # v2.26.0
library(shiny) # v1.0.0
library(pathview) # v1.14
library(biomaRt)
# Bioconductor 3.4

##################
### Fonctions ####
##################

## Proteins Domains ##

# Nettoyer les ids vides renvoys par BIOCONDUCTOR, prend un entre un objet getBM()   
cleanBlankLines<- function(datatoclean){ 
    count=1   
    datacleaned=""   
    for(i in which(as.vector(datatoclean)!="")) {
        datacleaned[count]=as.matrix(datatoclean)[i]
        count=count+1
    }
    return(datacleaned)
}

#retourne metadonnes d'une database (arg bioCdatabase) mapps sur ids panther family de l'echantillon
#Prend en entre un rsultat de cleanBlankLines()
add_metadata <- function(bioCdatabase) {
    metadata=""
    for(i in 1:length(Panther_Family_ID)) {
        #if domain id exists : do id retrieving
        if(Panther_Family_ID[i]!="") {
            metadata[i]=getBM(attributes=c(bioCdatabase),
                              filters ='family',
                              values = Panther_Family_ID[i],
                              mart = ensembl) 
        }
        #if domain id is blank (bug happened somewhere) : add blank line to metadata[i]
        #ce else n'est normalement plus utile depuis l'implementation de la fonction cleanBlankLines qui nettoie les ids vides
        else {
            metadata[i]="-"
        }
    }
    return(metadata)
}

## Pour Goenrichment ##

# DropGO
drop <- function(enrichment, drop1, drop2) {
    dropped <- NULL
    for(i in 1:drop1) {
        dropped <- dropGO(enrichment,level=i)
    }
    for(i in drop2:14) {
        dropped <- dropGO(dropped,level=i)
    }
    return(dropped)
}

##############
###  Page ###
##############

ui <- fluidPage(
    titlePanel("Gene Set Enrichment Analysis"),
    wellPanel(fileInput("data", "Choose a txt file", multiple=FALSE, accept = NULL)),
  
    mainPanel(
        
        ##### Onglets d'Input #####
      
        tabsetPanel(
            
            ### Onglet Data
            tabPanel("Data", id="data", dataTableOutput("contents")),
        
            ### Onglet GO
            tabPanel(
                "GO analysis", id="fgo",
                radioButtons("dogroup", label="DO Go classification ?", choices=c("Yes"="y", "No"="n"), selected="n"),
                numericInput(inputId="golvl", label="Level for go classification", value=7, min=1, max=14, step=1),
                sliderInput(inputId="gop", label="P-value cut-off", value=0.01, min=0.001, max=0.05 ),
                sliderInput(inputId="goq", label="q-value cut-off", value=0.01, min=0.001, max=0.05 ),
                sliderInput(inputId="gocat", label="categories to show", value=12, min=8, max=18),
                sliderInput(inputId="drop", label="Keep GO levels between", min = 1, max = 14, value = c(1,14)),
                selectInput(inputId="goont", label="Ontology", c("CC", "BP", "MF")),
                selectInput(inputId="adjm", label="p-value adjustment method", c("holm", "hochberg", "hommel", "bonferroni", "BH", "BY", "fdr", "none")),
                actionButton("startgo", "Start")
            ),
                   
            ### Onglet Pathways
            tabPanel("Pathways analysis", id="fkegg",
                sliderInput(inputId="keggp", label="P-value cut-off", value=0.01, min=0.001, max=0.05 ),
                sliderInput(inputId="keggq", label="q-value cut-off", value=0.01, min=0.001, max=0.05 ),
                selectInput(inputId="adjmk", label="p-value adjustment method", c("holm", "hochberg", "hommel", "bonferroni", "BH", "BY", "fdr", "none")),
                actionButton("startkegg", "Start")
            ),
        
            ### Onglet Protein Domains
            tabPanel("ProteinDomains", id="fpd",
                sliderInput(inputId="pdp", label="P-value cut-off", value=0.01, min=0.001, max=0.05 ),
                sliderInput(inputId="pdq", label="q-value cut-off", value=0.01, min=0.001, max=0.05 ),
                actionButton("startpd", "Start")
            )
        ),
        
        ##### Onglets d'Outputs  #####
        
        tabsetPanel(
            tabPanel("Plots", id="plots", plotOutput("vulcano")),
            tabPanel("GO results", id="rgo", 
                dataTableOutput("groupGo"), 
                plotOutput("enrichGo"), 
                plotOutput("enrichGodot"),
                plotOutput("gograph")
            ),
            tabPanel("Pathways results", id="rkegg", 
                uiOutput("pathway1"),
                uiOutput("pathway2"),
                uiOutput("pathway3")
            ),
            tabPanel("Protein Domains results", id="rpd", "Protein Domains")
        )
    ) # fin main panel
) # fin fluid page

###############
###  SERVER ###
###############

server <- function(input, output, session) {
    ### Data ###################
    data <- reactive({
        if (is.null(input$data)) {return(NULL)}
        read.table(input$data$datapath, header = TRUE, dec=",")
    })
  
    output$contents <- renderDataTable(data(), options=list(pageLength=10))
  
    ### Vulcano ################
    output$vulcano <- renderPlot({
        if (is.null(data())) {return(NULL)}
        ggplot(data(), aes(x=log2FC, y=-log10(adj_pvalue))) +
            geom_point(size = 1.5, alpha = 0.7, na.rm = T) +
            scale_x_continuous(name="log2FC") +
            scale_y_continuous(name="-log10(p_value)",limits=c(0,200))
    })
  
    genes <- reactive({as.character(data()$Ensembl_Gene_ID)})
    entrezIDs <- reactive({ bitr(genes(), fromType="ENSEMBL" , toType="ENTREZID", OrgDb="org.Hs.eg.db")$ENTREZID})
  
    ### GO ########################
  
    observe({
        if(input$startgo != 0) {
            # groupgo
            if(input$dogroup=="y") {
                gogroup <- reactive({ groupGO(entrezIDs(), OrgDb='org.Hs.eg.db', ont=input$goont, level=input$golvl , readable = TRUE) })
                output$groupGo <- renderDataTable(as.data.frame(gogroup()), options=list(pageLength=10))
            }
            # enrichGo
            goenrich <- reactive({ 
                enrichGO(entrezIDs(), OrgDb="org.Hs.eg.db", ont=input$goont, pAdjustMethod=input$adjm, pvalueCutoff=input$gop, qvalueCutoff=input$goq, readable= TRUE) 
            })
            enrichDrop <-reactive({ 
                drop(goenrich(), input$drop[1], input$drop[2]) 
            })
            output$enrichGo <- renderPlot({barplot(enrichDrop(), drop = TRUE, showCategory=input$gocat, title="Go Enrichment")}) #goenrich
            output$enrichGodot <- renderPlot({ dotplot(enrichDrop(), showCategory=input$gocat, title="GO Enrichment") }) #goenrich
      
            # svg file plotgograph
            output$gograph <- renderImage({
                go = enrichDrop() #goenrich
                width  <- session$clientData$output_gograph_width
                height <- session$clientData$output_gograph_height
                outfile <- tempfile(fileext=".svg")
                svg(outfile)
                plotGOgraph(go)
                dev.off()
        
                list(src = outfile, contentType = 'image/svg+xml', width = width, height = height, alt = "plotGoGraph")
            }, deleteFile = TRUE)
        }
    })
  
    ### Pathways ###################
  
    observe({
        if (input$startkegg != 0) {
      
            kegg <- reactive({ enrichKEGG(entrezIDs(), pAdjustMethod=input$adjmk, pvalueCutoff=input$keggp, qvalueCutoff=input$keggq)})
            pathview(gene.data=entrezIDs(), pathway.id=as.matrix(summary(kegg())[1])[1:3],species = "hsa", limit=3)
            system("copy hsa04060.pathview.png www")
            system("copy hsa04064.pathview.png www")
            system("copy hsa04668.pathview.png www")
            #Essayer de mettre la commande pathview dans un tempoutfile : un temp par pathview
            
            output$pathway1 <- renderUI({ tags$img(src="hsa04060.pathview.png", width="400px", height="247px") })
            output$pathway2 <- renderUI({ tags$img(src="hsa04064.pathview.png", width="400px", height="318px") })
            output$pathway3 <- renderUI({ tags$img(src="hsa04668.pathview.png", width="400px", height="273px") })
            
        }
    })
  
  
    ### Protein Domains ############
    
    observe({
        if (input$startpd != 0) {
            ensembl <- reactive({ useEnsembl(biomart="ensembl", dataset="hsapiens_gene_ensembl") })
            univers <- reactive({ getBM(attributes=c('ensembl_gene_id','family'), mart = ensembl()) })
            echantillon <- reactive({ getBM(attributes=c('ensembl_gene_id','family'), filters ='ensembl_gene_id', values = genes(), mart = ensembl()) })
            Univers_uniq_domains <- reactive({ as.matrix(as.matrix(univers()[2])[!duplicated(as.matrix(univers()[2])),]) })
            Univers_domains <- reactive ({ as.matrix(as.matrix(univers()[2])) })
            echantillon_uniq_domains <- reactive({ as.matrix(as.matrix(echantillon()[2])[!duplicated(as.matrix(echantillon()[2])),]) })
            echantillon_domains <- reactive({ as.matrix(echantillon()[2]) })
            Panther_Family_ID <- reactive({ cleanBlankLines(echantillon_uniq_domains()) })
        }
    })
}

shinyApp(ui = ui, server = server)
