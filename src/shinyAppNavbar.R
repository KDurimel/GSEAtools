library(shiny) # v1.0.0
library(ggplot2) # v2.2.1
library(clusterProfiler)# v3.2.11
library(GO.db) # v2.9
library(org.Hs.eg.db) # v3.4.0
library(topGO) # v2.26.0
library(pathview) # v1.14
library(DT)
library(biomaRt)
library(MASS)

ui <- navbarPage(theme=shinytheme("sandstone"), "Gene Set Enrichment Analysis",
           
  tabPanel("Data",
    sidebarLayout(
      sidebarPanel(
        fileInput("data", "Choose a txt file", multiple=FALSE, accept = NULL)
      ),
      mainPanel(
        dataTableOutput("contents")
      )
    )
  ),
  
  tabPanel("Vulcano Plot",
    plotOutput("vulcano")
  ),
  
  tabPanel("GO enrichment",
    sidebarLayout(
      sidebarPanel(
        radioButtons("dogroup", label="DO Go classification ?", choices=c("Yes"="y", "No"="n"), selected="n"),
        numericInput(inputId="golvl", label="Level for go classification", value=7, min=1, max=14, step=1),
        sliderInput(inputId="gop", label="P-value cut-off", value=0.01, min=0.001, max=0.05 ),
        sliderInput(inputId="goq", label="q-value cut-off", value=0.01, min=0.001, max=0.05 ),
        sliderInput(inputId="gocat", label="categories to show", value=12, min=8, max=18),
        sliderInput(inputId="drop", label="Keep GO levels between", min = 1, max = 14, value = c(1,14)),
        selectInput(inputId="goont", label="Ontology", c("CC", "BP", "MF")),
        selectInput(inputId="adjm", label="p-value adjustment method", 
                    c("holm", "hochberg", "hommel", "bonferroni", "BH", "BY", "fdr", "none")),
        actionButton("startgo", "Start")
      ),
      mainPanel(
        dataTableOutput("groupGo"), 
        plotOutput("enrichGo"), 
        plotOutput("enrichGodot"),
        plotOutput("gograph")
      )
    )
  ),
  
  tabPanel("Pathways",
    sidebarLayout(
      sidebarPanel(
        sliderInput(inputId="keggp", label="P-value cut-off", value=0.01, min=0.001, max=0.05 ),
        sliderInput(inputId="keggq", label="q-value cut-off", value=0.01, min=0.001, max=0.05 ),
        selectInput(inputId="adjmk", label="p-value adjustment method", 
                    c("holm", "hochberg", "hommel", "bonferroni", "BH", "BY", "fdr", "none")),
        actionButton("startkegg", "Start")
      ),
      mainPanel(
        uiOutput("pathway1"),
        uiOutput("pathway2"),
        uiOutput("pathway3")
      )
    )
  ),
  
  tabPanel("Protein domains",
    sidebarLayout(
      sidebarPanel(
        sliderInput(inputId="pdp", label="P-value cut-off", value=0.01, min=0.001, max=0.05 ),
        selectInput(inputId="adjmpd", label="p-value adjustment method", 
                    c("holm", "hochberg", "hommel", "bonferroni", "BH", "BY", "fdr", "none")),
        actionButton("startpd", "Start")
      ),
      mainPanel(
        dataTableOutput("Domains")
      )
    )
  )
)

server <- function(input, output, session) {
  ### Data ###################
  data <- reactive({
    if (is.null(input$data)) {return(NULL)}
    read.table(input$data$datapath, header = TRUE, dec=",", sep="\t")
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
  
  ### IDs ####################
  genes <- reactive({
    as.character(data()$Ensembl_Gene_ID)
  })
  entrezIDs <- reactive({
    bitr(genes(), fromType="ENSEMBL" , toType="ENTREZID", OrgDb="org.Hs.eg.db")$ENTREZID
  })
  
  ### GO #####################
  observe({
    if(input$startgo != 0) {
      print("Go Analysis started...")
      # groupgo
      if(input$dogroup=="y") {
        print("Go classification ...")
        gogroup <- reactive({ 
          groupGO(entrezIDs(), OrgDb='org.Hs.eg.db', ont=input$goont, level=input$golvl , readable = TRUE) 
        })
        output$groupGo <- renderDataTable(as.data.frame(gogroup()), options=list(pageLength=10))
      }
      
      # enrichGo
      goenrich <- reactive({
        print("GO enrichment ...")
        enrichGO(entrezIDs(), OrgDb="org.Hs.eg.db", ont=input$goont, pAdjustMethod=input$adjm, pvalueCutoff=input$gop, qvalueCutoff=input$goq, readable= TRUE) 
      })
      enrichDrop <-reactive({ 
        print("Dropping specified GO levels ...")
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
  
  ### Pathways ###############
  observe({
    if (input$startkegg != 0) {
      print("Pathway analysis Started ...")
      kegg <- reactive({ enrichKEGG(entrezIDs(), pAdjustMethod=input$adjmk, pvalueCutoff=input$keggp, qvalueCutoff=input$keggq) })
      pathview(gene.data=genes(), pathway.id=as.matrix(as.data.frame(kegg())[1])[1:3],species = "hsa", limit=c(1,3))
      
      shell("copy *.png www")
      
      output$pathway1 <- renderUI({ tags$img(src="hsa04060.pathview.png", width="400px", height="247px") })
      output$pathway2 <- renderUI({ tags$img(src="hsa04064.pathview.png", width="400px", height="318px") })
      output$pathway3 <- renderUI({ tags$img(src="hsa04668.pathview.png", width="400px", height="273px") })
    }
  })
  
  ### Protein Domains ########
  
  observe({
    if (input$startpd != 0) {
      output$Domains <- DT::renderDataTable({
        withProgress(message = 'Calculation in progress',
          detail = 'This may take a while...', value = 0, {
                       
            incProgress(1/7)
            print("Protein Domains analysis Started ...")
            print("Building organism and universe from Panther IDs ...")
            # Variables qui ne bougent pas : organism, univers
            ensembl <- useEnsembl(biomart="ensembl", dataset="hsapiens_gene_ensembl")
            univers <- getBM(attributes=c('ensembl_gene_id','family'), mart = ensembl)
            univers_domains <- as.matrix(as.matrix(univers[2])) # Liste Panther ID domaines de l'univers
            univers_uniq_domains <- as.matrix(as.matrix(univers[2])[!duplicated(as.matrix(univers[2])),]) # Pareil, sans doublons
                       
            incProgress(2/7)
            # Create sample data
            print("Building sample data ...")
            echantillon <- reactive({ 
              getBM(attributes=c('ensembl_gene_id','family'), filters ='ensembl_gene_id', values = genes(), mart = ensembl)
            }) 
            echantillon_uniq_domains <- reactive({ as.matrix(as.matrix(echantillon()[2])[!duplicated(as.matrix(echantillon()[2])),]) })
            echantillon_domains <- reactive({ as.matrix(echantillon()[2]) }) # pris en entree par test stat
                       
            incProgress(3/7)
            # Nettoyage ids vides
            print("Cleaning empty IDs from sample data ...")
            Panther_Family_ID <- reactive({ cleanBlankLines(echantillon_uniq_domains()) })
                       
            incProgress(4/7)
            # Test hypergeometrique
            print("Hypergeometric test ...")
            test_values <- reactive({ pvals_and_hits(Panther_Family_ID(), echantillon_domains(), univers_domains) })
            pvalue_ajustee <- reactive({ p.adjust(test_values()[[1]], method=input$adjmpd, n=length(univers_domains)) })
                       
            incProgress(5/7)
            # Metadata
            print("Adding metadata ...")
            Panther_Family_Description <- reactive ({ add_metadata('family_description', Panther_Family_ID(), ensembl) })
                       
            incProgress(6/7)
            # Final table v1
            Domain_enrich_ORA_results <- cbind(Panther_Family_ID(), test_values()[[2]], test_values()[[3]], pvalue_ajustee(), Panther_Family_Description())
                       
            # Filter p-values
            print("Filtering pvalues ...")
            flags <- reactive({ which(as.numeric(Domain_enrich_ORA_results[,4]) > input$pdp) })
            Domain_enrich_ORA_results_clean <- reactive ({ as.data.frame(as.matrix(Domain_enrich_ORA_results)[-flags(),]) })
                       
            incProgress(7/7)
            # Writing output file
            print("Writing output file ...")
            write.matrix(Domain_enrich_ORA_results_clean(), file="DOMAINS_ORA_results.txt", sep="\t")
                       
            print("Done ! Graphical output is coming ...")
          })
        
        # Output
        dataPD <- reactive({ read.table("DOMAINS_ORA_results.txt", header=TRUE, quote="\"", dec=".", sep="\t") })
        
        datatable(dataPD(),
                  colnames = c("Panther_Family_ID","Expected_Nb_of_Hits","Pvalue","Adjusted_Pvalue","Family_Description"),
                  options=list(pageLength=20)) #, autoWidth=TRUE
      })
    } # end if input$startpd != 0
  }) # end Protein Domains
}

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
add_metadata <- function(bioCdatabase, Panther_Family_ID, ensembl) {
  metadata=""
  for(i in 1:length(Panther_Family_ID)) { ## CA NE PEUT PAS MARCHER COMME CA, PASSER EN PARAMETRE
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

# Test hypergeometrique
pvals_and_hits <- function(Panther_Family_ID, echantillon_domains, Univers_domains) {
  pvals = ""
  nbhits = ""
  expected_nb_hits = ""
  for(i in 1:length(Panther_Family_ID)) {
    
    targeted_domain = Panther_Family_ID[i] #id recherch
    
    #enrichment tests
    pvals[i] <- phyper(length(which(echantillon_domains==targeted_domain)),  #nb hits de cet id dans l'echantillon
                       length(which(Univers_domains==targeted_domain)), #nb hits de cet id dans l'univers
                       length(Univers_domains),
                       length(echantillon_domains),
                       lower.tail = FALSE) 
    nbhits[i]=length(which(echantillon_domains==targeted_domain))
    expected_nb_hits[i]=(((length(which(Univers_domains==targeted_domain)))/(length(Univers_domains)))*length(echantillon_domains))
  }
  pvhits <- list(pvals, nbhits, expected_nb_hits)
  return(pvhits)
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

shinyApp(ui = ui, server = server)