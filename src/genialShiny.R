library(ggplot2) # v2.2.1
library(clusterProfiler)# v3.2.11
library(GO.db) # v2.9
library(org.Hs.eg.db) # v3.4
library(shiny)

##############
###  Page ###
##############

ui <- fluidPage(
  titlePanel("Gene Set Enrichment Analysis"),
  
    fileInput("data", "Choose a txt file", multiple=FALSE,
              accept = NULL
    ),
    mainPanel(
      
      tabsetPanel(
        # Tableau de données
        tabPanel("Data", id="data", dataTableOutput("contents")),
        
  ####################################
        
        # Onglet formulaire GO
        tabPanel("GO analysis", id="fgo",
                 
          checkboxGroupInput("fctsgo", "Choose what to run", 
                            c("GO classification"="groupGO", 
                              "Over-representation test"="enrichGO", 
                              "GO Gene Set Enrichment Analysis"="GOGSEA")),
          
          sliderInput(inputId="gop", label="P-value cut-off", value=0.01, min=0.001, max=0.05 ),
          
          sliderInput(inputId="goq", label="q-value cut-off", value=0.01, min=0.001, max=0.05 ),
          
          #sliderInput(IinputId="gocat", label="categories to show", value=12, min=8, max=18),
          
          selectInput(inputId="goont", label="Ontology", c("CC", "BP", "MF"))),
        
  ####################################
                   
        # Onglet formulaire KEGG
        tabPanel("Pathways analysis", id="fkegg",
                 
          checkboxGroupInput("fctskegg", "Choose what to run",
                              c("KEGG over-representation test"="enrichkegg",
                                "KEGG gene set enrichment"="kegggesea")),
          
          sliderInput(inputId="keggp", label="P-value cut-off", value=0.01, min=0.001, max=0.05 ),
          
          sliderInput(inputId="keggq", label="q-value cut-off", value=0.01, min=0.001, max=0.05 ),
          
          #sliderInput(IinputId="keggcat", label="categories to show", value=12, min=8, max=18),
          
          selectInput(inputId="goont", label="Ontology", c("CC", "BP", "MF"))),
        
  ####################################
        
        #Onglet formulaire Protein Domains
        tabPanel("ProteinDomains", id="fpd", inputPanel("formpd"))),
        
  ####################################
      tabsetPanel(
        tabPanel("Plots", id="plots", plotOutput("vulcano")),
        tabPanel("GO results", id="rgo", plotOutput("enrichGo")),
        tabPanel("Pathways results", id="rkegg", plotOutput("enrichKegg")),
        tabPanel("Protein Domains results", id="rpd", "Protein Domains")
      )
    )
  )

###############
###  SERVER ###
###############


server <- function(input, output) {
  # Lit le jeu de donnees
  data <- reactive({
    read.table(input$data$datapath, header = TRUE, dec=",")
  })
  # Affiche le jeu de donnees
  output$contents <- renderDataTable(data(), options=list(pageLength=10))
  
  # Vulcano plot
  output$vulcano <- renderPlot({
    ggplot(data(), aes(x=log2FC, y=-log10(adj_pvalue))) +
      geom_point(size = 1.5, alpha = 0.7, na.rm = T) +
      scale_x_continuous(name="log2FC") +
      scale_y_continuous(name="-log10(p_value)",limits=c(0,200))
  })
  
  genes <- reactive({ genes <- data()$Ensembl_Gene_ID })
  entrezIDs <- reactive({ bitr(genes(), fromType= IdsType , toType="ENTREZID", OrgDb="org.Hs.eg.db")$ENTREZID})
  
  # enrichGo
  goenrich <- reactive({
    enrichGO(entrezIDs(), OrgDb="org.Hs.eg.db", ont=input$goont, pAdjustMethod="fdr", pvalueCutoff=input$gop, qvalueCutoff=input$goq, readable= TRUE)
  })
  output$enrichgo <- renderPlot({barplot(goenrich(), drop = TRUE, showCategory = 12)}) # ajouter un param showcategory
}

shinyApp(ui = ui, server = server)
