library(shiny) # v1.0.1
library(shinythemes) # v 1.1.1
library(ggplot2) # v2.2.1
library(clusterProfiler)# v3.2.11
library(GO.db) # v3.2.0
library(org.Hs.eg.db) # v3.4.0
library(topGO) # v2.26.0
library(DT) # v0.2
library(biomaRt) # v2.30.0
library(MASS) # v 7.3-45

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
                conditionalPanel(
                    condition = "input.startgo !=0",
                    downloadButton('Table_GOclassif', 'Download Table')
                ),
                h4("Barplot from Go classification"),
                plotOutput("groupGoBar"),
                h4("Go enrichment"),
                plotOutput("enrichGo"),
                h4("Dotplot from GO enrichment"),
                plotOutput("enrichGodot"),
                h4("GO graph"),
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
                textInput("browskegg", "ID for browser KEGG :", "Please Enter an ID from the table below"),
                conditionalPanel(
                    condition="input.startkegg != 0",
                    actionButton("startkeggbrows", "Launch browserKEGG"),
                    downloadButton('Table_Kegg', 'Download Table')
                ),
                dataTableOutput("keggtab")
            )
        )
    ),
                 
    tabPanel("Protein domains",
        sidebarLayout(
            sidebarPanel(
                sliderInput(inputId="pdp", label="P-value cut-off", value=0.01, min=0.001, max=0.05 ),
                selectInput(inputId="adjmpd", label="p-value adjustment method", 
                    c("holm", "hochberg", "bonferroni", "BH", "BY", "fdr", "none")),
                radioButtons("metafam", label="Add Families descriptions", choices=c("Yes"="y", "No"="n"), selected="n" ),
                radioButtons("metatransid", label="Add Transcript IDs", choices=c("Yes"="y", "No"="n"), selected="n" ),
                radioButtons("metaensid", label="Add Ensembl Peptides IDs", choices=c("Yes"="y", "No"="n"), selected="n" ),
                    actionButton("startpd", "Start")
                ),
            mainPanel(
                dataTableOutput("Domains"),
                downloadButton('Table_Protein_Domains', 'Download Table'),
                plotOutput("dotplotpd")
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
        ggplot(data(), aes(x=log2FC, y=-log10(adj_pvalue), color=-log10(adj_pvalue) > 4)) +
            geom_point(size = 1.5, alpha = 1, na.rm = T) +
            scale_x_continuous(name="log2FC") +
            scale_y_continuous(name="-log10(adj_pvalue)",limits=c(0,200))
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
          
            # groupGO
          
            gogroup <- reactive({
                print("Go classification ...")
                groupGO(entrezIDs(), OrgDb='org.Hs.eg.db', ont=input$goont, level=input$golvl , readable = TRUE) 
            })
            output$groupGo <- renderDataTable({
                isolate({
                    as.data.frame(gogroup())
                })
            },options=list(pageLength=10), caption="Go classification :")
          
            output$groupGoBar <- renderPlot({
                isolate ({
                    barplot(gogroup(), drop = TRUE, showCategory=input$gocat, title="Barplot from Go classification")
                })
            }) #goenrich
              
            output$Table_GOclassif <- downloadHandler(
                filename = "GO_classification.txt",
                content = function(file) {
                    write.table(as.data.frame(gogroup()), file, sep="\t")
                }
            )
          
            # enrichGo
          
            goenrich <- reactive({
                print("GO enrichment ...")
                enrichGO(entrezIDs(), OrgDb="org.Hs.eg.db", ont=input$goont, pAdjustMethod=input$adjm, pvalueCutoff=input$gop, qvalueCutoff=input$goq, readable= TRUE) 
            })
            enrichDrop <-reactive({ 
                print("Dropping specified GO levels ...")
                drop(goenrich(), input$drop[1], input$drop[2]) 
            })
          
            # Sans isolate, les show category s actualisent tout de suite sans refaire tout le calcul.
            # avec isolate on diminue la reactivite mais tout reprend a chaque fois au debut. reflechir la dessus
            output$enrichGo <- renderPlot({
                isolate ({
                    barplot(enrichDrop(), drop = TRUE, showCategory=input$gocat, title="Go Enrichment")
                })
            }) #goenrich
            output$enrichGodot <- renderPlot({ 
                isolate ({ 
                    dotplot(enrichDrop(), showCategory=input$gocat, title="GO Enrichment")
                })
            }) #goenrich
          
            # svg file plotgograph
            output$gograph <- renderImage({
                isolate({
                    go = enrichDrop() #goenrich
                    width  <- session$clientData$output_gograph_width
                    height <- session$clientData$output_gograph_height
                    outfile <- tempfile(fileext=".svg")
                    svg(outfile)
                    plotGOgraph(go)
                    dev.off()
                    list(src = outfile, contentType = 'image/svg+xml', width = width, height = height, alt = "plotGoGraph")
                })
            }, deleteFile = TRUE)
        }
    })
  
  ### Pathways ###############
    observe({
        if (input$startkegg != 0) {
      
            print("Pathway analysis Started ...")
      
            kegg <- reactive({ enrichKEGG(entrezIDs(), pAdjustMethod=input$adjmk, pvalueCutoff=input$keggp, qvalueCutoff=input$keggq) })
            output$keggtab <- renderDataTable({ 
                as.matrix(as.data.frame(kegg()))
            })
      
            output$Table_Kegg <- downloadHandler(
                filename = "KEGG_enrichment.txt",
                content = function(file) {
                    write.table(as.matrix(as.data.frame(kegg())), file, sep="\t")
                }
            )
      
            if (input$startkeggbrows == 0) {
                return()
            }
            isolate(browseKEGG(keggs, input$browskegg))
        }
    })
  
  ### Protein Domains ########
  
    observe({
        if (input$startpd != 0) {
            # Variables qui ne bougent pas : organism, univers
            print("Building organism and universe from Panther IDs ...")
            ensembl <- useEnsembl(biomart="ensembl", dataset="hsapiens_gene_ensembl")
            univers <- getBM(attributes=c('ensembl_gene_id','family'), mart = ensembl)
            univers_domains <- as.matrix(as.matrix(univers[2])) # Liste Panther ID domaines de l'univers
            univers_uniq_domains <- as.matrix(as.matrix(univers[2])[!duplicated(as.matrix(univers[2])),]) # Pareil, sans doublons
      
            output$Domains <- DT::renderDataTable({
                isolate({
                    withProgress(message = 'Calculation in progress',
                        detail = 'This may take a while...', value = 0, {
                             
                        incProgress(1/9)
                        print("Protein Domains analysis Started ...")
                             
                        incProgress(2/9)
                        # Create sample data
                        print("Building sample data ...")
                        echantillon <- reactive({ 
                            getBM(attributes=c('ensembl_gene_id','family'), filters ='ensembl_gene_id', values = genes(), mart = ensembl)
                        }) 
                        echantillon_uniq_domains <- reactive({ as.matrix(as.matrix(echantillon()[2])[!duplicated(as.matrix(echantillon()[2])),]) })
                        echantillon_domains <- reactive({ as.matrix(echantillon()[2]) }) # pris en entree par test stat
                             
                        incProgress(3/9)
                        # Nettoyage ids vides
                        print("Cleaning empty IDs from sample data ...")
                        Panther_Family_ID <- reactive({ cleanBlankLines(echantillon_uniq_domains()) })
                             
                        incProgress(4/9)
                        # Test hypergeometrique
                        print("Hypergeometric test ...")
                        test_values <- reactive({ pvals_and_hits(Panther_Family_ID(), echantillon_domains(), univers_domains) })
                        pvalue_ajustee <- reactive({ p.adjust(test_values()[[1]], method=input$adjmpd, n=length(univers_domains)) })
                             
                        Domain_enrich_ORA_results <- cbind(Panther_Family_ID(), test_values()[[2]], test_values()[[3]], pvalue_ajustee())
                        colnames(Domain_enrich_ORA_results) <- c("Panther_Family_ID","Nb_Hits","Pvalue","Adj_Pvalue")
                             
                        incProgress(5/9)
                        # Metadata
                        if (input$metafam == "y") {
                            print("Adding Panther families descriptions ...")
                            Panther_Family_Description <- reactive ({ add_metadata('family_description', Panther_Family_ID(), ensembl) })
                            Domain_enrich_ORA_results <- cbind(Domain_enrich_ORA_results, Panther_Family_Description())
                            colnames(Domain_enrich_ORA_results)[length(colnames(Domain_enrich_ORA_results))] <- "Family_description"
                        }
                        incProgress(6/9)
                        if (input$metatransid == "y") {
                            print("Adding Ensemble TRanscripts IDs ...")
                            EnsemblTranscriptID <- reactive ({ add_metadata('ensembl_transcript_id',Panther_Family_ID(), ensembl) })
                            Domain_enrich_ORA_results <- cbind(Domain_enrich_ORA_results, EnsemblTranscriptID())
                                colnames(Domain_enrich_ORA_results)[length(colnames(Domain_enrich_ORA_results))] <- "Ensembl_transcript_ID"
                            }
                        incProgress(7/9)
                        if (input$metaensid == "y") {
                            print("Adding Ensemble Peptides IDs ...")
                            Ensembl_PeptideID <- reactive ({ add_metadata('ensembl_peptide_id',Panther_Family_ID(), ensembl) })
                            Domain_enrich_ORA_results <- cbind(Domain_enrich_ORA_results, Ensembl_PeptideID())
                            colnames(Domain_enrich_ORA_results)[length(colnames(Domain_enrich_ORA_results))] <- "Ensembl_peptide_ID"
                        }
                             
                        incProgress(8/9)
                        # Filter p-values
                        print("Filtering pvalues ...")
                        flags <- reactive({ which(as.numeric(Domain_enrich_ORA_results[,4]) > input$pdp) })
                        Domain_enrich_ORA_results_clean <- reactive ({ as.data.frame(as.matrix(Domain_enrich_ORA_results)[-flags(),]) })
                             
                        incProgress(9/9)
                        # Writing output file
                        print("Writing output file ...")
                        write.matrix(Domain_enrich_ORA_results_clean(), file="DOMAINS_ORA_results.txt", sep="\t")
                        print("Done ! Graphical output is coming ...")
                    })
                })
                # Output
                dataPD <- reactive({ 
                    a <- read.table("DOMAINS_ORA_results.txt", header=TRUE, quote="\"", dec=".", sep="\t", fill=TRUE)
                })
            
                datatable(dataPD(),
                #colnames = c("Panther_Family_ID","Expected_Nb_of_Hits","Pvalue","Adj_Pvalue","Family_Description"),
                options=list(pageLength=20)) #, autoWidth=TRUE
            })
          
            # truc bizarre, a corriger.
            output$Table_Protein_Domains <- downloadHandler(
                filename = "Protein_Domains.txt",
                content = function(file) {
                    a <- read.table("DOMAINS_ORA_results.txt", header=TRUE, quote="\"", dec=".", sep="\t", fill=TRUE)
                    write.table(a, file, sep="\t")
                }
            )
          
            output$dotplotpd <- renderPlot({
                table <- reactive({ read.table("DOMAINS_ORA_results.txt", header=TRUE, quote="\"", dec=".", sep="\t", fill=TRUE) })
                dotplot.domains(table(), genes(), "DotPlot Protein Domains")
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
        pvals[i] <- phyper(
            length(which(echantillon_domains==targeted_domain)),  #nb hits de cet id dans l'echantillon
            length(which(Univers_domains==targeted_domain)), #nb hits de cet id dans l'univers
            length(Univers_domains),
            length(echantillon_domains),
            lower.tail = FALSE
        ) 
        nbhits[i]=length(which(echantillon_domains==targeted_domain))
        expected_nb_hits[i]=(((length(which(Univers_domains==targeted_domain)))/(length(Univers_domains)))*length(echantillon_domains))
    }
    pvhits <- list(pvals, nbhits, expected_nb_hits)
    return(pvhits)
}

# dotplot Protein Domains
dotplot.domains <- function(results, listGenes, title) {
  data<-results[1:10,]
  geneRatio<-data$Pvalue/length(listGenes)
  pal<-colorRampPalette(c("red", "blue"))
  col<-pal(length(unique(results$Pvalue)))[seq(from=1, to = min (length(data), 20))]
  
  par(mar=c(5,max(sapply(rownames(data), nchar))/1.5, 4,6), xpd=TRUE)
  plot(geneRatio, 1:nrow(data)+1, col=col[1], pch=19, yaxt="n", ylab="", xlab="Gene Ratio", main=title)
  axis(2, at=1:nrow(data)+1, label=rownames(data), las=2)
  
  legend("topright", inset=c(-0.2, -0.04), legend =seq(0,1, 0.3), pch = 19, col = pal(length(unique(results$Pvalue)))[sapply(seq(0,1, 0.1), function(level){
    result<-which(round(unique(results$Pvalue), 1)==level)
    if(length(result)==0)return(0)
    else return(result[1])
  })],
  
  title="P value", bty="n")
  legend("topright", inset=c(-0.2, 0.5), legend=unique(data$Nb_Hits), pch=19, col="gray", pt.cex=log10(unique(data$Nb_Hits))+1, bty="n", title="Hit\n number")
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