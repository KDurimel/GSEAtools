#################################### WORKFLOW INITIALISATION : START ####################################

#Set Working directory to source file location  
#this.dir <- dirname(parent.frame(2)$ofile)
#setwd(this.dir)
#source("premierScriptKev.R", chdir = TRUE)
####### First Lib requierements
library(ggplot2) # v2.2.1
###############################


#Raise number of printable strings in order to be able to capture large outputs
options(max.print=99999)


####### LOG FILE AND REPORT FILE INITIALISATION ####### 

## Date ##
date=as.character.Date(date())
date=gsub(" ", "_", date, fixed = TRUE)
date=gsub(":", "-", date, fixed = TRUE)



####### Create GSEA Analysis directory and log/report files #######
dirname = paste("../results/",date,sep="") #from shiny
dir.create(dirname)

## Log file ##
logfilename=paste("../results/",date,"/Genial_",date,".log",sep="") #file name format Genial_Day_Month_Hour_min_sec_year.log
#Initialize log file 
write(date,file=logfilename,append=FALSE) 

 ## Report file ##
reportfilename=logfilename=paste("../results/",date,"/Genial_",date,".txt",sep="") #file name format Genial_Day_Month_Hour_min_sec_year.log
#Initialize report file
write(date,file=reportfilename,append=FALSE)



            ################################################################
######################           SHINY INITIALISATION               ######################
            ################################################################


# Shiny Server initialisation

   #Shiny Victor#

#Progress bar
pb<-winProgressBar(title="Step 1/x - Initialising workflow", label="0% done", min=0, max=100, initial=0)
getWinProgressBar(pb)

## Set Workflow parameters : FROM SHINY STDINPUT ##

inputfile="../raw_data/dataGO/human_retine.txt" # Fichier d'entree
Log2FC_column = 3 # Colonne du log2FC dans fichier d'entree
Log2FC_cutoff = 2 # Cutoff du log2FC
RnaPvalue_column = 5 # Colonne p-valeur ajustee
Qvalue_cutoff = 0.05 # Qvalue max pour ajustement p-valeur enrichissement GO
#RnaPvalue_cutoff = 0.05 # Cutoff p-valeur ajustee : #maj : pas besoin car fichier deja filtré
IdsType="ENSEMBL" #Type d'identifiants du fichier d'entrée (prendre dans COLONNE deroulante shiny)
MinGOlevel=1 #Niveau minimal a garder dans l'enrichissement GO
MaxGOlevel=6 #Niveau maximal a garder dans l'enrichissement GO

setWinProgressBar(pb,50,title="Step 1/x - Initialising workflow",label="50% done")

#Input file
data <- read.table(inputfile, header=TRUE, dec=",") 
setWinProgressBar(pb,100,title="Step 1/x - Initialising workflow",label="100% done")
close(pb)

#################################### WORKFLOW INITIALISATION : END ####################################
#=====================================================================================================#
#=====================================================================================================#
#######################################   GSEA-GO : START   ###########################################

pb<-winProgressBar(title="Step 2/x - GSEA analysis", label="0% done \t Generating vulcano plots ...",
                   min=0, max=100, initial=0)
getWinProgressBar(pb)


#P_value : variable temporaire pour Coloration conditionelle selon p-value
Log10_P_value <- cut(as.matrix(data[RnaPvalue_column]),
                     breaks = c(-Inf,1.410000e-23, Inf), # threshold x=1.410000e-23 soit valeur de log(x)
                     labels = c("<=25", ">25"))

#Plot -log10(p-value)=f(Log2FC)
png(file = paste("../results/",date,"/VulcanoPlot.png",sep=""), bg = "transparent", width = 800 , height= 600) # Img output settings
fig1<-ggplot(data, aes(x=data[Log2FC_column], y=-log10(data[RnaPvalue_column]),color=Log10_P_value)) + 
  geom_point(size = 1.5, alpha = 0.7, na.rm = T) +
  scale_color_manual(values = c("black", "cyan")) +
  scale_x_continuous(name="Log2FC") + 
  scale_y_continuous(name="-Log10(Pvalue)",limits=c(0,200))

plot(fig1)

setWinProgressBar(pb,50,
                  title="Step 2/x - GSEA analysis",label="50% done \t Generating vulcano plots ...")
dev.off() # close stdout device and export output 

plot(fig1) #rmarkdown



#Plot p-value=f(Log2FC) for eDge
png(file = paste("../results/",date,"/Pval_adjEDGE.png",sep=""), bg = "transparent", width = 800 , height= 600) # Img output settings
fig2<-ggplot(data, aes(x=log2FC_edge, y=-log10(data[Pvalue_adjcolumnEDGE]),color=P_value_Edge)) + 
  geom_point(size = 1.5, alpha = 0.7, na.rm = T) +
  scale_color_manual(values = c("black", "cyan")) +
  scale_x_continuous(name="log2FC_edge") + 
  scale_y_continuous(name="p_edge",limits=c(0,200))
plot(fig2)
dev.off() # close stdout device and export output
plot(fig2) #rmarkdown

setWinProgressBar(pb,100,
                  title="Step 2/x - GSEA analysis",label="100% done")
close(pb)


          ################################################################
######################            GO-ENRICHMENT               ######################
          ################################################################

pb<-winProgressBar(title="Step 3/x - GO enrichment", label="0% done \t Loading required libraries",
                   min=0, max=100, initial=0)
getWinProgressBar(pb)



################## LOAD Second lib Requirements ################
library(clusterProfiler) # v3.2.11
library(GO.db) # v2.9
library(AnnotationDbi) # v 1.34.4
library(org.Hs.eg.db) # v3.4
library(shiny) # v 1.0.0
################################################################

setWinProgressBar(pb,10,
                  title="Step 3/x - GO enrichment",label="10% done \t Generating vulcano plots ...")

#Extract all ENSEMBL ids from input : p-values are already filtered
#genes <- data$Ensembl_Gene_ID[abs(data$log2FC_DeSeq)>Log2FC_cutoff] filtrer log2FC
#genes <- as.list(data[Idscolumn])
setWinProgressBar(pb,15,
                  title="Step 3/x - GO enrichment",label="15% done \t Mapping dataset ...")
genes <- data$Ensembl_Gene_ID
#genes <- as.character(genes)


#Human ENTREZID ids retrieving from ENSEMBL ids : use bitr (biological id translator)
setWinProgressBar(pb,20,
                  title="Step 3/x - GO enrichment",label="20% done \t Retrieving ids ...")
EntrezIDs = bitr(genes, fromType= IdsType , toType="ENTREZID", OrgDb="org.Hs.eg.db")$ENTREZID # $ENTREZID in order to keep entrezID column

## Retrieve GO-TERMS from ENTREZID ##
setWinProgressBar(pb,25,
                  title="Step 3/x - GO enrichment",label="25% done \t Retrieving GO-terms ...")
#BP Ontology
BP <- groupGO(gene = EntrezIDs, 'org.Hs.eg.db', ont = "BP", level = 7)
#CC Ontology
#CC <- groupGO(gene = EntrezIDs, 'org.Hs.eg.db', ont = "CC", level = 7) 
#MF Ontology
#MF <- groupGO(gene = EntrezIDs, 'org.Hs.eg.db', ont = "MF", level = 7) 

setWinProgressBar(pb,100,
                  title="Step 3/x - GO enrichment",label="100% done \t")
close(pb)

            ################################################################
######################      OVER-REPRESENTATION ANALYSIS       ######################
            ################################################################

pb<-winProgressBar(title="Step 4/x - GO ORA Analysis", label="0% done \t BP Enrichment ...",
                   min=0, max=100, initial=0)
getWinProgressBar(pb)

#genes<-data$Ensembl_Gene_ID[abs(data$log2FC_DeSeq)>2]
#genes <- as.character(genes)
#, qvalueCutoff=0.5
#pvalue adjust : si 0.01 --> correspond a que j'accepte jusqu'a 1% de faux positifs dans ma liste

###################################
### ORA analysis
###################################

#BP Ontology
BPenrich <- enrichGO(EntrezIDs , OrgDb = "org.Hs.eg.db", ont = "BP", 
                     pAdjustMethod="fdr", pvalueCutoff=0.05, qvalueCutoff = 0.05, 
                     readable= TRUE)

setWinProgressBar(pb,35,
                  title="Step 4/x - GO ORA Analysis",label="35% done \t CC Enrichment ...")

#CC Ontology
CCenrich <- enrichGO(EntrezIDs , OrgDb = "org.Hs.eg.db", ont = "CC",
                     pAdjustMethod="fdr", pvalueCutoff=0.05, qvalueCutoff = 0.05,
                     readable= TRUE)

setWinProgressBar(pb,60,
                  title="Step 4/x - GO ORA Analysis",label="60% done \t MF Enrichment ...")

#MF Ontology
MFenrich <- enrichGO(EntrezIDs , OrgDb = "org.Hs.eg.db", ont = "MF", 
                     pAdjustMethod="fdr", pvalueCutoff=0.05, qvalueCutoff = 0.05,
                     readable= TRUE)

setWinProgressBar(pb,85,
                  title="Step 4/x - GO ORA Analysis",label="85% done \t Trimming terms ...")

## Drop too mainstream terms ##

  # future variables without dropped terms 
BPenrichDrop<-BPenrich 
CCenrichDrop<-CCenrich
MFenrichDrop<-MFenrich

  #Drop enriched terms between level min and max thresholds
for(i in MinGOlevel:MaxGOlevel)
{
  BPenrichDrop <- dropGO(BPenrichDrop,level=i)
  CCenrichDrop <- dropGO(CCenrichDrop,level=i)
  MFenrichDrop <- dropGO(MFenrichDrop,level=i)
}

setWinProgressBar(pb,95,
                  title="Step 4/x - GO ORA Analysis",label="95% done \t save results in tabulated file")


################################################
### Save trimmed ORA results in tabulated file
################################################
#BP
#capture.output(
#   {
#     summary(BPenrichDrop)[2:9]
#   },file=(paste("../results/",date,"/BPresults.txt",sep="")))
# #CC
# capture.output(
#   {
#     summary(CCenrichDrop)[2:9]
#   },file=(paste("../results/",date,"/CCresults.txt",sep="")))
# #MF
# capture.output(
#   {
#     summary(MFenrichDrop)[2:9]
#   },file=(paste("../results/",date,"/MFresults.txt",sep="")))



#VISUALISER LES RESULTATS D'ENRICHISSEMENT DE LA VARIABLE

  #View(summary(BPenrichDrop)) <---INTERFACER DANS SHINY!!!



#SORTIR LES RESULTATS D'ENRICHISSEMENT DANS UN FICHIER TABULE

header=paste("GO:id","description","Ontology","GeneRatio","BackgroundRatio","Pvalue","Adjusted_Pvalue","Qvalue","Hits","GeneIds",sep = "\t")
write(header,file=paste("../results/",date,"/GOresults.txt",sep=""),append=FALSE) 

Goline="" #Go:id , description, etc ... for each line
capture.output(
{
  for(i in 1:length(as.matrix(summary(BPenrichDrop)[1]))) #Number of enriched term : nb of summary lines for a criteria
  {
    Goline[i]=paste(as.matrix(summary(BPenrichDrop)[1])[i],as.matrix(summary(BPenrichDrop)[2])[i],"BP",
          as.matrix(summary(BPenrichDrop)[3])[i],as.matrix(summary(BPenrichDrop)[4])[i],
          as.matrix(summary(BPenrichDrop)[5])[i],as.matrix(summary(BPenrichDrop)[6])[i],
          as.matrix(summary(BPenrichDrop)[7])[i],as.matrix(summary(BPenrichDrop)[9])[i],
          as.matrix(summary(BPenrichDrop)[8])[i],sep="\t")
    cat(Goline[i])
    cat("\n")
  }
  
  for(i in 1:length(as.matrix(summary(CCenrichDrop)[1]))) #Number of enriched term : nb of summary lines for a criteria
  {
    Goline[i]=paste(as.matrix(summary(CCenrichDrop)[1])[i],as.matrix(summary(CCenrichDrop)[2])[i],"CC",
                    as.matrix(summary(CCenrichDrop)[3])[i],as.matrix(summary(CCenrichDrop)[4])[i],
                    as.matrix(summary(CCenrichDrop)[5])[i],as.matrix(summary(CCenrichDrop)[6])[i],
                    as.matrix(summary(CCenrichDrop)[7])[i],as.matrix(summary(CCenrichDrop)[9])[i],
                    as.matrix(summary(CCenrichDrop)[8])[i],sep="\t")
    cat(Goline[i])
    cat("\n")
  }
  
  for(i in 1:length(as.matrix(summary(MFenrichDrop)[1]))) #Number of enriched term : nb of summary lines for a criteria
  {
    Goline[i]=paste(as.matrix(summary(MFenrichDrop)[1])[i],as.matrix(summary(MFenrichDrop)[2])[i],"MF",
                    as.matrix(summary(MFenrichDrop)[3])[i],as.matrix(summary(MFenrichDrop)[4])[i],
                    as.matrix(summary(MFenrichDrop)[5])[i],as.matrix(summary(MFenrichDrop)[6])[i],
                    as.matrix(summary(MFenrichDrop)[7])[i],as.matrix(summary(MFenrichDrop)[9])[i],
                    as.matrix(summary(MFenrichDrop)[8])[i],sep="\t")
    cat(Goline[i])
    cat("\n")
  }
},file=(paste("../results/",date,"/GOresults.txt",sep="")),append=TRUE)




setWinProgressBar(pb,100,
                  title="Step 4/x - GO ORA Analysis",label="100% done \t ok")
close(pb)

           ##################################################################
################# OVER REPRESENTATION ANALYSYS : graphical outputs  ####################
           ##################################################################


pb<-winProgressBar(title="Step 5/x - GO ORA Analysis charts", label="0% done \t BP Barplot and dotplots...",
                   min=0, max=100, initial=0)
getWinProgressBar(pb)


## BARplots and dotplots ##

#BP
png(file = paste("../results/",date,"/barplot_BP.png",sep=""), bg = "transparent", width = 800 , height= 600) # Img output settings
barplot(BPenrich, drop = TRUE,showCategory = 12) #Show the 12 most enriched terms
dev.off()
barplot(BPenrich, drop = TRUE, showCategory = 12) #rmarkdown

png(file = paste("../results/",date,"/dotplot_BP.png",sep=""), bg = "transparent", width = 800 , height= 600) # Img output settings
dotplot(BPenrich, showCategory = 12) #Show the 12 most enriched terms
dev.off()
dotplot(BPenrich, showCategory = 12) #rmarkdown

setWinProgressBar(pb,20,
                  title="Step 5/x - GO ORA Analysis charts", label="20% done \t CC Barplot and dotplots...")

#CC
png(file = paste("../results/",date,"/barplot_CC.png",sep=""), bg = "transparent", width = 800 , height= 600) # Img output settings
barplot(CCenrich, drop = TRUE, showCategory = 12) #Show the 12 most enriched terms
dev.off()
barplot(CCenrich, drop = TRUE, showCategory = 12) #Show the 12 most enriched terms

png(file = paste("../results/",date,"/dotplot_CC.png",sep=""), bg = "transparent", width = 800 , height= 600) # Img output settings
dotplot(CCenrich, showCategory = 12) #Show the 12 most enriched terms
dev.off()
dotplot(CCenrich, showCategory = 12) #Show the 12 most enriched terms

setWinProgressBar(pb,35,
                  title="Step 5/x - GO ORA Analysis charts", label="35% done \t MF Barplot and dotplots...")


#MF
png(file = paste("../results/",date,"/barplot_MF.png",sep=""), bg = "transparent", width = 800 , height= 600) # Img output settings
barplot(BPenrich, drop = TRUE, showCategory = 12) #Show the 12 most enriched terms
dev.off()
barplot(BPenrich, drop = TRUE, showCategory = 12) #Show the 12 most enriched terms

png(file = paste("../results/",date,"/dotplot_MF.png",sep=""), bg = "transparent", width = 800 , height= 600) # Img output settings
dotplot(BPenrich, showCategory = 12) #Show the 12 most enriched terms
dev.off()
dotplot(BPenrich, showCategory = 12) #Show the 12 most enriched terms


## Tree plot ##

setWinProgressBar(pb,50,
                  title="Step 5/x - GO ORA Analysis charts", label="50% done \t BP Tree plot...")

#BP : 12 most significant nodes
svg(file = paste("../results/",date,"/Treeplot_BP.svg",sep=""))
#plotGOgraph(BPenrich, firstSigNodes = 12, useInfo = "all", sigForAll = TRUE, useFullNames = TRUE)
plotGOgraph(BPenrich, useInfo = "all", sigForAll = TRUE, useFullNames = TRUE)
dev.off()
#plotGOgraph(BPenrich, firstSigNodes = 12, useInfo = "all", sigForAll = TRUE, useFullNames = TRUE)
plotGOgraph(BPenrich, useInfo = "all", sigForAll = TRUE, useFullNames = TRUE)

setWinProgressBar(pb,55,
                  title="Step 5/x - GO ORA Analysis charts", label="55% done \t CC Tree plot...")

#CC
svg(file = paste("../results/",date,"/Treeplot_CC.svg",sep=""))
plotGOgraph(CCenrich, useInfo = "all", sigForAll = TRUE, useFullNames = TRUE)
dev.off()
plotGOgraph(CCenrich, useInfo = "all", sigForAll = TRUE, useFullNames = TRUE)

setWinProgressBar(pb,60,
                  title="Step 5/x - GO ORA Analysis charts", label="60% done \t MF Tree plot...")

#MF
svg(file = paste("../results/",date,"/Treeplot_MF.svg",sep=""))
plotGOgraph(MFenrich, useInfo = "all", sigForAll = TRUE, useFullNames = TRUE)
dev.off()
plotGOgraph(MFenrich, useInfo = "all", sigForAll = TRUE, useFullNames = TRUE)

## Enrichment map plot ##

setWinProgressBar(pb,65,
                  title="Step 5/x - GO ORA Analysis charts", label="65% done \t BP enrichment plot...")

#BP
png(file = paste("../results/",date,"/Enrichmentplot_BP.png",sep=""), bg = "transparent", width = 800 , height= 600)
enrichMap(BPenrich)
dev.off()
enrichMap(BPenrich)

setWinProgressBar(pb,70,
                  title="Step 5/x - GO ORA Analysis charts", label="70% done \t CC enrichment plot...")
#CC
png(file = paste("../results/",date,"/Enrichmentplot_CC.png",sep=""), bg = "transparent", width = 800 , height= 600)
enrichMap(CCenrich)
dev.off()
enrichMap(CCenrich)


setWinProgressBar(pb,75,
                  title="Step 5/x - GO ORA Analysis charts", label="75% done \t MF enrichment plot...")

#MF
png(file = paste("../results/",date,"/Enrichmentplot_MF.png",sep=""), bg = "transparent", width = 800 , height= 600)
enrichMap(MFenrich)
dev.off()
enrichMap(MFenrich)


## Cnet plot : reseau des relations entre GO-termes enrichis et les genes associes a celui-ci ##

setWinProgressBar(pb,80,
                  title="Step 5/x - GO ORA Analysis charts", label="80% done \t BP cnet plot...")

#BP
png(file = paste("../results/",date,"/Relations_BP.png",sep=""), bg = "transparent", width = 800 , height= 600)
cnetplot(BPenrich, categorySize = "pvalue") # categorySize can be scaled by 'pvalue' or 'geneNum'
dev.off()
cnetplot(BPenrich, categorySize = "pvalue") # categorySize can be scaled by 'pvalue' or 'geneNum'

setWinProgressBar(pb,85,
                  title="Step 5/x - GO ORA Analysis charts", label="85% done \t CC cnet plot...")

#CC
png(file = paste("../results/",date,"/Relations_CC.png",sep=""), bg = "transparent", width = 800 , height= 600)
cnetplot(CCenrich, categorySize = "pvalue")
dev.off()
cnetplot(CCenrich, categorySize = "pvalue")

setWinProgressBar(pb,90,
                  title="Step 5/x - GO ORA Analysis charts", label="90% done \t MF cnet plot...")

#MF
png(file = paste("../results/",date,"/Relations_MF.png",sep=""), bg = "transparent", width = 800 , height= 600)
cnetplot(MFenrich, categorySize = "pvalue")
dev.off()
cnetplot(MFenrich, categorySize = "pvalue")

## Dynamic HTML chart ##

setWinProgressBar(pb,95,
                  title="Step 5/x - GO ORA Analysis charts", label="95% done \t Dynamic chart...")
  
  
  
  MFdescriptions <- as.matrix(summary(MFenrichDrop)[2])
  ##convertir pvalues en coeffs exploitable : avec log(1/pvalue) : plus p-val faible (significatif), plus coeff sera élevé
  MFcoeffs <- as.matrix(summary(MFenrichDrop)[5]) #p-value contenu dans colonne[5]
  MFcoeffs <- log(1/MFcoeffs) #convertir p-value en coefficient exploitable dans highcharts
  
  ##highcharts##
  #generer code higchart
  for(i in 1:length(MFdescriptions))
  {
    #generer code highchart pour mfdescr
  }


#######################################   GSEA-GO : END   #############################################
#=====================================================================================================#
#=====================================================================================================#
#####################################   GSEA-KEGG : START   ###########################################


#Close progesss bar
close(pb)
```

library(pathview)

keggs=enrichKEGG(EntrezIDs ,pAdjustMethod="fdr", pvalueCutoff=0.05, qvalueCutoff = 0.05)

viewKegg=pathview(gene.data = EntrezIDs, pathway.id = as.matrix(summary(keggs)[1])[1:5],species = "hsa", limit = list(gene = length(EntrezIDs)))



#######################################   GSEA-KEGG : END   ###########################################
#=====================================================================================================#
#=====================================================================================================#
#####################################   GSEA-PROTFAM : START  #########################################

pvals=""
for(i in 1:length(Go_uniq))
{
  targeted_goterm=toString(Go_uniq[i]) #tostring in order to request value
  #Bro-tip= lower.tail=TRUE because (lower.tail=TRUE)=FALSE in next loop, for obscure reasons
  pvals[i]<-phyper(length(which(Go_ech==targeted_goterm)),length(which(Go_univ==targeted_goterm)),length(Go_univ[,1]),length(Go_ech[,1]),lower.tail = FALSE) 
}
pvals2<-as.vector(pvals)
pvals_adjusted<-p.adjust(pvals2, method="bonferroni", n=length(pvals))
a=1

######PHYPER OK, NOW ENHANCE PHYPER RESULTS#####
capture.output(
  #FILE Header
  cat("GO:ID","Go term","Number of hits","Expected number of hits","Go level","P-value","Corrected p-value",sep=";"),
  cat("\n"),
  
  for(i in 1:length(Go_uniq))
  {
    targeted_goterm=toString(Go_uniq[i]) #tostring in order to request value
    #Targeted_goterm: GOterm currently processed##Number of hits##Expected number of hits##phyper(count of this GO in the batch,count of this go in the universe, universe length,batch length)
    #With:number of hits: count of this goterm in the batch
    #expected number of hits: count of the same goterm in the universe
    expectednumberofhits=(((length(which(Go_univ==targeted_goterm)))/(length(Go_univ[,1])))*length(Go_ech[,1]))
    trm=Term(targeted_goterm)
    ##########################################################################
    #RETRIEVE GO-TERM LEVEL : by parents "walking" method (is-a , part-of)
    ##########################################################################
    #First request
    parents=get(targeted_goterm,GOMFPARENTS) 
    #If parents founded-->next DAG level exists.So repeat request with the parents of the next DAG level
    if(parents[1]!=parent)
    {
      level=1 #level in the DAG
      tryCatch({
        for(i in 1:1000) #Arbitrary.
        {
          parents=get(parents[1],GOMFPARENTS) #[1] because prioritize for "is-a"
          level=level+1
          if(parents[1]==parent) stop()#If node ="all" , all go three processed,no more DAG levels.Stop the requests loop
        }
      }, error=function(e){cat()})#Nothing in message error in order to dont print msgerros in results file
    }
    #If no parents founded-->term is on highest level. This is Aspect term-->lvl
    else
    {
      level=1
    }
    #Print results foreach go-term
    if(length(which(Go_ech==targeted_goterm))!=0) #if go-term is present in sample
    {
      col3=as.character(length(which(Go_ech==targeted_goterm)))
      col6=as.character(phyper(length(which(Go_ech==targeted_goterm)),length(which(Go_univ==targeted_goterm)),length(Go_univ[,1]),length(Go_ech[,1]),lower.tail = FALSE))
      cat(targeted_goterm,trm,col3,expectednumberofhits,level,col6,pvals_adjusted[a],sep=";")
      cat("\n")
    }
    rm(targeted_goterm)
    a=a+1
  }
  ,file=paste(compvalue,"/hyperesults_mf.txt",sep=""))
