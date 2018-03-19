#################################### WORKFLOW INITIALISATION : START ####################################

# Remember : this is a draft script --> Best practices level ~ 0%

#Set Working directory to source file location  
#this.dir <- dirname(parent.frame(2)$ofile)
#setwd(this.dir)

####### External dependencies
library(ggplot2)
library(clusterProfiler) 
library(GO.db) 
library(AnnotationDbi) 
library(org.Hs.eg.db) 
library(biomaRt)
library(shiny)
library(pathview) 
library(MASS)

start.time <- Sys.time() #Init basic counter

#Raise number of printable strings in order to be able to catch large outputs
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

#Progress bar
pb<-winProgressBar(title="Step 1/7 - Initialising workflow", label="0% done", min=0, max=100, initial=0)
getWinProgressBar(pb)

## Set Workflow parameters
inputfile=file.choose()# Fichier d'entree
Log2FC_column = 3 # Colonne du log2FC dans fichier d'entree
Log2FC_cutoff = 2 # Cutoff du log2FC
RnaPvalue_column = 5 # Colonne p-valeur ajustee
Qvalue_cutoff = 0.05 # Qvalue max pour ajustement p-valeur enrichissement GO
#RnaPvalue_cutoff = 0.05 # Cutoff p-valeur ajustee : #maj : pas besoin car fichier deja filtr
IdsType="ENSEMBL" #Type d'identifiants du fichier d'entre (prendre dans COLONNE deroulante shiny)
MinGOlevel=7 #Niveau minimal a garder dans l'enrichissement GO
MaxGOlevel=14 #Niveau maximal a garder dans l'enrichissement GO
TopKEGG = 3 # Nombre de Pathways max a plotter : default = 3
Padjustmethod="fdr"
pvaluecutoff=0.05 #pvalue utilisee dans les tests

setWinProgressBar(pb,50,title="Step 1/7 - Initialising workflow",label="50% done")

#Input file
data <- read.table(inputfile, header=TRUE, dec=",") 
setWinProgressBar(pb,100,title="Step 1/7 - Initialising workflow",label="100% done")
close(pb)

#################################### WORKFLOW INITIALISATION : END ####################################
#=====================================================================================================#
#=====================================================================================================#
#######################################   GSEA-GO : START   ###########################################

pb<-winProgressBar(title="Step 2/7 - GSEA analysis", label="0% done \t Generating vulcano plots ...",
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
                  title="Step 2/7 - GSEA analysis",label="50% done \t Generating vulcano plots ...")
dev.off() # close stdout device and export output 

setWinProgressBar(pb,100,
                  title="Step 2/7 - GSEA analysis",label="100% done")
close(pb)


          ################################################################
######################            GO-ENRICHMENT               ######################
          ################################################################

pb<-winProgressBar(title="Step 3/7 - GO enrichment", label="0% done \t Loading required libraries",
                   min=0, max=100, initial=0)

getWinProgressBar(pb)

setWinProgressBar(pb,10,
                  title="Step 3/7 - GO enrichment",label="10% done \t Generating vulcano plots ...")

#Extract all ENSEMBL ids from input : p-values are already filtered
setWinProgressBar(pb,15,
                  title="Step 3/7 - GO enrichment",label="15% done \t Mapping dataset ...")
genes <- data$Ensembl_Gene_ID

#Human ENTREZID ids retrieving from ENSEMBL ids : use bitr (biological id translator)
setWinProgressBar(pb,20,
                  title="Step 3/7 - GO enrichment",label="20% done \t Retrieving ids ...")
EntrezIDs = bitr(genes, fromType= IdsType , toType="ENTREZID", OrgDb="org.Hs.eg.db")$ENTREZID # $ENTREZID in order to keep entrezID column

## Retrieve GO-TERMS from ENTREZID ##
setWinProgressBar(pb,25,
                  title="Step 3/7 - GO enrichment",label="25% done \t Retrieving GO-terms ...")

setWinProgressBar(pb,100,
                  title="Step 3/7 - GO enrichment",label="100% done \t")
close(pb)

            ################################################################
######################      OVER-REPRESENTATION ANALYSIS       ######################
            ################################################################

pb<-winProgressBar(title="Step 4/7 - GO ORA Analysis", label="0% done \t BP Enrichment ...",
                   min=0, max=100, initial=0)
getWinProgressBar(pb)


###################################
### ORA analysis
###################################

#BP Ontology
BPenrich <- enrichGO(EntrezIDs , OrgDb = "org.Hs.eg.db", ont = "BP", 
                     pAdjustMethod="fdr", pvalueCutoff=0.05, qvalueCutoff = 0.05, 
                     readable= TRUE)

setWinProgressBar(pb,35,
                  title="Step 4/7 - GO ORA Analysis",label="35% done \t CC Enrichment ...")

#CC Ontology
CCenrich <- enrichGO(EntrezIDs , OrgDb = "org.Hs.eg.db", ont = "CC",
                     pAdjustMethod="fdr", pvalueCutoff=0.05, qvalueCutoff = 0.05,
                     readable= TRUE)

setWinProgressBar(pb,60,
                  title="Step 4/7 - GO ORA Analysis",label="60% done \t MF Enrichment ...")

#MF Ontology
MFenrich <- enrichGO(EntrezIDs , OrgDb = "org.Hs.eg.db", ont = "MF", 
                     pAdjustMethod="fdr", pvalueCutoff=0.05, qvalueCutoff = 0.05,
                     readable= TRUE)

setWinProgressBar(pb,85,
                  title="Step 4/7 - GO ORA Analysis",label="85% done \t Trimming terms ...")

## Drop too generic terms according some  thresholds ##


BPenrichDrop<-BPenrich # without dropped (too generic) terms 
CCenrichDrop<-CCenrich # without dropped (too generic) terms 
MFenrichDrop<-MFenrich # without dropped (too generic) terms 


#Drop terms below left threshold
for(i in 1:MinGOlevel)
{
  BPenrichDrop <- dropGO(BPenrichDrop,level=i)
  CCenrichDrop <- dropGO(CCenrichDrop,level=i)
  MFenrichDrop <- dropGO(MFenrichDrop,level=i)
}
#Drop terms behind right threshold
for(i in MaxGOlevel:14)
{
  BPenrichDrop <- dropGO(BPenrichDrop,level=i)
  CCenrichDrop <- dropGO(CCenrichDrop,level=i)
  MFenrichDrop <- dropGO(MFenrichDrop,level=i)
}

setWinProgressBar(pb,95,
                  title="Step 4/7 - GO ORA Analysis",label="95% done \t save results in tabulated file")


#Save results in a TSV file

header=paste("GO:id","description","Ontology","GeneRatio","BackgroundRatio","Pvalue","Adjusted_Pvalue","Qvalue","Hits","GeneIds",sep = "\t")
write(header,file=paste("../results/",date,"/GOresults.tsv",sep=""),append=FALSE) 

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
},file=(paste("../results/",date,"/GOresults.tsv",sep="")),append=TRUE)




setWinProgressBar(pb,100,
                  title="Step 4/7 - GO ORA Analysis",label="100% done \t ok")
close(pb)

           ##################################################################
################# OVER REPRESENTATION ANALYSYS : graphical outputs  ####################
           ##################################################################


pb<-winProgressBar(title="Step 5/7 - GO ORA Analysis charts", label="0% done \t BP Barplot and dotplots...",
                   min=0, max=100, initial=0)
getWinProgressBar(pb)


## BARplots and dotplots ##

#BP
png(file = paste("../results/",date,"/barplot_BP.png",sep=""), bg = "transparent", width = 800 , height= 600) # Img output settings
barplot(BPenrich, drop = TRUE,showCategory = 12) #Show the 12 most enriched terms
dev.off()

png(file = paste("../results/",date,"/dotplot_BP.png",sep=""), bg = "transparent", width = 800 , height= 600) # Img output settings
dotplot(BPenrich, showCategory = 12) #Show the 12 most enriched terms
dev.off()

setWinProgressBar(pb,20,
                  title="Step 5/7 - GO ORA Analysis charts", label="20% done \t CC Barplot and dotplots...")

#CC
png(file = paste("../results/",date,"/barplot_CC.png",sep=""), bg = "transparent", width = 800 , height= 600) # Img output settings
barplot(CCenrich, drop = TRUE, showCategory = 12) #Show the 12 most enriched terms
dev.off()

png(file = paste("../results/",date,"/dotplot_CC.png",sep=""), bg = "transparent", width = 800 , height= 600) # Img output settings
dotplot(CCenrich, showCategory = 12) #Show the 12 most enriched terms
dev.off()


setWinProgressBar(pb,35,
                  title="Step 5/7 - GO ORA Analysis charts", label="35% done \t MF Barplot and dotplots...")


#MF
png(file = paste("../results/",date,"/barplot_MF.png",sep=""), bg = "transparent", width = 800 , height= 600) # Img output settings
barplot(BPenrich, drop = TRUE, showCategory = 12) #Show the 12 most enriched terms
dev.off()


png(file = paste("../results/",date,"/dotplot_MF.png",sep=""), bg = "transparent", width = 800 , height= 600) # Img output settings
dotplot(BPenrich, showCategory = 12) #Show the 12 most enriched terms
dev.off()



## Tree plot ##

setWinProgressBar(pb,50,
                  title="Step 5/7 - GO ORA Analysis charts", label="50% done \t BP Tree plot...")

#BP : 12 most significant nodes
svg(file = paste("../results/",date,"/Treeplot_BP.svg",sep=""))
plotGOgraph(BPenrich, useInfo = "all", sigForAll = TRUE, useFullNames = TRUE)
dev.off()


setWinProgressBar(pb,55,
                  title="Step 5/7 - GO ORA Analysis charts", label="55% done \t CC Tree plot...")

#CC
svg(file = paste("../results/",date,"/Treeplot_CC.svg",sep=""))
plotGOgraph(CCenrich, useInfo = "all", sigForAll = TRUE, useFullNames = TRUE)
dev.off()

setWinProgressBar(pb,60,
                  title="Step 5/7 - GO ORA Analysis charts", label="60% done \t MF Tree plot...")

#MF
svg(file = paste("../results/",date,"/Treeplot_MF.svg",sep=""))
plotGOgraph(MFenrich, useInfo = "all", sigForAll = TRUE, useFullNames = TRUE)
dev.off()


## Enrichment map plot ##

setWinProgressBar(pb,65,
                  title="Step 5/7 - GO ORA Analysis charts", label="65% done \t BP enrichment plot...")

#BP
png(file = paste("../results/",date,"/Enrichmentplot_BP.png",sep=""), bg = "transparent", width = 800 , height= 600)
enrichMap(BPenrich)
dev.off()

setWinProgressBar(pb,70,
                  title="Step 5/7 - GO ORA Analysis charts", label="70% done \t CC enrichment plot...")
#CC
png(file = paste("../results/",date,"/Enrichmentplot_CC.png",sep=""), bg = "transparent", width = 800 , height= 600)
enrichMap(CCenrich)
dev.off()


setWinProgressBar(pb,75,
                  title="Step 5/7 - GO ORA Analysis charts", label="75% done \t MF enrichment plot...")

#MF
png(file = paste("../results/",date,"/Enrichmentplot_MF.png",sep=""), bg = "transparent", width = 800 , height= 600)
enrichMap(MFenrich)
dev.off()


## Cnet plot : reseau des relations entre GO-termes enrichis et les genes associes a celui-ci ##

setWinProgressBar(pb,80,
                  title="Step 5/7 - GO ORA Analysis charts", label="80% done \t BP cnet plot...")

#BP
png(file = paste("../results/",date,"/Relations_BP.png",sep=""), bg = "transparent", width = 800 , height= 600)
cnetplot(BPenrich, categorySize = "pvalue") # categorySize can be scaled by 'pvalue' or 'geneNum'
dev.off()

setWinProgressBar(pb,85,
                  title="Step 5/7 - GO ORA Analysis charts", label="85% done \t CC cnet plot...")

#CC
png(file = paste("../results/",date,"/Relations_CC.png",sep=""), bg = "transparent", width = 800 , height= 600)
cnetplot(CCenrich, categorySize = "pvalue")
dev.off()

setWinProgressBar(pb,90,
                  title="Step 5/7 - GO ORA Analysis charts", label="90% done \t MF cnet plot...")

#MF
png(file = paste("../results/",date,"/Relations_MF.png",sep=""), bg = "transparent", width = 800 , height= 600)
cnetplot(MFenrich, categorySize = "pvalue")
dev.off()

## Dynamic HTML chart ##

setWinProgressBar(pb,95,
                  title="Step 5/7 - GO ORA Analysis charts", label="95% done \t Dynamic chart...")
  
  
  
  MFdescriptions <- as.matrix(summary(MFenrichDrop)[2])
  ##convertir pvalues en coeffs exploitable : avec log(1/pvalue) : plus p-val faible (significatif), plus coeff sera lev
  MFcoeffs <- as.matrix(summary(MFenrichDrop)[5]) #p-value contenu dans colonne[5]
  MFcoeffs <- log(1/MFcoeffs) #convertir p-value en coefficient exploitable dans highcharts
  
  ##highcharts##
  for(i in 1:length(MFdescriptions))
  {
    #todo
  }
  
  #Close progesss bar
  close(pb)

#######################################   GSEA-GO : END   #############################################
#=====================================================================================================#
#=====================================================================================================#
#####################################   GSEA-KEGG : START   ###########################################



pb<-winProgressBar(title="Step 6/7 -  KEGG enrichment", label="0% done \t Loading required libraries", min=0, max=100, initial=0)
getWinProgressBar(pb)
  


setWinProgressBar(pb,15,
                  title="Step 6/7 -  KEGG enrichment", label="15% done \t EnrichKEGG...")

keggs=enrichKEGG(EntrezIDs ,pAdjustMethod="fdr", pvalueCutoff=pvaluecutoff, qvalueCutoff = Qvalue_cutoff)

setWinProgressBar(pb,50,
                  title="Step 6/7 -  KEGG enrichment", label="50% done \t Grapical representation")


setWinProgressBar(pb,99,
                  title="Step 6/7 -  KEGG enrichment", label="95% done \t Write results")

write.matrix(summary(keggs),file=paste("../results/",date,"/KEGGresults.tsv",sep=""),sep="\t")

close(pb)

#######################################   GSEA-KEGG : END   ###########################################
#=====================================================================================================#
#=====================================================================================================#
#####################################   GSEA-PROTFAM : START  #########################################

pb<-winProgressBar(title="Step 7/7 -  Proteins Domains enrichment", label="0% done \t Loading required libraries", min=0, max=100, initial=0)
getWinProgressBar(pb)


##############################################################################
### Generate sample and universe datasets before hypergeometric test
##############################################################################

#Set organism
ensembl = useEnsembl(biomart="ensembl", dataset="hsapiens_gene_ensembl")

## Donwload universe
setWinProgressBar(pb,5,
                  title="Step 7/7 -  Proteins Domains enrichment", label="5% done \t Download universe...")

univers <- getBM(attributes=c('ensembl_gene_id','family'), mart = ensembl)

setWinProgressBar(pb,15,
                  title="Step 7/7 -  Proteins Domains enrichment", label="15% done \t Downlaod sample metadata...")

## Create Sample data
echantillon <- getBM(attributes=c('ensembl_gene_id','family'), filters ='ensembl_gene_id',
                                      values = genes, mart = ensembl)


setWinProgressBar(pb,20,
                  title="Step 7/7 -  Proteins Domains enrichment", label="20% done \t Setting and cleaning data...")

## Liste de tous les Panther ID domaines sans doublons
Univers_uniq_domains <- as.matrix(as.matrix(univers[2])[!duplicated(as.matrix(univers[2])),]) 
## Liste de tous les Panther ID domaines de l'univers (= avec doublons)
Univers_domains <- as.matrix(as.matrix(univers[2]))
## Liste de tous les Panther ID domaines de l'echantillon sans doublons (utile pour generation tableau resultats et test stat non redondant sur un id)
echantillon_uniq_domains <- as.matrix(as.matrix(echantillon[2])[!duplicated(as.matrix(echantillon[2])),]) # liste des domaines a parcourir
## Liste de tous les Panther ID domaines de l'echatillon avec doublons (pris en entree par test stat)
echantillon_domains=as.matrix(echantillon[2])




cleanBlankLines<- function(datatoclean) #nettoyer les ids vides renvoys par BIOCONDUCTOR
{                                       #prend un entre un objet getBM() 
  count=1
  datacleaned=""
  for(i in which(as.vector(datatoclean)!=""))
  {
    datacleaned[count]=as.matrix(datatoclean)[i]
    count=count+1
  }
  return(datacleaned)
}

Panther_Family_ID=cleanBlankLines(echantillon_uniq_domains) #nettoyage ids vides

##############################################################################
### TEST STATISTIQUE : test hypergeometrique
##############################################################################

setWinProgressBar(pb,25,
                  title="Step 7/7 -  Proteins Domains enrichment", label="25% done \t Hypergeometric test...")

pvals="" #stockage p-values
nbhits=""
expectednumberofhits=""


#calcul p-values et nbhits
for(i in 1:length(Panther_Family_ID))
{
  targeted_domain = Panther_Family_ID[i] #id recherch
  
  #enrichment tests
  pvals[i] <- phyper(length(which(echantillon_domains==targeted_domain)),  #nb hits de cet id dans l'echantillon
                     length(which(Univers_domains==targeted_domain)), #nb hits de cet id dans l'univers
                     length(Univers_domains),
                     length(echantillon_domains),
                     lower.tail = FALSE) 
  nbhits[i]=length(which(echantillon_domains==targeted_domain))
  expectednumberofhits[i]=(((length(which(Univers_domains==targeted_domain)))/(length(Univers_domains)))*length(echantillon_domains))
}


##############################################################################
### TEST STATISTIQUE : correction des p-values
##############################################################################

pValue_ajustee <- p.adjust(pvals, method=Padjustmethod, n=length(Univers_domains))  # valeur du n a verifier


##############################################################################
### AJOUT DE METADONNES AU TEST STATS pour generation fichier de resultats
##############################################################################

#objectif : remplir metadatas[x] avec les id pfam correspondant aux id panther family a une position x
#remplir de blank lines si aucun id ne mappe, afin de conserver la mme dimension que echantillon_uniq_domains[x]
#pourquoi? l'objectif est de fusionner ces variables afin d'obtenir un tableau 

setWinProgressBar(pb,30,
                  title="Step 7/7 -  Proteins Domains enrichment", label="30% done \t Add supplementary metadata...")

add_metadata<- function(bioCdatabase) #retourne metadonnes d'une database (arg bioCdatabase) mapps sur ids panther family de l'echantillon
{                                     #Prend en entre un rsultat de cleanBlankLines()
  metadata=""
  for(i in 1:length(Panther_Family_ID))
  {
    #if domain id exists : do id retrieving
    if(Panther_Family_ID[i]!="")
    {
      metadata[i]=getBM(attributes=c(bioCdatabase),
                        filters ='family',
                        values = Panther_Family_ID[i],
                        mart = ensembl) 
    }
    #if domain id is blank (bug happened somewhere) : add blank line to metadata[i] : useless since cleanblanklines impl.
    else
    {
      metadata[i]="-"
    }
  }
  return(metadata)
}

# Metadonnee 1 : description de la famille proteique
setWinProgressBar(pb,45,
                  title="Step 7/7 -  Proteins Domains enrichment", label="45% done \t Add supplementary metadata...family description")

Panther_Family_Description=add_metadata('family_description')


# Metadonnee 2 : ensembl transcript ID
setWinProgressBar(pb,55,
                  title="Step 7/7 -  Proteins Domains enrichment", label="55% done \t Add supplementary metadata...Ensembl transcript ID")

EnsemblTranscriptID=add_metadata('ensembl_transcript_id')


# Metadonnee 3 : ensembl peptide ID
setWinProgressBar(pb,65,
                  title="Step 7/7 -  Proteins Domains enrichment", label="65% done \t Add supplementary metadata...Ensembl peptide ID")

EnsemblPeptideID=add_metadata('ensembl_peptide_id')


# Metadonnee 4 : pfam id
setWinProgressBar(pb,75,
                  title="Step 7/7 -  Proteins Domains enrichment", label="75% done \t Add supplementary metadata...Pfam")

PfamID=add_metadata('pfam')


# Metadonnee 5 : pfam start
setWinProgressBar(pb,85,
                  title="Step 7/7 -  Proteins Domains enrichment", label="85% done \t Add supplementary metadata...Pfam S")

PfamSTART=add_metadata('pfam_start')


# Metadonnee 6 : pfam end
setWinProgressBar(pb,95,
                  title="Step 7/7 -  Proteins Domains enrichment", label="95% done \t Add supplementary metadata...Pfam E")

PfamEND=add_metadata('pfam_end')

setWinProgressBar(pb,98,
                  title="Step 7/7 -  Proteins Domains enrichment", label="98% done \t Settings...")


# Merge results and metatadas inside a Matrix
Domain_enrich_results = cbind(Panther_Family_ID,
                              Panther_Family_Description,
                              EnsemblTranscriptID,
                              EnsemblPeptideID,
                              PfamID,
                              PfamSTART,
                              PfamEND)

# Merge ORA results and metatadas inside a Matrix
Domain_enrich_ORA_results = cbind(Panther_Family_ID,
                              nbhits,
                              expectednumberofhits,
                              pValue_ajustee,
                              Panther_Family_Description,
                              EnsemblTranscriptID,
                              EnsemblPeptideID,
                              PfamID,
                              PfamSTART,
                              PfamEND)


# flag p-value > pvaluecutoff
flags=which(as.numeric(Domain_enrich_ORA_results[,4]) > pvaluecutoff) 

# new variable without lines with p-value > pvaluecutoff
Domain_enrich_ORA_results_clean=as.matrix(Domain_enrich_ORA_results)[-flags,] 

setWinProgressBar(pb,99,
                  title="Step 7/7 -  Proteins Domains enrichment", label="98% done \t Write results...")

#Afficher l'enrichissement en domaines
write.matrix(Domain_enrich_results,file=paste("../results/",date,"/DOMAINSresults.tsv",sep=""),sep="\t") #MASS write.matrix used

write.matrix(Domain_enrich_ORA_results_clean,file=paste("../results/",date,"/DOMAINS_ORA_results.tsv",sep=""),sep="\t") #MASS write.matrix used

close(pb)

end.time <- Sys.time()
time.taken <- end.time - start.time

winDialog(type = c("ok"), paste("Analysis finished with success.\n","Execution time : ",time.taken,sep=""))
