# Essai avec un des jeux de données sur l'human qu'on avait trouvé.
# Les reference a tous les articles qu'on a trouvé sont restés dans les favoris
# firefox de mon PC à la fac ._.

# Y'a des trucs ci dessous qui marchent, mais je sais toujours pas si on part du bon
# jeu de données (genre celui là, il y a des p-value et des log2FC, mais dans les data
# du papier conseillé par Dauchel, y'en a pas ...)

# Requirements
library(ggplot2) # v2.2.1
library(clusterProfiler)# v3.2.11
library(GO.db) # v2.9
library(org.Hs.eg.db) # v3.4

# Lecture des donnees
data <- read.table("human_retine.txt", header=TRUE, dec=",")
head(data)

# Plots :
  # Graduations des axes
gradu_x <- c(-10,-9,-8,-7,-6,-5,-4,-3,-2,-1,0,1,2,3,4,5,6,7,8,9,10)
gradu_y <- c(0,0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1)

#dev.off() : j'avais un bug avec les graphes - dev.off() réintialise tout
ggplot(data, aes(x=log2FC_DeSeq, y=p_DeSeq)) + 
  geom_point(shape=1, color="blue", size=1) +
  scale_x_continuous(name="log2FC_DeSeq", breaks=gradu_x) + 
  scale_y_continuous(name="p_DeSeq")
ggplot(data, aes(x=log2FC_edge, y=p_edge)) + 
  geom_point(shape=16, color="blue", size=1) + 
  scale_x_continuous(name="log2FC_edge", breaks=gradu_x) + 
  scale_y_continuous(name="p_edge")

# Gene Ontology classification
genes <- data$Ensembl_Gene_ID
genes <- as.character(genes)
head(genes)
eg = bitr(genes, fromType="ENSEMBL", toType="ENTREZID", OrgDb="org.Hs.eg.db")
eg <- eg$ENTREZID
ggo <- groupGO(gene = eg, 'org.Hs.eg.db', ont = "BP", level = # Conversion des ID ENSEMBL au type ENTREZID, 
               readable = TRUE)

### Pas encore re-testé :

# Enrich GO
genes<-data$Ensembl_Gene_ID[abs(data$log2FC_DeSeq)>2]
genes <- as.character(genes)
ego <- enrichGO(genes , organism= "human", ont= "CC", pAdjustMethod="BH", pvalueCutoff=0.1, qvalueCutoff=0.5, readable= TRUE)
# ego est vide -> RTFW

# A voir
barplot(ggo, drop = TRUE, showCategory = 12)
barplot(ego, showCategory = 8)
enrichMap(ego)
cnetplot(ego, categorySize = "pvalue", foldChange = geneList)
cnetplot(kk, categorySize = "geneNum", foldChange = geneList)