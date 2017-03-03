# Essai avec un des jeux de données sur l'human qu'on avait trouvé.
# Les reference a tous les articles qu'on a trouvé sont restés dans les favoris
# firefox de mon PC à la fac ._.

# Requirements
library(ggplot2)
library(clusterProfiler)
library(GO.db)
library(org.Hs.eg.db)

# Lecture des donnees
data <- read.table("human_retine.txt", header=TRUE, dec=",")
head(data)

# Plots :
gradu_x <- c(-10,-9,-8,-7,-6,-5,-4,-3,-2,-1,0,1,2,3,4,5,6,7,8,9,10)
gradu_y <- c(0,0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1)

#dev.off() : j'avais un bug avec les graphes - dev.off() réintialise tout
ggplot(data, aes(x=log2FC_DeSeq, y=p_DeSeq)) + geom_point(shape=1, color="blue", size=1) +
  scale_x_continuous(name="log2FC_DeSeq", breaks=gradu_x) + 
  scale_y_continuous(name="p_DeSeq")
ggplot(data, aes(x=log2FC_edge, y=p_edge)) + geom_point(shape=16, color="blue", size=1) + 
  scale_x_continuous(name="log2FC_edge", breaks=gradu_x) + 
  scale_y_continuous(name="p_edge")

# Enrich GO
genes<-data$Ensembl_Gene_ID[abs(data$log2FC_DeSeq)>2]
ego <- enrichGO(genes , organism= "human", ont= "CC", pAdjustMethod="BH", pvalueCutoff=0.1, qvalueCutoff=0.5, readable= TRUE)
# ego est vide -> RTFW