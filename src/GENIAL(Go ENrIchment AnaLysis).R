library(biomaRt)

####### DRAFT - START

ensembl = useMart("ensembl",dataset="hsapiens_gene_ensembl") #uses human ensembl annotations
#gets gene symbol, transcript_id and go_id for all genes annotated with GO:000750

library(biomaRt)
ensembl = useMart("ensembl",dataset="hsapiens_gene_ensembl") #uses human ensembl annotations
#gets gene symbol, transcript_id and go_id for all genes annotated with GO:0007507
gene.data <- getBM(attributes=c('hgnc_symbol', 'ensembl_transcript_id', 'go_id'),
                   filters = 'go_id', values = 'GO:0007507', mart = ensembl)

capture.output(
  #GO terms retrieving
  #gene.data <- getBM(attributes=c('hgnc_symbol', 'ensembl_transcript_id', 'go_id'),mart = ensembl)
  #gene.data
  
,file="genedata.txt")

ensembl = useMart("ensembl",dataset="hsapiens_gene_ensembl")

results <- getBM(attributes = c("hgnc_symbol","ensembl_transcript_id", "go_id"), filters = "refseq_dna",
                 values = c("ENSG00000137634"), mart = ensembl)

####### DRAFT - END


###### REAL

library(biomaRt)

#use hsapiens DB
ensembl = useMart("ensembl",dataset="hsapiens_gene_ensembl") #uses human ensembl annotations

#map ENSG ids in order to retrieve associated GO ids (add ENSG ids in values = c(ids,ids...) to retrieve these ids)
getBM(attributes = c("hgnc_symbol","ensembl_transcript_id", "go_id"), filters = "ensembl_gene_id",
      values = c("ENSG00000137634"), mart = ensembl)