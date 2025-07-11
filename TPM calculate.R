library(ggplot2)
library(DESeq2)
getwd()
setwd("E:/RNAseqResult/")
txdb <- makeTxDbFromGFF("FACHB-2_genomic.gff", format="gff")
gene <- exonsBy(txdb,by="gene")
gene_lens <- lapply(gene, function(x){sum(width(reduce(x)))})
# check
class(gene_lens)
length(gene_lens)
# change to data frame
gene_lens1 <- as.data.frame(gene_lens)
class(gene_lens1)
dim(gene_lens1)

gene_lens1 <- t(gene_lens1)
dim(gene_lens1)

## change names
colnames(gene_lens1) <- c("Geneid","Length")
write.csv(exons_gene_lens1, file="gene_Length.csv")

readscounts <- read.table("E:/featurecount.txt", header=T)
readscounts <- readcounts+1

## merge
mycounts <- merge(exons_gene_lens2, readscounts, by="Geneid", all=F)
dim(mycounts)

## set row1 as row name
rownames(mycounts)<-mycounts[,1]
mycounts<-mycounts[,-1]

## kb calculate
kb <- mycounts$Length / 1000
write.csv(kb, file="kb.csv")

countdata <- mycounts[,2:7]
## rpk
rpk <- countdata / kb

## TPM
tpm <- t(t(rpk)/colSums(rpk)*1000000)
head(tpm)
write.csv(tpm, file="transcript.tpm.matrix.csv")
## FPKM
fpkm <- t(t(rpk)/colSums(countdata)*10^6) 
head(fpkm)
write.csv(fpkm, file="mRNA_fpkm.csv")


