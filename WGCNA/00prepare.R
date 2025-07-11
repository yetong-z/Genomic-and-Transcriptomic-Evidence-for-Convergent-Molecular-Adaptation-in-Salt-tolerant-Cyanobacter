rm(list = ls())

install.packages("WGCNA")
install.packages("FactoMineR")
install.packages("factoextra")
BiocManager::install("GO.db") 
BiocManager::install("impute")
BiocManager::install("preprocessCore")

options(stringsAsFactors = F)

Sys.setenv(LANGUAGE = "en") ## Set language
library(WGCNA)
library(FactoMineR)
library(factoextra)
library(tidyverse) # ggplot2 stringer dplyr tidyr readr purrr  tibble forcats
library(data.table) 
library(GO.db)
library(dplyr)


enableWGCNAThreads(nThreads = 0.75*parallel::detectCores())

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

################################# 0.prepare input data ################################
# get data
tpm00 <- read.csv("transcript.tpm.matrix.csv", header = T)
rownames(tpm00) <- tpm00[,1]
tpm00 <- tpm00[,-1] 
tpm <- tpm00
data <- log2(tpm+1)
## MADtop5000
keep_data <- data[order(apply(data,1,mad), decreasing = T)[1:5000],]
datExpr0 <- as.data.frame(t(keep_data))
datExpr <- datExpr0
## create datTraits，including information of group, traits, etc
datTraits <- data.frame(row.names = colnames(tpm),group=colnames(tpm))
fix(datTraits)
grouptype <- data.frame(group=sort(unique(datTraits$group)),
                        groupNo=1:length(unique(datTraits$group)))
fix(grouptype)
datTraits$groupNo = "NA"
for(i in 1:nrow(grouptype)){
  datTraits[which(datTraits$group == grouptype$group[i]),'groupNo'] <- grouptype$groupNo[i]
}
datTraits
