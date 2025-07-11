############################## 1.Data quality control ################################

### detect missing values
gsg <- goodSamplesGenes(datExpr0,verbose = 3)
gsg$allOK   # If return "True", no missing values
## If have missing values, remove them
if (!gsg$allOK){
  # Optionally, print the gene and sample names that were removed:
  if (sum(!gsg$goodGenes)>0)
    printFlush(paste("Removing genes:", paste(names(datExpr)[!gsg$goodGenes],
                                              collapse = ", ")));
  if (sum(!gsg$goodSamples)>0)
    printFlush(paste("Removing samples:",
                     paste(rownames(datExpr)[!gsg$goodSamples], collapse = ", ")));
  # Remove the offending genes and samples from the data:
  datExpr = datExpr[gsg$goodSamples, gsg$goodGenes]
}
gsg <- goodSamplesGenes(datExpr,verbose = 3)
gsg$allOK

### the systematic clustering tree
if(T){
  sampleTree <- hclust(dist(datExpr), method = "average")
  par(mar = c(0,5,2,0))
  plot(sampleTree, main = "Sample clustering", sub="", xlab="", cex.lab = 1, 
       cex.axis = 1, cex.main = 2)
  sample_colors <- numbers2colors(as.numeric(factor(datTraits$group)), 
                                  colors = rainbow(length(table(datTraits$group))), 
                                  signed = FALSE)   
  par(mar = c(1,4,3,1),cex=0.8) 
  pdf("step1_Sample dendrogram and trait.pdf",width = 8,height = 6)
  plotDendroAndColors(sampleTree, sample_colors,
                      groupLabels = "trait",
                      cex.dendroLabels = 0.8,
                      marAll = c(1, 4, 3, 1),
                      cex.rowText = 0.01,
                      main = "Sample dendrogram and trait" )
  dev.off()
}


if(F){
  clust <- cutreeStatic(sampleTree, cutHeight = 23500, minSize = 10) 
  table(clust)
  keepSamples <- (clust==1)
  datExpr <- datExpr[keepSamples, ]
  datTraits <- datTraits[keepSamples,]
  dim(datExpr) 
}

### PCA
# rm(list = ls())  
# load("step1_input.Rdata")
group_list <- datTraits$group
dat.pca <- PCA(datExpr, graph = F)
pca <- fviz_pca_ind(dat.pca,
                    title = "Principal Component Analysis",
                    legend.title = "Groups",
                    geom.ind = c("point","text"),
                    pointsize = 2,
                    labelsize = 4,
                    repel = TRUE,   # Avoid text overlapping (slow if many points)
                    col.ind = group_list,   # color by groups
                    axes.linetype=NA,    # remove axeslines,
                    mean.point=F
) +
  theme(legend.position = "none")+  # "none": REMOVE legend
  coord_fixed(ratio = 1)
pca
ggsave(pca,filename= "step1_Sample PCA analysis.pdf", width = 8, height = 8)

## save data
nGenes <- ncol(datExpr)
nSamples <- nrow(datExpr)
save(nGenes,nSamples,datExpr,datTraits,file="step1_input.Rdata")
