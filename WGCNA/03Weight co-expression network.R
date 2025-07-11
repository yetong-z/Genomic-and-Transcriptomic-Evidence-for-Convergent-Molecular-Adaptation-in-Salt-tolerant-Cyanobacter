##################### 3.Construct a weighted co-expression network to identify gene modules ####################
#load(file = "step1_input.Rdata")
#load(file = "step2_power_value.Rdata")

cor = WGCNA::cor
if(T){
  net <- blockwiseModules(
    datExpr,
    power = 12, 
    maxBlockSize = ncol(datExpr),
    corType = "pearson", 
    networkType = "unsigned",
    TOMType = "unsigned",
    reassignThreshold = 0,
    minModuleSize = 25, 
    mergeCutHeight = 0.25,
    numericLabels = TRUE, 
    pamRespectsDendro = FALSE,
    saveTOMs = F,
    verbose = 3
  )
  table(net$colors) 
}
cor = stats::cor


if(T){
  # Convert labels to colors for plotting
  moduleColors <- labels2colors(net$colors)
  table(moduleColors)
  # Plot the dendrogram and the module colors underneath
  pdf("step3_genes-modules_ClusterDendrogram.pdf",width = 16,height = 12)
  plotDendroAndColors(net$dendrograms[[1]], moduleColors[net$blockGenes[[1]]],
                      "Module colors",
                      dendroLabels = FALSE, hang = 0.03,
                      addGuide = TRUE, guideHang = 0.05)
  dev.off()
}
save(net, moduleColors, file = "step3_genes_modules.Rdata")
