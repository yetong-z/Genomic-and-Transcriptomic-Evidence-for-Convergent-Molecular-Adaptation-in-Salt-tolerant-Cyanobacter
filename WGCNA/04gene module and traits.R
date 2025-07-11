####################### 4.Gene modules and phenotypes #####################################
# rm(list = ls())  
# load(file = "step1_input.Rdata")
# load(file = "step2_power_value.Rdata")
# load(file = "step3_genes_modules.Rdata")

if(T){
  datTraits$group <- relevel(as.factor(datTraits$group), "Normal")
  design <- model.matrix(~0+datTraits$group)  
  colnames(design) <- levels(datTraits$group) 
  rownames(design) <- rownames(datTraits)
  MES0 <- moduleEigengenes(datExpr, moduleColors)$eigengenes  #Calculate module eigengenes
  MEs <- orderMEs(MES0)  #Put close eigenvectors next to each other
  moduleTraitCor <- cor(MEs, design, use = "p")
  moduleTraitPvalue <- corPvalueStudent(moduleTraitCor, nSamples)
  sizeGrWindow(10,6)
  textMatrix <- paste0(signif(moduleTraitCor,2),"\n(",
                       signif(moduleTraitPvalue,1),")", sep = "")
  dim(textMatrix) <- dim(moduleTraitCor)
  pdf("step4_Module-trait-relationship_heatmap.pdf",
      width = 2*length(colnames(design)), 
      height = 0.6*length(names(MEs)) )
  par(mar=c(5, 9, 3, 3))
  labeledHeatmap(Matrix = moduleTraitCor,
                 xLabels = colnames(design),
                 yLabels = names(MEs),
                 ySymbols = names(MEs),
                 colorLabels = F,
                 colors = blueWhiteRed(50),
                 textMatrix = textMatrix,
                 setStdMargins = F,
                 cex.text = 0.5,
                 zlim = c(-1,1), 
                 main = "Module-trait relationships")
  dev.off()
}

save(design, file = "step4_design.Rdata")

### boxplot
if(T){
  mes_group <- merge(MEs,datTraits,by="row.names") 
  library(gplots)
  library(ggpubr)
  library(grid)
  library(gridExtra) 
  draw_ggboxplot <- function(data,Module="Module",group="group"){
    ggboxplot(data,x=group, y=Module,
              ylab = paste0(Module),
              xlab = group,
              fill = group,
              palette = "jco",
              #add="jitter",
              legend = "") +
      stat_compare_means()+
      theme(axis.text.x = element_text(angle = 45, hjust = 1))
  }
  colorNames <- names(MEs)
  pdf("step4_Module-trait-relationship_boxplot.pdf", width = 7.5,height = 2.1*ncol(MEs))
  p <- lapply(colorNames,function(x) {
    draw_ggboxplot(mes_group, Module = x, group = "group")
  })
  do.call(grid.arrange,c(p,ncol=2)) 
  dev.off()
}


levels(datTraits$group)
choose_group <- "Recover"  

if(T){
  modNames <- substring(names(MEs), 3)
  geneModuleMembership <- as.data.frame(cor(datExpr, MEs, use = "p"))
  MMPvalue <- as.data.frame(corPvalueStudent(as.matrix(geneModuleMembership), nSamples))
  names(geneModuleMembership) <- paste0("MM", modNames)
  names(MMPvalue) <- paste0("p.MM", modNames)
  trait <- as.data.frame(design[,choose_group])
  geneTraitSignificance <- as.data.frame(cor(datExpr, trait, use = "p"))
  GSPvalue <- as.data.frame(corPvalueStudent(as.matrix(geneTraitSignificance), nSamples))
  names(geneTraitSignificance) <- paste0("GS")
  names(GSPvalue) <- paste0("GS")
  selectModule <- modNames  
  pdf("step4_gene-Module-trait-significance.pdf", width = 7, height = 2.0*ncol(MEs))
  par(mfrow = c(ceiling(length(selectModule)/2), 2)) 
  for(module in selectModule){
    column <- match(module, selectModule)
    print(module)
    moduleGenes <- moduleColors==module
    verboseScatterplot(abs(geneModuleMembership[moduleGenes, column]),
                       abs(geneTraitSignificance[moduleGenes, 1]),
                       xlab = paste("Module Membership in", module, "module"),
                       ylab = "Gene significance for trait",
                       main = paste("Module membership vs. gene significance\n"),
                       cex.main = 1.2, cex.lab = 1.2, cex.axis = 1.2, col = module)
  }
  dev.off()
}

## select brown module genes
selectModule <- "brown"
column <- match(selectModule, modNames)
moduleGenes <- moduleColors==selectModule
MM_brown <- abs(geneModuleMembership[moduleGenes, column])
GS_brown <- abs(geneTraitSignificance[moduleGenes, 1])
c <- as.data.frame(cbind(MM_brown,GS_brown))
write.csv(c, "MMGS_brown.csv")