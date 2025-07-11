
BiocManager::install("clusterProfiler")
BiocManager::install("AnnotationHub")
library("dplyr")
library(ggplot2)
library(ggrepel)
library("AnnotationHub")
library("clusterProfiler") 

#使用clusterProfiler进行KEGG富集分析
#载入背景基因集数据,第一列ko_num,第二列基因ID,顺序一定不能调换
term2gene <- read.csv("././term2gene_all.csv",header=T,sep = ",")

#载入背景基因集数据,第一列基因ko_num,第二列pathway_description
term2name <- read.csv("././term2name_all.csv",header = T,sep = ",")
#指定目的基因向量
gene <- read.csv("CK_vs_Y0_75_2d.deseq2_down_geneid.csv",header = F,sep = ",")
gene <- as.factor(gene$V1)
#富集分析
x <- enricher(gene,TERM2GENE = term2gene,TERM2NAME = term2name,pvalueCutoff = 0.05,
              pAdjustMethod = "none",qvalueCutoff = 1)
head(x, 10)

#输出结果
write.csv(x, "KEGG_CK_vs_Y0_75_2d.deseq2_down.csv", quote = F, row.names = F)
#绘制条形图
barplot(x, color="pvalue")
#绘制气泡图
pvalue_ <- pvalue<0.05
dotplot(x, x="Count", color="pvalue")
cnetplot(x,
         node_label = 'all',
         colorEdge = TRUE
         )
#将富集到的基因标在通路上
browseKEGG(x, "ko03010")
#富集到的pathways的基因之间重叠关系
emapplot(x)


#有对应注释的KEGG富集分析
gene <- read.csv("892_rapid_gene_KEGGanno.csv",header=T,sep = ",")
gene <- gene[,2]
egg <- enrichKEGG(gene = gene,keyType = "kegg",organism = "nfl",pvalueCutoff = 0.05,
                  pAdjustMethod = "BH",qvalueCutoff = 0.2)

x <- GSEA(geneList = gene,TERM2GENE = term2gene,TERM2NAME = term2name)

##########################绘制气泡图#####################

KEGG_dataset <- read.table(file ="WGCNA for KEGG enrichment/KEGG_yellow.csv",
                           header = TRUE, sep = ",")

#按照Pvalue从低到高排序[升序]
KEGG_dataset <- arrange(KEGG_dataset,KEGG_dataset[,5])
#Pathway列最好转化成因子型，否则作图时ggplot2会将所有Pathway按字母顺序重排序
#将Pathway列转化为因子型
KEGG_dataset$Description <- factor(KEGG_dataset$Description,levels = rev(KEGG_dataset$Description))

#绘制KEGG气泡图
p <- ggplot(KEGG_dataset,aes(x=GeneRatio,y=Description,colour=pvalue,size=Count))+
  geom_point()+
  scale_size(range=c(2, 8))+
  scale_colour_gradient(low = "blue",high = "red")+
  theme_bw()+
  ylab("KEGG Pathway Terms")+
  xlab("GeneRatio")+
  labs(color = expression(pvalue))+
  theme(legend.title=element_text(size=14), legend.text = element_text(size = 14))+
  theme(axis.title.y = element_text(margin = margin(r = 50)),axis.title.x = element_text(margin = margin(t = 20)))+
  theme(axis.text.x = element_text(face ="bold",color="black",angle=0,vjust=1))
p


