rm(list = ls())
library(ggplot2)
library(ggrepel)
library(viridis)
library(hrbrthemes)
library(ggthemes)

sholl <- read.csv(file = '/Users/sanachintamen/Desktop/sholl.csv', header = TRUE)
colnames(sholl) <- c('Region', 'Crossings')
sholl <- as.data.frame(sholl)

p <- ggplot(sholl, aes(x = Region, y= Crossings, fill = Region)) + 
  geom_boxplot() +
  geom_jitter(shape=1, position=position_jitter(0.2)) +
  ggtitle('Sholl Analysis') +
  xlab('Region') + ylab('Total Intersections') +
  scale_colour_hc()+
  theme_bw() +
  NoLegend()
p <- p +   theme( 
  plot.title = element_text(hjust = 0.5, size =12),
  axis.title.x = element_text(size=10),
  axis.title.y = element_text(size=10))

se <- function(x) sqrt(var(x)/length(x))
se(ctx)
se(sgz)
sgz <- sholl[1:5, 2]
ctx <- sholl[6:9,2]

SE_c<-sd(ctx)/sqrt(length(ctx))

sem(ctx, na.rm = FALSE)

# 
# summary(sholl[1:5,])
# Region            Crossings 
# Length:5           Min.   : 6  
# Class :character   1st Qu.:10  
# Mode  :character   Median :13  
# Mean   :12  
# 3rd Qu.:15  
# Max.   :16  
# > summary(sholl[6:9,])
# Region            Crossings    
# Length:4           Min.   :30.00  
# Class :character   1st Qu.:32.25  
# Mode  :character   Median :35.50  
# Mean   :36.50  
# 3rd Qu.:39.75  
# Max.   :45.00  

ggplot_build(p)$data

ggsave(
  '/Users/sanachintamen/Desktop/sholl.jpeg',
  plot = p,
  device = 'jpeg',
  width = 65.2,
  height = 80,
  units = 'mm',
  dpi = 300
)



p <- ggplot() + 
  geom_col(data= sholl, aes(x = Region, y = Crossings, fill)) +
  scale_fill_viridis(discrete = TRUE, alpha = 0.6, option = 'A') +
  theme_classic() 
p  
\



library(ggplot2)
M1M2markers <- c('Il6', 'Mertk', 'Tnf', 'Cd86', 'Stat1', 'Nos2', 'Il1a','Cd115','Ifng', 'Arg1', 'Cd163', 'Tgfb1','Tgfb2','Tgfb3' ,'Stat6', 'Retnla')

DotPlot(clusters, features = M1M2markers , cols = 'RdBu') +RotatedAxis()

m2 <- c('Ccl2', 'Gas6', 'Tgfbi', 'Arg1', 'Cd163', 'Mertk', 'Csf1r', 'Mrc1', 'Stat6', 'Retnla', 'Il10', 'Irf4')
DotPlot(clusters, features = m2 , cols = 'RdBu') +RotatedAxis()

#cd16/32 is Fcgr2b
  
m1 <- c('Il1a', 'Il6', 'Ifng', 'Tnf', 'Cd68','Nos2', 'Stat1',  'Fcgr2b', 'Il12a', 'Il23a', 'Irf5')

DotPlot(clusters, features = m1 , cols = 'RdBu') +RotatedAxis()


DotPlot(clusters, features = c(m1,m2 ), cols = 'RdBu') +RotatedAxis() + xlab('Genes') + ylab('Cluster ID')
#####

library(VennDiagram)
library(ggvenn)
library(dplyr)


SGZvsMG <- FindMarkers(clusters, ident.1 = 8, ident.2 = c(13,14), logfc.threshold = 0.25)
SGZup <- filter(SGZvsMG, avg_log2FC > 0)
SGZup <- rownames(SGZup)
SGZdown <- filter(SGZvsMG, avg_log2FC < 0)
SGZdown <- rownames(SGZdown)

DAM <- read.csv('/Users/sanachintamen/Documents/kerenshaul1vs3.csv')
DAM <- as.data.frame(DAM)
DAMup <- filter(DAM, up.down == 1)
DAMup <- DAMup$X
DAMdown <- filter(DAM, up.down == -1)
DAMdown <- DAMdown$X
rm(DAM)
i <- intersect(SGZup, DAMup)


PAM <- read.csv(file = '/Users/sanachintamen/Desktop/PAM.csv')
PAM <- as.data.frame(PAM)
PAMup <- filter(PAM, avg_logFC > 0)
PAMup <- PAMup$gene
PAMdown <- filter(PAM, avg_logFC <0)
PAMdown <- PAMdown$gene
rm(PAM)



plot <- venn.diagram(list(SGZup, DAMup, PAMup),
                     filename = '/Users/sanachintamen/Desktop/upvenn2.tiff',
                     cex = 1,
                     fill = c(combined[10], combined[14], combined[12]),
                     main = 'Upregulated Genes',
                     main.cex = 2,
                     fontface = 'plain',
                     fontfamily = 'Arial',
                     height = 3100,
                     width = 4000,
                     cat.cex = 1.5,
                     margin =0.15,
                     output = TRUE,
                     category.names = c('SGZ', 'DAM ', 'PAM'))



plot <- venn.diagram(list(SGZup, DAMup),
                     filename = '/Users/sanachintamen/Analysis/ThesisPlots/VENNup.tiff',
                     cex = 2,
                     fill = c(combined[10], combined[14]),
                     main = 'Upregulated Genes',
                     main.cex = 1,
                     fontface = 'plain',
                     fontfamily = 'Arial',
                     height = 3100,
                     width = 4000,
                     cat.cex = 1.5,
                     cat.pos = 3,
                     margin =0.15,
                     cat.dist = 0.5,
                     output = TRUE,
                     category.names = c('Cluster 8', 'DAM '))




plot <- venn.diagram(list(SGZup, DAMup),
                     filename = '/Users/sanachintamen/Analysis/ThesisPlots/VENNup.tiff',
                     cex = 2,
                     fill = c(combined[10], combined[14]),
                     main.cex = 1,
                     fontface = 'plain',
                     fontfamily = 'Arial',
                     height = 3100,
                     width = 4000,
                     cat.cex = 1.5,
                     cat.pos =c( -5,3),
                     margin =0.15,
                     output = TRUE,
                     category.names = c('Cluster 8', 'DAM '))


plot <- venn.diagram(list(SGZdown, DAMdown),
                     filename = '/Users/sanachintamen/Analysis/ThesisPlots/VENNdown.tiff',
                     cex = 2,
                     fill = c(combined[10], combined[14]),
                     main.cex = 1,
                     fontface = 'plain',
                     fontfamily = 'Arial',
                     height = 3100,
                     width = 4000,
                     cat.cex = 1.5,
                     cat.pos =c( -5,3),
                     margin =0.15,
                     output = TRUE,
                     category.names = c('Cluster 8', 'DAM '))


plot <- venn.diagram(list(SGZdown, DAMdown, PAMdown),
                     filename = '/Users/sanachintamen/Desktop/downvenn2.tiff',
                     cex = 1,
                     fill =c(combined[10], combined[14], combined[12]),
                     main = 'Downregulated Genes',
                     main.cex = 2,
                     fontface = 'plain',
                     fontfamily = 'Arial',
                     height = 3100,
                     width = 4000,
                     cat.cex = 3,
                     margin =0.15,
                     category.names = c('SGZ', 'DAM ', 'PAM'))

#####
DAM_twostep <- c('Cx3cr1', 'P2ry12','Tmem119','Tyrobp','Ctsb','Ctsd','Apoe', 'B2m', 'Fth1', 'Lyz2','Trem2', 'Axl', 'Cst7', 'Ctsl', 'Lpl', 'Cd9', 'Csf1', 'Ccl6', 'Itgax', 'Lilrb4', 'Timp2')
p6 <- DotPlot(clusters, features = c('Cx3cr1', 'P2ry12','Tmem119','Tyrobp','Ctsb','Ctsd','Apoe', 'B2m', 'Fth1', 'Lyz2','Trem2', 'Axl', 'Cst7', 'Ctsl', 'Lpl', 'Cd9', 'Csf1', 'Ccl6', 'Itgax', 'Lilrb4', 'Timp2'), cols = 'RdBu')+ RotatedAxis() + labs(y="Cluster ID", x = "Gene") 

pdf(file = '/Users/sanachintamen/Analysis/ThesisPlots/DAM_subcluster.pdf',
    width = 9,
    height = 5.5)
DotPlot(SGZ_sub, features = DAM_twostep, cols = 'RdBu')+ RotatedAxis() + labs(y="Cluster ID", x = "Gene") 
dev.off()

ggsave(plot = p6,
       filename = '/Users/sanachintamen/Desktop/DAMdotplot.pdf',
       device = pdf,
       width = 200,
       height = 130,
       units = 'mm')

DAMinSGZ <- intersect(SGZup, DAMup)

library(viridis)


pdf(file = '/Users/sanachintamen/Analysis/ThesisPlots/DAMheatmap.pdf',
    width = 16,
    height = 14)
DoHeatmap(SGZ, features = DAMinSGZ) + scale_fill_viridis(option="plasma") + NoLegend()
dev.off()





DoHeatmap(SGZ, features = DAMinSGZ) + scale_fill_gradientn(colors = c("blue", "white", "red"))

DoHeatmap(SGZ, features = DAMinSGZ) + scale_fill_viridis_c(option="cividis")
DoHeatmap(SGZ, features = DAMup) + NoLegend() 


scale_fill_viridis(option="magma")




library(Seurat)
library(topGO)
library("org.Mm.eg.db")
library(genefilter)
library(dplyr)
library(tidyr)
library(tidyselect)
library(ggplot2)
library(enrichplot)

clusters <- readRDS(file = '/Users/sanachintamen/Desktop/Cx3cr1_clusters.Rds', refhook = NULL)
DefaultAssay(clusters) <- 'RNA'
clusters <- NormalizeData(clusters)
clusters <- FindVariableFeatures(clusters)
all.genes <- rownames(clusters)
clusters <- ScaleData(clusters, features = all.genes)

SGZvsMG <- FindMarkers(clusters, ident.1 = 8, ident.2 = c(13,14), logfc.threshold = 0.25)
gene <- rownames(SGZvsMG)
SGZvsMG <- cbind(gene, SGZvsMG)
SGZvsMG <- as_data_frame(SGZvsMG)
#write.csv(SGZvsMG, file = '/Users/sanachintamen/Desktop/SGZvsMG.csv')

topDiffGenes <- function (genelist) 
{
  return(genelist < 0.01)
}


UpGenes <- function (genelist) 
{
  genes_up <- filter(genelist, avg_log2FC > 0)
  gene <- genes_up$gene
  genes_up <- as.numeric(genes_up$p_val_adj)
  names(genes_up) <- gene
  rm(gene)
  return(genes_up)
}


DownGenes <- function (genelist) 
{
  genes_down <- filter(genelist, avg_log2FC < 0)
  gene <- genes_down$gene
  genes_down <- as.numeric(genes_down$p_val_adj)
  names(genes_down) <- gene
  rm(gene)
  return(genes_down)
}



up <- UpGenes(SGZvsMG)
down <- DownGenes(SGZvsMG)


annotation <- function(genelist)
{
  GOdata <- new("topGOdata",
                description = 'test run', ontology = "BP",
                allGenes = genelist, geneSel = topDiffGenes,
                annot = annFUN.org,
                nodeSize = 10,
                mapping = 'org.Mm.eg.db',
                ID = 'symbol')
  resultFisher <- runTest(GOdata, 
                          algorithm = "classic", 
                          statistic = "fisher")
  resultKS <- runTest(GOdata, 
                      algorithm = "classic", 
                      statistic = "ks")
  resultKS.elim <- runTest(GOdata, 
                           algorithm = "elim", 
                           statistic = "ks")
  allRes <- GenTable(GOdata, classicFisher = resultFisher,
                     classicKS = resultKS, 
                     elimKS = resultKS.elim,
                     orderBy = "elimKS", 
                     ranksOf = "classicFisher", 
                     topNodes = 200,
                     numChar = 200)  
  return(allRes)
}

downGO <- annotation(down)

upGO <- annotation(up)





neuronup <- upGO[c(11,23,43,56,57,86,95,92,107,122,129,134,151,177,178),]
neurondown <- downGO[c(4,5,6,7,8,10,13,14,27,39,42,53,59,74,98),]

neuronup <- as_data_frame(neuronup)
neuronup <- neuronup %>% arrange(desc(Annotated))



immuneUP <- upGO[c(25, 137,1, 2,14,22, 24,114,181, 30, 51, 49,90),]


pdf(file = '/Users/sanachintamen/Analysis/ThesisPlots/immuneUPGO.pdf',
    width = 12)

ggplot(immuneUP, aes(x =Significant,y = reorder(Term, Significant))) +
  theme_classic()+
  geom_col(fill = 'firebrick')  +
  theme(text = element_text(size = 16), axis.title=element_text(size=18,face="bold")) +  
  scale_y_discrete(position = 'left') +
  xlab('Number of genes') +
  ylab('GO term')

dev.off()




