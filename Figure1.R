library(Seurat)
library(wesanderson)
library(RColorBrewer)
library(ggplot2)
library(patchwork)
library(dplyr)


clusters <- readRDS(file = '/Users/sanachintamen/Desktop/Cx3cr1_clusters.Rds', refhook = NULL)
DefaultAssay(clusters) <- 'integrated'
MG_tree <- BuildClusterTree(clusters, reorder.numeric = TRUE, reorder = TRUE)

PlotClusterTree(MG_tree)

library(ape)
phylo = plot.phylo(MG_tree@tools[["BuildClusterTree"]], edge.width = 0.5, label.offset = 1,cex = 1,font = 1,main = "Unsupervised clustering")
phylo = plot.phylo(MG_tree@tools[["BuildClusterTree"]], 
                   edge.width = 0.5, 
                   label.offset = 1,
                   cex = 1,font = 1,
                   direction = 'downwards',
                   node.pos = 1,
                   main = "Unsupervised clustering") 
library(clustree)
clustree(MG_tree, assay = 'integrated')


d <-  wes_palette('Darjeeling1')
r <- wes_palette('Royal2')
g <- wes_palette('GrandBudapest2')
f <- wes_palette('FantasticFox1')
c <- wes_palette('Cavalcanti1')
combined <- c(g,r,d,c,f)

#####  UMAP PLOT
pC<-UMAPPlot(MG_tree, label = T, cols=combined)
pC <- pC + NoAxes() + NoLegend()

ggsave(pC,
        file = '/Users/sanachintamen/Desktop/UMAP.pdf',
       device = pdf,
       width = 150,
       height = 90,
       units = 'mm'
       )

##### 
#Cell type classification
DefaultAssay(MG_tree) <- 'RNA'
clusters <- NormalizeData(MG_tree)
clusters <- FindVariableFeatures(clusters)
all.genes <- rownames(clusters)
clusters <- ScaleData(clusters, features = all.genes)
DGE <- FindAllMarkers(clusters, logfc.threshold = 0.25, min.pct = 0.7, only.pos = TRUE, assay = 'RNA')

#Rbfox3 = NeuN
markergenes <- c('Cx3cr1', 'Itgam', 'Ptprc', 'Cd4', 'Cd8a','Cd19','Cd80', 'Rbfox3', 'Cspg4', 'Nes', 'Gfap', 'Mbp', 'Olig1')

p1 <- DotPlot(clusters, features  = markergenes, cols = 'RdBu')+ RotatedAxis() + labs(y="Cluster ID", x = "Gene")
p1
ggsave(plot = p1,
       filename = '/Users/sanachintamen/Desktop/markergenes.pdf',
       device = pdf,
       width = 250,
       height = 110,
       units = 'mm')

######
p3 <- VlnPlot(clusters, features = 'nFeature_RNA', pt.size = 0, cols = combined) + 
  NoLegend() + 
  ggtitle('Number of Unique Genes') + 
  xlab('ClusterID')
p4 <- VlnPlot(clusters, features = 'nCount_RNA', pt.size = 0, cols = combined) + 
  NoLegend() +
  ggtitle('Number of Transcripts') + 
  xlab('ClusterID')


p5 <- p4/p3


ggsave(plot = p5,
       filename = '/Users/sanachintamen/Desktop/QC.pdf',
       width = 180,
       height = 188,
       device = pdf,
       units = 'mm')




#####

library(dplyr)
library(tidyr)

meta <- clusters@meta.data
meta <- as.data.frame(meta)

run1 <- filter(meta, orig.ident == 'SS01')
run1 <- as.data.frame(table(run1$tree.ident))
run <- '1'
length(run) <- 14
run <- replace_na(run, '1')
run1 <- cbind(run, run1)
colnames(run1) <-c('RunID', 'ClusterID', 'CellCount')
totalcells <- sum(run1$CellCount)
Percent <- ((run1$CellCount)/totalcells) * 100
Percent <- round(Percent, digits=2)
run1 <- cbind(run1, Percent)


run2<- filter(meta, orig.ident == 'SS02')
run2 <- as.data.frame(table(run2$tree.ident))
run <- '2'
length(run) <- 14
run <- replace_na(run, '2')
run2 <- cbind(run, run2)
colnames(run2) <-c('RunID', 'ClusterID', 'CellCount')
totalcells <- sum(run2$CellCount)
Percent <- ((run2$CellCount)/totalcells) * 100
Percent <- round(Percent, digits=2)
run2 <- cbind(run2, Percent)


run3<- filter(meta, orig.ident == 'SS03')
run3 <- as.data.frame(table(run3$tree.ident))
run <- '3'
length(run) <- 14
run <- replace_na(run, '3')
run3 <- cbind(run, run3)
colnames(run3) <-c('RunID', 'ClusterID', 'CellCount')
totalcells <- sum(run3$CellCount)
Percent <- ((run3$CellCount)/totalcells) * 100
Percent <- round(Percent, digits=2)
run3 <- cbind(run3, Percent)

runs <- rbind(run1, run2, run3)

p <- ggplot(runs, aes(fill= ClusterID, y= Percent, x= RunID)) +
  geom_bar( stat="identity", width = 0.4) +
  scale_fill_manual(values = combined)+
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        panel.background = element_blank(), 
        axis.line = element_line(colour = "black")) 
p


##### 

p2 <- RidgePlot(clusters, features = 'Xist', cols = combined) + NoLegend() + labs(y="Cluster ID", x = 'Expression Level') + 
  theme(axis.text = element_text(hjust = 0.5))
p2
VlnPlot(clusters, features = 'Xist', cols = combined, pt.size = 0)



#####

FvM <- FindMarkers(clusters, ident.1 = 14, ident.2 = c(12,13), min.pct = 0.25)


FvM<- as.data.frame(FvM)
genes <- filter(FvM, p_val_adj < 0.001)
genes <- rownames(genes)
genes 

VlnPlot(clusters, idents = c(12,13,14), features = genes, pt.size = 0,  cols = c(combined[12], combined[13], combined[14]), ncol = 5) 
RidgePlot(clusters, idents = c(0,1), features = genes,cols = c(combined[12], combined[13], combined[14]), ncol = 7)



pdf(file = '/Users/sanachintamen/Analysis/ThesisPlots/XistUMAP.pdf',
    width = 9)
FeaturePlot(clusters, features = 'Xist', cols = scheme)
dev.off()
