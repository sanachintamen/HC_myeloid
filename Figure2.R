library(patchwork)
library(Seurat)
library(wesanderson)
scheme <- c("#d5d0bc", "#ff5e62")

clusters <- readRDS(file = '/Users/sanachintamen/Desktop/Cx3cr1_clusters.Rds', refhook = NULL)
DefaultAssay(clusters) <- 'integrated'
MG_tree <- BuildClusterTree(clusters, reorder.numeric = TRUE, reorder = TRUE)


d <-  wes_palette('Darjeeling1')
r <- wes_palette('Royal2')
g <- wes_palette('GrandBudapest2')
f <- wes_palette('FantasticFox1')
c <- wes_palette('Cavalcanti1')
combined <- c(g,r,d,c,f)

######

DefaultAssay(MG_tree) <- 'RNA'

clusters <- NormalizeData(MG_tree)
clusters <- FindVariableFeatures(clusters)
all.genes <- rownames(clusters)
clusters <- ScaleData(clusters, features = all.genes)
DGE <- FindAllMarkers(clusters, logfc.threshold = 0.25, min.pct = 0.7, only.pos = TRUE, assay = 'RNA')


write.csv(DGE, file = '/Volumes/Critical\ Care/Kernie\ Lab/Sana/DGE.csv')

top_markers <- DGE %>% group_by(cluster) %>% top_n(n = 5, wt = avg_log2FC)
genelist <- unique(top_markers$gene)

# check http://oshlacklab.com/combes-organoid-paper/07_Combined_Clustering.html for better heat map


cols <- viridis(100)[c(1, 50, 100)]


DoHeatmap(subset(clusters, downsample =30), features = genelist, group.colors = combined) +  scale_fill_viridis()

FeaturePlot(clusters, features = 'Rtp4')
#####

scheme <- c("#FFFF33", "#ff5e62")
scheme <- c("#CCFFE5", "#ff5e62")
CCFFE5

FeaturePlot(clusters, features = 'Ccr1', cols = scheme)+ 
  theme_void() + ggtitle('Ccr1' ) +
  theme(plot.title = element_text(hjust = 0.5, face = 'bold'))

FeaturePlot(clusters, features = 'Gas6', cols = scheme)+ 
  theme(plot.title = element_text(hjust = 0.5, face = 'bold'))
FeaturePlot(clusters, features = 'Ccr1', cols = scheme)+ 
  theme(plot.title = element_text(hjust = 0.5, face = 'bold'))
FeaturePlot(clusters, features = 'Lpl', cols = scheme)+ 
  theme(plot.title = element_text(hjust = 0.5, face = 'bold'))
FeaturePlot(clusters, features = 'Ifitm3', cols = scheme)+ 
  theme(plot.title = element_text(hjust = 0.5, face = 'bold'))

FeaturePlot(clusters, features = 'Mt1', cols = scheme)+ 
  theme(plot.title = element_text(hjust = 0.5, face = 'bold'))
FeaturePlot(clusters, features = 'Gpr34', cols = scheme)+ 
  theme_void() + ggtitle('Cd9' ) +
  theme(plot.title = element_text(hjust = 0.5, face = 'bold'))

FeaturePlot(clusters, features = 'Mrc1', cols = scheme)+ 
  theme_void() + ggtitle('Ccl4' ) +
  theme(plot.title = element_text(hjust = 0.5, face = 'bold'))




p0 <- FeaturePlot(clusters, features = 'Cx3cr1', cols = scheme)+ 
  theme_void() + ggtitle('Cx3cr1' ) +
  theme(plot.title = element_text(hjust = 0.5, face = 'bold'))


p1 <- FeaturePlot(clusters, features = 'Tmem119', cols = scheme)+ 
  theme_void() + ggtitle('Tmem119' ) +
  theme(plot.title = element_text(hjust = 0.5, face = 'bold'))

p2 <- FeaturePlot(clusters, features = 'Aif1',cols = scheme) + 
  theme_void() + ggtitle('Iba1' ) +
  theme(plot.title = element_text(hjust = 0.5, face = 'bold'))

p3 <- FeaturePlot(clusters, features = 'Ptprc', cols = scheme)+ 
  theme_void() + ggtitle('Cd45' ) +
  theme(plot.title = element_text(hjust = 0.5, face = 'bold'))

p4 <-FeaturePlot(clusters, features = 'Cd68',cols = scheme)+
  theme_void() + ggtitle('Cd68' ) +
  theme(plot.title = element_text(hjust = 0.5, face = 'bold'))

p5 <- FeaturePlot(clusters, features = 'Cd63',cols = scheme)+ 
  theme_void() + ggtitle('Cd63' ) +
  theme(plot.title = element_text(hjust = 0.5, face = 'bold'))

p6 <- FeaturePlot(clusters, features = 'Apoe',cols = scheme)+
  theme_void() + ggtitle('Apoe' ) +
  theme(plot.title = element_text(hjust = 0.5, face = 'bold'))

p7 <- FeaturePlot(clusters, features = 'Egr1', cols = scheme)+
  theme_void() + ggtitle('Egr1' ) +
  theme(plot.title = element_text(hjust = 0.5, face = 'bold'))

p8 <- FeaturePlot(clusters, features = 'Exosc2', cols =scheme)+
  theme_void() + ggtitle('Exosc2' ) +
  theme(plot.title = element_text(hjust = 0.5, face = 'bold'))

p9 <- FeaturePlot(clusters, features = 'Rtp4', cols =scheme)+
  theme_void() + ggtitle('Rtp4' ) +
  theme(plot.title = element_text(hjust = 0.5, face = 'bold'))

p10 <- FeaturePlot(clusters, features = 'Ccr1',cols = scheme)+
  theme_void() + ggtitle('Ccr1' ) +
  theme(plot.title = element_text(hjust = 0.5, face = 'bold'))

p11 <- FeaturePlot(clusters, features = 'Ifitm3', cols = scheme)+ 
  theme_void() + ggtitle('Ifitm3' ) +
  theme(plot.title = element_text(hjust = 0.5, face = 'bold'))

p12 <- FeaturePlot(clusters, features = 'Top2a', cols = scheme)+ 
  theme_void() + ggtitle('Top2a' ) +
  theme(plot.title = element_text(hjust = 0.5, face = 'bold'))


p13 <- FeaturePlot(clusters, features = 'Mrc1', cols = scheme)+ 
  theme_void() + ggtitle('Mrc1' ) +
  theme(plot.title = element_text(hjust = 0.5, face = 'bold'))


p14 <- FeaturePlot(clusters, features = 'Cd74',cols = scheme)+ 
  theme_void() + ggtitle('Cd74' ) +
  theme(plot.title = element_text(hjust = 0.5, face = 'bold'))

p15 <- FeaturePlot(clusters, features = 'Ttr',cols = scheme)+ 
  theme_void() + ggtitle('Ttr' ) +
  theme(plot.title = element_text(hjust = 0.5, face = 'bold'))

pdf(file = '/Users/sanachintamen/Analysis/ThesisPlots/FeaturePlots.pdf',
    height = 8,
    width = 12)
(p0/p1/p2/p3) | (p4/p5/p6/p7) | (p8/p9/p10/p11) | (p12/p13/p14/p15)

dev.off()


DGE <- FindAllMarkers(clusters, logfc.threshold = 0.25, min.pct = 0.7, only.pos = TRUE, assay = 'RNA')

top_markers <- DGE %>% group_by(cluster) %>% top_n(n = 4, wt = avg_log2FC)
genelist <- unique(top_markers$gene)


DoHeatmap(subset(clusters, downsample =30), features = genelist, group.colors = combined) + NoLegend()

FeaturePlot(clusters, features = 'Rtp4')

#####
library(Seurat)
library(topGO)
library("org.Mm.eg.db")
library(genefilter)
library(dplyr)
library(tidyr)
library(tidyselect)
library(ggplot2)
library(enrichplot)

clusterGO <- function(DGE, clusterID, GOquant)
{
  dge <- as.data.frame(DGE)
  dge <- filter(dge, cluster == clusterID)
  dge <- dplyr::select(dge, gene, p_val_adj)
  genes <- dge$gene
  dge <- as.numeric(dge$p_val_adj)
  names(dge) <- genes
  dge <- annotation(dge, GOquant)
  return(dge)
}


cluster1 <- clusterGO(DGE, 1,30)

cluster2 <- clusterGO(DGE_all, 2, 30)

cluster3 <- clusterGO(DGE_all, 3, 30)

cluster4 <- clusterGO(DGE, 4, 30)

#none GO found
cluster5 <- clusterGO(DGE_all, 5, 30)

cluster6 <- clusterGO(DGE, 6, 30)

cluster7 <- clusterGO(DGE, 7, 30)

#no nodes
cluster8 <- clusterGO(DGE_all, 8, 200)

cluster9 <- clusterGO(DGE, 9, 30)

cluster10 <- clusterGO(DGE, 10, 50)

cluster11 <- clusterGO(DGE, 11, 30)

cluster12 <- clusterGO(DGE_all, 12, 30)

cluster13 <- clusterGO(DGE_all, 13, 30)

cluster14 <- clusterGO(DGE_all, 14, 30)

cluster2 <- as.data.frame(cluster2)
up2 <- cluster2[c(4,7,24,29,22,28,25), ]
ggplot(up2, aes(Significant, Term)) +
  geom_col(fill = combined[2]) +
  scale_y_discrete(position = 'left') +
  xlab('Number of significant genes')

cluster9 <- as.data.frame(cluster9)
up9 <- cluster9[c(2,3,19,1,9,7,29), ]
ggplot(up9, aes(Significant, Term)) +
  geom_col(fill = combined[9]) +
  scale_y_discrete(position = 'left') +
  xlab('Number of significant genes')

cluster1 <- as.data.frame(cluster1)
up1 <- cluster1[c(1,25,20,18,28,17,8), ]
ggplot(up1, aes(Significant, Term)) +
  geom_col(fill = combined[1]) +
  scale_y_discrete(position = 'left') +
  xlab('Number of significant genes')

cluster3 <- as.data.frame(cluster3)
up3 <- cluster3[c(2,21,3,8,15,4,19), ]
ggplot(up3, aes(Significant, Term)) +
  geom_col(fill = combined[3]) +
  scale_y_discrete(position = 'left') +
  xlab('Number of significant genes')

cluster10 <- as.data.frame(cluster10)
up10 <- cluster10[c(2,3,4,5,29,25,13), ]
ggplot(up10, aes(Significant, Term)) +
  geom_col(fill = combined[10]) +
  scale_y_discrete(position = 'left') +
  xlab('Number of significant genes')


#####
FeaturePlot(clusters, idents = 10, features = 'Apoe')
