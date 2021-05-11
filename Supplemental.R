library(Seurat)
library(wesanderson)

clusters <- readRDS(file = '/Users/sanachintamen/Desktop/Cx3cr1_clusters.Rds', refhook = NULL)
DefaultAssay(clusters) <- 'integrated'

######
p0 <- UMAPPlot(clusters, label = TRUE) + NoLegend() + ggtitle('Resolution = 0.5')
pdf(file = '/Users/sanachintamen/Analysis/ThesisPlots/DefaultClustersUMAP.pdf',
    width = 8)
p0
dev.off()
DefaultAssay(clusters) <- 'integrated'

clusters <- FindClusters(clusters, resolution = 0.1)
p1 <- UMAPPlot(clusters, label = TRUE) + NoLegend() + ggtitle('Resolution = 0.1')
pdf(file = '/Users/sanachintamen/Analysis/ThesisPlots/DefaultClustersUMAP_0.1.pdf',
    width = 8)
p1
dev.off()



clusters <- FindClusters(clusters, resolution = 0.3)
p2 <- UMAPPlot(clusters, label = TRUE) + NoLegend() + ggtitle('Resolution = 0.3')
pdf(file = '/Users/sanachintamen/Analysis/ThesisPlots/DefaultClustersUMAP_0.3.pdf',
    width = 8)
p2
dev.off()

clusters <- FindClusters(clusters, resolution = 0.8)
p3 <- UMAPPlot(clusters, label = TRUE) + NoLegend() + ggtitle('Resolution = 0.8')
pdf(file = '/Users/sanachintamen/Analysis/ThesisPlots/DefaultClustersUMAP_0.8.pdf',
    width = 8)
p3
dev.off()


clusters <- FindClusters(clusters, resolution = 1.0)
p4 <- UMAPPlot(clusters, label = TRUE) + NoLegend() + ggtitle('Resolution = 1.0')
pdf(file = '/Users/sanachintamen/Analysis/ThesisPlots/DefaultClustersUMAP_1.0.pdf',
    width = 8)
p4
dev.off()

clusters <- FindClusters(clusters, resolution = 1.2)
p5 <- UMAPPlot(clusters, label = TRUE) + NoLegend() + ggtitle('Resolution = 1.2')
pdf(file = '/Users/sanachintamen/Analysis/ThesisPlots/DefaultClustersUMAP_1.2.pdf',
    width = 8)
p5
dev.off()


#####
library(clustree)
clusters <- FindClusters(
  object = clusters,
  reduction.type = "pca",
  resolution = c(0.1, 0.3, 0.5, 0.8, 1.0, 1.2),
  dims.use = 1:40,
  save.SNN = TRUE
)


clustree(clusters)

pdf(file = '/Users/sanachintamen/Analysis/ThesisPlots/Clustree.pdf',
    width = 10,
    height = 8.5)

clustree(clusters)

dev.off()

#####

clusters <- readRDS(file = '/Users/sanachintamen/Desktop/Cx3cr1_clusters.Rds', refhook = NULL)
DefaultAssay(clusters) <- 'integrated'

MG_tree <- BuildClusterTree(clusters, reorder.numeric = TRUE, reorder = TRUE)

p0 <- UMAPPlot(MG_tree, label = TRUE) + NoLegend() +ggtitle('Reordered')

pdf(file = '/Users/sanachintamen/Analysis/ThesisPlots/DendroUMAP.pdf',
    width = 8)
p0
dev.off()

library(wesanderson)
d <-  wes_palette('Darjeeling1')
r <- wes_palette('Royal2')
g <- wes_palette('GrandBudapest2')
f <- wes_palette('FantasticFox1')
c <- wes_palette('Cavalcanti1')
combined <- c(g,r,d,c,f)

p1 <- UMAPPlot(MG_tree, label = TRUE, cols= combined) + NoLegend() +ggtitle('Reordered')

pdf(file = '/Users/sanachintamen/Analysis/ThesisPlots/DendroUMAP_recolored.pdf',
    width = 8)
p1
dev.off()


#####
library(dplyr)
all.markers <- FindAllMarkers(object = clusters)
top20 <- all.markers %>% group_by(cluster) %>% top_n(20, avg_log2FC)
DoHeatmap(object = clusters, features =top20$gene)


library(dplyr)
DefaultAssay(MG_tree) <- 'RNA'

clusters <- NormalizeData(MG_tree)
clusters <- FindVariableFeatures(clusters)
all.genes <- rownames(clusters)
clusters <- ScaleData(clusters, features = all.genes)

pdf(file = '/Users/sanachintamen/Analysis/ThesisPlots/Apoe_Trem2_Dap12_unscaled.pdf',
    width = 6.5,
    height = 4)

VlnPlot(clusters, features = c('Apoe', 'Trem2', 'Tyrobp'), 
        pt.size = 0, 
        cols = combined, 
        slot = 'counts', 
        stack = TRUE, 
        fill.by = 'ident', 
        flip = TRUE)  + 
  geom_boxplot(outlier.shape = NA) + 
  RotatedAxis()

dev.off()


pdf(file = '/Users/sanachintamen/Analysis/ThesisPlots/Apoe_Trem2_Dap12_scaled.pdf',
    width = 6.5,
    height = 4)

VlnPlot(clusters, features = c('Apoe', 'Trem2', 'Tyrobp'), 
        pt.size = 0, 
        cols = combined,  
        stack = TRUE, 
        fill.by = 'ident', 
        flip = TRUE)  + 
  geom_boxplot(outlier.shape = NA) + 
  RotatedAxis()

dev.off()




#####
DGE <- FindAllMarkers(clusters, logfc.threshold = 0.25, min.pct = 0.7, only.pos = TRUE, assay = 'RNA')

clusters <- RunPCA(object = clusters)
clusters <- FindNeighbors(object = clusters, dims = 1:40)
clusters <- FindClusters(object = clusters, resolution = 0.5)
clusters <- RunUMAP(object = clusters)
DimPlot(object = clusters, reduction = "umap", label = TRUE)


#####
SGZ_sub <- FindSubCluster(
  clusters,
   '8',
  graph.name = "integrated_snn",
  subcluster.name = "sub.cluster",
  resolution = 0.1,
  algorithm = 1
)

BAM_sub <- FindSubCluster(
  SGZ_sub,
  '2',
  graph.name = "integrated_snn",
  subcluster.name = "sub.cluster",
  resolution = 0.1,
  algorithm = 1
)

DimPlot(BAM_sub, group.by = 'sub.cluster')

SGZ_sub <- SetIdent(SGZ_sub, value = 'sub.cluster')
DimPlot(SGZ_sub, label = TRUE)

DimPlot(SGZ_sub, group.by = 'sub.cluster', label = TRUE)

DimPlot(SGZ_sub, group.by = 'sub.cluster', label = TRUE)

SGZ <- subset(SGZ_sub, idents = 8)

VlnPlot(SGZ, features = 'Pcna', split.by = 'sub.cluster')

DimPlot(SGZ, group.by = 'sub.cluster', cols = combined)


SGZ <- SetIdent(SGZ, value = 'sub.cluster')
SGZ_sub_DGE <- FindAllMarkers(SGZ, only.pos = TRUE, logfc.threshold = 0.5)


SGZ_sub <- SetIdent(SGZ_sub, value = 'sub.cluster')

p10 <- DimPlot(SGZ, group.by = 'sub.cluster', cols = combined, label = TRUE) + NoLegend() + ggtitle('Sub Clusters')

pdf(file = '/Users/sanachintamen/Analysis/ThesisPlots/subclusterUMAP.pdf',
    height = 8)
p10
dev.off()

write.csv(SGZ_sub_DGE, file = '/Users/sanachintamen/Analysis/ThesisTables/SGZ_DGE.csv')

FeaturePlot(SGZ, features = c('Cd68', 'Lpl'), blend = TRUE )
FeaturePlot(SGZ, features = c('Cd9', 'Cd63'), blend = TRUE )
FeaturePlot(SGZ, features = c('Cd9', 'Cd83'), blend = TRUE )
FeaturePlot(SGZ, features = c('Ccl3', 'Ccl4'), blend = TRUE )
pdf(file  = '/Users/sanachintamen/Analysis/ThesisPlots/Feature_subcluster_c3ar1_c5ar1.pdf',
    height = 7,
    width = 21)
FeaturePlot(SGZ, features = c('C3ar1', 'C5ar1'), blend = TRUE )
dev.off()

pdf(file  = '/Users/sanachintamen/Analysis/ThesisPlots/Feature_subcluster_c3ar1_c5ar1.pdf',
    height = 7,
    width = 21)
FeaturePlot(SGZ, features = c('C3ar1', 'C5ar1'), blend = TRUE )
dev.off()

pdf(file  = '/Users/sanachintamen/Analysis/ThesisPlots/Feature_subcluster_cd83_ctsz.pdf',
    height = 7,
    width = 21)
FeaturePlot(SGZ, features = c('Cd83', 'Ctsz'), blend = TRUE)
dev.off()
FeaturePlot(SGZ, features = c('Mt1', 'Csf1'), blend = TRUE)
newIdent <- c('8_0', '8_1')
SGZ <- RenameIdents(SGZ,c('8_1', '8_2'))


SGZ_sub_DGE <- FindAllMarkers(SGZ_sub, logfc.threshold = 0.25, only.pos = TRUE, min.pct = 0.7)


library(dplyr)
library(magrittr)

stats <- table(MG_tree@meta.data[["tree.ident"]])
stats <- as.data.frame(stats)
count <- stats[[2]]
percentage <- 100 * (count/sum(stats[[2]]))

SGZvsBAM <- FindMarkers(clusters, ident.1 = 8, ident.2 = 2, min.pct = 0.3, only.pos = TRUE )

library(ggplot2)




pdf(file = '/Users/sanachintamen/Analysis/ThesisPlots/Lgals3.pdf',
    height = 5,
    width = 13)
VlnPlot(clusters, cols = combined, features = 'Lgals3', pt.size = 0.1, assay = 'RNA') + ylim(c(0.1,4))

dev.off()












