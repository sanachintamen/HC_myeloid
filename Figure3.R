library(Seurat)
library(wesanderson)
library(RColorBrewer)
library(ggplot2)
library(patchwork)
library(dplyr)


clusters <- readRDS(file = '/Users/sanachintamen/Desktop/Cx3cr1_clusters.Rds', refhook = NULL)
DefaultAssay(clusters) <- 'integrated'
MG_tree <- BuildClusterTree(clusters, reorder.numeric = TRUE, reorder = TRUE)
DefaultAssay(MG_tree) <- 'RNA'
VlnPlot(MG_tree, features = 'Cd68', pt.size = 0, split.by = 'orig.ident', slot = 'counts') + 
  geom_boxplot(outlier.size = 0 ,position = position_dodge(width = 0.9)) +
  xlab('ClusterID')

pdf('/Users/sanachintamen/Desktop/Fig3a.pdf', width = 16, height = 3.5)
VlnPlot(MG_tree, features = 'Cd68', pt.size = 0, split.by = 'orig.ident', slot = 'counts') + 
  geom_boxplot(outlier.size = 0 ,
  position = position_dodge(width = 0.9)) +
  xlab('Cluster ID') 
dev.off()

VlnPlot(MG_tree, features = 'Cd68', pt.size = 0, split.by = 'orig.ident', slot = 'counts')

#with range cut
VlnPlot(MG_tree, features = 'Cd68', pt.size = 0, split.by = 'orig.ident', slot = 'counts', y.max = 15) + geom_boxplot()


