library(Seurat)
library(ggplot2)
library(RColorBrewer)
library(patchwork)
library(gdata)
saveRDS(clusters, file = '/Volumes/Critical\ Care/Kernie\ Lab/Sana/Analysis/R_objects/2021\ files/clusters.Rds')


MvF <- FindMarkers(clusters, ident.1 = 0, ident.2 = 1, only.pos = TRUE)
FvM <- FindMarkers(clusters, ident.1 = 1, ident.2 = 0, only.pos = TRUE)
MvF <- FindMarkers(clusters, ident.1 = 0, ident.2 = 1)

split

genes <- rownames(MandF)
DotPlot(clusters, features = genes)

MvF<- as.data.frame(MvF)
genes <- filter(MvF, p_val_adj < 0.001)
genes <- rownames(genes)

p0 <- VlnPlot(clusters, idents = c(1,0), features = genes, pt.size = 0, cols = c('dodgerblue2', 'tomato2'))
RidgePlot(clusters, idents = c(0,1), features = genes,cols = c('dodgerblue2', 'tomato2'))


FeaturePlot(clusters, features = 'Xist')


SGZ_MG <- subset(clusters, idents = 4)
female <- FetchData(SGZ_MG, vars = 'Xist')

MG_F <- SGZ_MG[, which(x = female > 0)]
Fcells <- colnames(MG_F@assays[["RNA"]]@counts)
tag_f <- 'F'
length(tag_f) <- length(Fcells)
tag_f[is.na(tag_f)] = 'F'
names(tag_f) <- Fcells
MG_F <- AddMetaData(
  object = MG_F,
  metadata = tag_f,
  col.name= "Sex")


MG_M <- SGZ_MG[, which(x = female == 0)]
Mcells <- colnames(MG_M@assays[["RNA"]]@counts)
tag_m <- 'M'
length(tag_m) <- length(Mcells)
tag_m[is.na(tag_m)]= 'M'
names(tag_m) <- Mcells
MG_M <- AddMetaData(
  object = MG_M,
  metadata = tag_m,
  col.name= "Sex")

cells <- rbind(tag_f, tag_m)

MG_sex <- AddMetaData(SGZ_MG, cells, col.name = 'Sex')

SGZ <- merge(MG_F, MG_M)

DefaultAssay(SGZ) <- 'integrated'
Idents(object=SGZ) <- "Sex"
DimPlot(SGZ)
FeaturePlot(clusters, features = 'Xist', cols = c('aliceblue', 'tomato1'), label = TRUE)

library(dplyr)
DefaultAssay(SGZ) <- 'RNA'
SGZ_MvF <- FindMarkers(SGZ, ident.1 = 'M', ident.2 = 'F')
SGZ_MvF <- as.data.frame(SGZ_MvF)
genes <- filter(SGZ_MvF, p_val_adj < 0.001)
genes <- rownames(genes)
p1 <- VlnPlot(SGZ, idents = c('M','F'), features = genes, pt.size = 0, cols = c( 'tomato2','dodgerblue2'))
p1


SGZ_FvM <- FindMarkers(SGZ, ident.1 = 'M', ident.2 = 'F', only.pos = TRUE)




female <- FetchData(clusters, vars = 'Xist')
MG_F <- clusters[, which(x = female > 0)]
Fcells <- colnames(MG_F@assays[["RNA"]]@counts)
tag_f <- 'F'
length(tag_f) <- length(Fcells)
tag_f[is.na(tag_f)] = 'F'
names(tag_f) <- Fcells
MG_F <- AddMetaData(
  object = MG_F,
  metadata = tag_f,
  col.name= "Sex")


MG_M <- clusters[, which(x = female == 0)]
Mcells <- colnames(MG_M@assays[["RNA"]]@counts)
tag_m <- 'M'
length(tag_m) <- length(Mcells)
tag_m[is.na(tag_m)]= 'M'
names(tag_m) <- Mcells
MG_M <- AddMetaData(
  object = MG_M,
  metadata = tag_m,
  col.name= "Sex")

library(dplyr)
cells <- rbind(tag_f, tag_m)

MG_sex <- AddMetaData(clusters, cells, col.name = 'Sex')



Idents(object=MG_sex) <- "Sex"
MG_FvM <- FindMarkers(MG_sex, ident.1 = 'M', ident.2 = 'F', only.pos = TRUE)




MvF<- as.data.frame(MG_FvM)
genes <- filter(MvF, p_val_adj < 0.001)
genes <- rownames(genes)

p1 <- VlnPlot(MG_sex, idents = c('M', 'F'), features = genes, pt.size = 0, cols = c('tomato2','dodgerblue2'))
