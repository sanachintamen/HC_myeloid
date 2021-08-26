
library(Seurat)
library(harmony)
library(Matrix)
library(ggplot2)
library(dittoSeq)
library(SeuratWrappers)
library(pheatmap)
library(dplyr)

## Load Hippocampal data comprising only SGZ and homeostatic clusters (8,13 and 14)
load("chintamen_MG_raw_seurat_81314.rda",verbose=T) 
chintamen_MG_raw_seurat_81314 #31040 features across 11433 cells; mt removed
DefaultAssay(chintamen_MG_raw_seurat_81314) <- "RNA";

## Load Karen-Shaul et. al.'s data comprising only microglial cells

load("shaul_MG_raw_seurat_feb.fin.rda")
shaul_MG_raw_seurat_ #33979 features across 10337 cells; mt removed
DefaultAssay(shaul_MG_raw_seurat_) <- "RNA"

## Common genes in both datasets
chintamen_shaul.genes<- intersect(rownames(shaul_MG_raw_seurat_),rownames(chintamen_MG_raw_seurat_81314))

## Integration with Hamrony
merge_chintamen_shaul_intersect<-merge(shaul_MG_raw_seurat_[chintamen_shaul.genes,],
                                       chintamen_MG_raw_seurat_81314[chintamen_shaul.genes,], 
                                       add.cell.ids = c("shaul","Chintamen"),
                                       project = "Chintamen_shaul")
merge_chintamen_shaul.norm_intersect <- NormalizeData(merge_chintamen_shaul_intersect) %>% 
                                        FindVariableFeatures() %>%
                                        ScaleData() %>% RunPCA(verbose = FALSE)
#ElbowPlot(merge_chintamen_shaul.norm_intersect)
chintamen_shaul.runharmony_intersect <- RunHarmony(merge_chintamen_shaul.norm_intersect,
                                                   group.by.vars = "project_orig",
                                                   epsilon.cluster = -Inf, epsilon.harmony = -Inf,
                                                   max.iter.harmony=10,theta = 2)
chintamen_shaul.harmony.rumumap_intersect <- RunUMAP(chintamen_shaul.runharmony_intersect, 
                                                     reduction = "harmony", dims = 1:10,seed.use = 45)

chintamen_shaul.harmony.neigh_intersect <- FindNeighbors(chintamen_shaul.harmony.rumumap_intersect,
                                                         reduction = "harmony", dims = 1:10) %>%
                                                          FindClusters(resolution = c(0.5, 0.7, 0.9))

#save(chintamen_shaul.harmony.neigh_intersect,file="chintamen_shaul.harmony.neigh_intersect.rda")

## UMAP
DimPlot(chintamen_shaul.harmony.neigh_intersect, 
          group.by = "seurat_clusters",label = T)+CTCN_chroma("mixed3")# customized color palette
D1<-DimPlot(chintamen_shaul.harmony.neigh_intersect, group.by = "project_orig",cols = c("#B00B1E","#ffbb91"),label = T)#original identity
D2<-DimPlot(chintamen_shaul.harmony.neigh_intersect, group.by = "seurat_clusters",label = T)+CTCN_chroma("mixed3")#integrated clusters
D3<-DimPlot(chintamen_shaul.harmony.neigh_intersect, group.by = "tree.ident_sana",label = T)+CTCN_chroma("mixed4") #original dataset clusters
D4<-DimPlot(chintamen_shaul.harmony.neigh_intersect, group.by = "seurat_clusters_shaul",label = T)+CTCN_chroma("mixed5")#original dataset clusters
D1|D2
D3|D4

## Dotplot of up or down-regulated stage-1 TREM2-independent and Step-2 Trem2-dependent marker genes
DotPlot(chintamen_shaul.harmony.clust_intersect, features = c(DAM_down, DAMstep1, DAMstep2),
        cols = "RdBu",scale.max = 100)+RotatedAxis()






