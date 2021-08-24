library(Seurat)
library(sctransform)
library(magrittr)
library(dplyr)

# Load count matrices and filter
matrix_dir_1 = '/Users/sanachintamen/Analysis/Data/Run1/'
MG_1 <- Read10X(data.dir = matrix_dir) %>% CreateSeuratObject( , pproject = 'SS01')
MG_1[["percent.mt"]] <- PercentageFeatureSet(MG_1, pattern = "^mt-")
MG_1 <- subset(MG_1, subset = nFeature_RNA > 1000 & nCount_RNA >2500 &percent.mt < 20) %>% SCTransform( ,vars.to.regress = "percent.mt", verbose = FALSE)
DefaultAssay(MG_1) <- 'SCT'

matrix_dir_2 = '/Users/sanachintamen/Analysis/Data/Run2/'
MG_2 <- Read10X(data.dir = matrix_dir_2) %>% CreateSeuratObject(, project = 'SS02')
MG_2[["percent.mt"]] <- PercentageFeatureSet(MG_2, pattern = "^mt-")
MG_2 <- subset(MG_2, subset = nFeature_RNA > 1000 & nCount_RNA >2500 &percent.mt < 20) %>% SCTransform( , vars.to.regress = "percent.mt", verbose = FALSE)
DefaultAssay(MG_2) <- 'SCT'

matrix_dir_3 = '/Users/sanachintamen/Analysis/Data/Run3/'
MG_3 <- Read10X(data.dir = matrix_dir_3) %>% CreateSeuratObject( , project = 'SS03')
MG_3[["percent.mt"]] <- PercentageFeatureSet(MG_3, pattern = "^mt-")
MG_3 <- subset(MG_3, subset = nFeature_RNA > 1000 & nCount_RNA >2500 &percent.mt < 20) %>% SCTransform(, vars.to.regress = "percent.mt", verbose = FALSE)
DefaultAssay(MG_3) <- 'SCT'

mg_list <- c(MG_1, MG_2, MG_3)
# Integrate runs
features <- SelectIntegrationFeatures(object.list = mg_list, 
                                      nfeatures = 3000)
samples.combined <- PrepSCTIntegration(object.list = mg_list , anchor.features = features) %>%
  FindIntegrationAnchors( ,normalization.method = 'SCT', anchor.features = features)%>%
  IntegrateData( , normalization.method = 'SCT')
#Cluster 
clusters <- RunPCA(samples.combined, verbose = FALSE) %>% 
  RunUMAP( , reduction = 'pca', dims = 1:40) %>% 
  FindNeighbors( ,dims = 1:40) %>% FindClusters(, resolution = 0.5) %>%  
  BuildClusterTree( , reorder.numeric = TRUE, reorder = TRUE)
#Scale RNA
DefaultAssay(clusters) <- 'RNA'
all.genes <- rownames(clusters)
clusters <- NormalizeData(clusters) %>% FindVariableFeatures() %>% ScaleData() 
# Differential gene expression 
DGE <- FindAllMarkers(clusters, logfc.threshold = 0.25, min.pct = 0.7, assay = 'RNA')
DGE <- filter(DGE, p_val_adj<0.05)
SGZvsMG <- FindMarkers(clusters, ident.1 = 8, ident.2 = c(12,13,14), logfc.threshold = 0.25)
SGZvsMG_filt <- filter(SGZvsMG, p_val_adj<0.05)
write.csv(SGZvsMG_filt, file = '/Users/sanachintamen/Analysis/ThesisTables/Table2.csv')
