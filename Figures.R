library(wesanderson)
library(Seurat)
library(ggplot2)
library(patchwork)
library(dplyr)
library(tidyr)
library(magrittr)
library(EnhancedVolcano)
library(topGO)
library("org.Mm.eg.db")

combined <- c(wes_palette('GrandBudapest2'), wes_palette('Royal2'),wes_palette('Darjeeling1'), wes_palette('Cavalcanti1'), wes_palette('FantasticFox1'))

#UMAP plot for initial clusters + bar plot Fig 1D
pC<-UMAPPlot(clusters, label = T, cols=combined, label.size =6)
clust <-as.data.frame(table(clusters$tree.ident))
run <- 'Integrated'
length(run) <- 14
run <- replace_na(run, 'Integrated')
clust <- cbind(run, clust)
colnames(clust) <-c('RunID', 'ClusterID', 'CellCount')
totalcells <- sum(clust$CellCount)
Percent <- ((clust$CellCount)/totalcells) * 100
Percent <- round(Percent, digits=0)
clust <- cbind(clust, Percent)


ggplot(clust, aes(fill= ClusterID, y= CellCount, x= RunID)) +
  geom_bar( stat="identity", width = 0.1) +
  scale_fill_manual(values = combined)+
  theme( panel.grid.major = element_blank(), 
         panel.grid.minor = element_blank(),
         panel.background = element_blank(),
         axis.text.y = element_blank(),
         axis.title.y = element_blank(),
         axis.line = element_line(colour = "black")) +
  labs(y = 'Cell Count') + 
  scale_x_discrete(expand = c(0,0))+
  scale_y_continuous(expand = c(0,0))+
  geom_text(aes(label = Percent), position = position_stack(vjust = 0.5)) + 
  NoLegend() + coord_flip()


#Figure 1E
DotPlot(clusters, features  =c('Cx3cr1', 'Itgam', 'Ptprc', 'Cd4', 'Cd8a','Cd19','Cd80', 'Rbfox3', 'Cspg4', 'Nes', 'Gfap', 'Mbp', 'Olig1'), cols = 'RdBu')+ RotatedAxis() + labs(y="Cluster ID", x = "Gene")

#Figure 1F Feature PLOTS


p0 <- FeaturePlot(clusters, features = 'Cx3cr1', cols = c("#d5d0bc", "#ff5e62"))+ 
  theme_void() + ggtitle('Cx3cr1' ) +
  theme(plot.title = element_text(hjust = 0.5, face = 'bold'))

p1 <- FeaturePlot(clusters, features = 'Tmem119', cols = c("#d5d0bc", "#ff5e62"))+ 
  theme_void() + ggtitle('Tmem119' ) +
  theme(plot.title = element_text(hjust = 0.5, face = 'bold'))

p2 <- FeaturePlot(clusters, features = 'Aif1',cols = c("#d5d0bc", "#ff5e62")) + 
  theme_void() + ggtitle('Iba1' ) +
  theme(plot.title = element_text(hjust = 0.5, face = 'bold'))

p3 <- FeaturePlot(clusters, features = 'Ptprc', cols = c("#d5d0bc", "#ff5e62"))+ 
  theme_void() + ggtitle('Cd45' ) +
  theme(plot.title = element_text(hjust = 0.5, face = 'bold'))

p4 <-FeaturePlot(clusters, features = 'Cd68',cols = c("#d5d0bc", "#ff5e62"))+
  theme_void() + ggtitle('Cd68' ) +
  theme(plot.title = element_text(hjust = 0.5, face = 'bold'))

p5 <- FeaturePlot(clusters, features = 'Cd63',cols = c("#d5d0bc", "#ff5e62"))+ 
  theme_void() + ggtitle('Cd63' ) +
  theme(plot.title = element_text(hjust = 0.5, face = 'bold'))

p6 <- FeaturePlot(clusters, features = 'Apoe',cols = c("#d5d0bc", "#ff5e62"))+
  theme_void() + ggtitle('Apoe' ) +
  theme(plot.title = element_text(hjust = 0.5, face = 'bold'))

p7 <- FeaturePlot(clusters, features = 'Egr1', cols = c("#d5d0bc", "#ff5e62"))+
  theme_void() + ggtitle('Egr1' ) +
  theme(plot.title = element_text(hjust = 0.5, face = 'bold'))

p8 <- FeaturePlot(clusters, features = 'Exosc2', cols =c("#d5d0bc", "#ff5e62"))+
  theme_void() + ggtitle('Exosc2' ) +
  theme(plot.title = element_text(hjust = 0.5, face = 'bold'))

p9 <- FeaturePlot(clusters, features = 'Rtp4', cols =c("#d5d0bc", "#ff5e62"))+
  theme_void() + ggtitle('Rtp4' ) +
  theme(plot.title = element_text(hjust = 0.5, face = 'bold'))

p10 <- FeaturePlot(clusters, features = 'Ccr1',cols = c("#d5d0bc", "#ff5e62"))+
  theme_void() + ggtitle('Ccr1' ) +
  theme(plot.title = element_text(hjust = 0.5, face = 'bold'))

p11 <- FeaturePlot(clusters, features = 'Ifitm3', cols = c("#d5d0bc", "#ff5e62"))+ 
  theme_void() + ggtitle('Ifitm3' ) +
  theme(plot.title = element_text(hjust = 0.5, face = 'bold'))

p12 <- FeaturePlot(clusters, features = 'Top2a', cols = c("#d5d0bc", "#ff5e62"))+ 
  theme_void() + ggtitle('Top2a' ) +
  theme(plot.title = element_text(hjust = 0.5, face = 'bold'))


p13 <- FeaturePlot(clusters, features = 'Mrc1', cols = c("#d5d0bc", "#ff5e62"))+ 
  theme_void() + ggtitle('Mrc1' ) +
  theme(plot.title = element_text(hjust = 0.5, face = 'bold'))


p14 <- FeaturePlot(clusters, features = 'Cd74',cols = c("#d5d0bc", "#ff5e62"))+ 
  theme_void() + ggtitle('Cd74' ) +
  theme(plot.title = element_text(hjust = 0.5, face = 'bold'))


(p0|p1|p2|p3|p4)/(p5|p6|p7|p8|p9)/(p10|p11|p12|p13|p14)


#Figure 2H

VlnPlot(clusters, features = 'Cd68', pt.size = 0, split.by = 'orig.ident', slot = 'counts') + 
  geom_boxplot(outlier.size = 0 ,
               position = position_dodge(width = 0.9)) +
  xlab('Cluster ID') 

#Figure 3A Volcano Plot for cluster 8

SGZ <- FindMarkers(clusters, ident.1 = 8, ident.2 = c(13,14), logfc.threshold = 0.001)
EnhancedVolcano(SGZ,
                lab = rownames(SGZ),
                x = 'avg_log2FC',
                y = 'p_val_adj',
                pCutoff = 10e-32,
                FCcutoff = 0.5,
                pointSize = 3.0,
                labSize = 6.0,
                colAlpha = 1,
                legendLabels=c('n.s.','Log2FC','adj p-val',
                               'adj p-val + log2FC'),
                legendPosition = 'bottom',
                legendLabSize = 16,
                legendIconSize = 5.0,
                gridlines.major = FALSE,
                gridlines.minor = FALSE,
                title = NULL,
                subtitle =  NULL)

#Figure 3B Violin plots of downregulated/upregulated cluster 8
VlnPlot(clusters, features = c('Tmem119', 'P2ry12','Fcrls',  'Olfml3', 'Hexb', 'Tgfbr1', 'Gpr34', 'Sall1', 'Mertk', 'P2ry13', 'Csf1r', 'Selplg','C1qa', 'C1qb','Cd34'), idents = c(8,13,14), pt.size = 0, ncol = 5, cols =combined[c(8,13,14)])
VlnPlot(clusters, features = c('Cd9', 'Cd63', 'Ftl1', 'Ctsz', 'Ctsb','Ctsd', 'Mif', 'C3ar1', 'C5ar1', 'Apoe', 'Mt1', 'Lyz2', 'Cd83','Lpl', 'Csf1' ), idents = c(8,13,14), pt.size = 0, ncol = 5,cols = combined[c(8,13,14)])


#Figure 3C GO terms

DGE_all <- FindAllMarkers(clusters, logfc.threshold = 0.25, only.pos = TRUE,assay = 'RNA')
annotation <- function(genelist, GOquant)
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
                     topNodes = GOquant,
                     numChar = 200)  
  return(allRes)
}


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

topDiffGenes <- function (genelist) 
{
  return(genelist < 0.01)
}

cluster8 <- clusterGO(DGE_all, 8, 500) %>% filter(elimKS < 0.05)
neuronup <- c8[c(126,80,142,71, 125, 21, 123,89,122, 88,43, 55, 138, 40, 132), ]

ggplot(neuronup, aes(x =Annotated,y = reorder(Term, Annotated))) +
  theme_classic()+
  geom_col(fill = 'forestgreen')  +
  theme(text = element_text(size = 20), axis.title=element_text(size=18,face="bold")) +  
  scale_y_discrete(position = 'left') +
  xlab('Number of genes') +
  ylab('GO term')


#Figure 4D Immune GO terms

immuneUP <- c8[c(34,22,139,2,19,5,13,62,10.150,100,72,64),]

ggplot(immuneUP, aes(x =Significant,y = reorder(Term, Significant))) +
  theme_classic()+
  geom_col(fill = 'firebrick')  +
  theme(text = element_text(size = 20), axis.title=element_text(size=18,face="bold")) +  
  scale_y_discrete(position = 'left') +
  xlab('Number of genes') +
  ylab('GO term')

#Figure 4E UMAP Plot
load('/Volumes/Critical\ Care/Kernie\ Lab/Sana/Sana_pal/fin/chintamen_shaul.rda')
DimPlot(chintamen_shaul.harmony.clust_intersect, pt.size = 0.6,label.size = 6,split.by = 'project_orig', label = TRUE, cols = c(g,r,c)) + NoLegend() 
#Figure 4F DotPlot for DAM genes
DotPlot(chintamen_shaul.harmony.clust_intersect, features = c(c('Cd33', 'Cx3cr1', 'P2ry12', 'P2ry13', 'Tgfbr1', 'Txnip', 'Glul', 'Tmem119'), c('Tyrobp', 'Ctsb', 'Cstb', 'Ctsd', 'Apoe', 'B2m', 'Fth1', 'Timp2', 'H2-D1', 'Lyz2'),c('Trem2', 'Ank', 'Cd63', 'Cd9','Serpine2', 'Ctsz', 'Cd68', 'Cadm1', 'Spp1','Cd52', 'Ctsa', 'Clec7a', 'Axl', 'Ctsl','Lpl','Ccl6', 'Csf1', 'Hif1a', 'Gusb','Cst7', 'Itgax')), cols = 'RdBu') + RotatedAxis() + xlab('')

#Figure 4G Venn Diagram of upregulated genes
shaul <- FetchData(chintamen_shaul.harmony.clust_intersect, vars = c('project_orig', "seurat_clusters"))

Shaul_homeostatic <- filter(shaul, project_orig == 'Shaul') %>% filter(seurat_clusters %in% homeostatic) %>% rownames()

Shaul_DAM <- filter(shaul, project_orig == 'Shaul') %>%filter(seurat_clusters %in% DAM) %>% rownames()

Chintamen_homeostatic <- filter(shaul, project_orig == 'Chintamen') %>% filter(seurat_clusters %in% homeostatic)%>% rownames()

Chintamen_DAM <- filter(shaul, project_orig == 'Chintamen') %>% filter(seurat_clusters %in% DAM) %>% rownames()



# Cluster breakdown
clust <-as.data.frame(table(clusters$tree.ident))
run <- 'Integrated'
length(run) <- 14
run <- replace_na(run, 'Integrated')
clust <- cbind(run, clust)
colnames(clust) <-c('RunID', 'ClusterID', 'CellCount')
totalcells <- sum(clust$CellCount)
Percent <- ((clust$CellCount)/totalcells) * 100
Percent <- round(Percent, digits=0)
clust <- cbind(clust, Percent)

ggplot(clust, aes(fill= ClusterID, y= CellCount, x= RunID)) +
  geom_bar( stat="identity", width = 0.1) +
  scale_fill_manual(values = combined)+
  theme( panel.grid.major = element_blank(), 
         panel.grid.minor = element_blank(),
         panel.background = element_blank(),
         axis.text.y = element_blank(),
         axis.title.y = element_blank(),
         axis.line = element_line(colour = "black")) +
  labs(y = 'Cell Count') + 
  scale_x_discrete(expand = c(0,0))+
  scale_y_continuous(expand = c(0,0))+
  geom_text(aes(label = Percent), position = position_stack(vjust = 0.5)) + 
  NoLegend()

#Myeloid markers
DotPlot(clusters, assay = 'RNA', features =c('Cx3cr1','Tmem119', 'P2ry12', 'Hexb','Sall1', 'Aif1','Csf1r','Siglech','Itgam','Ptprc','Cd68','Itgax','Ccr2','Pf4','Mrc1','Ms4a7','Cd163', 'Tlr2','Tlr4'), cols = 'RdBu')+ RotatedAxis() + labs(y="Cluster ID", x = "Gene")

#Sex specific differences
FeaturePlot(clusters, features = 'Xist', cols = c('aliceblue', 'tomato1'), label = TRUE)

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

MG_sex <- merge(MG_M, MG_F)
Idents(object=MG_sex) <- "Sex"
MG_FvM <- FindMarkers(MG_sex, ident.1 = 'M', ident.2 = 'F')
MvF<- as.data.frame(MG_FvM)
genes <- filter(MvF, p_val_adj < 0.001) %>% rownames()

levels(MG_sex) <- c('F', 'M')
VlnPlot(MG_sex, idents = c('F', 'M'), features = genes, pt.size = 0, cols = c('tomato1','cornflowerblue'))

clust_m <-as.data.frame(table(MG_M$tree.ident))
colnames(clust_m) <-c('ClusterID', 'CellCount')
totalcells <- sum(clust_m$CellCount)
Percent <- ((clust_m$CellCount)/totalcells) * 100
clust_m <- cbind(clust_m, Percent)
Sex <- 'M'
length(Sex) <- 14
Sex <- replace_na(Sex, 'M')
clust_m <- cbind(clust_m, Sex)

clust_f <-as.data.frame(table(MG_F$tree.ident))
colnames(clust_f) <-c('ClusterID', 'CellCount')
totalcells <- sum(clust_f$CellCount)
Percent <- ((clust_f$CellCount)/totalcells) * 100
clust_f <- cbind(clust_f, Percent)
Sex <- 'F'
length(Sex) <- 14
Sex <- replace_na(Sex, 'F')
clust_f <- cbind(clust_f, Sex)

clust_breakdown <- rbind (clust_m, clust_f) %>%as.data.frame()

ggplot(clust_breakdown, aes(fill= Sex, y= Percent, x= ClusterID)) +
  geom_bar( stat="identity", position = "dodge") +
  theme(text = element_text(size=18),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        panel.background = element_blank(),
        axis.line = element_line(colour = "black")) +
  scale_fill_manual(values =c('tomato1','cornflowerblue'))


SGZ_MG <- subset(clusters, idents = 8)
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

SGZ <- merge(MG_F, MG_M)

DefaultAssay(SGZ) <- 'RNA'
Idents(object=SGZ) <- "Sex"
SGZ_MvF <- FindMarkers(SGZ, ident.1 = 'M', ident.2 = 'F')
SGZ_MvF <- as.data.frame(SGZ_MvF)
genes <- filter(SGZ_MvF, p_val_adj < 0.001)
genes <- rownames(genes)
VlnPlot(SGZ, idents = c('M','F'), features = genes, pt.size = 0, cols = c( 'tomato1','cornflowerblue'))


#DGE between clusters
top_markers <- DGE %>% group_by(cluster) %>% top_n(n = 5, wt = avg_log2FC)
genelist <- unique(top_markers$gene)
DoHeatmap(subset(clusters, downsample =30), features = genelist, group.colors = combined) + scale_fill_gradientn(colors = c("blue", "white", "red"))+ theme(text = element_text(size = 20))
