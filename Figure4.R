library(dplyr)


cluster8 <- clusterGO(DGE_all, 8, 500)
neuronup <- cluster8[c(11,23,43,56,57,86,95,92,107,122,129,134,143,151,177,178, 123,69),]
neuronup <- cluster8[c(69,23,43,56,57,86,95,92,107,122,129,134,143,151,177,178, 123,11),]
ggplot(neuronup, aes(Significant, Term)) +
  geom_col(fill = combined[8]) +
  scale_y_discrete(position = 'left') +
  xlab('Number of genes')


ggplot(neuronup, aes(x = Annotated, y = reorder(Term, Annotated))) +
  geom_col(fill = combined[8]) +
  scale_y_discrete(position = 'left') +
  xlab('Number of genes') +
  ylab('GO term')

VlnPlot(MG_tree, assay = 'RNA', slot = 'counts', sort = 'decreasing',features = c(MG_markers, MG_up), ncol = 5, cols = pal, idents = c(8,13,14), pt.size = 0)

pdf(file = '/Users/sanachintamen/Analysis/ThesisPlots/VlnPlots_RNAunscaled.pdf',
    width = 9,
    height = 14)
VlnPlot(MG_tree, assay = 'RNA', slot = 'counts', features = c(MG_markers, MG_up), ncol = 5, cols = pal, idents = c(8,13,14), pt.size = 0) 
dev.off()

pal <- combined[c(8,13,14)]

MG_markers <- c('Tmem119', 'P2ry12','Fcrls',  'Olfml3', 'Hexb', 'Tgfbr1', 'Gpr34', 'Sall1', 'Mertk', 'P2ry13', 'Csf1r', 'Selplg','C1qa', 'C1qb','Cd34')

MG_up <- c('Cd9', 'Cd63', 'Ftl1', 'Ctsz', 'Ctsb','Ctsd', 'Mif', 'C3ar1', 'C5ar1', 'Apoe', 'Mt1', 'Lyz2', 'Cd83','Lpl', 'Csf1' )
p0 <- VlnPlot(clusters, features = MG_markers, idents = c(13,14,8), pt.size = 0, ncol = 5, cols = pal)
p1 <- VlnPlot(clusters, features = MG_up, idents = c(13,14,8), pt.size = 0, ncol = 5,cols = pal)

p0/p1

pdf(file = '/Users/sanachintamen/Analysis/ThesisPlots/VlnPlots.pdf',
    width = 9,
    height = 14)
p0/p1

dev.off()

VlnPlot(clusters, features = 'Cx3cr1', slot = 'counts')+ geom_boxplot()




#####


library(ggrepel)
library(EnhancedVolcano)


SGZ <- FindMarkers(clusters, ident.1 = 8, logfc.threshold = 0.25)
SGZvsMG <- FindMarkers(clusters, ident.1 = 8, ident.2 = c(13,14), logfc.threshold = 0.25)
#write.csv(SGZ, file = '/Users/sanachintamen/Analysis/ThesisTables/SGZ_Fig4volcano_DGE.csv')


p2 <- EnhancedVolcano(SGZ,
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


pdf(file = "/Users/sanachintamen/Analysis/ThesisPlots/SGZvsMGVolcano.pdf",   # The directory you want to save the file in
    width = 11, # The width of the plot in inches
    height = 15.5) # The height of the plot in inches

p3

dev.off()


pdf(file = "/Users/sanachintamen/Analysis/ThesisPlots/SGZvsMGVolcano2.pdf",   # The directory you want to save the file in
    width = 8, # The width of the plot in inches
    height = 13) # The height of the plot in inches

p3

dev.off()



p3 <- EnhancedVolcano(SGZvsMG,
                      lab = rownames(SGZvsMG),
                      x = 'avg_log2FC',
                      y = 'p_val_adj',
                      xlim = c(-1.5,3),
                      ylim = c(0,300),
                      pCutoff = 0.001,
                      FCcutoff = 0.4,
                      pointSize = 1.0,
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


p3
