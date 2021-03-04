#To filter for p values
topDiffGenes <- function (genelist) 
{
  return(genelist < 0.01)
}

#To filter for upregulated genes and provide appropriate input for TopGO
UpGenes <- function (genelist) 
{
  genes_up <- filter(genelist, avg_log2FC > 0)
  gene <- genes_up$gene
  genes_up <- as.numeric(genes_up$p_val_adj)
  names(genes_up) <- gene
  rm(gene)
  return(genes_up)
}


#To filter for downregulated genes and provide appropriate input for TopGO
DownGenes <- function (genelist) 
{
  genes_down <- filter(genelist, avg_log2FC < 0)
  gene <- genes_down$gene
  genes_down <- as.numeric(genes_down$p_val_adj)
  names(genes_down) <- gene
  rm(gene)
  return(genes_down)
}



# To perform GO with TopGOO
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
