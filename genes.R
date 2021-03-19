
#Ccl6, Ccl9 
Chemup <- c('Cxcl3', 'Cxcl1', 'Cxcl2', 'Ccl5', 'Pf4', 'Ccl6', 'Ccl9', 'Ppbp', 'Cxcl11', 'Ccl2', 'Ccl25', 'Cmtm3', 'Ybx1', 'Ccl7', 'Ccl12', 'Cxcl10', 'Cx3cl1')

#Cmtm7
Chemdown <- c('Cxcl5', 'Ccl3', 'Cmtm4', 'Cmtm7', 'Cklf', 'Cmtm2a', 'Ccl4', 'Cmtm8', 'Ccl8', 'Ccl27a', 'Cxcl14', 'Cxcl12', 'Cmtm5')

ECMup <- c('Thbs1', 'Lamc1', 'Lamb2', 'Mmp9', 'Nog', 'Tnfaip6', 'Adam17', 'Sepp1')
ECMdown <- c('Mmp3', 'Mmp10', 'Ecm1', 'Lepre1', 'Spp1', 'Gas6', 'Mgp')

EnzIn<- c('Camp', 'Serpinf2', 'Serpine1', 'Timp2', 'Cst3', 'Spink2')

#B4galt1
Memup <- c('Dll4', 'Sema3b', 'Cd70', 'Sema3f', 'Sema3a', 'Jag1', 'Sdc1', 'Dll3', 'B4galt1', 'Anxa1')
Memdown <- c('Sema3c', 'Sema4d', 'Efna1')

Cytoup <- c('Csf3', 'Ptn', 'Il1f6', 'Il1b','Il6','Il1a', 'OSsm', 'Lif', 'Il19', 'Tnfsf13b', 'Clcf1', 'Il27', 'Mdk', 'C3', 'Il24', 'Il1rn', 'Il18', 'Ebi3', 'Ifnb1', 'Ltb', 'Tnfsf14', 'Il10')
#Mif, Csf1, Aimp1 upregulated, Il6 downregulated
Cytodown <- c('Il1f9', 'Tnfsf9', 'Nampt', 'Mif', 'Tnf', 'Il15', 'Twsg1', 'Aimp1', 'Lta', 'Metrn', 'Kitl', 'Ctf1', 'Il16', 'Il4', 'Tnfsf15', 'Csf1', 'Tnfsf10')

TFup <- c('Prl3c1', 'Ereg', 'Areg', 'Inha', 'Artn', 'Bmp1', 'Fgf23', 'Nrtn', 'Fgf2', 'Prl2a1', 'Ins1', 'Vegfa','Ins2', 'Inhba', 'Prl3d1', 'Fgf14', 'Cntf', 'Gal','Fgf8', 'Flt3l','Fbrs', 'Fgf18', 'Gdf15', 'Tgfb1', 'Nenf', 'Hbegf', 'Pdgfa', 'Tgfb3', 'Pdgfb', 'Cyr61')
TFdown <- c('Nrg4', 'Nodal', 'Lefty2', 'Gpi', 'Hdgf', 'Prl13a1', 'Gdf9', 'Vegfb', 'Pdgfc', 'Egfl7', 'Grn', 'Fgf9', 'Lefty1', 'Bmp2', 'Cdnf', 'Nov', 'Grem1', 'Igf', 'Gdf3', 'Vegfc', 'Ctgf')

PepHup <- c('Vgf', 'Scg2', 'Amh', 'Cartpt', 'Nmb', 'Edn2', 'Npy', 'Sct', 'Retn', 'Fndc5', 'Cort', 'Nts', 'Npff', 'Rln3', 'Calca', 'Gnrh1', 'Prok2')
PepHdown <- c('Adm', 'Gcg', 'Stc2', 'Nppa', 'Pomc', 'Insl6', 'Edn1', 'Hgf', 'Gnl3')


Chemup <- c('Cxcl3', 'Cxcl1', 'Cxcl2', 'Ccl5', 'Pf4', 'Ccl6', 'Ccl9', 'Ppbp', 'Cxcl11', 'Ccl2', 'Ccl25', 'Cmtm3', 'Ybx1', 'Ccl7', 'Ccl12', 'Cxcl10', 'Cx3cl1')
Chemdown <- c('Cxcl5', 'Ccl3', 'Cmtm4', 'Cmtm7', 'Cklf', 'Cmtm2a', 'Ccl4', 'Cmtm8', 'Ccl8', 'Ccl27a', 'Cxcl14', 'Cxcl12', 'Cmtm5')
MISCup <- c('Dmkn', 'Hilpda', 'Hrg', 'Pecam1', 'Uts2b')
MISCdown <- c('Copa', 'Atxn10', 'Aggf1', 'Ngrn')

DotPlot(clusters, features = Chemup)
DotPlot(clusters, features = Chemdown)
DotPlot(clusters, features = ECMdown)
DotPlot(clusters, features = Memup)
DotPlot(clusters, features = Memdown, cols = 'Spectral')
DotPlot(clusters, features = Cytoup, cols = 'RdBu')
DotPlot(clusters, features = Cytodown, cols = 'RdBu')

DotPlot(clusters, features = TFdown, cols = 'RdBu')
DotPlot(clusters, features = MISCdown, cols = 'RdBu')


