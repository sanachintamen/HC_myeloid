

######
library(dplyr)
library(tidyr)

meta <- clusters@meta.data
meta <- as.data.frame(meta)

run1 <- filter(meta, orig.ident == 'SS01')
run1 <- as.data.frame(table(run1$tree.ident))
run <- '1'
length(run) <- 14
run <- replace_na(run, '1')
run1 <- cbind(run, run1)
colnames(run1) <-c('RunID', 'ClusterID', 'CellCount')
totalcells <- sum(run1$CellCount)
Percent <- ((run1$CellCount)/totalcells) * 100
Percent <- round(Percent, digits=2)
run1 <- cbind(run1, Percent)


run2<- filter(meta, orig.ident == 'SS02')
run2 <- as.data.frame(table(run2$tree.ident))
run <- '2'
length(run) <- 14
run <- replace_na(run, '2')
run2 <- cbind(run, run2)
colnames(run2) <-c('RunID', 'ClusterID', 'CellCount')
totalcells <- sum(run2$CellCount)
Percent <- ((run2$CellCount)/totalcells) * 100
Percent <- round(Percent, digits=2)
run2 <- cbind(run2, Percent)


run3<- filter(meta, orig.ident == 'SS03')
run3 <- as.data.frame(table(run3$tree.ident))
run <- '3'
length(run) <- 14
run <- replace_na(run, '3')
run3 <- cbind(run, run3)
colnames(run3) <-c('RunID', 'ClusterID', 'CellCount')
totalcells <- sum(run3$CellCount)
Percent <- ((run3$CellCount)/totalcells) * 100
Percent <- round(Percent, digits=2)
run3 <- cbind(run3, Percent)

runs <- rbind(run1, run2, run3)

p <- ggplot(runs, aes(fill= ClusterID, y= Percent, x= RunID)) +
  geom_bar( stat="identity", width = 0.4) +
  scale_fill_manual(values = combined)+
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        panel.background = element_blank(), 
        axis.line = element_line(colour = "black")) 
p
