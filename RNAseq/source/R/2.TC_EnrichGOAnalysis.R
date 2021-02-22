##### Code for Mouse Stomach Time-course RNA-seq data
##### Author: Somi Kim (ksomi47@dgist.ac.kr)
##### Last Update: 2021/02/17

### Time-course RNA-seq data
### 2. Functional enrichment analysis of DEGs at 1 dpi


library(DESeq2)
library(topGO)
library(plyr)
library(dplyr)
library(pheatmap)
library(RColorBrewer)


##### load data


load("tca.RData")

day1_up_genes = c(levels(p[[2]]$data$group), levels(p[[4]]$data$group))
day1_down_genes = c(levels(p[[5]]$data$group), levels(p[[6]]$data$group))


### GO analysis

backgroundGenes <- names(tca@clusterRes@cluster)


allGene = factor(as.integer(backgroundGenes %in% day1_up_genes))
names(allGene) = backgroundGenes

tgd = new("topGOdata", ontology="BP", allGenes=allGene, nodeSize=5,
          annot=annFUN.org, mapping="org.Mm.eg.db", ID="Ensembl")
resultTopGO.elim = runTest(tgd, algorithm="elim", statistic="Fisher")
goTable_up = GenTable(tgd, Fisher.elim=resultTopGO.elim,
                      orderBy="Fisher.classic", topNodes = topnodes)

goTable_up$Fisher.elim = as.numeric(goTable_up$Fisher.elim)
goTable_up$Fisher.elim[is.na(goTable_up$Fisher.elim)] = 0
goTable_up$log10P = -log10(goTable_up$Fisher.elim)


allGene = factor(as.integer(backgroundGenes %in% day1_down_genes))
names(allGene) = backgroundGenes

tgd = new("topGOdata", ontology="BP", allGenes=allGene, nodeSize=5,
          annot=annFUN.org, mapping="org.Mm.eg.db", ID="Ensembl")
resultTopGO.elim = runTest(tgd, algorithm="elim", statistic="Fisher")
goTable_down = GenTable(tgd, Fisher.elim=resultTopGO.elim,
                        orderBy="Fisher.classic", topNodes = topnodes)

goTable_down$Fisher.elim = as.numeric(goTable_down$Fisher.elim)
goTable_down$Fisher.elim[is.na(goTable_down$Fisher.elim)] = 0
goTable_down$log10P = -log10(goTable_down$Fisher.elim)

