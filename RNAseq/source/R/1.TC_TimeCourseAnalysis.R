##### Code for Mouse Stomach Time-course RNA-seq data
##### Author: Somi Kim (ksomi47@dgist.ac.kr)
##### Last Update: 2021/02/17

### Time-course RNA-seq data
### 1. Time course sequencing data analysis


library(TCseq)
library(RColorBrewer)
library(biomaRt)
library(fpc)


##### load data

load("ensemblGenes_mmusculus_2019-12-12.RData")


### load raw count table
counts <- read.csv(file="../FILE/troy+bulkdata/raw_readcounts_SN6_SN8.csv", header=TRUE,
                   row.names=1, check.names = FALSE)

use_colnames = grep("Troy", colnames(counts), value=TRUE)

counts = counts[,use_colnames]
counts = counts[-c((nrow(counts)-4):nrow(counts)),]

rownames(counts) = gsub("\\..*", "", rownames(counts))


counts = counts[rownames(counts) %in% rownames(ensemblGenes),]


### make colData
sampleNames <- c(
  paste0("d0_", 1:6),
  paste0("d1_", 1:3),
  paste0("d3_", 1:3),
  paste0("d7_", 1:3),
  paste0("d14_", 1:2)
)
sampleDays <- c(rep("d0", 6), rep("d1", 3), rep("d3", 3), rep("d7", 3), rep("d14", 2))
sampleGroups <- c(rep(1,6), rep(2:4,each=3), rep(5,2))

coldata <- data.frame(names = sampleNames,
                      days = sampleDays,
                      group = sampleGroups)


### make genomic reference
ensembl <- useMart(biomart = "ensembl", 
                   dataset="mmusculus_gene_ensembl", 
                   host="http://sep2019.archive.ensembl.org")

ensemblRef <- getBM(attributes=c('ensembl_gene_id', 'external_gene_name',  'chromosome_name',
                                   'start_position', 'end_position', 'gene_biotype'), 
                    mart=ensembl)
rownames(ensemblRef) <- ensemblRef[,1]

TC_ref = ensemblRef[,c('ensembl_gene_id',  'chromosome_name', 'start_position', 'end_position')]
head(TC_ref)
colnames(TC_ref) = c("id", "chr", "start", "end")

rm(ensemblRef)


### run TCseq
rownames(counts) = gsub("[.][0-9]+", "", rownames(counts))
TC_ref = TC_ref[rownames(counts),]

tca <- TCA(design = coldata,
           counts = counts,
           genomicFeature = TC_ref)

tca <- DBanalysis(tca)

tca <- timecourseTable(tca, value = "FC", norm.method = FALSE, filter = TRUE, abs.fold=1, pvalue.threshold = 0.05)

i=6
set.seed(i)
tca <- timeclust(tca, algo="cm", k=i, standardize=TRUE)

p <- timeclustplot(tca, value = "z-score(RPKM)", cols=3,
                   membership.color = brewer.pal(9, "YlGnBu"))

day1_up_genes = c(levels(p[[2]]$data$group), levels(p[[4]]$data$group))
day1_down_genes = c(levels(p[[5]]$data$group), levels(p[[6]]$data$group))


#save(tca, "tca.RData")