.libPaths(Sys.getenv("R_LIBS_STOMACH"))

library(SingleCellExperiment)

load("D:/somi_function/ensemblGenes_mmusculus_2019-12-12.RData")

rdatadir = "rdata3/"
plotdir = "plots3/"

samples = paste0("s00", 4:9)
samples

cbs <- read.table(file="bc_celseq2.tsv", sep = "\t")
cb.control = cbs$V1[c(357:360, 381:384)]

for(sample in samples){
  
  count = readRDS(paste0("D:/Project/ISTHMUS/gif_sortseq/zumis/", sample, "/GIF_lineage_tracing.dgecounts.rds"))
  sce = SingleCellExperiment(assays = list(counts = count$umicount$exon$all))
  
  cells.control = cb.control
  is.control = grepl(paste(cells.control, collapse="|"), colnames(sce))
  
  sce = sce[, !is.control]
  
  colnames(sce) = paste0(sample, "_", colnames(sce))
  
  assign(paste0("gif_sce.", sample), sce)
  
}
rm(sce)
rm(count)

for(sample in samples){
  save(list = paste0("gif_sce.", sample), file=paste0(rdatadir, "sce/gif_sce.", sample))
}
