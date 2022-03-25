load("D:/somi_function/ensemblGenes_mmusculus_2019-12-12.RData")

rdatadir = "rdata3/"
plotdir = "plots3/"

library(scran)
library(scMerge)

# load sce
samples = paste0("s00", 4:9)
for(i in samples){
  objectname = paste0("sce.qc.", i)
  load(paste0(rdatadir, "sce.qc/", objectname, ".RData"))
}

# merge singlecellexperiment objects
for(sample in samples){
  sce = eval(parse(text = paste0("gif_sce.qc.", sample)))
  counts(sce) = as.matrix(counts(sce))
  
  assign(paste0("gif_sce.qc.", sample), sce)
}

sce_list = list(
  gif_sce.qc.s004,
  gif_sce.qc.s005,
  gif_sce.qc.s006,
  gif_sce.qc.s007,
  gif_sce.qc.s008,
  gif_sce.qc.s009
)

coldataNames = colnames(colData(gif_sce.qc.s004))
gif_sce = sce_cbind(sce_list, method="union", cut_off_batch = 0, cut_off_overall = 0, exprs = "counts",
                    colData_names = coldataNames, batch_names = c(paste("NT", 1:3, sep="-"), paste("D2", 1:3, sep="-")))
gif_sce$batch = gsub("s004", "control-1", gif_sce$batch)
gif_sce$batch = gsub("s005", "control-2", gif_sce$batch)
gif_sce$batch = gsub("s006", "control-3", gif_sce$batch)

gif_sce$batch = gsub("s007", "injury-1", gif_sce$batch)
gif_sce$batch = gsub("s008", "injury-2", gif_sce$batch)
gif_sce$batch = gsub("s009", "injury-3", gif_sce$batch)

gif_sce$condition = gsub("-[0-9]", "", gif_sce$batch)

norm_function <- function(sce){
  
  clusters <- quickCluster(sce, use.ranks = TRUE)
  sce <- computeSumFactors(sce, cluster=clusters)
  
  rm(clusters)
  
  sce <- normalize(sce)
  
  rownames(sce@assays@data$logcounts) = rownames(sce)
  colnames(sce@assays@data$logcounts) = colnames(sce)
  
  return(sce)
}

gif_sce.norm <- norm_function(gif_sce)

hvg_function <- function(sce){
  var.fit <- trendVar(sce, method="spline", parametric=FALSE, use.spikes=FALSE)
  var.out <- decomposeVar(sce, var.fit)
  var.out.order = var.out[order(var.out$FDR),]
  
  #hvg <- var.out[1:500,]
  hvg <- var.out[which(var.out$FDR < 0.05 & var.out$bio > 0.05),]
  print(nrow(hvg))
  
  return(hvg)
}

gif_hvg.norm <- hvg_function(gif_sce.norm) # 239

