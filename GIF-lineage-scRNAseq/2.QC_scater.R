load("D:/somi_function/ensemblGenes_mmusculus_2019-12-12.RData")

rdatadir = "rdata3/"
plotdir = "plots3/"

library(scater)
library(dplyr)
library(RColorBrewer)

samples = paste0("s00", 4:9)
samples

# S 4,5,6 = NT
# S 7,8,9 = D2

qc_hist_function <- function(sce, v=0, column, name, path){
  
  hist_column = sce[[column]]
  
  png(paste0(path, "hist_", column, ".png"))
  hist(
    hist_column,
    breaks = 100,
    xlab = name,
    main="")
  if(v != 0){
    abline(v = v, col="red")
  }
  dev.off()
}
qc_dimplot_function <- function(df, column, name, path){
  
  df = df[order(df[,column]),]
  dd <- ggplot(df) +
    geom_point(aes(x=PC1, y=PC2, color = eval(parse(text = column)))) +
    scale_colour_gradientn(colours = brewer.pal(9, "Reds"),
                           limits = c(min(df[,column]), max(df[,column]))) +
    ggtitle(name) +
    theme_bw() +
    theme(plot.title = element_text(size = 15, face="italic", hjust = 0.5, margin = margin(0,0,4,0,"mm")),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          panel.border = element_blank(),
          axis.line = element_blank(),
          axis.ticks = element_blank(),
          axis.text = element_blank(),
          axis.title = element_text(size=15, face="italic"),
          legend.title = element_blank(),
          legend.key = element_blank())
  ggsave(paste0(path, "dimplot_", column, ".png"), plot=dd, width=6, height=5)
  
}
qc_function <- function(sce, total_counts_cutoff=0, total_features_cutoff=0, pct_counts_ERCC_cutoff=0,
                        save = FALSE, objectname=NULL){
  
  if(sum(is.null(objectname)) == 1){
    objectname = deparse(substitute(sce))
    objectname = gsub("raw", "", objectname)
  }
  
  path = paste0(plotdir, strsplit(objectname, "_")[[1]][1], "_qc/")
  if(sum(dir.exists(path)) == 0){
    dir.create(path)
  }
  
  # mtGenes = ensemblGenes[ensemblGenes[,3] == "MT",]
  # is.mito = rownames(sce) %in% mtGenes[,2]
  # print(sum(is.mito))
  
  is.ERCC = grepl("ERCC", rownames(sce))
  print(sum(is.ERCC))
  
  mt.genes = ensemblGenes[ensemblGenes$chromosome_name == "MT", ]$ensembl_gene_id
  is.MT = rownames(sce) %in% mt.genes
  print(sum(is.MT))
  
  #cells.control = c(357:360, 381:384)
  #is.control = grepl(paste(cells.control, collapse="|"), colnames(sce))
  #print(sum(is.control))
  
  sce <- calculateQCMetrics(sce,
                            feature_controls = list(ERCC = is.ERCC,
                                                    MT = is.MT))
  sce$log10_sum = log10(sce$total_counts + 1)
  sce$log10_detected = log10(sce$total_features_by_counts + 1)
  
  ### plot function
  
  coldata = c("log10_sum", "log10_detected", "pct_counts_ERCC")
  columnName = c("log10 total umi counts", "log10 total feature counts", "pct counts ERCC")
  cutoffs = c(total_counts_cutoff, total_features_cutoff, pct_counts_ERCC_cutoff)
  
  # run pca
  # vector <-  c(unique(colnames(sce@colData)))
  vector <- c("total_counts", "total_features_by_counts", "pct_counts_ERCC")
  sce <- runColDataPCA(sce, ncomponents=5, variables = vector)
  
  # hist
  for(i in 1:length(coldata)){
    qc_hist_function(sce, v=cutoffs[i], coldata[i], columnName[i], path=path)
  }
  
  df <- as.data.frame(sce@int_colData$reducedDims$PCA_coldata)
  df$log10_sum <- sce$log10_sum
  df$log10_detected <- sce$log10_detected
  df$pct_counts_ERCC <- sce$pct_counts_ERCC
  
  # pca plot
  for(i in 1:length(coldata)){
    qc_dimplot_function(df, coldata[i], columnName[i], path=path)
  }
  
  filter_by_total_counts = sce$log10_sum > total_counts_cutoff
  filter_by_feature_counts = sce$log10_detected > total_features_cutoff
  filter_by_pct_counts_ERCC = sce$pct_counts_ERCC < pct_counts_ERCC_cutoff
  
  sce$use <- (
    filter_by_feature_counts &
      filter_by_total_counts &
      filter_by_pct_counts_ERCC
  )
  print(table(sce$use))
  
  ### save function
  if(sum(save) == 1){
    
    sce.qc <- sce[, sce$use]
    
    pp <- plotReducedDim(sce, "PCA_coldata",
                   colour_by = "use",
                   size_by = "total_features_by_counts") +
      theme(panel.grid.major=element_blank(),
            panel.grid.minor=element_blank(),
            axis.line=element_blank(),
            axis.ticks=element_blank(),
            axis.title = element_blank(),
            panel.border = element_blank(),
            legend.text=element_text(size=13),
            legend.key=element_blank(),
            axis.text = element_blank())
    ggsave(paste0(path, "dimplot_use.png"), plot=pp, width=6, height=5)
    
    print(paste0("sce.qc is generated. There are ", ncol(sce.qc), " cells finally."))
    return(sce.qc)
    
  }else{
    return(sce)
  }
}

for(i in 1:6){
  objectname = paste0("gif_sce.", samples[i])
  rawsce <- eval(parse(text=objectname))
  
  sce1 <- qc_function(rawsce, total_counts_cutoff = 2.5, pct_counts_ERCC_cutoff = 20, save=FALSE, objectname = samples[i])
  sce.qc <- qc_function(rawsce, total_counts_cutoff = 2.5, pct_counts_ERCC_cutoff = 20, save=TRUE, objectname = samples[i])
  
  assign(paste0("gif_sce.", samples[i]), sce1)
  assign(paste0("gif_sce.qc.", samples[i]), sce.qc)
}
