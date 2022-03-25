.libPaths(Sys.getenv("R_LIBS_STOMACH"))

load("D:/somi_function/ensemblGenes_mmusculus_2019-12-12.RData")

rdatadir = "rdata3/"
plotdir = "plots3/"

library(KernelKnn)
library(ggplot2)
library(Seurat)
library(dplyr)
library(LaCroixColoR)

load(paste0(rdatadir, "seurat/gif_seurat.combined.3.RData"))
load("D:/Project/ISTHMUS/pgc_scRNA/data/Seurat/pgc_seurat_cc2.RData")

run_kernelknn <- function(hvg, ref_mat, test_mat, ref_embeddings, k, method="pearson_correlation", scale=FALSE, size=1.5){
  
  hvg = hvg[hvg %in% rownames(test_mat)]
  hvg = hvg[hvg %in% rownames(ref_mat)]
  
  ref_mat = as.matrix(ref_mat[hvg,])
  test_mat = as.matrix(test_mat[hvg,])
  
  if(sum(scale) == 1){
    ref_mat = t(scale(t(ref_mat)))
    test_mat = t(scale(t(test_mat)))
  }
  
  indexN <- knn.index.dist(t(ref_mat), t(test_mat), k=k, threads=4, method=method)
  
  iN2 <- indexN$test_knn_idx
  rownames(iN2) <- colnames(test_mat)
  
  iN3 <- apply(iN2, 2, function(x) colnames(ref_mat)[x])
  idx <- apply(iN3, 2, function(x) ref_embeddings[x,1])
  idy <- apply(iN3, 2, function(x) ref_embeddings[x,2])
  
  prjx <- rowMeans(idx)
  prjy <- rowMeans(idy)
  names(prjx) <- rownames(iN2)
  names(prjy) <- rownames(iN2)
  
  df=data.frame(x=c(ref_embeddings[,1],prjx), 
                y=c(ref_embeddings[,2],prjy), 
                group = c(rep("ref",nrow(ref_embeddings)),
                          rep("test",length(prjx))))
  
  g <- ggplot(df, aes(x=x, y=y, color = group)) +
    geom_point(size=size) +
    guides(colour = guide_legend(override.aes = list(size=10)))+
    theme_bw() +
    theme(text = element_text(size = 20),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          panel.border = element_blank(),
          axis.line = element_blank(),
          axis.ticks = element_blank(),
          axis.title = element_blank(),
          axis.text = element_blank(),
          legend.title = element_blank(),
          legend.key = element_blank())
  print(g)
  return(df)
}

# gif_seurat
DefaultAssay(gif_seurat.combined.3) = "integrated"
hvg.gif = VariableFeatures(gif_seurat.combined.3)
mat.gif <- gif_seurat.combined.3@assays$RNA@data

#hvg.gif = ensemblGenes[match(hvg.gif, ensemblGenes$external_gene_name),]$ensembl_gene_id
#hvg.gif = hvg.gif[!is.na(hvg.gif)]
#rownames(mat.gif) = ensemblGenes[match(rownames(mat.gif), ensemblGenes$external_gene_name),]$ensembl_gene_id
#mat.gif = mat.gif[!is.na(rownames(mat.gif)),]

# pgc_seurat_cc2
hvg.pgc = VariableFeatures(pgc_seurat_cc2)
mat.pgc = pgc_seurat_cc2@assays$RNA@data

umap.pgc = pgc_seurat_cc2@reductions$umap@cell.embeddings

# perform mapping
df.gif_to_pgc <- run_kernelknn(hvg.gif, mat.pgc, mat.gif, umap.pgc, k=10, scale=TRUE)

# df.gif_to_pgc$celltype = "ref"
# df.gif_to_pgc[colnames(mat.gif),]$celltype = gif_seurat.combined.3$seurat_clusters %>% as.character()
# df.gif_to_pgc$celltype = factor(df.gif_to_pgc$celltype, levels = c(0:(clustern-1), "ref"))
df.gif_to_pgc$celltype = "ref"
df.gif_to_pgc[colnames(mat.gif),]$celltype = gif_seurat.combined.3$gif_clusters %>% as.character()
df.gif_to_pgc$celltype = factor(df.gif_to_pgc$celltype, levels = c("Gif_high", "Gif_low", "ref"))
clustern = 2

df.gif_to_pgc$cond = ""
df.gif_to_pgc[colnames(mat.gif),]$cond = gif_seurat.combined.3$condition %>% as.character()

df.gif_to_pgc$celltype_cond = paste(df.gif_to_pgc$celltype, df.gif_to_pgc$cond, sep="_")
df.gif_to_pgc$celltype_cond = gsub("ref_", "ref", df.gif_to_pgc$celltype_cond)

df.gif_to_pgc$phase = "ref"
df.gif_to_pgc[colnames(mat.gif),]$phase = as.character(gif_seurat.combined.3$Phase)
df.gif_to_pgc$phase = factor(df.gif_to_pgc$phase, levels = c("G1", "S", "G2M", "Undecided", "ref"))

geneid = ensemblGenes[ensemblGenes$external_gene_name == "Cdkn1c", ]$ensembl_gene_id
df.gif_to_pgc$p57 = NA
df.gif_to_pgc[colnames(mat.gif), ]$p57 = gif_seurat.combined.3@assays$RNA@data[geneid,]
df.gif_to_pgc$p57 = as.numeric(df.gif_to_pgc$p57)

geneid = ensemblGenes[ensemblGenes$external_gene_name == "Mki67", ]$ensembl_gene_id
df.gif_to_pgc$ki67 = NA
df.gif_to_pgc[colnames(mat.gif), ]$ki67 = gif_seurat.combined.3@assays$RNA@data[geneid,]
df.gif_to_pgc$ki67 = as.numeric(df.gif_to_pgc$ki67)

geneid = ensemblGenes[ensemblGenes$external_gene_name == "Cblif", ]$ensembl_gene_id
df.gif_to_pgc$gif = NA
df.gif_to_pgc[colnames(mat.gif), ]$gif = gif_seurat.combined.3@assays$RNA@data[geneid,]
df.gif_to_pgc$gif = as.numeric(df.gif_to_pgc$gif)

ggg <- ggplot(df.gif_to_pgc[order(df.gif_to_pgc$gif, na.last = FALSE),]) +
  geom_point(aes(x=x, y=y, color=gif, size=celltype))+
  #guides(colour = guide_legend(override.aes = list(size=10)))+
  scale_size_manual(values = c(rep(2, clustern),1.5)) + 
  scale_color_gradientn(colors=col, na.value = "grey50") +
  theme_bw() +
  theme(text = element_text(size = 20),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        axis.line = element_blank(),
        axis.ticks = element_blank(),
        axis.title = element_blank(),
        axis.text = element_blank(),
        legend.title = element_blank(),
        legend.key = element_blank())
ggsave(paste0(plotdir, "Seurat.3/projection_gif.png"), plot=ggg,  width=8, height=6)

gg_color_hue <- function(n) {
  hues = seq(15, 375, length = n + 1)
  hcl(h = hues, l = 65, c = 100)[1:n]
}
g <- ggplot(df.gif_to_pgc) +
  geom_point(aes(x=x, y=y, color=celltype, size=celltype))+
  guides(colour = guide_legend(override.aes = list(size=10)))+
  scale_size_manual(values = c(rep(2, clustern),1.5)) + 
  scale_color_manual(values = c(gg_color_hue(clustern), "lightgray")) +
  theme_bw() +
  theme(text = element_text(size = 20),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        axis.line = element_blank(),
        axis.ticks = element_blank(),
        axis.title = element_blank(),
        axis.text = element_blank(),
        legend.title = element_blank(),
        legend.key = element_blank())
g
ggsave(paste0(plotdir, "Seurat.3/projection_gif_to_pgc.pdf"), plot=g, height=6, width=7.5)

for(i in 0:(clustern-1)){
  df = df.gif_to_pgc[df.gif_to_pgc$celltype %in% c(as.character(i), "ref"),]
  g <- ggplot(df) +
    geom_point(aes(x=x, y=y, color=celltype, size=celltype))+
    guides(colour = guide_legend(override.aes = list(size=10)))+
    scale_size_manual(values = c(2,1.5)) + 
    scale_color_manual(values = c(gg_color_hue(clustern)[i+1], "lightgray")) +
    #scale_color_manual(values = c("lightgray", as.character(lacroix_palette("CranRaspberry"))[1], as.character(lacroix_palette("Pure"))[2])) +
    theme_bw() +
    theme(text = element_text(size = 20),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          panel.border = element_blank(),
          axis.line = element_blank(),
          axis.ticks = element_blank(),
          axis.title = element_blank(),
          axis.text = element_blank(),
          legend.title = element_blank(),
          legend.key = element_blank())
  ggsave(paste0(plotdir, "Seurat.3/projection_gif_to_pgc_cl", i, ".pdf"), plot=g, height=6, width=7.5)
}

celltype = c("Gif_high", "Gif_low")
for(i in 0:(clustern-1)){
  df = df.gif_to_pgc[df.gif_to_pgc$celltype %in% c(celltype[i+1], "ref"),]
  g <- ggplot(df) +
    geom_point(aes(x=x, y=y, color=cond, size=celltype))+
    guides(colour = guide_legend(override.aes = list(size=10)))+
    scale_size_manual(values = c(2,1.5)) + 
    #scale_color_manual(values = c(gg_color_hue(7)[i+1], "lightgray")) +
    scale_color_manual(values = c("lightgray", as.character(lacroix_palette("CranRaspberry"))[1], as.character(lacroix_palette("Pure"))[2])) +
    theme_bw() +
    guides(size="none") + 
    theme(text = element_text(size = 20),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          panel.border = element_blank(),
          axis.line = element_blank(),
          axis.ticks = element_blank(),
          axis.title = element_blank(),
          axis.text = element_blank(),
          legend.title = element_blank(),
          legend.key = element_blank())
  ggsave(paste0(plotdir, "Seurat.3/projection_gif_to_pgc_gifcl_cond", celltype[i+1], ".pdf"), plot=g, height=6, width=7.5)
}

for(i in 0:(clustern-1)){
  df = df.gif_to_pgc[df.gif_to_pgc$celltype %in% c(as.character(i), "ref"),]
  g <- ggplot(df) +
    geom_point(aes(x=x, y=y, color=phase, size=celltype))+
    guides(colour = guide_legend(override.aes = list(size=10)))+
    scale_size_manual(values = c(2,1.5)) + 
    #scale_color_manual(values = c(gg_color_hue(7)[i+1], "lightgray")) +
    scale_color_manual(values = c(colors[[4]], "lightgray")) +
    theme_bw() +
    theme(text = element_text(size = 20),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          panel.border = element_blank(),
          axis.line = element_blank(),
          axis.ticks = element_blank(),
          axis.title = element_blank(),
          axis.text = element_blank(),
          legend.title = element_blank(),
          legend.key = element_blank())
  ggsave(paste0(plotdir, "Seurat.3/projection_gif_to_pgc_cl_phase", i, ".pdf"), plot=g, height=6, width=7.5)
}
