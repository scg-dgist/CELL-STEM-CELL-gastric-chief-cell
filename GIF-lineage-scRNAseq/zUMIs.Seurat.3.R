# cl 4 remove
Idents(gif_seurat.combined.2) = gif_seurat.combined.2$seurat_clusters
gif_seurat.3 <- subset(gif_seurat.combined.2, idents = c(4), invert=TRUE)

gif_sce.3 <- as.SingleCellExperiment(gif_seurat.3, assay="RNA")
gif_sce.norm.3 <- norm_function(gif_sce.3)

gif_hvg.norm.3 <- hvg_function(gif_sce.norm.3)

gif_seurat.3 <- make_seurat_somi_function(gif_sce.norm.3, gif_hvg.norm.3)
gif_seurat.3 <- seurat_cluster_somi_function(gif_seurat.3, 20, 0.8)

DimPlot(gif_seurat.3)
DimPlot(gif_seurat.3, group.by = "Phase")

gif_seurat.combined.3 <- seurat_integ_function(gif_seurat.3, PCA=30, resolution=0.8)

gif_seurat.combined.3 <- FindClusters(gif_seurat.combined.3, resolution = 1)
DimPlot(gif_seurat.combined.3, pt.size=2, label=TRUE, label.size=8)

somi_featureplot(gif_seurat.combined.3, "Amy2a3")

DefaultAssay(gif_seurat.combined.3) = "RNA"
somi_featureplot(gif_seurat.combined.3, "Stmn1", col=col, size=2)

markers <- FindAllMarkers(gif_seurat.combined.3, logfc.threshold = 0)
markers$symbol = ensemblGenes[markers$gene, "external_gene_name"]
head(subset(markers, cluster==0 & avg_logFC > 0), n=20)

FeaturePlot(gif_seurat.combined.3, "log10_detected", pt.size=2)
VlnPlot(gif_seurat.combined.3, "pct_counts_ERCC")
VlnPlot(gif_seurat.combined.3, "pct_counts_MT")

gif_seurat.combined.3$gif_clusters = gif_seurat.combined.3$seurat_clusters
gif_seurat.combined.3$gif_clusters = gsub("3", "Gif_high", gif_seurat.combined.3$gif_clusters)
gif_seurat.combined.3$gif_clusters = gsub("0", "Gif_high", gif_seurat.combined.3$gif_clusters)
gif_seurat.combined.3$gif_clusters = gsub("1", "Gif_high", gif_seurat.combined.3$gif_clusters)
gif_seurat.combined.3$gif_clusters = gsub("2", "Gif_low", gif_seurat.combined.3$gif_clusters)

DimPlot(gif_seurat.combined.3, group.by = "gif_clusters")

save(gif_seurat.3, file=paste0(rdatadir, "Seurat/gif_seurat.3.RData"))
save(gif_seurat.combined.3, file=paste0(rdatadir, "Seurat/gif_seurat.combined.3.RData"))

#
# dimplot
clustern = length(unique(gif_seurat.combined.3$seurat_clusters))

library(LaCroixColoR)
gg_color_hue <- function(n) {
  hues = seq(15, 375, length = n + 1)
  hcl(h = hues, l = 65, c = 100)[1:n]
}
colors = list(
  gg_color_hue(clustern),
  c(as.character(lacroix_palette("CranRaspberry"))[1:3], as.character(lacroix_palette("Pure"))[1:3]),
  c(as.character(lacroix_palette("CranRaspberry"))[1], as.character(lacroix_palette("Pure"))[1]),
  c(lacroix_palette("Pure")[1], lacroix_palette("KiwiSandia")[3], lacroix_palette("KiwiSandia")[2], "lightgray")
)

for(i in 1:length(List)){
  d <- somi_dimplot2(gif_seurat.combined.3, group.by = List[i], size=3, cols = colors[[i]])
  ggsave(paste0(plotdir, "Seurat.3/dimplot/umap_", List[i], ".png"), plot=d, height=8, width=10)
}

# NT, D2
#
NT.cells = colnames(gif_seurat.combined.3[, gif_seurat.combined.3$condition == "NT"])
d <- somi_dimplot2(gif_seurat.combined.3, group.by = "seurat_clusters", cols = colors[[1]], cells=NT.cells, size=3)
ggsave(paste0(plotdir, "Seurat.3/dimplot/umap_seurat_clusters_NT.pdf"), plot=d, height=8, width=10)

d <- somi_dimplot2(gif_seurat.combined.3, group.by = "Phase", cols = colors[[4]], cells=NT.cells, size=3)
ggsave(paste0(plotdir, "Seurat.3/dimplot/umap_phase_NT.png"), plot=d, height=8, width=10)

D2.cells = colnames(gif_seurat.combined.3[, gif_seurat.combined.3$condition == "D2"])
d <- somi_dimplot2(gif_seurat.combined.3, group.by = "seurat_clusters", cols = colors[[1]], cells=D2.cells, size=3)
ggsave(paste0(plotdir, "Seurat.3/dimplot/umap_seurat_clusters_D2.pdf"), plot=d, height=8, width=10)

d <- somi_dimplot2(gif_seurat.combined.3, group.by = "Phase", cols = colors[[4]], cells=D2.cells, size=3)
ggsave(paste0(plotdir, "Seurat.3/dimplot/umap_phase_D2.png"), plot=d, height=8, width=10)


# featureplot
genes = c("Cblif", "GFP", "Gkn3", "Muc6", "Muc5ac", "Cd44", "Stmn1", "Mki67", "Lgr5", "Tnfrsf19", "Cdkn1c", "Tff2", "Atp4a", "Atp4b", "Kcne2", "Clu", "Sox9", "Amy2b", "Amy2a3")
genes = c("Pnliprp2", "Clps", "Furin")

DefaultAssay(gif_seurat.combined.3) = "RNA"
for(i in genes){
  f <- somi_featureplot(gif_seurat.combined.3, i, reduction="umap", col=col, size=3)
  ggsave(paste0(plotdir, "Seurat.3/featureplot/umap_expr_", i, ".pdf"), plot=f, height=8, width=10)
}

somi_featureplot_together <- function(seurat, feature, reduction="umap", col=NULL, size=1){
  
  feature_id = rownames(ensemblGenes[ensemblGenes$external_gene_name == feature,])
  print(c(feature, feature_id))
  
  if(is.null(col)){
    col = c("#EEEEEE", brewer.pal(9, "Reds"))
  }
  if(sum(feature_id %in% rownames(seurat)) == 0){
    print(paste0("no ", feature))
    
  }else{
    expr = GetAssayData(seurat, slot="data", assay="RNA")
    max = max(expr[feature_id,])
    min = min(expr[feature_id,])
    
    Idents(seurat) = seurat$condition
    
    f1 <- FeaturePlot(subset(seurat, idents = "NT"),
                      features = feature_id,
                      reduction = reduction,
                      order=TRUE, pt.size=size) +
      scale_color_gradientn(colours = col, breaks=c(min,(max-min)/2 + min,max), limits = c(min, max),
                            labels=c(round(min, digit=1),round((max-min)/2 + min, digit=1),round(max, digit=1))) +
      labs(title = feature) +
      theme(plot.title = element_text(size=50, face="bold.italic"),
            panel.grid.major=element_blank(),
            panel.grid.minor=element_blank(), 
            axis.line=element_blank(),
            axis.ticks=element_blank(),
            axis.title = element_blank(),
            panel.border = element_blank(),
            #legend.position = "none",
            legend.text=element_text(size=25), 
            legend.title=element_blank(),
            legend.key=element_blank(),
            axis.text = element_blank())
    
    ggsave(paste0(plotdir, "Seurat.3/NT.featureplot/umap_expr_NT_", feature, ".pdf"), plot=f1, height=8, width=10)
    
    f2 <- FeaturePlot(subset(seurat, idents = "D2"),
                      features = feature_id,
                      reduction = reduction,
                      order=TRUE, pt.size=size) +
      scale_color_gradientn(colours = col, breaks=c(min,(max-min)/2 + min,max), limits = c(min, max),
                            labels=c(round(min, digit=1),round((max-min)/2 + min, digit=1),round(max, digit=1))) +
      labs(title = feature) +
      theme(plot.title = element_text(size=50, face="bold.italic"),
            panel.grid.major=element_blank(),
            panel.grid.minor=element_blank(), 
            axis.line=element_blank(),
            axis.ticks=element_blank(),
            axis.title = element_blank(),
            panel.border = element_blank(),
            #legend.position = "none",
            legend.text=element_text(size=25), 
            legend.title=element_blank(),
            legend.key=element_blank(),
            axis.text = element_blank())
    ggsave(paste0(plotdir, "Seurat.3/D2.featureplot/umap_expr_D2_", feature, ".pdf"), plot=f2, height=8, width=10)
  }
}
for(g in genes){
  somi_featureplot_together(gif_seurat.combined.3, g, reduction = "umap", col=col, size=3)
}

# # dotplot
Genes = c("Cblif", "Gper1", "Bhlha15", "Pga5", "Cdkn1c", "Muc6", "Gkn3", "Mki67", "Stmn1",
          "Tff2", "Muc5ac", "Tff1", "Gkn1", "Gkn2")
DimPlot(gif_seurat.combined.3, group.by = "seurat_clusters", label=TRUE, label.size=8)
cl_order = c(3,0,1,2) %>% as.character

s <- somi_dotplot(gif_seurat.combined.3, rev(genes), cl_order, "seurat_clusters")
ggsave(paste0(plotdir, "Seurat.3/gif_dotplot.pdf"), plot=s, width=6, height=5)

Idents(gif_seurat.combined.3) = gif_seurat.combined.3$condition
s <- somi_dotplot(subset(gif_seurat.combined.3, idents = "NT"), rev(genes), cl_order, "seurat_clusters")
ggsave(paste0(plotdir, "Seurat.3/gif_NT_dotplot.pdf"), plot=s, width=6, height=5)
s <- somi_dotplot(subset(gif_seurat.combined.3, idents = "D2"), rev(genes), cl_order, "seurat_clusters")
ggsave(paste0(plotdir, "Seurat.3/gif_D2_dotplot.pdf"), plot=s, width=6, height=5)

cl_order = c("Gif_high", "Gif_low") %>% as.character

s <- somi_dotplot(gif_seurat.combined.3, rev(genes), cl_order, "gif_clusters")
ggsave(paste0(plotdir, "Seurat.3/gif_dotplot.pdf"), plot=s, width=5, height=5)

Idents(gif_seurat.combined.3) = gif_seurat.combined.3$condition
s <- somi_dotplot(subset(gif_seurat.combined.3, idents = "NT"), rev(genes), cl_order, "gif_clusters")
ggsave(paste0(plotdir, "Seurat.3/gif_NT_dotplot.pdf"), plot=s, width=6, height=5)
s <- somi_dotplot(subset(gif_seurat.combined.3, idents = "D2"), rev(genes), cl_order, "gif_clusters")
ggsave(paste0(plotdir, "Seurat.3/gif_D2_dotplot.pdf"), plot=s, width=6, height=5)

#####
genes = c("Cblif", "Gper1", "Bhlha15", "Pga5", "Cdkn1c", "Muc6", "Gkn3", "Mki67", "Stmn1",
          "Tff2", "Muc5ac", "Tff1", "Gkn1", "Gkn2")
genes = c("Hmgb2", "Top2a", "Smc2", "Ccne1", "Ccne2", "Cdt1", "Pcna", "Cdc45", "Mcm10", "Rad51", "Ccnb1", "Ccnb2", "Plk1", "Aurka")
genes = c("Ccnd1", "Ccnd2", "Ccnd3")

gif_high.cells <- colnames(gif_seurat.combined.3[, gif_seurat.combined.3$gif_clusters == "Gif_high"])
gif_low.cells <- colnames(gif_seurat.combined.3[, gif_seurat.combined.3$gif_clusters == "Gif_low"])

selected_geneids = ensemblGenes[match(genes, ensemblGenes$external_gene_name), ]$ensembl_gene_id
data <- gif_seurat.combined.3@assays$RNA@data[selected_geneids, ]
scale.data <- t(scale(t(data)))

mat <- matrix(nrow = nrow(scale.data),
              ncol = 2,
              dimnames = list(selected_geneids, c("gif_high", "gif_low")))
for(i in selected_geneids){
  mat[i, "gif_high"] = mean(scale.data[i, gif_high.cells])
  mat[i, "gif_low"] = mean(scale.data[i, gif_low.cells])
}
rownames(mat) = genes

mat[mat > 0.5] = 0.5
mat[mat < -0.5] = -0.5

prop.mat <- matrix(nrow = nrow(scale.data),
                   ncol = 2,
                   dimnames = list(selected_geneids, c("gif_high", "gif_low")))
for(i in selected_geneids){
  prop.mat[i, "gif_high"] = sum(data[i, gif_high.cells] > 0) / length(gif_high.cells) * 100
  prop.mat[i, "gif_low"] = sum(data[i, gif_low.cells] > 0) / length(gif_low.cells) * 100
}

dot.df <- data.frame(row.names = c("gif_high", "gif_low"))
dot.df <- cbind.data.frame(dot.df, t(mat))

data.df <- data.frame(genes = rep(genes, 2))
data.df$celltype = c(rep("gif_high", length(genes)), rep("gif_low", length(genes)))
data.df$expr = c(mat[,1], mat[,2])
data.df$prop = c(prop.mat[,1], prop.mat[,2])

min.value = min(data.df$expr)
max.value = max(data.df$expr)
abs.value = max(abs(min.value), max.value)
gd <- ggplot(data.df) + 
  geom_point(aes(x=celltype, y=genes, color=expr, size=prop)) +
  scale_y_discrete(limits = rev(genes)) + 
  scale_color_gradientn(colors = rev(brewer.pal(11, "RdBu"))[c(1:4, 6, 8:11)], limits = c(-abs.value, abs.value)) + 
  scale_size(range = c(0,9), limits = c(0, 100)) +
  theme_classic() +
  theme(text = element_text(size=20)) +
  ylab("") + xlab("")
gd
ggsave(paste0(plotdir, "Seurat.3/gif_dotplot_rescaled_representative_genes.pdf"), plot=gd, width=5, height=6)
ggsave(paste0(plotdir, "Seurat.3/gif_dotplot_rescaled_phase_cyclegenes.pdf"), plot=gd, width=5, height=6)

#####
genes = c("Hmgb2", "Top2a", "Smc2", "Ccne1", "Ccne2", "Cdt1", "Pcna", "Cdc45", "Mcm10", "Rad51", "Ccnb1", "Ccnb2", "Plk1", "Aurka")
genes = c("Ccnd1", "Ccnd2", "Ccnd3")

D2.gif_high.cells <- colnames(gif_seurat.combined.3[, gif_seurat.combined.3$gif_clusters == "Gif_high" & gif_seurat.combined.3$condition == "D2"])
D2.gif_low.cells <- colnames(gif_seurat.combined.3[, gif_seurat.combined.3$gif_clusters == "Gif_low" & gif_seurat.combined.3$condition == "D2"])
NT.gif_high.cells <- colnames(gif_seurat.combined.3[, gif_seurat.combined.3$gif_clusters == "Gif_high" & gif_seurat.combined.3$condition == "NT"])
NT.gif_low.cells <- colnames(gif_seurat.combined.3[, gif_seurat.combined.3$gif_clusters == "Gif_low" & gif_seurat.combined.3$condition == "NT"])

selected_geneids = ensemblGenes[match(genes, ensemblGenes$external_gene_name), ]$ensembl_gene_id
data <- gif_seurat.combined.3@assays$RNA@data[selected_geneids, ]
scale.data <- t(scale(t(data)))

mat <- matrix(nrow = nrow(scale.data),
              ncol = 4,
              dimnames = list(selected_geneids, c("NT.gif_high", "NT.gif_low", "D2.gif_high", "D2.gif_low")))
for(i in selected_geneids){
  mat[i, "NT.gif_high"] = mean(scale.data[i, NT.gif_high.cells])
  mat[i, "NT.gif_low"] = mean(scale.data[i, NT.gif_low.cells])
  mat[i, "D2.gif_high"] = mean(scale.data[i, D2.gif_high.cells])
  mat[i, "D2.gif_low"] = mean(scale.data[i, D2.gif_low.cells])
}
rownames(mat) = genes

mat[mat > 0.5] = 0.5
mat[mat < -0.5] = -0.5

prop.mat <- matrix(nrow = nrow(scale.data),
                   ncol = 4,
                   dimnames = list(selected_geneids, c("NT.gif_high", "NT.gif_low", "D2.gif_high", "D2.gif_low")))
for(i in selected_geneids){
  prop.mat[i, "NT.gif_high"] = sum(data[i, NT.gif_high.cells] > 0) / length(NT.gif_high.cells) * 100
  prop.mat[i, "NT.gif_low"] = sum(data[i, NT.gif_low.cells] > 0) / length(NT.gif_low.cells) * 100
  prop.mat[i, "D2.gif_high"] = sum(data[i, D2.gif_high.cells] > 0) / length(D2.gif_high.cells) * 100
  prop.mat[i, "D2.gif_low"] = sum(data[i, D2.gif_low.cells] > 0) / length(D2.gif_low.cells) * 100
}

dot.df <- data.frame(row.names = c("NT.gif_high", "NT.gif_low", "D2.gif_high", "D2.gif_low"))
dot.df <- cbind.data.frame(dot.df, t(mat))

data.df <- data.frame(genes = rep(genes, 4))
data.df$celltype = c(rep("NT.gif_high", length(genes)), rep("NT.gif_low", length(genes)), rep("D2.gif_high", length(genes)), rep("D2.gif_low", length(genes)))
data.df$expr = c(mat[,1], mat[,2], mat[,3], mat[,4])
data.df$prop = c(prop.mat[,1], prop.mat[,2], prop.mat[,3], prop.mat[,4])

min.value = min(data.df$expr)
max.value = max(data.df$expr)
abs.value = max(abs(min.value), max.value)

data.df$celltype = factor(data.df$celltype, levels = c("NT.gif_high","NT.gif_low", "D2.gif_high",  "D2.gif_low"))
gd <- ggplot(data.df) + 
  geom_point(aes(x=celltype, y=genes, color=expr, size=prop)) +
  scale_y_discrete(limits = rev(genes)) + 
  scale_color_gradientn(colors = rev(brewer.pal(11, "RdBu"))[c(1:4, 6, 8:11)], limits = c(-abs.value, abs.value)) + 
  scale_size(range = c(0,15), limits = c(0, 100)) +
  theme_classic() +
  theme(text = element_text(size=20),
        axis.text.x = element_text(angle=45, hjust=1)) +
  ylab("") + xlab("")
gd
ggsave(paste0(plotdir, "Seurat.3/gif_dotplot_rescaled_phase_4_order1.pdf"), plot=gd, width=7, height=5)

somi_featureplot(gif_seurat.combined.3, "Ccnd2")
#
DimPlot(gif_seurat.combined.3, pt.size=2, group.by = "seurat_clusters")
df <- data.frame(row.names = colnames(gif_seurat.combined.3),
                 #cluster = gif_seurat.combined.3$seurat_clusters,
                 cluster = gif_seurat.combined.3$gif_clusters,
                 condition = gif_seurat.combined.3$condition,
                 phase = gif_seurat.combined.3$Phase)
#df$cluster = factor(df$cluster, levels = c(3,0,1,2))
df$condition = factor(df$condition, levels = c("NT", "D2"))
v <- ggplot(df) + geom_bar(aes(x=condition, fill=phase), position = "fill") +
  scale_fill_manual(values = c(colors[[4]])) +
  theme_classic() #+
#facet_wrap(~cluster, ncol=5)
v
ggsave(paste0(plotdir, "Seurat.3/phase_condition_barplot.png"), plot=v, width=3, height=3)

v <- ggplot(df) + geom_bar(aes(x=condition, fill=phase), position = "fill") +
  scale_fill_manual(values = c(colors[[4]])) +
  theme_classic() +
  facet_wrap(~cluster, ncol=4)
v
ggsave(paste0(plotdir, "Seurat.3/phase_cl_condition_barplot.png"), plot=v, width=4, height=3)

v <- ggplot(df) + geom_bar(aes(x=cluster, fill=phase), position = "fill") +
  scale_fill_manual(values = c(colors[[4]])) +
  theme_classic()
v
ggsave(paste0(plotdir, "Seurat.3/phase_cl_barplot.png"), plot=v, width=3, height=3)

# barcode files
unique(gif_seurat.combined.3$batch)
batch = unique(gif_seurat.combined.3$batch)
samples = paste0("s00", 4:9)
samples

# 
write.csv(gif_seurat.combined.3@meta.data, file=paste0(rdatadir, "Seurat.3/gif_seurat.combined.3.metadata.csv"))
for(i in 1:6){
  cell.barcodes = colnames(gif_seurat.combined.3[, gif_seurat.combined.3$batch == batch[i]])
  cell.barcodes = gsub("s00.*_", "", cell.barcodes)
  write.table(cell.barcodes, file=paste0(rdatadir, "Seurat.3/cellbarcodes_", samples[i], ".tsv"), sep="\t", quote=FALSE, row.names=FALSE, col.names=FALSE)
}

#
Idents(gif_seurat.combined.3) = gif_seurat.combined.3$condition
gif.NT = subset(gif_seurat.combined.3, idents = "NT")
gif.D2 = subset(gif_seurat.combined.3, idents = "D2")

table(gif.NT$seurat_clusters)
table(gif.D2$seurat_clusters)

genes = ensemblGenes[ensemblGenes$external_gene_name %in% c("Mki67", "Stmn1"),]$ensembl_gene_id
seurat = gif.NT

cells = colnames(seurat[, seurat$seurat_clusters == 2])

mat = seurat@assays$RNA@counts[genes, cells]
table(colSums(mat) > 0)
