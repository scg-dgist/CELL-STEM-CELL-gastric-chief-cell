
gif_seurat.combined

Idents(gif_seurat.combined) = gif_seurat.combined$seurat_clusters
gif_seurat.2 <- subset(gif_seurat.combined, idents = c(5), invert=TRUE)
gif_seurat.2

gif_sce.2 <- as.SingleCellExperiment(gif_seurat.2, assay="RNA")
gif_sce.norm.2 <- norm_function(gif_sce.2)

gif_hvg.norm.2 <- hvg_function(gif_sce.norm.2)

gif_seurat.2 <- make_seurat_somi_function(gif_sce.norm.2, gif_hvg.norm.2)
gif_seurat.2 <- seurat_cluster_somi_function(gif_seurat.2, 20, 0.8)

DimPlot(gif_seurat.2)
DimPlot(gif_seurat.2, group.by = "Phase")

gif_seurat.combined.2 <- seurat_integ_function(gif_seurat.2, PCA=20, resolution=0.8)

DimPlot(gif_seurat.combined.2)

DefaultAssay(gif_seurat.combined.2) = "RNA"
somi_featureplot(gif_seurat.combined.2, "Amy2b", col=col, size=2)

markers.2 <- FindAllMarkers(gif_seurat.combined.2, logfc.threshold = 0)
markers.2$symbol = ensemblGenes[markers.2$gene, "external_gene_name"]
head(subset(markers.2, cluster==4), n=20)

save(gif_seurat.2, file=paste0(rdatadir, "Seurat/gif_seurat.2.RData"))
save(gif_seurat.combined.2, file=paste0(rdatadir, "Seurat/gif_seurat.combined.2.RData"))

# dimplot
clustern = length(unique(gif_seurat.combined.2$seurat_clusters))

library(LaCroixColoR)
gg_color_hue <- function(n) {
  hues = seq(15, 375, length = n + 1)
  hcl(h = hues, l = 65, c = 100)[1:n]
}
colors = list(
  gg_color_hue(clustern),
  c(as.character(lacroix_palette("CranRaspberry"))[1:3], as.character(lacroix_palette("Pure"))[1:3]),
  c(as.character(lacroix_palette("CranRaspberry"))[1], as.character(lacroix_palette("Pure"))[1]),
  c(lacroix_palette("Pure")[1], lacroix_palette("KiwiSandia")[c(2,3)], "lightgray")
)

for(i in 1:length(List)){
  d <- somi_dimplot2(gif_seurat.combined.2, group.by = List[i], size=3, cols = colors[[i]])
  ggsave(paste0(plotdir, "Seurat.2/dimplot/umap_", List[i], ".png"), plot=d, height=8, width=10)
}

# featureplot
genes = c("Cblif", "GFP", "Gkn3", "Muc6", "Muc5ac", "Cd44", "Stmn1", "Mki67", "Lgr5", "Tnfrsf19", "Cdkn1c", "Tff2", "Atp4a", "Atp4b", "Kcne2", "Clu", "Sox9", "Amy2b", "Amy2a3")

DefaultAssay(gif_seurat.combined.2) = "RNA"
for(i in genes){
  f <- somi_featureplot(gif_seurat.combined.2, i, reduction="umap", col=col, size=3)
  ggsave(paste0(plotdir, "Seurat.2/featureplot/umap_expr_", i, ".png"), plot=f, height=8, width=10)
}

# NT, D2
#
NT.cells = colnames(gif_seurat.combined.2[, gif_seurat.combined.2$condition == "NT"])
d <- somi_dimplot2(gif_seurat.combined.2, group.by = "seurat_clusters", cols = colors[[1]], cells=NT.cells, size=3)
ggsave(paste0(plotdir, "Seurat.2/dimplot/umap_seurat_clusters_NT.png"), plot=d, height=8, width=10)

D2.cells = colnames(gif_seurat.combined.2[, gif_seurat.combined.2$condition == "D2"])
d <- somi_dimplot2(gif_seurat.combined.2, group.by = "seurat_clusters", cols = colors[[1]], cells=D2.cells, size=3)
ggsave(paste0(plotdir, "Seurat.2/dimplot/umap_seurat_clusters_D2.png"), plot=d, height=8, width=10)

#
# proportion bar plot
Idents(gif_seurat.combined.2) = gif_seurat.combined.2$condition
gif_seurat.NT = subset(gif_seurat.combined.2, idents = "NT")
gif_seurat.D2 = subset(gif_seurat.combined.2, idents = "D2")

# control
gif_seurat.NT

gif_sce.NT <- as.SingleCellExperiment(gif_seurat.NT, assay = "RNA")
gif_sce.NT <- norm_function(gif_sce.NT)
gif_hvg.NT <- hvg_function(gif_sce.NT)

gif_seurat.NT <- make_seurat_somi_function(gif_sce.NT, gif_hvg.NT)
gif_seurat.NT <- seurat_cluster_somi_function(gif_seurat.NT, 20, resolution = 0.6)
DimPlot(gif_seurat.NT, label=TRUE, label.size = 8)
DimPlot(gif_seurat.NT, group.by = "batch")

gif_seurat.NT <- seurat_integ_function(gif_seurat.NT)
DimPlot(gif_seurat.NT, group.by = "batch")

gif_seurat.NT$orig.seurat_cluster <- gif_seurat.combined.2$seurat_clusters[colnames(gif_seurat.NT)]
DimPlot(gif_seurat.NT, group.by = "orig.seurat_cluster")

#save(gif_seurat.NT, file=paste0(rdatadir, "Seurat/gif_seurat.NT.RData"))

# comparison in injury
gif_seurat.D2

gif_sce.D2 <- as.SingleCellExperiment(gif_seurat.D2, assay = "RNA")
gif_sce.D2 <- norm_function(gif_sce.D2)
gif_hvg.D2 <- hvg_function(gif_sce.D2)

gif_seurat.D2 <- make_seurat_somi_function(gif_sce.D2, gif_hvg.D2)
gif_seurat.D2 <- seurat_cluster_somi_function(gif_seurat.D2, 20, resolution = 0.6)
DimPlot(gif_seurat.D2, label=TRUE, label.size = 8)
DimPlot(gif_seurat.D2, group.by = "batch")

gif_seurat.D2 <- seurat_integ_function(gif_seurat.D2)
DimPlot(gif_seurat.D2, group.by = "batch")

#save(gif_seurat.D2, file=paste0(rdatadir, "Seurat/gif_seurat.D2.RData"))

gif_seurat.D2$orig.seurat_cluster <- gif_seurat.combined.2$seurat_clusters[colnames(gif_seurat.D2)]
DimPlot(gif_seurat.D2, group.by = "orig.seurat_cluster")



