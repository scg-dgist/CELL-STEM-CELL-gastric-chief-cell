.libPaths(Sys.getenv("R_LIBS_STOMACH"))

load("D:/somi_function/ensemblGenes_mmusculus_2019-12-12.RData")
head(ensemblGenes)
GFP.df <- data.frame(row.names = "GFP", ensembl_gene_id="GFP", external_gene_name = "GFP", chromosome_name="ext", gene_biotype="GFP")
ensemblGenes = rbind(ensemblGenes, GFP.df)
tail(ensemblGenes)

rdatadir = "rdata3/"
plotdir = "plots3/"

library(Seurat)
library(ggplot2)
library(scran)

load(paste0(rdatadir, "sce/gif_sce.norm.RData"))
load(paste0(rdatadir, "sce/gif_hvg.norm.RData"))

#
gif_seurat.save <- gif_seurat # hvg 0.05, pc 30 k.param 30

#
make_seurat_somi_function <- function(sce, hvg){
  
  if(is.null(reducedDims(sce)) == FALSE){
    reducedDims(sce) = NULL
  }
  
  Seurat <- as.Seurat(sce,
                      counts = "counts",
                      data = "logcounts",
                      assay = "RNA")
  VariableFeatures(Seurat) <- rownames(hvg)
  
  all.genes <- rownames(Seurat)
  Seurat <- ScaleData(Seurat, features = all.genes, vars.to.regress = "pct_counts_MT")
  
  Seurat <- RunPCA(Seurat,
                   features = VariableFeatures(Seurat),
                   reduction.key = "pca_",
                   verbose = FALSE)
  print(plot(Seurat@reductions$pca@stdev,
             xlab = "PC",
             ylab = "Eigenvalue"))
  
  return(Seurat)
}
seurat_cluster_somi_function <- function(Seurat, PCA, resolution=0.8){
  
  Seurat <- FindNeighbors(Seurat, dims=1:PCA,
                          features = VariableFeatures(Seurat))
  Seurat <- FindClusters(Seurat,
                         resolution = resolution)
  
  Seurat <- RunTSNE(Seurat,
                    dims = 1:PCA,
                    features = VariableFeatures(Seurat),
                    check_duplicates = FALSE)
  Seurat <- RunUMAP(Seurat,
                    dims = 1:PCA)
  
  return(Seurat)
}

gif_seurat <- make_seurat_somi_function(gif_sce.norm, gif_hvg.norm)
gif_seurat <- seurat_cluster_somi_function(gif_seurat, 20)

DimPlot(gif_seurat, group.by = "seurat_clusters", label=TRUE, label.size=8)
DimPlot(gif_seurat, group.by = "batch", pt.size=1.5)
FeaturePlot(gif_seurat, "pct_counts_MT")

# cell cycle scoring
s.genes = cc.genes$s.genes
g2m.genes = cc.genes$g2m.genes

convertHumanGeneList <- function(x){
  
  require("biomaRt")
  human = useMart("ensembl", dataset = "hsapiens_gene_ensembl")
  mouse = useMart("ensembl", dataset = "mmusculus_gene_ensembl")
  
  genesV2 = getLDS(attributes = c("hgnc_symbol"), filters = "hgnc_symbol", values = x , mart = human, attributesL = c("mgi_symbol"), martL = mouse, uniqueRows=T)
  
  humanx <- unique(genesV2[, 2])
  
  # Print the first 6 genes found to the screen
  print(head(humanx))
  return(humanx)
}

s.genes <- convertHumanGeneList(s.genes)
g2m.genes <- convertHumanGeneList(g2m.genes)

s.geneids = ensemblGenes[ensemblGenes$external_gene_name %in% s.genes, ]$ensembl_gene_id
g2m.geneids = ensemblGenes[ensemblGenes$external_gene_name %in% g2m.genes, ]$ensembl_gene_id

DefaultAssay(gif_seurat) = "RNA"
gif_seurat <- CellCycleScoring(gif_seurat, s.geneids, g2m.geneids)

DimPlot(gif_seurat, group.by = "Phase", pt.size=2)
DimPlot(gif_seurat, group.by = "seurat_clusters", label=TRUE)

gif_seurat$CC.Difference = gif_seurat$S.Score - gif_seurat$G2M.Score

#
seurat_integ_function <- function(seurat, split.by="batch", PCA=10, resolution=0.5, min.size=100){
  
  DefaultAssay(seurat) = "RNA"
  seurat.list <- SplitObject(seurat, split.by)
  
  for(i in 1:length(seurat.list)){
    counts <- seurat.list[[i]]@assays$RNA@counts
    sce <- SingleCellExperiment(assays = list(counts = counts))
    if(ncol(sce) < 100){
      min.size=ncol(sce)
    }
    sce <- norm_function(sce, min.size=min.size)
    data <- sce@assays@data$logcounts
    seurat.list[[i]]@assays$RNA@data <- data
    print("hvg")
    
    hvg <- hvg_function(sce)
    print(nrow(hvg))
    VariableFeatures(seurat.list[[i]]) <- rownames(hvg)
  }
  
  features <- SelectIntegrationFeatures(object.list = seurat.list)
  
  anchors <- FindIntegrationAnchors(object.list = seurat.list, anchor.features = features)
  
  combined <- IntegrateData(anchorset = anchors)
  
  DefaultAssay(combined) <- "integrated"
  
  # Run the standard workflow for visualization and clustering
  combined <- ScaleData(combined, verbose = FALSE, vars.to.regress=c("pct_counts_MT", "CC.Difference"))
  #combined <- ScaleData(combined, verbose=FASE)
  combined <- RunPCA(combined, npcs = PCA, verbose = FALSE)
  combined <- RunUMAP(combined, reduction = "pca", dims = 1:PCA)
  combined <- FindNeighbors(combined, reduction = "pca", dims = 1:PCA, k.param=40)
  combined <- FindClusters(combined, resolution = resolution)
  
  return(combined)
}
norm_function <- function(sce, min.size=100){
  
  clusters <- quickCluster(sce, use.ranks = TRUE, min.size=min.size)
  sce <- computeSumFactors(sce, cluster=clusters)
  
  rm(clusters)
  
  sce <- normalize(sce)
  
  rownames(sce@assays@data$logcounts) = rownames(sce)
  colnames(sce@assays@data$logcounts) = colnames(sce)
  
  return(sce)
}
gif_seurat.combined <- seurat_integ_function(gif_seurat, PCA=30, resolution=1.0)

gif_seurat.combined <- FindClusters(gif_seurat.combined, resolution = 1.0)
DimPlot(gif_seurat.combined, pt.size = 2)
DimPlot(gif_seurat.combined, group.by = "Phase")

VlnPlot(gif_seurat.combined, "pct_counts_MT")

#
DefaultAssay(gif_seurat.combined) = "RNA"

DimPlot(gif_seurat.combined, pt.size=1.5, label=TRUE, label.size=8)
DimPlot(gif_seurat.combined, group.by = "batch", pt.size=1.5)

somi_featureplot <- function(seurat, feature, reduction="umap", col=NULL, size=1){

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
    
    FeaturePlot(seurat,
                features = feature_id,
                reduction = reduction,
                order=TRUE, pt.size=size) +
      scale_color_gradientn(colours = col, breaks=c(min,(max-min)/2 + min,max), 
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
    
  }
}
DefaultAssay(gif_seurat.combined) = "RNA"
somi_featureplot(gif_seurat.combined, "Cblif")
gif_seurat.combined

FeaturePlot(gif_seurat.combined, "pct_counts_MT")

markers <- FindAllMarkers(gif_seurat.combined, logfc.threshold = 0)
markers$symbol = ensemblGenes[markers$gene, "external_gene_name"]
head(subset(markers, cluster==3))

write.csv(markers, file=paste0(rdatadir, "Seurat/cluster_markers.csv"))

DimPlot(gif_seurat.combined, label=TRUE, label.size=8)

somi_featureplot(gif_seurat.combined, "Amy2a3")

clustern = length(unique(gif_seurat.combined$seurat_clusters))

save(gif_seurat, file=paste0(rdatadir, "Seurat/gif_seurat.RData"))
save(gif_seurat.combined, file=paste0(rdatadir, "Seurat/gif_seurat.combined.RData"))

# color
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

somi_dimplot2 <- function(seurat, reduction="umap", group.by, size=1, cols=NULL, cells=NULL){
  DimPlot(seurat, cells=cells,
          reduction = reduction,
          group.by = group.by,
          pt.size = size,
          label = FALSE,
          label.size = 6) +
    labs(title = paste0("by ", group.by)) +
    scale_colour_manual(values= cols) +
    theme(plot.title = element_text(size=40, face="italic", hjust=0.5, margin=margin(0,0,15,0,"mm")),
          axis.line=element_blank(),
          axis.ticks=element_blank(),
          axis.title = element_blank(),
          axis.text = element_blank(),
          legend.text = element_text(size=20, face="italic")
    )
}
List = c("seurat_clusters", "batch", "condition", "Phase")
for(i in 1:length(List)){
  d <- somi_dimplot2(gif_seurat.combined, group.by = List[i], size=3, cols = colors[[i]])
  ggsave(paste0(plotdir, "Seurat/dimplot/umap_", List[i], ".png"), plot=d, height=8, width=10)
}

FeaturePlot(gif_seurat.combined, "log10_sum", pt.size = 1.5)
FeaturePlot(gif_seurat.combined, "log10_detected", pt.size = 1.5)
FeaturePlot(gif_seurat.combined, "pct_counts_ERCC", pt.size = 1.5)

somi_featureplot <- function(seurat, feature, reduction="umap", col=NULL, size=1){
  
  feature_id = rownames(ensemblGenes[ensemblGenes$external_gene_name == feature,])
  print(c(feature, feature_id))
  
  if(feature == "GFP"){
    feature_id = feature
  }
  
  if(is.null(col)){
    col = c("#EEEEEE", brewer.pal(9, "Reds"))
  }
  if(sum(feature_id %in% rownames(seurat)) == 0){
    print(paste0("no ", feature))
    
  }else{
    expr = GetAssayData(seurat, slot="data", assay="RNA")
    max = max(expr[feature_id,])
    min = min(expr[feature_id,])
    
    FeaturePlot(seurat,
                features = feature_id,
                reduction = reduction,
                order=TRUE, pt.size=size) +
      scale_color_gradientn(colours = col, breaks=c(min,(max-min)/2 + min,max), 
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
    
  }
}
library(inlmisc)
col = GetColors(256, start=0.3, end=1)
col=c("lightgray", GetColors(256, start=0.3, end=1))

genes = c("Cblif", "GFP", "Gkn3", "Muc6", "Muc5ac", "Cd44", "Stmn1", "Mki67", "Lgr5", "Tnfrsf19", "Cdkn1c", "Tff2", "Atp4a", "Atp4b", "Kcne2", "Clu", "Sox9", "Amy2a3", "Amy2b")
genes = c("Cldn18", "Anxa10")
genes = c("Pnliprp2", "Clps", "Furin")
DefaultAssay(gif_seurat.combined) = "RNA"
for(i in genes){
  f <- somi_featureplot(gif_seurat.combined, i, reduction="umap", col=col, size=3)
  ggsave(paste0(plotdir, "Seurat/featureplot/umap_expr_", i, ".png"), plot=f, height=8, width=10)
}

#
NT.cells = colnames(gif_seurat.combined[, gif_seurat.combined$condition == "NT"])
d <- somi_dimplot2(gif_seurat.combined, group.by = "seurat_clusters", cols = colors[[1]], cells=NT.cells, size=3)
ggsave(paste0(plotdir, "Seurat/dimplot/umap_seurat_clusters_NT.png"), plot=d, height=8, width=10)

D2.cells = colnames(gif_seurat.combined[, gif_seurat.combined$condition == "D2"])
d <- somi_dimplot2(gif_seurat.combined, group.by = "seurat_clusters", cols = colors[[1]], cells=D2.cells, size=3)
ggsave(paste0(plotdir, "Seurat/dimplot/umap_seurat_clusters_D2.png"), plot=d, height=8, width=10)

# proportion bar plot
Idents(gif_seurat.combined) = gif_seurat.combined$condition
gif_seurat.NT = subset(gif_seurat.combined, idents = "NT")
gif_seurat.D2 = subset(gif_seurat.combined, idents = "D2")

cl_order = c(0,1,5,3,4,2) %>% as.character

df = data.frame(pct = c((table(gif_seurat.D2$seurat_clusters) / ncol(gif_seurat.D2)) %>% as.character(),
                        (table(gif_seurat.NT$seurat_clusters) / ncol(gif_seurat.NT)) %>% as.character()) %>% as.numeric(),
                cluster = rep(c(0:(clustern-1)), 2) %>% as.character(),
                condition = rep(c("D2", "NT"), each=clustern))
df$condition = factor(df$condition, levels = c("NT", "D2"))

gg_color_hue <- function(n) {
  hues = seq(15, 375, length = n + 1)
  hcl(h = hues, l = 65, c = 100)[1:n]
}
condition_col = rev(c(as.character(lacroix_palette("CranRaspberry"))[1], as.character(lacroix_palette("Pure"))[1]))
g <- ggplot(df, aes(x=cluster, y=pct, fill=condition)) +
  geom_bar(stat = "identity", position="dodge") +
  scale_y_continuous(name = "percentage (%)", 
                     #breaks = c(0.05 * c(0:18)),
                     #labels = c(5 * c(0:18)),
                     #limits = c(0, 0.225),
                     expand=c(0,0)) +
  scale_x_discrete(breaks = as.character(cl_order), limits = as.character(cl_order)) +
  scale_fill_manual(values = condition_col) +
  #geom_text(aes(label = paste(round(pct * 100, digit=2), "%", sep="")), position = position_dodge(0.9),
  #          vjust = -0.5, color = "black", size=4, fontface="italic") +
  labs(x="cluster") +
  theme_classic() + 
  theme(plot.title = element_text(size=20, face="italic", hjust = 0.5, margin=margin(1,0,5,0,"mm")),
        axis.title.y = element_text(size=15, margin=margin(0,5,0,0,"mm")),
        axis.title.x = element_text(size=15, margin=margin(3,0,0,0,"mm")),
        axis.text = element_text(size=15),
        panel.grid.major.y = element_line(size=0.5))
g
ggsave(paste0(plotdir, "Seurat/barplot_proportion_seurat_clusters.png"), plot=g, width=6, height=4)

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

gif_seurat.NT$orig.seurat_cluster <- gif_seurat.combined$seurat_clusters[colnames(gif_seurat.NT)]
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

gif_seurat.D2$orig.seurat_cluster <- gif_seurat.combined$seurat_clusters[colnames(gif_seurat.D2)]
DimPlot(gif_seurat.D2, group.by = "orig.seurat_cluster")

#
somi_dotplot <- function(seurat, features, cluster_order=NULL, ident=NULL){
  
  features = unlist(features)
  df = ensemblGenes[ensemblGenes$external_gene_name %in% features,]
  df = df[match(features, df$external_gene_name),]
  geneids = df$ensembl_gene_id
  
  cluster_order = as.character(cluster_order)
  
  DotPlot(seurat, assay = "RNA", group.by = ident, features = geneids, col.min = -1.5, col.max = 1.5) + ylab("cluster") + xlab("") +
    scale_y_discrete(limits = cluster_order) +
    scale_x_discrete(limits = geneids, labels=features) +
    scale_color_gradientn(colors = c("#2166AC", "#F7F7F7", "#B2182B"), limits = c(-1.5,  1.5), breaks = c(-1.5, 0, 1.5)) +
    scale_radius(range = c(1, 8)) +
    theme(plot.title = element_text(size=50, face="bold.italic"),
          panel.grid.major=element_blank(),
          panel.grid.minor=element_blank(), 
          #axis.line=element_blank(),
          #axis.ticks=element_blank(),
          #axis.title = element_blank(),
          panel.border = element_blank(),
          #legend.position = "none",
          legend.text=element_text(size=18), 
          legend.title=element_text(size=18),
          #legend.key=element_blank(),
          axis.text.x = element_text(size=20, angle = 0, hjust=0.5, vjust=1),
          axis.text.y = element_text(size=20),
          axis.title.y = element_blank()
    ) + coord_flip()
}

genes = c("Cblif", "Gper1", "Bhlha15", "Pga5", "Cdkn1c", "Muc6", "Gkn3", "Mki67", "Stmn1",
          "Tff2", "Muc5ac", "Tff1", "Gkn1", "Gkn2")
DimPlot(gif_seurat.combined, group.by = "seurat_clusters", label=TRUE, label.size=8)


s <- somi_dotplot(gif_seurat.combined, rev(genes), cl_order, "seurat_clusters")
ggsave(paste0(plotdir, "Seurat/gif_dotplot.png"), plot=s, width=6, height=5)
s <- somi_dotplot(gif_seurat.NT, rev(genes), cl_order, "orig.seurat_cluster")
ggsave(paste0(plotdir, "Seurat/gif_NT_dotplot.png"), plot=s, width=6, height=5)
s <- somi_dotplot(gif_seurat.D2, rev(genes), cl_order, "orig.seurat_cluster")
ggsave(paste0(plotdir, "Seurat/gif_D2_dotplot.png"), plot=s, width=6, height=5)


#
somi_featureplot(gif_seurat.D2, "Atp4b")

List = c("seurat_clusters", "batch", "condition")
for(i in 1:3){
  d <- somi_dimplot2(gif_seurat.D2, group.by = List[i], size=3, cols = colors[[i]])
  ggsave(paste0(plotdir, "Seurat.D2/dimplot/umap_", List[i], ".png"), plot=d, height=8, width=10)
}

markers.D2 <- FindAllMarkers(gif_seurat.D2, logfc.threshold = 0)
markers.D2$symbol = ensemblGenes[markers.D2$gene, "external_gene_name"]
head(subset(markers.D2, cluster == 0), n=30)

DimPlot(gif_seurat.D2, label=TRUE)

somi_featureplot(gif_seurat.D2, "Atp4b")

DimPlot(gif_seurat.D2, group.by = "seurat_clusters", label=TRUE, label.size=8)
markers.23 <- FindMarkers(gif_seurat.D2, ident.1 = c(2,3), logfc.threshold = 0)
markers.23$symbol = ensemblGenes[rownames(markers.23), "external_gene_name"]

head(subset(markers.23, avg_logFC > 0), n=30)

write.csv(markers.23, file=paste0(rdatadir, "Seurat.D2/markers.23.csv"))

up.genes.df <- subset(markers.23, avg_logFC > 0.1 & p_val_adj < 0.05 & !is.na(symbol))
head(up.genes.df, n=20)
down.genes.df <- subset(markers.23, avg_logFC < -0.1 & p_val_adj < 0.05 & !is.na(symbol))
head(down.genes.df)

library(topGO)
library(dplyr)
library(org.Mm.eg.db)
library(plyr)

BGgenes = rownames(gif_seurat.D2)
topGO_function <- function(sublist, onts=c("BP", "MF", "CC"), topnodes=1000){
  
  targetGenes = rownames(sublist)
  
  backgroundGenes = BGgenes
  
  allGene = factor(as.integer(backgroundGenes %in% targetGenes))
  names(allGene) = backgroundGenes
  
  tab = as.list(onts)
  names(tab) = onts
  
  for(i in 1:length(onts)){
    tgd = new("topGOdata", ontology=onts[i], allGenes=allGene, nodeSize=5,
              annot=annFUN.org, mapping="org.Mm.eg.db", ID="Ensembl")
    resultTopGO.elim = runTest(tgd, algorithm="elim", statistic="Fisher")
    tab[[i]] = GenTable(tgd, Fisher.elim=resultTopGO.elim,
                        orderBy="Fisher.classic", topNodes=topnodes)
    tab[[i]]$onts = onts[i]
    
    genes_in_term = genesInTerm(tgd, tab[[i]]$GO.ID)
    sigGenesList = lapply(genes_in_term, function(x){ intersect(x, targetGenes) })
    sigSymsList = lapply(sigGenesList, function(x){ ensemblGenes[x, "external_gene_name"] })  
    
    tab[[i]]$sig_geneids = sapply(sigSymsList, function(x){paste(x, collapse=" ")})
    tab[[i]]$sig_symbols = sapply(sigSymsList, function(x){paste(x, collapse=" ")})
    
    print(head(tab[[i]]))
  }
  
  topGOResults = rbind.fill(tab)
  topGOResults$Fisher.elim = as.numeric(topGOResults$Fisher.elim)
  if(sum(is.na(topGOResults$Fisher.elim)) > 0){
    topGOResults[is.na(topGOResults$Fisher.elim),]$Fisher.elim = 0
  }
  topGOResults = topGOResults[order(topGOResults$Fisher.elim),]
  #topGOResults = subset(topGOResults, Fisher.elim < 0.05)
  topGOResults[["log10P"]] = -log10(topGOResults[["Fisher.elim"]])
  
  return(topGOResults)
  
}

topgoresult.up <- topGO_function(up.genes.df, onts="BP")
head(subset(topgoresult.up, onts=="BP"))

topgoresult.down <- topGO_function(down.genes.df, onts="BP")
head(subset(topgoresult.down, onts=="BP"))

write.csv(topgoresult.up, file=paste0(rdatadir, "Seurat.D2/Cl23.goUP.csv"))
write.csv(topgoresult.down, file=paste0(rdatadir, "Seurat.D2/Cl23.goDOWN.csv"))
