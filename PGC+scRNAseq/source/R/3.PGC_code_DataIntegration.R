##### Code for Mouse Stomach Pgc+ Single-cell RNA-seq data
##### Author: Somi Kim (ksomi47@dgist.ac.kr)
##### Last Update: 2021/02/17

### Pgc+ scRNA-seq data
### 3. Data Integration


library(scran)
library(Seurat)
library(plyr)


##### load data
load("chief_sce.norm.RData")
load("stomach_sce.norm.RData")
load("stomach_hvg.RData")
load("chief_hvg.RData")
load("ensemblGenes_mmusculus_2019-12-12.RData")

### Combine two data

chief_counts = chief_sce.norm@assays@data$counts
stomach_counts = stomach_sce.norm@assays$data$counts

pgc_counts = rbind.fill(as.data.frame(t(chief_counts)), as.data.frame(t(stomach_counts)))
pgc_counts = t(pgc_counts)

# select metadata to keep
keep_coldata = c("cell", "sample", "condition", "total_features_by_counts",
                 "log10_total_features_by_counts", "total_counts", "log10_total_counts",
                 "pct_counts_MT")
chief_coldata = colData(chief_sce.norm)[,keep_coldata]
stomach_coldata = colData(stomach_sce.norm)[,keep_coldata]

pgc_coldata = rbind(chief_coldata, stomach_coldata)
pgc_coldata$data = c(rep("data1", nrow(chief_coldata)), rep("data2", nrow(stomach_coldata)))
pgc_coldata = pgc_coldata[,c("cell", "sample", "condition", "data", "total_features_by_counts",
                             "log10_total_features_by_counts", "total_counts", "log10_total_counts",
                             "pct_counts_MT")]

rm(keep_coldata)
rm(chief_coldata)
rm(stomach_coldata)


### Create SingleCellExperiment Object
pgc_sce <- SingleCellExperiment(assays = list(counts = pgc_counts),
                                colData = pgc_coldata)

pgc_sce$sample = as.factor(pgc_sce$sample)
pgc_sce$condition = as.factor(pgc_sce$condition)
pgc_sce$data = as.factor(pgc_sce$data)


### Integrate two normalized data for further integration

chief_normTable = chief_sce.norm@assays@data$logcounts
stomach_normTable = stomach_sce.norm@assays$data$logcounts

normTable = rbind.fill(as.data.frame(t(chief_normTable)), as.data.frame(t(stomach_normTable)))
normTable = t(normTable)
normTable[is.na(normTable)] = 0
colnames(normTable) = c(colnames(chief_normTable), colnames(stomach_normTable))

chief_normTable = normTable[,colnames(chief_normTable)]
stomach_normTable = normTable[,colnames(stomach_normTable)]

normTableList = list(chief_normTable, stomach_normTable)

# hvgenes list
hvgList = list(rownames(chief_hvg), rownames(stomach_hvg))


### seurat integration
pgc_seurat = as.Seurat(pgc_sce,
                       counts = "counts",
                       data = "logcounts",
                       assay = "RNA")
all.genes = rownames(pgc_seurat)

seurat.list <- SplitObject(seurat, split.by = "data")

for(i in 1:length(seurat.list)){
  seurat.list[[i]]@assays$RNA@data <- normTableList[[i]]
  seurat.list[[i]]@assays$RNA@var.features = hvgList[[i]]
}

reference.list <- seurat.list[names(seurat.list)]
anchors <- FindIntegrationAnchors(object.list = reference.list)

pgc_seurat <- IntegrateData(anchorset = anchors, 
                            features.to.integrate = all.genes)

DefaultAssay(pgc_seurat) = "integrated"

pgc_seurat <- ScaleData(pgc_seurat, features = all.genes)
pgc_seurat <- RunPCA(pgc_seurat,
                     features = VariableFeatures(pgc_seurat),
                     reduction.key = "pca_",
                     verbose=FALSE)

pgc_seurat$data_condition <- as.factor(paste(pgc_seurat$data, pgc_seurat$condition, sep="_"))


### CellCycleScoring

s.genes = paste0(substr(cc.genes$s.genes, 1, 1),
                 tolower(substring(cc.genes$s.genes, 2)))
g2m.genes = paste0(substr(cc.genes$g2m.genes, 1, 1),
                   tolower(substring(cc.genes$g2m.genes, 2)))

s.genes.id = ensemblGenes[ensemblGenes$external_gene_name %in% s.genes,]$ensembl_gene_id
g2m.genes.id = ensemblGenes[ensemblGenes$external_gene_name %in% g2m.genes,]$ensembl_gene_id

pgc_seurat <- CellCycleScoring(pgc_seurat,
                           s.features = s.genes.id,
                           g2m.features = g2m.genes.id)
pgc_seurat$CC.Difference <- pgc_seurat$S.Score - pgc_seurat$G2M.Score

pgc_seurat <- ScaleData(pgc_seurat, 
                        vars.to.regress = "CC.Difference", 
                        features = rownames(pgc_seurat))


### Define highly variable genes after data integration
sce = as.SingleCellExperiment(pgc_seurat,
                              assay = "RNA")

var.fit <- trendVar(sce, method="spline", parametric=TRUE, use.spikes=FALSE)
var.out <- decomposeVar(sce, var.fit)

fdr_cutoff = 0.01
bio_cutoff = 0.5

hvg <- var.out[which(var.out$FDR < fdr_cutoff & var.out$bio > bio_cutoff),]

plot(x=var.out$mean, y=var.out$total, pch=16, cex=0.3,
     ylab="Variance of log-expression", xlab="Mean log-expression",
     ylim=range(0, 20), main = "initial")
o <- order(var.out$mean)
lines(x=var.out$mean[o], y=var.out$tech[o], col="dodgerblue", lwd=2)
points(y = var.out$total[var.out$FDR < fdr_cutoff & var.out$bio > bio_cutoff],
       x = var.out$mean[var.out$FDR < fdr_cutoff & var.out$bio > bio_cutoff],
       pch=16, cex=0.3, col="red")

VariableFeature(pgc_seurat) <- rownames(hvg)



pgc_seurat <- RunPCA(pgc_seurat,
                     assay = "integrated",
                     npcs=300,
                     features = VariableFeatures(pgc_seurat),
                     verbose=FALSE)

PCA=20
pgc_seurat <- RunTSNE(pgc_seurat,
                      dims = 1:PCA,
                      features = VariableFeatures(pgc_seurat))
pgc_seurat <- RunUMAP(pgc_seurat,
                      dims = 1:PCA)

pgc_seurat <- FindNeighbors(pgc_seurat, 
                            dims=1:PCA, 
                            reduction="pca", 
                            assay="integrated")
pgc_seurat <- FindClusters(pgc_seurat, resolution=1.5)
