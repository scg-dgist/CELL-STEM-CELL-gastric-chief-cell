##### Code for Mouse Stomach Pgc+ Single-cell RNA-seq data
##### Author: Somi Kim (ksomi47@dgist.ac.kr)
##### Last Update: 2021/02/17

### Pgc+ scRNA-seq data
### QC



library(SingleCellExperiment)
library(scater)
library(scran)



##### load data
countData = read.csv("/stomach_injury_count_table.csv", row.names=1, header=TRUE, check.names=FALSE)

# remove low quality sample
countData = countData[,!grepl("24146.*", colnames(countData))]
countData = countData[rowSums(countData) > 0, ] 

### create SingleCellExperiment (SCE) object
cellInfo <- data.frame(cell = colnames(countData),
                       sample = sapply(colnames(countData), function(x) substr(x, 1, 5)),
                       condition = c(rep("Injury", 384), rep("Control", 384), rep("Injury", 384)))
cd <- as.matrix(countData)
logcd <- log2(cd + 1)
stomach_sce <- SingleCellExperiment(assays = list(counts = cd, logcounts = logcd), colData = cellInfo)
rm(cd)
rm(logcd)
rm(cellInfo)
rm(countData)

# remove control wells
control_ind = grepl("^.*_1$|^.*_381$|^.*_382$|^.*_383$|^.*_384$", colnames(stomach_sce))
stomach_sce = stomach_sce[, !control_ind]
rm(control_ind)

keep_feature <- rowSums(counts(stomach_sce) >0)>0

stomach_sce = stomach_sce[keep_feature,]

rm(keep_feature)

mtGenes = ensemblGenes[ensemblGenes[,3] == "MT",]
is.mito = rownames(stomach_sce) %in% mtGenes[,1]

rm(mtGenes)

stomach_sce <- calculateQCMetrics(stomach_sce,
                                  feature_controls = list(MT = is.mito))
rm(is.mito)

stomach_sce = runPCA(stomach_sce)
plotPCA(stomach_sce,
        size_by = "total_features_by_counts",
        colour_by = "total_counts")


### Set filters to exclude poor quality cells
## total counts filter
hist(stomach_sce$total_counts, breaks=100, xlab = "Total Counts", main="")
abline(v=300000, col="red")
filter_by_total_counts = stomach_sce$total_counts > 300000

## total features by counts
hist(stomach_sce$total_features_by_counts, breaks=100, xlab = "Total Feature Counts", main="")
abline(v=4000, col="red")
filter_by_total_features = stomach_sce$total_features_by_counts > 4000

## mt percents
hist(stomach_sce$pct_counts_MT, breaks=100, xlab = "Mitochondrial %", main="")
abline(v=10, col="red")
filter_by_MT = stomach_sce$pct_counts_MT < 10


### Quality Control
stomach_sce$use <- (filter_by_total_counts &
                      filter_by_total_features &
                      filter_by_MT)
rm(filter_by_MT)
rm(filter_by_total_counts)
rm(filter_by_total_features)

plotPCA(stomach_sce,
        size_by="total_features_by_counts",
        colour_by="use")

stomach_sce.qc <- stomach_sce[, colData(stomach_sce)$use]


### Normalization
keep <- rowSums(counts(stomach_sce.qc) > 0) > 0
stomach_sce.qc = stomach_sce.qc[keep,] # 29905 genes
rm(keep)

clusters <- quickCluster(stomach_sce.qc, use.ranks=TRUE)
stomach_sce.qc <- computeSumFactors(stomach_sce.qc, cluster=clusters)
rm(clusters)

plot(sizeFactors(stomach_sce.qc), stomach_sce.qc$total_counts/1e3, log="xy",
     ylab="Library size (kilos)", xlab="Size factor")

stomach_sce.norm <- normalize(stomach_sce.qc)

rownames(stomach_sce.norm@assays$data$logcounts) = rownames(stomach_sce.norm)
colnames(stomach_sce.norm@assays$data$logcounts) = colnames(stomach_sce.norm)


### Identify highly variable genes
var.fit <- trendVar(stomach_sce.norm, method="spline", parametric=TRUE, use.spikes=FALSE)
var.out <- decomposeVar(stomach_sce.norm, var.fit)
stomach_hvg <- var.out[which(var.out$FDR < 0.01 & var.out$bio>0.2),]

plot(x=var.out$mean, y=var.out$total, pch=16, cex=0.3,
     ylab="Variance of log-expression", xlab="Mean log-expression",
     ylim=range(0, 20), main = "initial")
o <- order(var.out$mean)
lines(x=var.out$mean[o], y=var.out$tech[o], col="dodgerblue", lwd=2)
points(y = var.out$total[var.out$FDR < 0.01 & var.out$bio > .2],
       x = var.out$mean[var.out$FDR < 0.01 & var.out$bio > .2],
       pch=16, cex=0.3, col="red")

rm(var.fit)
rm(var.out)
rm(o)
