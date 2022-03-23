##### Code for Mouse Stomach Pgc+ Single-cell RNA-seq data
##### Author: Somi Kim (ksomi47@dgist.ac.kr)
##### Last Update: 2021/02/17

### Pgc+ scRNA-seq data - published
### QC



library(SingleCellExperiment)
library(scater)
library(scran)



##### load data
chiefCd = read.csv("stomach_chief_count_table.csv", row.names=1, header=TRUE, check.names=FALSE)

keep = rowSums(chiefCd) > 0
chiefCd = chiefCd[keep, ]

### Create SingleCellExperiment(SCE) object
cellInfo <- data.frame(cell = colnames(chiefCd),
                       sample = sapply(colnames(chiefCd), function(x) substr(x, 1, 5)),
                       condition = c(rep("Control", ncol(chiefCd))))

cd <- as.matrix(chiefCd)
logcd <- log2(cd+1)
chief_sce <- SingleCellExperiment(assays = list(counts = cd, logcounts = logcd),
                                  colData = cellInfo)
rm(cd)
rm(logcd)
rm(cellInfo)
rm(chiefCd)

# remove control wells
control_ind = grepl("(^.*_1$|^.*_381$|^.*_382$|^.*_383$|^.*_384$)",colnames(chief_sce))
chief_sce = chief_sce[, !control_ind]

### Set filters to exclude poor quality cells
keep_feature <- rowSums(counts(chief_sce) >0 ) >0
chief_sce = chief_sce[keep_feature,]

rm(keep_feature)

load("E:/somi_function/ensemblGenes_mmusculus_2019-12-12.RData")
mtGenes = ensemblGenes[ensemblGenes[,3] == "MT",]
is.mito = rownames(chief_sce) %in% mtGenes[,1]

rm(mtGenes)

chief_sce <- calculateQCMetrics(chief_sce,
                                feature_controls = list(MT = is.mito),
                                detection_limit = 5)

rm(is.mito)

chief_sce = runPCA(chief_sce)
plotPCA(chief_sce,
        size_by = "total_features_by_counts",
        colour_by = "total_counts")

hist(chief_sce$pct_counts_MT, breaks=100)
abline(v=10, col="red")
filter_by_MT = chief_sce$pct_counts_MT < 10
table(filter_by_MT)

filter_by_feature_counts = chief_sce$total_features_by_counts / nrow(chief_sce) * 100 > 3

chief_sce$use <- (filter_by_MT & 
                    filter_by_feature_counts)
table(chief_sce$use) # 743 cells use

rm(filter_by_MT)
rm(filter_by_feature_counts)

plotPCA(chief_sce,
        size_by="total_features_by_counts",
        colour_by="use")

chief_sce.qc <- chief_sce[, colData(chief_sce)$use]
dim(chief_sce.qc) # 27820 743

rm(chief_sce)

chief_sce.qc = runPCA(chief_sce.qc)

save(chief_sce.qc, file="data/Sce/chief_sce.qc.RData")


### Normalization
keep <- rowSums(counts(chief_sce.qc) > 0 ) > 0
table(keep)
chief_sce.qc = chief_sce.qc[keep,]
chief_sce.qc # 26855 743

rm(keep)

clusters <- quickCluster(chief_sce.qc, use.ranks = TRUE)
chief_sce.qc <- computeSumFactors(chief_sce.qc, cluster=clusters)

rm(clusters)

summary(sizeFactors(chief_sce.qc))

plot(sizeFactors(chief_sce.qc), chief_sce.qc$total_counts/1e3, log="xy",
     ylab="Library size (kilos)", xlab="Size factor")

chief_sce.norm <- normalize(chief_sce.qc)

rownames(chief_sce.norm@assays@data$logcounts) = rownames(chief_sce.norm)
colnames(chief_sce.norm@assays@data$logcounts) = colnames(chief_sce.norm)

save(chief_sce.norm, file="data/chief_data/chief_sce.norm.RData")

rm(chief_sce.qc)


### Identify highly variable genes
var.fit <- trendVar(chief_sce.norm, method="spline", parametric=TRUE, use.spikes=FALSE)
var.out <- decomposeVar(chief_sce.norm, var.fit)
chief_hvg <- var.out[which(var.out$FDR < 0.01 & var.out$bio > 0.5),]

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

save(chief_hvg, file="data/chief_data/chief_hvg.RData")
