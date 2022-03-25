load(paste0(rdatadir, "Seurat/gif_seurat.combined.3.RData"))

library(slingshot)

gif.cells = colnames(gif_seurat.combined.3)
df.gif = df.gif_to_pgc[gif.cells,]

embeddings <- data.frame(x=df.gif$x,
                         y=df.gif$y)
clusters <- gif_seurat.combined.3$seurat_clusters

slg <- slingshot(embeddings, clusters, start.clus=3)

plot(reducedDim(slg), pch = 16, cex = 0.5)
lines(slg, lwd = 2, col = 'yellow')

pt <- slingPseudotime(slg)
pt
pt = pt / max(pt)

gif_seurat.combined.3$Pseudotime = pt

curve.embeddings = slg@curves$curve1$s

df.gif_to_pgc$curve.x = NA
df.gif_to_pgc$curve.y = NA
df.gif_to_pgc[gif.cells,]$curve.x = curve.embeddings[,1]
df.gif_to_pgc[gif.cells,]$curve.y = curve.embeddings[,2]

df.gif_to_pgc$pseudotime = NA
df.gif_to_pgc[colnames(mat.gif), ]$pseudotime = gif_seurat.combined.3$Pseudotime
df.gif_to_pgc$pseudotime = as.numeric(df.gif_to_pgc$pseudotime)

library(inlmisc)
col=c("lightgray", GetColors(256, start=0.3, end=1))

g0 <- ggplot(df.gif_to_pgc[order(df.gif_to_pgc$pseudotime, na.last = FALSE),]) +
  geom_point(aes(x=x, y=y, color=pseudotime, size=celltype))+
  #guides(colour = guide_legend(override.aes = list(size=10)))+
  scale_size_manual(values = c(rep(2.5, clustern),2)) + 
  #scale_color_gradientn(colors=col, na.value = "grey50") +
  scale_color_viridis_c() +
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
g0 + geom_point(aes(x=curve.x, y=curve.y), color="black") + guides(size="none")
ggsave(paste0(plotdir, "Seurat.3/pseudotime_projection_trajectory.pdf"), plot=g0, height=6, width=7.5)

df.gif_to_pgc$celltype = factor(df.gif_to_pgc$celltype, levels = c("Gif_low", "Gif_high", "ref"))
color.gif = c(as.character(lacroix_palette("Orange"))[1], as.character(lacroix_palette("Pure"))[1], "lightgray")
gp <- ggplot(df.gif_to_pgc[rev(order(df.gif_to_pgc$celltype)),]) +
  geom_point(aes(x=x, y=y, color=celltype, size=celltype))+
  #guides(colour = guide_legend(override.aes = list(size=10)))+
  scale_size_manual(values = c(rep(2.5, clustern),2)) + 
  scale_color_manual(values = color.gif) +
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
ggsave(paste0(plotdir, "Seurat.3/projection_gif_hl.png"), plot=gp, width=6, height=5)

df.gif_to_pgc$celltype_cond = paste(df.gif_to_pgc$celltype, df.gif_to_pgc$cond, sep=".")
df.gif_to_pgc$celltype_cond = gsub("ref.", "ref", df.gif_to_pgc$celltype_cond)
df.gif_to_pgc$celltype_cond = gsub("NT", "control", df.gif_to_pgc$celltype_cond)
df.gif_to_pgc$celltype_cond = gsub("D2", "injury", df.gif_to_pgc$celltype_cond)
df.gif_to_pgc$celltype_cond = factor(df.gif_to_pgc$celltype_cond, levels = c("Gif_low.injury", "Gif_high.injury", "Gif_low.control", "Gif_high.control", "ref"))
color.1 = c(lacroix_palette("PeachPear", 6)[1],
            lacroix_palette("Orange", 6)[2],
            lacroix_palette("Pure", 6)[1],
            lacroix_palette("Pure", 6)[4],
            "lightgray")
gc <- ggplot(df.gif_to_pgc[rev(order(df.gif_to_pgc$celltype_cond)),]) +
  geom_point(aes(x=x, y=y, color=celltype_cond , size=celltype))+
  guides(colour = guide_legend(override.aes = list(size=10)))+
  scale_size_manual(values = c(rep(2.5, clustern),2)) + 
  scale_color_manual(values = color.1) +
  guides(size="none") +
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
ggsave(paste0(plotdir, "Seurat.3/projection_4groups.png"), plot=gc, width=8.5, height=6)

df.gif_to_pgc$celltype_cond2 = df.gif_to_pgc$celltype_cond
df.gif_to_pgc$celltype_cond2 = gsub("ref.", "ref", df.gif_to_pgc$celltype_cond2)
df.gif_to_pgc$celltype_cond2 = gsub("Gif_high.NT", "control", df.gif_to_pgc$celltype_cond2)
df.gif_to_pgc$celltype_cond2 = gsub("Gif_low.NT", "control", df.gif_to_pgc$celltype_cond2)
df.gif_to_pgc$celltype_cond2 = gsub("Gif_low.D2", "Gif_low.injury", df.gif_to_pgc$celltype_cond2)
df.gif_to_pgc$celltype_cond2 = gsub("Gif_high.D2", "Gif_high.injury", df.gif_to_pgc$celltype_cond2)
df.gif_to_pgc$celltype_cond2 = factor(df.gif_to_pgc$celltype_cond2, levels = c("Gif_low.injury", "Gif_high.injury", "control", "ref"))
color.2 = c(lacroix_palette("PeachPear", 6)[1],
            lacroix_palette("Orange", 6)[2],
            lacroix_palette("Pure", 6)[1],
            "lightgray")
ga <- ggplot(df.gif_to_pgc[rev(order(df.gif_to_pgc$celltype_cond2)),]) +
  geom_point(aes(x=x, y=y, color=celltype_cond2, size=celltype))+
  guides(colour = guide_legend(override.aes = list(size=10)))+
  scale_size_manual(values = c(rep(2.5, clustern),2)) + 
  scale_color_manual(values = color.2) +
  guides(size="none") +
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
ggsave(paste0(), plot=ga, width=8, height=6)

df.2 <- df.gif_to_pgc[!df.gif_to_pgc$celltype_cond %in% c("Gif_low.control"),]
df.2$celltype_cond = factor(df.2$celltype_cond, levels = c("Gif_low.injury", "Gif_high.injury", "Gif_high.control", "ref"))
color.3 = c(lacroix_palette("PeachPear", 6)[1],
            lacroix_palette("Orange", 6)[2],
            lacroix_palette("Pure", 6)[4],
            "lightgray")

gb <- ggplot(df.2[rev(order(df.2$celltype_cond)),]) +
  geom_point(aes(x=x, y=y, color=celltype_cond , size=celltype))+
  guides(colour = guide_legend(override.aes = list(size=10)))+
  scale_size_manual(values = c(rep(2.5, clustern),2)) + 
  scale_color_manual(values = color.3) +
  guides(size="none") +
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
ggsave(paste0(plotdir, "Seurat.3/projection_3groups.png"), plot=gb, width=8.5, height=6)

# draw line graph

df.gif <- df.gif_to_pgc[gif.cells,]
#df.gif$pseudotime = df.gif$pseudotime / max(df.gif$pseudotime)

g1 <- ggplot(df.gif, aes(x=pseudotime, y=gif, color=pseudotime)) + 
  geom_point() + 
  scale_color_viridis_c() +
  expand_limits(x = 0, y = 0) +
  geom_smooth(method = "loess", col="black", lwd=1.5) +
  theme_classic()
ggsave(paste0(plotdir, "Seurat.3/pseudotime_gif.pdf"), plot=g1, height=6, width=7.5)

g2 <- ggplot(df.gif, aes(x=pseudotime, y=ki67, color=pseudotime)) + 
  geom_point() + 
  scale_color_viridis_c() +
  expand_limits(x = 0, y = 0) +
  geom_smooth(method = "loess", col="black", lwd=1.5) +
  theme_classic()
ggsave(paste0(plotdir, "Seurat.3/pseudotime_ki67.pdf"), plot=g2, height=6, width=7.5)
g3 <- ggplot(df.gif, aes(x=pseudotime, y=p57, color=pseudotime)) + 
  geom_point() + 
  scale_color_viridis_c() +
  expand_limits(x = 0, y = 0) +
  geom_smooth(method = "loess", col="black", lwd=1.5) +
  theme_classic()
ggsave(paste0(plotdir, "Seurat.3/pseudotime_p57.pdf"), plot=g3, height=6, width=7.5)

# draw line graph - z-score
load(paste0(rdatadir, "Seurat/gif_seurat.3.RData"))
normCounts <- gif_seurat.3@assays$RNA@data
scaledCounts <- t(scale(t(normCounts)))

df.gif = df.gif_to_pgc[colnames(gif_seurat.combined.3),]
cell.order = rownames(df.gif[order(df.gif$pseudotime),])

df.gif = df.gif[cell.order,]

scaledCounts = scaledCounts[, cell.order]

library(zoo)
windowsize = ncol(gif_seurat.combined.3) * 0.05
windowsize=100

genes = c("Cdkn1c", "Mki67", "Cblif")
geneids = ensemblGenes[ensemblGenes$external_gene_name %in% genes, ]$ensembl_gene_id

TSmat = matrix(data=0, nrow=length(geneids),
               ncol=ncol(scaledCounts), dimnames=list(genes, cell.order))

for(i in 1:length(geneids)){
  y = scaledCounts[geneids[i],]
  Y_hat = rollapply(y, windowsize, mean, aline='center', fill='extend')
  TSmat[genes[i],] = Y_hat
}

identical(colnames(TSmat), rownames(df.gif))

df.gif$p57 = TSmat["Cdkn1c",]
df.gif$ki67 = TSmat["Mki67",]
df.gif$gif = TSmat["Cblif",]

library(ghibli)
p57.color = ghibli_palette("KikiMedium", 7)[3]
gif.color = ghibli_palette("MononokeMedium", 7)[6]
ki67.color = ghibli_palette("KikiMedium", 7)[4]

gt <- ggplot(df.gif) +
  geom_line(aes(x = pseudotime, y=p57), col=p57.color) +
  geom_line(aes(x = pseudotime, y=ki67), col=ki67.color) +
  geom_line(aes(x = pseudotime, y=gif), col=gif.color) +
  ylab("") + 
  theme_classic()
ggsave(paste0(plotdir, "Seurat.3/auto.expr.png"), plot=gt, height=2, width=3)

gt1 <- ggplot(df.gif) +
  geom_line(aes(x = pseudotime, y=p57), col=p57.color) +
  ylab("") + 
  theme_classic()
ggsave(paste0(plotdir, "Seurat.3/auto.expr.p57.png"), plot=gt1, height=2, width=3)

gt2 <- ggplot(df.gif) +
  geom_line(aes(x = pseudotime, y=ki67), col=ki67.color) +
  ylab("") + 
  theme_classic()
ggsave(paste0(plotdir, "Seurat.3/auto.expr.ki67.png"), plot=gt2, height=2, width=3)

gt3 <- ggplot(df.gif) +
  geom_line(aes(x = pseudotime, y=gif), col=gif.color) +
  ylab("") + 
  theme_classic()
ggsave(paste0(plotdir, "Seurat.3/auto.expr.gif.png"), plot=gt3, height=2, width=3)

par(mfrow=c(3,1))
gt1
gt2

library(ggpubr)
ggarrange(gt1 + ylab("p57") + xlab(""), 
          gt3 + ylab("Gif") + xlab(""), 
          gt2 + ylab("Ki67"),
          #labels = c("p57", "Gif", "Ki67"),
          ncol = 1, nrow = 3)
