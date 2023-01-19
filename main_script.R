library(scLCA)
# read into raw count data
load('mycount.RData')
# load meta information
load('mymeta.RData')
# Run single cell clustering algorithm scLCA
myclust.res <- myscLCA(datmatrix=count, clust.max=100, datBatch=sample$sampleID, topN=1)
# Cosine distance between cells
myDist.best <- cosDist(myclust.res$pca$rotation[,myclust.res$bestfactors]) 

# Create Seurat object.
library(Seurat)
sample <- Seurat::CreateSeuratObject(counts = count)
# Add meta information
sample <- AddMetaData(sample, metadata = mymeta.df)
# Add cluster information
sample$clustID <- myclust.res[[2]][[1]]$membership
# Collect QC metrics and compute percentage of various type of genes.
sample <- PercentageFeatureSet(sample, "^mt-", col.name = "percent_mito")
#sample <- PercentageFeatureSet(sample, "^Hb[^(p)]", col.name = "percent_hb")
#sample <- PercentageFeatureSet(sample, "^Rp[sl]", col.name = "percent_ribo")
#sample <- PercentageFeatureSet(sample, "Pecam1|Pf4", col.name = "percent_plat")
# Generate QC filter.
mask1 <- sample$nCount_RNA >= 500
mask2 <- sample$percent.mt <= 20
mask <- mask1 & mask2
# Apply filter
sample <- sample[, mask]
# Remove cluster 5 with more than half cells with low QC metric
sample <- sample[,sample$clustID != 5]
# Update distance matrix
myDist.best <-  myDist.best[ match(colnames(sample),colnames(count)), match(colnames(sample),colnames(count)) ]

# Single cell Differential Expression analysis
source("svaBasedAnalysis.R")
method = "SVA_auto_max20_UMI_scran-B10"
keepAll = T
imbalanceRatio = 1e8
outputRoot = "output"
if (!file.exists(outputRoot)) {
     dir.create(outputRoot, recursive = T)
}

mymeta.df <- sample@meta.data
count <- sample@assays$RNA@counts

mymeta.df$idx = 1:nrow(mymeta.df)
source("DEAdjustForBatch.R")
# Normal vs Tumor
g1 = mymeta.df[ mymeta.df$treatment == "Tumor", "idx"]
g2 = mymeta.df[ mymeta.df$treatment == "Normal", "idx"]
group =  c(rep("g1", length(g1)), rep("g2", length(g2)))
index = c(g1, g2)
countUsed = count[, index, drop = F]
srcUsed = src.all[index]
plateID = paste0(group, "_", srcUsed)
#print(table(group, srcUsed))
resultFile = paste0(outputRoot, "/tumor_vs_normal.", method, ".DEResult.tsv")
SVFile = paste0(outputRoot, "/tumor_vs_normal.", method,  ".usedSV.tsv")
svaBasedAnalysis(countUsed, group, plateID, as.character(srcUsed), method, resultFile, SVFile, imbalanceRatio)
# Cluster 16 vs other tumor associated clusters
g1 = mymeta.df[ mymeta.df$clustID %in% c(16),"idx"]
g2 = mymeta.df[ mymeta.df$clustID %in% c(2,3,7,9,12),"idx"]
group =  c(rep("g1", length(g1)), rep("g2", length(g2)))
index = c(g1, g2)
countUsed = count[, index, drop = F]
srcUsed = src.all[index]
plateID = paste0(group, "_", srcUsed)
resultFile = paste0(outputRoot, "/C16_vs_others.", method, ".DEResult.tsv")
SVFile = paste0(outputRoot, "/C16_vs_others.", method,  ".usedSV.tsv")
svaBasedAnalysis(countUsed, group, plateID, as.character(srcUsed), method, resultFile, SVFile, imbalanceRatio)
# Cluster 7 vs other tumor associated clusters
g1 = mymeta.df[ mymeta.df$clustID %in% c(7),"idx"]
g2 = mymeta.df[ mymeta.df$clustID %in% c(2,3,9,12,16),"idx"]
group =  c(rep("g1", length(g1)), rep("g2", length(g2)))
index = c(g1, g2)
countUsed = count[, index, drop = F]
srcUsed = src.all[index]
plateID = paste0(group, "_", srcUsed)
resultFile = paste0(outputRoot, "/C7_vs_others.", method, ".DEResult.tsv")
SVFile = paste0(outputRoot, "/C7_vs_others.", method,  ".usedSV.tsv")
svaBasedAnalysis(countUsed, group, plateID, as.character(srcUsed), method, resultFile, SVFile, imbalanceRatio)

# TSNE plot
library(Rtsne)
tSNE <- Rtsne(myDist.best, is_distance = TRUE, dims = 2)
myTSNE.df <- data.frame(tsneX=tSNE$Y[,1],tsneY=tSNE$Y[,2], cellID <- rownames(myDist.best))
myTSNE.df$clustID <- sample$clustID
myTSNE.df$treatment <- sample$treatment
# Plot
library(ggplot2)
library(RColorBrewer)
# Color code
load(mycol.df)
myg <- ggplot(myTSNE.df,aes(tsneX,tsneY,col=factor(clustID)))
myg + geom_point(size=1.2)+xlab("tsne-X")+ylab('tsne-Y')+ scale_color_manual( values=(mycol.df$col)) + theme_classic(base_size = 24) + theme(legend.position = "none")
ggsave('TSNE_16_clusters.pdf',width = 14,height = 14)
myg <- ggplot(myTSNE.df,aes(tsneX,tsneY,col=factor(treatment)))
myg + geom_point(size=1.2)+xlab("tsne-X")+ylab('tsne-Y')+ theme_classic(base_size = 24) + theme(legend.position = "none")
ggsave('TSNE_treatment.pdf',width = 14,height = 14)

# Dotplot for top genes
mydotplot <- SCpubr::do_DotPlot(sample = sample, features = mytopgenes, group.by = "clustID", flip = T, cluster.idents = F,rotate_x_axis_labels = 0, dot.scale = 12,colors.use = c("skyblue","orange", "red"),scale.by = "size")
mydotplot
ggsave("Dotplot_by_clusters.pdf",width = 10,height = 7)

mydotplot2 <- SCpubr::do_DotPlot(sample = sample,features = mytopgenes,group.by = "sampleID", flip = T, cluster.idents = F,rotate_x_axis_labels = 0, dot.scale = 12, colors.use = c("skyblue","orange", "red"),scale.by = 'size')
mydotplot2
ggsave("Dotplot_by_sample.pdf",width = 8,height = 7)

# SingleR and celldex mapping
library(SingleR)
myref <- celldex::MouseRNAseqData()
mypred <- SingleR(test=as.SingleCellExperiment(sample), ref=ref, labels=ref$label.main)
mytab <- table(Assigned=mypred$pruned.labels, Cluster=sample$clustID)
# Heatmap  
library(pheatmap)
pheatmap(log2(mytab+10), color=colorRampPalette(c("white", "blue"))(101))

# Spatial data
library(hdf5r)
library(Seurat)
library(ggplot2)
library(patchwork)
library(dplyr)
myspatial <- Load10X_Spatial("outs/",filename = "raw_feature_bc_matrix.h5")
#=== spatial data here ===#
myspatial <- Load10X_Spatial("outs/",filename = "raw_feature_bc_matrix.h5",filter.matrix = TRUE)
myspatial <- PercentageFeatureSet(m, "^mt-", col.name = "percent_mito")
myspatial <- SCTransform(myspatial, assay = "Spatial", verbose = FALSE)
SpatialFeaturePlot(m, features = c("nCount_Spatial", "nFeature_Spatial", "percent_mito"))

# Combine spatial data
mycombined <- merge(mycombined,myspatial,add.cell.ids=c("","spatialID"))
mycombined <- PercentageFeatureSet(mycombined, "^mt-", col.name = "percent_mito")
mycombined <- SCTransform(mycombined, assay = "Spatial", verbose = FALSE)
myspx.ls <- SplitObject(mycombined, split.by = "sampleID")

# Spatial plots
library(stringr)
myfeat <- c("MKI67", "PCNA", "TDO2", "IDO2", "CYP2E1", "Hpx", "VEGFA", "Pecam1", "CD4", "CD8A", "CD14", "CD68", "CD3G", "COL1A1")
myfeat <- str_to_sentence(myfeat)
for(jj in 1:length(myfeat)){
mytmpfeat.ls <- list()
for(j in c(1:7)){
m = subset(myspx.ls[[j]], subset = nCount_Spatial > 500  & percent_mito < 20);
mytmpfeat.ls[[j]] <- data.frame(value=colMeans(m@assays$SCT@data[rownames(m@assays$SCT@data) %in% c(myfeat[jj]),,drop=F]), sampleID=j)
}
mytmpfeat <- do.call("rbind",mytmpfeat.ls)
mytmpfeat$value <- scale(mytmpfeat$value)
mytmpfeat.ls <- split(mytmpfeat$value,mytmpfeat$sampleID)
for(ii in 1:length(mytmpfeat.ls)){thistmp <- mytmpfeat.ls[[ii]]; thistmp[thistmp > 2] =2; thistmp[thistmp < -2] = -2; if(min(thistmp) > -2){thistmp[which.min(thistmp)] = -2}; if(max(thistmp) < 2){thistmp[which.max(thistmp)] = 2}; mytmpfeat.ls[[ii]] <- thistmp;  }
mytmpfeat$value2 <- unlist(mytmpfeat.ls)
myplot.ls <- list()
k <- 1
for(j in c(1,5,6,7,2,3,4)){
m = subset(myspx.ls[[j]], subset = nCount_Spatial > 1000  & percent_mito < 20);
m <- AddMetaData(
object = m,
metadata = mytmpfeat$value2[mytmpfeat$sampleID == j],
col.name = "mytmp"
)
myplot.ls[[k]] <- SpatialFeaturePlot(m, features = "mytmp",image.alpha = .6)[[j]]+ ggtitle(sampleInfo$MouseID[j])+scale_fill_distiller(palette = "Spectral", rescaler = ~ scales::rescale_mid(.x, mid = 0) )
k <- k+1;
}
thisplot <- ggarrange(myplot.ls[[1]], myplot.ls[[2]], myplot.ls[[3]], myplot.ls[[4]], myplot.ls[[5]], myplot.ls[[6]], myplot.ls[[7]],nrow=2,ncol=4, common.legend = T,legend = "none")
annotate_figure(thisplot, top = text_grob(myfeat[jj], color = "black", face = "bold", size = 24))
ggsave(paste(myfeat[jj],"spatial_plot_dec19_2022round1.pdf",sep="-"),width = 10,height = 7)
}