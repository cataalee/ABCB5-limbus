#' ---
#' title: "FACS of ABCB5+ limbal cells"
#' author: "Catherine A. A. Lee"
#' date: "June 21, 2021"
#' ---

########################################################################################
#################################### Load packages ####################################
########################################################################################

# install.packages('BiocManager')
# install.packages('Seurat')
# install.packages('devtools')
# install.packages('RcppArmadillo')
# install.packages('harmony')

library(Seurat)
library(devtools)
library(harmony)
library(dplyr)
library(cowplot)
library(ggplot2)
library(RColorBrewer)

########################################################################################
############################ Load data and replicate info #############################
########################################################################################

limbus.data <- Read10X(data.dir 
                       = "agg-78-108-109/outs/filtered_gene_bc_matrices_mex/GRCh38/")

cellcodes <- as.data.frame(limbus.data@Dimnames[[2]])
colnames(cellcodes) <- "barcodes"
rownames(cellcodes) <- cellcodes$barcodes

cellcodes$libcodes <- as.factor(gsub(pattern=".+-", replacement="", cellcodes$barcodes))

sampleidentity <- cellcodes["libcodes"]

# Create Seurat object
limbus <- CreateSeuratObject(counts = limbus.data, project = "limbus3reps", min.cells 
                             = 3, min.features = 200, meta.data = sampleidentity)

########################################################################################
########################## Standard pre-processing workflow ###########################
########################################################################################

# The [[ operator can add columns to object metadata. 
limbus[["percent.mt"]] <- PercentageFeatureSet(limbus, pattern = "^MT-")

limbus <- subset(limbus, subset = nFeature_RNA > 400 & nFeature_RNA < 4000 & 
                   percent.mt < 10)

########################################################################################
################### Remove contaminating melanocytes and stromal cells #################
########################################################################################

# Before removal
#VlnPlot(limbus, features = c("MLANA", "TYRP1","KERA", "LUM"))

# Cutoff for these markers >2
limbus <- subset(x = limbus, subset = MLANA < 2 & TYRP1 < 2 & KERA < 2 & LUM < 2)
#VlnPlot(limbus, features = c("MLANA", "TYRP1","KERA", "LUM"))


########################################################################################
#################################### Norm and Cluster ##################################
########################################################################################

limbus.cleaned.norm <- NormalizeData(limbus, normalization.method = "LogNormalize", 
                                     scale.factor = 10000)

all.genes <- rownames(limbus.cleaned.norm)
limbus.cleaned.norm <- ScaleData(limbus.cleaned.norm, features = all.genes)

limbus.cleaned.norm <- FindVariableFeatures(limbus.cleaned.norm, selection.method 
                                            = "vst", nfeatures = 2000)

limbus.cleaned.norm <- RunPCA(limbus.cleaned.norm, features = VariableFeatures(object 
                                            = limbus.cleaned.norm))

######################################## Harmony #######################################

# options(repr.plot.height = 5, repr.plot.width = 12)
# p1 <- DimPlot(object = limbus.cleaned.norm, reduction = "pca", pt.size = .1, 
#               group.by = "libcodes")
# p2 <- VlnPlot(object = limbus.cleaned.norm, features = "PC_1", group.by = 
#                 "libcodes", pt.size = .1)
# plot_grid(p1,p2)

#options(repr.plot.height = 2.5, repr.plot.width = 6)
limbus.cleaned.norm <- limbus.cleaned.norm %>% 
  RunHarmony("libcodes", plot_convergence = TRUE)

# options(repr.plot.height = 5, repr.plot.width = 12)
# p1 <- DimPlot(object = limbus.cleaned.norm, reduction = "harmony", pt.size = .1, 
#               group.by = "libcodes")
# p2 <- VlnPlot(object = limbus.cleaned.norm, features = "harmony_1", group.by 
#               = "libcodes", pt.size = .1)
# plot_grid(p1,p2)

################################# Downstream Analysis ###################################
limbus.harmony <- limbus.cleaned.norm %>% 
  RunUMAP(reduction = "harmony", dims = 1:20) %>%
  FindNeighbors(reduction = "harmony", dims = 1:20) %>% 
  FindClusters(resolution = 0.01) %>% 
  identity()

# options(repr.plot.height = 4, repr.plot.width = 10)
# DimPlot(limbus.harmony, reduction = "umap", group.by = "libcodes", pt.size = .1, 
#                              split.by = 'libcodes')
# 
# options(repr.plot.height = 4, repr.plot.width = 6)
# DimPlot(limbus.harmony, reduction = "umap", group.by = "libcodes", pt.size = .1)

options(repr.plot.height = 4, repr.plot.width = 6) 
DimPlot(limbus.harmony, reduction = "umap", label = FALSE, pt.size = .1)

# find markers for every cluster compared to all remaining cells, 
#report both positive and negative
limbus.markers <- FindAllMarkers(limbus.harmony, only.pos = TRUE, min.pct = 0.25, 
                                 logfc.threshold = 0.25)
limbus.markers %>% group_by(cluster) %>% top_n(n = 100, wt = avg_log2FC)

# Write cluster markers to file
#write.csv(limbus.markers, "limbus.markers.csv")


################################# Enhanced UMAP Plot ####################################

# Renumber clusters to start at 1 instead of 0
current.cluster.ids <- c(0, 1)
new.cluster.ids <- c("qLSC", "aLSC")
limbus.harmony@meta.data$seurat_clusters <- plyr::mapvalues(x = limbus.harmony
                                  @meta.data$seurat_clusters, from = 
                                  current.cluster.ids, to = new.cluster.ids)

colors = c('cornflowerblue', 'darkgoldenrod2')

options(repr.plot.height = 4, repr.plot.width = 6) 
DimPlot(limbus.harmony, reduction = "umap", group.by = "seurat_clusters", 
        cols = colors, label = TRUE, pt.size = .3, label.size = 8)


################################## Stacked Violin Plots ################################

## remove the x-axis text and tick
## plot.margin to adjust the white space between each plot.
## ... pass any arguments to VlnPlot in Seurat
modify_vlnplot<- function(obj, 
                          feature, 
                          pt.size = 0, 
                          plot.margin = unit(c(-0.75, 0, -0.75, 0), "cm"),
                          ...) {
  p<- VlnPlot(obj, features = feature, pt.size = pt.size, ... )  + 
    xlab("") + ylab(feature) + ggtitle("") + 
    theme(legend.position = "none", 
          axis.text.x = element_blank(), 
          axis.ticks.x = element_blank(), 
          axis.title.y = element_text(size = rel(1), angle = 0), 
          axis.text.y = element_text(size = rel(1)), 
          plot.margin = plot.margin ) 
  return(p)
}

## extract the max value of the y axis
extract_max<- function(p){
  ymax<- max(ggplot_build(p)$layout$panel_scales_y[[1]]$range$range)
  return(ceiling(ymax))
}


## main function
StackedVlnPlot<- function(obj, features,
                          pt.size = 0, 
                          plot.margin = unit(c(-0.75, 0, -0.75, 0), "cm"),
                          ...) {
  
  plot_list<- purrr::map(features, function(x) modify_vlnplot(obj = obj,feature
                                                              = x, ...))
  
  # Add back x-axis title to bottom plot. patchwork is going to support this?
  plot_list[[length(plot_list)]]<- plot_list[[length(plot_list)]] +
    theme(axis.text.x=element_text(), axis.ticks.x = element_line())
  
  # change the y-axis tick to only max value 
  ymaxs<- purrr::map_dbl(plot_list, extract_max)
  plot_list<- purrr::map2(plot_list, ymaxs, function(x,y) x + 
                            scale_y_continuous(breaks = c(y)) + 
                            expand_limits(y = y))
  
  p<- patchwork::wrap_plots(plotlist = plot_list, ncol = 1)
  return(p)
}

# Stacked violin plots
features<- c("PAX6", "KRT15", "TP63")
features<- c("BCAM", "MYC", "CAV1", "ITGA6", "ITGB1", "ITGB4")
features<- c("CDKN2B", "CCND2")
features<- c("LAMA3", "LAMA5")
features<- c("CYR61")

StackedVlnPlot(obj = limbus.harmony, features = features, group.by 
               = "seurat_clusters", cols = colors)

FeaturePlot(object = limbus.harmony, features = features, combine 
               = TRUE) & NoAxes() & NoLegend()

########################################################################################
####################################### Slingshot ######################################
########################################################################################

# BiocManager::install("slingshot")

library(slingshot)
library(scales)

# run without star.clus
limbus.slingshot <- slingshot(Embeddings(limbus.harmony, "umap"), clusterLabels 
                              = limbus.harmony@meta.data$seurat_clusters)

current.cluster.ids <- c("qLSC", "aLSC")
new.cluster.ids <- c(1, 2)
limbus.harmony@meta.data$seurat_clusters <- plyr::mapvalues(x 
                              = limbus.harmony@meta.data$seurat_clusters, 
                              from = current.cluster.ids, to = colors)

colors = as.vector(limbus.harmony@meta.data$seurat_clusters)

plot(reducedDim(limbus.slingshot), col = colors, pch = 16, cex = 0.4)
lines(limbus.slingshot, lwd = 2, type = 'lineages', col = 'black')
