# 20230922_rat_retina_with_GEO.R

#### SETUP ####
  #### Samples and path variables ####
# Rat retina treated with human cells
# WT rat retina samples from GEO GSE209872
# ALIGNED to RAT ref rnor
# samples used:
#
# rnor_08182_P90_Treated_rat_human_cells_outs
# rnor_13933_P60_Treated_rat_cells_outs
# rnor_13933_P60_Untreated_rat_cells_outs
# rnor_16573_P90_Treated_rat_cells_outs
# rnor_16573_P90_Untreated_rat_cells_outs
# rnor_18482_P60_Treated_rat_cells_outs
# rnor_18482_P60_Untreated_rat_cells_outs
# rnor_19328_P60_Treated_rat_cells_outs
# rnor_19328_P60_Untreated_rat_cells_outs
# rnor_SRR20570827_outs
# rnor_SRR20570828_outs

date <- "20230922"
project <- "Rat_retina_with_GEO"
datadir <- "/Users/bells/Library/CloudStorage/Box-Box/20230922_rat_retina_with_GEO"
sourcedir1 <- "/Users/bells/Library/CloudStorage/Box-Box/20230317_rnor_outs_copy"
sourcedir2 <- "/Users/bells/Library/CloudStorage/Box-Box/20230926_rnor_wang_geo_outs_copy"

  #### Load required packages ####

library(Seurat)
library(tidyverse)
library(Matrix)
library(cowplot)
library(viridis)
library(patchwork)
library(pheatmap)
library(irlba)
library(beepr)
library(SeuratDisk)
library(EnhancedVolcano)
library(Vennerable)

  #### Set working dir and save session info ####
setwd(datadir)

sessionInfo()

sink(paste0(date,"_",project,"_devtools_sessionInfo3.txt"))
devtools::session_info()
sink()

sink(paste0(date,"_",project,"_sessionInfo3.txt"))
sessionInfo()
sink()

  #### Functions for saving images in hires ####

hirestiff <- function(saveas){
  ggsave(
    saveas,
    plot = last_plot(),
    device = "tiff",
    scale = 1,
    dpi = 600,
    limitsize = TRUE
  )
}

lowrestiff <- function(saveas){
  ggsave(
    saveas,
    plot = last_plot(),
    device = "tiff",
    scale = 1,
    dpi = 100,
    limitsize = TRUE
  )
}

hirestiffsquare <- function(saveas){
  ggsave(
    saveas,
    plot = last_plot(),
    device = "tiff",
    scale = 1,
    dpi = 600,
    limitsize = TRUE,
    bg = "white",
    width = 9,
    height = 9,
    units = c("in")
  )
}

lowrestiffsquare <- function(saveas){
  ggsave(
    saveas,
    plot = last_plot(),
    device = "tiff",
    scale = 1,
    dpi = 100,
    limitsize = TRUE,
    bg = "white",
    width = 9,
    height = 9,
    units = c("in")
  )
}

savepdf <- function(saveas){
  ggsave(
    saveas,
    plot = last_plot(),
    device = "pdf",
    scale = 1,
    dpi = 600,
    limitsize = TRUE
  )
}

savepdfsquare <- function(saveas){
  ggsave(
    saveas,
    plot = last_plot(),
    device = "pdf",
    scale = 1,
    dpi = 600,
    width = 9,
    height = 9,
    limitsize = TRUE
  )
}

  #### restore prefixes ####

# only save after determining cluster resolution
# must enter at minimum the date variable and datadir

setwd(datadir)

# for sct data:

# save(datadir, date, ident_res, meaningful.PCs,
#      prefixPC, prefixPCres, project, res, sourcedir,
#      hirestiff, hirestiffsquare, 
#      lowrestiff, lowrestiffsquare,
#      savepdf, savepdfsquare,
#      file = paste0(datadir,"/",date,"_",project,"_","prefixes.rdata"))
 
# load(paste0(datadir,"/",date,"_",project,"_","prefixes.rdata"))

# for integrated data:

# save(datadir, date, ident_res, meaningful.PCs, prefixPC, prefixPCres,
#      project, res, sourcedir1, sourcedir2,
#      hirestiff, hirestiffsquare,
#      lowrestiff, lowrestiffsquare,
#      savepdf, savepdfsquare,
#      palopal, sorted_palopal,
#      file = paste0(datadir,"/",date,"_",project,"_","prefixes.rdata"))

project <- "Rat_retina_with_GEO_int"

load(paste0(datadir,"/",date,"_",project,"_","prefixes.rdata"))

seurat_sct <- readRDS(file = paste0(prefixPCres,"_seurat_after_RNAnorm_CELLTYPES.rds"))

#### LOAD DATA ####
  #### Load data from Cellranger output ####

sample_list <- c(
  "rnor_08182_P90_Treated_rat_human_cells_outs",
  "rnor_13933_P60_Treated_rat_cells_outs",
  "rnor_13933_P60_Untreated_rat_cells_outs",
  "rnor_16573_P90_Treated_rat_cells_outs",
  "rnor_16573_P90_Untreated_rat_cells_outs",
  "rnor_18482_P60_Treated_rat_cells_outs",
  "rnor_18482_P60_Untreated_rat_cells_outs",
  "rnor_19328_P60_Treated_rat_cells_outs",
  "rnor_19328_P60_Untreated_rat_cells_outs",
  "rnor_SRR20570827_outs",
  "rnor_SRR20570828_outs"
)

# since we have data from two separate source directories
# we make two filelists and pathlists
filelist1 <- list.files(path = sourcedir1)
filelist1 <- as.data.frame(filelist1, files = list.files(path = sourcedir1))
path_list1 <- paste0(sourcedir1,"/",filelist1[,1],"/filtered_feature_bc_matrix")

filelist2 <- list.files(path = sourcedir2)
filelist2 <- as.data.frame(filelist2, files = list.files(path = sourcedir2))
path_list2 <- paste0(sourcedir2,"/",filelist2[,1],"/filtered_feature_bc_matrix")

path_list <- c(path_list1,path_list2)
path_list

# check that all files exist
all(file.exists(path_list))

# check that the sample list and path list are in the same order
for (i in 1:length(sample_list)) {
  TF <- grepl(sample_list[i],path_list[i])
  print(paste0(sample_list[i]," ",TF))
}

# make an empty list to store data into
data_list <- vector("list")

# make a list of all the 10X data
for (i in 1:length(sample_list)) {
  # read in cellranger output
  cellranger <- Read10X(data.dir = path_list[i])
  # append data to the list
  data_list <- c(data_list,cellranger)
}

updated_samps <- c("08182.p90.TR",
                   "13933.p60.TR",
                   "13933.p60.UT",
                   "16573.p90.TR",
                   "16573.p90.UT",
                   "18482.p60.TR",
                   "18482.p60.UT",
                   "19328.p60.TR",
                   "19328.p60.UT",
                   "GEO27.p00.WT",
                   "GEO28.p00.WT")

names(data_list) <- updated_samps

  #### Seurat Object ####

# Initialize the Seurat object with the raw (non-normalized data).
# Keep all genes expressed in >= 1 cell

# make an empty list to store data into
seurat_list <- vector("list")

# make a list of all the 10X data
for (i in 1:length(data_list)) {
  # create seurat object
  s.data <- CreateSeuratObject(counts = data_list[[i]], 
                               project = updated_samps[i], 
                               min.cells = 1, 
                               min.features = 0)
  # append data to the list
  seurat_list <- c(seurat_list,s.data)
}

names(seurat_list) <- updated_samps

seurat_list

  #### Add metadata ####

# Create lists for metadata attributes
time_point <- c(
  "08182.p90.TR" = "p90",
  "13933.p60.TR" = "p60",
  "13933.p60.UT" = "p60",
  "16573.p90.TR" = "p90",
  "16573.p90.UT" = "p90",
  "18482.p60.TR" = "p60",
  "18482.p60.UT" = "p60",
  "19328.p60.TR" = "p60",
  "19328.p60.UT" = "p60",
  "GEO27.p00.WT" = "WT",
  "GEO28.p00.WT" = "WT"
)

treatment <- c(
  "08182.p90.TR" = "Cell_Treated",
  "13933.p60.TR" = "Cell_Treated",
  "13933.p60.UT" = "Untreated",
  "16573.p90.TR" = "Cell_Treated",
  "16573.p90.UT" = "Untreated",
  "18482.p60.TR" = "Cell_Treated",
  "18482.p60.UT" = "Untreated",
  "19328.p60.TR" = "Cell_Treated",
  "19328.p60.UT" = "Untreated",
  "GEO27.p00.WT" = "Untreated",
  "GEO28.p00.WT" = "Untreated"
)

batch <- c(
  "08182.p90.TR" = "08182",
  "13933.p60.TR" = "13933",
  "13933.p60.UT" = "13933",
  "16573.p90.TR" = "16573",
  "16573.p90.UT" = "16573",
  "18482.p60.TR" = "18482",
  "18482.p60.UT" = "18482",
  "19328.p60.TR" = "19328",
  "19328.p60.UT" = "19328",
  "GEO27.p00.WT" = "GSE209872",
  "GEO28.p00.WT" = "GSE209872"
)

orig.ident <- c(
  "08182.p90.TR" = "08182.p90.TR",
  "13933.p60.TR" = "13933.p60.TR",
  "13933.p60.UT" = "13933.p60.UT",
  "16573.p90.TR" = "16573.p90.TR",
  "16573.p90.UT" = "16573.p90.UT",
  "18482.p60.TR" = "18482.p60.TR",
  "18482.p60.UT" = "18482.p60.UT",
  "19328.p60.TR" = "19328.p60.TR",
  "19328.p60.UT" = "19328.p60.UT",
  "GEO27.p00.WT" = "GEO27.WT.UT",
  "GEO28.p00.WT" = "GEO28.WT.UT"
)

# Iterate through the Seurat objects and assign metadata
for (key in names(time_point)) {
  seurat_list[[key]][["time_point"]] <- time_point[key]
  seurat_list[[key]][["treatment"]] <- treatment[key]
  seurat_list[[key]][["batch"]] <- batch[key]
  seurat_list[[key]][["time_point_treatment"]] <- paste(time_point[key], treatment[key], sep = "_")
  seurat_list[[key]][["time_point_batch"]] <- paste(time_point[key], batch[key], sep = "_")
  seurat_list[[key]][["treatment_batch"]] <- paste(treatment[key], batch[key], sep = "_")
  seurat_list[[key]][["time_point_treatment_batch"]] <- paste(time_point[key], treatment[key], batch[key], sep = "_")
  seurat_list[[key]][["orig.ident"]] <- orig.ident[key]
}

  #### Merge all data into a single seurat object ####

tenx1 <- merge(seurat_list$"08182.p90.TR", 
               y = c(
                 seurat_list$"13933.p60.TR",
                 seurat_list$"13933.p60.UT",
                 seurat_list$"16573.p90.TR",
                 seurat_list$"16573.p90.UT",
                 seurat_list$"18482.p60.TR",
                 seurat_list$"18482.p60.UT",
                 seurat_list$"19328.p60.TR",
                 seurat_list$"19328.p60.UT",
                 seurat_list$"GEO27.p00.WT",
                 seurat_list$"GEO28.p00.WT"
                 ),
               add.cell.ids = c(
                 "08182.p90.TR",
                 "13933.p60.TR",
                 "13933.p60.UT",
                 "16573.p90.TR",
                 "16573.p90.UT",
                 "18482.p60.TR",
                 "18482.p60.UT",
                 "19328.p60.TR",
                 "19328.p60.UT",
                 "GEO27.p00.WT",
                 "GEO28.p00.WT"
                 ),
               project = project)

tenx1

rm(cellranger,
   data_list,
   filelist1,
   filelist2,
   s.data,
   seurat_list,
   batch,
   i,
   key,
   path_list,
   path_list1,
   path_list2,
   sample_list,
   TF,
   time_point,
   treatment,
   updated_samps)

# save merged seurat object 
# saveRDS(tenx1, file = paste0("./",date,"_",project,"_merged_seurat_prefilter.rds"))

# tenx1 <- read_rds(file = paste0("./",date,"_",project,"_merged_seurat_prefilter.rds"))

#### QC FILTERING ####
  #### QC Setup ####
# copy tenx1 to another object data just in case 
data <- tenx1

# check to make sure pattern is correct and not grepping other genes
grep("^Mt-",rownames(data),value = T) # mitochondrial genes
grep("^Rp[sl]",rownames(data),value = T) # ribosomal genes
grep("^Mrp[sl]",rownames(data),value = T) # mitochondrial-ribosomal genes
grep("^Hb[abq]",rownames(data),value = T) # hemoglobin genes

# Add PercentageFeatureSet of various factors to metadata
data[["percent.mt"]] <- PercentageFeatureSet(data, pattern = "^Mt-")
data[["percent.ribo"]] <- PercentageFeatureSet(data, pattern = "^Rp[sl]")
data[["percent.mt.ribo"]] <- PercentageFeatureSet(data, pattern = "^Mrp[sl]")
data[["percent.hb"]] <- PercentageFeatureSet(data, pattern = "^Hb[abq]")

  #### Prefiltered Plots ####
vp <- VlnPlot(data, features = "nCount_RNA", pt.size = 0) +
  NoLegend() +
  ggtitle("nCount_RNA prefiltered")
vp
pdf(paste0("./",date,"_",project,"_prefilter_","nCount_RNA",".pdf"))
print(vp)
dev.off()

vp <- VlnPlot(data, features = "nFeature_RNA", pt.size = 0) +
  NoLegend() +
  ggtitle("nFeature_RNA prefiltered")
vp
pdf(paste0("./",date,"_",project,"_prefilter_","nFeature_RNA",".pdf"))
print(vp)
dev.off()

vp <- VlnPlot(data, features = "percent.mt", pt.size = 0) +
  NoLegend() +
  ggtitle("percent.mt prefiltered")
vp
pdf(paste0("./",date,"_",project,"_prefilter_","percent.mt",".pdf"))
print(vp)
dev.off()

vp <- VlnPlot(data, features = "percent.ribo", pt.size = 0) +
  NoLegend() +
  ggtitle("percent.ribo prefiltered")
vp
pdf(paste0("./",date,"_",project,"_prefilter_","percent.ribo",".pdf"))
print(vp)
dev.off()

vp <- VlnPlot(data, features = "percent.mt.ribo", pt.size = 0) +
  NoLegend() +
  ggtitle("percent.mt.ribo prefiltered")
vp
pdf(paste0("./",date,"_",project,"_prefilter_","percent.mt.ribo",".pdf"))
print(vp)
dev.off()

vp <- VlnPlot(data, features = "percent.hb", pt.size = 0) +
  NoLegend() +
  ggtitle("percent.hb prefiltered")
vp
pdf(paste0("./",date,"_",project,"_prefilter_","percent.hb",".pdf"))
print(vp)
dev.off()

p1 <- FeatureScatter(data, feature1 = "nCount_RNA", feature2 = "percent.mt") 
p2 <- FeatureScatter(data, feature1 = "nCount_RNA", feature2 = "nFeature_RNA") 
p1 + p2

pdf(paste0("./",date,"_",project,"_prefilter_","scatterplots",".pdf"))
print(p1 + p2)
dev.off()
  #### Additional Prefiltered QC metrics and plots ####
qc.metrics <- as_tibble(data[[c("nCount_RNA",
                                "nFeature_RNA",
                                "percent.mt",
                                "percent.ribo",
                                "percent.mt.ribo",
                                "percent.hb")]],
                        rownames="Cell.Barcode")

sp <- qc.metrics %>%
  arrange(percent.mt) %>%
  ggplot(aes(nCount_RNA,nFeature_RNA,colour=percent.mt)) + 
  geom_point() + 
  scale_color_gradientn(colors=c("black","blue","green2","red","yellow")) +
  ggtitle("QC metrics - percent Mitochondial DNA Prefiltered")
sp
pdf(paste0("./",date,"_",project,"_QC_Metrics_pct.mt_prefilter",".pdf"))
print(sp)
dev.off()

sp <- qc.metrics %>%
  arrange(percent.ribo) %>%
  ggplot(aes(nCount_RNA,nFeature_RNA,colour=percent.ribo)) + 
  geom_point() + 
  scale_color_gradientn(colors=c("black","blue","green2","red","yellow")) +
  ggtitle("QC metrics - percent Ribosomal DNA Prefiltered")
sp
pdf(paste0("./",date,"_",project,"_QC_Metrics_pct.ribo_prefilter",".pdf"))
print(sp)
dev.off()

sp <- qc.metrics %>%
  arrange(percent.mt.ribo) %>%
  ggplot(aes(nCount_RNA,nFeature_RNA,colour=percent.mt.ribo)) + 
  geom_point() + 
  scale_color_gradientn(colors=c("black","blue","green2","red","yellow")) +
  ggtitle("QC metrics - percent Mito-Ribosomal DNA Prefiltered")
sp
pdf(paste0("./",date,"_",project,"_QC_Metrics_pct.mt.ribo_prefilter",".pdf"))
print(sp)
dev.off()

sp <- qc.metrics %>%
  arrange(percent.hb) %>%
  ggplot(aes(nCount_RNA,nFeature_RNA,colour=percent.hb)) + 
  geom_point() + 
  scale_color_gradientn(colors=c("black","blue","green2","red","yellow")) +
  ggtitle("QC metrics - percent Hemoglobin DNA Prefiltered")
sp
pdf(paste0("./",date,"_",project,"_QC_Metrics_pct.hb_prefilter",".pdf"))
print(sp)
dev.off()

sp <- qc.metrics %>%
  ggplot(aes(percent.mt)) + 
  geom_histogram(binwidth = 0.5, fill="yellow", colour="black") +
  ggtitle("Distribution of Percentage of reads from Mitochondria Prefiltered")
sp
pdf(paste0("./",date,"_",project,"_Dist_mito_prefilter",".pdf"))
print(sp)
dev.off()

sp <- qc.metrics %>%
  ggplot(aes(percent.ribo)) + 
  geom_histogram(binwidth = 0.5, fill="yellow", colour="black") +
  ggtitle("Distribution of Percentage of reads from Ribosomes Prefiltered")
sp
pdf(paste0("./",date,"_",project,"_Dist_ribo_prefilter",".pdf"))
print(sp)
dev.off()

sp <- qc.metrics %>%
  ggplot(aes(percent.mt.ribo)) + 
  geom_histogram(binwidth = 0.05, fill="yellow", colour="black") +
  ggtitle("Distribution of Percentage of reads from Mito-Ribosomes Prefiltered")
sp
pdf(paste0("./",date,"_",project,"_Dist_mt_ribo_prefilter",".pdf"))
print(sp)
dev.off()

sp <- qc.metrics %>%
  ggplot(aes(percent.hb)) + 
  geom_histogram(binwidth = 0.5, fill="yellow", colour="black") +
  ggtitle("Distribution of Percentage of reads from Hemoglobin Prefiltered")
sp
pdf(paste0("./",date,"_",project,"_Dist_hb_prefilter",".pdf"))
print(sp)
dev.off()

  #### Add Z-scores  ####
data[["nUMI.z"]] <- scale(data$nCount_RNA)
data[["nGene.z"]] <- scale(data$nFeature_RNA)
data[["percent.mt.z"]] <- scale(data$percent.mt)
data[["percent.ribo.z"]] <- scale(data$percent.ribo)
data[["percent.mt.ribo.z"]] <- scale(data$percent.mt.ribo)
data[["percent.hb.z"]] <- scale(data$percent.hb)

  #### Prefiltered QC data  ####
length(data@meta.data$orig.ident)
mean(data@meta.data$percent.mt)
mean(data@meta.data$percent.ribo)
mean(data@meta.data$percent.mt.ribo)
mean(data@meta.data$percent.hb)
mean(data@meta.data$nCount_RNA)
median(data@meta.data$nCount_RNA)
mean(data@meta.data$nFeature_RNA)
median(data@meta.data$nFeature_RNA)
max(data@meta.data$nCount_RNA)

# save prefilter QC data
sink(paste0("./",date,"_",project,"_prefilter_QC_metrics.txt"))
cat("length")
length(data@meta.data$orig.ident)
cat("mean percent mito")
mean(data@meta.data$percent.mt)
cat("mean percent ribo")
mean(data@meta.data$percent.ribo)
cat("mean percent mito-ribo")
mean(data@meta.data$percent.mt.ribo)
cat("mean hemoglobin")
mean(data@meta.data$percent.hb)
cat("mean percent counts")
mean(data@meta.data$nCount_RNA)
cat("median counts")
median(data@meta.data$nCount_RNA)
cat("mean features")
mean(data@meta.data$nFeature_RNA)
cat("median features")
median(data@meta.data$nFeature_RNA)
cat("max counts")
max(data@meta.data$nCount_RNA)
sink()

  #### Filter cells based on Z-score  ####
data <- subset(data, subset = percent.mt.z < 3)
data <- subset(data, subset = percent.ribo.z < 3)
data <- subset(data, subset = percent.mt.ribo.z < 3)
data <- subset(data, subset = percent.hb.z < 3)

# data <- subset(data, subset = nGene.z < 3)
# data <- subset(data, subset = nUMI.z < 3)

  #### Postfiltered QC data ####
length(data@meta.data$orig.ident)
mean(data@meta.data$percent.mt)
mean(data@meta.data$percent.ribo)
mean(data@meta.data$percent.hb)
mean(data@meta.data$nCount_RNA)
median(data@meta.data$nCount_RNA)
mean(data@meta.data$nFeature_RNA)
median(data@meta.data$nFeature_RNA)
max(data@meta.data$nCount_RNA)

# save postfilter QC data
sink(paste0("./",date,"_",project,"_postfilter_QC_metrics.txt"))
cat("length")
length(data@meta.data$orig.ident)
cat("mean percent mito")
mean(data@meta.data$percent.mt)
cat("mean percent ribo")
mean(data@meta.data$percent.ribo)
cat("mean percent mito-ribo")
mean(data@meta.data$percent.mt.ribo)
cat("mean hemoglobin")
mean(data@meta.data$percent.hb)
cat("mean percent counts")
mean(data@meta.data$nCount_RNA)
cat("median counts")
median(data@meta.data$nCount_RNA)
cat("mean features")
mean(data@meta.data$nFeature_RNA)
cat("median features")
median(data@meta.data$nFeature_RNA)
cat("max counts")
max(data@meta.data$nCount_RNA)
sink()

  #### Postfiltered Plots ####
vp <- VlnPlot(data, features = "nCount_RNA", pt.size = 0, group.by = "orig.ident")  + 
  NoLegend() + 
  ggtitle("nCount_RNA filtered")
vp
pdf(paste0("./",date,"_",project,"_filtered_","nCount_RNA",".pdf"))
print(vp)
dev.off()

vp <- VlnPlot(data, features = "nFeature_RNA", pt.size = 0, group.by = "orig.ident") + 
  NoLegend() + 
  ggtitle("nFeature_RNA filtered")
vp
pdf(paste0("./",date,"_",project,"_filtered_","nFeature_RNA",".pdf"))
print(vp)
dev.off()

vp <- VlnPlot(data, features = "percent.mt", pt.size = 0, group.by = "orig.ident") +
  NoLegend() +
  ggtitle("percent.mt filtered")
vp
pdf(paste0("./",date,"_",project,"_filtered_","percent.mt",".pdf"))
print(vp)
dev.off()

vp <- VlnPlot(data, features = "percent.ribo", pt.size = 0, group.by = "orig.ident") +
  NoLegend() +
  ggtitle("percent.ribo filtered")
vp
pdf(paste0("./",date,"_",project,"_filtered_","percent.ribo",".pdf"))
print(vp)
dev.off()

vp <- VlnPlot(data, features = "percent.mt.ribo", pt.size = 0) +
  NoLegend() +
  ggtitle("percent.mt.ribo filtered")
vp
pdf(paste0("./",date,"_",project,"_filtered_","percent.mt.ribo",".pdf"))
print(vp)
dev.off()

vp <- VlnPlot(data, features = "percent.hb", pt.size = 0) +
  NoLegend() +
  ggtitle("percent.hb filtered")
vp
pdf(paste0("./",date,"_",project,"_filtered_","percent.hb",".pdf"))
print(vp)
dev.off()

p1 <- FeatureScatter(data, feature1 = "nCount_RNA", feature2 = "percent.mt") 
p2 <- FeatureScatter(data, feature1 = "nCount_RNA", feature2 = "nFeature_RNA") 
p1 + p2

pdf(paste0("./",date,"_",project,"_filtered_","scatterplots",".pdf"))
print(p1 + p2)
dev.off()

  #### Additional postfiltered QC metrics and plots ####
qc.metrics <- as_tibble(data[[c("nCount_RNA",
                                "nFeature_RNA",
                                "percent.mt",
                                "percent.ribo",
                                "percent.mt.ribo",
                                "percent.hb")]],
                        rownames="Cell.Barcode")

sp <- qc.metrics %>%
  arrange(percent.mt) %>%
  ggplot(aes(nCount_RNA,nFeature_RNA,colour=percent.mt)) + 
  geom_point() + 
  scale_color_gradientn(colors=c("black","blue","green2","red","yellow")) +
  ggtitle("QC metrics - percent Mitochondial DNA Postfilter")
sp
pdf(paste0("./",date,"_",project,"_QC_Metrics_pct.mt_postfilter",".pdf"))
print(sp)
dev.off()

sp <- qc.metrics %>%
  arrange(percent.ribo) %>%
  ggplot(aes(nCount_RNA,nFeature_RNA,colour=percent.ribo)) + 
  geom_point() + 
  scale_color_gradientn(colors=c("black","blue","green2","red","yellow")) +
  ggtitle("QC metrics - percent Ribosomal DNA Postfilter")
sp
pdf(paste0("./",date,"_",project,"_QC_Metrics_pct.ribo_postfilter",".pdf"))
print(sp)
dev.off()

sp <- qc.metrics %>%
  arrange(percent.mt.ribo) %>%
  ggplot(aes(nCount_RNA,nFeature_RNA,colour=percent.mt.ribo)) + 
  geom_point() + 
  scale_color_gradientn(colors=c("black","blue","green2","red","yellow")) +
  ggtitle("QC metrics - percent Mito-Ribosomal DNA Postfilter")
sp
pdf(paste0("./",date,"_",project,"_QC_Metrics_pct.mt.ribo_postfilter",".pdf"))
print(sp)
dev.off()

sp <- qc.metrics %>%
  arrange(percent.hb) %>%
  ggplot(aes(nCount_RNA,nFeature_RNA,colour=percent.hb)) + 
  geom_point() + 
  scale_color_gradientn(colors=c("black","blue","green2","red","yellow")) +
  ggtitle("QC metrics - percent Hemoglobin DNA Postfilter")
sp
pdf(paste0("./",date,"_",project,"_QC_Metrics_pct.hb_postfilter",".pdf"))
print(sp)
dev.off()

sp <- qc.metrics %>%
  ggplot(aes(percent.mt)) + 
  geom_histogram(binwidth = 0.25, fill="yellow", colour="black") +
  ggtitle("Distribution of Percentage of reads from Mitochondria Postfilter")
sp
pdf(paste0("./",date,"_",project,"_Dist_mito_postfilter",".pdf"))
print(sp)
dev.off()

sp <- qc.metrics %>%
  ggplot(aes(percent.ribo)) + 
  geom_histogram(binwidth = 0.5, fill="yellow", colour="black") +
  ggtitle("Distribution of Percentage of reads from Ribosomes Postfilter")
sp
pdf(paste0("./",date,"_",project,"_Dist_ribo_postfilter",".pdf"))
print(sp)
dev.off()

sp <- qc.metrics %>%
  ggplot(aes(percent.mt.ribo)) + 
  geom_histogram(binwidth = 0.05, fill="yellow", colour="black") +
  ggtitle("Distribution of Percentage of reads from Mito-Ribosomes Postfilter")
sp
pdf(paste0("./",date,"_",project,"_Dist_mt_ribo_postfilter",".pdf"))
print(sp)
dev.off()

sp <- qc.metrics %>%
  ggplot(aes(percent.hb)) + 
  geom_histogram(binwidth = 0.05, fill="yellow", colour="black") +
  ggtitle("Distribution of Percentage of reads from Hemoglobin Postfilter")
sp
pdf(paste0("./",date,"_",project,"_Dist_hb_postfilter",".pdf"))
print(sp)
dev.off()

  #### post QC cleanup ####

tenx1
data

tenx1 <- data

rm(vp,sp,p1,p2,
   qc.metrics,data)

  #### set factor levels ####
seurat_sct <- tenx1

seurat_sct$orig.ident <- factor(
  x = seurat_sct$orig.ident,
  levels = c(
    "GEO27.WT.UT",
    "GEO28.WT.UT",
    "08182.p90.TR",
    "13933.p60.TR",
    "13933.p60.UT",
    "16573.p90.TR",
    "16573.p90.UT",
    "18482.p60.TR",
    "18482.p60.UT",
    "19328.p60.TR",
    "19328.p60.UT"
  )
)

seurat_sct$treatment <- factor(
  x = seurat_sct$treatment,
  levels = c(
    "Untreated",
    "Cell_Treated"
  )
)

seurat_sct$batch <- factor(
  x = seurat_sct$batch,
  levels = c(
    "GSE209872",
    "08182",
    "13933",
    "16573",
    "18482",
    "19328"
  )
)

seurat_sct$time_point <- factor(
  x = seurat_sct$time_point,
  levels = c(
    "WT",
    "p60",
    "p90"
  )
)

seurat_sct$time_point_treatment <- factor(
  x = seurat_sct$time_point_treatment,
  levels = c(
    "WT_Untreated",
    "p60_Cell_Treated",
    "p60_Untreated",
    "p90_Cell_Treated",
    "p90_Untreated"
  )
)

seurat_sct$time_point_batch <- factor(
  x = seurat_sct$time_point_batch,
  levels = c(
    "WT_GSE209872",
    "p60_13933",
    "p60_18482",
    "p60_19328",
    "p90_08182",
    "p90_16573"
  )
)

seurat_sct$treatment_batch <- factor(
  x = seurat_sct$treatment_batch,
  levels = c(
    "Untreated_GSE209872",
    "Untreated_13933",
    "Untreated_18482",
    "Untreated_19328",
    "Untreated_16573",
    "Cell_Treated_13933",
    "Cell_Treated_18482",
    "Cell_Treated_19328",
    "Cell_Treated_08182",
    "Cell_Treated_16573"
  )
)

seurat_sct$time_point_treatment_batch <- factor(
  x = seurat_sct$time_point_treatment_batch,
  levels = c(
    "WT_Untreated_GSE209872",
    "p60_Untreated_13933",
    "p60_Untreated_18482",
    "p60_Untreated_19328",
    "p90_Untreated_16573",
    "p60_Cell_Treated_13933",
    "p60_Cell_Treated_18482",
    "p60_Cell_Treated_19328",
    "p90_Cell_Treated_08182",
    "p90_Cell_Treated_16573"
  )
)

rm(tenx1)

#### SCTransform ####
seurat_sct

seurat_sct <- SCTransform(seurat_sct,
                          ncells = 11848, #total cells / 10 (118478/10 = 11848)
                          vars.to.regress = c("percent.mt","percent.ribo"), 
                          verbose = TRUE)

  #### PCA ####
# this performs PCA on the seurat object
seurat_sct <- RunPCA(seurat_sct, npcs = 50, verbose = TRUE)

# make PC coordinate object a data frame
xx.coord <- as.data.frame(seurat_sct@reductions$pca@cell.embeddings)

# make PC feature loadings object a data frame
xx.gload <- as.data.frame(seurat_sct@reductions$pca@feature.loadings)

# calculate eigenvalues for arrays

# generate squares of all sample coordinates
sq.xx.coord <- as.data.frame(xx.coord^2)
# create empty list for eigenvalues first
eig <- c()
# calculate the eigenvalue for each PC in sq.xx.coord by taking the sqrt of the sum of squares
for(i in 1:ncol(sq.xx.coord))
  eig[i] = sqrt(sum(sq.xx.coord[,i]))
# calculate the total variance by adding up all the eigenvalues
sum.eig <- sum(eig)
# calculate the expected contribution of all PCs if they all contribute equally to the total variance
expected.contribution <- sum.eig/(length(xx.coord)-1)
# return the number of principal components with an eigenvalue greater than expected by equal variance
meaningful.PCs <- sum(eig > expected.contribution)

# create empty list for eigenvalue percentage
eig.percent <- c()
# calculate the percentage of the total variance by each PC eigenvalue
for(i in 1:length(eig))
  eig.percent[i] = 100*eig[i]/sum.eig
# sum of all eig.percent should total to 100
sum(eig.percent)
# create empty list for scree values
scree <- c()
# calculate a running total of variance contribution
for(i in 1:length(eig))
  if(i == 1) scree[i] = eig.percent[i] else scree[i] = scree[i-1] + eig.percent[i]

# create data frame for eigenvalue summaries
eigenvalues <- data.frame("PC" = colnames(xx.coord), "eig" = eig, "percent" = eig.percent, "scree" = scree)

# write csv for eigenvalues
# write.csv(eigenvalues, file = paste0("./",date,"_",project,"_PCA_eigenvalues.csv"), row.names = F)

# plot scree values
plot(eigenvalues$percent, ylim = c(0,100), type = "S", xlab = "PC", ylab = "Percent of variance",
     main = paste0(date,"_",project," scree plot all samples PCA"))
points(eigenvalues$scree, ylim = c(0,100), type = "p", pch = 16)
lines(eigenvalues$scree)
# add red line to indicate cut-off
cut.off <- 100/(length(eig)-1)
abline(h = cut.off, col = "red")
# add blue line to indicate which PCs are meaningful and kept
abline(v = meaningful.PCs, col = "blue")
text(meaningful.PCs, cut.off, label = paste("cutoff PC",meaningful.PCs),
     adj = c(-0.1, -0.5))

dev.copy(pdf, paste0("./",date,"_",project,"_scree_plot.pdf"))
dev.off()

rm(eigenvalues,sq.xx.coord,xx.coord,xx.gload,cut.off,
   eig,eig.percent,expected.contribution,i,scree,sum.eig)

# meaningful.PCs <- 15

  #### Run UMAP and look at UMAP plots ####

seurat_sct <- RunUMAP(seurat_sct, reduction = "pca", dims = 1:meaningful.PCs, verbose = TRUE)

# update prefixed variable
prefixPC <- paste0("./",date,"_",project,"_",meaningful.PCs,"PCs")

## UMAP plot by sample name ("orig.ident")
DimPlot(seurat_sct,
        reduction = "umap",
        label = FALSE, 
        pt.size = .25, 
        group.by = "orig.ident", 
        raster = F
        )
hirestiff(paste0(prefixPC,"_UMAP_by_sample_hires.tiff"))
lowrestiff(paste0(prefixPC,"_UMAP_by_sample_lowres.tiff"))

DimPlot(seurat_sct, 
        reduction = "umap", 
        label = FALSE, 
        pt.size = .25, 
        group.by = "treatment",
        raster = F
        )
hirestiff(paste0(prefixPC,"_UMAP_by_treatment_hires.tiff"))
lowrestiff(paste0(prefixPC,"_UMAP_by_treatment_lowres.tiff"))

DimPlot(seurat_sct, 
        reduction = "umap", 
        label = FALSE, 
        pt.size = .25, 
        split.by = "treatment", 
        group.by = "treatment",
        raster = F
        )
hirestiff(paste0(prefixPC,"_UMAP_split_by_treatment_hires.tiff"))
lowrestiff(paste0(prefixPC,"_UMAP_split_by_treatment_lowres.tiff"))

DimPlot(seurat_sct, 
        reduction = "umap", 
        label = FALSE, 
        pt.size = .25, 
        group.by = "batch",
        raster = F
        )
hirestiff(paste0(prefixPC,"_UMAP_by_batch_hires.tiff"))
lowrestiff(paste0(prefixPC,"_UMAP_by_batch_lowres.tiff"))

DimPlot(seurat_sct, 
        reduction = "umap", 
        label = FALSE, 
        pt.size = .25, 
        group.by = "batch", 
        split.by = "batch",
        raster = F
        )
hirestiff(paste0(prefixPC,"_UMAP_split_by_batch_hires.tiff"))
lowrestiff(paste0(prefixPC,"_UMAP_split_by_batch_lowres.tiff"))

DimPlot(seurat_sct, 
        reduction = "umap", 
        label = FALSE, 
        pt.size = .25, 
        group.by = "time_point",
        raster = F
        )
hirestiff(paste0(prefixPC,"_UMAP_by_time_point_hires.tiff"))
lowrestiff(paste0(prefixPC,"_UMAP_by_time_point_lowres.tiff"))

DimPlot(seurat_sct, 
        reduction = "umap", 
        label = FALSE, 
        pt.size = .25, 
        split.by = "time_point", 
        group.by = "time_point",
        raster = F
        )
hirestiff(paste0(prefixPC,"_UMAP_split_by_time_point_hires.tiff"))
lowrestiff(paste0(prefixPC,"_UMAP_split_by_time_point_lowres.tiff"))

# saveRDS(seurat_sct, file = paste0(prefixPC,"_seurat_integrated_preclustering.rds"))

# seurat_sct <- read_rds(file = paste0(prefixPC,"_seurat_integrated_preclustering.rds"))
#### (skip to integration) ####
  #### SCT Clustering and Resolution ####
# DefaultAssay(seurat_sct) <- "SCT"

# Determine the K-nearest neighbor graph
seurat_sct <- FindNeighbors(object = seurat_sct, reduction = "pca", dims = 1:meaningful.PCs)

# Determine the clusters                              

seurat_sct <- FindClusters(object = seurat_sct,
                           resolution = c(0.1,0.2,0.3,0.4,0.5))

res <- "_res.0.1"
ident_res <- paste0("SCT_snn",res)

Idents(seurat_sct) <- ident_res
DimPlot(seurat_sct, reduction = "umap", label = TRUE, label.size = 5, pt.size = 0.8) + 
  NoLegend() + 
  ggtitle(paste0(ident_res))

#update prefix
prefixPCres <- paste0(prefixPC,res)

# Plot the UMAP
DimPlot(seurat_sct, reduction = "umap", label = TRUE, label.size = 5, pt.size = 0.8)  + 
  NoLegend() + 
  ggtitle(paste0(ident_res))
hirestiff(paste0(prefixPCres,"_UMAP","_by_","cluster","_hires.tiff"))
lowrestiff(paste0(prefixPCres,"_UMAP","_by_","cluster","_lowres.tiff"))

# UMAP of cells in each cluster by treatment without cluster labels
DimPlot(seurat_sct, reduction = "umap", label = FALSE, split.by = "treatment", pt.size = 0.8)  + NoLegend()
hirestiff(paste0(prefixPCres,"_UMAP","_by_","treatment","_no_labels","_hires.tiff"))
lowrestiff(paste0(prefixPCres,"_UMAP","_by_","treatment","_no_labels","_lowres.tiff"))

# UMAP of cells in each cluster by treatment with cluster labels
DimPlot(seurat_sct, reduction = "umap", label = TRUE, 
        split.by = "treatment", pt.size = 0.8, label.size = 5)  + NoLegend()
hirestiff(paste0(prefixPCres,"_UMAP","_by_","treatment","_with_clusters_labels","_hires.tiff"))
lowrestiff(paste0(prefixPCres,"_UMAP","_by_","treatment","_with_clusters_labels","_lowres.tiff"))

# UMAP of cells in each cluster by batch without cluster labels
DimPlot(seurat_sct, reduction = "umap", label = FALSE, split.by = "batch", pt.size = 0.8)  + NoLegend()
hirestiff(paste0(prefixPCres,"_UMAP","_by_","batch","_no_labels","_hires.tiff"))
lowrestiff(paste0(prefixPCres,"_UMAP","_by_","batch","_no_labels","_lowres.tiff"))

# UMAP of cells in each cluster by batch with cluster labels
DimPlot(seurat_sct, reduction = "umap", label = TRUE, 
        split.by = "batch", pt.size = 0.8, label.size = 5)  + NoLegend()
hirestiff(paste0(prefixPCres,"_UMAP","_by_","batch","_with_clusters_labels","_hires.tiff"))
lowrestiff(paste0(prefixPCres,"_UMAP","_by_","batch","_with_clusters_labels","_lowres.tiff"))

# UMAP of cells in each cluster by time_point without cluster labels
DimPlot(seurat_sct, reduction = "umap", label = FALSE, split.by = "time_point", pt.size = 0.8)  + NoLegend()
hirestiff(paste0(prefixPCres,"_UMAP","_by_","time_point","_no_labels","_hires.tiff"))
lowrestiff(paste0(prefixPCres,"_UMAP","_by_","time_point","_no_labels","_lowres.tiff"))

# UMAP of cells in each cluster by time_point with cluster labels
DimPlot(seurat_sct, reduction = "umap", label = TRUE, 
        split.by = "time_point", pt.size = 0.8, label.size = 5)  + NoLegend()
hirestiff(paste0(prefixPCres,"_UMAP","_by_","time_point","_with_clusters_labels","_hires.tiff"))
lowrestiff(paste0(prefixPCres,"_UMAP","_by_","time_point","_with_clusters_labels","_lowres.tiff"))

  #### Extract number of cells per cluster per orig.ident ####
n_cells <- FetchData(seurat_sct, vars = c("ident", "orig.ident")) %>%
  dplyr::count(ident, orig.ident) %>%
  tidyr::spread(ident, n)

write.csv(n_cells, file = paste0(prefixPCres,"_cells_per_cluster.csv"))

rm(n_cells)

  #### save RDS containing reduction and cluster idents ####
saveRDS(seurat_sct, paste0(prefixPCres,"_seurat_after_clustering.rds"))

# seurat_sct <- readRDS(file = paste0(prefixPCres,"_seurat_after_clustering.rds"))

  #### normalize rna slot ####

# Select the RNA counts slot to be the default assay for visualization purposes
DefaultAssay(seurat_sct) <- "RNA"

# Normalize, find variable features, scale data 
seurat_sct <- NormalizeData(seurat_sct)
seurat_sct <- FindVariableFeatures(seurat_sct)
all.genes <- rownames(seurat_sct)
seurat_sct <- ScaleData(seurat_sct, features = all.genes)

# save object containing RNA normalized data
saveRDS(seurat_sct, paste0(prefixPCres,"_seurat_after_RNAnorm.rds"))

seurat_sct <- readRDS(file = paste0(prefixPCres,"_seurat_after_RNAnorm.rds"))

# SaveH5Seurat(seurat_sct, 
#              filename = paste0(prefixPCres,"_seurat_after_RNAnorm.h5Seurat"))
# 
# seurat_sct_test <- LoadH5Seurat(file = paste0(prefixPCres,"_seurat_after_RNAnorm.h5Seurat"))

#### ####
  #### Run integration (RPCA method) ####
DefaultAssay(seurat_sct) <- "RNA"

s.list <- SplitObject(seurat_sct, split.by = "orig.ident")

# normalize data and find variable features
s.list <- lapply(X = s.list, FUN = function(x) {
  x <- NormalizeData(x, verbose = FALSE)
  x <- FindVariableFeatures(x, verbose = FALSE)
})

# Next, select features for downstream integration, and run PCA on each 
# object in the list, which is required for running the alternative 
# reciprocal PCA workflow.

features <- SelectIntegrationFeatures(object.list = s.list)

s.list <- lapply(X = s.list, FUN = function(x) {
  x <- ScaleData(x, features = features, verbose = FALSE)
  x <- RunPCA(x, features = features, verbose = FALSE)
})

anchors <- FindIntegrationAnchors(object.list = s.list, 
                                  reduction = "rpca",
                                  k.anchor = 10) # default is 5

s.integrated <- IntegrateData(anchorset = anchors, dims = 1:50)
s.integrated <- ScaleData(s.integrated, verbose = FALSE)

# saveRDS(s.integrated, file = paste0("./",date,"_",project,"_seurat_integrated.rds"))

# s.integrated <- read_rds(file = paste0("./",date,"_",project,"_seurat_integrated.rds"))

seurat_sct <- s.integrated

rm(s.list,anchors,s.integrated)

  #### set factor levels ####
seurat_sct$orig.ident <- factor(
  x = seurat_sct$orig.ident,
  levels = c(
    "GEO27.WT.UT",
    "GEO28.WT.UT",
    "08182.p90.TR",
    "13933.p60.TR",
    "13933.p60.UT",
    "16573.p90.TR",
    "16573.p90.UT",
    "18482.p60.TR",
    "18482.p60.UT",
    "19328.p60.TR",
    "19328.p60.UT"
  )
)

seurat_sct$treatment <- factor(
  x = seurat_sct$treatment,
  levels = c(
    "Untreated",
    "Cell_Treated"
  )
)

seurat_sct$batch <- factor(
  x = seurat_sct$batch,
  levels = c(
    "GSE209872",
    "08182",
    "13933",
    "16573",
    "18482",
    "19328"
  )
)

seurat_sct$time_point <- factor(
  x = seurat_sct$time_point,
  levels = c(
    "WT",
    "p60",
    "p90"
  )
)

seurat_sct$time_point_treatment <- factor(
  x = seurat_sct$time_point_treatment,
  levels = c(
    "WT_Untreated",
    "p60_Cell_Treated",
    "p60_Untreated",
    "p90_Cell_Treated",
    "p90_Untreated"
  )
)

seurat_sct$time_point_batch <- factor(
  x = seurat_sct$time_point_batch,
  levels = c(
    "WT_GSE209872",
    "p60_13933",
    "p60_18482",
    "p60_19328",
    "p90_08182",
    "p90_16573"
  )
)

seurat_sct$treatment_batch <- factor(
  x = seurat_sct$treatment_batch,
  levels = c(
    "Untreated_GSE209872",
    "Untreated_13933",
    "Untreated_18482",
    "Untreated_19328",
    "Untreated_16573",
    "Cell_Treated_13933",
    "Cell_Treated_18482",
    "Cell_Treated_19328",
    "Cell_Treated_08182",
    "Cell_Treated_16573"
  )
)

seurat_sct$time_point_treatment_batch <- factor(
  x = seurat_sct$time_point_treatment_batch,
  levels = c(
    "WT_Untreated_GSE209872",
    "p60_Untreated_13933",
    "p60_Untreated_18482",
    "p60_Untreated_19328",
    "p90_Untreated_16573",
    "p60_Cell_Treated_13933",
    "p60_Cell_Treated_18482",
    "p60_Cell_Treated_19328",
    "p90_Cell_Treated_08182",
    "p90_Cell_Treated_16573"
  )
)

  #### update project variable to include integration ####
project <- paste0(project,"_int")

  #### INT PCA ####
# this performs PCA on the seurat object
seurat_sct <- RunPCA(seurat_sct, npcs = 50, verbose = TRUE)

# make PC coordinate object a data frame
xx.coord <- as.data.frame(seurat_sct@reductions$pca@cell.embeddings)

# make PC feature loadings object a data frame
xx.gload <- as.data.frame(seurat_sct@reductions$pca@feature.loadings)

# calculate eigenvalues for arrays

# generate squares of all sample coordinates
sq.xx.coord <- as.data.frame(xx.coord^2)
# create empty list for eigenvalues first
eig <- c()
# calculate the eigenvalue for each PC in sq.xx.coord by taking the sqrt of the sum of squares
for(i in 1:ncol(sq.xx.coord))
  eig[i] = sqrt(sum(sq.xx.coord[,i]))
# calculate the total variance by adding up all the eigenvalues
sum.eig <- sum(eig)
# calculate the expected contribution of all PCs if they all contribute equally to the total variance
expected.contribution <- sum.eig/(length(xx.coord)-1)
# return the number of principal components with an eigenvalue greater than expected by equal variance
meaningful.PCs <- sum(eig > expected.contribution)

# create empty list for eigenvalue percentage
eig.percent <- c()
# calculate the percentage of the total variance by each PC eigenvalue
for(i in 1:length(eig))
  eig.percent[i] = 100*eig[i]/sum.eig
# sum of all eig.percent should total to 100
sum(eig.percent)
# create empty list for scree values
scree <- c()
# calculate a running total of variance contribution
for(i in 1:length(eig))
  if(i == 1) scree[i] = eig.percent[i] else scree[i] = scree[i-1] + eig.percent[i]

# create data frame for eigenvalue summaries
eigenvalues <- data.frame("PC" = colnames(xx.coord), "eig" = eig, "percent" = eig.percent, "scree" = scree)

# write csv for eigenvalues
# write.csv(eigenvalues, file = paste0("./",date,"_",project,"_PCA_eigenvalues.csv"), row.names = F)

# plot scree values
plot(eigenvalues$percent, ylim = c(0,100), type = "S", xlab = "PC", ylab = "Percent of variance",
     main = paste0(date,"_",project," scree plot all samples PCA"))
points(eigenvalues$scree, ylim = c(0,100), type = "p", pch = 16)
lines(eigenvalues$scree)
# add red line to indicate cut-off
cut.off <- 100/(length(eig)-1)
abline(h = cut.off, col = "red")
# add blue line to indicate which PCs are meaningful and kept
abline(v = meaningful.PCs, col = "blue")
text(meaningful.PCs, cut.off, label = paste("cutoff PC",meaningful.PCs),
     adj = c(-0.1, -0.5))

dev.copy(pdf, paste0("./",date,"_",project,"_scree_plot.pdf"))
dev.off()

rm(eigenvalues,sq.xx.coord,xx.coord,xx.gload,cut.off,
   eig,eig.percent,expected.contribution,i,scree,sum.eig)

# meaningful.PCs <- 12

  #### INT Run UMAP and look at UMAP plots ####

seurat_sct <- RunUMAP(seurat_sct, 
                      reduction = "pca", 
                      dims = 1:meaningful.PCs, 
                      verbose = TRUE)

# update prefixed variable
prefixPC <- paste0("./",date,"_",project,"_",meaningful.PCs,"PCs")

## UMAP plot by sample name ("orig.ident")
DimPlot(seurat_sct, 
        reduction = "umap", 
        label = FALSE, 
        pt.size = .25, 
        group.by = "orig.ident",
        raster = F
        )
hirestiff(paste0(prefixPC,"_UMAP_by_sample_hires.tiff"))
lowrestiff(paste0(prefixPC,"_UMAP_by_sample_lowres.tiff"))

DimPlot(seurat_sct, 
        reduction = "umap", 
        label = FALSE, 
        pt.size = .25, 
        group.by = "treatment",
        raster = F
        )
hirestiff(paste0(prefixPC,"_UMAP_by_treatment_hires.tiff"))
lowrestiff(paste0(prefixPC,"_UMAP_by_treatment_lowres.tiff"))

DimPlot(seurat_sct, 
        reduction = "umap", 
        label = FALSE, 
        pt.size = .25, 
        split.by = "treatment", 
        group.by = "treatment",
        raster = F
        )
hirestiff(paste0(prefixPC,"_UMAP_split_by_treatment_hires.tiff"))
lowrestiff(paste0(prefixPC,"_UMAP_split_by_treatment_lowres.tiff"))

DimPlot(seurat_sct, 
        reduction = "umap", 
        label = FALSE, 
        pt.size = .25, 
        group.by = "batch",
        raster = F
        )
hirestiff(paste0(prefixPC,"_UMAP_by_batch_hires.tiff"))
lowrestiff(paste0(prefixPC,"_UMAP_by_batch_lowres.tiff"))

DimPlot(seurat_sct, 
        reduction = "umap", 
        label = FALSE, 
        pt.size = .25, 
        group.by = "batch", 
        split.by = "batch",
        raster = F
        )
hirestiff(paste0(prefixPC,"_UMAP_split_by_batch_hires.tiff"))
lowrestiff(paste0(prefixPC,"_UMAP_split_by_batch_lowres.tiff"))

DimPlot(seurat_sct, 
        reduction = "umap", 
        label = FALSE, 
        pt.size = .25, 
        group.by = "time_point",
        raster = F
        )
hirestiff(paste0(prefixPC,"_UMAP_by_time_point_hires.tiff"))
lowrestiff(paste0(prefixPC,"_UMAP_by_time_point_lowres.tiff"))

DimPlot(seurat_sct, 
        reduction = "umap", 
        label = FALSE, 
        pt.size = .25, 
        split.by = "time_point", 
        group.by = "time_point",
        raster = F
        )
hirestiff(paste0(prefixPC,"_UMAP_split_by_time_point_hires.tiff"))
lowrestiff(paste0(prefixPC,"_UMAP_split_by_time_point_lowres.tiff"))

DimPlot(seurat_sct, 
        reduction = "umap", 
        label = FALSE, 
        pt.size = .25, 
        group.by = "time_point_treatment",
        raster = F
        )
hirestiff(paste0(prefixPC,"_UMAP_by_time_point_treatment_hires.tiff"))
lowrestiff(paste0(prefixPC,"_UMAP_by_time_point_treatment_lowres.tiff"))

DimPlot(seurat_sct, 
        reduction = "umap", 
        label = FALSE, 
        pt.size = .25, 
        split.by = "time_point_treatment", 
        group.by = "time_point_treatment",
        raster = F
        )
hirestiff(paste0(prefixPC,"_UMAP_split_by_time_point_treatment_hires.tiff"))
lowrestiff(paste0(prefixPC,"_UMAP_split_by_time_point_treatment_lowres.tiff"))

DimPlot(seurat_sct, 
        reduction = "umap", 
        label = FALSE, 
        pt.size = .25, 
        group.by = "time_point_batch",
        raster = F
        )
hirestiff(paste0(prefixPC,"_UMAP_by_time_point_batch_hires.tiff"))
lowrestiff(paste0(prefixPC,"_UMAP_by_time_point_batch_lowres.tiff"))

DimPlot(seurat_sct, 
        reduction = "umap", 
        label = FALSE, 
        pt.size = .25, 
        split.by = "time_point_batch", 
        group.by = "time_point_batch",
        raster = F
        )
hirestiff(paste0(prefixPC,"_UMAP_split_by_time_point_batch_hires.tiff"))
lowrestiff(paste0(prefixPC,"_UMAP_split_by_time_point_batch_lowres.tiff"))

DimPlot(seurat_sct, 
        reduction = "umap", 
        label = FALSE, 
        pt.size = .25, 
        group.by = "treatment_batch",
        raster = F
        )
hirestiff(paste0(prefixPC,"_UMAP_by_treatment_batch_hires.tiff"))
lowrestiff(paste0(prefixPC,"_UMAP_by_treatment_batch_lowres.tiff"))

DimPlot(seurat_sct, 
        reduction = "umap", 
        label = FALSE, 
        pt.size = .25, 
        split.by = "treatment_batch", 
        group.by = "treatment_batch",
        raster = F
        )
hirestiff(paste0(prefixPC,"_UMAP_split_by_treatment_batch_hires.tiff"))
lowrestiff(paste0(prefixPC,"_UMAP_split_by_treatment_batch_lowres.tiff"))

DimPlot(seurat_sct, 
        reduction = "umap", 
        label = FALSE, 
        pt.size = .25, 
        group.by = "time_point_treatment_batch",
        raster = F
        )
hirestiff(paste0(prefixPC,"_UMAP_by_time_point_treatment_batch_hires.tiff"))
lowrestiff(paste0(prefixPC,"_UMAP_by_time_point_treatment_batch_lowres.tiff"))

DimPlot(seurat_sct, 
        reduction = "umap", 
        label = FALSE, 
        pt.size = .25, 
        split.by = "time_point_treatment_batch", 
        group.by = "time_point_treatment_batch",
        raster = F
        )
hirestiff(paste0(prefixPC,"_UMAP_split_by_time_point_treatment_batch_hires.tiff"))
lowrestiff(paste0(prefixPC,"_UMAP_split_by_time_point_treatment_batch_lowres.tiff"))

# saveRDS(seurat_sct, file = paste0(prefixPC,"_seurat_integrated_preclustering.rds"))

# seurat_sct <- read_rds(file = paste0(prefixPC,"_seurat_integrated_preclustering.rds"))

  #### INT Clustering and Resolution ####
# DefaultAssay(seurat_sct) <- "integrated"

# Determine the K-nearest neighbor graph
seurat_sct <- FindNeighbors(object = seurat_sct, 
                            reduction = "pca", 
                            dims = 1:meaningful.PCs)

# Determine the clusters                              
# seurat_sct <- FindClusters(object = seurat_sct,
                           # resolution = c(0.1,0.2,0.3,0.4,0.5))

seurat_sct <- FindClusters(object = seurat_sct,
                           resolution = c(0.5))

res <- "_res.0.5"
ident_res <- paste0("integrated_snn",res)

Idents(seurat_sct) <- ident_res

DimPlot(seurat_sct, 
        reduction = "umap", 
        label = TRUE, 
        label.size = 5, 
        pt.size = 0.8,
        raster = F
        ) + 
  NoLegend() + 
  ggtitle(paste0(ident_res))

# set cluster colors with Palo (makes sure adjacent colors are different)
palopal <- Palo::Palo(seurat_sct[["umap"]]@cell.embeddings,
                      as.character(Idents(seurat_sct)),
                      scales::hue_pal()(length(levels(seurat_sct))))

# sort these so they are in numeric order and can be transferred to new idents
sorted_palopal <- palopal[order(as.numeric(names(palopal)))]

# test resolution with Palo colors
DimPlot(seurat_sct, 
        reduction = "umap", 
        label = TRUE, 
        label.size = 5, 
        pt.size = 0.8,
        cols = sorted_palopal,
        raster = F
        ) + 
  NoLegend() + 
  ggtitle(paste0(ident_res))

#update prefix
prefixPCres <- paste0(prefixPC,res)

# Plot the UMAP
DimPlot(seurat_sct, 
        reduction = "umap", 
        label = TRUE, 
        label.size = 5, 
        pt.size = 0.8,
        cols = sorted_palopal,
        raster = F
        ) + 
  NoLegend() + 
  ggtitle(paste0(ident_res))
hirestiff(paste0(prefixPCres,"_UMAP","_by_","cluster","_hires.tiff"))
lowrestiff(paste0(prefixPCres,"_UMAP","_by_","cluster","_lowres.tiff"))

# UMAP of cells in each cluster by treatment without cluster labels
DimPlot(seurat_sct, 
        reduction = "umap", 
        label = FALSE, 
        split.by = "treatment", 
        pt.size = 0.8,
        cols = sorted_palopal,
        raster = F
        ) + 
  NoLegend()
hirestiff(paste0(prefixPCres,"_UMAP","_by_","treatment",
                 "_no_labels","_hires.tiff"))
lowrestiff(paste0(prefixPCres,"_UMAP","_by_","treatment",
                  "_no_labels","_lowres.tiff"))

# UMAP of cells in each cluster by treatment with cluster labels
DimPlot(seurat_sct, 
        reduction = "umap", 
        label = TRUE, 
        split.by = "treatment", 
        pt.size = 0.8, 
        label.size = 5,
        cols = sorted_palopal,
        raster = F
        ) + 
  NoLegend()
hirestiff(paste0(prefixPCres,"_UMAP","_by_","treatment",
                 "_with_clusters_labels","_hires.tiff"))
lowrestiff(paste0(prefixPCres,"_UMAP","_by_","treatment",
                  "_with_clusters_labels","_lowres.tiff"))

# UMAP of cells in each cluster by batch without cluster labels
DimPlot(seurat_sct, 
        reduction = "umap", 
        label = FALSE, 
        split.by = "batch", 
        pt.size = 0.8,
        cols = sorted_palopal,
        raster = F
        ) + 
  NoLegend()
hirestiff(paste0(prefixPCres,"_UMAP","_by_","batch",
                 "_no_labels","_hires.tiff"))
lowrestiff(paste0(prefixPCres,"_UMAP","_by_","batch",
                  "_no_labels","_lowres.tiff"))

# UMAP of cells in each cluster by batch with cluster labels
DimPlot(seurat_sct, 
        reduction = "umap", 
        label = TRUE, 
        split.by = "batch", 
        pt.size = 0.8, 
        label.size = 3,
        cols = sorted_palopal,
        raster = F
        ) + 
  NoLegend()
hirestiff(paste0(prefixPCres,"_UMAP","_by_","batch",
                 "_with_clusters_labels","_hires.tiff"))
lowrestiff(paste0(prefixPCres,"_UMAP","_by_","batch",
                  "_with_clusters_labels","_lowres.tiff"))

# UMAP of cells in each cluster by time_point without cluster labels
DimPlot(seurat_sct, 
        reduction = "umap", 
        label = FALSE, 
        split.by = "time_point", 
        pt.size = 0.8,
        cols = sorted_palopal,
        raster = F
        ) + 
  NoLegend()
hirestiff(paste0(prefixPCres,"_UMAP","_by_","time_point",
                 "_no_labels","_hires.tiff"))
lowrestiff(paste0(prefixPCres,"_UMAP","_by_","time_point",
                  "_no_labels","_lowres.tiff"))

# UMAP of cells in each cluster by time_point with cluster labels
DimPlot(seurat_sct, 
        reduction = "umap", 
        label = TRUE, 
        split.by = "time_point", 
        pt.size = 0.8, 
        label.size = 5,
        cols = sorted_palopal,
        raster = F
        ) + 
  NoLegend()
hirestiff(paste0(prefixPCres,"_UMAP","_by_","time_point",
                 "_with_clusters_labels","_hires.tiff"))
lowrestiff(paste0(prefixPCres,"_UMAP","_by_","time_point",
                  "_with_clusters_labels","_lowres.tiff"))

# UMAP of cells in each cluster by time_point without cluster labels
DimPlot(seurat_sct, 
        reduction = "umap", 
        label = FALSE, 
        split.by = "time_point_treatment", 
        pt.size = 0.8,
        cols = sorted_palopal,
        raster = F
        ) + 
  NoLegend()
hirestiff(paste0(prefixPCres,"_UMAP","_by_","time_point_treatment",
                 "_no_labels","_hires.tiff"))
lowrestiff(paste0(prefixPCres,"_UMAP","_by_","time_point_treatment",
                  "_no_labels","_lowres.tiff"))

# UMAP of cells in each cluster by time_point with cluster labels
DimPlot(seurat_sct, 
        reduction = "umap", 
        label = TRUE, 
        split.by = "time_point_treatment", 
        pt.size = 0.8, 
        label.size = 5,
        cols = sorted_palopal,
        raster = F
        ) + 
  NoLegend()
hirestiff(paste0(prefixPCres,"_UMAP","_by_","time_point_treatment",
                 "_with_clusters_labels","_hires.tiff"))
lowrestiff(paste0(prefixPCres,"_UMAP","_by_","time_point_treatment",
                  "_with_clusters_labels","_lowres.tiff"))

  #### INT Extract number of cells per cluster per orig.ident ####
n_cells <- FetchData(seurat_sct, vars = c("ident", "orig.ident")) %>%
  dplyr::count(ident, orig.ident) %>%
  tidyr::spread(ident, n)

write.csv(n_cells, file = paste0(prefixPCres,"_cells_per_cluster.csv"))

rm(n_cells)

  #### INT save RDS containing reduction and cluster idents ####
# saveRDS(seurat_sct, paste0(prefixPCres,"_seurat_after_clustering.rds"))

# seurat_sct <- readRDS(file = paste0(prefixPCres,"_seurat_after_clustering.rds"))

  #### INT normalize rna slot ####

# Select the RNA counts slot to be the default assay for visualization purposes
DefaultAssay(seurat_sct) <- "RNA"

# Normalize, find variable features, scale data 
seurat_sct <- NormalizeData(seurat_sct)
seurat_sct <- FindVariableFeatures(seurat_sct)
all.genes <- rownames(seurat_sct)
seurat_sct <- ScaleData(seurat_sct, features = all.genes)

# Export normalized counts
norm_counts <- GetAssayData(seurat_sct, slot = "data")

save(norm_counts, file = paste0(prefixPCres,"_normalized_counts_RNA_sparse_matrix.rdata"))

rm(norm_counts)

# save RDS containing RNA normalized data
# saveRDS(seurat_sct, paste0(prefixPCres,"_seurat_after_RNAnorm.rds"))
# seurat_sct <- readRDS(file = paste0(prefixPCres,"_seurat_after_RNAnorm.rds"))

# save RDS containing RNA normalized data AND CELLTYPE LABELS
# saveRDS(seurat_sct, paste0(prefixPCres,"_seurat_after_RNAnorm_CELLTYPES.rds"))
# seurat_sct <- readRDS(file = paste0(prefixPCres,"_seurat_after_RNAnorm_CELLTYPES.rds"))

  #### INT Cluster QC #### 

# Look at QC metrics for clustering by cell quality

VlnPlot(seurat_sct,
        features = "nCount_RNA",
        raster = F) + 
  NoLegend()

hirestiff(paste(prefixPCres,"vln","plot","QC","nCount_RNA","hires.tiff", sep = "_"))
lowrestiff(paste(prefixPCres,"vln","plot","QC","nCount_RNA","lowres.tiff", sep = "_"))

# Low count clusters:  nothing obvious

VlnPlot(seurat_sct,
        features = "nFeature_RNA",
        raster = F) + 
  NoLegend()

hirestiff(paste(prefixPCres,"vln","plot","QC","nFeature_RNA","hires.tiff", sep = "_"))
lowrestiff(paste(prefixPCres,"vln","plot","QC","nFeature_RNA","lowres.tiff", sep = "_"))

# Low feature clusters:  nothing obvious

VlnPlot(seurat_sct,
        features = "percent.mt",
        raster = F) + 
  NoLegend()

hirestiff(paste(prefixPCres,"vln","plot","QC","percent.mt","hires.tiff", sep = "_"))
lowrestiff(paste(prefixPCres,"vln","plot","QC","percent.mt","lowres.tiff", sep = "_"))

# High percent mito clusters:  various

VlnPlot(seurat_sct,
        features = "percent.ribo",
        raster = F) + 
  NoLegend()

hirestiff(paste(prefixPCres,"vln","plot","QC","percent.ribo","hires.tiff", sep = "_"))
lowrestiff(paste(prefixPCres,"vln","plot","QC","percent.ribo","lowres.tiff", sep = "_"))

# High percent ribo clusters:  varies

#### PLOTS ####
  #### Test Plots ####
# apoptosis genes (bax and casp3 up, other two down)
FeaturePlot(seurat_sct, 
            reduction = "umap", 
            features = c("Bax","Casp3","Bnip3","Bcl2"), 
            order = TRUE,
            min.cutoff = 'q10',
            label = TRUE,
            label.size = 3,
            pt.size = 0.5,
            raster = F)

# Plot the UMAP
DimPlot(seurat_sct, 
        reduction = "umap", 
        label = TRUE, 
        label.size = 5, 
        pt.size = 0.8) + 
  NoLegend()

# Plot the UMAP split by Treatment
DimPlot(seurat_sct, 
        reduction = "umap", 
        label = TRUE, 
        split.by = "treatment", 
        pt.size = 0.8) + 
  NoLegend()


Rods <- c("Rho","Gnat1","Nrl","Pde6a","Pde6b","Cngb1")
Rods1 <- c("Rho","Gnat1","Nrl")
Rods2 <- c("Pde6a","Pde6b","Cngb1")
Cones <- c("Arr3","Opn1sw","Gnat2","Pde6h")
Mueller_Glia <- c("Glul","Rlbp1")
Astrocytes <- c("Gfap","Aqp4","Gpr37")
Microglia <- c("Tyrobp","Cx3cr1","Tmem119","Clec7a")
Bipolar <- c("Scgn","Vsx1","Vsx2")
Retinal_Ganglion <- c("Rbpms","Pou4f1","Pou4f2")
Amacrine <- c("Calb1","Gad1","Gad2")
Horizontal <- c("Onecut2")
Pericytes <- c("Pecam1","Acta2")
Reticulocytes <- c("Hbb-bs")

default_FP <- function(obj,features){
  FeaturePlot(obj, 
              reduction = "umap", 
              features = features, 
              order = TRUE,
              # min.cutoff = 'q10',
              label = TRUE,
              label.size = 3,
              pt.size = 0.5)
}

  #### Genes from Dr. Wang ####

FeaturePlot(seurat_sct, 
            reduction = "umap", 
            features = c("Ifitm3","Socs3","Slpi","Ch25h"), 
            order = TRUE,
            min.cutoff = 'q10',
            label = TRUE,
            label.size = 3,
            pt.size = 0.5)

FeaturePlot(seurat_sct, 
            reduction = "umap", 
            features = c("Rlbp1", # Rlbp1 = Crlbp1
                         "Anxa1","Anxa2","Nfkbia"), 
            order = TRUE,
            min.cutoff = 'q10',
            label = TRUE,
            label.size = 3,
            pt.size = 0.5)

FeaturePlot(seurat_sct, 
            reduction = "umap", 
            features = c("Aif1", # Aif1 = Iba1
                         "Adgre1", # Adgre1 = F4/80
                         "Cd68"), 
            order = TRUE,
            min.cutoff = 'q10',
            label = TRUE,
            label.size = 3,
            pt.size = 0.5)

FeaturePlot(seurat_sct, 
            reduction = "umap", 
            features = c("Aifm1","Aif1"), 
            order = TRUE,
            min.cutoff = 'q10',
            label = TRUE,
            label.size = 3,
            pt.size = 0.5)

FeaturePlot(seurat_sct, 
            reduction = "umap", 
            features = c("Mydgf","Kdelr1","Kdelr2","Ern1"), 
            order = TRUE,
            min.cutoff = 'q10',
            label = TRUE,
            label.size = 3,
            pt.size = 0.5)

FeaturePlot(seurat_sct, 
            reduction = "umap", 
            features = c("Mydgf"), 
            order = TRUE,
            min.cutoff = 'q10',
            label = TRUE,
            label.size = 3,
            pt.size = 0.5,
            split.by = "time_point_treatment")

FeaturePlot(seurat_sct, 
            reduction = "umap", 
            features = c("Kdelr1"), 
            order = TRUE,
            min.cutoff = 'q10',
            label = TRUE,
            label.size = 3,
            pt.size = 0.5,
            split.by = "time_point_treatment")

FeaturePlot(seurat_sct, 
            reduction = "umap", 
            features = c("Kdelr2"), 
            order = TRUE,
            min.cutoff = 'q10',
            label = TRUE,
            label.size = 3,
            pt.size = 0.5,
            split.by = "time_point_treatment")

FeaturePlot(seurat_sct, 
            reduction = "umap", 
            features = c("Ern1"), 
            order = TRUE,
            min.cutoff = 'q10',
            label = TRUE,
            label.size = 3,
            pt.size = 0.5,
            split.by = "time_point_treatment")

FeaturePlot(seurat_sct, 
            reduction = "umap", 
            features = c("Manf"), 
            order = TRUE,
            min.cutoff = 'q10',
            label = TRUE,
            label.size = 3,
            pt.size = 0.5,
            split.by = "time_point_treatment")

VlnPlot(seurat_sct,
        idents = c("ROD","CONE"),
        features = c("Mydgf","Kdelr1","Kdelr2","Ern1","Manf"),
        split.by = "time_point_treatment",
        stack = T,
        flip = T)

FeaturePlot(seurat_sct, 
            reduction = "umap", 
            features = c("Creb1", 
                         "Crem"), 
            order = TRUE,
            min.cutoff = 'q10',
            label = TRUE,
            label.size = 3,
            pt.size = 0.5)

FeaturePlot(seurat_sct, 
            reduction = "umap", 
            features = c("Tmem119","P2ry13","P2ry12","Siglech"), 
            order = TRUE,
            min.cutoff = 'q10',
            label = TRUE,
            label.size = 3,
            pt.size = 0.5)

  #### retina genes from Saba SCT ####

# Rod Photoreceptors
# Cluster 7
VlnPlot(seurat_sct,
        features = c("Rho","Gnat1","Nrl",
                     "Pde6a","Pde6b","Cngb1"),
        stack = TRUE,
        flip = TRUE) +
  NoLegend() +
  ggtitle("Rod Photoreceptors")

hirestiff(paste(prefixPCres,"VLN","Rod","Photoreceptors","hires.tiff", sep = "_"))
lowrestiff(paste(prefixPCres,"VLN","Rod","Photoreceptors","lowres.tiff", sep = "_"))

# Cone Photoreceptors
# Cluster 8
VlnPlot(seurat_sct,
        features = c("Arr3","Opn1sw","Gnat2","Pde6h"),
        stack = TRUE,
        flip = TRUE) +
  NoLegend() +
  ggtitle("Cone Photoreceptors")

hirestiff(paste(prefixPCres,"VLN","Cone","Photoreceptors","hires.tiff", sep = "_"))
lowrestiff(paste(prefixPCres,"VLN","Cone","Photoreceptors","lowres.tiff", sep = "_"))

# Muller Glia
# Clusters 4, 10, 12 
VlnPlot(seurat_sct,
        features = c("Glul","Rlbp1","Nes","Vim"),
        stack = TRUE,
        flip = TRUE) +
  NoLegend() +
  ggtitle("Muller Glia")

hirestiff(paste(prefixPCres,"VLN","Muller","Glia","hires.tiff", sep = "_"))
lowrestiff(paste(prefixPCres,"VLN","Muller","Glia","lowres.tiff", sep = "_"))

# Retinal Astrocytes
# Clusters 4, 10, 12 
VlnPlot(seurat_sct,
        features = c("Gfap","Aqp4","Gpr37"),
        stack = TRUE,
        flip = TRUE) +
  NoLegend() +
  ggtitle("Retinal Astrocytes")

hirestiff(paste(prefixPCres,"VLN","Retinal","Astrocytes","hires.tiff", sep = "_"))
lowrestiff(paste(prefixPCres,"VLN","Retinal","Astrocytes","lowres.tiff", sep = "_"))

# Microglia 
# Clusters 3, 13, 16
VlnPlot(seurat_sct,
        features = c("Tyrobp","Cx3cr1","Tmem119"),
        stack = TRUE,
        flip = TRUE) +
  NoLegend() +
  ggtitle("Microglia")

hirestiff(paste(prefixPCres,"VLN","Microglia","hires.tiff", sep = "_"))
lowrestiff(paste(prefixPCres,"VLN","Microglia","lowres.tiff", sep = "_"))

# Bipolar Cells 
# Clusters 0, 1, (4), 6, 9, (10), ?11 ,(12,14)
VlnPlot(seurat_sct,
        features = c("Scgn","Vsx1","Vsx2"),
        stack = TRUE,
        flip = TRUE) +
  NoLegend() +
  ggtitle("Bipolar Cells")

hirestiff(paste(prefixPCres,"VLN","Bipolar","Cells","hires.tiff", sep = "_"))
lowrestiff(paste(prefixPCres,"VLN","Bipolar","Cells","lowres.tiff", sep = "_"))

# Retinal Ganglion Cells 
# Clusters 2, 5, 13, 14, 16
VlnPlot(seurat_sct,
        features = c("Rbpms","Pou4f1","Pou4f2"),
        stack = TRUE,
        flip = TRUE) +
  NoLegend() +
  ggtitle("Retinal Ganglion Cells")

hirestiff(paste(prefixPCres,"VLN","Retinal","Ganglion","Cells","hires.tiff", sep = "_"))
lowrestiff(paste(prefixPCres,"VLN","Retinal","Ganglion","Cells","lowres.tiff", sep = "_"))

# Amacrine Cells 
# Part of Cluster 6
VlnPlot(seurat_sct,
        features = c("Calb1","Gad1","Gad2"),
        stack = TRUE,
        flip = TRUE) +
  NoLegend() +
  ggtitle("Amacrine Cells")

hirestiff(paste(prefixPCres,"VLN","Amacrine","Cells","hires.tiff", sep = "_"))
lowrestiff(paste(prefixPCres,"VLN","Amacrine","Cells","lowres.tiff", sep = "_"))

# Horizontal Cells 
# Cluster
VlnPlot(seurat_sct,
        features = c("Onecut2"),
        # stack = TRUE,
        flip = TRUE,
        pt.size = 0) +
  NoLegend() +
  ggtitle("Horizontal Cells")

hirestiff(paste(prefixPCres,"VLN","Horizontal","Cells","hires.tiff", sep = "_"))
lowrestiff(paste(prefixPCres,"VLN","Horizontal","Cells","lowres.tiff", sep = "_"))

# Pericytes
# Clusters 2,16
VlnPlot(seurat_sct,
        features = c(#"Pecam1",
                     "Acta2"),
        # stack = TRUE,
        flip = TRUE,
        pt.size = 0) +
  NoLegend() #+ ggtitle("Pericytes")

hirestiff(paste(prefixPCres,"VLN","Pericytes","hires.tiff", sep = "_"))
lowrestiff(paste(prefixPCres,"VLN","Pericytes","lowres.tiff", sep = "_"))

# Reticulocytes
# Cluster 15
VlnPlot(seurat_sct,
        features = c(#"Pecam1",
          "Hbb-bs"),
        # stack = TRUE,
        flip = TRUE,
        pt.size = 0) +
  NoLegend() #+ ggtitle("Pericytes")

hirestiff(paste(prefixPCres,"VLN","Pericytes","hires.tiff", sep = "_"))
lowrestiff(paste(prefixPCres,"VLN","Pericytes","lowres.tiff", sep = "_"))

  #### retina genes from Saba Integrated ####

VlnPlot(seurat_sct,
        features = c("Nfe2l2"),
        # stack = TRUE,
        flip = TRUE) +
  NoLegend()

default_FP(seurat_sct,c("Nfe2l2"))

VlnPlot(seurat_sct,
        features = c(RPEpi,"Mertk"),
        stack = TRUE,
        flip = TRUE) +
  NoLegend()

# Rod Photoreceptors
# Cluster 8
VlnPlot(seurat_sct,
        features = Rods,
        stack = TRUE,
        flip = TRUE) +
  NoLegend() +
  ggtitle("Rod Photoreceptors")

hirestiff(paste(prefixPCres,"VLN","Rod","Photoreceptors","hires.tiff", sep = "_"))
lowrestiff(paste(prefixPCres,"VLN","Rod","Photoreceptors","lowres.tiff", sep = "_"))

default_FP(seurat_sct,Rods1)

hirestiff(paste(prefixPCres,"UMAP","Rod","genes1","hires.tiff", sep = "_"))
lowrestiff(paste(prefixPCres,"UMAP","Rod","genes1","lowres.tiff", sep = "_"))

default_FP(seurat_sct,Rods2)

hirestiff(paste(prefixPCres,"UMAP","Rod","genes2","hires.tiff", sep = "_"))
lowrestiff(paste(prefixPCres,"UMAP","Rod","genes2","lowres.tiff", sep = "_"))


default_FP(seurat_sct,"Grm6")

default_FP(seurat_sct,c("Pou4f1","Pou4f2","Pou4f3"))

VlnPlot(seurat_sct,
        features = c("Onecut1","Onecut2"),
        stack = TRUE,
        flip = TRUE) +
  NoLegend()

# Cone Photoreceptors
# Cluster 9
VlnPlot(seurat_sct,
        features = Cones,
        stack = TRUE,
        flip = TRUE) +
  NoLegend() +
  ggtitle("Cone Photoreceptors")

hirestiff(paste(prefixPCres,"VLN","Cone","Photoreceptors","hires.tiff", sep = "_"))
lowrestiff(paste(prefixPCres,"VLN","Cone","Photoreceptors","lowres.tiff", sep = "_"))

default_FP(seurat_sct,Cones)

hirestiff(paste(prefixPCres,"UMAP","Cone","genes","hires.tiff", sep = "_"))
lowrestiff(paste(prefixPCres,"UMAP","Cone","genes","lowres.tiff", sep = "_"))

FeaturePlot(seurat_sct, 
            reduction = "umap", 
            features = Cones, 
            order = TRUE,
            min.cutoff = 'q10',
            label = TRUE,
            label.size = 3,
            pt.size = 0.5)

# Muller Glia
# Clusters 10,11,14,15,16,20,24,27,29
VlnPlot(seurat_sct,
        features = c(Mueller_Glia,"Gfap","Aqp4","S100b"),
        stack = TRUE,
        flip = TRUE) +
  NoLegend() +
  ggtitle("Muller Glia")

hirestiff(paste(prefixPCres,"VLN","Muller","Glia","hires.tiff", sep = "_"))
lowrestiff(paste(prefixPCres,"VLN","Muller","Glia","lowres.tiff", sep = "_"))

default_FP(seurat_sct,Mueller_Glia)

hirestiff(paste(prefixPCres,"UMAP","Mueller_Glia","genes","hires.tiff", sep = "_"))
lowrestiff(paste(prefixPCres,"UMAP","Mueller_Glia","genes","lowres.tiff", sep = "_"))

Macroglia <- c("Apoe","Slc1a3","Trpm3","Glul","Clu" )


FeaturePlot(seurat_sct, 
            reduction = "umap", 
            features = Macroglia, 
            order = TRUE,
            min.cutoff = 'q10',
            label = TRUE,
            label.size = 3,
            pt.size = 0.5)

FeaturePlot(seurat_sct, 
            reduction = "umap", 
            features = Mueller_Glia, 
            order = TRUE,
            min.cutoff = 'q10',
            label = TRUE,
            label.size = 3,
            pt.size = 0.5)

FeaturePlot(seurat_sct, 
            reduction = "umap", 
            features = c("Glul"), 
            order = TRUE,
            min.cutoff = 'q10',
            label = TRUE,
            label.size = 3,
            pt.size = 0.5,
            split.by = "time_point")

FeaturePlot(TRseurat, 
            reduction = "umap", 
            features = c("Glul"), 
            order = TRUE,
            min.cutoff = 'q10',
            label = TRUE,
            label.size = 3,
            pt.size = 0.5,
            split.by = "time_point")

FeaturePlot(UTseurat, 
            reduction = "umap", 
            features = c("Glul"), 
            order = TRUE,
            min.cutoff = 'q10',
            label = TRUE,
            label.size = 3,
            pt.size = 0.5,
            split.by = "time_point")

# Retinal Astrocytes
# Clusters 10,11,14,15,16,20,24,27,29
VlnPlot(seurat_sct,
        features = Astrocytes,
        stack = TRUE,
        flip = TRUE) +
  NoLegend() +
  ggtitle("Retinal Astrocytes")

hirestiff(paste(prefixPCres,"VLN","Retinal","Astrocytes","hires.tiff", sep = "_"))
lowrestiff(paste(prefixPCres,"VLN","Retinal","Astrocytes","lowres.tiff", sep = "_"))

default_FP(seurat_sct,Astrocytes)

hirestiff(paste(prefixPCres,"UMAP","Astrocytes","genes","hires.tiff", sep = "_"))
lowrestiff(paste(prefixPCres,"UMAP","Astrocytes","genes","lowres.tiff", sep = "_"))

# Microglia 
# Clusters 5,7,13,21,22,26,28,29
VlnPlot(seurat_sct,
        features = Microglia,
        stack = TRUE,
        flip = TRUE) +
  NoLegend() +
  ggtitle("Microglia")

hirestiff(paste(prefixPCres,"VLN","Microglia","hires.tiff", sep = "_"))
lowrestiff(paste(prefixPCres,"VLN","Microglia","lowres.tiff", sep = "_"))

default_FP(seurat_sct,Microglia)

hirestiff(paste(prefixPCres,"UMAP","Microglia","genes","hires.tiff", sep = "_"))
lowrestiff(paste(prefixPCres,"UMAP","Microglia","genes","lowres.tiff", sep = "_"))

# Bipolar Cells 
# Clusters (0),1,2,3
VlnPlot(seurat_sct,
        features = c(Bipolar,"Snap25","Neto1"),
        stack = TRUE,
        flip = TRUE) +
  NoLegend() +
  ggtitle("Bipolar Cells")

hirestiff(paste(prefixPCres,"VLN","Bipolar","Cells","hires.tiff", sep = "_"))
lowrestiff(paste(prefixPCres,"VLN","Bipolar","Cells","lowres.tiff", sep = "_"))

default_FP(seurat_sct,Bipolar)

hirestiff(paste(prefixPCres,"UMAP","Bipolar","genes","hires.tiff", sep = "_"))
lowrestiff(paste(prefixPCres,"UMAP","Bipolar","genes","lowres.tiff", sep = "_"))

# Retinal Ganglion Cells 
# Clusters 4,6,12,18,19,21,24,25,26,28
VlnPlot(seurat_sct,
        features = c(Retinal_Ganglion,
                     "Slc17a6","Pou4f3","Cacna1c"),
        stack = TRUE,
        flip = TRUE,
        raster = F) +
  NoLegend() +
  ggtitle("Retinal Ganglion Cells")

hirestiff(paste(prefixPCres,"VLN","Retinal","Ganglion","Cells","hires.tiff", sep = "_"))
lowrestiff(paste(prefixPCres,"VLN","Retinal","Ganglion","Cells","lowres.tiff", sep = "_"))

default_FP(seurat_sct,Retinal_Ganglion)

hirestiff(paste(prefixPCres,"UMAP","Retinal_Ganglion","genes","hires.tiff", sep = "_"))
lowrestiff(paste(prefixPCres,"UMAP","Retinal_Ganglion","genes","lowres.tiff", sep = "_"))

FeaturePlot(seurat_sct, 
            reduction = "umap", 
            features = c(Retinal_Ganglion,
                         "Slc17a6","Pou4f3"), 
            order = TRUE,
            min.cutoff = 'q10',
            label = TRUE,
            label.size = 3,
            pt.size = 0.5,
            raster = F)

# Amacrine Cells 
# Cluster 17
VlnPlot(seurat_sct,
        features = Amacrine,
        stack = TRUE,
        flip = TRUE) +
  NoLegend() +
  ggtitle("Amacrine Cells")

hirestiff(paste(prefixPCres,"VLN","Amacrine","Cells","hires.tiff", sep = "_"))
lowrestiff(paste(prefixPCres,"VLN","Amacrine","Cells","lowres.tiff", sep = "_"))

default_FP(seurat_sct,Amacrine)

hirestiff(paste(prefixPCres,"UMAP","Amacrine","genes","hires.tiff", sep = "_"))
lowrestiff(paste(prefixPCres,"UMAP","Amacrine","genes","lowres.tiff", sep = "_"))

# Horizontal Cells 
# Clusters ?
VlnPlot(seurat_sct,
        features = c(Horizontal,"Onecut1","Gfap"),
        stack = TRUE,
        flip = TRUE,
        pt.size = 0) +
  NoLegend() +
  ggtitle("Horizontal Cells")

hirestiff(paste(prefixPCres,"VLN","Horizontal","Cells","hires.tiff", sep = "_"))
lowrestiff(paste(prefixPCres,"VLN","Horizontal","Cells","lowres.tiff", sep = "_"))

default_FP(seurat_sct,Horizontal)

hirestiff(paste(prefixPCres,"UMAP","Horizontal","genes","hires.tiff", sep = "_"))
lowrestiff(paste(prefixPCres,"UMAP","Horizontal","genes","lowres.tiff", sep = "_"))

# Pericytes
# Cluster 4,12,18,24,25,26
VlnPlot(seurat_sct,
        features = c(Pericytes,"Cacna1c"),
        stack = TRUE,
        flip = TRUE,
        pt.size = 0) +
  NoLegend() + 
  ggtitle("Pericytes")

hirestiff(paste(prefixPCres,"VLN","Pericytes","hires.tiff", sep = "_"))
lowrestiff(paste(prefixPCres,"VLN","Pericytes","lowres.tiff", sep = "_"))

default_FP(seurat_sct,Pericytes)

hirestiff(paste(prefixPCres,"UMAP","Pericytes","genes","hires.tiff", sep = "_"))
lowrestiff(paste(prefixPCres,"UMAP","Pericytes","genes","lowres.tiff", sep = "_"))

# Reticulocytes 23
VlnPlot(seurat_sct,
        features = c(Reticulocytes,"Gfap"),
        stack = TRUE,
        flip = TRUE,
        pt.size = 0) +
  NoLegend() +
  ggtitle("Reticulocytes")

hirestiff(paste(prefixPCres,"VLN","Reticulocytes","hires.tiff", sep = "_"))
lowrestiff(paste(prefixPCres,"VLN","Reticulocytes","lowres.tiff", sep = "_"))

default_FP(seurat_sct,Reticulocytes)

hirestiff(paste(prefixPCres,"UMAP","Reticulocytes","genes","hires.tiff", sep = "_"))
lowrestiff(paste(prefixPCres,"UMAP","Reticulocytes","genes","lowres.tiff", sep = "_"))

  #### retinal genes ####

Rods <- c("Rho","Gnat1","Nrl","Pde6a","Pde6b","Cngb1")
Rods1 <- c("Rho","Gnat1","Nrl")
Rods2 <- c("Pde6a","Pde6b","Cngb1")
Cones <- c("Arr3","Opn1sw","Gnat2","Pde6h")
Mueller_Glia <- c("Glul","Rlbp1")
Astrocytes <- c("Gfap","Aqp4","Gpr37")
Microglia <- c("Tyrobp","Cx3cr1","Tmem119")
Bipolar <- c("Scgn","Vsx1","Vsx2")
Retinal_Ganglion <- c("Rbpms","Pou4f1","Pou4f2")
Amacrine <- c("Calb1","Gad1","Gad2")
Horizontal <- c("Onecut2")
Pericytes <- c("Pecam1","Acta2")
Reticulocytes <- c("Hbb-bs")
RPEpi <- c("Mitf","Tjp1","Rpe65","Rlbp1","Best1")

retina_genes <- c(Rods,Cones,Mueller_Glia,Astrocytes,
                     Microglia,Bipolar,Retinal_Ganglion,
                     Amacrine,Horizontal,Reticulocytes)

FeaturePlot(seurat_sct, 
            reduction = "umap", 
            features = c(Horizontal,"Lhx1","Isl1","Prox1"), 
            order = TRUE,
            min.cutoff = 'q10',
            label = TRUE,
            label.size = 3,
            pt.size = 0.5)

FeaturePlot(seurat_sct, 
            reduction = "umap", 
            features = c("Ntrk2"), 
            order = TRUE,
            min.cutoff = 'q10',
            label = TRUE,
            label.size = 3,
            pt.size = 0.5)

  #### retinal genes from Chen et al 2022 ####

Rods <- c("Rho","Nr2e3","Nrl","Pdc","Rp1")
Cones <- c("Opn1sw","Opn1mw","Pde6h","Arr3","Gnat2")
Mueller_Glia <- c("Zfp36l1","Dbi","Apoe","Slc1a3","Sparc")
Bipolar <- c("Pcp2","Isl1","Grm6","Trnp1")
Amacrine <- c("C1ql1","Snhg11","Tfap2b","Pcsk1n")
Horizontal <- c("Slc4a3","Calb1","Septin4","Tpm3")
Microglia <- c("Ctsd","Ccl4","C1qb","C1qc")
Vascular <- c("Trpm1","Igfbp7")

# Astrocytes <- c("Gfap","Aqp4","Gpr37")
# Retinal_Ganglion <- c("Rbpms","Pou4f1","Pou4f2")
# Pericytes <- c("Pecam1","Acta2")
# Reticulocytes <- c("Hbb-bs")

retina_genes <- c(Rods,Cones,Mueller_Glia,Bipolar,
                  Amacrine,Horizontal,Microglia,Vascular)

Idents(seurat_sct) <- ident_res

Idents(seurat_sct) <- "Cell.Type"

Idents(seurat_sct) <- "Cluster.Type"

VlnPlot(seurat_sct,
        features = c("Olig2","Apc","Cnp",
                     "Nes","Fabp7","Pax6"),
        stack = T,
        flip = T) + NoLegend()

VlnPlot(seurat_sct,
        features = c("Casp3","Casp6","Bax",
                     "Bcl2","Fas","Faslg","Tnfrsf1a"),
        stack = T,
        flip = T) + NoLegend()

VlnPlot(seurat_sct,
        features = Rods,
        stack = T,
        flip = T) + NoLegend()

VlnPlot(seurat_sct,
        features = Cones,
        stack = T,
        flip = T) + NoLegend()

VlnPlot(seurat_sct,
        features = Mueller_Glia,
        stack = T,
        flip = T) + NoLegend()

VlnPlot(seurat_sct,
        features = Retinal_Ganglion,
        stack = T,
        flip = T) + NoLegend()

VlnPlot(seurat_sct,
        features = c("Glul","Rlbp1","Gfap","S100b"),
        stack = T,
        flip = T) + NoLegend()

VlnPlot(seurat_sct,
        features = Bipolar,
        stack = T,
        flip = T) + NoLegend()

VlnPlot(seurat_sct,
        features = c("Isl1","Grm6","Scgn","Vsx2","Prdm8",
                     "Prkca","Pcp2","Bhlhe23","Cabp5","Car8",
                     "Trpm1","Slc5a8","Cacna2d3","Adamts5"),
        stack = T,
        flip = T) + NoLegend()

VlnPlot(seurat_sct,
        features = c("Gsg1","Tmem215","Trnp1"),
        stack = T,
        flip = T) + NoLegend()

VlnPlot(seurat_sct,
        features = Amacrine,
        stack = T,
        flip = T) + NoLegend()

VlnPlot(seurat_sct,
        features = c("Snhg11","Tfap2b","Gad1","Gad2"),
        stack = T,
        flip = T) + NoLegend()

VlnPlot(seurat_sct,
        features = Horizontal,
        stack = T,
        flip = T) + NoLegend()

VlnPlot(seurat_sct,
        features = Microglia,
        stack = T,
        flip = T) + NoLegend()

VlnPlot(seurat_sct,
        features = c("Ccl4","C1qb","C1qc","Tyrobp","Cx3cr1","Tmem119"),
        stack = T,
        flip = T) + NoLegend()

VlnPlot(seurat_sct,
        features = Vascular,
        stack = T,
        flip = T) + NoLegend()

VlnPlot(seurat_sct,
        features = Retinal_Ganglion,
        stack = T,
        flip = T) + NoLegend()

VlnPlot(seurat_sct,
        features = Reticulocytes,
        # stack = T,
        flip = T) + NoLegend()

VlnPlot(seurat_sct,
        features = c("Ruvbl1","Tfrc","Ucp3","Folr1",
                     "Ube2o","Cd36","Itga4","Slc11a2"),
        stack = T,
        flip = T) + NoLegend() # reticulocyte markers

VlnPlot(seurat_sct,
        features = c("Tfrc","Cga","Wrn","Kit","Piezo1",
                     "Kmt5a","Tal1","Gata1","Lmo2","Eng",
                     "Cfp","Itgb1","Use1","Grsf1","Lrp8",
                     "Lmnb1","Gypa","Tfr2","Igbp1","Epor",
                     "Ets1","Bpgm","Cd44","Cdh1"),
        stack = T,
        flip = T) + NoLegend() # Erythroblast markers

DotPlot(seurat_sct,
        features = c("Tfrc","Cga","Wrn","Kit","Piezo1",
        "Kmt5a","Tal1","Gata1","Lmo2","Eng",
        "Cfp","Itgb1","Use1","Grsf1","Lrp8",
        "Lmnb1","Gypa","Tfr2","Igbp1","Epor",
        "Ets1","Bpgm","Cd44","Cdh1"),
        cluster.idents = T)



# Erythroid-like and erythroid precursor cells
DotPlot(seurat_sct,
        features = c("Ahsp","Hba-a1","Hba-a2","Hbb-bs"),
        cluster.idents = T)

# Muller
DotPlot(seurat_sct,
        features = c("S100a16","Dbi","Gnai2","Crb1",
                     "Rdh10","Abca8a","Crot","Dapl1",
                     "Itm2b"),
        cluster.idents = T)

default_FP(seurat_sct,Rods)
default_FP(seurat_sct,Cones)
default_FP(seurat_sct,Mueller_Glia)
default_FP(seurat_sct,Bipolar)
default_FP(seurat_sct,Amacrine)
default_FP(seurat_sct,Horizontal)
default_FP(seurat_sct,Microglia)
default_FP(seurat_sct,Vascular)


default_FP(seurat_sct,c("Trpm1","Prkca","Pcp2"))

default_FP(seurat_sct,c("Tp53"))

VlnPlot(seurat_sct,
        features = c("Tp53")) + NoLegend()

FeaturePlot(seurat_sct, 
            reduction = "umap", 
            features = c("Trpm1","Prkca","Pcp2","Bhlhe23"), 
            order = TRUE,
            min.cutoff = 'q10',
            label = TRUE,
            label.size = 3,
            pt.size = 0.5)

FeaturePlot(seurat_sct, 
            reduction = "umap", 
            features = c("Prdm8","Prkca","Pcp2","Bhlhe23"), 
            order = TRUE,
            min.cutoff = 'q10',
            label = TRUE,
            label.size = 3,
            pt.size = 0.5)

  #### dotplots ####

DotPlot(seurat_sct,
        features = c("Rho",    # 7 rods 
                     "Arr3",   # 8 cones
                     "Rlbp1",  # 4,10,12 muller glia
                     "Scgn",   # 6 BC
                     "Sebox",  # 0,1 RBCs
                     "Slc6a9", # 8,7,2,6 Gly-AC
                     "Gad1",   # 2,6 GABA-AC
                     "Thy1",   # 2,6 GABA-AC / GLY -AC
                     "Onecut2",# 6 HC
                     "Cx3cr1", # 16,3 MG
                     # "Pecam1", # NONE Pericytes
                     "Acta2",  # 2,16 Pericytes
                     "Hbb-bs"),# 15 Reticulocytes
        cluster.idents = TRUE)

DotPlot(seurat_sct,
        features = c("Rho","Gnat1","Nrl","Pde6a","Pde6b","Cngb1",    # rods 
                     "Arr3","Opn1sw","Gnat2","Pde6h",                # cones
                     "Glul","Rlbp1",                                 # mueller glia
                     "Gfap","Aqp4","Gpr37",                          # astrocytes
                     "Tyrobp","Cx3cr1","Tmem119",                    # Microglia
                     "Scgn","Vsx1","Vsx2",                           # Bipolar
                     "Rbpms","Pou4f1","Pou4f2",                      # Retinal Ganglion
                     "Calb1","Gad1","Gad2",                          # Amacrine
                     "Onecut2",                                      # Horizontal
                     "Pecam1","Acta2",                               # Pericytes
                     "Hbb-bs"),                                      # Reticulocytes
        cluster.idents = TRUE) + 
  theme(axis.text.x = element_text(angle = 45, vjust = 0.5))

Idents(p60seurat) <- ident_res

DotPlot(p60seurat,
        features = c("Rho","Gnat1","Nrl","Pde6a","Pde6b","Cngb1",    # rods 
                     "Arr3","Opn1sw","Gnat2","Pde6h",                # cones
                     "Glul","Rlbp1",                                 # mueller glia
                     "Gfap","Aqp4","Gpr37",                          # astrocytes
                     "Tyrobp","Cx3cr1","Tmem119",                    # Microglia
                     "Scgn","Vsx1","Vsx2",                           # Bipolar
                     "Rbpms","Pou4f1","Pou4f2",                      # Retinal Ganglion
                     "Calb1","Gad1","Gad2",                          # Amacrine
                     "Onecut2",                                      # Horizontal
                     "Pecam1","Acta2",                               # Pericytes
                     "Hbb-bs"),                                      # Reticulocytes
        cluster.idents = TRUE) + 
  theme(axis.text.x = element_text(angle = 45, vjust = 0.5))

Idents(p90seurat) <- ident_res

DotPlot(p90seurat,
        features = c("Rho","Gnat1","Nrl","Pde6a","Pde6b","Cngb1",    # rods 
                     "Arr3","Opn1sw","Gnat2","Pde6h",                # cones
                     "Glul","Rlbp1",                                 # mueller glia
                     "Gfap","Aqp4","Gpr37",                          # astrocytes
                     "Tyrobp","Cx3cr1","Tmem119",                    # Microglia
                     "Scgn","Vsx1","Vsx2",                           # Bipolar
                     "Rbpms","Pou4f1","Pou4f2",                      # Retinal Ganglion
                     "Calb1","Gad1","Gad2",                          # Amacrine
                     "Onecut2",                                      # Horizontal
                     "Pecam1","Acta2",                               # Pericytes
                     "Hbb-bs"),                                      # Reticulocytes
        cluster.idents = TRUE) + 
  theme(axis.text.x = element_text(angle = 45, vjust = 0.5))

DotPlot(seurat_sct,
        features = c(Rods,Cones,Mueller_Glia,Astrocytes,
                     Microglia,Bipolar,Retinal_Ganglion,
                     Amacrine,Horizontal,Pericytes,Reticulocytes),
        cluster.idents = TRUE) + 
  theme(axis.text.x = element_text(angle = 45, vjust = 0.5))

VlnPlot(seurat_sct,
        features = Retinal_Ganglion,
        stack = TRUE,
        flip = TRUE,
        sort = TRUE) + NoLegend()

# rods 7
VlnPlot(seurat_sct,
        features = Rods,
        stack = TRUE,
        flip = TRUE,
        sort = FALSE) + NoLegend()

FeaturePlot()

# cones 8
VlnPlot(seurat_sct,
        features = Cones,
        stack = TRUE,
        flip = TRUE,
        sort = TRUE) + NoLegend()

# Mueller_Glia 
VlnPlot(seurat_sct,
        features = Mueller_Glia,
        stack = TRUE,
        flip = TRUE,
        sort = TRUE) + NoLegend()




  #### Identify cells with non-zero expression for each gene ####
cells_with_expression <- WhichCells(seurat_sct, expression = )
all.genes <- rownames(seurat_sct)

# Create a list of all genes with expression greater than zero
genes_with_expression <- names(which(Sums(cells_with_expression) > 0))

which(rowSums(seurat_sct) > 0)

rowSums(seurat_sct) > 0
namesgenes <- names(rowSums(seurat_sct) > 0)
FeaturePlot(seurat_sct, 
            reduction = "umap", 
            features = c("Opn1sw"), 
            order = TRUE,
            min.cutoff = 'q10',
            label = TRUE,
            label.size = 3,
            pt.size = 0.5)

FeaturePlot(seurat_sct, 
            reduction = "umap", 
            features = c("Arr3"), 
            order = TRUE,
            min.cutoff = 'q10',
            label = TRUE,
            label.size = 3,
            pt.size = 0.5)

FeaturePlot(seurat_sct, 
            reduction = "umap", 
            features = c("Pax6","Tfap2b","Gad1","Slc6a9"), 
            order = TRUE,
            min.cutoff = 'q10',
            label = TRUE,
            label.size = 3,
            pt.size = 0.5)

FeaturePlot(seurat_sct, 
            reduction = "umap", 
            features = c("Tfap2b"), 
            order = TRUE,
            min.cutoff = 'q10',
            label = FALSE,
            label.size = 3,
            pt.size = 0.5,
            split.by = "treatment")

FeaturePlot(seurat_sct, 
            reduction = "umap", 
            features = c("Tmem119","Adgre1","Tyrobp"), 
            order = TRUE,
            min.cutoff = 'q10',
            label = TRUE,
            label.size = 3,
            pt.size = 0.5)

FeaturePlot(seurat_sct, 
            reduction = "umap", 
            features = c("Glul"), 
            order = TRUE,
            min.cutoff = 'q10',
            label = TRUE,
            label.size = 3,
            pt.size = 0.5)

FeaturePlot(seurat_sct, 
            reduction = "umap", 
            features = c("Gfap"), 
            order = TRUE,
            min.cutoff = 'q10',
            label = TRUE,
            label.size = 3,
            pt.size = 0.5)

FeaturePlot(seurat_sct, 
            reduction = "umap", 
            features = c("Gsg1","Vsx1"), 
            order = TRUE,
            min.cutoff = 'q10',
            label = TRUE,
            label.size = 3,
            pt.size = 0.5)

FeaturePlot(seurat_sct, 
            reduction = "umap", 
            features = c("Onecut1"), 
            order = TRUE,
            min.cutoff = 'q10',
            label = TRUE,
            label.size = 3,
            pt.size = 0.5)

FeaturePlot(seurat_sct, 
            reduction = "umap", 
            features = c("Grm6"), 
            order = TRUE,
            min.cutoff = 'q10',
            label = TRUE,
            label.size = 3,
            pt.size = 0.5)

FeaturePlot(seurat_sct, 
            reduction = "umap", 
            features = c("Rbpms"), 
            order = TRUE,
            min.cutoff = 'q10',
            label = TRUE,
            label.size = 3,
            pt.size = 0.5)

VlnPlot(seurat_sct,
        features = "Grm6",
        flip = TRUE,
        pt.size = 0)



  #### Wang et al 2022 ####

DimPlot(seurat_sct, 
        reduction = "umap", 
        label = TRUE,
        label.size = 5,
        pt.size = 1,
        raster = F,
        cols = sorted_palopal) + 
  NoLegend() 

DimPlot(seurat_sct, 
        reduction = "umap", 
        label = TRUE,
        label.size = 5,
        pt.size = 1,
        raster = F,
        cols = sorted_palopal,
        split.by = "time_point") + 
  NoLegend() 

# Rod -	Pde6a, Rho, Sag, Gnat1, Nrl
VlnPlot(seurat_sct,
        features = c("Pde6a", "Rho", "Sag", "Gnat1", "Nrl",
                     "Map2","Rom1"),
        stack = T,
        flip = T,
        raster = F) + NoLegend()

FeaturePlot(seurat_sct, 
            reduction = "umap", 
            features = c("Pde6a", "Rho", "Sag", "Gnat1", "Nrl",
                         "Map2"), 
            order = TRUE,
            min.cutoff = 'q10',
            label = TRUE,
            label.size = 3,
            pt.size = 0.5,
            raster = F)

# S Cone - Pde6h, Arr3, Gnat2, Opn1sw, Ccdc136, Ttr
VlnPlot(seurat_sct,
        features = c("Pde6h", "Arr3", "Gnat2", "Opn1sw", "Ccdc136" ,"Ttr"),
        stack = T,
        flip = T,
        raster = F) + NoLegend()

FeaturePlot(seurat_sct, 
            reduction = "umap", 
            features = c("Pde6h", "Arr3", "Gnat2", "Opn1mw"), 
            order = TRUE,
            min.cutoff = 'q10',
            label = TRUE,
            label.size = 3,
            pt.size = 0.5,
            raster = F)

FeaturePlot(seurat_sct, 
            reduction = "umap", 
            features = c("Ccdc136" ,"Ttr", "Vopp1", "Lmo4"), 
            order = TRUE,
            min.cutoff = 'q10',
            label = TRUE,
            label.size = 3,
            pt.size = 0.5,
            raster = F)

# M/L Cone - Pde6h, Arr3, Gnat2, Opn1mw/Opn1lw, Vopp1, Lmo4
VlnPlot(seurat_sct,
        features = c("Pde6h", "Arr3", "Gnat2", "Opn1mw", 
                     "Opn1lw", "Vopp1", "Lmo4"),
        stack = T,
        flip = T,
        raster = F) + NoLegend()

# HC - Pvalb, Lhx1, Pax6, Onecut1
VlnPlot(seurat_sct,
        features = c("Pvalb", "Lhx1", "Pax6", "Onecut1"),
        stack = T,
        flip = T,
        raster = F) + NoLegend()

FeaturePlot(seurat_sct, 
            reduction = "umap", 
            features = c("Pvalb", "Lhx1", "Pax6", "Onecut1"), 
            order = TRUE,
            min.cutoff = 'q10',
            label = TRUE,
            label.size = 3,
            pt.size = 0.5,
            raster = F)

# AC - Snhg11, Pax6, Slc32a1, Crabp1, Nrxn2, Gad1/Gad2, Maf, Tfap2a, Slc6a9, Ebf3
VlnPlot(seurat_sct,
        features = c("Snhg11", "Pax6", "Slc32a1", "Crabp1", 
                     "Nrxn2", "Gad1", "Gad2", "Maf", 
                     "Tfap2a", "Slc6a9", "Ebf3"),
        stack = T,
        flip = T) + NoLegend()

hirestiff(paste(prefixPCres,"UMAP","Amacrine","markers","hires.tiff", sep = "_"))
lowrestiff(paste(prefixPCres,"UMAP","Amacrine","markers","lowres.tiff", sep = "_"))

# RBC - Gsg1, Pax6, Slc6a9, Vsx2, Otx2, Prkca, Isl1, Grm6, Cabp5, Vstm2b, Casp7, Rpa1
VlnPlot(seurat_sct,
        features = c("Gsg1", "Pax6", "Slc6a9", "Vsx2", "Otx2", 
                     "Prkca", "Isl1", "Grm6", "Cabp5", "Vstm2b", 
                     "Casp7", "Rpa1"),
        stack = T,
        flip = T) + NoLegend()

FeaturePlot(seurat_sct, 
            reduction = "umap", 
            features = c("Gsg1", "Pax6", "Slc6a9", "Vsx2", "Otx2",
                         "Prkca", "Isl1", "Grm6", "Cabp5", "Vstm2b", 
                         "Casp7", "Rpa1"), 
            order = TRUE,
            min.cutoff = 'q10',
            label = TRUE,
            label.size = 3,
            pt.size = 0.5,
            raster = F)

# ON CBC - Gsg1, Pax6, Slc6a9, Vsx2, Otx2, App, Scgn, Vsx1, Isl1, Grm6
VlnPlot(seurat_sct,
        features = c("Gsg1", "Pax6", "Slc6a9", "Vsx2", "Otx2", 
                     "App", "Scgn", "Vsx1", "Isl1", "Grm6"),
        stack = T,
        flip = T) + NoLegend()

FeaturePlot(seurat_sct, 
            reduction = "umap", 
            features = c("App", "Scgn", "Vsx1", "Isl1", "Grm6"), 
            order = TRUE,
            min.cutoff = 'q10',
            label = TRUE,
            label.size = 3,
            pt.size = 0.5,
            raster = F)

# OFF CBC - Gsg1, Pax6, Slc6a9, Vsx2, Otx2, App, Scgn, Vsx1, Grin2b, Grik1
VlnPlot(seurat_sct,
        features = c("Gsg1", "Pax6", "Slc6a9", "Vsx2", "Otx2", 
                     "App", "Scgn", "Vsx1", "Grin2b", "Grik1"),
        stack = T,
        flip = T) + NoLegend()

FeaturePlot(seurat_sct, 
            reduction = "umap", 
            features = c("App", "Scgn", "Vsx1", "Grin2b", "Grik1"), 
            order = TRUE,
            min.cutoff = 'q10',
            label = TRUE,
            label.size = 3,
            pt.size = 0.5,
            raster = F)

# Müller - Gsta1, Pax6, Rlbp1, Aqp4, Slc1a3, Apoe, Dkk3, Gpr37, Rax, Hes1, Notch1, Glul
VlnPlot(seurat_sct,
        features = c("Gsta1", "Pax6", "Rlbp1", "Aqp4", "Slc1a3", 
                     "Apoe", "Dkk3", "Gpr37", "Rax", "Hes1", 
                     "Notch1", "Glul"),
        stack = T,
        flip = T,
        raster = F) + NoLegend()

FeaturePlot(seurat_sct, 
            reduction = "umap", 
            features = c("Gfap", "Aqp4", "S100b"), 
            order = TRUE,
            min.cutoff = 'q10',
            label = TRUE,
            label.size = 3,
            pt.size = 0.5,
            raster = F)

# Microglia - C1qa, Tmem119, Cx3cr1, Apoe, Pax6, P2ry12, Aif1
VlnPlot(seurat_sct,
        features = c("C1qa", "Tmem119", "Cx3cr1", "Apoe", "Pax6", 
                     "P2ry12", "Aif1"),
        stack = T,
        flip = T,
        raster = F) + NoLegend()

FeaturePlot(seurat_sct, 
            reduction = "umap", 
            features = c("C1qa", "Tmem119", "Cx3cr1", 
                         "P2ry12", "Aif1"), 
            order = TRUE,
            min.cutoff = 'q10',
            label = TRUE,
            label.size = 3,
            pt.size = 0.5,
            raster = F)

# EC - Vwf, Pecam1, Cldn5, Cdh5, Tek, Kdr, Flt1
VlnPlot(seurat_sct,
        features = c("Vwf", "Pecam1", "Cldn5", "Cdh5", "Tek", 
                     "Kdr", "Flt1", "Rbpms"),
        stack = T,
        flip = T,
        raster = F) + NoLegend()

# Pericyte - Rgs5, Kcnj8, Myl9, Cspg4, Pdgfrb, Myh11, Acta2
VlnPlot(seurat_sct,
        features = c("Rgs5","Kcnj8","Myl9","Cspg4","Pdgfrb", 
                     "Myh11","Acta2","Mgp", "Crip1", "Tpm2",
                     "Rbpms","Pgf","Igf1","Igf2r","Vegfa",
                     "Axl"),
        stack = T,
        flip = T,
        raster = F) + NoLegend()

VlnPlot(seurat_sct,
        features = c("Dnmbp","Tubb3","Nefl","Nefm","Nefh","Rbfox3","Map2","Syp","Gap43"),
        stack = T,
        flip = T,
        raster = F) + NoLegend()

FeaturePlot(seurat_sct, 
            reduction = "umap", 
            features = c("Dnmbp","Tubb3","Nefl","Nefm","Nefh",
                         "Rbfox3","Map2","Syp","Gap43"), 
            order = TRUE,
            min.cutoff = 'q10',
            label = TRUE,
            label.size = 3,
            pt.size = 0.5,
            raster = F)

VlnPlot(seurat_sct,
        features = c("Mgp", "Crip1", "Tpm2"),
        stack = T,
        flip = T,
        raster = F) + NoLegend()


# Macrophage - Cxcr4, Cd53, Ptprc
VlnPlot(seurat_sct,
        features = c("Cxcr4", "Cd53", "Ptprc"),
        stack = T,
        flip = T,
        raster = F) + NoLegend()

FeaturePlot(seurat_sct, 
            reduction = "umap", 
            features = c("Cxcr4", "Cd53", "Ptprc"), 
            order = TRUE,
            min.cutoff = 'q10',
            label = TRUE,
            label.size = 3,
            pt.size = 0.5,
            raster = F)

# RGC - Rbpms
VlnPlot(seurat_sct,
        features = c("Rbpms","Pou4f1","Pou4f2","Pou4f3","Thy1",
                     "Map2","Opn4","Slc17a6","Jam2","Tbr1"),
        stack = T,
        flip = T,
        raster = F) + NoLegend()

VlnPlot(seurat_sct,
        features = c("Rbpms","Pou4f1","Pou4f2",
                     "Jam2","Foxp1","Foxp2"),
        stack = T,
        flip = T,
        raster = F) + NoLegend()

# Ctxn3+Müller - Cd9, Penk, Ctxn3
VlnPlot(seurat_sct,
        features = c("Cd9", "Penk", "Ctxn3"),
        stack = T,
        flip = T,
        raster = F) + NoLegend()

# Ctxn3-Müller - Cd9, Penk
VlnPlot(seurat_sct,
        features = c("Cd9", "Penk"),
        stack = T,
        flip = T,
        raster = F) + NoLegend()

# EC - Pltp, Csrp2, Cxcl12, Id1, Ramp2, Slc2a1
VlnPlot(seurat_sct,
        features = c("Pltp", "Csrp2", "Cxcl12", "Id1", 
                     "Ramp2", "Slc2a1"),
        stack = T,
        flip = T,
        raster = F) + NoLegend()

# Pericyte - Mgp, Crip1, Tpm2
VlnPlot(seurat_sct,
        features = c("Mgp", "Crip1", "Tpm2"),
        stack = T,
        flip = T,
        raster = F) + NoLegend()

  #### Genes from 9/18 ####

# Egr11 in rods and cones - probably Egr1
# Amacrine: Snhg11, Pax6, Slc32a1, Crabp1, Nrxn2, 
#           Gad1/Gad2, Maf, Tfap2a, Slc6a9, Ebf3 - done above
# Muller glia and astrocytes: 
#           glutamine synthetase (GLUL), clusterin (CLU), 
#           and apolipoprotein E (APOE).
# NxnL1: metabolic dysfunction-glucosse metabolism.
# Hk2: a key enzyme for aerobic glycolysis and required for survival of photoreceptors during aging.
# Cebpd: upregulation is related to neuroprotection.

VlnPlot(seurat_sct,
        features = c("Egr1","Nxnl1","Hk2","Cebpd"),
        stack = T,
        flip = T) + NoLegend()

hirestiff(paste(prefixPCres,"VLN","Egr1","Nxnl1","Hk2","Cebpd","hires.tiff", sep = "_"))
lowrestiff(paste(prefixPCres,"VLN","Egr1","Nxnl1","Hk2","Cebpd","lowres.tiff", sep = "_"))

VlnPlot(seurat_sct,
        features = c("Egr1"),
        flip = T,
        idents = c(8,9),
        split.by = "time_point_treatment")

VlnPlot(seurat_sct,
        features = c("Egr1"),
        flip = T,
        idents = c("ROD","CONE"),
        split.by = "time_point_treatment")

hirestiff(paste(prefixPCres,"VLN","Egr1","rod-cone","hires.tiff", sep = "_"))
lowrestiff(paste(prefixPCres,"VLN","Egr1","rod-cone","lowres.tiff", sep = "_"))

VlnPlot(seurat_sct,
        features = c("Nxnl1"),
        flip = T,
        idents = c("ROD","CONE"),
        split.by = "time_point_treatment")

hirestiff(paste(prefixPCres,"VLN","Nxnl1","rod-cone","hires.tiff", sep = "_"))
lowrestiff(paste(prefixPCres,"VLN","Nxnl1","rod-cone","lowres.tiff", sep = "_"))

VlnPlot(seurat_sct,
        features = c("Hk2"),
        flip = T,
        idents = c("ROD","CONE"),
        split.by = "time_point_treatment")

hirestiff(paste(prefixPCres,"VLN","Hk2","rod-cone","hires.tiff", sep = "_"))
lowrestiff(paste(prefixPCres,"VLN","Hk2","rod-cone","lowres.tiff", sep = "_"))

VlnPlot(seurat_sct,
        features = c("Cebpd"),
        flip = T,
        idents = c("GLN1","GLN2","MULG"),
        split.by = "time_point_treatment")

hirestiff(paste(prefixPCres,"VLN","Cebpd","gln_mulg","hires.tiff", sep = "_"))
lowrestiff(paste(prefixPCres,"VLN","Cebpd","gln_mulg","lowres.tiff", sep = "_"))

VlnPlot(seurat_sct,
        features = c("Glul", "Clu", "Apoe"),
        stack = T,
        flip = T) + NoLegend()

hirestiff(paste(prefixPCres,"VLN","Muller","markers","09_18","hires.tiff", sep = "_"))
lowrestiff(paste(prefixPCres,"VLN","Muller","markers","09_18","lowres.tiff", sep = "_"))


#### Genes from 9/28 ####
# VlnPlot(seurat_sct,
#         features = c("Pou4f1", "Pou4f2",  "Pou4f3",  "Sncg", "Isl1", "Nefl", "Nefm",
#                      "Elavl4", "Eomes", "Calb1", "Calb2", "Opn4", "Fes", "Tpbg", "Irx3", 
#                      "Irx4", "Trb1", "Neurod2", "Penk", "Jam2", "Calca", "Spp1", "Plpp4", 
#                      "Coch", "Foxp1", "Foxp2", "Tusc5", "Satb2", "Cdk15", "Anxa3", "Slc17a6", 
#                      "Tbx20", "Serpine2", "Nmb", "Adcyap1", "Nmb", "Cartpt", "Mmp17", "Gpr88", 
#                      "Fam19a4", "Ebf3", "Col25a1", "A230065H16Rik", "Dcx"),
#         stack = T,
#         flip = T) + NoLegend()

VlnPlot(seurat_sct,
        features = c("Pou4f1","Pou4f2","Pou4f3","Sncg","Isl1", 
                     "Nefl","Nefm","Elavl4",#"Eomes",
                     "Calb1",
                     "Calb2","Opn4","Fes","Tpbg","Irx3"),
        stack = T,
        flip = T,
        idents = c(6,15,20)) + NoLegend() + ggtitle("RGC set 1")

hirestiff(paste(prefixPCres,"VLN","RGC","set","1","sept_28","hires.tiff", sep = "_"))
lowrestiff(paste(prefixPCres,"VLN","RGC","set","1","sept_28","lowres.tiff", sep = "_"))

VlnPlot(seurat_sct,
        features = c("Irx4","Trib1","Neurod2","Penk","Jam2",
                     "Calca","Spp1","Plpp4","Coch","Foxp1",
                     "Foxp2","Trarg1","Satb2","Cdk15","Anxa3"),
        stack = T,
        flip = T,
        idents = c(6,15,20)) + NoLegend() + ggtitle("RGC set 2")

hirestiff(paste(prefixPCres,"VLN","RGC","set","2","sept_28","hires.tiff", sep = "_"))
lowrestiff(paste(prefixPCres,"VLN","RGC","set","2","sept_28","lowres.tiff", sep = "_"))

VlnPlot(seurat_sct,
        features = c("Slc17a6","Tbx20","Serpine2","Nmb","Adcyap1",
                     "Nmb","Cartpt","Mmp17","Gpr88","Tafa4",
                     "Ebf3","Col25a1","Lbhd2","Dcx"),
        stack = T,
        flip = T,
        idents = c(6,15,20)) + NoLegend() + ggtitle("RGC set 3")

hirestiff(paste(prefixPCres,"VLN","RGC","set","3","sept_28","hires.tiff", sep = "_"))
lowrestiff(paste(prefixPCres,"VLN","RGC","set","3","sept_28","lowres.tiff", sep = "_"))

# VlnPlot(seurat_sct,
#         features = c("Gldn","Chl1","Qpct","Dgkg","Timp2",
#                      "Junb","Ifitm10","Fgf1","Sema5a","Ramp1",
#                      "Ostf1","Esrrg","Etl4","Kbtbd11","Ctxn3",
#                      "Igf11","Sdk1","Sdk2","Man1a","Ifi27"),
#         stack = T,
#         flip = T) + NoLegend()

VlnPlot(seurat_sct,
        features = c("Gldn","Chl1","Qpct","Dgkg","Timp2",
                     "Junb","Ifitm10","Fgf1","Sema5a","Ramp1"),
        stack = T,
        flip = T,
        idents = c(6,15,20)) + NoLegend() + ggtitle("RGC set 4")

hirestiff(paste(prefixPCres,"VLN","RGC","set","4","sept_28","hires.tiff", sep = "_"))
lowrestiff(paste(prefixPCres,"VLN","RGC","set","4","sept_28","lowres.tiff", sep = "_"))

VlnPlot(seurat_sct,
        features = c("Ostf1","Esrrg","Etl4","Kbtbd11","Ctxn3",
                     "Igf2","Sdk1","Sdk2","Man1a1","Ifi27"),
        stack = T,
        flip = T,
        idents = c(6,15,20)) + NoLegend() + ggtitle("RGC set 5")

hirestiff(paste(prefixPCres,"VLN","RGC","set","5","sept_28","hires.tiff", sep = "_"))
lowrestiff(paste(prefixPCres,"VLN","RGC","set","5","sept_28","lowres.tiff", sep = "_"))

VlnPlot(seurat_sct,
        features = c("Pdlim1","Col1a1","Pdgfrb","Rbpms"),
        stack = T,
        flip = T,
        idents = c(6,15,20)) + NoLegend() + ggtitle("Pericytes set 1")

hirestiff(paste(prefixPCres,"VLN","Pericytes","set","1","sept_28","hires.tiff", sep = "_"))
lowrestiff(paste(prefixPCres,"VLN","Pericytes","set","1","sept_28","lowres.tiff", sep = "_"))

VlnPlot(seurat_sct,
        features = c("Rgs5","Kcnj8","Myl9","Cspg4","Pdgfrb", 
                     "Myh11","Acta2","Mgp","Crip1","Tpm2",
                     "Pdlim1","Col1a1"),
        stack = T,
        flip = T,
        idents = c(6,15,20)) + NoLegend() + ggtitle("Pericytes set 1 with Wang et al markers")

hirestiff(paste(prefixPCres,"VLN","Pericytes","with","wang","markers","sept_28","hires.tiff", sep = "_"))
lowrestiff(paste(prefixPCres,"VLN","Pericytes","with","wang","markers","sept_28","lowres.tiff", sep = "_"))


VlnPlot(seurat_sct,
        features = c("Prph","Ebf3","Atp13a5",
                     "Cend1","Esrrb","Uchl1"),
        stack = T,
        flip = T,
        idents = c(6,15,20)) + 
  NoLegend() + 
  ggtitle("New markers 10-6-23")

hirestiff(paste(prefixPCres,"VLN","new","markers","Oct6","hires.tiff", sep = "_"))
lowrestiff(paste(prefixPCres,"VLN","new","markers","Oct6","lowres.tiff", sep = "_"))

Prph, Ebf3
Pericyte: Atp13a5

BM88(Cend1),  ERRbeta (Esrrb), PGP9.5 (Uchl1)


FeaturePlot(seurat_sct, 
            reduction = "umap", 
            features = c("Rbpms"), 
            order = TRUE,
            min.cutoff = 'q10',
            label = TRUE,
            label.size = 3,
            pt.size = 0.5,
            raster = F,
            split.by = "time_point")

FeaturePlot(seurat_sct, 
            reduction = "umap", 
            features = c("Sncg","Rbpms","Pou4f1","Kcnj8"), 
            order = TRUE,
            min.cutoff = 'q10',
            label = TRUE,
            label.size = 3,
            pt.size = 0.5,
            raster = F)

FeaturePlot(seurat_sct, 
            reduction = "umap", 
            features = c("Sncg","Rbpms","Isl1","Rbfox3"), 
            order = TRUE,
            min.cutoff = 'q10',
            label = TRUE,
            label.size = 3,
            pt.size = 0.5,
            raster = F)

FeaturePlot(seurat_sct, 
            reduction = "umap", 
            features = c("Sncg","Rbpms","Pou4f1","Pou4f2"), 
            order = TRUE,
            min.cutoff = 'q10',
            label = TRUE,
            label.size = 3,
            pt.size = 0.5,
            raster = F)

  #### Find Cluster Markers for all cells ####

# find all markers
cluster_markers <- FindAllMarkers(seurat_sct, only.pos = TRUE, min.pct = 0.50, logfc.threshold = 0.4)

write.csv(cluster_markers, 
          file = paste(prefixPCres,"clustermarkers_pos_only_minpct0.50_logfcthresh0.4.csv", sep = "_"))

  #### find DEGs between treatments ####

Idents(object = seurat_sct) <- "treatment"

cluster_markers <- FindAllMarkers(seurat_sct, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.4)

write.csv(cluster_markers, file = 
            paste(prefixPCres,"degs_treatment_all_cells_minpct0.25_logfcthresh0.4.csv", sep = "_"))

Idents(object = seurat_sct) <- ident_res

  #### find DEGs between time points ####

Idents(object = seurat_sct) <- "time_point"

cluster_markers <- FindAllMarkers(seurat_sct, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.4)

write.csv(cluster_markers, file = 
            paste(prefixPCres,"degs_timepoint_all_cells_minpct0.25_logfcthresh0.4.csv", sep = "_"))

Idents(object = seurat_sct) <- ident_res

  #### subsets for further DEG analysis ####

# subset seurat obj by time point for treated vs untreated DEGs
Idents(object = seurat_sct) <- "time_point"

p60seurat <- subset(seurat_sct, idents = "p60") 
p90seurat <- subset(seurat_sct, idents = "p90") 

Idents(object = p60seurat) <- "treatment"
cluster_markers <- FindAllMarkers(p60seurat, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.4)

write.csv(cluster_markers, file = 
            paste(prefixPCres,"degs_treatment_p60_minpct0.25_logfcthresh0.4.csv", sep = "_"))

Idents(object = p90seurat) <- "treatment"
cluster_markers <- FindAllMarkers(p90seurat, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.4)

write.csv(cluster_markers, file = 
            paste(prefixPCres,"degs_treatment_p90_minpct0.25_logfcthresh0.4.csv", sep = "_"))

# subset seurat obj by treatment for p60 vs p90 DEGs
Idents(object = seurat_sct) <- "treatment"

TRseurat <- subset(seurat_sct, idents = "Cell_Treated") 
UTseurat <- subset(seurat_sct, idents = "Untreated") 

Idents(object = TRseurat) <- "time_point"
cluster_markers <- FindAllMarkers(TRseurat, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.4)

write.csv(cluster_markers, file = 
            paste(prefixPCres,"degs_treated_p60_vs_p90_minpct0.25_logfcthresh0.4.csv", sep = "_"))

Idents(object = UTseurat) <- "time_point"
cluster_markers <- FindAllMarkers(UTseurat, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.4)

write.csv(cluster_markers, file = 
            paste(prefixPCres,"degs_untreated_p60_vs_p90_minpct0.25_logfcthresh0.4.csv", sep = "_"))

Idents(object = seurat_sct) <- ident_res
Idents(object = p60seurat) <- ident_res
Idents(object = p90seurat) <- ident_res
Idents(object = TRseurat) <- ident_res
Idents(object = UTseurat) <- ident_res

  #### save subset objects ####
saveRDS(p60seurat, paste0(prefixPCres,"_p60_seurat.rds"))

# p60seurat <- readRDS(file = paste0(prefixPCres,"_p60_seurat.rds"))

saveRDS(p90seurat, paste0(prefixPCres,"_p90_seurat.rds"))

p90seurat <- readRDS(file = paste0(prefixPCres,"_p90_seurat.rds"))

saveRDS(TRseurat, paste0(prefixPCres,"_cell_treated_seurat.rds"))

# TRseurat <- readRDS(file = paste0(prefixPCres,"_cell_treated_seurat.rds"))

saveRDS(UTseurat, paste0(prefixPCres,"_untreated_seurat.rds"))

# UTseurat <- readRDS(file = paste0(prefixPCres,"_untreated_seurat.rds"))

  #### for loop for DEGs within clusters  ####

clusterlist <- as.list(0:(max(as.integer(seurat_sct@meta.data$integrated_snn_res.0.2))-1))
clusterlist

for (k in 1:length(clusterlist)){
  
  nam <- paste("clusterDEG_treatment", clusterlist[k], sep = "_")
  assign(nam, FindMarkers(seurat_sct, ident.1 = "Cell_Treated", group.by = "treatment", subset.ident = clusterlist[k]))
  
  write.csv(get(nam), file = paste(prefixPCres, nam, "list.csv", sep = "_"), row.names = TRUE)
  
} 

# positive LFC means upregulation in Cell_Treated cells

for (k in 1:length(clusterlist)){
  
  nam <- paste("clusterDEG_timepoint", clusterlist[k], sep = "_")
  assign(nam, FindMarkers(seurat_sct, ident.1 = "p60", group.by = "time_point", subset.ident = clusterlist[k]))
  
  write.csv(get(nam), file = paste(prefixPCres, nam, "list.csv", sep = "_"), row.names = TRUE)
  
} 

# positive LFC means upregulation in p60 cells

#### LABEL CELL TYPES ####
  #### Label cell types (short name) ####
seurat_sct <- RenameIdents(object = seurat_sct, 
                           '0' = "RBC",
                           '1' = "RBC", 
                           '2' = "ROD", 
                           '3' = "MG",
                           '4' = "CBC",
                           '5' = "RBC",
                           '6' = "RGC",
                           '7' = "EC",
                           '8' = "ROD",
                           '9' = "MULG",
                           '10' = "MULG",
                           '11' = "CONE",
                           '12' = "MG",
                           '13' = "MG",
                           '14' = "CBC",
                           '15' = "RGC",
                           '16' = "AC",
                           '17' = "MULG",
                           '18' = "MULG",
                           '19' = "MULG",
                           '20' = "RGC",
                           '21' = "MG",
                           '22' = "EC",
                           '23' = "M0",
                           '24' = "MULG",
                           '25' = "ROD",
                           '26' = "EC",
                           '27' = "MULG",
                           '28' = "MG")

  #### save cell type (short name) as a metadata column ####
seurat_sct[["Cell.Type"]] <- seurat_sct@active.ident

  #### dimplots of cell type (short name) labeled clusters ####

DimPlot(seurat_sct, 
        reduction = "umap", 
        label = TRUE,
        label.size = 6, 
        pt.size = 1,
        raster = F) + 
  NoLegend()
hirestiff(paste(prefixPCres,"UMAP","celltypes","hires.tiff",sep = "_"))
lowrestiff(paste(prefixPCres,"UMAP","celltypes","lowres.tiff",sep = "_"))

DimPlot(seurat_sct, 
        reduction = "umap", 
        label = TRUE,
        label.size = 5, 
        pt.size = 1,
        raster = F,
        split.by = "time_point") + 
  NoLegend()
hirestiff(paste(prefixPCres,"UMAP","celltypes","by","time_point","hires.tiff",sep = "_"))
lowrestiff(paste(prefixPCres,"UMAP","celltypes","by","time_point","lowres.tiff",sep = "_"))

DimPlot(seurat_sct, 
        reduction = "umap", 
        label = TRUE,
        label.size = 3, 
        pt.size = 1,
        raster = F,
        split.by = "time_point_treatment") + 
  NoLegend()
hirestiff(paste(prefixPCres,"UMAP","celltypes","by","time_point_treatment","hires.tiff",sep = "_"))
lowrestiff(paste(prefixPCres,"UMAP","celltypes","by","time_point_treatment","lowres.tiff",sep = "_"))

for (i in 1:length(levels(seurat_sct$time_point_treatment))){
  p <- (DimPlot(seurat_sct,
                reduction = "umap",
                label = TRUE,
                label.size = 8,
                pt.size = 1,
                cells = c(WhichCells(seurat_sct, 
                                     expression = 
                                       time_point_treatment == 
                                       levels(seurat_sct$time_point_treatment)[i]))) +
          NoLegend() +
          ggtitle(paste0(levels(seurat_sct$time_point_treatment)[i])))
  
  print(p)
  
  hirestiff(paste(prefixPCres,"UMAP","split_by","time_point_treatment","with_cluster_labels",
                  levels(seurat_sct$time_point_treatment)[i],"hires.tiff",sep = "_"))
  lowrestiff(paste(prefixPCres,"UMAP","split_by","time_point_treatment","with_cluster_labels",
                   levels(seurat_sct$time_point_treatment)[i],"lowres.tiff",sep = "_"))
  hirestiffsquare(paste(prefixPCres,"UMAP","split_by","time_point_treatment","with_cluster_labels",
                        levels(seurat_sct$time_point_treatment)[i],"hires_square.tiff",sep = "_"))
  lowrestiffsquare(paste(prefixPCres,"UMAP","split_by","time_point_treatment","with_cluster_labels",
                         levels(seurat_sct$time_point_treatment)[i],"lowres_square.tiff",sep = "_"))
}

rm(p)
    #x#### Label cell types (long name) ####

seurat_sct <- RenameIdents(object = seurat_sct, 
                           'BIP1' = "Bipolar_1",
                           'BIP2' = "Bipolar_2", 
                           'GLN1' = "Ganglion_1",
                           'GLN2' = "Ganglion_2",
                           'MG' = "Microglia",
                           'ROD' = "Rods",
                           'CONE' = "Cones",
                           'MULG' = "Müller_Glia",
                           'AMCR' = "Amacrine",
                           'RTIC' = "Reticulocytes")

    #x#### save cell type (long name) as a metadata column ####
seurat_sct[["Cell.Type.Long"]] <- seurat_sct@active.ident

    #x#### Label cell types plus cluster number ####
seurat_sct <- RenameIdents(object = seurat_sct, 
                           '0' = "RGC0",
                           '1' = "BP1", 
                           '2' = "BP2", 
                           '3' = "BP_ROD3",
                           '4' = "FIB4",
                           '5' = "MG5",
                           '6' = "FIB6",
                           '7' = "MG7",
                           '8' = "ROD8",
                           '9' = "CONE9",
                           '10' = "MULG10",
                           '11' = "FIB_MULG_PGMT11",
                           '12' = "FIB12",
                           '13' = "MG13",
                           '14' = "MULG14",
                           '15' = "FIB_ASTR15",
                           '16' = "FIB_PGMT16",
                           '17' = "RGC_AMCR17",
                           '18' = "FIB18",
                           '19' = "FIB19",
                           '20' = "FIB_MEL20",
                           '21' = "MG21",
                           '22' = "MG22",
                           '23' = "BLOOD23",
                           '24' = "FIB24",
                           '25' = "FIB25",
                           '26' = "MG26",
                           '27' = "MULG27",
                           '28' = "MG28",
                           '29' = "MG29")

#### Label cell types
seurat_sct <- RenameIdents(object = seurat_sct, 
                           '0' = "BPC_0",
                           '1' = "BPC_1", 
                           '2' = "BPC_2", 
                           '3' = "BPC_3",
                           '4' = "RGC_4",
                           '5' = "MG_5",
                           '6' = "RGC_6",
                           '7' = "MG_7",
                           '8' = "ROD_8",
                           '9' = "CONE_9",
                           '10' = "MULG_10",
                           '11' = "MULG_11",
                           '12' = "RGC_12",
                           '13' = "MG_13",
                           '14' = "MULG_14",
                           '15' = "MULG_15",
                           '16' = "PGMT_16",
                           '17' = "AMCR_17",
                           '18' = "RGC_18",
                           '19' = "RGC_19",
                           '20' = "MEL_20",
                           '21' = "MG_21",
                           '22' = "MG_22",
                           '23' = "RTIC_23",
                           '24' = "RGC_24",
                           '25' = "RGC_25",
                           '26' = "MG_26",
                           '27' = "MULG_27",
                           '28' = "MG_28",
                           '29' = "MG_29")

    #x#### save cell type plus cluster numbers as a metadata column ####
seurat_sct[["Cluster.Type"]] <- seurat_sct@active.ident
seurat_sct[["Type.Cluster"]] <- seurat_sct@active.ident

saveRDS(seurat_sct, paste0(prefixPCres,"_seurat_after_RNAnorm_NEW_CELLTYPES.rds"))
seurat_sct <- readRDS(file = paste0(prefixPCres,"_seurat_after_RNAnorm_NEW_CELLTYPES.rds"))

    #x#### dimplots of cell type (long name) labeled clusters ####

DimPlot(seurat_sct, 
        reduction = "umap", 
        label = TRUE,
        label.size = 8, 
        pt.size = 1) + 
  NoLegend()
hirestiff(paste(prefixPCres,"UMAP","celltypesLN","hires.tiff",sep = "_"))
lowrestiff(paste(prefixPCres,"UMAP","celltypesLN","lowres.tiff",sep = "_"))

DimPlot(seurat_sct, 
        reduction = "umap", 
        label = TRUE,
        label.size = 6, 
        pt.size = 1,
        split.by = "time_point") + 
  NoLegend()
hirestiff(paste(prefixPCres,"UMAP","celltypesLN","by","time_point","hires.tiff",sep = "_"))
lowrestiff(paste(prefixPCres,"UMAP","celltypesLN","by","time_point","lowres.tiff",sep = "_"))

DimPlot(seurat_sct, 
        reduction = "umap", 
        label = TRUE,
        label.size = 6, 
        pt.size = 1,
        split.by = "treatment") + 
  NoLegend()
hirestiff(paste(prefixPCres,"UMAP","celltypesLN","by","treatment","hires.tiff",sep = "_"))
lowrestiff(paste(prefixPCres,"UMAP","celltypesLN","by","treatment","lowres.tiff",sep = "_"))

DimPlot(seurat_sct, 
        reduction = "umap", 
        label = TRUE,
        label.size = 4, 
        pt.size = 1,
        split.by = "time_point_treatment") + 
  NoLegend()
hirestiff(paste(prefixPCres,"UMAP","celltypesLN","by","time_point_treatment","hires.tiff",sep = "_"))
lowrestiff(paste(prefixPCres,"UMAP","celltypesLN","by","time_point_treatment","lowres.tiff",sep = "_"))

# set cluster colors
clust_cols <- scales::hue_pal()(length(levels(seurat_sct)))
names(clust_cols) <- levels(seurat_sct)

levels(seurat_sct)

# UMAP of cells in each cluster by time_point_treatment without cluster labels
DimPlot(seurat_sct, 
        reduction = "umap", 
        label = FALSE,
        cols = clust_cols, 
        split.by = "time_point_treatment",
        pt.size = 1) + 
  NoLegend() +
  ggtitle(paste0(ident_res))
hirestiff(paste(prefixPCres,"UMAP","split_by","time_point_treatment","without_cluster_labels","hires.tiff",sep = "_"))
lowrestiff(paste(prefixPCres,"UMAP","split_by","time_point_treatment","without_cluster_labels","lowres.tiff",sep = "_"))

for (i in 1:length(levels(seurat_sct$time_point_treatment))){
  p <- (DimPlot(seurat_sct,
                reduction = "umap",
                label = FALSE,
                pt.size = 1,
                cols = clust_cols,
                cells = c(WhichCells(seurat_sct, expression = time_point_treatment == levels(seurat_sct$time_point_treatment)[i]))) +
          NoLegend() +
          ggtitle(paste0(levels(seurat_sct$time_point_treatment)[i])))
  
  print(p)
  
  hirestiff(paste(prefixPCres,"UMAP","split_by","time_point_treatment","no_cluster_labels",
                  levels(seurat_sct$time_point_treatment)[i],"hires.tiff",sep = "_"))
  lowrestiff(paste(prefixPCres,"UMAP","split_by","time_point_treatment","no_cluster_labels",
                   levels(seurat_sct$time_point_treatment)[i],"lowres.tiff",sep = "_"))
  hirestiffsquare(paste(prefixPCres,"UMAP","split_by","time_point_treatment","no_cluster_labels",
                        levels(seurat_sct$time_point_treatment)[i],"hires_square.tiff",sep = "_"))
  lowrestiffsquare(paste(prefixPCres,"UMAP","split_by","time_point_treatment","no_cluster_labels",
                         levels(seurat_sct$time_point_treatment)[i],"lowres_square.tiff",sep = "_"))
}

# UMAP of cells in each cluster by time_point_treatment with cluster labels
DimPlot(seurat_sct, 
        reduction = "umap", 
        label = TRUE, 
        cols = clust_cols,
        split.by = "time_point_treatment",
        label.size = 5,
        pt.size = 1) +
  NoLegend() +
  ggtitle(paste0(ident_res))
hirestiff(paste(prefixPCres,"UMAP","split_by","time_point_treatment","with_cluster_labels","hires.tiff",sep = "_"))
lowrestiff(paste(prefixPCres,"UMAP","split_by","time_point_treatment","with_cluster_labels","lowres.tiff",sep = "_"))

for (i in 1:length(levels(seurat_sct$time_point_treatment))){
  p <- (DimPlot(seurat_sct,
                reduction = "umap",
                label = TRUE,
                label.size = 8,
                pt.size = 1,
                cols = clust_cols,
                cells = c(WhichCells(seurat_sct, expression = time_point_treatment == levels(seurat_sct$time_point_treatment)[i]))) +
          NoLegend() +
          ggtitle(paste0(levels(seurat_sct$time_point_treatment)[i])))
  
  print(p)
  
  hirestiff(paste(prefixPCres,"UMAP","split_by","time_point_treatment","with_cluster_labels",
                  levels(seurat_sct$time_point_treatment)[i],"hires.tiff",sep = "_"))
  lowrestiff(paste(prefixPCres,"UMAP","split_by","time_point_treatment","with_cluster_labels",
                   levels(seurat_sct$time_point_treatment)[i],"lowres.tiff",sep = "_"))
  hirestiffsquare(paste(prefixPCres,"UMAP","split_by","time_point_treatment","with_cluster_labels",
                        levels(seurat_sct$time_point_treatment)[i],"hires_square.tiff",sep = "_"))
  lowrestiffsquare(paste(prefixPCres,"UMAP","split_by","time_point_treatment","with_cluster_labels",
                         levels(seurat_sct$time_point_treatment)[i],"lowres_square.tiff",sep = "_"))
}

# UMAP of cells in each cluster by time_point_treatment with no cluster labels plus legend
DimPlot(seurat_sct, 
        reduction = "umap", 
        label = FALSE, 
        cols = clust_cols,
        label.size = 5,
        pt.size = 1) +
  ggtitle(paste0(ident_res))
hirestiff(paste(prefixPCres,"UMAP","with_no_labels_plus_legend","hires.tiff",sep = "_"))
lowrestiff(paste(prefixPCres,"UMAP","with_no_labels_plus_legend","lowres.tiff",sep = "_"))

DimPlot(seurat_sct, 
        reduction = "umap", 
        label = FALSE, 
        cols = clust_cols,
        split.by = "time_point_treatment",
        label.size = 5,
        pt.size = 1) +
  ggtitle(paste0(ident_res))
hirestiff(paste(prefixPCres,"UMAP","split_by","time_point_treatment","with_no_labels_plus_legend","hires.tiff",sep = "_"))
lowrestiff(paste(prefixPCres,"UMAP","split_by","time_point_treatment","with_no_labels_plus_legend","lowres.tiff",sep = "_"))

for (i in 1:length(levels(seurat_sct$time_point_treatment))){
  p <- (DimPlot(seurat_sct,
                reduction = "umap",
                label = FALSE,
                label.size = 8,
                pt.size = 1,
                cols = clust_cols,
                cells = c(WhichCells(seurat_sct, expression = time_point_treatment == levels(seurat_sct$time_point_treatment)[i]))) +
          ggtitle(paste0(levels(seurat_sct$time_point_treatment)[i])))
  
  print(p)
  
  hirestiff(paste(prefixPCres,"UMAP","split_by","time_point_treatment","with_no_labels_plus_legend",
                  levels(seurat_sct$time_point_treatment)[i],"hires.tiff",sep = "_"))
  lowrestiff(paste(prefixPCres,"UMAP","split_by","time_point_treatment","with_no_labels_plus_legend",
                   levels(seurat_sct$time_point_treatment)[i],"lowres.tiff",sep = "_"))
  hirestiffsquare(paste(prefixPCres,"UMAP","split_by","time_point_treatment","with_no_labels_plus_legend",
                        levels(seurat_sct$time_point_treatment)[i],"hires_square.tiff",sep = "_"))
  lowrestiffsquare(paste(prefixPCres,"UMAP","split_by","time_point_treatment","with_no_labels_plus_legend",
                         levels(seurat_sct$time_point_treatment)[i],"lowres_square.tiff",sep = "_"))
}

rm(p,i)



  #### Extract number of cells per cluster per orig.ident ####
n_cells <- FetchData(seurat_sct, vars = c("ident", "orig.ident")) %>%
  dplyr::count(ident, orig.ident) %>%
  tidyr::spread(ident, n)
n_cells
write.csv(n_cells, file = paste0(prefixPCres,"_celltypes_per_sample.csv"))

rm(n_cells)
  
  #### switch between cluster numbers and cell type idents ####

Idents(seurat_sct) <- ident_res

Idents(seurat_sct) <- "Cell.Type"

# Idents(seurat_sct) <- "Cell.Type.Long"
# 
# Idents(seurat_sct) <- "Cluster.Type"
# 
# Idents(seurat_sct) <- "Type.Cluster"

  #### check which ident is active by dimplot ####

unique(sort(seurat_sct@active.ident))

view(seurat_sct@meta.data)

DimPlot(seurat_sct, 
        reduction = "umap", 
        label = TRUE,
        label.size = 5,
        pt.size = 1,
        raster = F) + 
  NoLegend() 

DimPlot(seurat_sct, 
        reduction = "umap", 
        label = TRUE,
        label.size = 5,
        pt.size = 1,
        cells.highlight = WhichCells(seurat_sct, idents = "Rods")) + 
  NoLegend() 

DimPlot(seurat_sct, 
        reduction = "umap", 
        label = TRUE,
        label.size = 5,
        pt.size = 1,
        cells.highlight = WhichCells(seurat_sct, idents = "Cones")) + 
  NoLegend() 

DimPlot(seurat_sct, 
        reduction = "umap", 
        label = TRUE,
        label.size = 5,
        pt.size = 1,
        cols = palopal) + 
  NoLegend() 

  #### save RDS containing RNA normalized data AND CELLTYPE LABELS ####

# saveRDS(seurat_sct, paste0(prefixPCres,"_seurat_after_RNAnorm_CELLTYPES.rds"))

seurat_sct <- readRDS(file = paste0(prefixPCres,"_seurat_after_RNAnorm_CELLTYPES.rds"))


seurat_sct <- readRDS("~/Desktop/DEG_and_Pathway_Analysis/20230922_rat_retina_with_GEO/20230922_Rat_retina_with_GEO_int_12PCs_res.0.5_seurat_after_RNAnorm_CELLTYPES.rds")

View(seurat_sct@meta.data)

DimPlot(seurat_sct, reduction = 'umap', group.by = "Cell.Type", split.by = "treatment",label = T, pt.size = 0.5, raster=FALSE)
#### save seurat object as an hdf5 file for import into python ####

SaveH5Seurat(seurat_sct, 
             filename = paste(prefixPCres,"H5seurat_obj.H5seurat",sep = "_"), 
             overwrite = FALSE, 
             verbose = TRUE)

#### Find Cluster Markers for all cells as celltypes ####

# find all markers
cluster_markers <- FindAllMarkers(seurat_sct, 
                                  only.pos = TRUE, 
                                  min.pct = 0.50, 
                                  logfc.threshold = 0.4)

write.csv(cluster_markers, 
          file = paste(prefixPCres,"clustermarkers_CELLTYPES_pos_only_minpct0.50_logfcthresh0.4.csv", sep = "_"))

  #### find DEGs between treatments regardless of cell type ####

Idents(object = seurat_sct) <- "treatment"

cluster_markers <- FindAllMarkers(seurat_sct, 
                                  only.pos = TRUE, 
                                  min.pct = 0.25, 
                                  logfc.threshold = 0.4)

write.csv(cluster_markers, file = 
            paste(prefixPCres,"degs_treatment_all_cells_minpct0.25_logfcthresh0.4.csv", sep = "_"))

Idents(object = seurat_sct) <- ident_res

  #### find DEGs between time points regardless of cell type ####

Idents(object = seurat_sct) <- "time_point"

cluster_markers <- FindAllMarkers(seurat_sct, 
                                  only.pos = TRUE, 
                                  min.pct = 0.25, 
                                  logfc.threshold = 0.4)

write.csv(cluster_markers, file = 
            paste(prefixPCres,"degs_timepoint_all_cells_minpct0.25_logfcthresh0.4.csv", sep = "_"))

Idents(object = seurat_sct) <- ident_res

  #### find DEGs between time points and treatment regardless of cell type ####

Idents(object = seurat_sct) <- "time_point_treatment"

cluster_markers <- FindAllMarkers(seurat_sct, 
                                  only.pos = TRUE, 
                                  min.pct = 0.25, 
                                  logfc.threshold = 0.4)

write.csv(cluster_markers, file = 
            paste(prefixPCres,"degs_timepoint_treatment_all_cells_minpct0.25_logfcthresh0.4.csv", sep = "_"))

Idents(object = seurat_sct) <- ident_res

Idents(object = seurat_sct) <- "Cell.Type"
  #### for loop for DEGs within clusters TR vs UT ####

clusterlist <- c(levels(seurat_sct@meta.data$Cell.Type))

clusterlist
length(clusterlist)

for (k in 1:length(clusterlist)){
  message(crayon::red(paste("Calculating DEGs within",
                            clusterlist[k],
                            sep = " ")))
  
  nam <- paste("clusterDEG_TR_vs_UT", clusterlist[k], sep = "_")
  
  assign(nam, FindMarkers(seurat_sct, 
                          ident.1 = "Cell_Treated", 
                          group.by = "treatment", 
                          subset.ident = clusterlist[k]))
  
  write.csv(get(nam), file = paste(prefixPCres, 
                                   nam, 
                                   "list.csv", 
                                   sep = "_"), 
            row.names = TRUE)
  
} 

# positive LFC means upregulation in Cell_Treated cells

  #### for loop for DEGs within clusters by time_point ####

# get the name of clusters
clusterlist <- c(levels(seurat_sct@meta.data$Cell.Type))
clusterlist

# set the groups to compare
grouplist <- c(levels(seurat_sct@meta.data$time_point))
grouplist <- c(grouplist, grouplist[1])
grouplist

for (i in 1:(length(grouplist)-1)){
  for (k in 1:length(clusterlist)){
    message(crayon::red(paste("Calculating DEGs between",
                              grouplist[i],
                              "and",
                              grouplist[i+1],
                              "within the",
                              clusterlist[k],
                              "cluster",
                              sep = " ")))
    
    nam <- paste("clusterDEG",
                 grouplist[i],
                 "vs",
                 grouplist[i+1],
                 "in",
                 clusterlist[k],
                 sep = "_")
    
    assign(nam, FindMarkers(seurat_sct, 
                            ident.1 = grouplist[i], 
                            ident.2 = grouplist[i+1],
                            group.by = "time_point", 
                            subset.ident = clusterlist[k]))
    
    write.csv(get(nam), file = paste(prefixPCres, 
                                     nam, 
                                     "list.csv", 
                                     sep = "_"), 
              row.names = TRUE)
    
  } 
}



#### for loop for DEGs within clusters by time_point ####

# get the name of clusters
clusterlist <- c(levels(seurat_sct@meta.data$Cell.Type))
clusterlist

# set the groups to compare
grouplist <- c(levels(seurat_sct@meta.data$time_point_treatment))
grouplist

# Loop through each pair of groups (i vs j)
for (i in 1:(length(grouplist)-1)){
  for (j in (i+1):length(grouplist)){
    group1 <- grouplist[i]
    group2 <- grouplist[j]
    
    # Run find markers on group1 and group2 for each element of clusterlist
    # print the comparison to the terminal
    for (k in 1:length(clusterlist)){
      message(crayon::red(paste("Calculating DEGs between",
                                group1,
                                "and",
                                group2,
                                "within the",
                                clusterlist[k],
                                "cluster",
                                sep = " ")))
      
      # create the results variable unique to each pair of groups and 
      # element of clusterlist
      nam <- paste("clusterDEG",
                   group1,
                   "vs",
                   group2,
                   "in",
                   clusterlist[k],
                   sep = "_")
      
      # run FindMarkers and save results in the "nam" variable
      assign(nam, FindMarkers(seurat_sct,
                              ident.1 = group1,
                              ident.2 = group2,
                              group.by = "time_point_treatment",
                              subset.ident = clusterlist[k]))
      
      # write out the results
      write.csv(get(nam), file = paste(prefixPCres,
                                       nam,
                                       "list.csv",
                                       sep = "_"),
                row.names = TRUE)
    }
  } 
}



  #### subsets for further DEG analysis ####

# subset seurat obj by time point for treated vs untreated DEGs
Idents(seurat_sct) <- "time_point"

p60seurat <- subset(seurat_sct, idents = "p60")
rm(seurat_sct)
p90seurat <- subset(seurat_sct, idents = "p90") 
Idents(p90seurat) <- "Cell.Type"

Cone_p90seurat <- subset(p90seurat, idents = "CONE") 
Cone <- subset(Cone_p90seurat, subset = Arr3 > .5)

CONE <- WhichCells(object=Cone_p90seurat, expression = Arr3 > .5)
Idents(Cone) <- "treatment"
Cone_DEG <- FindMarkers(Cone,
                                 logfc.threshold = 0.2,
                                 min.pct = 0.2,
                        ident.1 = "Cell_Treated",
                        ident.2 = "Untreated")
setwd("~/Desktop/DEG_and_Pathway_Analysis/Rat_DEGs/RAT_Pathway_analysis_V5")
write.csv(Cone_DEG,"Cone_p90seurat1.csv")

Cone_DEG1 <- FindMarkers(Cone,
                        ident.1 = "Cell_Treated",
                        ident.2 = "Untreated")
write.csv(Cone_DEG1,"Cone_p90seurat.csv")




VlnPlot(Cone, features = "Tyrobp", pt.size=0)

FeaturePlot(Cone, features = "Apoe", min.cutoff = "q9", split.by = "treatment", cols = c("grey", "blue"), pt.size = .8)

DimPlot(p90seurat, reduction = 'umap', group.by = "Cell.Type", split.by = "treatment", label = T, pt.size = 0.5, raster=FALSE)
DimPlot(Cone_p90seurat, reduction = 'umap',  split.by = "treatment", label = T, pt.size = 0.5, raster=FALSE)
View(p90seurat@meta.data)
View(Cone_p90seurat@meta.data)
Idents(p60seurat) <- "treatment"
cluster_markers <- FindAllMarkers(p60seurat, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.4)

write.csv(cluster_markers, file = 
            paste(prefixPCres,"degs_treatment_p60_minpct0.25_logfcthresh0.4.csv", sep = "_"))

Idents(p90seurat) <- "treatment"
cluster_markers <- FindAllMarkers(p90seurat, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.4)

write.csv(cluster_markers, file = 
            paste(prefixPCres,"degs_treatment_p90_minpct0.25_logfcthresh0.4.csv", sep = "_"))

# subset seurat obj by treatment for p60 vs p90 DEGs
Idents(seurat_sct) <- "treatment"

TRseurat <- subset(seurat_sct, idents = "Cell_Treated") 
TRseurat@meta.data <- droplevels(TRseurat@meta.data)

UTseurat <- subset(seurat_sct, idents = "Untreated") 
UTseurat@meta.data <- droplevels(UTseurat@meta.data)

Idents(UTseurat) <- "time_point_treatment"

UTseuratp60p90 <- subset(UTseurat, 
                         idents = "WT_Untreated",
                         invert = T)
UTseuratp60p90@meta.data <- droplevels(UTseuratp60p90@meta.data)

Idents(object = seurat_sct) <- "Cell.Type"
# Idents(object = p60seurat) <- "Cell.Type"
# Idents(object = p90seurat) <- "Cell.Type"
Idents(object = TRseurat) <- "Cell.Type"
Idents(object = UTseurat) <- "Cell.Type"
Idents(object = UTseuratp60p90) <- "Cell.Type"

DimPlot(seurat_sct, 
        reduction = "umap", 
        label = TRUE,
        label.size = 5,
        pt.size = 1,
        raster = F) + 
  NoLegend() 

DimPlot(TRseurat, 
        reduction = "umap", 
        label = TRUE,
        label.size = 5,
        pt.size = 1,
        raster = F) + 
  NoLegend() 

DimPlot(UTseurat, 
        reduction = "umap", 
        label = TRUE,
        label.size = 5,
        pt.size = 1,
        raster = F) + 
  NoLegend() 

DimPlot(UTseuratp60p90, 
        reduction = "umap", 
        label = TRUE,
        label.size = 5,
        pt.size = 1,
        raster = F) + 
  NoLegend() 

  #### DEGs in RODS Treated p60 vs p90 ####

rod_tr_tp_markers <- FindMarkers(TRseurat,
                                 logfc.threshold = 0.4,
                                 min.pct = 0.25,
                                 only.pos = F, 
                                 ident.1 = "p60", 
                                 group.by = "time_point", 
                                 subset.ident = "ROD")

write.csv(rod_tr_tp_markers, file = 
            paste(prefixPCres,"degs_ROD_p60_TR_vs_p90_TR",
                  "minpct0.25_logfcthresh0.4.csv", sep = "_"))

  #### DEGs in RODS Untreated p60 vs p90 ####

rod_ut_tp_markers <- FindMarkers(UTseuratp60p90,
                                 logfc.threshold = 0.4,
                                 min.pct = 0.25,
                                 only.pos = F, 
                                 ident.1 = "p60", 
                                 group.by = "time_point", 
                                 subset.ident = "ROD")

write.csv(rod_tr_tp_markers, file = 
            paste(prefixPCres,"degs_ROD_p60_UT_vs_p90_UT",
                  "minpct0.25_logfcthresh0.4.csv", sep = "_"))

  #### DEGs in CONES Treated p60 vs p90 ####

cone_tr_tp_markers <- FindMarkers(TRseurat,
                                 logfc.threshold = 0.4,
                                 min.pct = 0.25,
                                 only.pos = F, 
                                 ident.1 = "p60", 
                                 group.by = "time_point", 
                                 subset.ident = "CONE")

write.csv(cone_tr_tp_markers, file = 
            paste(prefixPCres,"degs_CONE_p60_TR_vs_p90_TR",
                  "minpct0.25_logfcthresh0.4.csv", sep = "_"))

  #### DEGs in CONES Untreated p60 vs p90 ####

cone_ut_tp_markers <- FindMarkers(UTseuratp60p90,
                                 logfc.threshold = 0.4,
                                 min.pct = 0.25,
                                 only.pos = F, 
                                 ident.1 = "p60", 
                                 group.by = "time_point", 
                                 subset.ident = "CONE")

write.csv(cone_tr_tp_markers, file = 
            paste(prefixPCres,"degs_CONE_p60_UT_vs_p90_UT",
                  "minpct0.25_logfcthresh0.4.csv", sep = "_"))

  #### DEGs in MULG Treated p60 vs p90 ####

mulg_tr_tp_markers <- FindMarkers(TRseurat,
                                 logfc.threshold = 0.4,
                                 min.pct = 0.25,
                                 only.pos = F, 
                                 ident.1 = "p60", 
                                 group.by = "time_point", 
                                 subset.ident = "MULG")

write.csv(mulg_tr_tp_markers, file = 
            paste(prefixPCres,"degs_MULG_p60_TR_vs_p90_TR",
                  "minpct0.25_logfcthresh0.4.csv", sep = "_"))

  #### DEGs in MULG Untreated p60 vs p90 ####

mulg_ut_tp_markers <- FindMarkers(UTseuratp60p90,
                                 logfc.threshold = 0.4,
                                 min.pct = 0.25,
                                 only.pos = F, 
                                 ident.1 = "p60", 
                                 group.by = "time_point", 
                                 subset.ident = "MULG")

write.csv(mulg_tr_tp_markers, file = 
            paste(prefixPCres,"degs_MULG_p60_UT_vs_p90_UT",
                  "minpct0.25_logfcthresh0.4.csv", sep = "_"))

  #### save subset objects ####
saveRDS(p60seurat, paste0(prefixPCres,"_p60_seurat.rds"))

# p60seurat <- readRDS(file = paste0(prefixPCres,"_p60_seurat.rds"))

saveRDS(p90seurat, paste0(prefixPCres,"_p90_seurat.rds"))

p90seurat <- readRDS(file = paste0(prefixPCres,"_p90_seurat.rds"))

saveRDS(TRseurat, paste0(prefixPCres,"_cell_treated_seurat.rds"))

# TRseurat <- readRDS(file = paste0(prefixPCres,"_cell_treated_seurat.rds"))

saveRDS(UTseurat, paste0(prefixPCres,"_untreated_seurat.rds"))

# UTseurat <- readRDS(file = paste0(prefixPCres,"_untreated_seurat.rds"))

#### p60 TR vs UT DEGs within cell types #### 
Idents(p60seurat) <- "Cell.Type"

clusterlist <- c(levels(p60seurat@meta.data$Cell.Type))

clusterlist

for (k in 1:length(clusterlist)){
  
  nam <- paste("p60",
               "clusterDEG_treatment", 
               clusterlist[k], sep = "_")
  assign(nam, FindMarkers(p60seurat, 
                          ident.1 = "Cell_Treated", 
                          group.by = "treatment", 
                          subset.ident = clusterlist[k]))
  
  write.csv(get(nam), 
            file = paste(prefixPCres, 
                         nam, 
                         "list.csv", 
                         sep = "_"), 
            row.names = TRUE)
  
} 

# positive LFC means upregulation in Cell_Treated cells

#### p90 TR vs UT DEGs within cell types #### 
Idents(p90seurat) <- "Cell.Type"

clusterlist <- c(levels(p90seurat@meta.data$Cell.Type))

clusterlist

for (k in 1:length(clusterlist)){
  
  nam <- paste("p90",
               "clusterDEG_treatment", 
               clusterlist[k], sep = "_")
  assign(nam, FindMarkers(p90seurat, 
                          ident.1 = "Cell_Treated", 
                          group.by = "treatment", 
                          subset.ident = clusterlist[k]))
  
  write.csv(get(nam), 
            file = paste(prefixPCres, 
                         nam, 
                         "list.csv", 
                         sep = "_"), 
            row.names = TRUE)
  
} 

# positive LFC means upregulation in Cell_Treated cells

#### ACTIONet ####
library(ACTIONet)
library(SingleCellExperiment)
# convert seurat object to a single cell experiment
sce <- Seurat::as.SingleCellExperiment(seurat_sct)

# convert single cell experiment to an ACTIONet experiment
ace <- as.ACTIONetExperiment(sce)
ace = normalize.ace(ace)
ace = reduce.ace(ace)
ace = runACTIONet(ace)

# Annotate cell-types
data("curatedMarkers_human")
# markers = c(curatedMarkers_human$Brain$PFC$Mohammadi2020$marker.genes,
#             curatedMarkers_human$Brain$PFC$Velmeshev2019$marker.genes,
#             curatedMarkers_human$Brain$PFC$Schirmer2019$marker.genes,
#             curatedMarkers_human$Brain$PFC$MathysDavila2019$marker.genes,
#             curatedMarkers_human$Brain$PFC$Wang2018$marker.genes,
#             curatedMarkers_human$Brain$PFC$Layers$marker.genes)

markers <- c(curatedMarkers_human$Retina$marker.genes)

# List of human gene names
markers

h1 <- markers$Rods
h2 <- markers$Cones
h3 <- markers$RGCs
h4 <- markers$BPs
h5 <- markers$ACs
h6 <- markers$HCs
h7 <- markers$Macroglia
h8 <- markers$Microglia
h9 <- markers$Endo

  #### convert human genes to rat genes ####
# Load the biomaRt library
library(biomaRt)

# convert human to rat gene names

# listEnsemblArchives()
# https://oct2022.archive.ensembl.org DIDN"T WORK
# https://jul2022.archive.ensembl.org DIDN"T WORK
# https://apr2022.archive.ensembl.org DIDN"T WORK
# https://dec2021.archive.ensembl.org THIS ONE WORKED
# https://may2021.archive.ensembl.org 

human <- useMart("ensembl",
                dataset = "hsapiens_gene_ensembl",
                host = "https://dec2021.archive.ensembl.org")

rat <- useMart("ensembl",
              dataset = "rnorvegicus_gene_ensembl",
              host = "https://dec2021.archive.ensembl.org")

r1 = getLDS(attributes = c("hgnc_symbol"), 
            filters = "hgnc_symbol", 
            values = h1, 
            mart = human, 
            attributesL = c("rgd_symbol"), 
            martL = rat, 
            uniqueRows=T)

r2 <- getLDS(attributes = c("hgnc_symbol"), 
             filters = "hgnc_symbol", 
             values = h2, 
             mart = human, 
             attributesL = c("rgd_symbol"), 
             martL = rat, 
             uniqueRows=T)

r3 <- getLDS(attributes = c("hgnc_symbol"), 
             filters = "hgnc_symbol", 
             values = h3, 
             mart = human, 
             attributesL = c("rgd_symbol"), 
             martL = rat, 
             uniqueRows=T)

r4 <- getLDS(attributes = c("hgnc_symbol"), 
             filters = "hgnc_symbol", 
             values = h4, 
             mart = human, 
             attributesL = c("rgd_symbol"), 
             martL = rat, 
             uniqueRows=T)

r5 <- getLDS(attributes = c("hgnc_symbol"), 
             filters = "hgnc_symbol", 
             values = h5, 
             mart = human, 
             attributesL = c("rgd_symbol"), 
             martL = rat, 
             uniqueRows=T)

r6 <- getLDS(attributes = c("hgnc_symbol"), 
             filters = "hgnc_symbol", 
             values = h6, 
             mart = human, 
             attributesL = c("rgd_symbol"), 
             martL = rat, 
             uniqueRows=T)

r7 <- getLDS(attributes = c("hgnc_symbol"), 
             filters = "hgnc_symbol", 
             values = h7, 
             mart = human, 
             attributesL = c("rgd_symbol"), 
             martL = rat, 
             uniqueRows=T)

r8 <- getLDS(attributes = c("hgnc_symbol"), 
             filters = "hgnc_symbol", 
             values = h8, 
             mart = human, 
             attributesL = c("rgd_symbol"), 
             martL = rat, 
             uniqueRows=T)

r9 <- getLDS(attributes = c("hgnc_symbol"), 
             filters = "hgnc_symbol", 
             values = h9, 
             mart = human, 
             attributesL = c("rgd_symbol"), 
             martL = rat, 
             uniqueRows=T)

  #### make rat reference list ####
rat_markers <- list()

rat_markers$Rods <- r1$RGD.symbol
rat_markers$Cones <- r2$RGD.symbol
rat_markers$RGCs <- r3$RGD.symbol
rat_markers$BPs <- r4$RGD.symbol
rat_markers$ACs <- r5$RGD.symbol
rat_markers$HCs <- r6$RGD.symbol
rat_markers$Macroglia <- r7$RGD.symbol
rat_markers$Microglia <- r8$RGD.symbol
rat_markers$Endo <- r9$RGD.symbol

rat_markers



rat_markers$Rods <- c("Ahi1","Rho","Gnat1","Nrl","Nr2e3",
                      "Pde6a","Pde6b","Cngb1","Pdc","Rp1",
                      "Guca1a","Scaper","Ppef2")
rat_markers$Cones <- c("Arr3","Opn1mw","Opn1sw","Gnat2","Pde6h")
rat_markers$Mueller <- c("Glul","Rlbp1","Zfp36l1","Dbi","Apoe",
                         "Slc1a3","Sparc","Gfap","Aqp4","Gpr37",
                         "S100b")
rat_markers$BPs <- c("Scgn","Vsx1","Vsx2","Pcp2","Isl1",
                     "Grm6","Trnp1","Tmem215","Prkca","Camk2b")
rat_markers$ACs <- c("Calb1","Gad1","Gad2","C1ql1","C1ql2",
                     "Snhg11","Tfap2b","Tfap2c","Pcsk1n","Slc6a9")
rat_markers$HCs <- c("Onecut1","Onecut2","Lhx1","Slc4a3","Calb1",
                     "Septin4","Tpm3")
rat_markers$Macroglia <- c("Apoe","Slc1a3","Trpm3","Glul","Clu" )
rat_markers$Microglia <- c("Aif1","Tyrobp","Cx3cr1","Tmem119","Ctsd",
                           "Ccl4","C1qa","C1qb","C1qc","Cd163")
rat_markers$RGCs <- c("Rbpms","Pou4f1","Pou4f2","Slc17a6","Nefm",
                      "Nefl","Sncg" )
# rat_markers$Endo <- c("Tm4sf1","Apold1","Adamts9","Igfbp7","Cd34",
#                       "Rgs5","Cldn5","Cdh5")
# rat_markers$RPEpi <- c("Mitf","Tjp1","Rpe65","Rlbp1","Best1")

# rat_markers$Vascular <- c("Trpm1","Igfbp7")
# rat_markers$Pericytes <- c("Pecam1","Acta2")
rat_markers$Ret <- c("Hbb-bs")

rat_markers

#### annotate cells and export results to original seurat object ####
annot.out <- annotate.cells.using.markers(ace, rat_markers)

ace$ACEcelltypes <- annot.out$Label

# Export results to the seurat object as a new metadata column

seurat_sct[["ACEcelltypes"]] <- ace$ACEcelltypes

# rm(h1,h2,h3,h4,h5,h6,h7,h8,h9,
#    r1,r2,r3,r4,r5,r6,r7,r8,r9,
#    markers,rat_markers,human,rat,sce,
#    curatedMarkers_human,annot.out,ace)

  #### switch to ACEcelltypes ####

Idents(seurat_sct) <- "ACEcelltypes"

DimPlot(seurat_sct, 
        reduction = "umap", 
        label = TRUE,
        label.size = 7, 
        pt.size = 1) 

#### Module Scores for custom gene sets ####

irm_1125 <- read.csv("/Users/bells/Library/CloudStorage/Box-Box/20201124_MG_scRNAseq(old)/gene_module_lists/IRM_11_25.csv", header = FALSE)
irm_1125$V1 <- as.character(irm_1125$V1)
irm_1125 <- list(irm_1125$V1)
seurat_sct <- AddModuleScore(object = seurat_sct, features = irm_1125, search = TRUE, name = "IRM_11_25")

arm_pos_1125 <- read.csv("/Users/bells/Library/CloudStorage/Box-Box/20201124_MG_scRNAseq(old)/gene_module_lists/ARM_POS_11_25.csv", header = FALSE)
arm_pos_1125$V1 <- as.character(arm_pos_1125$V1)
arm_pos_1125 <- list(arm_pos_1125$V1)
seurat_sct <- AddModuleScore(object = seurat_sct, features = arm_pos_1125, search = TRUE, name = "ARM_POS_11_25")

arm_neg_1125 <- read.csv("/Users/bells/Library/CloudStorage/Box-Box/20201124_MG_scRNAseq(old)/gene_module_lists/ARM_NEG_11_25.csv", header = FALSE)
arm_neg_1125$V1 <- as.character(arm_neg_1125$V1)
arm_neg_1125 <- list(arm_neg_1125$V1)
seurat_sct <- AddModuleScore(object = seurat_sct, features = arm_neg_1125, search = TRUE, name = "ARM_NEG_11_25")

default_FP(seurat_sct,"IRM_11_251")
default_FP(seurat_sct,"ARM_POS_11_251")
default_FP(seurat_sct,"ARM_NEG_11_251")


   #### expression values for genes of interest ####

expression_data <- FetchData(object = seurat_sct, 
                             vars = c("genotype",ident_res,"Cell.Type",
                                      "ISL1","NEFL","NEFM",
                                      "NEFH","NGRN","PRPH","ELAVL3",
                                      "PFKP","RCAN1","SELPLG","STMN2","TARDBP"))

write.csv(expression_data, 
          file = paste(prefixPCres,"expression","in_all_cells","by_genotype.csv", sep = "_"))

#### VENN compare DEG gene lists ####
setwd(datadir)

# load data and convert the output "list" to a tibble 
# with the rownames converted to a column with the name "symbol"
p60_tr_vs_ut <- read.csv(file = "./20230327_Rat_retina_int_12PCs_res.0.2_degs_treatment_p60_minpct0.25_logfcthresh0.4.csv",
                             header = TRUE)
colnames(p60_tr_vs_ut) <- c("symbol","p_val","ave_log2FC","pct.1","pct.2","p_val_adj","cluster","gene")

p90_tr_vs_ut <- read.csv(file = "./20230327_Rat_retina_int_12PCs_res.0.2_degs_treatment_p90_minpct0.25_logfcthresh0.4.csv",
                            header = TRUE)
colnames(p90_tr_vs_ut) <- c("symbol","p_val","ave_log2FC","pct.1","pct.2","p_val_adj","cluster","gene")

p60_tr_vs_ut <- as_tibble(p60_tr_vs_ut)
p90_tr_vs_ut <- as_tibble(p90_tr_vs_ut)

# split tibbles by "cluster" and p_val_adj < 0.05
p60_up_in_tr <- filter(p60_tr_vs_ut, cluster == "Cell_Treated" & p_val_adj < 0.05)
p60_up_in_ut <- filter(p60_tr_vs_ut, cluster == "Untreated" & p_val_adj < 0.05)
p90_up_in_tr <- filter(p90_tr_vs_ut, cluster == "Cell_Treated" & p_val_adj < 0.05)
p90_up_in_ut <- filter(p90_tr_vs_ut, cluster == "Untreated" & p_val_adj < 0.05)

# venn diagrams of overlapping genes

### Up in Treated condition ###
treated_up <- list(p60_up_in_tr$gene,p90_up_in_tr$gene)
names(treated_up) <- c("p60_Up_in_Cell_Treated","p90_Up_in_Cell_Treated")

venn_treated_up <- Venn(treated_up)

plot(venn_treated_up, doWeights = TRUE, type = "circles")

v1 <- compute.Venn(venn_treated_up, doWeights = TRUE, type = "circles")

# change default labels
setlabels <- VennGetSetLabels(v1)
setlabels
facelabels <- VennGetFaceLabels(v1)
facelabels
universe <- VennGetUniverseRange(v1)
universe
# setlabels[setlabels$Label=="p60_Up_in_Cell_Treated","hjust"] <- "right"
# setlabels[setlabels$Label=="p90_Up_in_Cell_Treated","hjust"] <- "left"
# setlabels[setlabels$Label=="YoungV_vs_iMACS_13","x"] <- -1.5
# setlabels
# v1 <- VennSetSetLabels(v1,setlabels)

grid::grid.newpage()
plot(v1)

pdf(paste(prefixPCres,"venn","p60","vs","p90","up","in","Cell_Treated","overlapping","DEGs.pdf", sep = "_"))
print(plot(v1))
dev.off()

# get lists of non-overlapping and overlapping genes
p60_only_tr_up <- v1@IntersectionSets$'10'
p60_only_tr_up
overlap_tr_up <- v1@IntersectionSets$'11'
overlap_tr_up
p90_only_tr_up <- v1@IntersectionSets$'01'
p90_only_tr_up

write.csv(p60_only_tr_up, file = paste(prefixPCres,"venn","genes","up","in","p60","cell_treated.csv",sep = "_"))
write.csv(overlap_tr_up, file = paste(prefixPCres,"venn","genes","up","in","p60_and_p90","cell_treated.csv",sep = "_"))
write.csv(p90_only_tr_up, file = paste(prefixPCres,"venn","genes","up","in","p90","cell_treated.csv",sep = "_"))

### Up in Untreated condition ###
untreated_up <- list(p60_up_in_ut$gene,p90_up_in_ut$gene)
names(untreated_up) <- c("p60_Up_in_Untreated","p90_Up_in_Untreated")

venn_untreated_up <- Venn(untreated_up)

plot(venn_untreated_up, doWeights = TRUE, type = "circles")

v2 <- compute.Venn(venn_untreated_up, doWeights = TRUE, type = "circles")

# change default labels
# setlabels <- VennGetSetLabels(v2)
# setlabels
# facelabels <- VennGetFaceLabels(v2)
# facelabels
# universe <- VennGetUniverseRange(v2)
# universe
# setlabels[setlabels$Label=="p60_Up_in_Cell_Treated","hjust"] <- "right"
# setlabels[setlabels$Label=="p90_Up_in_Cell_Treated","hjust"] <- "left"
# setlabels[setlabels$Label=="YoungV_vs_iMACS_13","x"] <- -1.5
# setlabels
# v1 <- VennSetSetLabels(v1,setlabels)

grid::grid.newpage()
plot(v2)

pdf(paste(prefixPCres,"venn","p60","vs","p90","up","in","Untreated","overlapping","DEGs.pdf", sep = "_"))
print(plot(v2))
dev.off()

# get lists of non-overlapping and overlapping genes
p60_only_ut_up <- v2@IntersectionSets$'10'
p60_only_ut_up
overlap_ut_up <- v2@IntersectionSets$'11'
overlap_ut_up
p90_only_ut_up <- v2@IntersectionSets$'01'
p90_only_ut_up

write.csv(p60_only_ut_up, file = paste(prefixPCres,"venn","genes","up","in","p60","untreated.csv",sep = "_"))
write.csv(overlap_ut_up, file = paste(prefixPCres,"venn","genes","up","in","p60_and_p90","untreated.csv",sep = "_"))
write.csv(p90_only_ut_up, file = paste(prefixPCres,"venn","genes","up","in","p90","untreated.csv",sep = "_"))
#### Venn in ROD cluster ####

celltype <- "RODs"

# function to export DEG lists to csv
exportListToCSV <- function(data, filename) {
  # Determine the maximum length among all columns
  max_length <- max(sapply(data, length))
  
  # Pad shorter columns with empty strings
  filled_data <- lapply(data, function(x) {
    c(x, rep("", max_length - length(x)))
  })
  
  # Convert the list to a data frame
  df <- as.data.frame(filled_data)
  
  # Set the column names
  colnames(df) <- names(data)
  
  # Write the data frame to a CSV file
  write.csv(df, filename, row.names = FALSE)
}

# load the data
treatment_DEGs <- read.csv(file = "DEGs_p60_p90_by_treatment_RODs.csv",
                         header = TRUE)

p60_TR_UP <- treatment_DEGs$p60_TR
p60_UT_UP <- treatment_DEGs$p60_UT
p90_TR_UP <- treatment_DEGs$p90_TR
p90_UT_UP <- treatment_DEGs$p90_UT

rm(treatment_DEGs)

  #### p60-p90 overlap in Treated condition ####
treated_up <- list(p60_TR_UP,p90_TR_UP)
names(treated_up) <- c("p60_Up_in_Cell_Treated","p90_Up_in_Cell_Treated")

venn_treated_up <- Venn(treated_up)

plot(venn_treated_up, doWeights = TRUE, type = "circles")

v1 <- compute.Venn(venn_treated_up, doWeights = TRUE, type = "circles")

grid::grid.newpage()
plot(v1)

pdf(paste(prefixPCres,"venn","p60","p90","up","in","Cell_Treated","overlapping","DEGs",celltype,".pdf", sep = "_"))
print(plot(v1))
dev.off()

# get lists of non-overlapping and overlapping genes
p60_only_tr_up <- v1@IntersectionSets$'10'
p60_only_tr_up
overlap_tr_up <- v1@IntersectionSets$'11'
overlap_tr_up
p90_only_tr_up <- v1@IntersectionSets$'01'
p90_only_tr_up

p60p90_TR_up <- list(p60_only_tr_up = p60_only_tr_up,
                     overlap_tr_up = overlap_tr_up,
                     p90_only_tr_up = p90_only_tr_up)
p60p90_TR_up

exportListToCSV(p60p90_TR_up, paste(prefixPCres,"venn",
                                    "p60","p90","cell_treated",
                                    "overlapping","DEGs",celltype,
                                    ".csv",sep = "_"))

  #### p60 TR p90 UT overlap  ####
treated_up <- list(p60_TR_UP,p90_UT_UP)
names(treated_up) <- c("p60_Up_in_Cell_Treated","p90_Up_in_Unreated")

venn_treated_up <- Venn(treated_up)

plot(venn_treated_up, doWeights = TRUE, type = "circles")

v1 <- compute.Venn(venn_treated_up, doWeights = TRUE, type = "circles")

grid::grid.newpage()
plot(v1)

pdf(paste(prefixPCres,"venn","p60","TR","p90","UT",
          "overlapping","DEGs",celltype,".pdf", sep = "_"))
print(plot(v1))
dev.off()

# get lists of non-overlapping and overlapping genes
p60_only_tr_up <- v1@IntersectionSets$'10'
p60_only_tr_up
overlap_tr_up <- v1@IntersectionSets$'11'
overlap_tr_up
p90_only_ut_up <- v1@IntersectionSets$'01'
p90_only_ut_up

p60_TR_p90_UT <- list(p60_only_tr_up = p60_only_tr_up,
                     overlap_tr_up = overlap_tr_up,
                     p90_only_ut_up = p90_only_ut_up)
p60_TR_p90_UT

exportListToCSV(p60_TR_p90_UT, paste(prefixPCres,"venn",
                                    "p60","TR","p90","UT",
                                    "overlapping","DEGs",celltype,
                                    ".csv",sep = "_"))

  #### p60 UT p90 TR overlap  ####
treated_up <- list(p60_UT_UP,p90_TR_UP)
names(treated_up) <- c("p60_Up_in_Unreated","p90_Up_in_Cell_Treated")

venn_treated_up <- Venn(treated_up)

plot(venn_treated_up, doWeights = TRUE, type = "circles")

v1 <- compute.Venn(venn_treated_up, doWeights = TRUE, type = "circles")

grid::grid.newpage()
plot(v1)

pdf(paste(prefixPCres,"venn","p60","UT","p90","TR",
          "overlapping","DEGs",celltype,".pdf", sep = "_"))
print(plot(v1))
dev.off()

# get lists of non-overlapping and overlapping genes
p60_only_ut_up <- v1@IntersectionSets$'10'
p60_only_ut_up
overlap_tr_up <- v1@IntersectionSets$'11'
overlap_tr_up
p90_only_tr_up <- v1@IntersectionSets$'01'
p90_only_tr_up

p60_UT_p90_TR <- list(p60_only_ut_up = p60_only_ut_up,
                      overlap_tr_up = overlap_tr_up,
                      p90_only_tr_up = p90_only_tr_up)
p60_UT_p90_TR

exportListToCSV(p60_UT_p90_TR, paste(prefixPCres,"venn",
                                     "p60","UT","p90","TR",
                                     "overlapping","DEGs",celltype,
                                     ".csv",sep = "_"))

  #### p60-p90 overlap in Untreated condition ####
untreated_up <- list(p60_UT_UP,p90_UT_UP)
names(untreated_up) <- c("p60_Up_in_Untreated","p90_Up_in_Untreated")

venn_untreated_up <- Venn(untreated_up)

plot(venn_untreated_up, doWeights = TRUE, type = "circles")

v1 <- compute.Venn(venn_untreated_up, doWeights = TRUE, type = "circles")

grid::grid.newpage()
plot(v1)

pdf(paste(prefixPCres,"venn",
          "p60","p90",
          "up","in",
          "Untreated",
          "overlapping","DEGs",
          celltype,".pdf", sep = "_"))
print(plot(v1))
dev.off()

# get lists of non-overlapping and overlapping genes
p60_only_ut_up <- v1@IntersectionSets$'10'
p60_only_ut_up
overlap_ut_up <- v1@IntersectionSets$'11'
overlap_ut_up
p90_only_ut_up <- v1@IntersectionSets$'01'
p90_only_ut_up

p60p90_UT_up <- list(p60_only_ut_up = p60_only_ut_up,
                     overlap_ut_up = overlap_ut_up,
                     p90_only_ut_up = p90_only_ut_up)
p60p90_UT_up

exportListToCSV(p60p90_UT_up, paste(prefixPCres,"venn",
                                    "p60","p90","untreated",
                                    "overlapping","DEGs",celltype,
                                    ".csv",sep = "_"))

#### Venn in MG cluster ####

celltype <- "MG"

# function to export DEG lists to csv
exportListToCSV <- function(data, filename) {
  # Determine the maximum length among all columns
  max_length <- max(sapply(data, length))
  
  # Pad shorter columns with empty strings
  filled_data <- lapply(data, function(x) {
    c(x, rep("", max_length - length(x)))
  })
  
  # Convert the list to a data frame
  df <- as.data.frame(filled_data)
  
  # Set the column names
  colnames(df) <- names(data)
  
  # Write the data frame to a CSV file
  write.csv(df, filename, row.names = FALSE)
}

# load the data
treatment_DEGs <- read.csv(file = "DEGs_p60_p90_by_treatment_MG.csv",
                           header = TRUE)

p60_TR_UP <- treatment_DEGs$p60_TR
p60_UT_UP <- treatment_DEGs$p60_UT
p90_TR_UP <- treatment_DEGs$p90_TR
p90_UT_UP <- treatment_DEGs$p90_UT

rm(treatment_DEGs)

  #### p60-p90 overlap in Treated condition ####
treated_up <- list(p60_TR_UP,p90_TR_UP)
names(treated_up) <- c("p60_Up_in_Cell_Treated","p90_Up_in_Cell_Treated")

venn_treated_up <- Venn(treated_up)

plot(venn_treated_up, doWeights = TRUE, type = "circles")

v1 <- compute.Venn(venn_treated_up, doWeights = TRUE, type = "circles")

grid::grid.newpage()
plot(v1)

pdf(paste(prefixPCres,"venn","p60","p90","up","in","Cell_Treated","overlapping","DEGs",celltype,".pdf", sep = "_"))
print(plot(v1))
dev.off()

# get lists of non-overlapping and overlapping genes
p60_only_tr_up <- v1@IntersectionSets$'10'
p60_only_tr_up
overlap_tr_up <- v1@IntersectionSets$'11'
overlap_tr_up
p90_only_tr_up <- v1@IntersectionSets$'01'
p90_only_tr_up

p60p90_TR_up <- list(p60_only_tr_up = p60_only_tr_up,
                     overlap_tr_up = overlap_tr_up,
                     p90_only_tr_up = p90_only_tr_up)
p60p90_TR_up

exportListToCSV(p60p90_TR_up, paste(prefixPCres,"venn",
                                    "p60","p90","cell_treated",
                                    "overlapping","DEGs",celltype,
                                    ".csv",sep = "_"))

  #### p60 TR p90 UT overlap  ####
treated_up <- list(p60_TR_UP,p90_UT_UP)
names(treated_up) <- c("p60_Up_in_Cell_Treated","p90_Up_in_Untreated")

venn_treated_up <- Venn(treated_up)

plot(venn_treated_up, doWeights = TRUE, type = "circles")

v1 <- compute.Venn(venn_treated_up, doWeights = TRUE, type = "circles")

grid::grid.newpage()
plot(v1)

pdf(paste(prefixPCres,"venn","p60","TR","p90","UT",
          "overlapping","DEGs",celltype,".pdf", sep = "_"))
print(plot(v1))
dev.off()

# get lists of non-overlapping and overlapping genes
p60_only_tr_up <- v1@IntersectionSets$'10'
p60_only_tr_up
overlap_tr_up <- v1@IntersectionSets$'11'
overlap_tr_up
p90_only_ut_up <- v1@IntersectionSets$'01'
p90_only_ut_up

p60_TR_p90_UT <- list(p60_only_tr_up = p60_only_tr_up,
                      overlap_tr_up = overlap_tr_up,
                      p90_only_ut_up = p90_only_ut_up)
p60_TR_p90_UT

exportListToCSV(p60_TR_p90_UT, paste(prefixPCres,"venn",
                                     "p60","TR","p90","UT",
                                     "overlapping","DEGs",celltype,
                                     ".csv",sep = "_"))

  #### p60 UT p90 TR overlap  ####
treated_up <- list(p60_UT_UP,p90_TR_UP)
names(treated_up) <- c("p60_Up_in_Untreated","p90_Up_in_Cell_Treated")

venn_treated_up <- Venn(treated_up)

plot(venn_treated_up, doWeights = TRUE, type = "circles")

v1 <- compute.Venn(venn_treated_up, doWeights = TRUE, type = "circles")

grid::grid.newpage()
plot(v1)

pdf(paste(prefixPCres,"venn","p60","UT","p90","TR",
          "overlapping","DEGs",celltype,".pdf", sep = "_"))
print(plot(v1))
dev.off()

# get lists of non-overlapping and overlapping genes
p60_only_ut_up <- v1@IntersectionSets$'10'
p60_only_ut_up
overlap_tr_up <- v1@IntersectionSets$'11'
overlap_tr_up
p90_only_tr_up <- v1@IntersectionSets$'01'
p90_only_tr_up

p60_UT_p90_TR <- list(p60_only_ut_up = p60_only_ut_up,
                      overlap_tr_up = overlap_tr_up,
                      p90_only_tr_up = p90_only_tr_up)
p60_UT_p90_TR

exportListToCSV(p60_UT_p90_TR, paste(prefixPCres,"venn",
                                     "p60","UT","p90","TR",
                                     "overlapping","DEGs",celltype,
                                     ".csv",sep = "_"))

  #### p60-p90 overlap in Untreated condition ####
untreated_up <- list(p60_UT_UP,p90_UT_UP)
names(untreated_up) <- c("p60_Up_in_Untreated","p90_Up_in_Untreated")

venn_untreated_up <- Venn(untreated_up)

plot(venn_untreated_up, doWeights = TRUE, type = "circles")

v1 <- compute.Venn(venn_untreated_up, doWeights = TRUE, type = "circles")

grid::grid.newpage()
plot(v1)

pdf(paste(prefixPCres,"venn",
          "p60","p90",
          "up","in",
          "Untreated",
          "overlapping","DEGs",
          celltype,".pdf", sep = "_"))
print(plot(v1))
dev.off()

# get lists of non-overlapping and overlapping genes
p60_only_ut_up <- v1@IntersectionSets$'10'
p60_only_ut_up
overlap_ut_up <- v1@IntersectionSets$'11'
overlap_ut_up
p90_only_ut_up <- v1@IntersectionSets$'01'
p90_only_ut_up

p60p90_UT_up <- list(p60_only_ut_up = p60_only_ut_up,
                     overlap_ut_up = overlap_ut_up,
                     p90_only_ut_up = p90_only_ut_up)
p60p90_UT_up

exportListToCSV(p60p90_UT_up, paste(prefixPCres,"venn",
                                    "p60","p90","untreated",
                                    "overlapping","DEGs",celltype,
                                    ".csv",sep = "_"))

#### Trying ggVennDiagram and ggvenn instead of Vennerable ####

# ggVennDiagram
if (!require(devtools)) install.packages("devtools")
devtools::install_github("gaospecial/ggVennDiagram")

library("ggVennDiagram")

ggVennDiagram(treated_up,label_alpha = 0) +
  ggplot2::scale_fill_gradient(low="blue",high = "yellow")

# ggvenn
if (!require(devtools)) install.packages("devtools")
devtools::install_github("yanlinlin82/ggvenn")

library("ggvenn")

ggvenn(treated_up,
       # fill_color = c("green", "red"),
       stroke_size = 1,
       # set_name_size = 9,
       text_size = 8,
       show_percentage = FALSE,
       auto_scale = FALSE
)

?ggvenn
?geom_venn


#### Heatmap ####

Idents(seurat_sct) <- "treatment"

view(seurat_sct@meta.data)

p60_treatment_DEGs <- read_csv(file = "20230327_Rat_retina_int_12PCs_res.0.2_degs_treatment_p60_minpct0.25_logfcthresh0.4.csv")

p60_treatment_DEGs %>% sort()

p60_treatment_DEGs <- p60_treatment_DEGs[order(p60_treatment_DEGs$p_val_adj), ]

DoHeatmap(seurat_sct, features = c(p60_treatment_DEGs$gene[1:30]), 
          cells = WhichCells(seurat_sct, 
                             expression = time_point == "p60"))

p90_treatment_DEGs <- read_csv(file = "20230327_Rat_retina_int_12PCs_res.0.2_degs_treatment_p90_minpct0.25_logfcthresh0.4.csv")

DoHeatmap(seurat_sct, features = p90_treatment_DEGs$gene)

#### volcano plots ####

p60_treatment_DEGs <- read_csv(file = "20230327_Rat_retina_int_12PCs_res.0.2_degs_treatment_p60_minpct0.25_logfcthresh0.4.csv")

TR_DEGS <- subset(p60_treatment_DEGs, subset = cluster == "Cell_Treated")
UT_DEGS <- subset(p60_treatment_DEGs, subset = cluster == "Untreated")

UT_DEGS$avg_log2FC <- -1 * UT_DEGS$avg_log2FC

DEGS <- rbind(TR_DEGS,UT_DEGS)

volc <- DEGS

for (i in 1:length(volc$p_val_adj)) {
  if (volc$p_val_adj[i] == 0) {
    volc$p_val_adj[i] <- 1e-293
  }
}

EnhancedVolcano(volc,
                lab = volc$gene,
                x = 'avg_log2FC',
                y = 'p_val_adj',
                xlim = c(min(volc$avg_log2FC, na.rm = TRUE) - 0.1, 
                         max(volc$avg_log2FC, na.rm = TRUE) + 0.1),
                ylim = c(0, max(-log10(volc$p_val_adj), na.rm = TRUE) + 100),
                # title = "Aging iMPs vs Aging Veh",
                titleLabSize = 30,
                caption = NULL,
                subtitle = NULL,
                xlab = bquote(bold(~Log["2"] ~ "fold change")),
                ylab = bquote(bold(~"-" ~Log["10"] ~ p_val_adj)),
                axisLabSize = 25,
                # col = volc_colors,
                colAlpha = 1/1,
                legendLabels = c("NS", 
                                 expression(Log[2] ~ FC), 
                                 "p_val_adj", 
                                 expression(p_val_adj ~ and ~ log[2] ~ FC)),
                legendPosition = 'bottom',
                legendLabSize = 25,
                legendIconSize = 7,
                drawConnectors = TRUE,
                min.segment.length = 1,
                labSize = 7,
                max.overlaps = Inf,
                pCutoff = 0.05,
                FCcutoff = 0.58,
                gridlines.major = FALSE,
                gridlines.minor = FALSE,
                arrowheads = FALSE,
                pointSize = 7
) 

p90_treatment_DEGs <- read_csv(file = "20230327_Rat_retina_int_12PCs_res.0.2_degs_treatment_p90_minpct0.25_logfcthresh0.4.csv")

TR_DEGS <- subset(p90_treatment_DEGs, subset = cluster == "Cell_Treated")
UT_DEGS <- subset(p90_treatment_DEGs, subset = cluster == "Untreated")

UT_DEGS$avg_log2FC <- -1 * UT_DEGS$avg_log2FC

DEGS <- rbind(TR_DEGS,UT_DEGS)

volc <- DEGS

for (i in 1:length(volc$p_val_adj)) {
  if (volc$p_val_adj[i] == 0) {
    volc$p_val_adj[i] <- 1e-293
  }
}

EnhancedVolcano(volc,
                lab = volc$gene,
                x = 'avg_log2FC',
                y = 'p_val_adj',
                xlim = c(min(volc$avg_log2FC, na.rm = TRUE) - 0.1, 
                         max(volc$avg_log2FC, na.rm = TRUE) + 0.1),
                ylim = c(0, max(-log10(volc$p_val_adj), na.rm = TRUE) + 100),
                # title = "Aging iMPs vs Aging Veh",
                titleLabSize = 30,
                caption = NULL,
                subtitle = NULL,
                xlab = bquote(bold(~Log["2"] ~ "fold change")),
                ylab = bquote(bold(~"-" ~Log["10"] ~ p_val_adj)),
                axisLabSize = 25,
                # col = volc_colors,
                colAlpha = 1/1,
                legendLabels = c("NS", 
                                 expression(Log[2] ~ FC), 
                                 "p_val_adj", 
                                 expression(p_val_adj ~ and ~ log[2] ~ FC)),
                legendPosition = 'bottom',
                legendLabSize = 25,
                legendIconSize = 7,
                drawConnectors = TRUE,
                min.segment.length = 1,
                labSize = 3,
                max.overlaps = Inf,
                pCutoff = 0.05,
                FCcutoff = 1,
                gridlines.major = FALSE,
                gridlines.minor = FALSE,
                arrowheads = FALSE,
                pointSize = 3
) 


#### ####
#### Rod subset recluster ####
s_rods <- subset(seurat_sct,
                 idents = c("ROD"))

  #### SCT subset ####
DefaultAssay(s_rods) <- "RNA"

s_rods <- SCTransform(s_rods,
                      vars.to.regress = c("percent.mt","percent.ribo"), 
                      verbose = TRUE)  


  #### update project variable to include subset ####
project_ROD <- "Rat_retina_RODS"

  #### sub PCA ####
# this performs PCA on the seurat object
s_rods <- RunPCA(s_rods, npcs = 50, verbose = TRUE)

# make PC coordinate object a data frame
xx.coord <- as.data.frame(s_rods@reductions$pca@cell.embeddings)

# make PC feature loadings object a data frame
xx.gload <- as.data.frame(s_rods@reductions$pca@feature.loadings)

# calculate eigenvalues for arrays

# generate squares of all sample coordinates
sq.xx.coord <- as.data.frame(xx.coord^2)
# create empty list for eigenvalues first
eig <- c()
# calculate the eigenvalue for each PC in sq.xx.coord by taking the sqrt of the sum of squares
for(i in 1:ncol(sq.xx.coord))
  eig[i] = sqrt(sum(sq.xx.coord[,i]))
# calculate the total variance by adding up all the eigenvalues
sum.eig <- sum(eig)
# calculate the expected contribution of all PCs if they all contribute equally to the total variance
expected.contribution <- sum.eig/(length(xx.coord)-1)
# return the number of principal components with an eigenvalue greater than expected by equal variance
meaningful.PCs <- sum(eig > expected.contribution)

# create empty list for eigenvalue percentage
eig.percent <- c()
# calculate the percentage of the total variance by each PC eigenvalue
for(i in 1:length(eig))
  eig.percent[i] = 100*eig[i]/sum.eig
# sum of all eig.percent should total to 100
sum(eig.percent)
# create empty list for scree values
scree <- c()
# calculate a running total of variance contribution
for(i in 1:length(eig))
  if(i == 1) scree[i] = eig.percent[i] else scree[i] = scree[i-1] + eig.percent[i]

# create data frame for eigenvalue summaries
eigenvalues <- data.frame("PC" = colnames(xx.coord), "eig" = eig, "percent" = eig.percent, "scree" = scree)

# write csv for eigenvalues
# write.csv(eigenvalues, file = paste0("./",date,"_",project,"_PCA_eigenvalues.csv"), row.names = F)

# plot scree values
plot(eigenvalues$percent, ylim = c(0,100), type = "S", xlab = "PC", ylab = "Percent of variance",
     main = paste0(date,"_",project," scree plot all samples PCA"))
points(eigenvalues$scree, ylim = c(0,100), type = "p", pch = 16)
lines(eigenvalues$scree)
# add red line to indicate cut-off
cut.off <- 100/(length(eig)-1)
abline(h = cut.off, col = "red")
# add blue line to indicate which PCs are meaningful and kept
abline(v = meaningful.PCs, col = "blue")
text(meaningful.PCs, cut.off, label = paste("cutoff PC",meaningful.PCs),
     adj = c(-0.1, -0.5))

dev.copy(pdf, paste0("./",date,"_",project_ROD,"_scree_plot.pdf"))
dev.off()

rm(eigenvalues,sq.xx.coord,xx.coord,xx.gload,cut.off,
   eig,eig.percent,expected.contribution,i,scree,sum.eig)

# meaningful.PCs <- 9

  #### sub Run UMAP and look at UMAP plots ####

s_rods <- RunUMAP(s_rods, 
                  reduction = "pca", 
                  dims = 1:meaningful.PCs, 
                  verbose = TRUE)

# update prefixed variable
ROD_prefixPC <- paste0("./",date,"_",project_ROD,"_",meaningful.PCs,"PCs")

## UMAP plot by sample name ("orig.ident")
DimPlot(s_rods, reduction = "umap", label = FALSE, 
        pt.size = .25, group.by = "orig.ident")
hirestiff(paste0(ROD_prefixPC,"_UMAP_by_sample_hires.tiff"))
lowrestiff(paste0(ROD_prefixPC,"_UMAP_by_sample_lowres.tiff"))

DimPlot(s_rods, reduction = "umap", label = FALSE, 
        pt.size = .25, group.by = "treatment")
hirestiff(paste0(ROD_prefixPC,"_UMAP_by_treatment_hires.tiff"))
lowrestiff(paste0(ROD_prefixPC,"_UMAP_by_treatment_lowres.tiff"))

DimPlot(s_rods, reduction = "umap", label = FALSE, 
        pt.size = .25, split.by = "treatment", group.by = "treatment")
hirestiff(paste0(ROD_prefixPC,"_UMAP_split_by_treatment_hires.tiff"))
lowrestiff(paste0(ROD_prefixPC,"_UMAP_split_by_treatment_lowres.tiff"))

DimPlot(s_rods, reduction = "umap", label = FALSE, 
        pt.size = .25, group.by = "batch")
hirestiff(paste0(ROD_prefixPC,"_UMAP_by_batch_hires.tiff"))
lowrestiff(paste0(ROD_prefixPC,"_UMAP_by_batch_lowres.tiff"))

DimPlot(s_rods, reduction = "umap", label = FALSE, 
        pt.size = .25, group.by = "batch", split.by = "batch")
hirestiff(paste0(ROD_prefixPC,"_UMAP_split_by_batch_hires.tiff"))
lowrestiff(paste0(ROD_prefixPC,"_UMAP_split_by_batch_lowres.tiff"))

DimPlot(s_rods, reduction = "umap", label = FALSE, 
        pt.size = .25, group.by = "time_point")
hirestiff(paste0(ROD_prefixPC,"_UMAP_by_time_point_hires.tiff"))
lowrestiff(paste0(ROD_prefixPC,"_UMAP_by_time_point_lowres.tiff"))

DimPlot(s_rods, reduction = "umap", label = FALSE, 
        pt.size = .25, split.by = "time_point", group.by = "time_point")
hirestiff(paste0(ROD_prefixPC,"_UMAP_split_by_time_point_hires.tiff"))
lowrestiff(paste0(ROD_prefixPC,"_UMAP_split_by_time_point_lowres.tiff"))

DimPlot(s_rods, reduction = "umap", label = FALSE, 
        pt.size = .25, group.by = "time_point_treatment")
hirestiff(paste0(ROD_prefixPC,"_UMAP_by_time_point_treatment_hires.tiff"))
lowrestiff(paste0(ROD_prefixPC,"_UMAP_by_time_point_treatment_lowres.tiff"))

DimPlot(s_rods, reduction = "umap", label = FALSE, 
        pt.size = 1, split.by = "time_point_treatment", 
        group.by = "time_point_treatment")
hirestiff(paste0(ROD_prefixPC,"_UMAP_split_by_time_point_treatment_hires.tiff"))
lowrestiff(paste0(ROD_prefixPC,"_UMAP_split_by_time_point_treatment_lowres.tiff"))

DimPlot(s_rods, reduction = "umap", label = FALSE, pt.size = .25, group.by = "time_point_batch")
hirestiff(paste0(ROD_prefixPC,"_UMAP_by_time_point_batch_hires.tiff"))
lowrestiff(paste0(ROD_prefixPC,"_UMAP_by_time_point_batch_lowres.tiff"))

DimPlot(s_rods, reduction = "umap", label = FALSE, pt.size = .25, split.by = "time_point_batch", group.by = "time_point_batch")
hirestiff(paste0(ROD_prefixPC,"_UMAP_split_by_time_point_batch_hires.tiff"))
lowrestiff(paste0(ROD_prefixPC,"_UMAP_split_by_time_point_batch_lowres.tiff"))

DimPlot(s_rods, reduction = "umap", label = FALSE, pt.size = .25, group.by = "treatment_batch")
hirestiff(paste0(ROD_prefixPC,"_UMAP_by_treatment_batch_hires.tiff"))
lowrestiff(paste0(ROD_prefixPC,"_UMAP_by_treatment_batch_lowres.tiff"))

DimPlot(s_rods, reduction = "umap", label = FALSE, pt.size = .25, split.by = "treatment_batch", group.by = "treatment_batch")
hirestiff(paste0(ROD_prefixPC,"_UMAP_split_by_treatment_batch_hires.tiff"))
lowrestiff(paste0(ROD_prefixPC,"_UMAP_split_by_treatment_batch_lowres.tiff"))

DimPlot(s_rods, reduction = "umap", label = FALSE, pt.size = .25, group.by = "time_point_treatment_batch")
hirestiff(paste0(ROD_prefixPC,"_UMAP_by_time_point_treatment_batch_hires.tiff"))
lowrestiff(paste0(ROD_prefixPC,"_UMAP_by_time_point_treatment_batch_lowres.tiff"))

DimPlot(s_rods, reduction = "umap", label = FALSE, pt.size = .25, split.by = "time_point_treatment_batch", group.by = "time_point_treatment_batch")
hirestiff(paste0(ROD_prefixPC,"_UMAP_split_by_time_point_treatment_batch_hires.tiff"))
lowrestiff(paste0(ROD_prefixPC,"_UMAP_split_by_time_point_treatment_batch_lowres.tiff"))

# saveRDS(s_rods, file = paste0(ROD_prefixPC,"_seurat_integrated_preclustering.rds"))

# s_rods <- read_rds(file = paste0(ROD_prefixPC,"_seurat_integrated_preclustering.rds"))

  #### Remove outlier batch 08182 ####

Idents(s_rods) <- "batch"

s_rods <- subset(s_rods,
                 idents = c("08182"),
                 invert = TRUE)
  #### (rerun) SCT subset ####
DefaultAssay(s_rods) <- "RNA"

s_rods <- SCTransform(s_rods,
                      vars.to.regress = c("percent.mt","percent.ribo"), 
                      verbose = TRUE)  

  #### (rerun) sub PCA ####
# this performs PCA on the seurat object
s_rods <- RunPCA(s_rods, npcs = 50, verbose = TRUE)
# make PC coordinate object a data frame
xx.coord <- as.data.frame(s_rods@reductions$pca@cell.embeddings)
# make PC feature loadings object a data frame
xx.gload <- as.data.frame(s_rods@reductions$pca@feature.loadings)
# calculate eigenvalues for arrays
# generate squares of all sample coordinates
sq.xx.coord <- as.data.frame(xx.coord^2)
# create empty list for eigenvalues first
eig <- c()
# calculate the eigenvalue for each PC in sq.xx.coord by taking the sqrt of the sum of squares
for(i in 1:ncol(sq.xx.coord))
  eig[i] = sqrt(sum(sq.xx.coord[,i]))
# calculate the total variance by adding up all the eigenvalues
sum.eig <- sum(eig)
# calculate the expected contribution of all PCs if they all contribute equally to the total variance
expected.contribution <- sum.eig/(length(xx.coord)-1)
# return the number of principal components with an eigenvalue greater than expected by equal variance
ROD_meaningful.PCs <- sum(eig > expected.contribution)
# create empty list for eigenvalue percentage
eig.percent <- c()
# calculate the percentage of the total variance by each PC eigenvalue
for(i in 1:length(eig))
  eig.percent[i] = 100*eig[i]/sum.eig
# sum of all eig.percent should total to 100
sum(eig.percent)
# create empty list for scree values
scree <- c()
# calculate a running total of variance contribution
for(i in 1:length(eig))
  if(i == 1) scree[i] = eig.percent[i] else scree[i] = scree[i-1] + eig.percent[i]
# create data frame for eigenvalue summaries
eigenvalues <- data.frame("PC" = colnames(xx.coord), "eig" = eig, "percent" = eig.percent, "scree" = scree)

# plot scree values
plot(eigenvalues$percent, ylim = c(0,100), type = "S", xlab = "PC", ylab = "Percent of variance",
     main = paste0(date,"_",project," scree plot all samples PCA"))
points(eigenvalues$scree, ylim = c(0,100), type = "p", pch = 16)
lines(eigenvalues$scree)
# add red line to indicate cut-off
cut.off <- 100/(length(eig)-1)
abline(h = cut.off, col = "red")
# add blue line to indicate which PCs are meaningful and kept
abline(v = ROD_meaningful.PCs, col = "blue")
text(meaningful.PCs, cut.off, label = paste("cutoff PC",ROD_meaningful.PCs),
     adj = c(-0.1, -0.5))

dev.copy(pdf, paste0("./",date,"_",project_ROD,"_scree_plot.pdf"))
dev.off()

rm(eigenvalues,sq.xx.coord,xx.coord,xx.gload,cut.off,
   eig,eig.percent,expected.contribution,i,scree,sum.eig)

# ROD_meaningful.PCs <- 11

  #### (rerun) sub Run UMAP and look at UMAP plots ####

s_rods <- RunUMAP(s_rods, 
                  reduction = "pca", 
                  dims = 1:ROD_meaningful.PCs, 
                  verbose = TRUE)

# update prefixed variable
ROD_prefixPC <- paste0("./",date,"_",project_ROD,"_",ROD_meaningful.PCs,"PCs")

## UMAP plot by sample name ("orig.ident")
DimPlot(s_rods, reduction = "umap", label = FALSE, 
        pt.size = .25, group.by = "orig.ident")
hirestiff(paste0(ROD_prefixPC,"_UMAP_by_sample_hires.tiff"))
lowrestiff(paste0(ROD_prefixPC,"_UMAP_by_sample_lowres.tiff"))

DimPlot(s_rods, reduction = "umap", label = FALSE, 
        pt.size = .25, group.by = "treatment")
hirestiff(paste0(ROD_prefixPC,"_UMAP_by_treatment_hires.tiff"))
lowrestiff(paste0(ROD_prefixPC,"_UMAP_by_treatment_lowres.tiff"))

DimPlot(s_rods, reduction = "umap", label = FALSE, 
        pt.size = .25, split.by = "treatment", group.by = "treatment")
hirestiff(paste0(ROD_prefixPC,"_UMAP_split_by_treatment_hires.tiff"))
lowrestiff(paste0(ROD_prefixPC,"_UMAP_split_by_treatment_lowres.tiff"))

DimPlot(s_rods, reduction = "umap", label = FALSE, 
        pt.size = .25, group.by = "batch")
hirestiff(paste0(ROD_prefixPC,"_UMAP_by_batch_hires.tiff"))
lowrestiff(paste0(ROD_prefixPC,"_UMAP_by_batch_lowres.tiff"))

DimPlot(s_rods, reduction = "umap", label = FALSE, 
        pt.size = .25, group.by = "batch", split.by = "batch")
hirestiff(paste0(ROD_prefixPC,"_UMAP_split_by_batch_hires.tiff"))
lowrestiff(paste0(ROD_prefixPC,"_UMAP_split_by_batch_lowres.tiff"))

DimPlot(s_rods, reduction = "umap", label = FALSE, 
        pt.size = .25, group.by = "time_point")
hirestiff(paste0(ROD_prefixPC,"_UMAP_by_time_point_hires.tiff"))
lowrestiff(paste0(ROD_prefixPC,"_UMAP_by_time_point_lowres.tiff"))

DimPlot(s_rods, reduction = "umap", label = FALSE, 
        pt.size = .25, split.by = "time_point", group.by = "time_point")
hirestiff(paste0(ROD_prefixPC,"_UMAP_split_by_time_point_hires.tiff"))
lowrestiff(paste0(ROD_prefixPC,"_UMAP_split_by_time_point_lowres.tiff"))

DimPlot(s_rods, reduction = "umap", label = FALSE, 
        pt.size = .25, group.by = "time_point_treatment")
hirestiff(paste0(ROD_prefixPC,"_UMAP_by_time_point_treatment_hires.tiff"))
lowrestiff(paste0(ROD_prefixPC,"_UMAP_by_time_point_treatment_lowres.tiff"))

DimPlot(s_rods, reduction = "umap", label = FALSE, 
        pt.size = 1, split.by = "time_point_treatment", 
        group.by = "time_point_treatment")
hirestiff(paste0(ROD_prefixPC,"_UMAP_split_by_time_point_treatment_hires.tiff"))
lowrestiff(paste0(ROD_prefixPC,"_UMAP_split_by_time_point_treatment_lowres.tiff"))

DimPlot(s_rods, reduction = "umap", label = FALSE, pt.size = .25, group.by = "time_point_batch")
hirestiff(paste0(ROD_prefixPC,"_UMAP_by_time_point_batch_hires.tiff"))
lowrestiff(paste0(ROD_prefixPC,"_UMAP_by_time_point_batch_lowres.tiff"))

DimPlot(s_rods, reduction = "umap", label = FALSE, pt.size = .25, split.by = "time_point_batch", group.by = "time_point_batch")
hirestiff(paste0(ROD_prefixPC,"_UMAP_split_by_time_point_batch_hires.tiff"))
lowrestiff(paste0(ROD_prefixPC,"_UMAP_split_by_time_point_batch_lowres.tiff"))

DimPlot(s_rods, reduction = "umap", label = FALSE, pt.size = .25, group.by = "treatment_batch")
hirestiff(paste0(ROD_prefixPC,"_UMAP_by_treatment_batch_hires.tiff"))
lowrestiff(paste0(ROD_prefixPC,"_UMAP_by_treatment_batch_lowres.tiff"))

DimPlot(s_rods, reduction = "umap", label = FALSE, pt.size = .25, split.by = "treatment_batch", group.by = "treatment_batch")
hirestiff(paste0(ROD_prefixPC,"_UMAP_split_by_treatment_batch_hires.tiff"))
lowrestiff(paste0(ROD_prefixPC,"_UMAP_split_by_treatment_batch_lowres.tiff"))

DimPlot(s_rods, reduction = "umap", label = FALSE, pt.size = .25, group.by = "time_point_treatment_batch")
hirestiff(paste0(ROD_prefixPC,"_UMAP_by_time_point_treatment_batch_hires.tiff"))
lowrestiff(paste0(ROD_prefixPC,"_UMAP_by_time_point_treatment_batch_lowres.tiff"))

DimPlot(s_rods, reduction = "umap", label = FALSE, pt.size = .25, split.by = "time_point_treatment_batch", group.by = "time_point_treatment_batch")
hirestiff(paste0(ROD_prefixPC,"_UMAP_split_by_time_point_treatment_batch_hires.tiff"))
lowrestiff(paste0(ROD_prefixPC,"_UMAP_split_by_time_point_treatment_batch_lowres.tiff"))

# saveRDS(s_rods, file = paste0(ROD_prefixPC,"_seurat_integrated_preclustering.rds"))

# s_rods <- read_rds(file = paste0(ROD_prefixPC,"_seurat_integrated_preclustering.rds"))

  #### sub Clustering and Resolution ####
# DefaultAssay(s_rods) <- "integrated"

# Determine the K-nearest neighbor graph
s_rods <- FindNeighbors(object = s_rods, reduction = "pca", dims = 1:ROD_meaningful.PCs)

# Determine the clusters                              
s_rods <- FindClusters(object = s_rods,
                           resolution = c(0.1,0.2,0.3,0.4,0.5))

ROD_res <- "_res.0.1"
ROD_ident_res <- paste0("SCT_snn",ROD_res)

Idents(s_rods) <- ROD_ident_res
DimPlot(s_rods, reduction = "umap", label = TRUE, label.size = 5, pt.size = 0.8) + 
  NoLegend() + 
  ggtitle(paste0(ROD_ident_res))

#update prefix
ROD_prefixPCres <- paste0(ROD_prefixPC,ROD_res)

# Plot the UMAP
DimPlot(s_rods, reduction = "umap", label = TRUE, label.size = 6, pt.size = 1)  + 
  NoLegend() + 
  ggtitle(paste0(ROD_ident_res))
hirestiff(paste0(ROD_prefixPCres,"_UMAP","_by_","cluster","_hires.tiff"))
lowrestiff(paste0(ROD_prefixPCres,"_UMAP","_by_","cluster","_lowres.tiff"))

# UMAP of cells in each cluster by treatment without cluster labels
DimPlot(s_rods, reduction = "umap", label = FALSE, split.by = "treatment", pt.size = 1)  + NoLegend()
hirestiff(paste0(ROD_prefixPCres,"_UMAP","_by_","treatment","_no_labels","_hires.tiff"))
lowrestiff(paste0(ROD_prefixPCres,"_UMAP","_by_","treatment","_no_labels","_lowres.tiff"))

# UMAP of cells in each cluster by treatment with cluster labels
DimPlot(s_rods, reduction = "umap", label = TRUE, 
        split.by = "treatment", pt.size = 1, label.size = 6)  + NoLegend()
hirestiff(paste0(ROD_prefixPCres,"_UMAP","_by_","treatment","_with_clusters_labels","_hires.tiff"))
lowrestiff(paste0(ROD_prefixPCres,"_UMAP","_by_","treatment","_with_clusters_labels","_lowres.tiff"))

# UMAP of cells in each cluster by batch without cluster labels
DimPlot(s_rods, reduction = "umap", label = FALSE, split.by = "batch", pt.size = 1) + 
  NoLegend() + 
  ggtitle("Batch")
hirestiff(paste0(ROD_prefixPCres,"_UMAP","_by_","batch","_no_labels","_hires.tiff"))
lowrestiff(paste0(ROD_prefixPCres,"_UMAP","_by_","batch","_no_labels","_lowres.tiff"))

# UMAP of cells in each cluster by batch with cluster labels
DimPlot(s_rods, reduction = "umap", label = TRUE, 
        split.by = "batch", pt.size = 1, label.size = 6) + 
  NoLegend() + 
  ggtitle("Batch")
hirestiff(paste0(ROD_prefixPCres,"_UMAP","_by_","batch","_with_clusters_labels","_hires.tiff"))
lowrestiff(paste0(ROD_prefixPCres,"_UMAP","_by_","batch","_with_clusters_labels","_lowres.tiff"))

# UMAP of cells in each cluster by time_point without cluster labels
DimPlot(s_rods, reduction = "umap", label = FALSE, split.by = "time_point", pt.size = 1)  + NoLegend()
hirestiff(paste0(ROD_prefixPCres,"_UMAP","_by_","time_point","_no_labels","_hires.tiff"))
lowrestiff(paste0(ROD_prefixPCres,"_UMAP","_by_","time_point","_no_labels","_lowres.tiff"))

# UMAP of cells in each cluster by time_point with cluster labels
DimPlot(s_rods, reduction = "umap", label = TRUE, 
        split.by = "time_point", pt.size = 1, label.size = 6)  + NoLegend()
hirestiff(paste0(ROD_prefixPCres,"_UMAP","_by_","time_point","_with_clusters_labels","_hires.tiff"))
lowrestiff(paste0(ROD_prefixPCres,"_UMAP","_by_","time_point","_with_clusters_labels","_lowres.tiff"))

# UMAP of cells in each cluster by time_point without cluster labels
DimPlot(s_rods, reduction = "umap", label = FALSE, split.by = "time_point_treatment", pt.size = 1)  + NoLegend()
hirestiff(paste0(ROD_prefixPCres,"_UMAP","_by_","time_point_treatment","_no_labels","_hires.tiff"))
lowrestiff(paste0(ROD_prefixPCres,"_UMAP","_by_","time_point_treatment","_no_labels","_lowres.tiff"))

# UMAP of cells in each cluster by time_point with cluster labels
DimPlot(s_rods, reduction = "umap", label = TRUE, 
        split.by = "time_point_treatment", pt.size = 1, label.size = 6)  + NoLegend()
hirestiff(paste0(ROD_prefixPCres,"_UMAP","_by_","time_point_treatment","_with_clusters_labels","_hires.tiff"))
lowrestiff(paste0(ROD_prefixPCres,"_UMAP","_by_","time_point_treatment","_with_clusters_labels","_lowres.tiff"))

  #### sub Extract number of cells per cluster per orig.ident ####
n_cells <- FetchData(s_rods, vars = c("ident", "orig.ident")) %>%
  dplyr::count(ident, orig.ident) %>%
  tidyr::spread(ident, n)
n_cells
write.csv(n_cells, file = paste0(ROD_prefixPCres,"_cells_per_cluster.csv"))

rm(n_cells)

  #### sub save RDS containing reduction and cluster idents ####
# saveRDS(s_rods, paste0(ROD_prefixPCres,"_seurat_after_clustering.rds"))

# s_rods <- readRDS(file = paste0(ROD_prefixPCres,"_seurat_after_clustering.rds"))

  #### sub normalize rna slot ####

# Select the RNA counts slot to be the default assay for visualization purposes
DefaultAssay(s_rods) <- "RNA"

# Normalize, find variable features, scale data 
s_rods <- NormalizeData(s_rods)
s_rods <- FindVariableFeatures(s_rods)
all.genes <- rownames(s_rods)
s_rods <- ScaleData(s_rods, features = all.genes)

# Export normalized counts
# norm_counts <- GetAssayData(s_rods, slot = "data")
# save(norm_counts, file = paste0(ROD_prefixPCres,"_normalized_counts_RNA_sparse_matrix.rdata"))

# write.csv(norm_counts, file = paste0(ROD_prefixPCres,"_normalized_counts_RNA.csv"))

# save RDS containing RNA normalized data
# saveRDS(s_rods, paste0(ROD_prefixPCres,"_seurat_after_RNAnorm.rds"))
# s_rods <- readRDS(file = paste0(ROD_prefixPCres,"_seurat_after_RNAnorm.rds"))

  #### Plots ####

# PLOTS

VlnPlot(s_rods,
        features = Rods,
        stack = TRUE,
        flip = TRUE) +
  NoLegend() +
  ggtitle("Rod Photoreceptors")

hirestiff(paste(ROD_prefixPCres,"VLN","Rod","genes","hires.tiff", sep = "_"))
lowrestiff(paste(ROD_prefixPCres,"VLN","Rod","genes","lowres.tiff", sep = "_"))

FeaturePlot(s_rods, 
            reduction = "umap", 
            features = c(Rods1), 
            order = TRUE,
            min.cutoff = 'q10',
            label = TRUE,
            label.size = 5,
            pt.size = 0.8)

hirestiff(paste(ROD_prefixPCres,"FTR","Rod","genes","1","hires.tiff", sep = "_"))
lowrestiff(paste(ROD_prefixPCres,"FTR","Rod","genes","1","lowres.tiff", sep = "_"))

FeaturePlot(s_rods, 
            reduction = "umap", 
            features = c(Rods2), 
            order = TRUE,
            min.cutoff = 'q10',
            label = TRUE,
            label.size = 5,
            pt.size = 0.8)

hirestiff(paste(ROD_prefixPCres,"FTR","Rod","genes","2","hires.tiff", sep = "_"))
lowrestiff(paste(ROD_prefixPCres,"FTR","Rod","genes","2","lowres.tiff", sep = "_"))

  #### DEGs ####

    #### Find Cluster Markers for ####

# find all markers
cluster_markers <- FindAllMarkers(s_rods, only.pos = TRUE, min.pct = 0.50, logfc.threshold = 0.4)

write.csv(cluster_markers, 
          file = paste(ROD_prefixPCres,"clustermarkers_pos_only_minpct0.50_logfcthresh0.4.csv", sep = "_"))

    #### find DEGs between treatments ####

Idents(object = s_rods) <- "treatment"

cluster_markers <- FindAllMarkers(s_rods, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.4)

write.csv(cluster_markers, file = 
            paste(ROD_prefixPCres,"degs_treatment_all_cells_minpct0.25_logfcthresh0.4.csv", sep = "_"))

Idents(object = s_rods) <- ident_res

    #### find DEGs between time points ####

Idents(object = s_rods) <- "time_point"

cluster_markers <- FindAllMarkers(s_rods, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.4)

write.csv(cluster_markers, file = 
            paste(ROD_prefixPCres,"degs_timepoint_all_cells_minpct0.25_logfcthresh0.4.csv", sep = "_"))

Idents(object = s_rods) <- ident_res

    #### find DEGs between time points and treatment ####

Idents(object = s_rods) <- "time_point_treatment"

cluster_markers <- FindAllMarkers(s_rods, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.4)

write.csv(cluster_markers, file = 
            paste(ROD_prefixPCres,"degs_timepoint_treatment_all_cells_minpct0.25_logfcthresh0.4.csv", sep = "_"))

Idents(object = s_rods) <- ident_res

#### Cone subset recluster ####
s_cones <- subset(seurat_sct,
                 idents = c("CONE9"))

  #### SCT subset ####
DefaultAssay(s_cones) <- "RNA"

s_cones <- SCTransform(s_cones,
                      vars.to.regress = c("percent.mt","percent.ribo"), 
                      verbose = TRUE)  


  #### update project variable to include subset ####
project <- "Rat_retina_CONES"

  #### sub PCA ####
# this performs PCA on the seurat object
s_cones <- RunPCA(s_cones, npcs = 50, verbose = TRUE)

# make PC coordinate object a data frame
xx.coord <- as.data.frame(s_cones@reductions$pca@cell.embeddings)

# make PC feature loadings object a data frame
xx.gload <- as.data.frame(s_cones@reductions$pca@feature.loadings)

# calculate eigenvalues for arrays

# generate squares of all sample coordinates
sq.xx.coord <- as.data.frame(xx.coord^2)
# create empty list for eigenvalues first
eig <- c()
# calculate the eigenvalue for each PC in sq.xx.coord by taking the sqrt of the sum of squares
for(i in 1:ncol(sq.xx.coord))
  eig[i] = sqrt(sum(sq.xx.coord[,i]))
# calculate the total variance by adding up all the eigenvalues
sum.eig <- sum(eig)
# calculate the expected contribution of all PCs if they all contribute equally to the total variance
expected.contribution <- sum.eig/(length(xx.coord)-1)
# return the number of principal components with an eigenvalue greater than expected by equal variance
meaningful.PCs <- sum(eig > expected.contribution)

# create empty list for eigenvalue percentage
eig.percent <- c()
# calculate the percentage of the total variance by each PC eigenvalue
for(i in 1:length(eig))
  eig.percent[i] = 100*eig[i]/sum.eig
# sum of all eig.percent should total to 100
sum(eig.percent)
# create empty list for scree values
scree <- c()
# calculate a running total of variance contribution
for(i in 1:length(eig))
  if(i == 1) scree[i] = eig.percent[i] else scree[i] = scree[i-1] + eig.percent[i]

# create data frame for eigenvalue summaries
eigenvalues <- data.frame("PC" = colnames(xx.coord), "eig" = eig, "percent" = eig.percent, "scree" = scree)

# write csv for eigenvalues
# write.csv(eigenvalues, file = paste0("./",date,"_",project,"_PCA_eigenvalues.csv"), row.names = F)

# plot scree values
plot(eigenvalues$percent, ylim = c(0,100), type = "S", xlab = "PC", ylab = "Percent of variance",
     main = paste0(date,"_",project," scree plot all samples PCA"))
points(eigenvalues$scree, ylim = c(0,100), type = "p", pch = 16)
lines(eigenvalues$scree)
# add red line to indicate cut-off
cut.off <- 100/(length(eig)-1)
abline(h = cut.off, col = "red")
# add blue line to indicate which PCs are meaningful and kept
abline(v = meaningful.PCs, col = "blue")
text(meaningful.PCs, cut.off, label = paste("cutoff PC",meaningful.PCs),
     adj = c(-0.1, -0.5))

dev.copy(pdf, paste0("./",date,"_",project,"_scree_plot.pdf"))
dev.off()

rm(eigenvalues,sq.xx.coord,xx.coord,xx.gload,cut.off,
   eig,eig.percent,expected.contribution,i,scree,sum.eig)

# meaningful.PCs <- 10

  #### sub Run UMAP and look at UMAP plots ####

s_cones <- RunUMAP(s_cones, 
                  reduction = "pca", 
                  dims = 1:meaningful.PCs, 
                  verbose = TRUE)

# update prefixed variable
prefixPC <- paste0("./",date,"_",project,"_",meaningful.PCs,"PCs")

## UMAP plot by sample name ("orig.ident")
DimPlot(s_cones, reduction = "umap", label = FALSE, 
        pt.size = .25, group.by = "orig.ident")
hirestiff(paste0(prefixPC,"_UMAP_by_sample_hires.tiff"))
lowrestiff(paste0(prefixPC,"_UMAP_by_sample_lowres.tiff"))

DimPlot(s_cones, reduction = "umap", label = FALSE, 
        pt.size = .25, group.by = "treatment")
hirestiff(paste0(prefixPC,"_UMAP_by_treatment_hires.tiff"))
lowrestiff(paste0(prefixPC,"_UMAP_by_treatment_lowres.tiff"))

DimPlot(s_cones, reduction = "umap", label = FALSE, 
        pt.size = .25, split.by = "treatment", group.by = "treatment")
hirestiff(paste0(prefixPC,"_UMAP_split_by_treatment_hires.tiff"))
lowrestiff(paste0(prefixPC,"_UMAP_split_by_treatment_lowres.tiff"))

DimPlot(s_cones, reduction = "umap", label = FALSE, 
        pt.size = .25, group.by = "batch")
hirestiff(paste0(prefixPC,"_UMAP_by_batch_hires.tiff"))
lowrestiff(paste0(prefixPC,"_UMAP_by_batch_lowres.tiff"))

DimPlot(s_cones, reduction = "umap", label = FALSE, 
        pt.size = .25, group.by = "batch", split.by = "batch")
hirestiff(paste0(prefixPC,"_UMAP_split_by_batch_hires.tiff"))
lowrestiff(paste0(prefixPC,"_UMAP_split_by_batch_lowres.tiff"))

DimPlot(s_cones, reduction = "umap", label = FALSE, 
        pt.size = .25, group.by = "time_point")
hirestiff(paste0(prefixPC,"_UMAP_by_time_point_hires.tiff"))
lowrestiff(paste0(prefixPC,"_UMAP_by_time_point_lowres.tiff"))

DimPlot(s_cones, reduction = "umap", label = FALSE, 
        pt.size = .25, split.by = "time_point", group.by = "time_point")
hirestiff(paste0(prefixPC,"_UMAP_split_by_time_point_hires.tiff"))
lowrestiff(paste0(prefixPC,"_UMAP_split_by_time_point_lowres.tiff"))

DimPlot(s_cones, reduction = "umap", label = FALSE, 
        pt.size = .25, group.by = "time_point_treatment")
hirestiff(paste0(prefixPC,"_UMAP_by_time_point_treatment_hires.tiff"))
lowrestiff(paste0(prefixPC,"_UMAP_by_time_point_treatment_lowres.tiff"))

DimPlot(s_cones, reduction = "umap", label = FALSE, 
        pt.size = 1, split.by = "time_point_treatment", 
        group.by = "time_point_treatment")
hirestiff(paste0(prefixPC,"_UMAP_split_by_time_point_treatment_hires.tiff"))
lowrestiff(paste0(prefixPC,"_UMAP_split_by_time_point_treatment_lowres.tiff"))

DimPlot(s_cones, reduction = "umap", label = FALSE, pt.size = .25, group.by = "time_point_batch")
hirestiff(paste0(prefixPC,"_UMAP_by_time_point_batch_hires.tiff"))
lowrestiff(paste0(prefixPC,"_UMAP_by_time_point_batch_lowres.tiff"))

DimPlot(s_cones, reduction = "umap", label = FALSE, pt.size = .25, split.by = "time_point_batch", group.by = "time_point_batch")
hirestiff(paste0(prefixPC,"_UMAP_split_by_time_point_batch_hires.tiff"))
lowrestiff(paste0(prefixPC,"_UMAP_split_by_time_point_batch_lowres.tiff"))

DimPlot(s_cones, reduction = "umap", label = FALSE, pt.size = .25, group.by = "treatment_batch")
hirestiff(paste0(prefixPC,"_UMAP_by_treatment_batch_hires.tiff"))
lowrestiff(paste0(prefixPC,"_UMAP_by_treatment_batch_lowres.tiff"))

DimPlot(s_cones, reduction = "umap", label = FALSE, pt.size = .25, split.by = "treatment_batch", group.by = "treatment_batch")
hirestiff(paste0(prefixPC,"_UMAP_split_by_treatment_batch_hires.tiff"))
lowrestiff(paste0(prefixPC,"_UMAP_split_by_treatment_batch_lowres.tiff"))

DimPlot(s_cones, reduction = "umap", label = FALSE, pt.size = .25, group.by = "time_point_treatment_batch")
hirestiff(paste0(prefixPC,"_UMAP_by_time_point_treatment_batch_hires.tiff"))
lowrestiff(paste0(prefixPC,"_UMAP_by_time_point_treatment_batch_lowres.tiff"))

DimPlot(s_cones, reduction = "umap", label = FALSE, pt.size = .25, split.by = "time_point_treatment_batch", group.by = "time_point_treatment_batch")
hirestiff(paste0(prefixPC,"_UMAP_split_by_time_point_treatment_batch_hires.tiff"))
lowrestiff(paste0(prefixPC,"_UMAP_split_by_time_point_treatment_batch_lowres.tiff"))

# saveRDS(s_cones, file = paste0(prefixPC,"_seurat_integrated_preclustering.rds"))

# s_cones <- read_rds(file = paste0(prefixPC,"_seurat_integrated_preclustering.rds"))

  #### sub Clustering and Resolution ####
# DefaultAssay(s_cones) <- "integrated"

# Determine the K-nearest neighbor graph
s_cones <- FindNeighbors(object = s_cones, reduction = "pca", dims = 1:meaningful.PCs)

# Determine the clusters                              
s_cones <- FindClusters(object = s_cones,
                       resolution = c(0.1,0.2,0.3,0.4,0.5))
view(s_cones@meta.data)
res <- "_res.0.1"
ident_res <- paste0("SCT_snn",res)

Idents(s_cones) <- ident_res
DimPlot(s_cones, reduction = "umap", label = TRUE, label.size = 5, pt.size = 0.8) + 
  NoLegend() + 
  ggtitle(paste0(ident_res))

#update prefix
prefixPCres <- paste0(prefixPC,res)

# Plot the UMAP
DimPlot(s_cones, reduction = "umap", label = TRUE, label.size = 6, pt.size = 0.8)  + 
  NoLegend() + 
  ggtitle(paste0(ident_res))
hirestiff(paste0(prefixPCres,"_UMAP","_by_","cluster","_hires.tiff"))
lowrestiff(paste0(prefixPCres,"_UMAP","_by_","cluster","_lowres.tiff"))

# UMAP of cells in each cluster by treatment without cluster labels
DimPlot(s_cones, reduction = "umap", label = FALSE, split.by = "treatment", pt.size = 0.8)  + NoLegend()
hirestiff(paste0(prefixPCres,"_UMAP","_by_","treatment","_no_labels","_hires.tiff"))
lowrestiff(paste0(prefixPCres,"_UMAP","_by_","treatment","_no_labels","_lowres.tiff"))

# UMAP of cells in each cluster by treatment with cluster labels
DimPlot(s_cones, reduction = "umap", label = TRUE, 
        split.by = "treatment", pt.size = 0.8, label.size = 5)  + NoLegend()
hirestiff(paste0(prefixPCres,"_UMAP","_by_","treatment","_with_clusters_labels","_hires.tiff"))
lowrestiff(paste0(prefixPCres,"_UMAP","_by_","treatment","_with_clusters_labels","_lowres.tiff"))

# UMAP of cells in each cluster by batch without cluster labels
DimPlot(s_cones, reduction = "umap", label = FALSE, split.by = "batch", pt.size = 0.8)  + NoLegend()
hirestiff(paste0(prefixPCres,"_UMAP","_by_","batch","_no_labels","_hires.tiff"))
lowrestiff(paste0(prefixPCres,"_UMAP","_by_","batch","_no_labels","_lowres.tiff"))

# UMAP of cells in each cluster by batch with cluster labels
DimPlot(s_cones, reduction = "umap", label = TRUE, 
        split.by = "batch", pt.size = 0.8, label.size = 5)  + NoLegend()
hirestiff(paste0(prefixPCres,"_UMAP","_by_","batch","_with_clusters_labels","_hires.tiff"))
lowrestiff(paste0(prefixPCres,"_UMAP","_by_","batch","_with_clusters_labels","_lowres.tiff"))

# UMAP of cells in each cluster by time_point without cluster labels
DimPlot(s_cones, reduction = "umap", label = FALSE, split.by = "time_point", pt.size = 0.8)  + NoLegend()
hirestiff(paste0(prefixPCres,"_UMAP","_by_","time_point","_no_labels","_hires.tiff"))
lowrestiff(paste0(prefixPCres,"_UMAP","_by_","time_point","_no_labels","_lowres.tiff"))

# UMAP of cells in each cluster by time_point with cluster labels
DimPlot(s_cones, reduction = "umap", label = TRUE, 
        split.by = "time_point", pt.size = 0.8, label.size = 5)  + NoLegend()
hirestiff(paste0(prefixPCres,"_UMAP","_by_","time_point","_with_clusters_labels","_hires.tiff"))
lowrestiff(paste0(prefixPCres,"_UMAP","_by_","time_point","_with_clusters_labels","_lowres.tiff"))

# UMAP of cells in each cluster by time_point without cluster labels
DimPlot(s_cones, reduction = "umap", label = FALSE, split.by = "time_point_treatment", pt.size = 0.8)  + NoLegend()
hirestiff(paste0(prefixPCres,"_UMAP","_by_","time_point_treatment","_no_labels","_hires.tiff"))
lowrestiff(paste0(prefixPCres,"_UMAP","_by_","time_point_treatment","_no_labels","_lowres.tiff"))

# UMAP of cells in each cluster by time_point with cluster labels
DimPlot(s_cones, reduction = "umap", label = TRUE, 
        split.by = "time_point_treatment", pt.size = 0.8, label.size = 5)  + NoLegend()
hirestiff(paste0(prefixPCres,"_UMAP","_by_","time_point_treatment","_with_clusters_labels","_hires.tiff"))
lowrestiff(paste0(prefixPCres,"_UMAP","_by_","time_point_treatment","_with_clusters_labels","_lowres.tiff"))

  #### sub Extract number of cells per cluster per orig.ident ####
n_cells <- FetchData(s_cones, vars = c("ident", "orig.ident")) %>%
  dplyr::count(ident, orig.ident) %>%
  tidyr::spread(ident, n)

write.csv(n_cells, file = paste0(prefixPCres,"_cells_per_cluster.csv"))

rm(n_cells)

  #### sub save RDS containing reduction and cluster idents ####
# saveRDS(s_cones, paste0(prefixPCres,"_seurat_after_clustering.rds"))

# s_cones <- readRDS(file = paste0(prefixPCres,"_seurat_after_clustering.rds"))

  #### sub normalize rna slot ####

# Select the RNA counts slot to be the default assay for visualization purposes
DefaultAssay(s_cones) <- "RNA"

# Normalize, find variable features, scale data 
s_cones <- NormalizeData(s_cones)
s_cones <- FindVariableFeatures(s_cones)
all.genes <- rownames(s_cones)
s_cones <- ScaleData(s_cones, features = all.genes)

# Export normalized counts
# norm_counts <- GetAssayData(s_cones, slot = "data")
# save(norm_counts, file = paste0(prefixPCres,"_normalized_counts_RNA_sparse_matrix.rdata"))

# write.csv(norm_counts, file = paste0(prefixPCres,"_normalized_counts_RNA.csv"))

# save RDS containing RNA normalized data
# saveRDS(s_cones, paste0(prefixPCres,"_seurat_after_RNAnorm.rds"))
# s_cones <- readRDS(file = paste0(prefixPCres,"_seurat_after_RNAnorm.rds"))

# save RDS containing RNA normalized data AND CELLTYPE LABELS
saveRDS(s_cones, paste0(prefixPCres,"_seurat_after_RNAnorm.rds"))
s_cones <- readRDS(file = paste0(prefixPCres,"_seurat_after_RNAnorm.rds"))


#### MG subset recluster ####
s_mg <- subset(seurat_sct,
                  idents = c("MG5","MG7","MG13","MG21",
                             "MG22","MG26","MG28","MG29"))

  #### SCT subset ####
DefaultAssay(s_mg) <- "RNA"

s_mg <- SCTransform(s_mg,
                       vars.to.regress = c("percent.mt","percent.ribo"), 
                       verbose = TRUE)  


  #### update project variable to include subset ####
project <- "Rat_retina_MG"

  #### sub PCA ####
# this performs PCA on the seurat object
s_mg <- RunPCA(s_mg, npcs = 50, verbose = TRUE)

# make PC coordinate object a data frame
xx.coord <- as.data.frame(s_mg@reductions$pca@cell.embeddings)

# make PC feature loadings object a data frame
xx.gload <- as.data.frame(s_mg@reductions$pca@feature.loadings)

# calculate eigenvalues for arrays

# generate squares of all sample coordinates
sq.xx.coord <- as.data.frame(xx.coord^2)
# create empty list for eigenvalues first
eig <- c()
# calculate the eigenvalue for each PC in sq.xx.coord by taking the sqrt of the sum of squares
for(i in 1:ncol(sq.xx.coord))
  eig[i] = sqrt(sum(sq.xx.coord[,i]))
# calculate the total variance by adding up all the eigenvalues
sum.eig <- sum(eig)
# calculate the expected contribution of all PCs if they all contribute equally to the total variance
expected.contribution <- sum.eig/(length(xx.coord)-1)
# return the number of principal components with an eigenvalue greater than expected by equal variance
meaningful.PCs <- sum(eig > expected.contribution)

# create empty list for eigenvalue percentage
eig.percent <- c()
# calculate the percentage of the total variance by each PC eigenvalue
for(i in 1:length(eig))
  eig.percent[i] = 100*eig[i]/sum.eig
# sum of all eig.percent should total to 100
sum(eig.percent)
# create empty list for scree values
scree <- c()
# calculate a running total of variance contribution
for(i in 1:length(eig))
  if(i == 1) scree[i] = eig.percent[i] else scree[i] = scree[i-1] + eig.percent[i]

# create data frame for eigenvalue summaries
eigenvalues <- data.frame("PC" = colnames(xx.coord), "eig" = eig, "percent" = eig.percent, "scree" = scree)

# write csv for eigenvalues
# write.csv(eigenvalues, file = paste0("./",date,"_",project,"_PCA_eigenvalues.csv"), row.names = F)

# plot scree values
plot(eigenvalues$percent, ylim = c(0,100), type = "S", xlab = "PC", ylab = "Percent of variance",
     main = paste0(date,"_",project," scree plot all samples PCA"))
points(eigenvalues$scree, ylim = c(0,100), type = "p", pch = 16)
lines(eigenvalues$scree)
# add red line to indicate cut-off
cut.off <- 100/(length(eig)-1)
abline(h = cut.off, col = "red")
# add blue line to indicate which PCs are meaningful and kept
abline(v = meaningful.PCs, col = "blue")
text(meaningful.PCs, cut.off, label = paste("cutoff PC",meaningful.PCs),
     adj = c(-0.1, -0.5))

dev.copy(pdf, paste0("./",date,"_",project,"_scree_plot.pdf"))
dev.off()

rm(eigenvalues,sq.xx.coord,xx.coord,xx.gload,cut.off,
   eig,eig.percent,expected.contribution,i,scree,sum.eig)

# meaningful.PCs <- 16

  #### sub Run UMAP and look at UMAP plots ####

s_mg <- RunUMAP(s_mg, 
                   reduction = "pca", 
                   dims = 1:meaningful.PCs, 
                   verbose = TRUE)

# update prefixed variable
prefixPC <- paste0("./",date,"_",project,"_",meaningful.PCs,"PCs")

## UMAP plot by sample name ("orig.ident")
DimPlot(s_mg, reduction = "umap", label = FALSE, 
        pt.size = .25, group.by = "orig.ident")
hirestiff(paste0(prefixPC,"_UMAP_by_sample_hires.tiff"))
lowrestiff(paste0(prefixPC,"_UMAP_by_sample_lowres.tiff"))

DimPlot(s_mg, reduction = "umap", label = FALSE, 
        pt.size = .25, group.by = "treatment")
hirestiff(paste0(prefixPC,"_UMAP_by_treatment_hires.tiff"))
lowrestiff(paste0(prefixPC,"_UMAP_by_treatment_lowres.tiff"))

DimPlot(s_mg, reduction = "umap", label = FALSE, 
        pt.size = .25, split.by = "treatment", group.by = "treatment")
hirestiff(paste0(prefixPC,"_UMAP_split_by_treatment_hires.tiff"))
lowrestiff(paste0(prefixPC,"_UMAP_split_by_treatment_lowres.tiff"))

DimPlot(s_mg, reduction = "umap", label = FALSE, 
        pt.size = .25, group.by = "batch")
hirestiff(paste0(prefixPC,"_UMAP_by_batch_hires.tiff"))
lowrestiff(paste0(prefixPC,"_UMAP_by_batch_lowres.tiff"))

DimPlot(s_mg, reduction = "umap", label = FALSE, 
        pt.size = .25, group.by = "batch", split.by = "batch")
hirestiff(paste0(prefixPC,"_UMAP_split_by_batch_hires.tiff"))
lowrestiff(paste0(prefixPC,"_UMAP_split_by_batch_lowres.tiff"))

DimPlot(s_mg, reduction = "umap", label = FALSE, 
        pt.size = .25, group.by = "time_point")
hirestiff(paste0(prefixPC,"_UMAP_by_time_point_hires.tiff"))
lowrestiff(paste0(prefixPC,"_UMAP_by_time_point_lowres.tiff"))

DimPlot(s_mg, reduction = "umap", label = FALSE, 
        pt.size = .25, split.by = "time_point", group.by = "time_point")
hirestiff(paste0(prefixPC,"_UMAP_split_by_time_point_hires.tiff"))
lowrestiff(paste0(prefixPC,"_UMAP_split_by_time_point_lowres.tiff"))

DimPlot(s_mg, reduction = "umap", label = FALSE, 
        pt.size = .25, group.by = "time_point_treatment")
hirestiff(paste0(prefixPC,"_UMAP_by_time_point_treatment_hires.tiff"))
lowrestiff(paste0(prefixPC,"_UMAP_by_time_point_treatment_lowres.tiff"))

DimPlot(s_mg, reduction = "umap", label = FALSE, 
        pt.size = 1, split.by = "time_point_treatment", 
        group.by = "time_point_treatment")
hirestiff(paste0(prefixPC,"_UMAP_split_by_time_point_treatment_hires.tiff"))
lowrestiff(paste0(prefixPC,"_UMAP_split_by_time_point_treatment_lowres.tiff"))

DimPlot(s_mg, reduction = "umap", label = FALSE, pt.size = .25, group.by = "time_point_batch")
hirestiff(paste0(prefixPC,"_UMAP_by_time_point_batch_hires.tiff"))
lowrestiff(paste0(prefixPC,"_UMAP_by_time_point_batch_lowres.tiff"))

DimPlot(s_mg, reduction = "umap", label = FALSE, pt.size = .25, split.by = "time_point_batch", group.by = "time_point_batch")
hirestiff(paste0(prefixPC,"_UMAP_split_by_time_point_batch_hires.tiff"))
lowrestiff(paste0(prefixPC,"_UMAP_split_by_time_point_batch_lowres.tiff"))

DimPlot(s_mg, reduction = "umap", label = FALSE, pt.size = .25, group.by = "treatment_batch")
hirestiff(paste0(prefixPC,"_UMAP_by_treatment_batch_hires.tiff"))
lowrestiff(paste0(prefixPC,"_UMAP_by_treatment_batch_lowres.tiff"))

DimPlot(s_mg, reduction = "umap", label = FALSE, pt.size = .25, split.by = "treatment_batch", group.by = "treatment_batch")
hirestiff(paste0(prefixPC,"_UMAP_split_by_treatment_batch_hires.tiff"))
lowrestiff(paste0(prefixPC,"_UMAP_split_by_treatment_batch_lowres.tiff"))

DimPlot(s_mg, reduction = "umap", label = FALSE, pt.size = .25, group.by = "time_point_treatment_batch")
hirestiff(paste0(prefixPC,"_UMAP_by_time_point_treatment_batch_hires.tiff"))
lowrestiff(paste0(prefixPC,"_UMAP_by_time_point_treatment_batch_lowres.tiff"))

DimPlot(s_mg, reduction = "umap", label = FALSE, pt.size = .25, split.by = "time_point_treatment_batch", group.by = "time_point_treatment_batch")
hirestiff(paste0(prefixPC,"_UMAP_split_by_time_point_treatment_batch_hires.tiff"))
lowrestiff(paste0(prefixPC,"_UMAP_split_by_time_point_treatment_batch_lowres.tiff"))

# saveRDS(s_mg, file = paste0(prefixPC,"_seurat_integrated_preclustering.rds"))

# s_mg <- read_rds(file = paste0(prefixPC,"_seurat_integrated_preclustering.rds"))

  #### sub Clustering and Resolution ####
# DefaultAssay(s_mg) <- "integrated"

# Determine the K-nearest neighbor graph
s_mg <- FindNeighbors(object = s_mg, reduction = "pca", dims = 1:meaningful.PCs)

# Determine the clusters                              
s_mg <- FindClusters(object = s_mg,
                        resolution = c(0.1,0.2,0.3,0.4,0.5))
view(s_mg@meta.data)
res <- "_res.0.4"
ident_res <- paste0("SCT_snn",res)

Idents(s_mg) <- ident_res
DimPlot(s_mg, reduction = "umap", label = TRUE, label.size = 5, pt.size = 0.8) + 
  NoLegend() + 
  ggtitle(paste0(ident_res))

#update prefix
prefixPCres <- paste0(prefixPC,res)

# Plot the UMAP
DimPlot(s_mg, reduction = "umap", label = TRUE, label.size = 6, pt.size = 0.8)  + 
  NoLegend() + 
  ggtitle(paste0(ident_res))
hirestiff(paste0(prefixPCres,"_UMAP","_by_","cluster","_hires.tiff"))
lowrestiff(paste0(prefixPCres,"_UMAP","_by_","cluster","_lowres.tiff"))

# UMAP of cells in each cluster by treatment without cluster labels
DimPlot(s_mg, reduction = "umap", label = FALSE, split.by = "treatment", pt.size = 0.8)  + NoLegend()
hirestiff(paste0(prefixPCres,"_UMAP","_by_","treatment","_no_labels","_hires.tiff"))
lowrestiff(paste0(prefixPCres,"_UMAP","_by_","treatment","_no_labels","_lowres.tiff"))

# UMAP of cells in each cluster by treatment with cluster labels
DimPlot(s_mg, reduction = "umap", label = TRUE, 
        split.by = "treatment", pt.size = 0.8, label.size = 5)  + NoLegend()
hirestiff(paste0(prefixPCres,"_UMAP","_by_","treatment","_with_clusters_labels","_hires.tiff"))
lowrestiff(paste0(prefixPCres,"_UMAP","_by_","treatment","_with_clusters_labels","_lowres.tiff"))

# UMAP of cells in each cluster by batch without cluster labels
DimPlot(s_mg, reduction = "umap", label = FALSE, split.by = "batch", pt.size = 0.8)  + NoLegend()
hirestiff(paste0(prefixPCres,"_UMAP","_by_","batch","_no_labels","_hires.tiff"))
lowrestiff(paste0(prefixPCres,"_UMAP","_by_","batch","_no_labels","_lowres.tiff"))

# UMAP of cells in each cluster by batch with cluster labels
DimPlot(s_mg, reduction = "umap", label = TRUE, 
        split.by = "batch", pt.size = 0.8, label.size = 5)  + NoLegend()
hirestiff(paste0(prefixPCres,"_UMAP","_by_","batch","_with_clusters_labels","_hires.tiff"))
lowrestiff(paste0(prefixPCres,"_UMAP","_by_","batch","_with_clusters_labels","_lowres.tiff"))

# UMAP of cells in each cluster by time_point without cluster labels
DimPlot(s_mg, reduction = "umap", label = FALSE, split.by = "time_point", pt.size = 0.8)  + NoLegend()
hirestiff(paste0(prefixPCres,"_UMAP","_by_","time_point","_no_labels","_hires.tiff"))
lowrestiff(paste0(prefixPCres,"_UMAP","_by_","time_point","_no_labels","_lowres.tiff"))

# UMAP of cells in each cluster by time_point with cluster labels
DimPlot(s_mg, reduction = "umap", label = TRUE, 
        split.by = "time_point", pt.size = 0.8, label.size = 5)  + NoLegend()
hirestiff(paste0(prefixPCres,"_UMAP","_by_","time_point","_with_clusters_labels","_hires.tiff"))
lowrestiff(paste0(prefixPCres,"_UMAP","_by_","time_point","_with_clusters_labels","_lowres.tiff"))

# UMAP of cells in each cluster by time_point without cluster labels
DimPlot(s_mg, reduction = "umap", label = FALSE, split.by = "time_point_treatment", pt.size = 0.8)  + NoLegend()
hirestiff(paste0(prefixPCres,"_UMAP","_by_","time_point_treatment","_no_labels","_hires.tiff"))
lowrestiff(paste0(prefixPCres,"_UMAP","_by_","time_point_treatment","_no_labels","_lowres.tiff"))

# UMAP of cells in each cluster by time_point with cluster labels
DimPlot(s_mg, reduction = "umap", label = TRUE, 
        split.by = "time_point_treatment", pt.size = 0.8, label.size = 5)  + NoLegend()
hirestiff(paste0(prefixPCres,"_UMAP","_by_","time_point_treatment","_with_clusters_labels","_hires.tiff"))
lowrestiff(paste0(prefixPCres,"_UMAP","_by_","time_point_treatment","_with_clusters_labels","_lowres.tiff"))

  #### sub Extract number of cells per cluster per orig.ident ####
n_cells <- FetchData(s_mg, vars = c("ident", "orig.ident")) %>%
  dplyr::count(ident, orig.ident) %>%
  tidyr::spread(ident, n)

write.csv(n_cells, file = paste0(prefixPCres,"_cells_per_cluster.csv"))

rm(n_cells)

  #### sub save RDS containing reduction and cluster idents ####
# saveRDS(s_mg, paste0(prefixPCres,"_seurat_after_clustering.rds"))

# s_mg <- readRDS(file = paste0(prefixPCres,"_seurat_after_clustering.rds"))

  #### sub normalize rna slot ####

# Select the RNA counts slot to be the default assay for visualization purposes
DefaultAssay(s_mg) <- "RNA"

# Normalize, find variable features, scale data 
s_mg <- NormalizeData(s_mg)
s_mg <- FindVariableFeatures(s_mg)
all.genes <- rownames(s_mg)
s_mg <- ScaleData(s_mg, features = all.genes)

# Export normalized counts
# norm_counts <- GetAssayData(s_mg, slot = "data")
# save(norm_counts, file = paste0(prefixPCres,"_normalized_counts_RNA_sparse_matrix.rdata"))

# write.csv(norm_counts, file = paste0(prefixPCres,"_normalized_counts_RNA.csv"))

# save RDS containing RNA normalized data
# saveRDS(s_mg, paste0(prefixPCres,"_seurat_after_RNAnorm.rds"))
# s_mg <- readRDS(file = paste0(prefixPCres,"_seurat_after_RNAnorm.rds"))

# save RDS containing RNA normalized data AND CELLTYPE LABELS
saveRDS(s_mg, paste0(prefixPCres,"_seurat_after_RNAnorm.rds"))

s_mg <- readRDS(file = paste0(prefixPCres,"_seurat_after_RNAnorm.rds"))


  #### VLN plots ####
# Complement
VlnPlot(s_mg, 
        features = c("C1qa","C1qb",
                     "C1qc","C3ar1","C3",
                     "Itgax"), 
        flip = TRUE, 
        stack = TRUE) +
  NoLegend() +
  ggtitle("Complement Genes")

hirestiff(paste(prefixPCres,"vln","complement","genes","hires.tiff", sep = "_"))
lowrestiff(paste(prefixPCres,"vln","complement","genes","lowres.tiff", sep = "_"))

VlnPlot(s_mg, 
        features = c("C1qa","C1qb",
                     "C1qc","C3ar1","C3",
                     "Itgax"), 
        flip = TRUE, 
        stack = TRUE,
        split.by = "treatment") +
  ggtitle("Complement Genes by Treatment")

hirestiff(paste(prefixPCres,"vln","complement","genes","by","treatment","hires.tiff", sep = "_"))
lowrestiff(paste(prefixPCres,"vln","complement","genes","by","treatment","lowres.tiff", sep = "_"))

# ARMS
VlnPlot(s_mg, 
        features = c("Apoe","Cd74","RT1-Bb",
                     "RT1-Ba","Ctsb","Ctsd"), 
        flip = TRUE, 
        stack = TRUE) +
  NoLegend() +
  ggtitle("ARMS")

hirestiff(paste(prefixPCres,"vln","ARMS","genes","hires.tiff", sep = "_"))
lowrestiff(paste(prefixPCres,"vln","ARMS","genes","lowres.tiff", sep = "_"))

VlnPlot(s_mg, 
        features = c("Apoe","Cd74","RT1-Bb",
                     "RT1-Ba","Ctsb","Ctsd"), 
        flip = TRUE, 
        stack = TRUE,
        split.by = "treatment") +
  ggtitle("ARMS by Treatment")

hirestiff(paste(prefixPCres,"vln","ARMS","genes","by","treatment","hires.tiff", sep = "_"))
lowrestiff(paste(prefixPCres,"vln","ARMS","genes","by","treatment","lowres.tiff", sep = "_"))

#IRMS
VlnPlot(s_mg, 
        features = c("Mnda","Ifi27l2b","Irf7",
                     "Isg15","Oasl2","Rnf213",
                     "Rtp4","Usp18","Xaf1"), 
        flip = TRUE, 
        stack = TRUE) +
  NoLegend() +
  ggtitle("IRMS")

hirestiff(paste(prefixPCres,"vln","IRMS","genes","hires.tiff", sep = "_"))
lowrestiff(paste(prefixPCres,"vln","IRMS","genes","lowres.tiff", sep = "_"))

VlnPlot(s_mg, 
        features = c("Mnda","Ifi27l2b","Irf7",
                     "Isg15","Oasl2","Rnf213",
                     "Rtp4","Usp18","Xaf1"), 
        flip = TRUE, 
        stack = TRUE,
        split.by = "treatment") +
  ggtitle("IRMS by Treatment")

hirestiff(paste(prefixPCres,"vln","IRMS","genes","by","treatment","hires.tiff", sep = "_"))
lowrestiff(paste(prefixPCres,"vln","IRMS","genes","by","treatment","lowres.tiff", sep = "_"))

# Homeostatic_MG Genes
VlnPlot(s_mg, 
        features = c("Cx3cr1","P2ry12","Tmem119","Hexb","Cst3"), 
        flip = TRUE, 
        stack = TRUE) +
  NoLegend() +
  ggtitle("Homeostatic MG Genes")

hirestiff(paste(prefixPCres,"vln","Homeostatic_MG","genes","hires.tiff", sep = "_"))
lowrestiff(paste(prefixPCres,"vln","Homeostatic_MG","genes","lowres.tiff", sep = "_"))

VlnPlot(s_mg, 
        features = c("Cx3cr1","P2ry12","Tmem119","Hexb","Cst3"), 
        flip = TRUE, 
        stack = TRUE,
        split.by = "treatment") +
  ggtitle("Homeostatic MG Genes by Treatment")

hirestiff(paste(prefixPCres,"vln","Homeostatic_MG","genes","by","treatment","hires.tiff", sep = "_"))
lowrestiff(paste(prefixPCres,"vln","Homeostatic_MG","genes","by","treatment","lowres.tiff", sep = "_"))

VlnPlot(s_mg, 
        features = c("P2ry12","Tmem119"), 
        flip = TRUE, 
        stack = TRUE,
        split.by = "treatment") +
  ggtitle("P2ry12 & Tmem119 MG Genes by Treatment")

hirestiff(paste(prefixPCres,"vln","P2ry12_Tmem119","genes","by","treatment","hires.tiff", sep = "_"))
lowrestiff(paste(prefixPCres,"vln","P2ry12_Tmem119","genes","by","treatment","lowres.tiff", sep = "_"))


# stage 1 DAM upreg
VlnPlot(s_mg, 
        features = c("Tyrobp","Apoe","B2m","Trem2"), 
        flip = TRUE, 
        stack = TRUE) +
  NoLegend() +
  ggtitle("DAM Stage 1, Upregulated Genes")

hirestiff(paste(prefixPCres,"vln","DAM_stage_1_up","genes","hires.tiff", sep = "_"))
lowrestiff(paste(prefixPCres,"vln","DAM_stage_1_up","genes","lowres.tiff", sep = "_"))

VlnPlot(s_mg, 
        features = c("Tyrobp","Apoe","B2m","Trem2"), 
        flip = TRUE, 
        stack = TRUE,
        split.by = "treatment") +
  ggtitle("DAM Stage 1, Upregulated Genes by Treatment")

hirestiff(paste(prefixPCres,"vln","DAM_stage_1_up","genes","by","treatment","hires.tiff", sep = "_"))
lowrestiff(paste(prefixPCres,"vln","DAM_stage_1_up","genes","by","treatment","lowres.tiff", sep = "_"))

# stage 2 DAM upreg
VlnPlot(s_mg, 
        features = c("Lpl","Cst7","Axl","Itgax",
                     "Spp1","Cd9","Ccl6","Csf1"), 
        flip = TRUE, 
        stack = TRUE) +
  NoLegend() +
  ggtitle("DAM Stage 2, Upregulated Genes")

hirestiff(paste(prefixPCres,"vln","DAM_stage_2_up","genes","hires.tiff", sep = "_"))
lowrestiff(paste(prefixPCres,"vln","DAM_stage_2_up","genes","lowres.tiff", sep = "_"))

VlnPlot(s_mg, 
        features = c("Lpl","Cst7","Axl","Itgax",
                     "Spp1","Cd9","Ccl6","Csf1"), 
        flip = TRUE, 
        stack = TRUE,
        split.by = "treatment") +
  ggtitle("DAM Stage 2, Upregulated Genes by Treatment")

hirestiff(paste(prefixPCres,"vln","DAM_stage_2_up","genes","by","treatment","hires.tiff", sep = "_"))
lowrestiff(paste(prefixPCres,"vln","DAM_stage_2_up","genes","by","treatment","lowres.tiff", sep = "_"))

#Aging MG
VlnPlot(s_mg, 
        features = c("Ccl4","Il1b","Ifitm3","Rtp4","Irf7",
                     "Isg15","Oasl2","F13a1","Mrc1","Pf4",
                     "Clec12a","Ms4a7"), 
        flip = TRUE, 
        stack = TRUE) +
  NoLegend() +
  ggtitle("Aging MG Genes")

hirestiff(paste(prefixPCres,"vln","Aging_MG","genes","hires.tiff", sep = "_"))
lowrestiff(paste(prefixPCres,"vln","Aging_MG","genes","lowres.tiff", sep = "_"))

VlnPlot(s_mg, 
        features = c("Ccl4","Il1b","Ifitm3","Rtp4","Irf7",
                     "Isg15","Oasl2","F13a1","Mrc1","Pf4",
                     "Clec12a","Ms4a7"), 
        flip = TRUE, 
        stack = TRUE,
        split.by = "treatment") + 
  ggtitle("Aging MG Genes by Treatment")

hirestiff(paste(prefixPCres,"vln","Aging_MG","genes","by","treatment","hires.tiff", sep = "_"))
lowrestiff(paste(prefixPCres,"vln","Aging_MG","genes","by","treatment","lowres.tiff", sep = "_"))

# Tlr4 related
VlnPlot(s_mg, 
        features = c("Tlr4","Irak1","Irak4","Traf6",
                     "Chuk","Mapk14","Jun","Junb","Jund"), 
        flip = TRUE, 
        stack = TRUE) +
  NoLegend() +
  ggtitle("TLR4-Related Genes")

hirestiff(paste(prefixPCres,"vln","TLR4_related","genes","hires.tiff", sep = "_"))
lowrestiff(paste(prefixPCres,"vln","TLR4_related","genes","lowres.tiff", sep = "_"))

VlnPlot(s_mg, 
        features = c("Tlr4","Irak1","Irak4","Traf6",
                     "Chuk","Mapk14","Jun","Junb","Jund"), 
        flip = TRUE, 
        stack = TRUE,
        split.by = "treatment") +
  ggtitle("TLR4-Related Genes by Treatment")

hirestiff(paste(prefixPCres,"vln","TLR4_related","genes","by","treatment","hires.tiff", sep = "_"))
lowrestiff(paste(prefixPCres,"vln","TLR4_related","genes","by","treatment","lowres.tiff", sep = "_"))

# M2
VlnPlot(s_mg, 
        features = c("Cd163","Msr1","Mrc1","Vegfa"), 
        flip = TRUE, 
        stack = TRUE) +
  NoLegend() +
  ggtitle("M2 Genes")

hirestiff(paste(prefixPCres,"vln","M2","genes","hires.tiff", sep = "_"))
lowrestiff(paste(prefixPCres,"vln","M2","genes","lowres.tiff", sep = "_"))

VlnPlot(s_mg, 
        features = c("Cd163","Msr1","Mrc1","Vegfa"), 
        flip = TRUE, 
        stack = TRUE,
        split.by = "treatment") +
  ggtitle("M2 Genes by Treatment")

hirestiff(paste(prefixPCres,"vln","M2","genes","by","treatment","hires.tiff", sep = "_"))
lowrestiff(paste(prefixPCres,"vln","M2","genes","by","treatment","lowres.tiff", sep = "_"))

# M1
VlnPlot(s_mg, 
        features = c("RT1-Db1","Stat1","Cd14",
                     "Fcgr2a","Fcgr2b","Cd40","Cd86",
                     "Il1b","Il6","Tnf",
                     "Ccl2","Nos1"), 
        flip = TRUE, 
        stack = TRUE) +
  NoLegend() +
  ggtitle("M1 Genes")

hirestiff(paste(prefixPCres,"vln","M1","genes","hires.tiff", sep = "_"))
lowrestiff(paste(prefixPCres,"vln","M1","genes","lowres.tiff", sep = "_"))

VlnPlot(s_mg, 
        features = c("RT1-Db1","Stat1","Cd14",
                     "Fcgr2a","Fcgr2b","Cd40","Cd86",
                     "Il1b","Il6","Tnf",
                     "Ccl2","Nos1"), 
        flip = TRUE, 
        stack = TRUE,
        split.by = "treatment") +
  ggtitle("M1 Genes by Treatment")

hirestiff(paste(prefixPCres,"vln","M1","genes","by","treatment","hires.tiff", sep = "_"))
lowrestiff(paste(prefixPCres,"vln","M1","genes","by","treatment","lowres.tiff", sep = "_"))

# Hammond-Stevens
VlnPlot(s_mg, 
        features = c("F13a1","Mgl2"), 
        flip = TRUE, 
        stack = TRUE) +
  NoLegend() +
  ggtitle("Stevens F13A1 & Mgl2 Expression")

hirestiff(paste(prefixPCres,"vln","Stevens_F13a1_Mgl2","genes","hires.tiff", sep = "_"))
lowrestiff(paste(prefixPCres,"vln","Stevens_F13a1_Mgl2","genes","lowres.tiff", sep = "_"))

VlnPlot(s_mg, 
        features = c("F13a1","Mgl2"), 
        flip = TRUE, 
        stack = TRUE,
        split.by = "Type",
        cols = aging_palette) +
  ggtitle("Stevens F13A1 & Mgl2 Expression by Type")

hirestiff(paste(prefixPCres,"vln","Stevens_F13a1_Mgl2","genes","by","type","hires.tiff", sep = "_"))
lowrestiff(paste(prefixPCres,"vln","Stevens_F13a1_Mgl2","genes","by","type","lowres.tiff", sep = "_"))

# # monomacs
# VlnPlot(s_mg, 
#         features = c("Cd44","Ptprc","Siglec1","Mrc1","Ccr2"), 
#         flip = TRUE, 
#         stack = TRUE,
#         split.by = "Type",
#         cols = aging_palette)
# 
# VlnPlot(s_mg, 
#         features = c("Cd44","Ptprc","Siglec1","Mrc1","Ccr2",
#                      "F10","Emilin2","F5","C3","Gda","Mki67",
#                      "Sell","Hp"), 
#         flip = TRUE, 
#         stack = TRUE,
#         split.by = "Type",
#         cols = aging_palette)

# Serum Amyloid Genes and downstream
VlnPlot(s_mg, 
        features = c("Saa3",
                     "Il6","Tnf"), 
        flip = TRUE, 
        stack = TRUE) +
  NoLegend() +
  ggtitle("Serum Amyloid Genes")

hirestiff(paste(prefixPCres,"vln","saa","genes","hires.tiff", sep = "_"))
lowrestiff(paste(prefixPCres,"vln","saa","genes","lowres.tiff", sep = "_"))

VlnPlot(s_mg, 
        features = c("Saa3",
                     "Il6","Tnf"), 
        flip = TRUE, 
        stack = TRUE,
        split.by = "Type",
        cols = aging_palette) +
  ggtitle("Serum Amyloid Genes by Type")

hirestiff(paste(prefixPCres,"vln","saa","genes","by","type","hires.tiff", sep = "_"))
lowrestiff(paste(prefixPCres,"vln","saa","genes","by","type","lowres.tiff", sep = "_"))