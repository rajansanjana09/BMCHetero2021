source("R:/RESRoberts/Bioinformatics/Analysis/Sanjana/scSeurat.R")
source("R:/RESRoberts/Bioinformatics/Analysis/Sanjana/Downstream.v2.R")
#loading libraries
library(Seurat)
library(future)
library(ggplot2)
# library(fgsea)
library(msigdbr)
library(clusterProfiler)
library(SingleR)
library(dplyr) # %<%
library(pheatmap)
library(RColorBrewer)
library(viridis) # inferno color palette
library(grid)

# Create Seurat objects and perform initial QC.  Label original source.
cx.raw <- tenXLoadQC("R:/RESRoberts/Bioinformatics/scRNAOuts/S0016-zymo/filtered_feature_bc_matrix/", spec = "mixHuman")
cx.raw <- subset(cx.raw, subset = nFeature_RNA >3500 & nCount_RNA <50000 & percent.mt <15)
cx.raw$src <- "Culture"

tib.raw <- tenXLoadQC("R:/RESRoberts/Bioinformatics/scRNAOuts/S0018xS0028/filtered_feature_bc_matrix/", spec = "mixHuman")
tib.raw <- subset(tib.raw, subset = nFeature_RNA >3000 & nCount_RNA <60000 & percent.mt <18)
tib.raw$src <- "Tibia"

lung.raw <- tenXLoadQC("R:/RESRoberts/Bioinformatics/scRNAOuts/S0024xS0029/filtered_feature_bc_matrix/", spec = "mixHuman")
lung.raw <- subset(lung.raw, subset = nFeature_RNA >1250 & nCount_RNA <60000 & percent.mt <25)
lung.raw$src <- "Lung"

# Add lineage tracing tags to the Seurat objects
cx.raw <- processLTBC(cx.raw,
                      lt.loc = "R:/RESRoberts/Bioinformatics/scRNAOuts/S0016-zymo/lt.fq",
                      cid.loc = "R:/RESRoberts/Bioinformatics/scRNAOuts/S0016-zymo/cid.fq",
                      histogram = T,
                      title = "Culture, Top 40 Clones",
                      ymax = 1.75,
                      col.fill = "#E64B35FF",
                      relative = T)
tib.raw <- processLTBC(tib.raw,
                       lt.loc = "R:/RESRoberts/Bioinformatics/scRNAOuts/S0018xS0028/lt.fq",
                       cid.loc = "R:/RESRoberts/Bioinformatics/scRNAOuts/S0018xS0028/cid.fq",
                       histogram = T,
                       title = "Tibia, Top 40 Clones",
                       ymax = 1.75,
                       col.fill = "#00A087FF",
                       relative = T)
lung.raw <- processLTBC(lung.raw,
                        lt.loc = "R:/RESRoberts/Bioinformatics/scRNAOuts/S0024xS0029/lt.fq",
                        cid.loc = "R:/RESRoberts/Bioinformatics/scRNAOuts/S0024XS0029/cid.fq",
                        histogram = T,
                        title = "Lineage Barcode Enrichment by Microenvironment",
                        ymax = 1.75,
                        relative = T)


#######Export Lineage Tag enrichment plots
pdf("Culture_LT.pdf", width = 7, height = 2)
cx.raw$fig
dev.off()

pdf("Tibia_LT.pdf", width = 7, height = 2)
tib.raw$fig
dev.off()

pdf("Lung_LT.pdf", width = 7, height = 2)
lung.raw$fig
dev.off()

cx.raw <- cx.raw$sobject
tib.raw <- tib.raw$sobject
lung.raw <- lung.raw$sobject

#Subset all to #2800 cells in each condition 
cx.raw <- subset(cx.raw, cells = sample(Cells(cx.raw), 2800))
tib.raw <- subset(tib.raw, cells = sample(Cells(tib.raw), 2800))
lung.raw <- subset(lung.raw, cells = sample(Cells(lung.raw), 2800))

# Merge into a single Seurat object
os17 <- merge(cx.raw, y = c(tib.raw, lung.raw),
              add.cell.ids = c("Culture", "Tibia", "Lung"),
              project = "LineageTracing")

# Process and cluster
os17 <- NormalizeData(os17) %>%
  FindVariableFeatures(selection.method = "vst") %>%
  ScaleData() %>%
  RunPCA(pc.genes = os.17@var.genes, npcs = 20) %>%
  RunUMAP(reduction = "pca", dims = 1:20) %>%
  FindNeighbors(reduction = "pca", dims = 1:20) %>%
  FindClusters(resolution = 0.3)

# CCR
# Attempt to regress out the effects of cell cycle on these tumor cells
s.genes <- cc.genes$s.genes
g2m.genes <- cc.genes$g2m.genes
os17 <- CellCycleScoring(object = os17, s.features = s.genes,
                         g2m.features = g2m.genes, set.ident = TRUE)

os17 <- ScaleData(object = os17, vars.to.regress = c("S.Score", "G2M.Score"),
                  features = rownames(x = os17)) 
os17 <- RunPCA(os17, pc.genes = os.17@var.genes, npcs = 20) %>%
  RunUMAP(reduction = "pca", dims = 1:20) %>%
  FindNeighbors(reduction = "pca", dims = 1:20) %>%
  FindClusters(resolution = 0.3)

# os17 <- RunPCA(os17, pc.genes = os17@var.genes, npcs = 20)
# os17 <- RunHarmony(os17, group.by.vars = "src", plot_convergence = T)
# os17 <- RunUMAP(os17, reduction = "harmony", dims = 1:20)
# os17 <- FindNeighbors(os17, reduction = "harmony", dims = 1:20)
# os17 <- FindClusters(os17, resolution = 0.3)

# pdf("Harmony.pdf", width = 7, height = 7)
# DimPlot(os17, reduction = "umap", group.by = "src", pt.size = 1, label = F, order = cell.ids) + 
#   coord_fixed() + 
#   ggtitle("OS17 by Source") + 
#   scale_color_npg()
# dev.off()

# Plot the data 
set.seed(100)
cell.ids <- sample(colnames(os17))
DimPlot(os17, reduction = "umap", group.by = "src", pt.size = 1, label = F, order = cell.ids) + 
  coord_fixed() + 
  ggtitle("OS17 by Source") + 
  scale_color_npg()
DimPlot(os17, reduction = "umap", pt.size = 1, label = T) + 
  coord_fixed() + 
  ggtitle("OS17 Clusters") + 
  scale_color_npg(alpha = 0.7)

##################Find the optimum resolution for clustering
os17.nC <- nRes(os17, 
                res = seq(from = 0.2, to = 0.3, by = 0.01))
plot <- pSil(os17.nC, 0.3)
os17.3 <- os17 %>% FindClusters(resolution = 0.3)
os17.2 <- os17 %>% FindClusters(resolution = 0.2)

# find markers for every cluster compared to all remaining cells, report only the positive ones
os17.markers <- FindAllMarkers(os17.2, only.pos = FALSE, min.pct = 0.25, logfc.threshold = 0.25)

#Heatmap
os17.markers %>%
  group_by(cluster) %>%
  top_n(n = 5, wt = avg_log2FC) -> top10

# pdf("All.pdf", width = 10, height = 7)
DoHeatmap(os17.2, features = top10$gene) + NoLegend()
# dev.off()

#Markers
VlnPlot(os17.3,
        features = c("COL1A1", "LUM", "RUNX2", "SOX9", "CTSK", "JUN", "CD74", "SFRP2"),
        slot = "counts",
        ncol = 4,
        log = TRUE)