####################subset to cells that have a lineage tag
memory.limit(size = 40000)
#Culture
Culture.LTcells <- vector(mode = "character")
for(i in 1:nrow(cx.lt.list)){
  temp <- WhichCells(cx.sub, expression = lt == cx.lt.list[[i,1]]) 
  Culture.LTcells <- append(Culture.LTcells, temp)
}
cx.tagged <- subset(cx.sub, cells = Culture.LTcells)

#Tibia
tib.LTcells <- vector(mode = "character")
for(i in 1:nrow(tib.lt.list)){
  temp <- WhichCells(tib.sub, expression = lt == tib.lt.list[[i,1]]) 
  tib.LTcells <- append(tib.LTcells, temp)
}
tib.tagged <- subset(tib.sub, cells = tib.LTcells)

#Lung
lung.LTcells <- vector(mode = "character")
for(i in 1:515){
  temp <- WhichCells(lung.sub, expression = lt == lung.lt.list[[i,1]]) 
  lung.LTcells <- append(lung.LTcells, temp)
}
lung.tagged <- subset(lung.sub, cells = lung.LTcells)

lung.sub <- subset(lung.sub, idents = c("0", "1", "2", "6")) ##Remove outlier cells in cluster 3 and 5

pdf("LungOutline.pdf", height = 5, width = 5)
DimPlot(lung.sub, reduction = "umap", pt.size = 1, label = TRUE) + 
  coord_fixed() + 
  ggtitle("Lung-derived") +
  xlim(0,9) + NoAxes()+ 
  scale_color_npg(alpha = 1)
dev.off()

#Lung
vplayout <- function(x, y) viewport(layout.pos.row = x, layout.pos.col = y)
pdf("grid_lung_v4.pdf", width = 10, height = 10)
grid.newpage()
pushViewport(viewport(layout = grid.layout(1, 4)))
for (i in 1:4) {
  p <- DimPlot(lung.tagged, 
               reduction = "umap",
               pt.size = 1,
               cells.highlight = WhichCells(lung.sub, expression = lt == lung.lt.list[[i,1]]), 
               cols.highlight = "#4DBBD5FF",
               sizes.highlight = 3) + 
    coord_fixed() +
    theme(legend.position = "none") +
    ggtitle(paste("Clone", lung.lt.list[[i,1]])) +
    xlim(0,9) + NoLegend() + NoAxes()
  print(p, vp=vplayout(ceiling(i/4), i))
}
dev.off()

#cell cycle distribution of top 4 lineages - change to Pie chart and percentage

vplayout <- function(x, y) viewport(layout.pos.row = x, layout.pos.col = y)
pdf("pie.pdf", width = 15, height = 15)
grid.newpage()
pushViewport(viewport(layout = grid.layout(1, 4)))
for(i in 1:4){
  temp <- WhichCells(lung.sub, expression = lt == lung.lt.list[[i,1]]) 
  LT.cells <- subset(lung.sub, cells = temp)
  df <- table(LT.cells$Phase)
  bp <- ggplot(as.data.frame(df), aes(x="", y=Freq, fill=Var1))+
    geom_bar(width = 1, stat = "identity")
  pie <- bp + coord_polar("y", start=0) + scale_fill_brewer(palette="Blues")+
    theme_minimal() + NoAxes() +
    ggtitle(paste("Clone", lung.lt.list[[i,1]]))
  print(pie, vp=vplayout(ceiling(i/4), i))
}
dev.off()

#barplot
bp <-  ggplot(as.data.frame(df), aes(x = Var1, y = Freq)) +
  geom_bar(fill = "#4DBBD5FF", stat = "identity") +
  ggtitle("Cell cycle distribution") +
  ylab("Count") +
  xlab("Phase")

#######################Pathways activated in top 10 lineages

lung.lt.list[1:10,1]
# 052-082154 
# 039-074146 
# 037-078792 
# 039-052640 
# 095-044696 
# 039-098388 
# 064-065939 
# 013-044371 
# 063-007191 
# 018-082025

data <- lung.tagged

#Subset cells of each lineage into a separate cluster in lung tagged
for(i in 1:10) {
  cells <- WhichCells(data, expression = lt == lung.lt.list[[i,1]])
  Idents(object = data, cells = cells) <- paste("Enriched Lineage", i, sep = " ")
}

library(forcats)
Idents(data) <- fct_rev(Idents(data))

# Idents(data) <- data$seurat_clusters
# Stash cell identity classes
data[["Lineage"]] <- Idents(data)

DimPlot(data)
#reorganize Lineage levels


#Top markers for each lineage

# find markers for every cluster compared to all remaining cells, report positive and negative ones
data.markers <- FindAllMarkers(data, only.pos = FALSE, min.pct = 0.25, logfc.threshold = 0.25)
data.markers <- FindAllMarkers(data, only.pos = TRUE)

#Run through function
em.hm <- DGEA(data)

#Remove non-LT clusters
em.hm.lt <- em.hm[-c(1:4)]
clust.ids = sort(unique(data@active.ident))
clust.ids <- clust.ids[-c(1:4)]
c <- rownames(em.hm.lt)
#remove "HALLMARK_"
c <- gsub("HALLMARK_", "", c)
#remove "_" by removing special characters
c <- gsub("_", " ", c)
rownames(em.hm.lt) <- c
#remove rows with low pvalues
row_names_df_to_remove<-c("FATTY ACID METABOLISM",
                          "ADIPOGENESIS",
                          "PEROXISOME",
                          "SPERMATOGENESIS",
                          "XENOBIOTIC METABOLISM",
                          "PI3K AKT MTOR SIGNALING",
                          "HEDGEHOG SIGNALING")
em.hm.lt <- em.hm.lt[!(row.names(em.hm.lt) %in% row_names_df_to_remove),]

col.breaks=seq(-log10(1),min(max(-log10(em.hm.lt))+1,20),by=0.5)
col=inferno(length(col.breaks)) # library(viridis)
col=c("white",colorRampPalette(brewer.pal(n = 7, name ="Reds"))(50))
pheatmap(-log10(em.hm.lt[,1:length(clust.ids)]),
         cluster_rows = TRUE,
         cluster_cols = FALSE,
         cellwidth = 20,
         cellheight = 10,
         treeheight_row = 0,
         treeheight_col=0,
         color = col,
         scale='none',
         breaks=col.breaks,
         fontsize = 11)

########################Run IPA on the DEGs from top ten lineages
#NOTE: After Enriched lineage 5, not more than 25 DEGs with P_adj_val < 0.01
#Think about if the enrichr data is valid and do we need to run a GSEA for Glycolysis alone
#How is DEG analysis done, do we filter genes based on pvalue, pr adj p value from FindMarkers output
#Do we have to re-do the analysis in CCR uncorrected space?
L1.markers <- FindMarkers(data, ident.1 = "Enriched Lineage 1", ident.2 = NULL, only.pos = FALSE, min.pct = 0.25)
L2.markers <- FindMarkers(data, ident.1 = "Enriched Lineage 2", ident.2 = NULL, only.pos = FALSE, min.pct = 0.25)
L3.markers <- FindMarkers(data, ident.1 = "Enriched Lineage 3", ident.2 = NULL, only.pos = FALSE, min.pct = 0.25)
L4.markers <- FindMarkers(data, ident.1 = "Enriched Lineage 4", ident.2 = NULL, only.pos = FALSE, min.pct = 0.25)
L5.markers <- FindMarkers(data, ident.1 = "Enriched Lineage 5", ident.2 = NULL, only.pos = FALSE, min.pct = 0.25)
L6.markers <- FindMarkers(data, ident.1 = "Enriched Lineage 6", ident.2 = NULL, only.pos = FALSE, min.pct = 0.25)
L7.markers <- FindMarkers(data, ident.1 = "Enriched Lineage 7", ident.2 = NULL, only.pos = FALSE, min.pct = 0.25)
L8.markers <- FindMarkers(data, ident.1 = "Enriched Lineage 8", ident.2 = NULL, only.pos = FALSE, min.pct = 0.25)
L9.markers <- FindMarkers(data, ident.1 = "Enriched Lineage 9", ident.2 = NULL, only.pos = FALSE, min.pct = 0.25)
L10.markers <- FindMarkers(data, ident.1 = "Enriched Lineage 10", ident.2 = NULL, only.pos = FALSE, min.pct = 0.25)

#remove unnecessary columns
myList <- list(L1.markers, L2.markers, L3.markers, L4.markers, L5.markers, L6.markers, L7.markers, L8.markers, L9.markers, L10.markers)
myList <- lapply(myList, function(x) { x["p_val"] <- NULL; x[, 2:3] <- NULL; x })

#write markers to excel sheet
names <- c("L1.markers", "L2.markers", "L3.markers", "L4.markers", "L5.markers", "L6.markers", "L7.markers", "L8.markers", "L9.markers", "L10.markers")
library(xlsx)
for (i in 1:10){
  location <- paste0("R:/RESRoberts/Bioinformatics/Analysis/Sanjana/2020/L",i,".markers.xlsx", sep = "")
  write.xlsx(myList[[i]], file = location, 
             col.names = TRUE, row.names = TRUE, append = FALSE)
} 