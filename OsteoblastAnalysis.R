############################Overlay osteoblast markers on Culture/Flank cells 
load("R:/RESRoberts/Bioinformatics/Analysis/Sanjana/2020/17.Harmony/Harmony in individual samples/Figure3/OB-S0031_CCR_v2.RData")

data <- OS.list[[1]] 
DimPlot(data, reduction = "umap", pt.size = 1, label = T) + 
  coord_fixed() + 
  scale_color_npg(alpha = 0.7)

#Set optimum resolution for clustering
data.nC <- nRes(data, 
                res = seq(from = 0.1, to = 0.3, by = 0.05))
plot <- pSil(data.nC, 0.15)
data <- data %>% FindClusters(resolution = 0.15)

###################Heatmap of OB markers
#Heatmap

ob.markers <- FindAllMarkers(OS.list[[1]], only.pos = FALSE, min.pct = 0.25)

ob.markers %>%
  group_by(cluster) %>%
  top_n(n = 10, wt = avg_log2FC) -> top10

# pdf("All.pdf", width = 15, height = 10)
p <- DoHeatmap(OS.list[[1]], features = top10$gene, size = 5) + NoLegend()
p + theme(text = element_text(size=16, color = "black"))
# dev.off()

#############################Generate phentoype module
cluster.0 <- FindMarkers(data, ident.1 = "0", only.pos = FALSE, min.pct = 0.25)
cluster.1 <- FindMarkers(data, ident.1 = "1", only.pos = FALSE, min.pct = 0.25)
cluster.2 <- FindMarkers(data, ident.1 = "2", only.pos = FALSE, min.pct = 0.25)

#Extract gene names
cluster0 <- list(rownames(cluster.0))
cluster1 <- list(rownames(cluster.1))
cluster2 <- list(rownames(cluster.2))

write.xlsx(cluster.0, file = "C:/Users/rssxr002/Downloads/cluster.0.xlsx")
write.xlsx(cluster.1, file = "C:/Users/rssxr002/Downloads/cluster.1.xlsx")
write.xlsx(cluster.2, file = "C:/Users/rssxr002/Downloads/cluster.2.xlsx")

##############Annotation of phenotype module
#Ran msigdb gene set enrichment on http://www.gsea-msigdb.org/gsea/index.jsp
#Visualization of p values usign barplot

SanjanaPlots2 <- read.delim("C:/Users/rssxr002/Downloads/SanjanaPlots2.txt")

ggplot(SanjanaPlots2, aes(x = -1 * Order,
                         y = FDR.log,
                         fill = FDR.log > 0)) +
  geom_bar(stat = "identity") +
  coord_flip() +
  facet_wrap(~ Cluster, ncol = 1) +
  scale_fill_brewer(type = "qual", palette = "Set1") +
  geom_text(aes(y = Label_y * 15, label = Pathway)) +
  theme_bw() +
  ylab("") 

#Addmodule score to OS cell culture/flank models
# OB <- OS.list[[1]]
# OS17.cx <- OS.list[[2]]
# OS17.Tibia <-  OS.list[[3]]
# OS17.Lung <-OS.list[[4]]
# t143B.cx <- OS.list[[5]]
# t143B.Tibia <-OS.list[[6]]
# t143B.Lung <- OS.list[[7]]
# OS2.Flank <- OS.list[[8]]
# OS2.Tibia <- OS.list[[9]]
# OS2.Lung <- OS.list[[10]]
# OS7.Flank <-  OS.list[[11]]
# OS7.Tibia <- OS.list[[12]]
# OS7.Lung <- OS.list[[13]]

for (i in 1:13) {
  OS.list[[i]] <- AddModuleScore(OS.list[[i]], features = cluster0, ctrl = 5, name = 'Red')
  OS.list[[i]] <- AddModuleScore(OS.list[[i]], features = cluster1, ctrl = 5, name = 'Blue')
  OS.list[[i]] <- AddModuleScore(OS.list[[i]], features = cluster2, ctrl = 5, name = 'Green')
}

plot1 <- DotPlot(OS.list[[2]], features = c("cluster01", "cluster11", "cluster21")) + coord_flip() 
plot2 <- DotPlot(OS.list[[3]], features = c("cluster01", "cluster11", "cluster21")) + coord_flip()
plot3 <- DotPlot(OS.list[[4]], features = c("cluster01", "cluster11", "cluster21")) + coord_flip()

# DotPlot(OS.list[[1]], features = c("cluster01", "cluster11", "cluster21")) + coord_flip() 

# plot1 <- RidgePlot(OS.list[[11]], features = c("cluster01", "cluster11", "cluster21"))
# plot2 <- RidgePlot(OS.list[[12]], features = c("cluster01", "cluster11", "cluster21"))
# plot3 <- RidgePlot(OS.list[[13]], features = c("cluster01", "cluster11", "cluster21"))

  CombinePlots(
    plots = list(plot1, plot2, plot3),
    legend = 'none', nrow = 1) 

  