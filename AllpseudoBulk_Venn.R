# Psudobulk to determine gene differnetially regulated on tissue colonization and 
# Venn diagram to visualize shared genes

# # Merge into a single Seurat object
# OS <- merge(OB, y = c(os17.cx.raw, os17.tib.raw, os17.lung.raw, 
#                       t143b.cx.raw, t143b.tib.raw, t143b.lung.raw, 
#                       OS2.cx.raw, OS2.tib.raw, OS2.lung.raw, 
#                       OS7.cx.raw, OS7.tib.raw, OS7.lung.raw),
#             add.cell.ids = c("OB", "OS17_Culture", "OS17_Tibia", "OS17_Lung", 
#                              "t143b_Culture", "t143b_Tibia", "t143b_Lung",
#                              "NCHOS2_Flank", "NCHOS2_Tibia", "NCHOS2_Lung", 
#                              "NCHOS7_Flank", "NCHOS7_Tibia", "NCHOS7_Lung"),
#             project = "Heterogeneity")

#OS17

os17 <- subset(OS, cells = WhichCells(OS, expression = cond == "OS17"))
t143B <- subset(OS, cells = WhichCells(OS, expression = cond == "t143b"))
NCHOS2 <- subset(OS, cells = WhichCells(OS, expression = cond == "NCHOS2"))
NCHOS7 <- subset(OS, cells = WhichCells(OS, expression = cond == "NCHOS7"))
  
#extract markers
data <- NCHOS7
DimPlot(data, group.by = "src")

#rename Idents to src type before DGE analysis
Idents(data, cells = WhichCells(data, expression = src == "NCHOS7_Flank")) <- "Culture"
Idents(data, cells = WhichCells(data, expression = src == "NCHOS7_Lung")) <- "Lung"
Idents(data, cells = WhichCells(data, expression = src == "NCHOS7_Tibia")) <- "Tibia"
data[["tissue"]] <- Idents(data)

NCHOS7 <- data
DimPlot(NCHOS7)

data <- list(
  os17 = os17,
  t143B = t143B,
  NCHOS2 = NCHOS2, 
  NCHOS7 = NCHOS7
)

#Lung marker - UP
Lung.up <- list(NULL)
for (i in 1:length(data)) {
  tmp <- data[[i]]
  x <- FindMarkers(tmp, ident.1 = "Lung", ident.2 = "Culture", only.pos = TRUE, min.pct = 0.25)
  Lung.up[[i]] <- rownames(x)
  }

#Tibia markers - UP 
Tibia.up <- list(NULL)
for (i in 1:length(data)) {
  tmp <- data[[i]]
  x <- FindMarkers(tmp, ident.1 = "Tibia", ident.2 = "Culture", only.pos = TRUE, min.pct = 0.25)
  Tibia.up[[i]] <- rownames(x)
}

#Lung markers - DOWN
Lung.down <- list(NULL)
for (i in 1:length(data)) {
  tmp <- data[[i]]
  x <- FindMarkers(tmp, ident.1 = "Lung", ident.2 = "Culture", only.pos = FALSE, min.pct = 0.25)
  x <- x[x$avg_log2FC<0,]
  Lung.down[[i]] <- rownames(x)
}

#Tibia markers - DOWN 
Tibia.down <- list(NULL)
for (i in 1:length(data)) {
  tmp <- data[[i]]
  x <- FindMarkers(tmp, ident.1 = "Tibia", ident.2 = "Culture", only.pos = FALSE, min.pct = 0.25)
  x <- x[x$avg_log2FC<0,]
  Tibia.down[[i]] <- rownames(x)
}

#UpDown genes
TC.markers <- FindMarkers(data, ident.1 = "Tibia", ident.2 = "Culture", only.pos = FALSE, min.pct = 0.25)
LC.markers <- FindMarkers(data, ident.1 = "Lung", ident.2 = "Culture", only.pos = FALSE, min.pct = 0.25)

#Separate pos and neg avg_log2FC
TC.dn.markers <- TC.markers[TC.markers$avg_log2FC<0,]
TC.up.markers <- TC.markers[TC.markers$avg_log2FC>0,]

LC.dn.markers <- LC.markers[LC.markers$avg_log2FC<0,]
LC.up.markers <- LC.markers[LC.markers$avg_log2FC>0,]

# view results
write.xlsx(TC.dn.markers, file = "R:/RESRoberts/Bioinformatics/Analysis/Sanjana/2020/2.TC.dn.markers.xlsx", 
           col.names = TRUE, row.names = TRUE, append = FALSE)
write.xlsx(TC.up.markers, file = "R:/RESRoberts/Bioinformatics/Analysis/Sanjana/2020/1.TC.up.markers.xlsx", 
           col.names = TRUE, row.names = TRUE, append = FALSE)

write.xlsx(LC.dn.markers, file = "R:/RESRoberts/Bioinformatics/Analysis/Sanjana/2020/4.LC.dn.markers.xlsx", 
           col.names = TRUE, row.names = TRUE, append = FALSE)
write.xlsx(LC.up.markers, file = "R:/RESRoberts/Bioinformatics/Analysis/Sanjana/2020/3.LC.up.markers.xlsx", 
           col.names = TRUE, row.names = TRUE, append = FALSE)


############################Example Venn Diagram
# if (!require(devtools)) install.packages("devtools")
# devtools::install_github("yanlinlin82/ggvenn")
# library(ggvenn)

# if (!require(devtools)) install.packages("devtools")
# devtools::install_github("gaospecial/ggVennDiagram")
# library("ggVennDiagram")
genes <- Tibia.up
x <- list(
  A = genes[[1]],
  B = genes[[2]],
  C = genes[[3]],
  D = genes[[4]]
)

pdf("Tibia.up.pdf", width = 5, height = 5)
plot <- ggvenn(
  x,
  fill_color = c("#0073C2FF", "#EFC000FF", "#868686FF", "#CD534CFF"),
  stroke_size = 0.9, set_name_size = 4,  show_percentage = TRUE
)
plot
dev.off()

# ggVennDiagram(x[1:4], label_alpha = 0.7, label = "count",
#               category.names = c("OS17",
#                                  "143B",
#                                  "NCHOS2",
#                                  "NCHOS7")) + 
#   scale_fill_gradient(low = "#F4FAFE", high = "#4981BF")+
#   theme(legend.title = element_text(color = "black"),
#         legend.position = "right")

###Extract genes that are shared in these datasets
genes <- Tibia.up
A = genes[[1]]
B = genes[[2]]
C = genes[[3]]
D = genes[[4]]

# # Preparing clusterProfiler to perform hypergeometric test on msigdb signatures
# m_t2g.c2 <- msigdbr(species = "Homo sapiens", category = "C2") %>%
#   dplyr::select(gs_name, human_gene_symbol)
# m_t2g.c6 <- msigdbr(species = "Homo sapiens", category = "C6") %>%
#   dplyr::select(gs_name, human_gene_symbol)
# m_t2g.h <- msigdbr(species = "Homo sapiens", category = "H") %>%
#   dplyr::select(gs_name, human_gene_symbol)
# m_t2n.h <- msigdbr(species = "Homo sapiens", category = "H") %>% 
#   dplyr::select(gs_id, gs_name)
# m_t2g=rbind(m_t2g.c2,m_t2g.c6)
# 
# # msigdb signature to use
# msig.gene.set = m_t2g.h
# msig.name = m_t2n.h

intersect <- c((intersect(intersect(intersect(A,B),C),D)), 
               (intersect(intersect(A,B),C)),
               (intersect(intersect(A,B),D)),
               (intersect(intersect(A,D),C)),
               (intersect(intersect(B,D),C)))

tmp <- enricher(intersect, TERM2GENE=msig.gene.set, TERM2NAME = msig.name)
em=tmp@result[,c("ID", "p.adjust")] #pvalue/p.adjust
rownames(em) <- NULL
em <- em[1:10,]
em[,2] <- -log10(em[,2])
# barplot(em)

# pathways <- list(NULL)
pathways[[1]] <- em

write.xlsx(pathways, file ="R:/RESRoberts/Bioinformatics/Analysis/Sanjana/2020/pathways.xlsx", col.names = TRUE, row.names = TRUE, append = FALSE)
