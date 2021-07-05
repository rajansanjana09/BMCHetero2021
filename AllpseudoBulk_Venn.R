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
os17 <- subset(OS, cells = OS[["cond"]] == "OS17")

os17 <- subset(OS, cells = WhichCells(OS, expression = cond == "OS17"))
t143B <- subset(OS, cells = WhichCells(OS, expression = cond == "t143b"))
NCHOS2 <- subset(OS, cells = WhichCells(OS, expression = cond == "NCHOS2"))
NCHOS7 <- subset(OS, cells = WhichCells(OS, expression = cond == "NCHOS7"))
 
#extract markers
data <- os17
DimPlot(data, group.by = "src")

#rename Idents to src type before DGE analysis
Idents(data, cells = WhichCells(data, expression = src == "OS17_Culture")) <- "Culture"
Idents(data, cells = WhichCells(data, expression = src == "OS17_Lung")) <- "Lung"
Idents(data, cells = WhichCells(data, expression = src == "OS17_Tibia")) <- "Tibia"

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