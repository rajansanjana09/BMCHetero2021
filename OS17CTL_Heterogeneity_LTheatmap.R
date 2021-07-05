#create heatmap 
B.list <- list(Culture = cx.sub,
               Tibia = tib.sub,
               Lung = lung.sub)

#Heatmap and em.hm files
em.hm.list <- list()
for (i in 1:length(B.list)) {
  em.hm.list[[i]] <- DGEA(B.list[[i]]) 
}

library(data.table)
for(i in 1:(length(em.hm.list))) {
  em.hm.list[[i]] <- setDT(em.hm.list[[i]], keep.rownames = TRUE)[]
}

temp1 <- em.hm.list[[1]]
temp2 <- em.hm.list[[2]]

cx.pd <- full_join(temp1, temp2, by = "rn")
cx.pd <- as.data.frame(cx.pd)

rownames(cx.pd)=cx.pd[,1]
cx.pd=cx.pd[,-1]
cx.pd[is.na(cx.pd)]=1
cx.pd <- cx.pd[!grepl('NA', rownames(cx.pd)), ] #remove rows with rownames 'NA'
# cx.pd <- cx.pd[grepl("^NA", rownames(cx.pd))==F,]

# cx.pd <- cx.pd[-(23:30),] #remove rows with names NA - figure out why we have NAs!!
cx.pd.log <- -log10(cx.pd) #log transform
library(pheatmap)   

cx.pd_transpose <- as.data.frame(t(cx.pd.log))
df <- as.matrix(cx.pd_transpose)

subj<-c("A0", "A1","A2", "A3", "A4",
        "B0", "B1","B2", "B3")
rownames(df)<-subj
aka2 = data.frame(ID = factor(c("OS17","OS17", "OS17","OS17", "OS17",
                                "t143B", "t143B", "t143B", "t143B")))
rownames(aka2)<-subj
aka3 = list(ID = c(OS17 = "red", t143B = "gold"))
# 
# df[] <- lapply(df, gsub, pattern='HALLMARK_', replacement='')
# df


pheatmap(t(scale(df)),
         annotation_col = aka2, 
         annotation_colors = aka3[1],
         annotation_legend = TRUE,
         gaps_col =  4,
         show_colnames = T, show_rownames = T, cluster_rows = T, 
         cluster_cols = T, legend = TRUE, 
         clustering_distance_rows = "euclidean", border_color = FALSE)