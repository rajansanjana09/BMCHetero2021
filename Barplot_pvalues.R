library(org.Hs.eg.db)
keytypes(org.Hs.eg.db)

data(geneList, package="DOSE")
gene <- names(geneList)[abs(geneList) > 2]

gmtfile <- system.file("extdata", "c5.cc.v5.0.entrez.gmt", package="clusterProfiler")
c5 <- read.gmt(gmtfile)

egmt <- enricher(gene, TERM2GENE=c5)
head(egmt)
egmt2 <- GSEA(geneList, TERM2GENE=c5, verbose=FALSE)
head(egmt2)

ggo <- groupGO(gene     = gene,
               OrgDb    = org.Hs.eg.db,
               ont      = "CC",
               level    = 3,
               readable = TRUE)

head(ggo)

barplot(tmp, drop=TRUE, showCategory=12, )
egmt@ontology


