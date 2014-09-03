#First, we try to find genes that predict treatment status. We use [_tspair_](http://www.bioconductor.org/packages/release/bioc/html/tspair.html) package. As we are comparing treatment/non-treatment groups and using normalized/ComBat-adjusted data, _tspair_ is particularly useful because it was designed for the binary phenotype, and is robust to normalization/processing steps.

tsp <- tspcalc(combat_edata, meta$Cohort)
tsp
# out <- tspsig(combat_edata, meta$Cohort, B=50, seed=1)

idx <- unique(c(tsp$index[,1], tsp$index[,2]))
genes.limma <- rownames(res)
genes.tsp <- rownames(combat_edata)[idx]
genes.all <- unique(c(genes.limma, genes.tsp))
res.venn <- matrix(0, nrow=length(genes.all), ncol=2)
rownames(res.venn) <- genes.all
colnames(res.venn) <- c("limma", "tsp")
res.venn[genes.limma, 1] <- 1
res.venn[genes.tsp, 2] <- 1
source("venn4.R")
vennDiagram(res.venn)

idx.tsp <- tsp$index[order(tsp$score, decreasing=T),]
head(apply(idx.tsp, 1, function(x) paste(annot.f[x[1], "GeneName"], annot.f[x[2], "GeneName"], sep=":")), n=20)
