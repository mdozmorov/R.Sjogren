# Collapse gene names by maximum expression
mtx<-read.table("clipboard", sep="\t", header=T, row.names=1) # 3 columns matrix of gene IDs, names, and FC
mtx.a<-aggregate(mtx$logFC, by=list(mtx$GeneName), max)
write.table(mtx.a, "clipboard", sep="\t", col.names=F, row.names=F)

# Collapse whole dataset
library(WGCNA)
combat_edata_collapsed<-collapseRows(combat_edata, as.character(annot.f[,"X.ProbeName."]), rownames(combat_edata))
write.table(combat_edata_collapsed$datETcollapsed, "results/combat_edata_collapsed.txt", sep="\t", col.names=NA)
write.table(meta[meta[,"Subject.ID"] %in% colnames(combat_edata),], "results/combat_edata_meta.txt", sep="\t", col.names=NA)
