# Create dataset collapsed by gene names
# After running Analysis.Rmd, we have batch-adjusted low variability filtered combat_edata
library("WGCNA")
datET <- as.matrix(combat_edata)
rowGroup <- annot.f[rownames(combat_edata), 2]# Gene names corresponding to row IDs
rowID <- rownames(combat_edata) # row IDs from combat_edata
combat_edata_collapsed <- collapseRows(datET, rowGroup, rowID, method="MaxMean")$datETcollapsed
write.table(combat_edata_collapsed,"results/combat_edata_collapsed.txt", sep="\t", col.names=NA)

# Heatmap of selected genes
query<-readLines("clipboard")
data<-exprs.q # Quantile normalized matrix, ordered by conditions, collapsed to gene names
data<-data[,grep("NC",colnames(data),invert=T)] # Remove NC samples
data.plot<-t(scale(t(data[annot$GeneName %in% query,])))
data.plot<-t(scale(t(data[query,])))

data<-data[,order(colnames(data))] # Sorted columns
head(data)
genenames<-annot$GeneName[annot$GeneName %in% query]
data<-exprs.q.av[annot$GeneName %in% query,] # Columns ordered as they should
data.collapsed<-collapseRows(data.plot,genenames,seq(1,dim(data.plot)[[1]]),"maxRowVariance") # Default method - MaxMean
data.plot<-as.matrix((data.collapsed$datETcollapsed))
data.plot<-data.plot[,c("AC.av", "SS.av")]

library(gplots) # install.packages("gplot")
library(RColorBrewer) # of source("http:/bioconductor.org/biocLite.R") biocLite("RColorBrewer")
opar<-par(no.readonly=T) #Save original settings
par(oma=c(5,0,0,15)) #Make right and bottom margins larger
color<-rev(brewer.pal(11,'RdYlGn')) #Red-yellow-green gradient
color<-brewer.pal(6,'Reds') #Only for overrepresentation, for positive numbers, low=white, high=red intensity
color<-greenred #Standard green-black-red palette
color<-colorRampPalette(c("blue", "white","red"))

dist.method<-"euclidean"  # "euclidean", "maximum", "manhattan", "canberra", "binary" or "minkowski"
hclust.method<-"average" # "ward", "single", "complete", "average", "mcquitty", "median" or "centroid"
heatmap.2(data.plot,trace="none",col=color, Colv=F, distfun=function(x){dist(x,method=dist.method)}, hclustfun=function(x){hclust(x,method=hclust.method)}, density.info="none",cexCol=1,cexRow=0.8,labRow=annot$GeneName[annot$GeneName %in% query])#, cellnote=formatC(1/10^abs(mtx), format="e", digits=2), notecol='darkgreen')

heatmap.2(t(scale(t(as.matrix((data.collapsed$datETcollapsed))))),trace="none",col=color, Colv=F, distfun=function(x){dist(x,method=dist.method)}, hclustfun=function(x){hclust(x,method=hclust.method)}, density.info="none",cexCol=1,cexRow=0.8,labRow=rownames(data.collapsed$datETcollapsed))#, cellnote=formatC(1/10^abs(mtx), format="e", digits=2), notecol='darkgreen')

# Expression of selected genes
rownames(annot.f["PDCD1" == annot.f[, "GeneName"], ])
dt1 <- combat_edata["gene23548", meta$Cohort == 1]
dt2 <- combat_edata["gene23548", meta$Cohort == 2]
