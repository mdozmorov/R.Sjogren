library("WGCNA")
tmp<-log2(exprs)
rownames(tmp)<-seq(1,dim(tmp)[[1]])
datET<-as.matrix(tmp)
rowGroup<-annot[,2]
rowID<-rownames(tmp)
tmp1<-collapseRows(datET,rowGroup,rowID,method="MaxMean")$datETcollapsed
write.table(tmp1,"F:/WorkOMRF/Databases/ImmGen/ShannonData1.txt",sep="\t")

mtx<-read.table("F:/WorkOMRF/Databases/ImmGen/matrix.txt",sep="\t",header=T,row.names=1)
IAC<-cor(mtx)
library(cluster)
cluster1=hclust(as.dist(1-IAC),method="ward") # "ward", "single", "complete", "average", "mcquitty", "median" or "centroid".
plot(cluster1, cex=0.1) #,cex=0.7,labels=colnames(tmp))

query<-read.table("clipboard",sep="\t")
a<-query[,1]
e<-query[,2]
write.table(aggregate(e,list(a),min),"clipboard",sep="\t",row.names=F,col.names=F)
write.table(aggregate(e,list(a),max),"clipboard",sep="\t",row.names=F,col.names=F)

# Prepare the data for barplot
query<-readLines("clipboard") # read gene names
unique(annot$GeneName[annot$GeneName %in% query])
setdiff(query,annot$GeneName[annot$GeneName %in% query])
data<-exprs.quantile[annot$GeneName %in% query,] # Subset quantile normalized data by these genes
genenames<-annot$GeneName[annot$GeneName %in% query] # Get gene names
data<-data[order(genenames),] # Organize by alphabetical order the data
genenames<-genenames[order(genenames)] # and the genes

#Dirty correlation between cytokines and genes
tmp2<-read.table("clipboard",sep="\t",header=T,row.names=1) # Read cytokine data from clipboard
tmp3<-cor(t(tmp2),t(exprs.q.collapsed[query,colnames(tmp2)])) # Correlation between columns, order kept
write.table(tmp3,"clipboard",sep="\t")
heatmap.2(tmp3,trace="none",col=color,distfun=function(x){dist(x,method=dist.method)}, hclustfun=function(x){hclust(x,method=hclust.method)}, density.info="none",cexCol=1,cexRow=1) # Visualization

for (i in 1:nrow(tmp2)){
  cor.minmax<-vector("numeric",nrow(exprs.q.collapsed))
  for (j in 1:nrow(exprs.q.collapsed)){
    cor.minmax[j]<-cor(as.numeric(tmp2[i,]),as.numeric(exprs.q.collapsed[j,colnames(tmp2)]))
  }
  print(paste(rownames(tmp2)[i],
              "best correlates with",rownames(exprs.q.collapsed)[cor.minmax == max(cor.minmax)],
              "at Pearson's",formatC(as.numeric(cor.minmax[cor.minmax == max(cor.minmax)])),
              "and best anticorrelates with",rownames(exprs.q.collapsed)[cor.minmax ==min(cor.minmax)],
              "at Pearson's",formatC(as.numeric(cor.minmax[cor.minmax == min(cor.minmax)]))
              ))
}




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
