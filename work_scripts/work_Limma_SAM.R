library(gplots) # install.packages("gplot")
library(RColorBrewer) # of source("http:/bioconductor.org/biocLite.R") biocLite("RColorBrewer")
opar<-par(no.readonly=T) #Save original settings
par(oma=c(10,0,0,10)) #Make right and bottom margins larger
color<-rev(brewer.pal(11,'RdYlGn')) #Red-yellow-green gradient
color<-greenred #Standard green-black-red palette
color<-bluered # Blue-red palette

# Normalization and limma analysis
source("http://bioconductor.org/biocLite.R")
biocLite(c("Biobase"))
biocLite("illuminaHumanv4.db")

library(Biobase)
library(vsn)
library("preprocessCore")
library("illuminaHumanv4.db")
library("genefilter")
library("arrayQualityMetrics")
library("WGCNA")
## loadData
# Full annotations. To avoid R errors on special characters, in Excel select "Description" column, use Ctrl+1, Custom category, \"@\" type
annot.f <- read.table("data//annot_full.txt",sep="\t", quote="\"", header=T, as.is=T) 
annot.f <- apply(annot.f, 2, function(x){sapply(x, function(y){gsub("\"", "", y)})}) # Remove \"
colnames(annot.f) <- c("ProbeName", "GeneName", "SystematicName", "Description") # Rename columns
rownames(annot.f) <- paste("gene", seq(1, nrow(annot.f)), sep="") # Dummy row names, to be later used for merging
# Expression data
exprs <- as.matrix(read.table("data//data.txt", sep="\t", header=T))
# Meta data
meta1 <- read.table("data//meta1.txt", sep="\t", header=T, row.names=1) # Has cohort
meta2 <- read.table("data//meta2.txt", sep="\t", header=T, row.names=1) # Has other clinical parameters
meta <- rbind(meta1[1:2, intersect(colnames(meta1), colnames(meta2))], 
              meta2[, intersect(colnames(meta1), colnames(meta2))]) # Need only the first two parameters from meta1. Join vertically by common names
patients <- colnames(meta)[meta["Microarray Class", ] == "AC" | meta["Microarray Class", ] == "PSS"] # Should be 34. This gives us 8 AC and 9 PSS for Cohort 1 and 8 AC and 9 PSS for Cohort 2
patients <- patients[!patients %in% c("p1033216.2", "p1033680.6...5.")]

meta <- as.data.frame(t(meta[, patients])); exprs <- exprs[, patients] # Subsetting
colnames(meta) <- sapply(colnames(meta), function(x) gsub(" ", "", x)) # Removing spaces from column names

# Excluding patients
meta <- meta[!(meta[,"Subject.ID"] %in% c("p1033216.2", "p1033680.6...5.")),]
exprs <- exprs[, !(colnames(exprs) %in% c("p1033216.2", "p1033680.6...5."))]

meta$Microarray.Class<-factor(meta$Microarray.Class) # Remove unused levels, like "Skip", "Partial"

# Read in expression data, and annotation in the same order# PCA
summary(prcomp(log10(exprs)))
pca<-prcomp(log10(exprs))$rotation
x = pca[,1]
y = pca[,2]
xadj<-0.1*(max(x)-min(x))
yadj<-0.1*(max(y)-min(y))
plot(x,y,xlab="PC1", ylab="PC2", main="PC analysis on cohorts",
     xlim=c(min(x) - xadj, max(x) + xadj), ylim=c(min(y) - yadj, max(y) + yadj),
     pch=ifelse(meta$Cohort == 1, 1, 2),
     col=ifelse(meta$MicroarrayClass == "AC", "red", "blue")) 
text(x,y+0.03,labels=rownames(pca), cex=0.7)
legend("bottomright", c("Cohort1/AC", "Cohort2/AC", "Cohort1/PCC", "Cohort2/PCC"), 
       col=c("red","red","blue","blue"),
       pch=c(1,2,1,2))


# Normalizing
library(limma)
exprs.n<-log2(normalizeQuantiles(exprs))
colnames(exprs.n)<-colnames(exprs)
rownames(exprs.n)<-paste("gene", seq(1,nrow(exprs.n)), sep="") # Dummy row names, to be later used for merging

#arrayQualityMetrics(new("ExpressionSet",exprs=exprs),outdir="./QC_quantile") # QC
# Calculating IACs for all pairs of samples and examining the distribution of IACs in the dataset
IAC=cor(exprs[rowMeans(log2(exprs))>8,],use="p")
library(cluster)
cluster1=hclust(as.dist(1-IAC),method="average")
plot(cluster1,cex=0.7)#,labels=paste(meta$Cohort,meta$MicroarrayClass,sep="-")) #dimnames(exprs)[[2]])

# Limma and SVA/Combat
eset.0<-new("ExpressionSet",exprs=as.matrix(exprs.n)) # Make ExpressionSet

# Dealing with batch effect
library(sva)
mod<-model.matrix(~as.factor(MicroarrayClass), data=meta) # Full model matrix
colnames(mod)[2]<-"outcome"
mod0<-model.matrix(~1, data=meta) # Null model
n.sv = num.sv(exprs.n,mod,method="leek")
n.sv # Number of confounding values
svobj<-sva(exprs.n, mod, mod0)#, n.sv=2) # Artificially set number of batches to 2
modSv<-cbind(mod, svobj$sv) # Modify model matrixes
mod0Sv<-cbind(mod0, svobj$sv)
pValues = f.pvalue(exprs.n,mod,mod0) # Unadjusted p-values
qValues = p.adjust(pValues,method="BH")
sum(qValues<0.1)
pValuesSv = f.pvalue(exprs.n,modSv,mod0Sv) # Batch effect adjusted
qValuesSv = p.adjust(pValuesSv,method="BH")
sum(qValuesSv<0.1)
writeLines(unique(as.character(annot.f[qValuesSv < 0.1, "GeneName"])), "clipboard")
modBatch<-model.matrix(~as.factor(Microarray.Class) + as.factor(Cohort), data=meta) # Specify batch
mod0Batch<-model.matrix(~as.factor(Cohort), data=meta)
pValuesBatch<-f.pvalue(exprs.n, modBatch, mod0Batch)
qValuesBatch<-p.adjust(pValuesBatch, method="BH")
sum(qValuesBatch<0.1)
unique(as.character(annot.f[qValuesBatch < 0.1, "GeneName"]))

# LIMMA begins
fit <- lmFit(eset.0, mod) # Check eset!!
fit2 <- eBayes(fit)
tmp<-topTable(fit2,coef="outcome",number=nrow(exprs(eset.0)),adjust.method="BH",p.value=0.1)
dim(tmp)
writeLines(paste("Total:",nrow(tmp),"Up:",sum(tmp$logFC>0),"Dn:",sum(tmp$logFC<0), sep=" "), "clipboard")
write.table(unique(annot.f[rownames(tmp), c(2,4)]), "clipboard-128", sep="\t", row.names=F)
write.table(merge(tmp, annot.f, by="row.names"), "f:/111.txt", sep="\t", row.names=F)

# ComBat
combat_edata<-ComBat(dat=exprs.n, batch=meta$Cohort, mod=mod, numCovs=NULL, par.prior=TRUE, prior.plots=F)
pValuesComBat = f.pvalue(combat_edata,mod,mod0)
qValuesComBat = p.adjust(pValuesComBat,method="BH")
sum(qValuesComBat < 0.1)
unique(as.character(annot.f[qValuesComBat < 0.1, "GeneName"]))
eset.0<-new("ExpressionSet",exprs=as.matrix(combat_edata)) # Make ExpressionSet from combat'ed data
fit <- lmFit(eset.0, mod)
fit2 <-eBayes(fit)

# Cohorts separately
cohort=2
eset.0<-new("ExpressionSet", exprs=exprs.n[,meta$Cohort==cohort])
mod<-model.matrix(~as.factor(Microarray.Class), data=meta[meta$Cohort==cohort,])

# Between cohorts
mod<-model.matrix(~as.factor(Cohort), data=meta)
# Rename column
colnames(mod)[2]<-"outcome"


# ber
library(ber)
combat_edata<-ber(t(exprs.quantile), pheno$batch)
combat_edata<-combat_p(t(exprs.quantile), as.factor(pheno$batch), covariates = NULL)

groups<-colnames(exprs.quantile)
groups[grep("AC",(groups))]<- "AC"
groups[grep("NC",(groups))]<- "NC"
groups[grep("SS",(groups))]<- "SS"
groups
# Averaging across conditions"
AC.av<-apply(exprs.quantile[,grep("AC",groups)],1,mean)
NC.av<-apply(exprs.quantile[,grep("NC",groups)],1,mean)
SS.av<-apply(exprs.quantile[,grep("SS",groups)],1,mean)
exprs.q.av<-cbind(NC.av,AC.av,SS.av)
tmp<-collapseRows(as.matrix(exprs.q.av),annot[,2],seq(1,length(annot[,2])),"maxRowVariance")
tmp$datETcollapsed['FLJ10661',]
write.table(tmp$datETcollapsed,"data.collapsed.txt",sep="\t")
#"Not averaging, jsut reordering sequentially"
exprs.q.AC<-exprs.quantile[,grep("AC",groups)]
dim(exprs.q.AC)
exprs.q.NC<-exprs.quantile[,grep("NC",groups)]
dim(exprs.q.NC)
exprs.q.SS<-exprs.quantile[,grep("SS",groups)]
dim(exprs.q.SS)
exprs.q<-cbind(exprs.q.NC,exprs.q.AC,exprs.q.SS)
colnames(exprs.q)<-c(rep("NC",5),rep("AC",9),rep("SS",10)) # 9 AC 07/20/2013☺
rownames(exprs.q)<-rownames(exprs.quantile)
tmp1<-collapseRows(as.matrix(exprs.q),annot[,2],seq(1,length(annot[,2])),"maxRowVariance")
head(tmp1$datETcollapsed)
exprs.q.collapsed<-tmp1$datETcollapsed # Store collapsed dataset
rownames(exprs.q.collapsed)<-annot[tmp1$selectedRow,'GeneName'] # Make row names as gene names

# Get all AC-NC, SS-NC gene indexes
query<-unique(as.numeric(readLines("clipboard"))) # Gene indexes
all.degs<-collapseRows(exprs.q.av[query,],annot[query,2],query,"maxRowVariance") # Average expression of indexes, collapsed on gene names, indexes as row names, highest variance
all.degs<-exprs.q[query,] # or Take all DEGs without collapsing
head(all.degs$datETcollapsed)
dim(all.degs$datETcollapsed)
# Heatmap of these data
dist.method<-"euclidean"  # "euclidean", "maximum", "manhattan", "canberra", "binary" or "minkowski"
hclust.method<-"average" # "ward", "single", "complete", "average", "mcquitty", "median" or "centroid"
temp<-as.matrix(t(scale(t(all.degs$datETcollapsed)))) # Scaled data, collapsed dataset
temp<-as.matrix(t(scale(t(all.degs))))
# Restart subcluster processing from here
hc.rows<-hclust(dist(temp,method=dist.method),method=hclust.method) # Cluster rows
plot(hc.rows) # Plot dendrogram, then use it for heatmap
h<-heatmap.2(temp,Rowv=as.dendrogram(hc.rows), Colv=F, dendrogram="row",col=color,density.info="none", trace="none")
table(cutree(hc.rows,k=4)) # How many genes in each cluster
writeLines(names(cutree(hc.rows,k=3)[cutree(hc.rows,k=3) == 3]),"clipboard") # Gene in a specific cluster
barplot(tmp$datETcollapsed['XBP1',])
tmp$datETcollapsed['PDE4DIP',] # Cluster 3
# Individual up- down-regulated clusters
cl1<-names(cutree(hc.rows,k=4)[cutree(hc.rows,k=4) == 1]) # Genes in cluster 1
cl2<-names(cutree(hc.rows,k=4)[cutree(hc.rows,k=4) == 4]) # Genes in cluster 1
temp<-as.matrix(t(scale(t(all.degs$datETcollapsed[cl2,])))) # Change cl1 to cl2 to get scaled data, subclusters
dim(temp)
# Then go through the above heatmap procedure
# Cluster 1, subcluster 3 - stronger downregulated in AC
# Cluster 1, subcluster 2 - stronger downregulated in SS

# PhenoData
pData<-read.table("phenoData.txt", row.names=1, header=T, sep='\t')
pData<-pData[-c(grep(72,rownames(pData))),]
dim(pData)
rownames(pData)
summary(pData)
rownames(pData)<-colnames(exprs)
all(rownames(pData) == colnames(exprs)) # Must be true
names(pData)
sapply(pData,class)
# # Metadata
# metadata<-data.frame(labelDescription = c("Patient ID", "Visit ID", "Sex", "Maturity", "Race", "Actual visit week", "white blood cells", "red blood cells", "platelets"), 
#                      row.names = c("patients", "visits", "gender", "age", "etnicity", "week", "WBC", "RBC", "PLT"))
phenoData<-new("AnnotatedDataFrame", data=pData) #, varMetadata=metadata)
annotation<-"illuminaHumanv4"

eset<-new("ExpressionSet",exprs=as.matrix(exprs.quantile), phenoData=phenoData, annotation="illuminaHumanv4")
# eset.filtered<-varFilter(eset)
# dim(exprs(eset.filtered)) # 23658 genes left)
# eset.filtered.avg<-cbind(rowMeans(exprs(eset.filtered[,1:15])),rowMeans(exprs(eset.filtered[,16:30])),rowMeans(exprs(eset.filtered[,31:45])))
library(betr)
prob<-betr(eset=eset, cond=pData(eset)$strain, timepoint=pData(eset)$time, replicate=pData(eset)$rep, alpha=0.05)
write.table(prob,"prob.txt",sep='\t')

heatmap.2(as.matrix(exprs.av[prob==1,]),Rowv=T,Colv=F,trace="none",col=color,distfun=function(x){dist(x,method=dist.method)}, hclustfun=function(x){hclust(x,method=hclust.method)}, density.info="none",cexCol=1.5,cexRow=1)#, cellnote=formatC(1/10^abs(mtx), format="e", digits=2), notecol='darkgreen')


m<-exprs.av[prob==1,]
m<-exprs.av
annot<-read.table("clipboard",sep='\t',header=F,row.names=1,as.is=T)

# Plotting time dynamics for Probe_ID
timeplot <- function(probe) {
  row <- as.numeric(m[probe,]) # Get the numerical data
  names(row) <- colnames(m)  # Column names
  rMax <- max(row) + (max(row) - min(row)) * 0.1 # Y axis Max
  rMin <- max(c(min(row) - (max(row) - min(row)) * 0.1, 0)) # Y axis Min
  # print(rMin)
  if (rMin<1) {rMin<-0} # Min can't be negative'
  KO <- grep("KO", names(row))     # Get the data for spacific condition
  WT <- grep("WT", names(row))
  
  
  plot(row[KO], type='o', col='red', lwd=3, lty=2, ylim=c(rMin,rMax), 
       ylab='Log2 Expression', xlab='Time Point',xaxt='n')
  par(xaxt='s') # Add x axis to the plot
  axis(1,at=1:4,labels=as.character(c(0,6,24,48))) # Add x axis labels
  lines(row[WT], type='o', col='green', lwd=3, lty=1)
  title(as.character(c((probe),"-",annot[probe,])),font.main=1)
  legend('topright', c('KO', 'WT'), 
         lty=c(2,1), lwd=c(rep(2)), col=c('red','green',cex=0.8,horiz=FALSE))
}

m<-read.table("clipboard",header=T,row.names=1,sep='\t',as.is=T)
dev.off()
pdf("Profiles.pdf")
par(mfrow=c(3,2))                           # Makes plot containing 2x2
for (i in 1:dim(m)[[1]]){
  timeplot(rownames(m)[i])
}
dev.off()

# Time plot for gene name
timeplot.gene <- function(gene) {
  #dev.off()
  # probes <- m[m$"ILMN_GENE"==gene,"PROBE_ID"]
  probes<-rownames(annot)[annot==gene]
  sapply(probes, function(probe) {timeplot(probe); par(ask=T)})
}

geneList<-readLines("clipboard")
dev.off()
pdf("Profiles-p38 MAPKinases.pdf")
par(mfrow=c(3,2))                           # Makes plot containing 2x2
for (i in 1:length(geneList)){
  timeplot.gene(geneList[i])
}
dev.off()






# SAM begins
library(samr)
# genenames<-readLines("clipboard") # Gene names only
genenames<-annot$GeneName
# y<-paste(c(rep(1,20),rep(2,20)),"Time",c(rep(0,5),rep(6,5),rep(24,5),rep(48,5)), sep="")
y=c(rep(1,3),rep(2,4))
groups
y<-groups[grep("AC|NC",groups)]
y[y =="AC"]<-2
y[y == "NC"]<-1
exprs.0<-exprs.quantile[,grep("AC|NC",groups)] # Get only N condition

data<-list(x=exprs.0,y=y, geneid=as.character(rownames(exprs.0)),genenames=genenames,logged2=T)
samr.obj<-samr(data,resp.type="Two class unpaired",nperms=100)
delta<-1
samr.plot(samr.obj,delta)
delta.table<-samr.compute.delta.table(samr.obj)
siggenes.table<-samr.compute.siggenes.table(samr.obj,delta,data,delta.table,min.foldchange=1)
nrow(siggenes.table$genes.up)
nrow(siggenes.table$genes.lo)
selected.IDs<-siggenes.table$genes.up[,2]

write.table(siggenes.table$genes.up,"clipboard-128",sep='\t')
write.table(siggenes.table$genes.lo,"clipboard-128",sep='\t')

tmp<-siggenes.table$genes.up
# tmp<- as.data.frame(tmp)
FC = as.numeric(tmp[,'Fold Change'])
genes = toupper(tmp[,2])
write.table(aggregate(FC,by=list(genes),min),"clipboard-128",sep="\t")#,row.names=F,col.names=F)


# Paired design
patients<-as.matrix(read.table("Paired_Patients-Visits.txt",header=F,sep='\t'))[,1]
patients.f<-factor(patients)
patients.f<-phenoData(eset)$patients

visits<-as.matrix(read.table("Paired_Patients-Visits.txt",header=F,sep='\t'))[,2]
visits.f<-factor(visits,levels=c("V1","V6","VF"))
visits.f<-phenoData(eset)$visits

design<-model.matrix(~patients.f+visits.f)
fit <- lmFit(eset.filtered, design)
fit <- eBayes(fit)
topTable(fit, coef="visits.fVF")
topTable(fit, coef="visits.fV6")
results <- decideTests(fit,adjust.method="BH",p.value=0.05,lfc=log2(1.5))
ID.fVF<-topTable(fit, coef="visits.fV6", number=2000, adjust.method="BH", p.value=0.05, lfc=log2(1.5))$ID
exprs.fVF.avg<-eset.filtered.avg[ID.fVF,]


write.table(exprs.fVF, "visits_fV6.txt", sep="\t")

dist.method<-"binary"  # "euclidean", "maximum", "manhattan", "canberra", "binary" or "minkowski"
hclust.method<-"average" # "ward", "single", "complete", "average", "mcquitty", "median" or "centroid"
heatmap.2(as.matrix(2^exprs.fVF.avg[4:40,-3]), distfun=function(x){dist(x,method=dist.method)}, hclustfun=function(x){hclust(x,method=hclust.method)}, col=greenred, key=T,  keysize=0.1, density.info="none", trace="none",  cexCol=1, cexRow=1) #lwid=c(1.5,3), lhei=c(1.5,4),



# Using design matrix as described in limma guide section 8.5
design.1<-read.table("design_(-33).txt", header=T, row.names=1, sep='\t')
fit.1<-lmFit(eset,design.1)
cont.matrix.1<-makeContrasts(V6vsV1=V6-V1, VFvsV1=VF-V1, VFvsV6=VF-V6, levels=design.1)
fit.1.2<-contrasts.fit(fit.1,cont.matrix.1)
fit.1.2<-eBayes(fit.1.2)
topTable(fit.1.2, coef="V6vsV1", adjust="BH") # None
topTable(fit.1.2, coef="VFvsV1", adjust="BH") # None
topTable(fit.1.2, coef="VFvsV6", adjust="BH") # None

# Using design matrix as in limma 8.6
design.2<-model.matrix(~0+visits.f)
colnames(design.2)<-c("V1", "V6", "VF")
contrast.matrix.2<-makeContrasts(V6-V1, VF-V1, VF-V6, levels=design.2)
fit.2<-lmFit(eset,design.2)
fit.2.2<-contrasts.fit(fit.2, contrast.matrix.2)
fit.2.2<-eBayes(fit.2.2)
topTable(fit.1.2, coef=1, adjust="BH") # None
topTable(fit.1.2, coef=2, adjust="BH") # None
topTable(fit.1.2, coef=3, adjust="BH") # None

# x<-read.table("clipboard", sep='\t')
# exp <- new("ExpressionSet", exprs=as.matrix(x))
exp<-new("ExpressionSet",exprs=exprs)
exp.norm<-justvsn(exp)
meanSdPlot(exp.norm)
write.table(exprs(exp.norm),"F:/111.txt",sep="\t",row.names=F,col.names=F)

design<-cbind(WT=1,MUvsWT=c(0,0,0,0,0,0,1,1,1,1,1,1))
fit<-lmFit(exp.norm,design)
fit<-eBayes(fit)
results<-topTable(fit,coef="MUvsWT",adjust="BH",number=2000)
write.table(results,"F:/111.txt",sep="\t")
