library(gplots) # install.packages("gplot")
library(RColorBrewer) # of source("http:/bioconductor.org/biocLite.R") biocLite("RColorBrewer")
opar<-par(no.readonly=T) #Save original settings
par(oma=c(10,0,0,10)) #Make right and bottom margins larger
color<-greenred #Standard green-black-red palette
dist.method<-"euclidean"  # "euclidean", "maximum", "manhattan", "canberra", "binary" or "minkowski"
hclust.method<-"ward" # "ward", "single", "complete", "average", "mcquitty", "median" or "centroid"

library(pamr)
response.data<-pamr.from.excel("BOLD2_baseline_pamr.txt",77, sample.labels=T)
# pamr.menu(response.data) # All automatic

tmp<-readLines("clipboard") # Read in batch labels
unique(tmp) # Check unique
par(mar=c(2,2,2,2))
hist(as.numeric(tmp[idx]),breaks=length(tmp[idx]))
# lines(density(as.numeric(tmp)))
sort(table(tmp1),decreasing=T) # Frequency of each
idx<-grep("African",tmp1) # Select what to keep
response.data$x<-response.data$x[,idx] # Reformat expression matrix, idx to keep
response.data$y<-(tmp[idx]) # Reformat batch labels, idx to keep
response.data$samplelabels<-response.data$samplelabels[idx] # Reformat patient IDs, idx to keep
tmp1<-pamr.makeclasses(response.data)
tmp1<-cut(as.numeric(tmp),quantile(as.numeric(tmp)),include.lowest=T,labels=1:4) # Cut by quantiles

# Random testing
# response.data$y<-sample((response.data$y),length(response.data$y),replace=F) # Random sampling, replace=T makes variable number of labels drawn
# All manual. http://www-stat.stanford.edu/~tibs/PAM/Rdist/doc/readme.html
response.train<-pamr.train(response.data) 
# response.train
response.cv<-pamr.cv(response.train,response.data)

# response.cv
pamr.plotcv(response.cv)
new.scales<-pamr.adaptthresh(response.train)
response.train2<-pamr.train(response.data,threshold.scale=new.scales)
response.cv2<-pamr.cv(response.train2,response.data)
pamr.plotcv(response.cv2)
t<-4
pamr.plotcen(response.train,response.data,t)
pamr.confusion(response.cv,t)
pamr.plotcvprob(response.train,response.data,t)
pamr.geneplot(response.train,response.data,t)
response.genes<-pamr.listgenes(response.train,response.data,t,genenames=T)
write.table(response.genes,"clipboard",sep='\t')

fdr.obj<-pamr.fdr(response.train,response.data)
pamr.plotfdr(fdr.obj)

response.train2<-pamr.train(response.data,hetero="F")
response.results2<-pamr.cv(response.train2,response.data)
response.scales<-pamr.adaptthresh(response.train)
response.train3<-pamr.train(response.data,threshold.scale=response.scales)
response.results3<-pamr.cv(response.train3,response.data)
response.results2


# Convert condition labels into class labels 1 and -1
v <- response.data$y
v[v=="M"] <- -1
v[v=="F"] <- 1
v <- as.numeric(v)
# Data for SVM should be transposed
response.svm<-t(response.data$x)
svm.model<-svm(response.svm, v, type="C-classification", kernel="linear")
predicted<-predict(svm.model, response.svm) # predict labels of training data
sum(predicted != v) # count differences
table(true=v, pred=predicted) # confusion matrix
svm.cross <- svm(response.svm,v, type="C-classification", kernel="linear", cross=10)


# http://cran.r-project.org/web/packages/ktspair/index.html
library(ktspair)
exprs<-as.matrix(read.table("data_raw.txt",header=T,sep='\t',row.names=1, as.is=T))
exprs<-response.data$x
colnames(exprs)<-response.data$samplelabels
rownames(exprs)<-response.data$geneid
grp<-response.data$y
ktsp<-ktspcalc(exprs,grp,3)
ktspplot(ktsp)
# plot.ktsp(ktsp)
# print.ktsp(ktsp)
write.table(ktsp$ktspdat,"clipboard",sep='\t')
write.table(cbind(rowMeans(ktsp$ktspdat[,1:22]),rowMeans(ktsp$ktspdat[,23:38])),"clipboard",sep='\t')
IDs<-ktsp$index
for (ii in 1:length(ktsp$index)){
IDs[[ii]]<-  rownames(exprs)[ktsp$index[[ii]]]
}
write.table(IDs,"clipboard",sep='\t')









names(ktsp)
# [1] "index"       "ktspscore"   "grp"         "ktspdat"     "k"          
# [6] "labels"      "rankscore"   "sensitivity" "specificity" "med"
heatmap.2(as.matrix(ktsp$ktspdat),Rowv=T, Colv=F,trace="none",col=color,distfun=function(x){dist(x,method=dist.method)}, hclustfun=function(x){hclust(x,method=hclust.method)}, density.info="none",cexCol=1.5,cexRow=1)#, cellnote=formatC(1/10^abs(mtx), format="e", digits=2), notecol='darkgreen')
