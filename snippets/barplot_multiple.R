# Plotting multiple histograms from a single matrix
# Guthridge BOLD data

# Function plotting error bars. http://monkeysuncle.stanford.edu/?p=485
error.bar <- function(x, y, upper, lower=upper, length=0.1,...){
  if(length(x) != length(y) | length(y) !=length(lower) | length(lower) != length(upper))
    stop("vectors must be same length")
  arrows(x,y+upper, x, y-lower, angle=90, code=3, length=length, ...)
}


gene<-readLines("clipboard")
data<-exprs.q[annot[,2] %in% gene,] # Full dataset
rowGroup<-annot[annot[,2] %in% gene,2]
rowID<-seq(1:length(rowGroup))
data<-collapseRows(data,rowGroup,rowID,method="MaxMean")$datETcollapsed

data<-read.table("clipboard",header=T,sep='\t',row.names=1) # Read in the matrix
# colnames(data)<-c("V1","VF") # ,"VF")
# rownames(data)<-data$SYMBOL # Assign row names to the column
# data$SYMBOL<-NULL # Remove that column

# Heme oxygenase
# PROBE_ID SYMBOL                                              DESCRIPTION
# 13050 ILMN_1800512  HMOX1 Homo sapiens heme oxygenase (decycling) 1 (HMOX1), mRNA.
# 13051 ILMN_1658807  HMOX2 Homo sapiens heme oxygenase (decycling) 2 (HMOX2), mRNA.
# data<-geneExpMy.filt.avg[annot[grep("HMOX",annot[,3]),1],]
# colnames(data)<-c("V1","V6","VF")
# rownames(data)<-annot[grep("HMOX",annot[,3]),2]

# Data as precomputed averages, no error bars possible
for (i in 1:nrow(data)) {
  data.min<-min(as.numeric(data[i,])) - 0.1*min(as.numeric(data[i,])) # Get the min -10%
  data.max<-max(as.numeric(data[i,])) + 0.1*max(as.numeric(data[i,])) # Get the max +10%
  barx<-barplot(as.numeric(data[i,]), main=row.names(data)[i], names.arg=colnames(data), 
          ylim=c(data.min,data.max), xpd=F) # Plot a hits for current data
}

# Data as full matrix, can calculate error bars
dev.off() # Clear the plot area
pdf("test.pdf")
par(mfrow=c(3,2)) # Num X Num plot area
indices <- lapply(c("NC","AC","SS"), function(p) {grep(p, colnames(data))}) # Group specific column indices
# Get all averages into one matrix
data.means <- matrix(unlist(lapply(indices, function(idxs) {
  rowMeans(data[,idxs])
})), nrow=nrow(data))
# Get all SDs into one matrix
data.sd <- matrix(unlist(lapply(indices, function(idxs) {
  apply(data[,idxs],1,sd)
})), nrow=nrow(data))
# Loop through it and plot
for (i in 1:dim(data)[[1]]) {
    data.min<-min(as.numeric(data[i,])) - 0.1*min(as.numeric(data[i,])) # Get the min -10%
    data.max<-max(as.numeric(data[i,])) + 0.1*max(as.numeric(data[i,])) # Get the max +10%
    barx <- barplot(data.means[i,], names.arg=c("NC","AC","SS") ,xlab=rownames(data)[i], ylim=c(data.min,data.max),xpd=F, ylab="Log2 expression")
    error.bar(barx, data.means[i,], 1.96*data.sd[i,]/sqrt(dim(data)[[2]]/3))
}
dev.off()

# More blunt but working method, processing one row at a time
dev.off() # Clear the plot area
par(mfrow=c(4,4)) # Num X Num plot area
names<-c("V1","V6","VF")
indices <- lapply(names, function(p) {grep(p, colnames(data))}) # Group specific column indices
for (i in 18:43) { # nrow(data)
  data.means<-c(rowMeans(data[i,indices[[1]]]),rowMeans(data[i,indices[[2]]]),rowMeans(data[i,indices[[3]]]))
  data.sd<-1.96*c(apply(data[i,indices[[1]]],1,sd)/sqrt(length(data[i,indices[[1]]])),apply(data[i,indices[[2]]],1,sd)/sqrt(length(data[i,indices[[2]]])),apply(data[i,indices[[3]]],1,sd)/sqrt(length(data[i,indices[[3]]])))
  data.min<-min(as.numeric(data[i,])) - 0.1*min(as.numeric(data[i,])) # Get the min -10%
  data.max<-max(as.numeric(data[i,])) + 0.1*max(as.numeric(data[i,])) # Get the max +10%
  barx<-barplot(data.means, main=row.names(data)[i], names.arg=names, 
           ylim=c(data.min,data.max), xpd=F) # Plot a hits for current data
  error.bar(barx,data.means,data.sd)
}

# Data as raw numbers
y1<-as.matrix(t(read.table("clipboard",sep='\t'))) # Separate visits
y2<-as.matrix(t(read.table("clipboard",sep='\t')))
y2<-na.omit(y2) # will remove two columns, as they have NAs
y3<-as.matrix(t(read.table("clipboard",sep='\t')))
y1.means <- apply(y1,2,mean)
y1.sd <- apply(y1,2,sd)
y2.means <- apply(y2,2,mean)
y2.sd <- apply(y2,2,sd)
y3.means <- apply(y3,2,mean)
y3.sd <- apply(y3,2,sd)


yy<-matrix(c(y1.means,y2.means,y3.means),3,32,byrow=T)
yy.names<-read.table("clipboard",sep='\t')
yyt.names<-t(yy.names)
ee<-matrix(c(y1.sd,y2.sd,y3.sd),3,32,byrow=T)*1.96/sqrt(32)

dev.off()
par(mfrow=c(4,4))
for (i in 1:ncol(yy)) {
  data.min<-min(yy[,i]) - 0.5*min(yy[,i]) # Get the min -10%
  data.max<-max(yy[,i]) + 0.5*max(yy[,i]) # Get the max +10%
  barx<-barplot(yy[,i], beside=TRUE,col=c("blue","magenta","red"),ylim=c(data.min,data.max),names.arg=yyt.names[i],xpd=F,axis.lty=1)
  error.bar(barx,yy[,i],ee[,i])
}

# Stacked barplot


barplot(as.table(as.matrix(data)),legend.text=T,col=rainbow(7))