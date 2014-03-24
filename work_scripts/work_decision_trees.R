
cols.nums <- c("Cohort", "MicroarrayClass", "Ageatbiopsy", "LGleft", "LGright", "OSSL", "OSSR", "SchirmerL", "SchirmerR", "WUSFvol", "SSFvolL", "SSFvolR", "FS", "LaBioRadvalue", "RoBioRadvalue", "RFvalue")

cols.text <- c("RaceNIH", "Ethnicity", "Sex", "AECGClass", "ACRClass", "LGcall", "LissamineGreenSumACR", "Schirmercall", "WUSFcall", "Biopsycall", "Lacall", "Rocall", "RFcall")

cols.excl <- c("DateClinic", "ANAcall", "HighestANA", "ReichlinANAvalue", "Secondaryforms", "Medications", "AECGClass", "ACRClass", "Biopsycall", "FS")

meta.ds <- read.table("data//meta.txt", sep="\t", row.names=1, header=T, as.is=T)[, cols.nums]#[, !(colnames(meta) %in% cols.excl)]

meta.target <- "MicroarrayClass"

(meta.numerics <- names(meta.ds)[sapply(meta.ds, mode) == "numeric"])
(meta.categorics <- names(meta.ds)[sapply(meta.ds, mode) == "character"])
(meta.form <- formula(paste(meta.target, "~ .")))

meta.model <- rpart(formula=meta.form, data=meta.ds)
fancyRpartPlot(meta.model)

summary(meta.model)

plot(meta.model, uniform=T)
text(meta.model, use.n=T, all=T, cex=0.8)


meta.model$cptable
meta.model$variable.importance

asRules.rpart <- function(model)
{
  if (!inherits(model, "rpart")) stop("Not a legitimate rpart tree")
  #
  # Get some information.
  #
  frm     <- model$frame
  names   <- row.names(frm)
  ylevels <- attr(model, "ylevels")
  ds.size <- model$frame[1,]$n
  #
  # Print each leaf node as a rule.
  #
  for (i in 1:nrow(frm))
  {
    if (frm[i,1] == "<leaf>")
    {
      # The following [,5] is hardwired - needs work!
      cat("\n")
      cat(sprintf(" Rule number: %s ", names[i]))
      cat(sprintf("[yval=%s cover=%d (%.0f%%) prob=%0.2f]\n",
                  ylevels[frm[i,]$yval], frm[i,]$n,
                  round(100*frm[i,]$n/ds.size), frm[i,]$yval2[,5]))
      pth <- path.rpart(model, nodes=as.numeric(names[i]), print.it=FALSE)
      cat(sprintf("   %s\n", unlist(pth)[-1]), sep="")
    }
  }
}
asRules(meta.model)

library(genefilter)

hist(apply(combat_edata, 1, sd), n=100)

ds <- data.frame(c(t(combat_edata), meta$MicroarrayClass))
colnames(ds)[length(colnames(ds))] <- "MicroarrayClass"
target <- "MicroarrayClass"

(form <- formula(paste(meta.target, "~ .")))

model <- rpart(formula=form, data=ds)
fancyRpartPlot(meta.model)
