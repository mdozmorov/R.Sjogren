Sjogren syndrome microarray data analysis
========================================================
Topics to be covered:
* Data preparation
* Using ComBat to account for batch effect
* TODO





```
## ==========================================================================
## *
## *  Package WGCNA 1.34 loaded.
## *
## *    Important note: It appears that your system supports multi-threading,
## *    but it is not enabled within WGCNA in R. 
## *    To allow multi-threading within WGCNA with all available cores, use 
## *
## *          allowWGCNAThreads()
## *
## *    within R. Use disableWGCNAThreads() to disable threading if necessary.
## *    Alternatively, set the following environment variable on your system:
## *
## *          ALLOW_WGCNA_THREADS=<number_of_processors>
## *
## *    for example 
## *
## *          ALLOW_WGCNA_THREADS=12
## *
## *    To set the environment variable in linux bash shell, type 
## *
## *           export ALLOW_WGCNA_THREADS=12
## *
## *     before running R. Other operating systems or shells will
## *     have a similar command to achieve the same aim.
## *
## ==========================================================================
```

















Data preparation
--------------------

`annot.txt` is taken from `Single-experiment raw data3//annot3.txt`

`data.txt` are combined from `01.Data.Analysis.xlsx` and `11.Data.xlsx`.

`meta1.txt` is taken from `11 Annotations final cohorts sheet 2 20DEC13 DF.xlsx` and contains cohort information. The data has been transposed  to have patients names as columns for compatibility with `data.txt`

`meta2.txt` is taken from `20140214FarrisMicroarraymgdb(1).xlsx` (accessed on 03-13-2014). The data has been transposed. ID `p1033680-6` has been renamed to `p1033680-6 (-5)` to be compatible with the `data.txt` header.




To identify the largest source of variability within the data, we perform principal component analysis and color the samples by Cohort/Treatment status.


```r
outliersRemove <- FALSE
meta <- loadMeta(outliersRemove)
exprs.n <- loadExprs(outliersRemove)
# arrayQualityMetrics(new('ExpressionSet', exprs=exprs.n), outdir='arrayQC_WithOutliers')
```



```r
prinComponents(log2(exprs.n))
```

<img src="img/PCAwithOutliers.png" title="plot of chunk PCAwithOutliers" alt="plot of chunk PCAwithOutliers" width="700" />


Clearly, the cohorts are very different. The two patients, p1033216.2 and p1033680.6...5., appear as outliers. They are also picked up by the _arrayQualityMetrics_ set of tests. 

We remove them and look at the principal components again.


```r
outliersRemove <- TRUE
meta <- loadMeta(outliersRemove)
exprs.n <- loadExprs(outliersRemove)
# arrayQualityMetrics(new('ExpressionSet', exprs=exprs.n), outdir='arrayQC_WithoutOutliers')
```



```r
prinComponents(log2(exprs.n))
```

<img src="img/PCAwithoutOutliers.png" title="plot of chunk PCAwithoutOutliers" alt="plot of chunk PCAwithoutOutliers" width="700" />


The data looks more homogeneous now, with the cohort effect still dominating.

To further investigate how the outliers affect differential gene expression, we will perform differential expression _limma_ analysis on the expression data with and without outlier, as well as check for the enriched pathways. 


```r
outliersRemove <- FALSE
meta <- loadMeta(outliersRemove)
exprs.n <- loadExprs(outliersRemove)
# arrayQualityMetrics(new('ExpressionSet', exprs=exprs.n), outdir='arrayQC_WithOutliers')
res <- limmaOnData(exprs.n, "MicroarrayClass", "limma_MicroarrayClass_WithOutliers_WithBatch.txt")
nrow(res)
```

```
## [1] 41
```

```r
res.pathway <- reactomeEnrichment(res)
res.pathway[[2]]
```

```
## [1] 0
```


Not many probes are differentially expressed out of 62976 total. No pathway enrichment either.


```r
outliersRemove <- TRUE
meta <- loadMeta(outliersRemove)
exprs.n <- loadExprs(outliersRemove)
# arrayQualityMetrics(new('ExpressionSet', exprs=exprs.n), outdir='arrayQC_WithoutOutliers')
res <- limmaOnData(exprs.n, "MicroarrayClass", "limma_MicroarrayClass_WithoutOutliers_WithBatch.txt")
nrow(res)
```

```
## [1] 169
```

```r
res.pathway <- reactomeEnrichment(res)
res.pathway[[2]]
```

```
## [1] 0
```


Removing the outliers definitely helps to improve sensitivity in detecting differentially expressed genes. But not enough to see any pathway enrichment.

Using ComBat to account for batch effect
-----------------------------------------

To investigate batch effect in raw data, we will test the differences between the cohorts of patients, without considering their treatment status. The outliers were removed.


```r
res <- limmaOnData(exprs.n, "Cohort", "limma_Cohorts_WithoutOutliers_WithBatch.txt")
nrow(res)
```

```
## [1] 19562
```

```r
res.pathway <- reactomeEnrichment(res)
res.pathway[[2]]
```

```
## [1] 0
```


Almost 1/3 of the probes are differentially expressed between cohorts. They show no pathway enrichment.

We use ComBat to manually adjust for the cohort effect and keep the treatment effect. If we look at the PCA plot after the adjustment, the data look more homogeneous, and the cohorts are now mixed together.


```r
combat_edata <- ComBat(dat = exprs.n, batch = meta$Cohort, mod = model.matrix(~as.factor(MicroarrayClass), data = meta), 
    numCovs = NULL, par.prior = TRUE, prior.plots = F)
```

```
## Found 2 batches
## Found 1  categorical covariate(s)
## Standardizing Data across genes
## Fitting L/S model and finding priors
## Finding parametric adjustments
## Adjusting the Data
```

```r
prinComponents(combat_edata)
```

<img src="img/combat.png" title="plot of chunk combat" alt="plot of chunk combat" width="700" />


And no differentially expressed genes can be detected after removing the batch effect.


```r
res <- limmaOnData(combat_edata, "Cohort", "limma_Cohort_WithoutOutliers_WithoutBatch.txt")
nrow(res)
```

```
## [1] 0
```


Enrichment analysis
===================
Reactome canonical pathway enrichment analysis
----------------------------------------------

So far we used the complete dataset with 62976 probes. However, probes constantly expressed across the conditions are less interesting to us. To increase power of detecting differentially expressed genes, we filter such constantly expressed genes.


```r
# Filtering low variability probes
combat_edata <- exprs(varFilter(new("ExpressionSet", exprs = combat_edata)))
nrow(combat_edata)
```

```
## [1] 31488
```


Almost half of the probes were removed. The number of DEGs expectedly increased (with less probes we have more power to detect DEGs) from 754 to 784.

After removing batch effect, we not only get more differentially expressed genes, but also enriched pathways. The total number, and the top 30 most tignificant pathways are shown.


```r
res <- limmaOnData(combat_edata, "MicroarrayClass", "limma_MicroarrayClass_WithoutOutliers_WithoutBatch.txt")
nrow(res)
```

```
## [1] 784
```

```r
res.pathway <- reactomeEnrichment(res)
res.pathway[[2]]
```

```
## [1] 6
```

```r
grid.table(res.pathway[[1]][1:30, ], gp = gpar(fontsize = 7))
```

<img src="img/limmaOnTreatmentCombat.png" title="plot of chunk limmaOnTreatmentCombat" alt="plot of chunk limmaOnTreatmentCombat" width="700" />


Let's have a look at top 50 differentially expressed genes. Multiple probes for the same gene are collapsed to one by maximum expression level.


```r
degs <- collapseRows(combat_edata[rownames(res), ], annot.f[rownames(res), "GeneName"], rownames(res))
# Sort the collapsed data by the largest absolute fold change
maxRatio <- rowMeans(degs$datETcollapsed[, meta$MicroarrayClass == "PSS"]) - rowMeans(degs$datETcollapsed[, meta$MicroarrayClass == 
    "AC"])
degs.sorted <- degs$datETcollapsed[order(abs(maxRatio), decreasing = T), order(meta$MicroarrayClass)]
# Or just use the most significant data, columns reordered degs.sorted <- degs$datETcollapsed[,
# order(meta$MicroarrayClass)]
color <- colorRampPalette(c("blue", "white", "red"))
if (HIDE) labCol <- "" else labCol <- colnames(degs.sorted)
if (HIDE) labRow <- "" else labRow <- rownames(degs.sorted)
heatmap.2(degs.sorted[1:50, ], Colv = F, Rowv = F, scale = "row", trace = "none", col = color, key = F, density.info = "none", 
    cexCol = 1, cexRow = 0.8, labCol = labCol, labRow = labRow[1:50])  #, cellnote=formatC(1/10^abs(mtx), format='e', digits=2), notecol='darkgreen')
```

<img src="img/limmaVisual.png" title="plot of chunk limmaVisual" alt="plot of chunk limmaVisual" width="700" />


We aggregate probe names summarizing expression of the same gene by maximum fold change (because Ingenuity does it wrong).


```r
# Merge limma results and annotations
tmp <- merge(res, annot.f, by = "row.names")
# Get unique 'gene name - max logFC' pairs'
degs <- aggregate(tmp[, "logFC"], by = list(tmp$GeneName), max)
write.table(degs, "results/limma_MicroarrayClass_WithoutOutliers_WithoutBatch_IPA.txt", sep = "\t", quote = F, col.names = NA)
```



Gene ontology enrichment analysis
-------------------------------------

First, look at the total number and the top 20 significant ontologies in the "Biological Process" space (Normally, most informative functions).


```r
res.go <- GOEnrichment(res, "BP")
res.go[[2]]
```

```
## [1] 640
```

```r
grid.table(res.go[[1]][1:20, ], gp = gpar(fontsize = 6))
```

<img src="img/GOBPEnrichment.png" title="plot of chunk GOBPEnrichment" alt="plot of chunk GOBPEnrichment" width="700" />


Second, look at the total number and the top 20 significant ontologies in the "Molecular Function" space (Second most informative).


```r
res.go <- GOEnrichment(res, "MF")
res.go[[2]]
```

```
## [1] 123
```

```r
grid.table(res.go[[1]][1:20, ], gp = gpar(fontsize = 6))
```

<img src="img/GOMFEnrichment.png" title="plot of chunk GOMFEnrichment" alt="plot of chunk GOMFEnrichment" width="700" />


Finally, look at the total number and the top 20 significant ontologies in the "Cellular Component" space (Least informative).


```r
res.go <- GOEnrichment(res, "CC")
res.go[[2]]
```

```
## [1] 87
```

```r
grid.table(res.go[[1]][1:20, ], gp = gpar(fontsize = 6))
```

<img src="img/GOCCEnrichment.png" title="plot of chunk GOCCEnrichment" alt="plot of chunk GOCCEnrichment" width="700" />



TODO:

* Differential expression on other categorical clinical parameters

* Simple correlation of genes with  numerical clinical parameters

* Machine learning on clinical meta-data
