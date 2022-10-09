Cooking miRNA-Seq data
================
Pilar González Marchante

- <a href="#loading" id="toc-loading">Loading</a>
- <a href="#filtering" id="toc-filtering">Filtering</a>
- <a href="#exploring" id="toc-exploring">Exploring</a>
- <a href="#normalizing" id="toc-normalizing">Normalizing</a>
  - <a href="#cqn" id="toc-cqn">cqn</a>
  - <a href="#edaseq" id="toc-edaseq">EDASeq</a>
- <a href="#analyzing-differential-expression"
  id="toc-analyzing-differential-expression">Analyzing differential
  expression</a>
  - <a href="#deseq2" id="toc-deseq2">DESeq2</a>
  - <a href="#limma" id="toc-limma">limma</a>
  - <a href="#limma-voom" id="toc-limma-voom">limma-voom</a>
  - <a href="#edger" id="toc-edger">edgeR</a>
  - <a href="#intersecting-degs" id="toc-intersecting-degs">Intersecting
    DEGs</a>

# Loading

We have information for 1,881 miRNAs, in rows, and 2,593 columns. For
each sample, there are 3 columns: `read_count` (raw counts),
`reads_per_million_miRNA_mapped` (counts normalized by dividing the read
counts of a miRNA by the total read counts of the sample) and
`cross-mapped` (which can be Y for YES or N for NO; if Y, it indicates
that a single read aligns to more than one miRNA). Let’s just use the
raw counts (`read_count`) for now.

``` r
load("data/raw/miRNA-Seq/miRNA.rda")

head(mirna)[, 1:4]
```

    ##       miRNA_ID read_count_TCGA-E2-A1L7-01A-11R-A143-13
    ## 1 hsa-let-7a-1                                   15576
    ## 2 hsa-let-7a-2                                   15513
    ## 3 hsa-let-7a-3                                   15781
    ## 4   hsa-let-7b                                   13032
    ## 5   hsa-let-7c                                    1225
    ## 6   hsa-let-7d                                     662
    ##   reads_per_million_miRNA_mapped_TCGA-E2-A1L7-01A-11R-A143-13
    ## 1                                                  10986.0037
    ## 2                                                  10941.5688
    ## 3                                                  11130.5935
    ## 4                                                   9191.6795
    ## 5                                                    864.0122
    ## 6                                                    466.9193
    ##   cross-mapped_TCGA-E2-A1L7-01A-11R-A143-13
    ## 1                                         N
    ## 2                                         N
    ## 3                                         N
    ## 4                                         N
    ## 5                                         N
    ## 6                                         N

``` r
dim(mirna)
```

    ## [1] 1881 2593

``` r
rownames(mirna) <- mirna$miRNA_ID
mirna <- mirna[, -1]

mirna.raw.counts <- colnames(mirna)[grep("count", colnames(mirna))]
mirna.raw.counts <- mirna[,mirna.raw.counts]
colnames(mirna.raw.counts) <- gsub("read_count_","", colnames(mirna.raw.counts))
rm(mirna)
```

Let’s look at the data.

``` r
boxplot(mirna.raw.counts[, 1:50] + 1, log = "y", outline = FALSE, las = 2)
```

![](images/cookingmiRNASeq/boxplot.raw.png)

# Filtering

We already have our miRNA data in `mirna.raw.counts`, but we still have
to define our factors as `condition`, `tss`, `plate`, `portion`and
`sample`, which we can easily extract with the `TCGAbiolinks` function
`get_IDs`.

``` r
library(NOISeq)
library(TCGAbiolinks)
```

``` r
barcodes <- get_IDs(mirna.raw.counts)
myfactors <- data.frame(barcodes$tss, barcodes$portion, barcodes$plate, barcodes$condition)
head(myfactors)
```

    ##   barcodes.tss barcodes.portion barcodes.plate barcodes.condition
    ## 1           E2              11R           A143             cancer
    ## 2           E2              33R           A143             normal
    ## 3           BH              11R           A22I             cancer
    ## 4           D8              11R           A14L             cancer
    ## 5           AC              11R           A36A             cancer
    ## 6           A2              11R           A085             cancer

Choosing a CPM threshold.

``` r
library(SummarizedExperiment)
library(edgeR)
library(limma)
library(ggplot2)

mean_log_cpm <- aveLogCPM(mirna.raw.counts)

filter_threshold <- log2(0.5)

ggplot() + aes(x=mean_log_cpm) +
    geom_histogram(binwidth=0.2) +
    geom_vline(xintercept=filter_threshold) +
    ggtitle("Histogram of logCPM before filtering")
```

![](images/cookingmiRNASeq/histogram.logCPM.png)

``` r
ggplot() + aes(x=mean_log_cpm) +
    geom_density() +
    geom_vline(xintercept=filter_threshold) +
    ggtitle("Density plot of logCPM before filtering") +
    xlim(-6.1, 13.5)
```

![](images/cookingmiRNASeq/density.logCPM.png)

So let’s try CPM filtering with a `CPM threshold = 0.2, 0.5 and 1` and a
`cv.cutoff = 500`, so that we remove those features with low expression
(but not with low variability). We will also apply Wilcoxon test
filtering and compare the results.

``` r
myfiltCPM01 <- filtered.data(mirna.raw.counts, factor = myfactors$barcodes.condition, norm = FALSE, depth = NULL, method = 1, cv.cutoff = 500, cpm = 0.1, p.adj = "fdr") # 799 features (35.2%) are to be kept for differential expression analysis with filtering method 1

myfiltCPM02 <- filtered.data(mirna.raw.counts, factor = myfactors$barcodes.condition, norm = FALSE, depth = NULL, method = 1, cv.cutoff = 500, cpm = 0.2, p.adj = "fdr") # 662 features (35.2%) are to be kept for differential expression analysis with filtering method 1

myfiltCPM05 <- filtered.data(mirna.raw.counts, factor = myfactors$barcodes.condition, norm = FALSE, depth = NULL, method = 1, cv.cutoff = 500, cpm = 0.5, p.adj = "fdr") # 501 features (26.6%) are to be kept for differential expression analysis with filtering method 1

myfiltCPM1 <- filtered.data(mirna.raw.counts, factor = myfactors$barcodes.condition, norm = FALSE, depth = NULL, method = 1, cv.cutoff = 500, cpm = 1, p.adj = "fdr") # 420 features (22.3%) are to be kept for differential expression analysis with filtering method 1

myfiltWilcoxon <- filtered.data(mirna.raw.counts, factor = myfactors$barcodes.condition, norm = FALSE, depth = NULL, method = 2, p.adj = "fdr") # 1434 (76.2%) features are to be kept for differential expression analysis with filtering method 2
```

``` r
boxplot(log10(myfiltCPM02[, 1:50])+1, outline = FALSE, las = 2)
```

![](images/cookingmiRNASeq/boxplot.filt.CPM.02.png)

``` r
boxplot(log10(myfiltCPM05[, 1:50])+1, outline = FALSE, las = 2)
```

![](images/cookingmiRNASeq/boxplot.filt.CPM.05.png)

``` r
boxplot(log10(myfiltCPM1[, 1:50])+1, outline = FALSE, las = 2)
```

![](images/cookingmiRNASeq/boxplot.filt.CPM.1.png)

``` r
boxplot(log10(myfiltWilcoxon[, 1:50])+1, outline = FALSE, las = 2)
```

![](images/cookingmiRNASeq/boxplot.filt.Wilcoxon.png)

We’ll choose method 1 with a CPM threshold of 0.5.

# Exploring

``` r
mirna.filt.counts <- myfiltCPM05

rm(myfiltCPM05)

save(mirna.filt.counts, file = "data/cooked/miRNA-Seq/miRNA.filt.rda")
```

And we’ll need GC content information.

``` r
# save ids to input into biomart
# write.csv(rownames(mirna.filt.counts), "results/preprocessing/cookingmiRNASeq/miRNA.IDs.txt", row.names = FALSE, col.names = FALSE, quote = FALSE)

GCcontent <- read.csv("results/preprocessing/cookingmiRNASeq/miRNA.biomart.txt", sep = "\t")

mygc = c(GCcontent$Gene...GC.content)
names(mygc) = c(GCcontent$miRBase.ID)
mygc = mygc[rownames(mirna.filt.counts)]
names(mygc) = rownames(mirna.filt.counts)

mygc <- na.omit(mygc) # 10 mirnas without information

mirna.filt.counts <- mirna.filt.counts[names(mygc),]

save(mygc, file = "results/preprocessing/cookingmiRNASeq/GC.miRNA.rda")

mymirnadata.filt <- NOISeq::readData(data = mirna.filt.counts, factors = myfactors, gc = mygc)
```

Since the number of miRNAs (491) was too low to be able to plot GC bias
(given that by default, `NOISeq` makes each bin have 200 features) I had
to change the source code of the functions so that each bin had 40
features instead.

``` r
source(file = "scripts/preprocessing/GCbias.R")
GC.plot(GC.dat(mymirnadata.filt, factor = "barcodes.condition"))
```

![](images/cookingmiRNASeq/gc.bias.filt.png)

There seems to be some GC content bias, so we’ll trying normalizing
those out with `cqn` and `EDASeq`.

``` r
library(NOISeq)
mymirnaPCA = dat(mymirnadata.filt, type = "PCA")
par(cex = 0.75)
explo.plot(mymirnaPCA, factor = "barcodes.condition", plottype = "scores")
```

![](images/cookingmiRNASeq/pca.scores.filt.png)

``` r
explo.plot(mymirnaPCA, factor = "barcodes.condition", plottype = "loadings")
```

![](images/cookingmiRNASeq/pca.loadings.filt.png)

``` r
explo.plot(mymirnaPCA, factor = "barcodes.tss")
```

![](images/cookingmiRNASeq/pca.tss.filt.png)

``` r
explo.plot(mymirnaPCA, factor = "barcodes.portion")
```

![](images/cookingmiRNASeq/pca.portion.filt.png)

``` r
explo.plot(mymirnaPCA, factor = "barcodes.plate")
```

![](images/cookingmiRNASeq/pca.plate.filt.png)

# Normalizing

## cqn

`cqn` with all length = 1000.

``` r
library(cqn)
load("results/preprocessing/cookingmiRNASeq/GC.miRNA.rda")
sizeFactors.mirna <- colSums(mirna.filt.counts)

mirna.cqn.norm <- cqn(mirna.filt.counts, lengthMethod = "fixed", lengths = rep(100, 491), x = mygc, sizeFactors = sizeFactors.mirna, verbose = TRUE)

save(mirna.cqn.norm, file = "reports/preprocessing/files/cookingmiRNASeq/miRNA.cqn.norm.rda")

# Extract normalized data to check for bias on NOISeq
mirna.cqn.norm.expression <- mirna.cqn.norm$y + mirna.cqn.norm$offset
mirna.cqn.norm.expression <- as.data.frame(mirna.cqn.norm.expression)
```

## EDASeq

Now `EDASeq`.

``` r
library(EDASeq)

feature <- data.frame(gc=mygc)

data <- newSeqExpressionSet(counts=as.matrix(mirna.filt.counts), featureData=feature, phenoData=data.frame(conditions=barcodes$condition, row.names=barcodes$barcode))

dataWithin <- withinLaneNormalization(data, "gc", which="full")
mirna.eda.norm <- betweenLaneNormalization(dataWithin, which="full")

save(mirna.eda.norm, file = "reports/preprocessing/files/cookingmiRNASeq/miRNA.eda.norm.rda")

load("reports/preprocessing/files/cookingmiRNASeq/miRNA.eda.norm.rda")
```

“Normalization factors should be on the scale of the counts, like size
factors, and unlike offsets which are typically on the scale of the
predictors (i.e. the logarithmic scale for the negative binomial GLM).
At the time of writing, the transformation from the matrices provided by
these packages should be”:

``` r
# Extract the offset, which will be input directly into DEseq2 to normalise the counts. 
normFactors <- withinLaneNormalization(data,"gc",
                                       which="full", offset=TRUE)
normFactors <- betweenLaneNormalization(normFactors,
                                        which="full", offset=TRUE)

normFactors <- exp(-1 * normFactors@assayData$offset)
save(normFactors, file = "reports/preprocessing/files/cookingmiRNASeq/miRNA.normFactors.rda")
```

``` r
# Extract normalized counts to check for bias on NOISeq
mirna.eda.counts <- mirna.eda.norm@assayData$normalizedCounts
mirna.eda.counts <- as.data.frame(mirna.eda.counts)
```

Check if it fixed bias.

``` r
mymirnadata.norm.cqn <- NOISeq::readData(data = mirna.cqn.norm.expression, factors = myfactors, gc = mygc)
GC.plot(GC.dat(mymirnadata.norm.cqn, factor = "barcodes.condition"))
```

![](images/cookingmiRNASeq/gc.bias.cqn.norm.png)

``` r
mymirnadata.norm.eda <- NOISeq::readData(data = mirna.eda.counts, factors = myfactors, gc = mygc)
GC.plot(GC.dat(mymirnadata.norm.eda, factor = "barcodes.condition"))
```

![](images/cookingmiRNASeq/gc.bias.eda.norm.png)

`cqn` normalization seems to increase mean expression by a lot, so we’ll
choose `EDASeq` as it also does a good job of reducing GC content bias.
In any case, we’ll check the separation of the samples in PCA plots to
help us make up our mind.

``` r
mymirnaPCA = dat(mymirnadata.norm.cqn, type = "PCA", logtransf = TRUE)
par(cex = 0.75)
explo.plot(mymirnaPCA, factor = "barcodes.condition", plottype = "scores")
```

![](images/cookingmiRNASeq/pca.scores.cqn.norm.png)

``` r
explo.plot(mymirnaPCA, factor = "barcodes.condition", plottype = "loadings")
```

![](images/cookingmiRNASeq/pca.loadings.cqn.norm.png)

``` r
mymirnaPCA = dat(mymirnadata.norm.eda, type = "PCA")
par(cex = 0.75)
explo.plot(mymirnaPCA, factor = "barcodes.condition", plottype = "scores")
```

![](images/cookingmiRNASeq/pca.scores.eda.norm.png)

``` r
explo.plot(mymirnaPCA, factor = "barcodes.condition", plottype = "loadings")
```

![](images/cookingmiRNASeq/pca.loadings.eda.norm.png)

`EDASeq` does a much better job at separating the sample groups by
scores, as they seem to be more spread out horizontally (even though the
variance explained by PC1 is the same in both cases: 12%). We’ll use
`EDASeq` normalization from now on.

# Analyzing differential expression

## DESeq2

``` r
library(DESeq2)
library(TCGAbiolinks)

load("data/cooked/miRNA-Seq/miRNA.filt.rda")

library(TCGAbiolinks)
barcodes <- get_IDs(mirna.filt.counts)
barcodes$condition <- as.factor(barcodes$condition)
barcodes$condition <- relevel(barcodes$condition, ref = "normal")
myfactors <- data.frame(barcodes$tss, barcodes$portion, barcodes$plate, barcodes$condition)

dds <- DESeqDataSetFromMatrix(countData = mirna.filt.counts,
                              colData = barcodes,
                              design = ~ condition)

load("data/cooked/miRNA-Seq/miRNA.normFactors.rda")
# Before inputing normalizationFactors into DESeq2, you should divide out the geometric mean
normFactors <- normFactors / exp(rowMeans(log(normFactors)))
normalizationFactors(dds) <- normFactors

dds <- DESeq(dds)
resultsNames(dds) # lists the coefficients
# [1] "Intercept"                          
# [2] "condition_cancer_vs_normal"
res.deseq2 <- results(dds, alpha = 0.05)
summary(res.deseq2)
```

``` r
out of 491 with nonzero total read count
adjusted p-value < 0.05
LFC > 0 (up)       : 244, 50%
LFC < 0 (down)     : 140, 29%
outliers [1]       : 0, 0%
low counts [2]     : 0, 0%
(mean count < 0)
```

We had to decrease our lFC threshold from 2 to 1, as otherwise we ended
up getting a ridiculous number of DEGs (9 in total, after intersecting
all DEA methods).

We’ll select as significant those genes with a p.adj \< 0.05 and a lFC
\> 1 or lFC \< -1. We get a total of 133 upregulated DEGs (a 27% of the
filtered genes) in cancer samples, compared to normal ones; and 81
downregulated DEGs (a 16.5% of the filtered genes), from a total of 491
filtered miRNAs (originally 1,881 miRNAs in our raw data).

``` r
log.fold.change <- res.deseq2$log2FoldChange
q.value <- res.deseq2$padj
genes.ids <- rownames(mirna.filt.counts)
names(log.fold.change) <- genes.ids
names(q.value) <- genes.ids
activated.genes.deseq2 <- genes.ids[log.fold.change > 1 & q.value < 0.05]
activated.genes.deseq2 <- activated.genes.deseq2[!is.na(activated.genes.deseq2)]
repressed.genes.deseq2 <- genes.ids[log.fold.change < - 1 & q.value < 0.05]
repressed.genes.deseq2 <- repressed.genes.deseq2[!is.na(repressed.genes.deseq2)]
length(activated.genes.deseq2) # 133
length(repressed.genes.deseq2) # 81

log.q.val <- -log10(q.value)
plot(log.fold.change,log.q.val,pch=19,col="grey",cex=0.8,
xlim=c(-8,8),ylim = c(0,150),
xlab="log2(Fold-change)",ylab="-log10(q-value)",cex.lab=1.5)
points(x = log.fold.change[activated.genes.deseq2],
y = log.q.val[activated.genes.deseq2],col="red",cex=0.8,pch=19)
points(x = log.fold.change[repressed.genes.deseq2],
y = log.q.val[repressed.genes.deseq2],col="blue",cex=0.8,pch=19)
```

![](images/cookingmiRNASeq/volcano.plot.deseq2.png)

We don’t have to convert our DEGs IDs, since they are using standard
miRNA IDs. Let’s just save them.

``` r
write.table(activated.genes.deseq2, file = "results/preprocessing/cookingmiRNASeq/DESeq2.up.txt", row.names = FALSE, col.names = FALSE, quote = FALSE)

write.table(repressed.genes.deseq2, file = "results/preprocessing/cookingmiRNASeq/DESeq2.down.txt", row.names = FALSE, col.names = FALSE, quote = FALSE)

resOrdered <- res.deseq2[order(res.deseq2$padj),]
head(resOrdered)
resOrderedDF <- as.data.frame(resOrdered)
resOrderedDF <- na.omit(resOrderedDF)
write.table(resOrderedDF, file = "results/preprocessing/cookingmiRNASeq/DESeq2.ordered.csv", row.names=TRUE, col.names=TRUE, sep="\t", quote=FALSE)
```

## limma

``` r
library(limma)
load("data/cooked/miRNA-Seq/miRNA.norm.rda")

mirna.norm.expression <- log(mirna.norm@assayData$normalizedCounts)

design <- model.matrix(~ barcodes$condition)

fit1 <- lmFit(mirna.norm.expression, design)

fit2 <- eBayes(fit1)

top.limma <- topTable(fit2, coef = 2, number = Inf)

log.fold.change <- top.limma$logFC
q.value <- top.limma$adj.P.Val
genes.ids <- rownames(top.limma)
names(log.fold.change) <- genes.ids
names(q.value) <- genes.ids

activated.genes.limma <- genes.ids[log.fold.change > 1 & q.value < 0.05]
repressed.genes.limma <- genes.ids[log.fold.change < -1 & q.value < 0.05]

length(activated.genes.limma) # 42
length(repressed.genes.limma) # 49

log.q.val <- -log10(q.value)
plot(log.fold.change,log.q.val,pch=19,col="grey",cex=0.8,
xlim=c(-4,4),ylim = c(0,130),
xlab="log2(Fold-change)",ylab="-log10(q-value)",cex.lab=1.5)
points(x = log.fold.change[activated.genes.limma],
y = log.q.val[activated.genes.limma],col="red",cex=0.8,pch=19)
points(x = log.fold.change[repressed.genes.limma],
y = log.q.val[repressed.genes.limma],col="blue",cex=0.8,pch=19)
```

![](images/cookingmiRNASeq/volcano.plot.limma.png)

``` r
write.table(activated.genes.limma, file = "results/preprocessing/cookingmiRNASeq/limma.up.txt", row.names = FALSE, col.names = FALSE, quote = FALSE)

write.table(repressed.genes.limma, file = "results/preprocessing/cookingmiRNASeq/limma.down.txt", row.names = FALSE, col.names = FALSE, quote = FALSE)

topOrdered <- top.limma[order(top.limma$adj.P.Val),]
topOrderedDF <- as.data.frame(topOrdered)
topOrderedDF <- na.omit(topOrderedDF)
write.table(topOrderedDF, file = "results/preprocessing/cookingmiRNASeq/limma.ordered.csv", row.names=TRUE, col.names=TRUE, sep="\t", quote=FALSE)
```

## limma-voom

Even though it says expression, they are (normalized) counts, so we can
apply `limma-voom`.

``` r
mirna.norm.expression <- mirna.norm@assayData$normalizedCounts
v <- voom(mirna.norm.expression, design, plot=TRUE)
fit <- lmFit(v, design)
fit2 <- eBayes(fit)
top.limma.voom <- topTable(fit2, coef=ncol(design), number = Inf)

log.fold.change <- top.limma.voom$logFC
q.value <- top.limma.voom$adj.P.Val
genes.ids <- rownames(top.limma.voom)
names(log.fold.change) <- genes.ids
names(q.value) <- genes.ids

activated.genes.limma.voom <- genes.ids[log.fold.change > 1 & q.value < 0.05]
repressed.genes.limma.voom <- genes.ids[log.fold.change < -1 & q.value < 0.05]

length(activated.genes.limma.voom) # 77
length(repressed.genes.limma.voom) # 93

log.q.val <- -log10(q.value)
plot(log.fold.change,log.q.val,pch=19,col="grey",cex=0.8,
xlim=c(-6,6),ylim = c(0,150),
xlab="log2(Fold-change)",ylab="-log10(q-value)",cex.lab=1.5)
points(x = log.fold.change[activated.genes.limma.voom],
y = log.q.val[activated.genes.limma.voom],col="red",cex=0.8,pch=19)
points(x = log.fold.change[repressed.genes.limma.voom],
y = log.q.val[repressed.genes.limma.voom],col="blue",cex=0.8,pch=19)
```

![](images/cookingmiRNASeq/volcano.plot.limma.voom.png) We get quite a
different result with `limma-voom` compared to `limma`. In order to end
up with a sizeable number of DEGs, we’ll use `limma-voom` instead of
`limma`.

``` r
write.table(activated.genes.limma.voom, file = "results/preprocessing/cookingmiRNASeq/limma.voom.up.txt", row.names = FALSE, col.names = FALSE, quote = FALSE)

write.table(repressed.genes.limma.voom, file = "results/preprocessing/cookingmiRNASeq/limma.voom.down.txt", row.names = FALSE, col.names = FALSE, quote = FALSE)

topOrdered <- top.limma.voom[order(top.limma.voom$adj.P.Val),]
topOrderedDF <- as.data.frame(topOrdered)
topOrderedDF <- na.omit(topOrderedDF)
write.table(topOrderedDF, file = "results/preprocessing/cookingmiRNASeq/limma.voom.ordered.csv", row.names=TRUE, col.names=TRUE, sep="\t", quote=FALSE)
```

## edgeR

``` r
library(edgeR)

y <- DGEList(counts = mirna.filt.counts, lib.size = colSums(mirna.filt.counts), group = barcodes$condition, genes = rownames(mirna.filt.counts))

y$offset <- normFactors
```

``` r
design <- model.matrix(~ barcodes$condition)

y <- estimateCommonDisp(y, design = design)
y <- estimateTagwiseDisp(y, design = design)
plotBCV(y)
```

![](images/cookingmiRNASeq/plotBCV.edger.png)

This time it’s also better to estimate dispersions tagwise instead of
using a common one for all genes.

``` r
et <- exactTest(y) # performs pair-wise tests for differential expression between two groups
top.edger <- topTags(et, n = Inf) # takes the output from exactTest(), adjusts the raw p-values using the False Discovery Rate (FDR) correction, and returns the top differentially expressed genes

topSig <- top.edger[top.edger$table$FDR < 0.05, ] # we select DEGs with alpha=0.05
dim(topSig)
topSig <- topSig[abs(top.edger$table$logFC) >= 1, ] # we filter the output of dataDEGs by abs(LogFC) >=1
dim(topSig)
```

``` r
# is equivalent to doing

summary(decideTestsDGE(et, lfc = 1, p.value = 0.05))

       cancer-normal
Down              69
NotSig           290
Up               132
```

``` r
activated.genes.edger <- topSig$table$genes[topSig$table$logFC > 0] # 132
repressed.genes.edger <- topSig$table$genes[topSig$table$logFC < 0] # 69

log.q.val <- -log10(q.value)
plot(log.fold.change,log.q.val,pch=19,col="grey",cex=0.8,
xlim=c(-6,6),ylim = c(0,120),
xlab="log2(Fold-change)",ylab="-log10(q-value)",cex.lab=1.5)
points(x = log.fold.change[activated.genes.edger],
y = log.q.val[activated.genes.edger],col="red",cex=0.8,pch=19)
points(x = log.fold.change[repressed.genes.edger],
y = log.q.val[repressed.genes.edger],col="blue",cex=0.8,pch=19)
```

![](images/cookingmiRNASeq/volcano.plot.edger.png)

``` r
write.table(activated.genes.edger, file = "results/preprocessing/cookingmiRNASeq/edgeR.up.txt", row.names = FALSE, col.names = FALSE, quote = FALSE)

write.table(repressed.genes.edger, file = "results/preprocessing/cookingmiRNASeq/edgeR.down.txt", row.names = FALSE, col.names = FALSE, quote = FALSE)

top.edger <- as.data.frame(top.edger)

topOrdered <- top.edger[order(top.edger$FDR),]
topOrderedDF <- as.data.frame(topOrdered)
topOrderedDF <- na.omit(topOrderedDF)
write.table(topOrderedDF, file = "results/preprocessing/cookingmiRNASeq/edgeR.ordered.csv", row.names=TRUE, col.names=TRUE, sep="\t", quote=FALSE)
```

## Intersecting DEGs

``` r
# with DESeq2, limma and edgeR
common.activated <- intersect(intersect(activated.genes.deseq2, activated.genes.edger), activated.genes.limma) # 37
common.repressed <- intersect(intersect(repressed.genes.deseq2, repressed.genes.edger), repressed.genes.limma) # 36

# with DESeq2, limma-voom and edgeR
common.activated <- intersect(intersect(activated.genes.deseq2, activated.genes.edger), activated.genes.limma.voom) # 66
common.repressed <- intersect(intersect(repressed.genes.deseq2, repressed.genes.edger), repressed.genes.limma.voom) # 53

# with DESeq2, limma, limma-voom and edgeR
common.activated <- intersect(intersect(intersect(activated.genes.deseq2, activated.genes.edger), activated.genes.limma.voom), activated.genes.limma) # 36
common.repressed <- intersect(intersect(intersect(repressed.genes.deseq2, repressed.genes.edger), repressed.genes.limma.voom), repressed.genes.limma) # 36
```

We’ll select as DEGs as those intersecting for `DESeq2`, `limma-voom`
and `edgeR`, since `limma` appears to miss out on quite a few of them.

``` r
# with DESeq2, limma-voom and edgeR
common.activated <- intersect(intersect(activated.genes.deseq2, activated.genes.edger), activated.genes.limma.voom) # 66
common.repressed <- intersect(intersect(repressed.genes.deseq2, repressed.genes.edger), repressed.genes.limma.voom) # 53

write.table(common.activated, file = "results/preprocessing/cookingmiRNASeq/common.up.txt", row.names = FALSE, col.names = FALSE, quote = FALSE)

write.table(common.repressed, file = "results/preprocessing/cookingmiRNASeq/common.down.txt", row.names = FALSE, col.names = FALSE, quote = FALSE)
```

|    DEGs     | DESeq2 | limma-voom | edgeR | Common |
|:-----------:|:------:|:----------:|:-----:|:------:|
| *Activated* |  133   |     77     |  132  |   66   |
| *Repressed* |   81   |     93     |  69   |   53   |
|   *Total*   | *214*  |   *170*    | *201* | *119*  |

We have 119 DEGs, a 24.2% of the filtered miRNAs (491) and a 6.3% of the
original miRNAs (1,881).
