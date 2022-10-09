#------------------------- loading required packages --------------------------#
library(TCGAbiolinks)
library(SummarizedExperiment)
library(NOISeq)
library(EDASeq)
library(edgeR)
library(limma)
library(DESeq2)
library(ggplot2)

#--------------- filtering with NOISeq with CPM threshold = 0.5 ---------------#
load("data/raw/miRNA-Seq/miRNA.rda")

rownames(mirna) <- mirna$miRNA_ID
mirna <- mirna[, -1]

mirna.raw.counts <- colnames(mirna)[grep("count", colnames(mirna))]
mirna.raw.counts <- mirna[,mirna.raw.counts]
colnames(mirna.raw.counts) <- gsub("read_count_","", colnames(mirna.raw.counts))
rm(mirna)

barcodes <- get_IDs(mirna.raw.counts)
myfactors <- data.frame(barcodes$tss, barcodes$portion, barcodes$plate, barcodes$condition)

mirna.filt.counts <- filtered.data(mirna.raw.counts, factor = myfactors$barcodes.condition, norm = FALSE, depth = NULL, method = 1, cv.cutoff = 500, cpm = 0.5, p.adj = "fdr") # 501 features (26.6%) are to be kept for differential expression analysis with filtering method 1

# save(mirna.filt.counts, file = "data/cooked/miRNA-Seq/miRNA.filt.rda")

#-------- checking for possible GC content and length bias with NOISeq --------#
GCcontent <- read.csv("results/preprocessing/cookingmiRNASeq/miRNA.biomart.txt", sep = "\t")

mygc = c(GCcontent$Gene...GC.content)
names(mygc) = c(GCcontent$miRBase.ID)
mygc = mygc[rownames(mirna.filt.counts)]
names(mygc) = rownames(mirna.filt.counts)

mygc <- na.omit(mygc) # 10 mirnas without information

mirna.filt.counts <- mirna.filt.counts[names(mygc),]

# save(mygc, file = "results/preprocessing/cookingmiRNASeq/GC.miRNA.rda")

mymirnadata.filt <- NOISeq::readData(data = mirna.filt.counts, factors = myfactors, gc = mygc)

source(file = "scripts/preprocessing/GCbias.R")
GC.plot(GC.dat(mymirnadata.filt, factor = "barcodes.condition"))

mymirnaPCA = dat(mymirnadata.filt, type = "PCA")
par(cex = 0.75)
explo.plot(mymirnaPCA, factor = "barcodes.condition", plottype = "scores")
explo.plot(mymirnaPCA, factor = "barcodes.condition", plottype = "loadings")
explo.plot(mymirnaPCA, factor = "barcodes.tss")
explo.plot(mymirnaPCA, factor = "barcodes.portion")
explo.plot(mymirnaPCA, factor = "barcodes.plate")

#--------------- normalizing for GC content with EDASeq ---------------#
load("results/preprocessing/cookingmiRNASeq/GC.miRNA.rda")
feature <- data.frame(gc=mygc)

data <- newSeqExpressionSet(counts=as.matrix(mirna.filt.counts), featureData=feature, phenoData=data.frame(conditions=barcodes$condition, row.names=barcodes$barcode))

dataWithin <- withinLaneNormalization(data, "gc", which="full")
mirna.norm <- betweenLaneNormalization(dataWithin, which="full")

# save(mirna.norm, file = "data/cooked/miRNA-Seq/miRNA.norm.rda")

# Extract the offset, which will be input directly into DEseq2 to normalise the counts. 
normFactors <- withinLaneNormalization(data,"gc",
                                       which="full", offset=TRUE)
normFactors <- betweenLaneNormalization(normFactors,
                                        which="full", offset=TRUE)

normFactors <- exp(-1 * normFactors@assayData$offset)
# save(normFactors, file = "results/preprocessing/cookingmiRNASeq/miRNA.normFactors.rda")

# Extract normalized counts to check for bias on NOISeq
mirna.counts <- mirna.norm@assayData$normalizedCounts
mirna.counts <- as.data.frame(mirna.counts)

#-------- checking if it fixed GC content and length bias with NOISeq ---------#
mymirnadata.norm<- NOISeq::readData(data = mirna.counts, factors = myfactors, gc = mygc)
GC.plot(GC.dat(mymirnadata.norm, factor = "barcodes.condition"))

mymirnaPCA = dat(mymirnadata.norm, type = "PCA", logtransf = TRUE)
par(cex = 0.75)
explo.plot(mymirnaPCA, factor = "barcodes.condition", plottype = "scores")
#------------------------------- DEA with DESeq2 ------------------------------#
barcodes <- get_IDs(mirna.filt.counts)
barcodes$condition <- as.factor(barcodes$condition)
barcodes$condition <- relevel(barcodes$condition, ref = "normal")
myfactors <- data.frame(barcodes$tss, barcodes$portion, barcodes$plate, barcodes$condition)

dds <- DESeqDataSetFromMatrix(countData = mirna.filt.counts,
                              colData = barcodes,
                              design = ~ condition)

load("data/cooked/miRNA-Seq/miRNA.normFactors.rda")
normFactors <- normFactors / exp(rowMeans(log(normFactors)))
normalizationFactors(dds) <- normFactors

dds <- DESeq(dds)
res.deseq2 <- results(dds, alpha = 0.05)

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

# write.table(activated.genes.deseq2, file = "results/preprocessing/cookingmiRNASeq/DESeq2.up.txt", row.names = FALSE, col.names = FALSE, quote = FALSE)

# write.table(repressed.genes.deseq2, file = "results/preprocessing/cookingmiRNASeq/DESeq2.down.txt", row.names = FALSE, col.names = FALSE, quote = FALSE)

resOrdered <- res.deseq2[order(res.deseq2$padj),]
head(resOrdered)
resOrderedDF <- as.data.frame(resOrdered)
resOrderedDF <- na.omit(resOrderedDF)
# write.table(resOrderedDF, file = "results/preprocessing/cookingmiRNASeq/DESeq2.ordered.csv", row.names=TRUE, col.names=TRUE, sep="\t", quote=FALSE)

#------------------------------- DEA with limma-voom -------------------------------#
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

# write.table(activated.genes.limma.voom, file = "results/preprocessing/cookingmiRNASeq/limma.voom.up.txt", row.names = FALSE, col.names = FALSE, quote = FALSE)

# write.table(repressed.genes.limma.voom, file = "results/preprocessing/cookingmiRNASeq/limma.voom.down.txt", row.names = FALSE, col.names = FALSE, quote = FALSE)

topOrdered <- top.limma.voom[order(top.limma.voom$adj.P.Val),]
topOrderedDF <- as.data.frame(topOrdered)
topOrderedDF <- na.omit(topOrderedDF)
# write.table(topOrderedDF, file = "results/preprocessing/cookingmiRNASeq/limma.voom.ordered.csv", row.names=TRUE, col.names=TRUE, sep="\t", quote=FALSE)

#------------------------------- DEA with edgeR -------------------------------#
y <- DGEList(counts = mirna.filt.counts, lib.size = colSums(mirna.filt.counts), group = barcodes$condition, genes = rownames(mirna.filt.counts))

y$offset <- normFactors

design <- model.matrix(~ barcodes$condition)

y <- estimateCommonDisp(y, design = design)
y <- estimateTagwiseDisp(y, design = design)
plotBCV(y)

et <- exactTest(y) # performs pair-wise tests for differential expression between two groups
top.edger <- topTags(et, n = Inf) # takes the output from exactTest(), adjusts the raw p-values using the False Discovery Rate (FDR) correction, and returns the top differentially expressed genes

topSig <- top.edger[top.edger$table$FDR < 0.05, ] # we select DEGs with alpha=0.05
dim(topSig)
topSig <- topSig[abs(top.edger$table$logFC) >= 1, ] # we filter the output of dataDEGs by abs(LogFC) >=1
dim(topSig)

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

# write.table(activated.genes.edger, file = "results/preprocessing/cookingmiRNASeq/edgeR.up.txt", row.names = FALSE, col.names = FALSE, quote = FALSE)

# write.table(repressed.genes.edger, file = "results/preprocessing/cookingmiRNASeq/edgeR.down.txt", row.names = FALSE, col.names = FALSE, quote = FALSE)

top.edger <- as.data.frame(top.edger)

topOrdered <- top.edger[order(top.edger$FDR),]
topOrderedDF <- as.data.frame(topOrdered)
topOrderedDF <- na.omit(topOrderedDF)
# write.table(topOrderedDF, file = "results/preprocessing/cookingmiRNASeq/edgeR.ordered.csv", row.names=TRUE, col.names=TRUE, sep="\t", quote=FALSE)

#----------------------------- intersecting DEGs ------------------------------#
# with DESeq2, limma-voom and edgeR
common.activated <- intersect(intersect(activated.genes.deseq2, activated.genes.edger), activated.genes.limma.voom) # 66
common.repressed <- intersect(intersect(repressed.genes.deseq2, repressed.genes.edger), repressed.genes.limma.voom) # 53

# write.table(common.activated, file = "results/preprocessing/cookingmiRNASeq/common.up.txt", row.names = FALSE, col.names = FALSE, quote = FALSE)

# write.table(common.repressed, file = "results/preprocessing/cookingmiRNASeq/common.down.txt", row.names = FALSE, col.names = FALSE, quote = FALSE)