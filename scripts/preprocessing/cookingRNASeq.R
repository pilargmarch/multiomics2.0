#------------------------- loading required packages --------------------------#
library(TCGAbiolinks)
library(SummarizedExperiment)
library(NOISeq)
library(cqn)
library(edgeR)
library(limma)
library(DESeq2)
library(AnnotationDbi)
library(org.Hs.eg.db)
library(ggplot2)

#--------------- filtering with NOISeq with CPM threshold = 0.5 ---------------#
load("data/raw/RNA-Seq/RNA.rda")
rna.raw.counts <- as.data.frame(assay(rna))
rna.genes.info <- as.data.frame(rowRanges(rna))
rna.sample.info <- as.data.frame(colData(rna))

barcodes <- get_IDs(rna)
myfactors <- data.frame(barcodes$tss, barcodes$portion, barcodes$plate, barcodes$condition)
head(myfactors)

GCcontent <- read.csv("results/preprocessing/cookingRNASeq/genes.biomart.txt", sep = "\t")
colnames(GCcontent) <- c("gene_id", "gc_content")
# problem: we have different versions so we need to merge by ID only
# so we create a new column in rna.genes.info with the gene ID without version number, and we do the same in GCcontent
rna.genes.info$gene_id_no_version <- sub('\\.[0-9]*$', '', rna.genes.info$gene_id)
GCcontent$gene_id_no_version <- sub('\\.[0-9]*$', '', GCcontent$gene_id)
# complete rna.genes.info with GC content
rna.genes.info <- merge(rna.genes.info, GCcontent, by = "gene_id_no_version") # merging turns gene_id into gene_id.x and the gene_id from GCcontent into gene_id.y, so we need to remove those
rna.genes.info$gene_id <- rna.genes.info$gene_id.x
rna.genes.info$gene_id.x <- NULL
rna.genes.info$gene_id.y <- NULL

mygc = c(rna.genes.info$gc_content)
names(mygc) = rna.genes.info$gene_id

mylength = c(rna.genes.info$width)
names(mylength) = rna.genes.info$gene_id

mybiotypes = c(rna.genes.info$gene_type)
names(mybiotypes) = rna.genes.info$gene_id

mychroms = data.frame(rna.genes.info$seqnames, rna.genes.info$start, rna.genes.info$end)
rownames(mychroms) = rna.genes.info$gene_id
colnames(mychroms) <- c("Chr", "GeneStart", "GeneEnd")

mydata <- NOISeq::readData(data = rna.raw.counts, factors = myfactors, length = mylength, gc = mygc, biotype = mybiotypes, chromosome = mychroms)

rna.filt.counts <- filtered.data(rna.raw.counts, factor = myfactors$barcodes.condition, norm = FALSE, depth = NULL, method = 1, cv.cutoff = 500, cpm = 0.5, p.adj = "fdr") # 19350 features (32%) are to be kept for differential expression analysis

# save(rna.filt.counts, file = "data/cooked/RNA-Seq/RNA.filt.rda")

#-------- checking for possible GC content and length bias with NOISeq --------#
myexpdata.filt <- NOISeq::readData(data = rna.filt.counts, factors = myfactors, length = mylength, gc = mygc, biotype = mybiotypes, chromosome = mychroms)

myexplengthbias.filt = dat(myexpdata.filt, factor = "barcodes.condition", type = "lengthbias")
explo.plot(myexplengthbias.filt, samples = NULL, toplot = "global")

myexpGCbias.filt = dat(myexpdata.filt, factor = "barcodes.condition", type = "GCbias")
explo.plot(myexpGCbias.filt, samples = NULL, toplot = "global")

myexpPCA = dat(myexpdata.filt, type = "PCA")
par(cex = 0.75)
explo.plot(myexpPCA, factor = "barcodes.condition", plottype = "scores")
explo.plot(myexpPCA, factor = "barcodes.condition", plottype = "loadings")

#--------------- normalizing for GC content and length with cqn ---------------#
load("results/preprocessing/cookingRNASeq/GC.length.RNA.rda")
gc.length.rna <- as.data.frame(gc.length.rna)
sizeFactors.rna <- colSums(rna.filt.counts)

rna.norm <- cqn(rna.filt.counts, lengths = gc.length.rna$length, x = gc.length.rna$gc, sizeFactors = sizeFactors.rna, verbose = TRUE)
# save(rna.norm, file = "data/cooked/RNA-Seq/RNA.norm.rda")

cqnOffset <- rna.norm$glm.offset
cqnNormFactors <- exp(cqnOffset)
# save(cqnNormFactors, file = "data/cooked/RNA-Seq/RNA.normFactors.rda")

rna.norm.expression <- rna.norm$y + rna.norm$offset
rna.norm.expression <- as.data.frame(rna.norm.expression)

#-------- checking if it fixed GC content and length bias with NOISeq ---------#
myexpdata.norm <- NOISeq::readData(data = rna.norm.expression, factors = myfactors, length = mylength.norm, gc = mygc.norm, biotype = mybiotypes.norm, chromosome = mychroms.norm)

myexplengthbias.norm = dat(myexpdata.norm, factor = "barcodes.condition", type = "lengthbias")
explo.plot(myexplengthbias.norm, samples = NULL, toplot = "global")

myexpGCbias.norm = dat(myexpdata.norm, factor = "barcodes.condition", type = "GCbias")
explo.plot(myexpGCbias.norm, samples = NULL, toplot = "global")

myexpPCA.norm = dat(myexpdata.norm, type = "PCA", norm = TRUE, logtransf = TRUE)
par(cex = 0.75)
explo.plot(myexpPCA.norm, factor = "barcodes.condition", plottype = "scores")
explo.plot(myexpPCA.norm, factor = "barcodes.condition", plottype = "loadings")

barcodes$condition <- as.factor(barcodes$condition)
barcodes$condition <- relevel(barcodes$condition, ref = "normal")

RPM <- sweep(log2(rna.filt.counts + 1), 2, log2(sizeFactors.rna/10^6))
RPKM.std <- sweep(RPM, 1, log2(gc.length.rna$length / 10^3))

whGenes <- which(rowMeans(RPKM.std) >= 2 & gc.length.rna$length >= 100)
M.std <- rowMeans(RPKM.std[whGenes, which(barcodes$condition == "cancer")]) - rowMeans(RPKM.std[whGenes, which(barcodes$condition == "normal")])
A.std <- rowMeans(RPKM.std[whGenes,])
M.cqn <- rowMeans(rna.norm.expression[whGenes, which(barcodes$condition == "cancer")]) - rowMeans(rna.norm.expression[whGenes, which(barcodes$condition == "normal")])
A.cqn <- rowMeans(rna.norm.expression[whGenes,])

par(mfrow = c(1,2))
plot(A.std, M.std, cex = 0.5, pch = 16, xlab = "A", ylab = "M", 
     main = "Standard RPKM", ylim = c(-4,4), xlim = c(0,12), 
     col = alpha("black", 0.25))
plot(A.cqn, M.cqn, cex = 0.5, pch = 16, xlab = "A", ylab = "M", 
     main = "CQN normalized RPKM", ylim = c(-4,4), xlim = c(0,12), 
     col = alpha("black", 0.25))

# We can plot the effect of GC and length
par(mfrow=c(1,2))
cqnplot(rna.norm, n = 1, xlab = "GC content", lty = 1, ylim = c(1,7))
cqnplot(rna.norm, n = 2, xlab = "length", lty = 1, ylim = c(1,7))

#------------------------------- DEA with DESeq2 ------------------------------#
dds <- DESeqDataSetFromMatrix(countData = rna.filt.counts,
                              colData = rna.sample.info,
                              design = ~ barcodes.condition)
dds$barcodes.condition <- relevel(dds$barcodes.condition, ref = "normal")

load("data/cooked/RNA-Seq/RNA.normFactors.rda")
normFactors <- normFactors / exp(rowMeans(log(normFactors)))
normalizationFactors(dds) <- normFactors
dds <- DESeq(dds)
res <- results(dds, alpha = 0.05)

log.fold.change <- res$log2FoldChange
q.value <- res$padj
genes.ids <- rownames(rna.filt.counts)
names(log.fold.change) <- genes.ids
names(q.value) <- genes.ids
activated.genes.deseq2 <- genes.ids[log.fold.change > 2 & q.value < 0.05]
activated.genes.deseq2 <- activated.genes.deseq2[!is.na(activated.genes.deseq2)]
repressed.genes.deseq2 <- genes.ids[log.fold.change < - 2 & q.value < 0.05]
repressed.genes.deseq2 <- repressed.genes.deseq2[!is.na(repressed.genes.deseq2)]
length(activated.genes.deseq2) # 670
length(repressed.genes.deseq2) # 1148

log.q.val <- -log10(q.value)
plot(log.fold.change,log.q.val,pch=19,col="grey",cex=0.8,
     xlim=c(-8,8),ylim = c(0,240),
     xlab="log2(Fold-change)",ylab="-log10(q-value)",cex.lab=1.5)
points(x = log.fold.change[activated.genes.deseq2],
       y = log.q.val[activated.genes.deseq2],col="red",cex=0.8,pch=19)
points(x = log.fold.change[repressed.genes.deseq2],
       y = log.q.val[repressed.genes.deseq2],col="blue",cex=0.8,pch=19)

ens.str <- substr(rownames(res), 1, 15)
res$symbol <- mapIds(org.Hs.eg.db,
                     keys=ens.str,
                     column="SYMBOL",
                     keytype="ENSEMBL",
                     multiVals="first")
res$entrez <- mapIds(org.Hs.eg.db,
                     keys=ens.str,
                     column="ENTREZID",
                     keytype="ENSEMBL",
                     multiVals="first")

activated.genes.deseq2 <- as.data.frame(activated.genes.deseq2)
ens.str <- substr(activated.genes.deseq2$activated.genes.deseq2, 1, 15)
activated.genes.deseq2$entrez <- mapIds(org.Hs.eg.db,
                                        keys=ens.str,
                                        column="ENTREZID",
                                        keytype="ENSEMBL",
                                        multiVals="first")
activated.genes.deseq2 <- na.omit(activated.genes.deseq2)

repressed.genes.deseq2 <- as.data.frame(repressed.genes.deseq2)
ens.str <- substr(repressed.genes.deseq2$repressed.genes.deseq2, 1, 15)
repressed.genes.deseq2$entrez <- mapIds(org.Hs.eg.db,
                                        keys=ens.str,
                                        column="ENTREZID",
                                        keytype="ENSEMBL",
                                        multiVals="first")
repressed.genes.deseq2 <- na.omit(repressed.genes.deseq2)

# write.table(activated.genes.deseq2$entrez, file = "results/preprocessing/cookingRNASeq/DESeq2.up.txt", row.names = FALSE, col.names = FALSE, quote = FALSE)

# write.table(repressed.genes.deseq2$entrez, file = "results/preprocessing/cookingRNASeq/DESeq2.down.txt", row.names = FALSE, col.names = FALSE, quote = FALSE)

resOrdered <- res[order(res$padj),]
head(resOrdered)
resOrderedDF <- as.data.frame(resOrdered)
resOrderedDF <- na.omit(resOrderedDF)
# write.table(resOrderedDF, file = "results/preprocessing/cookingRNASeq/DESeq2.ordered.csv", row.names=TRUE, col.names=TRUE, sep="\t", quote=FALSE)

#------------------------------- DEA with limma -------------------------------#
load("data/cooked/RNA-Seq/RNA.norm.rda")
rna.norm.expression <- rna.norm$y + rna.norm$offset

design <- model.matrix(~ barcodes$condition)

fit1 <- lmFit(rna.norm.expression, design)

fit2 <- eBayes(fit1)

top <- topTable(fit2, coef = 2, number = Inf)

log.fold.change <- top$logFC
q.value <- top$adj.P.Val
genes.ids <- rownames(top)
names(log.fold.change) <- genes.ids
names(q.value) <- genes.ids

activated.genes.limma <- genes.ids[log.fold.change > 2 & q.value < 0.05]
repressed.genes.limma <- genes.ids[log.fold.change < -2 & q.value < 0.05]

length(activated.genes.limma) # 571
length(repressed.genes.limma) # 933

log.q.val <- -log10(q.value)
plot(log.fold.change,log.q.val,pch=19,col="grey",cex=0.8,
     xlim=c(-8,8),ylim = c(0,160),
     xlab="log2(Fold-change)",ylab="-log10(q-value)",cex.lab=1.5)
points(x = log.fold.change[activated.genes.limma],
       y = log.q.val[activated.genes.limma],col="red",cex=0.8,pch=19)
points(x = log.fold.change[repressed.genes.limma],
       y = log.q.val[repressed.genes.limma],col="blue",cex=0.8,pch=19)

ens.str <- substr(rownames(top), 1, 15)
top$symbol <- mapIds(org.Hs.eg.db,
                     keys=ens.str,
                     column="SYMBOL",
                     keytype="ENSEMBL",
                     multiVals="first")
top$entrez <- mapIds(org.Hs.eg.db,
                     keys=ens.str,
                     column="ENTREZID",
                     keytype="ENSEMBL",
                     multiVals="first")

activated.genes.limma <- as.data.frame(activated.genes.limma)
ens.str <- substr(activated.genes.limma$activated.genes.limma, 1, 15)
activated.genes.limma$entrez <- mapIds(org.Hs.eg.db,
                                       keys=ens.str,
                                       column="ENTREZID",
                                       keytype="ENSEMBL",
                                       multiVals="first")
activated.genes.limma <- na.omit(activated.genes.limma)

repressed.genes.limma <- as.data.frame(repressed.genes.limma)
ens.str <- substr(repressed.genes.limma$repressed.genes.limma, 1, 15)
repressed.genes.limma$entrez <- mapIds(org.Hs.eg.db,
                                       keys=ens.str,
                                       column="ENTREZID",
                                       keytype="ENSEMBL",
                                       multiVals="first")
repressed.genes.limma <- na.omit(repressed.genes.limma)

# write.table(activated.genes.limma$entrez, file = "results/preprocessing/cookingRNASeq/limma.up.txt", row.names = FALSE, col.names = FALSE, quote = FALSE)

# write.table(repressed.genes.limma$entrez, file = "results/preprocessing/cookingRNASeq/limma.down.txt", row.names = FALSE, col.names = FALSE, quote = FALSE)

topOrdered <- top[order(top$adj.P.Val),]
topOrderedDF <- as.data.frame(topOrdered)
topOrderedDF <- na.omit(topOrderedDF)
# write.table(topOrderedDF, file = "results/preprocessing/cookingRNASeq/limma.ordered.csv", row.names=TRUE, col.names=TRUE, sep="\t", quote=FALSE)

#------------------------------- DEA with edgeR -------------------------------#
y <- DGEList(counts = rna.filt.counts, lib.size = colSums(rna.filt.counts), group = barcodes$condition, genes = rownames(rna.filt.counts))

y$offset <- normFactors

design <- model.matrix(~ barcodes$condition)

y <- estimateCommonDisp(y, design = design)
y <- estimateTagwiseDisp(y, design = design)
plotBCV(y)

et <- exactTest(y) # performs pair-wise tests for differential expression between two groups
top <- topTags(et, n = Inf) # takes the output from exactTest(), adjusts the raw p-values using the False Discovery Rate (FDR) correction, and returns the top differentially expressed genes

topSig <- top[top$table$FDR < 0.05, ] # we select DEGs with alpha=0.05
dim(topSig)
topSig <- topSig[abs(top$table$logFC) >= 2, ] # we filter the output of dataDEGs by abs(LogFC) >=2
dim(topSig)

# this is equivalent to doing
de <- (decideTestsDGE(et, lfc = 2, p.value = 0.05))
summary(de)
detags <- rownames(y)[as.logical(de)]
plotSmear(et, de.tags=detags, main="plotSmear")
abline(h=c(-2,2), col="blue")

activated.genes.edger <- topSig$table$genes[topSig$table$logFC > 0]
length(activated.genes.edger) # 930
repressed.genes.edger <- topSig$table$genes[topSig$table$logFC < 0]
length(repressed.genes.edger) # 661

top <- top$table
log.fold.change <- top$logFC
q.value <- top$FDR
genes.ids <- rownames(rna.filt.counts)
names(log.fold.change) <- genes.ids
names(q.value) <- genes.ids
activated.genes.edger <- genes.ids[log.fold.change > 2 & q.value < 0.05]
activated.genes.edger <- activated.genes.edger[!is.na(activated.genes.edger)]
repressed.genes.edger <- genes.ids[log.fold.change < - 2 & q.value < 0.05]
repressed.genes.edger <- repressed.genes.edger[!is.na(repressed.genes.edger)]

log.q.val <- -log10(q.value)
plot(log.fold.change,log.q.val,pch=19,col="grey",cex=0.8,
     xlim=c(-8,8),ylim = c(0,240),
     xlab="log2(Fold-change)",ylab="-log10(q-value)",cex.lab=1.5)
points(x = log.fold.change[activated.genes.edger],
       y = log.q.val[activated.genes.edger],col="red",cex=0.8,pch=19)
points(x = log.fold.change[repressed.genes.edger],
       y = log.q.val[repressed.genes.edger],col="blue",cex=0.8,pch=19)

top <- as.data.frame(top)
ens.str <- substr(rownames(top), 1, 15)
top$symbol <- mapIds(org.Hs.eg.db,
                     keys=ens.str,
                     column="SYMBOL",
                     keytype="ENSEMBL",
                     multiVals="first")
top$entrez <- mapIds(org.Hs.eg.db,
                     keys=ens.str,
                     column="ENTREZID",
                     keytype="ENSEMBL",
                     multiVals="first")

activated.genes.edger <- as.data.frame(activated.genes.edger)
ens.str <- substr(activated.genes.edger$activated.genes.edger, 1, 15)
activated.genes.edger$entrez <- mapIds(org.Hs.eg.db,
                                       keys=ens.str,
                                       column="ENTREZID",
                                       keytype="ENSEMBL",
                                       multiVals="first")
activated.genes.edger <- na.omit(activated.genes.edger)

repressed.genes.edger <- as.data.frame(repressed.genes.edger)
ens.str <- substr(repressed.genes.edger$repressed.genes.edger, 1, 15)
repressed.genes.edger$entrez <- mapIds(org.Hs.eg.db,
                                       keys=ens.str,
                                       column="ENTREZID",
                                       keytype="ENSEMBL",
                                       multiVals="first")
repressed.genes.edger <- na.omit(repressed.genes.edger)

# write.table(activated.genes.edger$entrez, file = "results/preprocessing/cookingRNASeq/edgeR.up.txt", row.names = FALSE, col.names = FALSE, quote = FALSE)

# write.table(repressed.genes.edger$entrez, file = "results/preprocessing/cookingRNASeq/edgeR.down.txt", row.names = FALSE, col.names = FALSE, quote = FALSE)

topOrdered <- top[order(top$FDR),]
topOrderedDF <- as.data.frame(topOrdered)
topOrderedDF <- na.omit(topOrderedDF)
# write.table(topOrderedDF, file = "results/preprocessing/cookingRNASeq/edgeR.ordered.csv", row.names=TRUE, col.names=TRUE, sep="\t", quote=FALSE)

#----------------------------- intersecting DEGs ------------------------------#
activated.genes.deseq2 <- read.table(file = "results/preprocessing/cookingRNASeq/DESeq2.up.txt")
activated.genes.deseq2 <- as.vector(activated.genes.deseq2$V1)

repressed.genes.deseq2 <- read.table(file = "results/preprocessing/cookingRNASeq/DESeq2.down.txt")
repressed.genes.deseq2 <- as.vector(repressed.genes.deseq2$V1)

activated.genes.limma <- read.table(file = "results/preprocessing/cookingRNASeq/limma.up.txt")
activated.genes.limma <- as.vector(activated.genes.limma$V1)

repressed.genes.limma <- read.table(file = "results/preprocessing/cookingRNASeq/limma.down.txt")
repressed.genes.limma <- as.vector(repressed.genes.limma$V1)

activated.genes.edger <- read.table(file = "results/preprocessing/cookingRNASeq/edgeR.up.txt")
activated.genes.edger <- as.vector(activated.genes.edger$V1)

repressed.genes.edger <- read.table(file = "results/preprocessing/cookingRNASeq/edgeR.down.txt")
repressed.genes.edger <- as.vector(repressed.genes.edger$V1)

common.activated <- intersect(intersect(activated.genes.deseq2, activated.genes.edger), activated.genes.limma) 
length(common.activated) # 452

common.repressed <- intersect(intersect(repressed.genes.deseq2, repressed.genes.edger), repressed.genes.limma) 
length(common.repressed) # 473

write.table(common.activated, file = "results/preprocessing/cookingRNASeq/all.up.txt", row.names = FALSE, col.names = FALSE, quote = FALSE)

write.table(common.repressed, file = "results/preprocessing/cookingRNASeq/all.down.txt", row.names = FALSE, col.names = FALSE, quote = FALSE)