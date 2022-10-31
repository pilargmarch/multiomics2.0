#------------------------- loading required packages --------------------------#
library(TCGAbiolinks)
library(SummarizedExperiment)
library(limma)
library(ggplot2)

load("data/raw/met/met.rda")

met.beta.values <- as.data.frame(assay(met))
met.genes.info <- as.data.frame(rowRanges(met))
met.sample.info <- as.data.frame(met@colData)

#---------------------- filtering NA's and all 0 probes -----------------------#

probe.na <- rowSums(is.na(met.beta.values))
table(probe.na == 0)

probe <- probe.na[probe.na == 0]
met.beta.values <- met.beta.values[row.names(met.beta.values) %in% names(probe), ]

probe.0 <- rowSums(met.beta.values)
table(probe.0 == 0)

met.beta.values[,c("TCGA-A7-A26F-01A-21D-A16A-05", "TCGA-A7-A26J-01A-11D-A16A-05", "TCGA-A7-A26J-01A-11D-A27B-05", "TCGA-A7-A13G-01A-11D-A13K-05", "TCGA-A7-A26E-01B-06D-A27B-05", "TCGA-A7-A26E-01A-11D-A16A-05", "TCGA-B6-A1KC-01A-11D-A13K-05", "TCGA-AC-A2QH-01B-04D-A22R-05", "TCGA-AC-A3OD-01A-11D-A21R-05")] <- NULL

#--------------------- converting beta-values to M-values ---------------------#
bval <- met.beta.values
saveRDS(bval, file = "data/cooked/met/bval.RDS", compress = FALSE)

mval <- t(apply(met.beta.values, 1, function(x) log2(x/(1-x))))
saveRDS(mval, file = "data/cooked/met/mval.RDS", compress = FALSE)

# to load
mval <- readRDS("data/cooked/met/mval.RDS")
bval <- readRDS("data/cooked/met/bval.RDS")

#-------------------- exploring with PCA and density plot ---------------------#
pc = prcomp(t(bval[floor(runif(100000, min=1, max=364019)),]), scale = FALSE)
bcodes <- substr(rownames(pc$x), 1, 15)
tumor <- substr(bcodes, 14, 15) == "01"
loads <- round(pc$sdev^2/sum(pc$sdev)*100, 1)
xlab <- c(paste("PC1", loads[1], "%"))
ylab <- c(paste("PC2", loads[2], "%"))
plot(pc$x[ , 1], pc$x[ , 2], xlab = xlab, ylab = ylab, type = "n",
     main = "PC1 vs . PC2 ", cex.axis = 1.5 , cex.lab = 1.3 , cex.main = 1.5)
points(pc$x[tumor, 1] , pc$x[tumor, 2], col = "red", pch = 16)
points(pc$x[!tumor, 1] , pc$x[!tumor, 2], col = "blue", pch = 16)
legend("topright", legend = c("Tumor" , "Normal" ), pch = 16, col = c("red", "blue"))

rand_smp = sample(rownames(bval), 10000, replace = F)
bcodes = substr(colnames(bval), 1, 15)
tumor <- substr(bcodes, 14, 15) == "01"

plot(density(bval[rand_smp, 1]), type = "n", ylim = c(0,3), main = "Density of B-values for tumor and normal samples", cex.axis = 1.5, xlab = "Beta-value", ylab = "Density", cex.main = 1.5)
for (i in 1:length(bcodes)) {
  if (tumor[i] == T) {
    points(density(bval[rand_smp, i]), col = "red", type = "l")
  }
  else {
    points(density(bval[rand_smp, i]), col = "blue", type = "l")
  }
}
legend("topright", legend = c("Tumor", "Normal"), col = c("red", "blue"), lty = 1)

#------------------------------- DEA with limma -------------------------------#

met.sample.info <- met.sample.info[-c(30, 188, 189, 192, 198, 199, 309, 314, 794), ]

group <- as.factor(met.sample.info$sample_type)

group <- relevel(group, ref = "Solid Tissue Normal")

design <- model.matrix(~ group)

fit <- lmFit(mval, design)
fit2 <- eBayes(fit)
top <- topTable(fit2, coef=ncol(design), number = Inf)

log.fold.change <- top$logFC
q.value <- top$adj.P.Val
genes.ids <- rownames(top)
names(log.fold.change) <- genes.ids
names(q.value) <- genes.ids

activated.genes.limma <- genes.ids[log.fold.change > 2 & q.value < 0.05]
repressed.genes.limma <- genes.ids[log.fold.change < -2 & q.value < 0.05]

length(activated.genes.limma) # 4400
length(repressed.genes.limma) # 8763

log.q.val <- -log10(q.value)
plot(log.fold.change,log.q.val,pch=19,col="grey",cex=0.8,
     xlim=c(-5,5),ylim = c(0,130),
     xlab="log2(Fold-change)",ylab="-log10(q-value)",cex.lab=1.5)
points(x = log.fold.change[activated.genes.limma],
       y = log.q.val[activated.genes.limma],col="red",cex=0.8,pch=19)
points(x = log.fold.change[repressed.genes.limma],
       y = log.q.val[repressed.genes.limma],col="blue",cex=0.8,pch=19)

# write.table(activated.genes.limma, file = "results/preprocessing/cookingMet/limma.up.txt", row.names = FALSE, col.names = FALSE, quote = FALSE)

# write.table(repressed.genes.limma, file = "results/preprocessing/cookingMet/limma.down.txt", row.names = FALSE, col.names = FALSE, quote = FALSE)

topOrdered <- top[order(top$adj.P.Val),]
topOrderedDF <- as.data.frame(topOrdered)
topOrderedDF <- na.omit(topOrderedDF)
# write.table(topOrderedDF, file = "results/preprocessing/cookingMet/limma.ordered.csv", row.names=TRUE, col.names=TRUE, sep="\t", quote=FALSE)