#------------------------- loading required packages --------------------------#
library(limma)
library(ggplot2)
library(dplyr)
library(NOISeq)

load("data/raw/prot/prot.rda")

#------------------------------- filtering NAs --------------------------------#

missing.per.patient <- as.data.frame(sapply(prot, function(x) sum(is.na(x))))
missing.per.patient$percentage <- missing.per.patient$`sapply(prot, function(x) sum(is.na(x)))`*100/487
missing.per.peptide <- as.data.frame(rowSums(is.na(prot)))
missing.per.peptide$percentage <- missing.per.peptide$`rowSums(is.na(prot))`*100/643

patients.to.remove <- rownames(missing.per.patient[which(missing.per.patient$percentage >= 20), ]) # 18

peptides.to.remove <- rownames(missing.per.peptide[which(missing.per.peptide$percentage >= 20), ]) # 30

prot[, patients.to.remove] <- NULL

prot <- prot[!(row.names(prot) %in% peptides.to.remove), ]

#-------------------- converting antibodies to Entrez IDs ---------------------#
library(dplyr)
antibodies <- read.table("data/raw/prot/RPPA_Antibodies.txt", sep = "\t", header = TRUE, row.names = 1)
rownames(antibodies)<-gsub("-","",as.character(rownames(antibodies)))
rownames(antibodies)<-gsub("_","",as.character(rownames(antibodies)))
rownames(antibodies)<-toupper(rownames(antibodies))
rownames(prot)<-gsub("-","",as.character(rownames(prot)))
rownames(prot)<-gsub("_","",as.character(rownames(prot)))
rownames(prot)<-toupper(rownames(prot))

existing.peptides <- intersect(rownames(prot), rownames(antibodies)) # 437

missing.peptides <- setdiff(rownames(prot), rownames(antibodies)) # 20

prot <- prot[!(row.names(prot) %in% missing.peptides), ]

prot$Entrez <- antibodies[existing.peptides, ]$NCBI_Entrez_Gene_ID

# write.table(prot, file = "data/cooked/prot/prot.txt", sep = "\t", quote = FALSE)

# after deleting duplicated Entrez IDs and keeping the most reliable peptide for it
load("data/cooked/prot/prot.filt.rda")

#---------------------------- exploring with PCA ------------------------------#
groups <- substr(colnames(prot), 13, 15)
groups <- gsub("-01", "tumor", groups)
groups <- gsub("-11", "normal", groups)

groups <- as.data.frame(groups)

mydata <- NOISeq::readData(data = prot, factors = groups)

myPCA = dat(mydata, type = "PCA", logtransf = TRUE, norm = TRUE)
par(cex = 0.75)
explo.plot(myPCA, factor = "groups", plottype = "scores")
explo.plot(myPCA, factor = "groups", plottype = "loadings")

#------------------------------- DEA with limma -------------------------------#
groups <- substr(colnames(prot), 13, 15)
groups <- gsub("-01", "tumor", groups)
groups <- gsub("-11", "normal", groups)

design <- model.matrix(~ groups)

fit1 <- lmFit(prot, design)

fit2 <- eBayes(fit1)

top.limma <- topTable(fit2, coef = 2, number = Inf)

log.fold.change <- top.limma$logFC
q.value <- top.limma$adj.P.Val
genes.ids <- rownames(top.limma)
names(log.fold.change) <- genes.ids
names(q.value) <- genes.ids

activated.genes.limma <- genes.ids[log.fold.change > 0.5 & q.value < 0.05]
repressed.genes.limma <- genes.ids[log.fold.change < -0.5 & q.value < 0.05]

length(activated.genes.limma) # 64
length(repressed.genes.limma) # 28

log.q.val <- -log10(q.value)
plot(log.fold.change,log.q.val,pch=19,col="grey",cex=0.8,
     xlim=c(-3,3),ylim = c(0,40),
     xlab="log2(Fold-change)",ylab="-log10(q-value)",cex.lab=1.5)
points(x = log.fold.change[activated.genes.limma],
       y = log.q.val[activated.genes.limma],col="red",cex=0.8,pch=19)
points(x = log.fold.change[repressed.genes.limma],
       y = log.q.val[repressed.genes.limma],col="blue",cex=0.8,pch=19)

# write.table(activated.genes.limma, file = "results/preprocessing/cookingProt/limma.up.txt", row.names = FALSE, col.names = FALSE, quote = FALSE)

# write.table(repressed.genes.limma, file = "results/preprocessing/cookingProt/limma.down.txt", row.names = FALSE, col.names = FALSE, quote = FALSE)

topOrdered <- top.limma[order(top.limma$adj.P.Val),]
topOrderedDF <- as.data.frame(topOrdered)
topOrderedDF <- na.omit(topOrderedDF)
# write.table(topOrderedDF, file = "results/preprocessing/cookingProt/limma.ordered.csv", row.names=TRUE, col.names=TRUE, sep="\t", quote=FALSE)