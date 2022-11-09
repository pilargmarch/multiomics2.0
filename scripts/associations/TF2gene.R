#------------------------- loading required packages --------------------------#
library(biomaRt)
library(tidyr)

hTFtarget <- read.table(file = "reports/associations/TF-gene/TF-Target-information.hTFtarget.txt", sep = "\t")
colnames(hTFtarget) <- c("TF", "target", "tissue")
hTFtarget <- hTFtarget[-1, ]
hTFtarget <- separate_rows(hTFtarget, TF, sep = "/")
hTFtarget <- hTFtarget[grep("breast", hTFtarget$tissue), ]  # keep only the associations experimentally validated on breast tissue
length(unique(hTFtarget$TF)) # 80

#----------------- converting Gene Symbols to ENSEMBL IDs ---------------------#
mart = useMart(biomart="ensembl", dataset="hsapiens_gene_ensembl")
TF.conversion <- getBM(filters="external_gene_name",
                       attributes=c("ensembl_gene_id", "external_gene_name"), 
                       values=unique(hTFtarget$TF),
                       mart=mart)

TF.conversion <- TF.conversion[order(TF.conversion$external_gene_name), ]
rownames(TF.conversion) <- c(1:length(TF.conversion$ensembl_gene_id))

duplicated.TFs <- which(duplicated(TF.conversion$external_gene_name))
duplicated.TFs <- duplicated.TFs-1 # we keep the last one, as this is the main sequence as opposed to alternative sequences (haplotypes/patches)
TF.conversion <- TF.conversion[-duplicated.TFs, ]

missing.TFs <- setdiff(unique(hTFtarget$TF), TF.conversion$external_gene_name) # 0

mart = useMart(biomart="ensembl", dataset="hsapiens_gene_ensembl")
target.conversion <- getBM(filters="external_gene_name",
                           attributes=c("ensembl_gene_id", "external_gene_name"), 
                           values=unique(hTFtarget$target),
                           mart=mart)

target.conversion <- target.conversion[order(target.conversion$external_gene_name), ]
rownames(target.conversion) <- c(1:length(target.conversion$ensembl_gene_id))

duplicated.targets <- which(duplicated(target.conversion$external_gene_name))
duplicated.targets <- duplicated.targets-1 # we keep the last one, as this is the main sequence as opposed to alternative sequences (haplotypes/patches)
target.conversion <- target.conversion[-duplicated.targets, ]

missing.targets <- setdiff(unique(hTFtarget$target), target.conversion$external_gene_name) # 11855

write.table(missing.targets, file = "reports/associations/TF-gene/missing.targets.txt", sep = "\t", quote = FALSE, row.names = FALSE, col.names = FALSE)

missing.targets <- read.table(file = "reports/associations/TF-gene/hgnc-symbol-check-target.csv", sep = ",")
colnames(missing.targets) <- missing.targets[1, ]
missing.targets <- missing.targets[-1, ]

missing.target.conversion <- getBM(filters="external_gene_name",
                                   attributes=c("ensembl_gene_id", "external_gene_name"), 
                                   values=unique(missing.targets$`Approved symbol`),
                                   mart=mart)

missing.target.conversion <- missing.target.conversion[order(missing.target.conversion$external_gene_name), ]
rownames(missing.target.conversion) <- c(1:length(missing.target.conversion$ensembl_gene_id))

duplicated.targets <- which(duplicated(missing.target.conversion$external_gene_name))
duplicated.targets <- duplicated.targets-1 # we keep the last one, as this is the main sequence as opposed to alternative sequences (haplotypes/patches)
missing.target.conversion <- missing.target.conversion[-duplicated.targets, ]

missing.target.conversion$external_gene_name <- missing.targets$Input[match(missing.target.conversion$external_gene_name, missing.targets$`Approved symbol`)] # replace new names by old names, so that it finds them in the associations

target.conversion <- rbind(target.conversion, missing.target.conversion)

#--------- keeping transcription factors and target genes in our data ---------#
rna.expression <- read.table(file = "results/preprocessing/cookingRNASeq/RNA.expression.txt")
TFs <- intersect(TF.conversion$ensembl_gene_id, rna.expression$V1)
TFs <- as.data.frame(TFs)
targets <- intersect(target.conversion$ensembl_gene_id, rna.expression$V1)
targets <- as.data.frame(targets)

TFs$Gene_Symbol <- TF.conversion$external_gene_name[match(TFs$TFs, TF.conversion$ensembl_gene_id)]
targets$Gene_Symbol <- target.conversion$external_gene_name[match(targets$targets, target.conversion$ensembl_gene_id)]

TFs.associations <- intersect(hTFtarget$TF, TFs$Gene_Symbol)
target.associations <- intersect(hTFtarget$target, targets$Gene_Symbol)

TF.target.associations <- subset(hTFtarget, TF %in% TFs.associations)
TF.target.associations <- subset(TF.target.associations, target %in% target.associations)

TF.target.associations <- TF.target.associations[, -3] # delete tissue column

write.table(TF.target.associations, file = "results/associations/TF-gene/TF.associations.gene.symbol.tab", sep = "\t", quote = FALSE, row.names = FALSE, col.names = FALSE)

#-------- saving transcription factor expression and relevant features --------#
TF.target.associations$TF <- TF.conversion$ensembl_gene_id[match(TF.target.associations$TF, TF.conversion$external_gene_name)]
TF.target.associations$target <- target.conversion$ensembl_gene_id[match(TF.target.associations$target, target.conversion$external_gene_name)]

write.table(TF.target.associations, file = "results/associations/TF-gene/TF.associations.tab", sep = "\t", quote = FALSE, row.names = FALSE, col.names = FALSE)

rna.expression <- read.table(file = "results/preprocessing/cookingRNASeq/RNA.expression.txt")

all.TFs.and.targets <- c(TF.target.associations$TF, TF.target.associations$target)
length(unique(all.TFs.and.targets)) # 15152

TF.expression <- subset(rna.expression, V1 %in% TF.target.associations$TF) # matrix with 80 transcription factors

write.table(TF.expression, file = "results/associations/TF-gene/TF.expression.tab", sep = "\t", quote = FALSE, row.names = FALSE, col.names = FALSE)

rna.DEGs <- read.table(file = "results/preprocessing/cookingRNASeq/common.RNA.DEGs.txt")
TF.DEGs <- subset(rna.DEGs, V1 %in% TF.target.associations$TF) # 7 DEGs are transcription factors
write.table(TF.DEGs, file = "results/associations/TF-gene/TF.DEGs.tab", sep = "\t", quote = FALSE, row.names = FALSE, col.names = FALSE)

TF.target.associations <- read.table(file = "results/associations/TF-gene/TF.associations.tab")
TF.DEGs <- read.table(file = "results/associations/TF-gene/TF.DEGs.tab")
relevant.TF.associations <- subset(TF.target.associations, V1 %in% TF.DEGs$V1) # 17941 TF-gene associations are relevant per this criterion
write.table(relevant.TF.associations, file = "results/associations/TF-gene/DEG.TF.associations.tab", sep = "\t", quote = FALSE, row.names = FALSE, col.names = FALSE)