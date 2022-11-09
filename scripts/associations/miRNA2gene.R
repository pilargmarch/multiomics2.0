#------------------------- loading required packages --------------------------#
library(biomaRt)

miRNA2gene <- read.table("reports/associations/miRNA-gene/all_miRWalk_miRNA_Targets.txt")
colnames(miRNA2gene) <- c("mirbase_id", "genesymbol")
miRNA2gene <- miRNA2gene[-1, ]
miRNA2gene <- unique(miRNA2gene)

#----------------- converting Gene Symbols to ENSEMBL IDs ---------------------#
gene.symbols <- unique(miRNA2gene$genesymbol) # we need to map 2293 genes
mart = useMart(biomart="ensembl", dataset="hsapiens_gene_ensembl")
gene.conversion <- getBM(filters="external_gene_name",
                         attributes=c("ensembl_gene_id", "external_gene_name"), 
                         values= gene.symbols,
                         mart=mart)
gene.conversion <- gene.conversion[order(gene.conversion$external_gene_name), ]
rownames(gene.conversion) <- c(1:length(gene.conversion$ensembl_gene_id))

duplicated.genes <- which(duplicated(gene.conversion$external_gene_name))
duplicated.genes <- duplicated.genes-1 # we keep the last one, as this is the main sequence as opposed to alternative sequences (haplotypes/patches)
gene.conversion <- gene.conversion[-duplicated.genes, ]

missing.genes <- setdiff(gene.symbols, gene.conversion$external_gene_name) # 4
write.table(missing.genes, file = "reports/associations/miRNA-gene/missing.genes.txt", sep = "\t", quote = FALSE, row.names = FALSE, col.names = FALSE)

missing.genes <- read.table(file = "reports/associations/miRNA-gene/hgnc-symbol-check.csv", sep = ",")
colnames(missing.genes) <- missing.genes[1, ]
missing.genes <- missing.genes[-1, ]

missing.gene.conversion <- getBM(filters="external_gene_name",
                                 attributes=c("ensembl_gene_id", "external_gene_name"),
                                 values=missing.genes$`Approved symbol`,
                                 mart=mart)

missing.gene.conversion$external_gene_name <- missing.genes$Input[match(missing.gene.conversion$external_gene_name, missing.genes$`Approved symbol`)] # replace new names by old names, so that it finds them in the associations

gene.conversion <- rbind(gene.conversion, missing.gene.conversion)

miRNA2gene$ensembl_gene_id <- gene.conversion$ensembl_gene_id[match(miRNA2gene$genesymbol, gene.conversion$external_gene_name)]

miRNA2gene[[1]] <- tolower(miRNA2gene[[1]]) # turn from uppercase to lowercase so that it can match DEGs, which are lowercase

write.table(miRNA2gene[, -2], file = "results/associations/miRNA-gene/all_miRWalk_miRNA_Targets.txt", sep = "\t", row.names = FALSE, quote = FALSE)

miRNA.associations <- read.table(file = "results/associations/miRNA-gene/all_miRWalk_miRNA_Targets.txt", header = TRUE)
miRNA.DEGs <- read.table(file = "results/preprocessing/cookingmiRNASeq/common.miRNA.DEGs.txt")
relevant.miRNA.associations <- subset(miRNA.associations, mirbase_id %in% miRNA.DEGs$V1) # 169 miRNA-gene associations are relevant per this criterion
write.table(relevant.miRNA.associations, file = "results/associations/miRNA-gene/DEG_miRWalk_miRNA_Targets.txt", sep = "\t", quote = FALSE, row.names = FALSE, col.names = FALSE)