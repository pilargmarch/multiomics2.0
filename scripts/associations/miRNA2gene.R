#------------------------- loading required packages --------------------------#
library(biomaRt)

miRNA2gene <- read.table("reports/associations/miRNA-gene/all_miRWalk_miRNA_Targets.txt")
colnames(miRNA2gene) <- c("mirbase_id", "genesymbol")
miRNA2gene <- miRNA2gene[-1, ]
miRNA2gene <- unique(miRNA2gene)

#----------------- converting Gene Symbols to ENSEMBL IDs ---------------------#
gene.symbols <- unique(miRNA2gene$genesymbol) # we need to map 231 genes
mart = useMart(biomart="ensembl", dataset="hsapiens_gene_ensembl")
gene.conversion <- getBM(filters="external_gene_name",
                         attributes=c("ensembl_gene_id", "external_gene_name"), 
                         values= gene.symbols,
                         mart=mart)
gene.conversion <- gene.conversion[order(gene.conversion$external_gene_name), ]
gene.conversion$external_gene_name[which(duplicated(gene.conversion$external_gene_name))]

miRNA2gene <- miRNA2gene[order(miRNA2gene$genesymbol), ]
gene.symbols <- unique(miRNA2gene$genesymbol)

unique.target.genes <- unique(gene.conversion$external_gene_name)

setdiff(gene.symbols, unique.target.genes) # we don't have info for "C5orf51" but this is because it's a synonym for RIMOC1, so we will have to manually add ENSG00000205765

gene.conversion <- rbind(gene.conversion, c("ENSG00000205765", "C5orf51"))

gene.conversion <- gene.conversion[order(gene.conversion$external_gene_name), ]
rownames(gene.conversion) <- c(1:length(gene.conversion$ensembl_gene_id))
unique.target.genes <- unique(gene.conversion$external_gene_name)

gene.conversion.no.duplicates <- gene.conversion[-c(179, 168, 205, 7, 88, 132, 228, 190), ]
gene.conversion.no.duplicates$external_gene_name == gene.symbols

rownames(gene.conversion.no.duplicates) <- gene.conversion.no.duplicates$external_gene_name
miRNA2gene$ENSEMBL <- gene.conversion.no.duplicates[miRNA2gene$genesymbol, ]$ensembl_gene_id
miRNA2gene$genesymbol <- NULL

# write.table(miRNA2gene, file = "results/associations/miRNA-gene/all_miRWalk_miRNA_Targets.txt", sep = "\t", row.names = FALSE, quote = FALSE)