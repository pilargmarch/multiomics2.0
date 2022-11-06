#------------------------- loading required packages --------------------------#
library(biomaRt)

dea.limma.prot <- read.table(file = "results/preprocessing/cookingProt/limma.ordered.tsv")


#--------------------- converting Entrez to ENSEMBL IDs -----------------------#
mart = useMart(biomart="ensembl", dataset="hsapiens_gene_ensembl")
gene.conversion <- getBM(filters="entrezgene_id",
                         attributes=c("ensembl_gene_id", "entrezgene_id", "entrezgene_accession"), 
                         values= rownames(dea.limma.prot),
                         mart=mart)

which(duplicated(gene.conversion$entrezgene_id))

gene.conversion <- gene.conversion[order(gene.conversion$entrezgene_id), ]

rownames(gene.conversion) <- gene.conversion$ensembl_gene_id

gene.conversion <- gene.conversion[-c(2, 32, 33, 34, 35, 36, 37, 38, 39, 66, 127, 128, 129, 130, 131, 132, 133, 134, 135, 136, 137, 138, 139, 143, 144, 145, 146, 172, 178, 179, 185, 197, 225, 226, 227, 228, 229, 230, 244, 250, 262, 270, 310, 317, 330, 363, 380, 415), ]

setdiff(rownames(dea.limma.prot), gene.conversion$entrezgene_id) # 780

missing.gene <- data.frame("ENSG00000204580", "780", "DDR1")
names(missing.gene) <- c("ensembl_gene_id", "entrezgene_id", "entrezgene_accession")
rownames(missing.gene) <- "ENSG00000204580"
gene.conversion <- rbind(gene.conversion, missing.gene)

write.table(gene.conversion, file = "results/associations/protein-gene/entrez2ensembl.txt", sep = "\t", quote = FALSE, row.names = FALSE)

#------------------------ replacing IDs in DEA table --------------------------#

dea.limma.prot <- dea.limma.prot[order(as.numeric(row.names(dea.limma.prot))), ]
gene.conversion <- gene.conversion[order(as.numeric(gene.conversion$entrezgene_id)), ]

rownames(dea.limma.prot) == gene.conversion$entrezgene_id

rownames(dea.limma.prot) <- gene.conversion$ensembl_gene_id

prot.expression <- as.data.frame(dea.limma.prot[, 1])
rownames(prot.expression) <- rownames(dea.limma.prot)
names(prot.expression) <- "logFC"

write.table(prot.expression, file = "results/associations/protein-gene/prot.expression.txt", sep = "\t", quote = FALSE, col.names = FALSE)

prot.DEGs <- read.table(file = "results/preprocessing/cookingProt/prot.DEGs.txt")
names(prot.DEGs) <- "entrezgene_id"

prot.DEGs <- merge(prot.DEGs, gene.conversion, by = "entrezgene_id")
prot.DEGs <- prot.DEGs[, 2]

write.table(prot.DEGs, file = "results/associations/protein-gene/prot.DEGs.txt", sep = "\t", quote = FALSE, col.names = FALSE)