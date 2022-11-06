From miRNAs to genes
================
Pilar González Marchante

-   <a href="#associating-mirbase-ids-to-gene-symbols-with-mirwalk"
    id="toc-associating-mirbase-ids-to-gene-symbols-with-mirwalk">Associating
    miRbase IDs to Gene Symbols with miRWalk</a>
-   <a href="#converting-gene-symbols-to-ensembl-ids"
    id="toc-converting-gene-symbols-to-ensembl-ids">Converting Gene Symbols
    to ENSEMBL IDs</a>

# Associating miRbase IDs to Gene Symbols with miRWalk

We discovered which miRNAs regulate which genes using online tool
miRWalk, which allows us to retrieve only those associations that have
been experimentally validated (not just predicted *in silico*) and that
appear on three widely used, independent databases: TargetScan, miRDB
and miRTarBase, using a combination of experimentation and prediction.

![](miRWalk.png)

``` r
miRNA2gene <- read.table("reports/associations/miRNA-gene/all_miRWalk_miRNA_Targets.txt")
colnames(miRNA2gene) <- c("mirbase_id", "genesymbol")
miRNA2gene <- miRNA2gene[-1, ]
miRNA2gene <- unique(miRNA2gene)
```

``` r
head(miRNA2gene)

   mirbase_id genesymbol
2 hsa-mir-326       FGF1
3 hsa-mir-326       FGF1
4 hsa-mir-326       FGF1
5 hsa-mir-326       FGF1
6 hsa-mir-326       FGF1
7 hsa-mir-326       FGF1
```

# Converting Gene Symbols to ENSEMBL IDs

However, there is one problem: our gene expression matrix comes in the
form of ENSEMBL IDs, yet this website gives genes as Gene Symbols. We
have one more task to complete: to convert Gene Symbols to ENSEMBL IDs,
with the help of `biomaRt`.

``` r
gene.symbols <- unique(miRNA2gene$genesymbol) # we need to map 231 genes
library(biomaRt)
mart = useMart(biomart="ensembl", dataset="hsapiens_gene_ensembl")
gene.conversion <- getBM(filters="external_gene_name",
              attributes=c("ensembl_gene_id", "external_gene_name"), 
              values= gene.symbols,
              mart=mart)
```

``` r
gene.conversion <- gene.conversion[order(gene.conversion$external_gene_name), ]
```

``` r
gene.conversion$external_gene_name[which(duplicated(gene.conversion$external_gene_name))]

[1] "SESN2"  "RBFOX2" "SURF4"  "AGPAT5" "HBP1"   "NIPA1"  "YTHDC1" "SMCR8" 
```

As expected, some gene symbols map to several ENSEMBL IDs. We decided to
map everything into ENSEMBL instead of going from ENSEMBL to Entrez/Gene
Symbol because some ENSEMBL transcripts are novel or putative so they
haven’t been given an external ID yet. If we converted from ENSEMBL to
Entrez, we would have to remove genes from the expression matrix; but if
we convert from Entrez/Symbol to ENSEMBL, we just have to add extra
associations to the associations matrix, which is fine.

``` r
miRNA2gene <- miRNA2gene[order(miRNA2gene$genesymbol), ]
gene.symbols <- unique(miRNA2gene$genesymbol)

unique.target.genes <- unique(gene.conversion$external_gene_name)

setdiff(gene.symbols, unique.target.genes) # we don't have info for "C5orf51" but this is because it's a synonym for RIMOC1, so we will have to manually add ENSG00000205765

gene.conversion <- rbind(gene.conversion, c("ENSG00000205765", "C5orf51"))

gene.conversion <- gene.conversion[order(gene.conversion$external_gene_name), ]
rownames(gene.conversion) <- c(1:length(gene.conversion$ensembl_gene_id))
unique.target.genes <- unique(gene.conversion$external_gene_name)
```

``` r
gene.conversion[gene.conversion$external_gene_name=="SESN2", ] # 179 (ENSG00000285069), 180 (ENSG00000130766)
gene.conversion[gene.conversion$external_gene_name=="RBFOX2", ] # 168 (ENSG00000277564), 169 (ENSG00000100320)
gene.conversion[gene.conversion$external_gene_name=="SURF4", ] # 205 (ENSG00000280951), 206 (ENSG00000148248)
gene.conversion[gene.conversion$external_gene_name=="AGPAT5", ] # 7 (ENSG00000284980), 8 (ENSG00000155189)
gene.conversion[gene.conversion$external_gene_name=="HBP1", ] # 88 (ENSG00000283847), 89 (ENSG00000105856)
gene.conversion[gene.conversion$external_gene_name=="NIPA1", ] # 132 (ENSG00000288478), 133 (ENSG00000170113)
gene.conversion[gene.conversion$external_gene_name=="YTHDC1", ] # 228 (ENSG00000275272), 229 (ENSG00000083896)
gene.conversion[gene.conversion$external_gene_name=="SMCR8", ] # 190 (ENSG00000283741), 191 (ENSG00000176994)

gene.conversion.no.duplicates <- gene.conversion[-c(179, 168, 205, 7, 88, 132, 228, 190), ]
gene.conversion.no.duplicates$external_gene_name == gene.symbols

rownames(gene.conversion.no.duplicates) <- gene.conversion.no.duplicates$external_gene_name
miRNA2gene$ENSEMBL <- gene.conversion.no.duplicates[miRNA2gene$genesymbol, ]$ensembl_gene_id
miRNA2gene$genesymbol <- NULL

write.table(miRNA2gene, file = "results/associations/miRNA-gene/all_miRWalk_miRNA_Targets.txt", sep = "\t", row.names = FALSE, quote = FALSE)
```

Then we manually added the 8 genes that had 2 ENSEMBL IDs, duplicating
those associations, leaving us with 269 miRNA-gene associations.

``` r
head(miRNA2gene)

    mirbase_id genesymbol
2  hsa-mir-326       FGF1
8  hsa-mir-107       CDK6
9  hsa-mir-107       CDK8
11 hsa-mir-107      CLOCK
12 hsa-mir-107    CSNK1G3
14 hsa-mir-107     DICER1 
```