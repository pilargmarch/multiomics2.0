From methylation sites to genes
================
Pilar González Marchante

-   <a href="#associating-cpg-sites-to-gene-symbols"
    id="toc-associating-cpg-sites-to-gene-symbols">Associating CpG sites to
    Gene Symbols</a>
-   <a href="#converting-gene-symbols-to-ensembl-ids"
    id="toc-converting-gene-symbols-to-ensembl-ids">Converting Gene Symbols
    to ENSEMBL IDs</a>

# Associating CpG sites to Gene Symbols

Methylation annotation was automatically retrieved from Methylation
Array Gene Annotation File (Gencode v36)
(<https://gdc.cancer.gov/about-data/gdc-data-processing/gdc-reference-files>)
when downloading through `TCGAbiolinks`. All the used methylation arrays
are Illumina Human Methylation 450.

![](met.genes.info.png)

``` r
library(TCGAbiolinks)
library(SummarizedExperiment)

load("data/raw/met/met.rda")

met.genes.info <- as.data.frame(rowRanges(met))
```

# Converting Gene Symbols to ENSEMBL IDs

We will now convert the filtered CpG sites (tested for differential
expression) and significantly methylated sites into ENSEMBL IDs.

``` r
dea.limma.met <- read.table("results/preprocessing/cookingMet/limma.ordered.tsv")

sites.to.convert <- rownames(dea.limma.met)

met.genes.info <- met.genes.info[sites.to.convert, ]
```

Some CpG probes map to several Gene Symbols, so we will have to split
them into several rows.

``` r
library(tidyr)
met.genes.info <- separate_rows(met.genes.info, Gene_Symbol, sep = ";")

met.genes.info <- met.genes.info[-which(met.genes.info$Gene_Symbol == ""), ]
```

``` r
length(unique(met.genes.info$Gene_Symbol)) # 20796 genes to convert

library(biomaRt)
mart = useMart(biomart="ensembl", dataset="hsapiens_gene_ensembl")
gene.conversion <- getBM(filters="hgnc_symbol",
              attributes=c("ensembl_gene_id", "hgnc_symbol"), 
              values=unique(met.genes.info$Gene_Symbol),
              mart=mart)

missing.genes <- setdiff(unique(met.genes.info$Gene_Symbol), gene.conversion$hgnc_symbol) # 3360

write.table(missing.genes, file = "reports/associations/met-gene/missing.genes.to.convert.txt", quote = FALSE, col.names = FALSE, row.names = FALSE)

duplicated.genes <- which(duplicated(gene.conversion$hgnc_symbol))
duplicated.genes <- duplicated.genes-1 # we keep the last one, as this is the main sequence as opposed to alternative sequences (haplotypes/patches)
gene.conversion <- gene.conversion[-duplicated.genes, ]
```

We have two problems to deal with: missing genes and genes with multiple
ENSEMBL IDs. For the missing genes part, we can input them into [this
multisymbol
checker](https://www.genenames.org/tools/multi-symbol-checker/) to get
the latest version of the gene name.

``` r
symbol.check <- read.table(file = "reports/associations/met-gene/hgnc-symbol-check.csv", sep = ",")
symbol.check <- symbol.check[order(symbol.check$V1), ]
which(duplicated(symbol.check$V1))

missing.gene.conversion <- getBM(filters="hgnc_symbol",
              attributes=c("ensembl_gene_id", "hgnc_symbol"), 
              values=symbol.check$V3,
              mart=mart)

missing.genes <- setdiff(unique(symbol.check$V3), missing.gene.conversion$hgnc_symbol) # 117

duplicated.genes <- which(duplicated(missing.gene.conversion$hgnc_symbol))
duplicated.genes <- duplicated.genes-1 # we keep the last one, as this is the main sequence as opposed to alternative sequences (haplotypes/patches)
missing.gene.conversion <- missing.gene.conversion[-duplicated.genes, ]

gene.conversion <- rbind(gene.conversion, missing.gene.conversion)

gene.conversion <- gene.conversion[order(gene.conversion$hgnc_symbol), ]

rownames(gene.conversion) <- c(1:length(gene.conversion$ensembl_gene_id))

duplicated.genes <- which(duplicated(gene.conversion$ensembl_gene_id))
gene.conversion <- gene.conversion[-duplicated.genes, ]

met.genes.info$Gene_Symbol <- gene.conversion$ensembl_gene_id[match(met.genes.info$Gene_Symbol, gene.conversion$hgnc_symbol)]

met.genes.info <- na.omit(met.genes.info)

met.associations <- met.genes.info[, c(6,7)]
```

Since I couldn’t get Paintomics to work with a methylation associations
file, I had to remap expression and relevant features to ENSEMBL IDs,
deleting those CpG islands without an association.

``` r
dea.limma.met <- read.table(file = "results/preprocessing/cookingMet/limma.ordered.tsv")
dea.limma.met <- cbind(rownames(dea.limma.met), dea.limma.met$logFC)
dea.limma.met <- as.data.frame(dea.limma.met)

dea.limma.met$V1 <- met.associations$Gene_Symbol[match(dea.limma.met$V1, met.associations$probeID)]
dea.limma.met <- na.omit(dea.limma.met)

write.table(dea.limma.met, file = "results/associations/met-gene/met.expression.tab", sep = "\t", quote = FALSE, row.names = FALSE, col.names = FALSE)
```

``` r
met.DEGs <- read.table(file = "results/preprocessing/cookingMet/met.DEGs.txt")
met.DEGs$V1 <- met.associations$Gene_Symbol[match(met.DEGs$V1, met.associations$probeID)]
met.DEGs <- na.omit(met.DEGs)
write.table(met.DEGs, file = "results/associations/met-gene/met.DEGs.tab", sep = "\t", quote = FALSE, row.names = FALSE, col.names = FALSE)
```
