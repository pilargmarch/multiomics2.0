From transcription factors to target genes
================
Pilar González Marchante

-   <a href="#finding-transcription-factors-and-their-target-genes"
    id="toc-finding-transcription-factors-and-their-target-genes">Finding
    transcription factors and their target genes</a>
-   <a
    href="#converting-transcription-factors-and-their-target-genes-from-gene-symbol-to-ensembl"
    id="toc-converting-transcription-factors-and-their-target-genes-from-gene-symbol-to-ensembl">Converting
    transcription factors and their target genes from Gene Symbol to
    ENSEMBL</a>
-   <a
    href="#seeing-which-transcription-factors-and-target-genes-are-in-our-gene-expression-data"
    id="toc-seeing-which-transcription-factors-and-target-genes-are-in-our-gene-expression-data">Seeing
    which transcription factors and target genes are in our gene expression
    data</a>
-   <a href="#saving-transcription-factor-expression-and-relevant-features"
    id="toc-saving-transcription-factor-expression-and-relevant-features">Saving
    transcription factor expression and relevant features</a>

# Finding transcription factors and their target genes

The TF-target regulations were retrieved from hTFtarget, which
integrates huge human TF target resources (7190 ChIP-seq samples of 659
TFs and high-confidence binding sites of 699 TFs) and epigenetic
modification information to predict accurate TF–target regulations. The
`TF-Target-information.txt` file was used, which indicates the
transcription factor, the target gene (both in Gene Symbol) and the
tissue where the relationship was discovered. If a TF-target regulation
was found in more than 30% of ChIP-Seq datasets (at least 3 datasets) in
given tissue, we defined the regulation was putative and listed in this
file.

# Converting transcription factors and their target genes from Gene Symbol to ENSEMBL

``` r
hTFtarget <- read.table(file = "reports/associations/TF-gene/TF-Target-information.hTFtarget.txt", sep = "\t")
colnames(hTFtarget) <- c("TF", "target", "tissue")
hTFtarget <- hTFtarget[-1, ]

library(tidyr)
hTFtarget <- separate_rows(hTFtarget, TF, sep = "/")

hTFtarget <- hTFtarget[grep("breast", hTFtarget$tissue), ]  # keep only the associations experimentally validated on breast tissue

length(unique(hTFtarget$TF)) # 80

library(biomaRt)
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
write.table(missing.TFs, file = "reports/associations/TF-gene/missing.TFs.txt", sep = "\t", quote = FALSE, row.names = FALSE, col.names = FALSE)
```

We have two problems to deal with: missing genes and genes with multiple
ENSEMBL IDs. For the missing genes part, we can input them into [this
multisymbol
checker](https://www.genenames.org/tools/multi-symbol-checker/) to get
the latest version of the gene name. For duplicated genes, what we did
was to keep the main sequence (which happens to always be the last
element) and delete the alternative ones (haplotypes/patches).

``` r
# missing.TFs <- read.table(file = "reports/associations/TF-gene/hgnc-symbol-check-TF.csv", sep = ",")
# colnames(missing.TFs) <- missing.TFs[1, ]
# missing.TFs <- missing.TFs[-1, ]
# 
# missing.TF.conversion <- getBM(filters="external_gene_name",
#               attributes=c("ensembl_gene_id", "external_gene_name"), 
#               values=missing.TFs$`Approved symbol`,
#               mart=mart)
# 
# TF.conversion <- rbind(TF.conversion, missing.TF.conversion)
```

Now let’s do the same, but for the target genes.

``` r
length(unique(hTFtarget$target)) # 32852

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
```

# Seeing which transcription factors and target genes are in our gene expression data

Now I need to see which of these ENSEMBL IDs are actually in our
expression matrix, both for TFs and target genes.

``` r
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
```

# Saving transcription factor expression and relevant features

After writing our associations file, we need to get the expression for
these transcription factors and target genes from our gene expression
matrix. For that, we need to convert them to ENSEMBL IDs.

``` r
TF.target.associations$TF <- TF.conversion$ensembl_gene_id[match(TF.target.associations$TF, TF.conversion$external_gene_name)]
TF.target.associations$target <- target.conversion$ensembl_gene_id[match(TF.target.associations$target, target.conversion$external_gene_name)]

write.table(TF.target.associations, file = "results/associations/TF-gene/TF.associations.tab", sep = "\t", quote = FALSE, row.names = FALSE, col.names = FALSE)

rna.expression <- read.table(file = "results/preprocessing/cookingRNASeq/RNA.expression.txt")

all.TFs.and.targets <- c(TF.target.associations$TF, TF.target.associations$target)
length(unique(all.TFs.and.targets)) # 15152

TF.expression <- subset(rna.expression, V1 %in% TF.target.associations$TF) # matrix with 80 transcription factors

write.table(TF.expression, file = "results/associations/TF-gene/TF.expression.tab", sep = "\t", quote = FALSE, row.names = FALSE, col.names = FALSE)
```

Same thing for our relevant features.

``` r
rna.DEGs <- read.table(file = "results/preprocessing/cookingRNASeq/common.RNA.DEGs.txt")
TF.DEGs <- subset(rna.DEGs, V1 %in% TF.target.associations$TF) # 7 DEGs are transcription factors
write.table(TF.DEGs, file = "results/associations/TF-gene/TF.DEGs.tab", sep = "\t", quote = FALSE, row.names = FALSE, col.names = FALSE)
```

As an approximation we’ll choose relevant associations as those where a
relevant TF appears.

``` r
TF.target.associations <- read.table(file = "results/associations/TF-gene/TF.associations.tab")
TF.DEGs <- read.table(file = "results/associations/TF-gene/TF.DEGs.tab")
relevant.TF.associations <- subset(TF.target.associations, V1 %in% TF.DEGs$V1) # 17941 TF-gene associations are relevant per this criterion
write.table(relevant.TF.associations, file = "results/associations/TF-gene/DEG.TF.associations.tab", sep = "\t", quote = FALSE, row.names = FALSE, col.names = FALSE)
```
