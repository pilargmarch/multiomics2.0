# Cooking proteomics data

## Filtering

We started with normalized (median-centered, log2 transformed) protein expression data for 487 peptide targets and 643 samples, of which 30 peptides and 18 samples had all NA values. After removing these, we ended up with a matrix with data for 457 peptides and 625 samples with no missing values.

## Converting peptide targets to gene IDs

Using [MD Anderson's current expanded antibody list](https://www.mdanderson.org/research/research-resources/core-facilities/functional-proteomics-rppa-core/antibody-information-and-protocols.html), peptide targets were converted to their Entrez ID equivalent. Some peptide targets had to be removed, as they were ambiguous, corresponded to multiple IDs or another target already had that ID. This way we went from 457 to 369 Entrez IDs.

## Analyzing differential expression

Differential expression analysis was performed using `limma`. Differentially expressed genes (DEGs) were selected as those with a q.value/p.adj/FDR < 0.05 and a logFoldChange > 1 or < 1, giving a total of [26 proteins](/results/preprocessing/cookingProt/prot.DEGs.txt) (a 7% of the filtered targets and a 5.7% of the raw peptide targets), 21 being upregulated in cancer samples when compared to normal samples and 5 being downregulated.