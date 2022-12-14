# Cooking RNA-Seq data

## Filtering

Low counts genes were filtered out with `NOISeq`'s `filtered.data` function using the first method, which removes those genes that have an average expression per condition less than `0.5` cpm (in our case) and a coefficient of variation higher than `500`. This filtering allowed us to go from 60,660 genes to 19,350, a 32% of the total.

## Normalizing

Sample-specific GC-content and length biases were observed using `NOIseq`, so we used `cqn` to normalize these effects (as well as the expected sequencing depth bias). This conditional quantile normalization method combines both within and between-lane (or sample) normalization and is based on a Poisson model for read counts. Lane-specific systematic biases, such as GC-content and length effects, are incorporated as smooth functions using natural cubic splines and estimated using robust quantile regression. In order to account for distributional differences between lanes, a full-quantile normalization procedure is adopted. GC-content and length biases were checked for again with `NOISeq`, where they seem to have been reduced.

## Analyzing differential expression

Differential expression analysis was performed using `DESeq2`, `limma` and `edgeR`. Differentially expressed genes (DEGs) were selected as those with a q.value/p.adj/FDR < 0.05 and a logFoldChange > 2 or < -2. Intersecting DEGs among the three methods were selected, giving a total of [1021 genes](/results/preprocessing/cookingRNASeq/common.RNA.DEGs.txt) (a 5.3% of the filtered genes and a 1.7% of the raw genes), 535 being upregulated in cancer samples when compared to normal samples and 486 being downregulated.

|    DEGs     | DESeq2 | limma  | edgeR  | Common |
|:-----------:|:------:|:------:|:------:|:------:|
| *Activated* |  1203   |  612   |  972   |  535   |
| *Repressed* |  635   |  914   |  586   |  486   |
|   *Total*   | *1838* | *1526* | *1558* | *1021*  |