# Cooking miRNA-Seq data

## Filtering

Low counts miRNAs were filtered out with `NOISeq`'s `filtered.data` function using the first method, which removes those genes that have an average expression per condition less than `0.5` cpm (in our case) and a coefficient of variation higher than `500`. This filtering allowed us to go from 1,881 miRNAs to 491, a 26% of the total.

## Normalizing

Sample-specific GC-content bias was observed using `NOIseq`, so we used `EDASeq` to normalize this effect (as well as the expected sequencing depth bias). 

This package considers two types of effects on gene-level counts: (1) within-lane gene especific effects, related to gene length and GC-content, and (2) effects related to between-lane distributional differences, e.g. sequencing depth. This two types of effects lead to a two-step normalization process: `withinLaneNormalization` accounts for GC-content or gene length biases, while `betweenLaneNormalization` focuses on normalization for sequencing depth.

`withinLaneNormalization` has four approaches to adjust for GC-content or length on each sample: loess robust local regression, global-scaling, using the median or the upper quantile, and full-quantile normalization. Loess regression performs a regression on the data, according to the gene effect of interest: it regresses log-scale gene counts to GC-content or length using the loess robust local regression method. The three approaches that use quantiles are a function of a defined number of equally sized bins. These bins divide the data according to GC content in several stratus: global scaling using the median will scale the data to have the same median for each bin, global scaling using the upper quantile scales the data to have the same upper quantile and full-quantile normalization (default) will take the several bins quantiles and pair them in order to obtain the median for every quantile, forcing each stratum to be the same. This is an approach similar to microarrays where for each lane the distribution of read counts is matched to a reference distribution that is defined according to the median counts of the sorted lane (or sample). `full` method was used in this `withinLaneNormalization`.

`betweenLaneNormalization` adjusts for lane sequencing depth, i.e. by the number of total read counts per lane i. This normalization aims at rendering lane differences, making the samples comparable. The authors of the package postulated three different types of normalization procedures, the same above referred: global-scaling normalization using upper quantile, global-scaling normalization using the median and full quantile normalization (which is the default). The way these normalizations process the quantiles is the same as described for `withinLaneNormalization` but applied to the lanes (or samples) in study. `full` method was also used in this `betweenLaneNormalization`.

GC bias was checked for again with `NOISeq`, where it seem to have been reduced.

## Analyzing differential expression

Differential expression analysis was performed using `DESeq2`, `limma-voom` and `edgeR`. Differentially expressed genes (DEGs) were selected as those with a q.value/p.adj/FDR < 0.05 and a logFoldChange > 1 or < -1. Intersecting DEGs among the three methods were selected, giving a total of [123 miRNAs](/results/preprocessing/cookingmiRNASeq/common.miRNA.DEGs.txt) (a 25% of the filtered miRNAs and a 6.5% of the raw miRNAs), 67 being upregulated in cancer samples when compared to normal samples and 56 being downregulated.

|    DEGs     | DESeq2 | limma-voom | edgeR | Common |
|:-----------:|:------:|:----------:|:-----:|:------:|
| *Activated* |  132   |     79     |  129  |   67   |
| *Repressed* |   83   |     93     |  76   |   56   |
|   *Total*   | *215*  |   *172*    | *205* | *123*  |