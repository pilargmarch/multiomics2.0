# Cooking methylation data

## Filtering

Any probes that had any NA values were removed, going from 485,577 CpG sites to 364,019 sites. There were no probes with all 0 values to remove.

## Analyzing differential expression

Differential expression analysis was performed using `limma`. Differentially expressed genes (DEGs) were selected as those with a q.value/p.adj/FDR < 0.05 and a logFoldChange > 2 or < -2, giving a total of [13,163 CpG sites](/results/preprocessing/cookingMet/met.DEGs.txt) (a 3.6% of the filtered probes and a 2.7% of the raw probes), 4,400 being upregulated in cancer samples when compared to normal samples and 8,763 being downregulated.