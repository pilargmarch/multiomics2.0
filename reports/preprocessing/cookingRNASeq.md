Cooking RNA-Seq data
================
Pilar González Marchante

-   <a href="#loading" id="toc-loading">Loading</a>
-   <a href="#filtering" id="toc-filtering">Filtering</a>
-   <a href="#exploring" id="toc-exploring">Exploring</a>
    -   <a href="#length-bias" id="toc-length-bias">Length bias</a>
    -   <a href="#gc-content-bias" id="toc-gc-content-bias">GC content bias</a>
    -   <a href="#rna-composition" id="toc-rna-composition">RNA composition</a>
    -   <a href="#pca-exploration" id="toc-pca-exploration">PCA exploration</a>
-   <a href="#normalizing" id="toc-normalizing">Normalizing</a>
    -   <a href="#cqn" id="toc-cqn">cqn</a>
    -   <a href="#edaseq" id="toc-edaseq">EDASeq</a>
-   <a href="#analyzing-differential-expression"
    id="toc-analyzing-differential-expression">Analyzing differential
    expression</a>
    -   <a href="#deseq2" id="toc-deseq2">DESeq2</a>
    -   <a href="#limma" id="toc-limma">limma</a>
    -   <a href="#edger" id="toc-edger">edgeR</a>
    -   <a href="#intersecting-degs" id="toc-intersecting-degs">Intersecting
        DEGs</a>

# Loading

The GDC harmonizes RNA-Seq data by aligning raw RNA reads to the GRCh38
reference genome build and calculating gene expression levels with
standardized protocols. We downloaded data of category
`Transcriptome Profiling` and type `Gene Expression Quantification`
where the experimental strategy is `RNA-Seq` and the workflow type is
`STAR - Counts`; and saved the `RangedSummarizedExperiment` object in an
`.RData` file named `RNA.rda`.

``` r
load("data/raw/RNA-Seq/RNA.rda")
```

There are 3 functions that allow us to access the most important data
present in `rna`: `colData()`, to access the clinical data;
`rowRanges()`, to access information about the genes; and `assay`, for
the raw counts.

``` r
library(TCGAbiolinks)
library(SummarizedExperiment)

rna.raw.counts <- as.data.frame(assay(rna))
rna.genes.info <- as.data.frame(rowRanges(rna))
rna.sample.info <- as.data.frame(colData(rna))
```

As mentioned, if we want to access the clinical data, we can use the
object created with `colData`:

``` r
colnames(rna.sample.info)
```

    ##  [1] "barcode"                                  
    ##  [2] "patient"                                  
    ##  [3] "sample"                                   
    ##  [4] "shortLetterCode"                          
    ##  [5] "definition"                               
    ##  [6] "sample_submitter_id"                      
    ##  [7] "sample_type_id"                           
    ##  [8] "sample_id"                                
    ##  [9] "sample_type"                              
    ## [10] "days_to_collection"                       
    ## [11] "state"                                    
    ## [12] "initial_weight"                           
    ## [13] "pathology_report_uuid"                    
    ## [14] "submitter_id"                             
    ## [15] "oct_embedded"                             
    ## [16] "is_ffpe"                                  
    ## [17] "tissue_type"                              
    ## [18] "synchronous_malignancy"                   
    ## [19] "ajcc_pathologic_stage"                    
    ## [20] "days_to_diagnosis"                        
    ## [21] "treatments"                               
    ## [22] "last_known_disease_status"                
    ## [23] "tissue_or_organ_of_origin"                
    ## [24] "days_to_last_follow_up"                   
    ## [25] "age_at_diagnosis"                         
    ## [26] "primary_diagnosis"                        
    ## [27] "prior_malignancy"                         
    ## [28] "year_of_diagnosis"                        
    ## [29] "prior_treatment"                          
    ## [30] "ajcc_staging_system_edition"              
    ## [31] "ajcc_pathologic_t"                        
    ## [32] "morphology"                               
    ## [33] "ajcc_pathologic_n"                        
    ## [34] "ajcc_pathologic_m"                        
    ## [35] "classification_of_tumor"                  
    ## [36] "diagnosis_id"                             
    ## [37] "icd_10_code"                              
    ## [38] "site_of_resection_or_biopsy"              
    ## [39] "tumor_grade"                              
    ## [40] "progression_or_recurrence"                
    ## [41] "alcohol_history"                          
    ## [42] "exposure_id"                              
    ## [43] "race"                                     
    ## [44] "gender"                                   
    ## [45] "ethnicity"                                
    ## [46] "vital_status"                             
    ## [47] "age_at_index"                             
    ## [48] "days_to_birth"                            
    ## [49] "year_of_birth"                            
    ## [50] "demographic_id"                           
    ## [51] "year_of_death"                            
    ## [52] "days_to_death"                            
    ## [53] "bcr_patient_barcode"                      
    ## [54] "primary_site"                             
    ## [55] "project_id"                               
    ## [56] "disease_type"                             
    ## [57] "name"                                     
    ## [58] "releasable"                               
    ## [59] "released"                                 
    ## [60] "preservation_method"                      
    ## [61] "days_to_sample_procurement"               
    ## [62] "paper_patient"                            
    ## [63] "paper_Tumor.Type"                         
    ## [64] "paper_Included_in_previous_marker_papers" 
    ## [65] "paper_vital_status"                       
    ## [66] "paper_days_to_birth"                      
    ## [67] "paper_days_to_death"                      
    ## [68] "paper_days_to_last_followup"              
    ## [69] "paper_age_at_initial_pathologic_diagnosis"
    ## [70] "paper_pathologic_stage"                   
    ## [71] "paper_Tumor_Grade"                        
    ## [72] "paper_BRCA_Pathology"                     
    ## [73] "paper_BRCA_Subtype_PAM50"                 
    ## [74] "paper_MSI_status"                         
    ## [75] "paper_HPV_Status"                         
    ## [76] "paper_tobacco_smoking_history"            
    ## [77] "paper_CNV.Clusters"                       
    ## [78] "paper_Mutation.Clusters"                  
    ## [79] "paper_DNA.Methylation.Clusters"           
    ## [80] "paper_mRNA.Clusters"                      
    ## [81] "paper_miRNA.Clusters"                     
    ## [82] "paper_lncRNA.Clusters"                    
    ## [83] "paper_Protein.Clusters"                   
    ## [84] "paper_PARADIGM.Clusters"                  
    ## [85] "paper_Pan.Gyn.Clusters"

Let’s look at some potentially interesting clinical variables.

``` r
table(rna.sample.info$vital_status)
```

    ## 
    ## Alive  Dead 
    ##   722   134

``` r
table(rna.sample.info$ajcc_pathologic_stage)
```

    ## 
    ##    Stage I   Stage IA   Stage IB   Stage II  Stage IIA  Stage IIB  Stage III 
    ##         66         66          4          6        272        203          2 
    ## Stage IIIA Stage IIIB Stage IIIC   Stage IV    Stage X 
    ##        138         21         55         12          5

``` r
table(rna.sample.info$days_to_death)
```

    ## 
    ##    0    1  116  160  172  197  227  239  255  266  295  302  320  322  348  365 
    ##    2    2    1    1    1    1    1    1    1    1    1    1    1    1    1    2 
    ##  377  385  446  524  538  558  571  573  577  584  614  616  639  678  723  785 
    ##    1    1    1    1    1    1    1    1    2    1    4    1    1    2    1    2 
    ##  786  792  821  825  860  879  904  912  959  967  976 1004 1009 1032 1034 1048 
    ##    2    2    1    1    1    1    1    1    2    1    2    1    1    1    3    1 
    ## 1072 1093 1104 1152 1174 1272 1275 1286 1324 1365 1388 1430 1439 1468 1508 1642 
    ##    2    1    1    1    2    1    1    2    1    1    2    1    1    1    1    2 
    ## 1649 1673 1688 1694 1759 1781 1793 1812 1900 1927 2097 2127 2192 2273 2361 2373 
    ##    1    1    1    2    2    1    1    1    1    2    1    2    2    2    1    1 
    ## 2417 2469 2520 2534 2551 2636 2712 2798 2854 2866 2965 3063 3126 3262 3461 3462 
    ##    1    1    2    2    1    1    1    2    1    1    2    1    1    1    1    1 
    ## 3472 3669 3959 6593 
    ##    2    2    2    1

One question we might ask ourselves is which was the tissue type that
was measured: primary tumor or solid tissue.

``` r
summary(factor(rna$sample_type))
```

    ##       Primary Tumor Solid Tissue Normal 
    ##                 781                  76

There are **76 controls** (`Solid Tissue Normal`) (note: these controls
are not healthy individuals, but normal tissue coming from those same
cancer patients) and **781 cancer samples** (`Primary Tumor`).

We’ll delete those features that are constant, redundant or that have no
clinical relevance.

``` r
names_toremove <- c("barcode", "patient", "sample", "sample_submitter_id", "sample_id", "sample_type_id", "state", "pathology_report_uuid", "submitter_id", "oct_embedded", "is_ffpe", "tissue_type", "synchronous_malignancy", "treatments", "last_known_disease_status", "tissue_or_organ_of_origin", "ajcc_staging_system_edition", "classification_of_tumor", "diagnosis_id", "site_of_resection_or_biopsy", "tumor_grade", "progression_or_recurrence", "alcohol_history", "exposure_id", "demographic_id", "bcr_patient_barcode", "primary_site", "project_id", "disease_type", "name", "releasable", "released", "preservation_method", "days_to_sample_procurement", "paper_patient", "paper_Tumor.Type", "paper_Included_in_previous_marker_papers", "paper_vital_status", "paper_days_to_birth", "paper_days_to_death", "paper_days_to_last_followup", "paper_age_at_initial_pathologic_diagnosis", "paper_Tumor_Grade", "paper_MSI_status", "paper_HPV_Status", "paper_tobacco_smoking_history", "paper_CNV Clusters", "paper_Mutation Clusters", "paper_DNA.Methylation Clusters", "paper_mRNA Clusters", "paper_miRNA Clusters", "paper_lncRNA Clusters", "paper_Protein Clusters", "paper_PARADIGM Clusters", "paper_Pan-Gyn Clusters")
names_toremain <- names(colData(rna))
names_toremain <- setdiff(names_toremain, names_toremove)
colData(rna) <- colData(rna)[, names_toremain]
rna.sample.info <- as.data.frame(colData(rna))

# load again because we changed the colnames earlier
load("data/raw/RNA-Seq/RNA.rda")
```

We were able to trim down our data from 85 to **30 clinical variables**.

``` r
head(rna.sample.info)
```

We can now load our RNA-Seq count matrix. We have 60,660 genes in rows
and 857 samples in columns.

``` r
head(rna.raw.counts[, 1:5])
```

    ##                    TCGA-E2-A1L7-01A-11R-A144-07 TCGA-E2-A1L7-11A-33R-A144-07
    ## ENSG00000000003.15                         1689                         4209
    ## ENSG00000000005.6                            16                           71
    ## ENSG00000000419.13                         1810                         1611
    ## ENSG00000000457.14                         1098                         1217
    ## ENSG00000000460.17                          715                          346
    ## ENSG00000000938.13                          624                          787
    ##                    TCGA-BH-A28O-01A-11R-A22K-07 TCGA-D8-A1XU-01A-11R-A14M-07
    ## ENSG00000000003.15                         4583                         5605
    ## ENSG00000000005.6                           135                            6
    ## ENSG00000000419.13                         1531                         4901
    ## ENSG00000000457.14                         1445                         1911
    ## ENSG00000000460.17                          298                          595
    ## ENSG00000000938.13                          515                          410
    ##                    TCGA-AC-A8OP-01A-11R-A36F-07
    ## ENSG00000000003.15                          786
    ## ENSG00000000005.6                            88
    ## ENSG00000000419.13                         1494
    ## ENSG00000000457.14                         1052
    ## ENSG00000000460.17                          229
    ## ENSG00000000938.13                          327

``` r
dim(rna.raw.counts)
```

    ## [1] 60660   857

We can also see the gene names associated with the Ensembl IDs in the
count matrix.

``` r
head(rna.genes.info)
```

    ##                    seqnames     start       end  width strand source type score
    ## ENSG00000000003.15     chrX 100627108 100639991  12884      - HAVANA gene    NA
    ## ENSG00000000005.6      chrX 100584936 100599885  14950      + HAVANA gene    NA
    ## ENSG00000000419.13    chr20  50934867  50958555  23689      - HAVANA gene    NA
    ## ENSG00000000457.14     chr1 169849631 169894267  44637      - HAVANA gene    NA
    ## ENSG00000000460.17     chr1 169662007 169854080 192074      + HAVANA gene    NA
    ## ENSG00000000938.13     chr1  27612064  27635185  23122      - HAVANA gene    NA
    ##                    phase            gene_id      gene_type gene_name level
    ## ENSG00000000003.15    NA ENSG00000000003.15 protein_coding    TSPAN6     2
    ## ENSG00000000005.6     NA  ENSG00000000005.6 protein_coding      TNMD     2
    ## ENSG00000000419.13    NA ENSG00000000419.13 protein_coding      DPM1     2
    ## ENSG00000000457.14    NA ENSG00000000457.14 protein_coding     SCYL3     2
    ## ENSG00000000460.17    NA ENSG00000000460.17 protein_coding  C1orf112     2
    ## ENSG00000000938.13    NA ENSG00000000938.13 protein_coding       FGR     2
    ##                       hgnc_id          havana_gene
    ## ENSG00000000003.15 HGNC:11858 OTTHUMG00000022002.2
    ## ENSG00000000005.6  HGNC:17757 OTTHUMG00000022001.2
    ## ENSG00000000419.13  HGNC:3005 OTTHUMG00000032742.2
    ## ENSG00000000457.14 HGNC:19285 OTTHUMG00000035941.6
    ## ENSG00000000460.17 HGNC:25565 OTTHUMG00000035821.9
    ## ENSG00000000938.13  HGNC:3697 OTTHUMG00000003516.3

RNA-Seq reads were aligned to the genome with STAR (Spliced Transcripts
Alignment to a Reference), a fast RNA-Seq read mapper with support for
splice-junction and fusion read detection. For <a
href="https://docs.gdc.cancer.gov/Data/Bioinformatics_Pipelines/Expression_mRNA_Pipeline/"
class="url">more information on the pipeline used for mRNA counts
generation</a>.

Note: our `rna.raw.counts` matrix contains **unnormalized, unstranded
raw counts**. The `rna` object has these raw `unstranded` counts, as
well as `stranded_first` and `stranded_second`, and `tpm_unstranded`,
`fpkm_unstranded` and `fpkm_uq_unstranded` normalized counts.

![](images/cookingRNASeq/STARcounts.png)

Before performing differential expression analysis, we should preprocess
the raw counts by filtering out low expression counts (and checking for
experimental bias, like batch effect). We also need to perform
normalization, as it’s likely that counts will vary depending on their
length and/or GC content, removing noise inherent to the experimental
technique.

# Filtering

We already have our expression data in `rna.raw.counts`, but we still
have to define our factors as `condition`, `tss`, `plate`, `portion`and
`sample`, which we can easily extract with the `TCGAbiolinks` function
`get_IDs`.

``` r
library(NOISeq)
library(TCGAbiolinks)
```

``` r
barcodes <- get_IDs(rna)
myfactors <- data.frame(barcodes$tss, barcodes$portion, barcodes$plate, barcodes$condition)
head(myfactors)
```

    ##   barcodes.tss barcodes.portion barcodes.plate barcodes.condition
    ## 1           E2              11R           A144             cancer
    ## 2           E2              33R           A144             normal
    ## 3           BH              11R           A22K             cancer
    ## 4           D8              11R           A14M             cancer
    ## 5           AC              11R           A36F             cancer
    ## 6           A2              11R           A084             cancer

And we’ll need additional biological information, such as feature
length, GC content, gene type and chromosome number, which we can get
through the Ensembl BioMart online interface (I tried to get them with
`getGeneLengthAndGCContent` from `EDASeq`, but it kept timing out). I
had to get rid of version numbers from the Ensembl gene IDs, since
otherwise I’d get information about fewer genes.

``` r
# the list used as input for BioMart
ids <- sub('\\.[0-9]*$', '', rna.genes.info$gene_id) # need to remove final digits after the dot (version numbers)
write.table(ids, file = "results/preprocessing/cookingRNASeq/genes.IDs.csv", row.names = FALSE, col.names = FALSE, quote = FALSE)
```

As a rough estimate of gene length, we can use the overlap of all exons
given by `rna.genes.info` as `width`. The same object gives us
`gene_type` and `seq.names` (chromosome number), so we only need to
download GC content from BioMart.

``` r
GCcontent <- read.csv("results/preprocessing/cookingRNASeq/genes.biomart.txt", sep = "\t")
length(GCcontent$Gene...GC.content) # we had 60660 genes but only have GC content info for 60513 of them, thus losing 147 genes
```

    ## [1] 60513

``` r
colnames(GCcontent) <- c("gene_id", "gc_content") # rename headers so we can merge

# problem: we have different versions so we need to merge by ID only
# so we create a new column in rna.genes.info with the gene ID without version number, and we do the same in GCcontent

rna.genes.info$gene_id_no_version <- sub('\\.[0-9]*$', '', rna.genes.info$gene_id)
GCcontent$gene_id_no_version <- sub('\\.[0-9]*$', '', GCcontent$gene_id)

# complete rna.genes.info with GC content
rna.genes.info <- merge(rna.genes.info, GCcontent, by = "gene_id_no_version") # merging turns gene_id into gene_id.x and the gene_id from GCcontent into gene_id.y, so we need to remove those

rna.genes.info$gene_id <- rna.genes.info$gene_id.x
rna.genes.info$gene_id.x <- NULL
rna.genes.info$gene_id.y <- NULL

mygc = c(rna.genes.info$gc_content)
names(mygc) = rna.genes.info$gene_id

mylength = c(rna.genes.info$width)
names(mylength) = rna.genes.info$gene_id

mybiotypes = c(rna.genes.info$gene_type)
names(mybiotypes) = rna.genes.info$gene_id

mychroms = data.frame(rna.genes.info$seqnames, rna.genes.info$start, rna.genes.info$end)
rownames(mychroms) = rna.genes.info$gene_id
colnames(mychroms) <- c("Chr", "GeneStart", "GeneEnd")

head(mygc)
```

    ## ENSG00000000003.15  ENSG00000000005.6 ENSG00000000419.13 ENSG00000000457.14 
    ##              40.40              40.78              40.20              40.14 
    ## ENSG00000000460.17 ENSG00000000938.13 
    ##              39.22              52.92

``` r
head(mylength)
```

    ## ENSG00000000003.15  ENSG00000000005.6 ENSG00000000419.13 ENSG00000000457.14 
    ##              12884              14950              23689              44637 
    ## ENSG00000000460.17 ENSG00000000938.13 
    ##             192074              23122

``` r
head(mybiotypes)
```

    ## ENSG00000000003.15  ENSG00000000005.6 ENSG00000000419.13 ENSG00000000457.14 
    ##   "protein_coding"   "protein_coding"   "protein_coding"   "protein_coding" 
    ## ENSG00000000460.17 ENSG00000000938.13 
    ##   "protein_coding"   "protein_coding"

``` r
head(mychroms)
```

    ##                      Chr GeneStart   GeneEnd
    ## ENSG00000000003.15  chrX 100627108 100639991
    ## ENSG00000000005.6   chrX 100584936 100599885
    ## ENSG00000000419.13 chr20  50934867  50958555
    ## ENSG00000000457.14  chr1 169849631 169894267
    ## ENSG00000000460.17  chr1 169662007 169854080
    ## ENSG00000000938.13  chr1  27612064  27635185

Once we have created the count data matrix, the data.frame for the
factors and the 4 biological annotation objects, we have to pack all
this information into a `NOISeq` object by using the `readData`
function.

``` r
mydata <- NOISeq::readData(data = rna.raw.counts, factors = myfactors, length = mylength, gc = mygc, biotype = mybiotypes, chromosome = mychroms)

myfirst50data <- NOISeq::readData(data = rna.raw.counts[, 1:50], factors = myfactors[1:50, ], length = mylength[1:50], gc = mygc[1:50], biotype = mybiotypes[1:50], chromosome = mychroms[1:50, ]) # for plots and tests that require a smaller sample size
```

Genes with very low counts in all samples provide little evidence for
differential expression. Often samples have many genes with zero or very
low counts. Testing for differential expression for many genes
simultaneously adds to the multiple testing burden, reducing the power
to detect DE genes. IT IS VERY IMPORTANT to filter out genes that have
all zero counts or very low counts. We filter using CPM values rather
than counts because they account for differences in sequencing depth
between samples.

Excluding features with low counts improves differential expression
results since the noise in the data is reduced. `NOISeq` includes three
methods to filter out these low count features: CPM, WIlcoxon test and
proportion test. We’ll try the first two.

But first, we’ll explore the raw counts a bit to later help us choose a
filtering method.

``` r
boxplot(log10(rna.raw.counts[, 1:50])+1, outline = FALSE, las = 2)
```

![](images/cookingRNASeq/boxplot.raw.png)

``` r
mybiodetection <- dat(mydata, k = 0, type = "biodetection", factor = NULL)
par(mfrow = c(1,2))
explo.plot(mybiodetection, samples = c(1, 2), toplot = "protein_coding", plottype = "comparison")
```

![](images/cookingRNASeq/mybiodetection.png)

``` r
[1] "Percentage of protein_coding biotype in each sample:"
TCGA-E2-A1L7-01A-11R-A144-07 
                     53.5015 
TCGA-E2-A1L7-11A-33R-A144-07 
                     51.4557 
[1] "Confidence interval at 95% for the difference of percentages: TCGA-E2-A1L7-01A-11R-A144-07 - TCGA-E2-A1L7-11A-33R-A144-07"
[1] 1.4823 2.6094
[1] "The percentage of this biotype is significantly DIFFERENT for these two samples (p-value = 1.012e-12 )."
```

How shall we choose a CPM threshold? The sensitivity plot can help us.
We can see that most features have between 0 and 1 CPM for these first
50 samples, so the threshold should be between 0.2 and 1, roughly
speaking.

``` r
mycountsbio = dat(myfirst50data, factor = NULL, type = "countsbio")
explo.plot(mycountsbio, toplot = 1, samples = NULL, plottype = "barplot")
```

![](images/cookingRNASeq/mycountsbio.png)

It can also help to visualize the log2(cpm) as histogram and density
plots. In them we see a bimodal distribution that can be split with a
filtering threshold of \~0.5. As a general rule, a good threshold can be
chosen by identifying the CPM that corresponds to a count of 10. Our
minimum library size is 19 millions, whereas median and mean are 57
millions and maximum is 114 million counts. An acceptable threshold
would be between 10/19 = 0.52 and 10/57 = 0.18, so we will try 0.2, 0.5
and 1 as possible CPM thresholds.

``` r
library(SummarizedExperiment)
library(edgeR)
library(limma)
library(ggplot2)

mean_log_cpm <- aveLogCPM(rna.raw.counts)

filter_threshold <- log2(0.5)

ggplot() + aes(x=mean_log_cpm) +
    geom_histogram(binwidth=0.2) +
    geom_vline(xintercept=filter_threshold) +
    ggtitle("Histogram of logCPM before filtering")
```

![](images/cookingRNASeq/histogram.logCPM.png)

``` r
ggplot() + aes(x=mean_log_cpm) +
    geom_density() +
    geom_vline(xintercept=filter_threshold) +
    ggtitle("Density plot of logCPM before filtering") +
    xlim(-6.1, 13.5)

summary(colSums(rna.raw.counts))
```

![](images/cookingRNASeq/density.logCPM.png)

So let’s try CPM filtering with a `CPM threshold = 0.2, 0.5 and 1` and a
`cv.cutoff = 500`, so that we remove those features with low expression
(but not with low variability). We will also apply Wilcoxon test
filtering and compare the results.

``` r
myfiltCPM02 <- filtered.data(rna.raw.counts, factor = myfactors$barcodes.condition, norm = FALSE, depth = NULL, method = 1, cv.cutoff = 500, cpm = 0.2, p.adj = "fdr") # 22797 features (37%) are to be kept for differential expression analysis with filtering method 1

myfiltCPM05 <- filtered.data(rna.raw.counts, factor = myfactors$barcodes.condition, norm = FALSE, depth = NULL, method = 1, cv.cutoff = 500, cpm = 0.5, p.adj = "fdr") # 19350 features (32%) are to be kept for differential expression analysis with filtering method 1

myfiltCPM1 <- filtered.data(rna.raw.counts, factor = myfactors$barcodes.condition, norm = FALSE, depth = NULL, method = 1, cv.cutoff = 500, cpm = 1, p.adj = "fdr") # 17276 features (28.5%) are to be kept for differential expression analysis with filtering method 1

myfiltWilcoxon <- filtered.data(rna.raw.counts, factor = myfactors$barcodes.condition, norm = FALSE, depth = NULL, method = 2, p.adj = "fdr") # 56401 (93%) features are to be kept for differential expression analysis with filtering method 2
```

``` r
boxplot(log10(myfiltCPM02[, 1:50])+1, outline = FALSE, las = 2)
```

![CPM filtering with threshold =
0.2](images/cookingRNASeq/boxplot.filt.CPM.02.png)

``` r
boxplot(log10(myfiltCPM05[, 1:50])+1, outline = FALSE, las = 2)
```

![CPM filtering with threshold =
0.5](images/cookingRNASeq/boxplot.filt.CPM.05.png)

``` r
boxplot(log10(myfiltCPM1[, 1:50])+1, outline = FALSE, las = 2)
```

![CPM filtering with threshold =
1](images/cookingRNASeq/boxplot.filt.CPM.1.png)

``` r
boxplot(log10(myfiltWilcoxon[, 1:50])+1, outline = FALSE, las = 2)
```

![Wilcoxon filtering](images/cookingRNASeq/boxplot.filt.Wilcoxon.png)

What kind of features are these methods filtering out?

``` r
myCPMdata02 <- NOISeq::readData(data = myfiltCPM02, factors = myfactors, length = mylength, gc = mygc, biotype = mybiotypes, chromosome = mychroms)

myCPMdata05 <- NOISeq::readData(data = myfiltCPM05, factors = myfactors, length = mylength, gc = mygc, biotype = mybiotypes, chromosome = mychroms)

myCPMdata1 <- NOISeq::readData(data = myfiltCPM1, factors = myfactors, length = mylength, gc = mygc, biotype = mybiotypes, chromosome = mychroms)

myWilcoxondata <- NOISeq::readData(data = myfiltWilcoxon, factors = myfactors, length = mylength, gc = mygc, biotype = mybiotypes, chromosome = mychroms)
```

``` r
mybiodetectionCPM02 <- dat(myCPMdata02, k = 0, type = "biodetection", factor = NULL)
par(mfrow = c(1,2))
explo.plot(mybiodetectionCPM02, samples = c(1, 2), toplot = "protein_coding", plottype = "comparison")
```

![CPM filtering with threshold =
0.2](images/cookingRNASeq/mybiodetection.CPM.02.png)

``` r
mybiodetectionCPM05 <- dat(myCPMdata05, k = 0, type = "biodetection", factor = NULL)
par(mfrow = c(1,2))
explo.plot(mybiodetectionCPM05, samples = c(1, 2), toplot = "protein_coding", plottype = "comparison")
```

![CPM filtering with threshold =
0.5](images/cookingRNASeq/mybiodetection.CPM.05.png)

``` r
mybiodetectionCPM1 <- dat(myCPMdata1, k = 0, type = "biodetection", factor = NULL)
par(mfrow = c(1,2))
explo.plot(mybiodetectionCPM1, samples = c(1, 2), toplot = "protein_coding", plottype = "comparison")
```

![CPM filtering with threshold =
1](images/cookingRNASeq/mybiodetection.CPM.1.png)

``` r
mybiodetectionWilcoxon <- dat(myWilcoxondata, k = 0, type = "biodetection", factor = NULL)
par(mfrow = c(1,2))
explo.plot(mybiodetectionWilcoxon, samples = c(1, 2), toplot = "protein_coding", plottype = "comparison")
```

![](images/cookingRNASeq/mybiodetection.Wilcoxon.png)

``` r
sum(mydata@featureData@data$Biotype=="protein_coding", na.rm=TRUE) # 19916 protein coding genes
sum(myCPMdata02@featureData@data$Biotype=="protein_coding", na.rm=TRUE) # 16021 protein coding genes
sum(myCPMdata05@featureData@data$Biotype=="protein_coding", na.rm=TRUE) # 15434 protein coding genes
sum(myCPMdata1@featureData@data$Biotype=="protein_coding", na.rm=TRUE) # 14772 protein coding genes
sum(myWilcoxondata@featureData@data$Biotype=="protein_coding", na.rm=TRUE) # 19391 protein coding genes
```

Method 2 (Wilcoxon test) barely does any filtering at all, so we’ll
stick to method 1 (CPM). Given the amount of protein coding features we
are left with with all the different thresholds, we’ll choose a CPM
threshold of 0.5 and will be left with only 32% of the original
features. This makes sense, as the prefiltered raw counts had a large
number of long non-coding RNAs and pseudogenes (which tend to have low
expression in an organism) that were removed with the CPM filtering,
enriching the counts in protein coding genes (going from \~70% to 85% of
total biotypes).

Let’s prep our data once again, this time with our filtered data!

``` r
# delete unnecessary objects
rm(mybiodetection, mybiodetectionCPM02, mybiodetectionCPM05, mybiodetectionCPM1, mybiodetectionWilcoxon, myCPMdata02, myCPMdata05, myCPMdata1, mydata, myfiltCPM02, myfiltCPM1, myfiltWilcoxon, myfirst50data, myWilcoxondata, mean_log_cpm, filter_threshold)

rna.filt.counts <- myfiltCPM05
rm(myfiltCPM05)

save(rna.filt.counts, file = "data/cooked/RNA-Seq/RNA.filt.rda")
```

``` r
load("data/cooked/RNA-Seq/RNA.filt.rda")

myexpdata.filt <- NOISeq::readData(data = rna.filt.counts, factors = myfactors, length = mylength, gc = mygc, biotype = mybiotypes, chromosome = mychroms)
```

# Exploring

What type of normalization should we use? In order to help us make a
decision, we will try to see if there are any length and/or GC content
biases that have to be corrected for.

## Length bias

This plot describes the relationship between the feature length and the
expression values.

``` r
myexplengthbias.filt = dat(myexpdata.filt, factor = "barcodes.condition", type = "lengthbias")
explo.plot(myexplengthbias.filt, samples = NULL, toplot = "global")
```

![](images/cookingRNASeq/length.bias.filt.png)

Since the p-values are significant (even though R<sup>2</sup>
coefficients aren’t higher than 95%) and we can see in the graph that
mean expression varies depending on feature length, we can conclude that
the expression depends on the feature length and as such, a length-based
normalization is required.

## GC content bias

This plot describes the relationship between the feature GC content and
the expression values.

``` r
myexpGCbias.filt = dat(myexpdata.filt, factor = "barcodes.condition", type = "GCbias")
explo.plot(myexpGCbias.filt, samples = NULL, toplot = "global")
```

![](images/cookingRNASeq/gc.bias.filt.png)

Although the p-values are significant, R² values are not high enough for
us to be able to confidently say there exists a GC content bias.
However, we will try a GC-content normalization as well, to try to
remove this potential bias.

## RNA composition

When two samples have different RNA composition, the distribution of
sequencing reads across the features is different in such a way that
although a feature had the same number of read counts in both samples,
it would not mean that it was equally expressed in both.

``` r
myexpcd = dat(myexpdata.filt, type = "cd", norm = FALSE, refColumn = 1)
# "Diagnostic test: FAILED. Normalization is required to correct this bias."
```

Since the test failed (median deviations of the samples with regard to
the reference sample are statistically significant) this means that a
normalization procedure should be used to correct this effect and make
the samples comparable before computing differential expression.

## PCA exploration

Are any of our factors producing some kind of batch effect?

``` r
# PCAs with NOISeq
myexpPCA = dat(myexpdata.filt, type = "PCA")
par(cex = 0.75)
explo.plot(myexpPCA, factor = "barcodes.condition", plottype = "scores")
```

![](images/cookingRNASeq/pca.scores.filt.png)

``` r
explo.plot(myexpPCA, factor = "barcodes.condition", plottype = "loadings")
```

![](images/cookingRNASeq/pca.loadings.filt.png)

``` r
explo.plot(myexpPCA, factor = "barcodes.tss")
```

![](images/cookingRNASeq/pca.tss.png)

``` r
explo.plot(myexpPCA, factor = "barcodes.portion")
```

![](images/cookingRNASeq/pca.portion.png)

``` r
explo.plot(myexpPCA, factor = "barcodes.plate")
```

![](images/cookingRNASeq/pca.plate.png)

We can appreciate in these plots that **none of the factors (TSS,
portion, plate) seem to be contributing towards batch effect** in our
RNA-Seq data, but we shall repeat these plots once our data is
normalized.

# Normalizing

## cqn

Prepare gene information for filtered data: GC content and length.

``` r
load("data/cooked/RNA-Seq/RNA.filt.rda")
library(cqn)

# GeneLengthAndGCContent <- getGeneLengthAndGCContent(sub('\\.[0-9]*$', '', rownames(rna.filt.counts)), "hsa")

# we need to delete the genes for which we have no information, after trying to download it using organism-based annotation packages from Bioconductor instead of the biomart package

# get missing information
# valores_NA <- GeneLengthAndGCContent[rowSums(is.na(GeneLengthAndGCContent)) > 0, ]
# NA_GeneLengthAndGCContent <- getGeneLengthAndGCContent(sub('\\.[0-9]*$', '', rownames(valores_NA)), "hsa", mode = "org.db")
# assembly: TxDb.Hsapiens.UCSC.hg38.knownGene
# genome assembly: BSgenome.Hsapiens.UCSC.hg38
# notna <- na.omit(NA_GeneLengthAndGCContent)

# replace 5 newly found genes :)
# GeneLengthAndGCContent[match(rownames(notna), rownames(GeneLengthAndGCContent)), ] <- notna
# valores_NA <- GeneLengthAndGCContent[rowSums(is.na(GeneLengthAndGCContent)) > 0, ]

# delete 32 missing genes :(
# GeneLengthAndGCContent <- GeneLengthAndGCContent[-(match(rownames(valores_NA), rownames(GeneLengthAndGCContent))), ]

# save information
# gc.length.rna <- GeneLengthAndGCContent 
# save(GeneLengthAndGCContent, file = "results/preprocessing/cookingRNASeq/GC.length.RNA.rda)

# also need to delete 32 missing genes from filtered counts
# rna.filt.counts <- rna.filt.counts[-match(rownames(valores_NA), sub('\\.[0-9]*$', '', rownames(rna.filt.counts))), ]

# and save those filtered complete counts
# save(rna.filt.counts, file = "data/cooked/RNA-Seq/RNA.filt.rda")

# load("results/preprocessing/cookingRNASeq/GC.length.RNA.rda")
```

Run normalization function. `cqn` requires an input of gene length, GC
content and the estimated library size per sample (which it will
estimate as the total sum of the counts if not provided by the user).

The hand off between the two packages is to use `DESeq2` with the
original counts and to supply the offset matrix calculated by `cqn` as a
`normalizationFactor` for the `dds` object.

``` r
load("data/cooked/RNA-Seq/RNA.filt.rda")
library(cqn)

sizeFactors.rna <- colSums(rna.filt.counts)

load("results/preprocessing/cookingRNASeq/GC.length.RNA.rda")
gc.length.rna <- as.data.frame(gc.length.rna)

rna.cqn.norm <- cqn(rna.filt.counts, lengths = gc.length.rna$length, x = gc.length.rna$gc, sizeFactors = sizeFactors.rna, verbose = TRUE)

save(rna.cqn.norm, file = "reports/preprocessing/files/cookingRNASeq/RNA.cqn.norm.rda")

rna.cqn.norm

# Extract the offset, which will be input directly into DEseq2 to normalise the counts
cqnOffset <- rna.cqn.norm$glm.offset
cqnNormFactors <- exp(cqnOffset)
save(cqnNormFactors, file = "reports/preprocessing/files/cookingRNASeq/RNA.cqn.normFactors.rda")

# Extract normalized data to check for bias on NOISeq
rna.cqn.norm.expression <- rna.cqn.norm$y + rna.cqn.norm$offset
rna.cqn.norm.expression <- as.data.frame(rna.cqn.norm.expression)
```

Did `cqn` normalization reduce our gene length and GC content biases?

``` r
load("reports/preprocessing/files/cookingRNASeq/RNA.cqn.norm.rda")
rna.cqn.norm.expression <- rna.cqn.norm$y + rna.cqn.norm$offset
rna.cqn.norm.expression <- as.data.frame(rna.cqn.norm.expression)

# need to prep gene information format
load("results/preprocessing/cookingRNASeq/GC.length.RNA.rda")
gc.length.rna <- as.data.frame(gc.length.rna)
mylength.norm <- as.integer(c(gc.length.rna$length))
names(mylength.norm) <- rownames(gc.length.rna)

mygc.norm <- c(gc.length.rna$gc)
names(mygc.norm) <- rownames(gc.length.rna)
  
mybiotypes.norm <- mybiotypes[match(rownames(gc.length.rna), sub('\\.[0-9]*$', '', names(mybiotypes)))]
mybiotypes.norm <- as.vector(mybiotypes.norm)
  
mychroms.norm <- mychroms[match(rownames(gc.length.rna), sub('\\.[0-9]*$', '', rownames(mychroms))), ]
rownames(mychroms.norm) <- sub('\\.[0-9]*$', '', rownames(mychroms.norm))

barcodes <- get_IDs(rna)
barcodes$condition <- as.factor(barcodes$condition)
barcodes$condition <- relevel(barcodes$condition, ref = "normal")
myfactors <- data.frame(barcodes$tss, barcodes$portion, barcodes$plate, barcodes$condition)

# need to drop the version number of the ENSEMBL IDs because there aren't any in the GC and length information
rownames(rna.cqn.norm.expression) <- sub('\\.[0-9]*$', '', rownames(rna.cqn.norm.expression))

myexpdata.norm <- NOISeq::readData(data = rna.cqn.norm.expression, factors = myfactors, length = mylength.norm, gc = mygc.norm, biotype = mybiotypes.norm, chromosome = mychroms.norm)
```

``` r
myexplengthbias.norm = dat(myexpdata.norm, factor = "barcodes.condition", type = "lengthbias")
explo.plot(myexplengthbias.norm, samples = NULL, toplot = "global")
```

![](images/cookingRNASeq/length.bias.cqn.norm.png)

``` r
myexpGCbias.norm = dat(myexpdata.norm, factor = "barcodes.condition", type = "GCbias")
explo.plot(myexpGCbias.norm, samples = NULL, toplot = "global")
```

![](images/cookingRNASeq/gc.bias.cqn.norm.png)

``` r
myexpPCA.norm = dat(myexpdata.norm, type = "PCA", norm = TRUE, logtransf = TRUE)
par(cex = 0.75)
explo.plot(myexpPCA.norm, factor = "barcodes.condition", plottype = "scores")
```

![](images/cookingRNASeq/pca.scores.cqn.norm.png)

``` r
explo.plot(myexpPCA.norm, factor = "barcodes.condition", plottype = "loadings")
```

![](images/cookingRNASeq/pca.loadings.cqn.norm.png)

``` r
boxplot(rna.cqn.norm.expression[, 1:50], outline = FALSE, las = 2)
```

![](images/cookingRNASeq/boxplot.cqn.norm.png)

GC content and length biases have significantly improved, to the point
where it’s almost gone. What about PCAs? Is there such good separation
between tumor and control samples still?

We can also assess the effect of normalization with some in-built `cqn`
functions.

``` r
library(ggplot2)
# we can compare the fold changes of this normalized data and standard RPKM
# first we compute standard RPKM on a log2 scale
RPM <- sweep(log2(rna.filt.counts + 1), 2, log2(sizeFactors.rna/10^6))
RPKM.std <- sweep(RPM, 1, log2(gc.length.rna$length / 10^3))

whGenes <- which(rowMeans(RPKM.std) >= 2 & gc.length.rna$length >= 100)
M.std <- rowMeans(RPKM.std[whGenes, which(barcodes$condition == "cancer")]) - rowMeans(RPKM.std[whGenes, which(barcodes$condition == "normal")])
A.std <- rowMeans(RPKM.std[whGenes,])
M.cqn <- rowMeans(rna.cqn.norm.expression[whGenes, which(barcodes$condition == "cancer")]) - rowMeans(rna.cqn.norm.expression[whGenes, which(barcodes$condition == "normal")])
A.cqn <- rowMeans(rna.cqn.norm.expression[whGenes,])

par(mfrow = c(1,2))
plot(A.std, M.std, cex = 0.5, pch = 16, xlab = "A", ylab = "M", 
     main = "Standard RPKM", ylim = c(-4,4), xlim = c(0,12), 
     col = alpha("black", 0.25))
plot(A.cqn, M.cqn, cex = 0.5, pch = 16, xlab = "A", ylab = "M", 
     main = "CQN normalized RPKM", ylim = c(-4,4), xlim = c(0,12), 
     col = alpha("black", 0.25))
```

![](images/cookingRNASeq/rpkm.cqn.norm.png)

``` r
# We can plot the effect of GC and length
par(mfrow=c(1,2))
cqnplot(rna.cqn.norm, n = 1, xlab = "GC content", lty = 1, ylim = c(1,7))
cqnplot(rna.cqn.norm, n = 2, xlab = "length", lty = 1, ylim = c(1,7))
```

![](images/cookingRNASeq/cqnplot.png)

## EDASeq

We’ll also try normalizing with `EDASeq`. Following (Risso et al. 2011),
we consider two main types of effects on gene-level counts: (1)
within-lane gene-specific (and possibly lane-specific) effects, e.g.,
related to gene length or GC-content, and (2) effects related to
between-lane distributional differences, e.g., sequencing depth.
Accordingly, `withinLaneNormalization` and `betweenLaneNormalization`
adjust for the first and second type of effects, respectively. We
recommend to normalize for within-lane effects prior to between-lane
normalization.

The `EDASeq` package provides the `SeqExpressionSet` class to store gene
counts, (lane-level) information on the sequenced libraries, and
(gene-level) feature information. We use the data frame met created in
Section secRead for the lane-level data. As for the feature data, we use
gene length and GC-content.

Since `EDASeq` can’t normalize for both GC content and length in one go,
we’ll try several configurations: (1) full GC, (2) full GC then length,
and (3) full length then GC.

``` r
library(EDASeq)

# need to drop the version number of the ENSEMBL IDs because there aren't any in the GC and length information
rownames(rna.filt.counts) <- sub('\\.[0-9]*$', '', rownames(rna.filt.counts))

feature <- data.frame(gc=mygc.norm,length=mylength.norm)
data <- newSeqExpressionSet(counts=as.matrix(rna.filt.counts),
                            featureData=feature,
                            phenoData=data.frame(
                              conditions=barcodes$condition,
                              row.names=barcodes$barcode))

rna.eda.norm.gc <- withinLaneNormalization(data,"gc", which="full")
rna.eda.norm.gc <- betweenLaneNormalization(rna.eda.norm.gc, which="full")

rna.eda.norm.gc.length <- withinLaneNormalization(data, "gc", which="full")
rna.eda.norm.gc.length <- withinLaneNormalization(rna.eda.norm.gc.length, "length", which="full")
rna.eda.norm.gc.length <- betweenLaneNormalization(rna.eda.norm.gc.length, which="full")

rna.eda.norm.length.gc <- withinLaneNormalization(data, "length", which="full")
rna.eda.norm.length.gc <- withinLaneNormalization(rna.eda.norm.length.gc, "gc", which="full")
rna.eda.norm.length.gc <- betweenLaneNormalization(rna.eda.norm.length.gc, which="full")

save(rna.eda.norm.gc, file = "reports/preprocessing/files/cookingRNASeq/RNA.eda.GC.norm.rda")

save(rna.eda.norm.gc.length, file = "reports/preprocessing/files/cookingRNASeq/RNA.eda.GC.length.norm.rda")

save(rna.eda.norm.length.gc, file = "reports/preprocessing/files/cookingRNASeq/RNA.eda.length.GC.norm.rda")
```

Did `EDASeq` normalization reduce our gene length and GC content biases?

``` r
load("reports/preprocessing/files/cookingRNASeq/RNA.eda.GC.norm.rda")
load("reports/preprocessing/files/cookingRNASeq/RNA.eda.GC.length.norm.rda")
load("reports/preprocessing/files/cookingRNASeq/RNA.eda.length.GC.norm.rda")

# Extract normalized counts to check for bias on NOISeq
rna.eda.counts.gc <- rna.eda.norm.gc@assayData$normalizedCounts
rna.eda.counts.gc <- as.data.frame(rna.eda.counts.gc)

rna.eda.counts.gc.length <- rna.eda.norm.gc.length@assayData$normalizedCounts
rna.eda.counts.gc.length <- as.data.frame(rna.eda.counts.gc.length)

rna.eda.counts.length.gc <- rna.eda.norm.length.gc@assayData$normalizedCounts
rna.eda.counts.length.gc <- as.data.frame(rna.eda.counts.length.gc)

library(NOISeq)

myexpdata.norm.eda.gc <- NOISeq::readData(data = rna.eda.counts.gc, factors = myfactors, length = mylength.norm, gc = mygc.norm)

myexpdata.norm.eda.gc.length <- NOISeq::readData(data = rna.eda.counts.gc.length, factors = myfactors, length = mylength.norm, gc = mygc.norm)

myexpdata.norm.eda.length.gc <- NOISeq::readData(data = rna.eda.counts.length.gc, factors = myfactors, length = mylength.norm, gc = mygc.norm)
```

``` r
myexplengthbias.norm.eda.gc = dat(myexpdata.norm.eda.gc, factor = "barcodes.condition", type = "lengthbias")
explo.plot(myexplengthbias.norm.eda.gc, samples = NULL, toplot = "global")
```

![](images/cookingRNASeq/length.bias.eda.gc.norm.png)

``` r
myexplengthbias.norm.eda.gc.length = dat(myexpdata.norm.eda.gc.length, factor = "barcodes.condition", type = "lengthbias")
explo.plot(myexplengthbias.norm.eda.gc.length, samples = NULL, toplot = "global")
```

![](images/cookingRNASeq/length.bias.eda.gc.length.norm.png)

``` r
myexplengthbias.norm.eda.length.gc = dat(myexpdata.norm.eda.length.gc, factor = "barcodes.condition", type = "lengthbias")
explo.plot(myexplengthbias.norm.eda.length.gc, samples = NULL, toplot = "global")
```

![](images/cookingRNASeq/length.bias.eda.length.gc.norm.png)

``` r
myexpGCbias.norm.eda.gc = dat(myexpdata.norm.eda.gc, factor = "barcodes.condition", type = "GCbias")
explo.plot(myexpGCbias.norm.eda.gc, samples = NULL, toplot = "global")
```

![](images/cookingRNASeq/gc.bias.eda.gc.norm.png)

``` r
myexpGCbias.norm.eda.gc.length = dat(myexpdata.norm.eda.gc.length, factor = "barcodes.condition", type = "GCbias")
explo.plot(myexpGCbias.norm.eda.gc.length, samples = NULL, toplot = "global")
```

![](images/cookingRNASeq/gc.bias.eda.gc.length.norm.png)

``` r
myexpGCbias.norm.eda.length.gc = dat(myexpdata.norm.eda.length.gc, factor = "barcodes.condition", type = "GCbias")
explo.plot(myexpGCbias.norm.eda.length.gc, samples = NULL, toplot = "global")
```

![](images/cookingRNASeq/gc.bias.eda.length.gc.norm.png)

Length then GC normalization has the same effect as only GC
normalization, leading me to think it’s not applying it properly. The
third option (length then GC normalization) is actually increasing
length bias, so that one is out of the question. The second option (GC
then length normalization) is the most interesting out of the three, but
`cqn` does a better job for GC bias (while being pretty comparable on
length bias), so we’ll normalize our filtered RNA-Seq counts with `cqn`
instead of with `EDASeq`.

What about PCAs? Is there such good separation between tumor and control
samples still?

``` r
myexpPCA.norm.eda.gc = dat(myexpdata.norm.eda.gc, type = "PCA", norm = TRUE, logtransf = FALSE)
par(cex = 0.75)
explo.plot(myexpPCA.norm.eda.gc, factor = "barcodes.condition", plottype = "scores")
```

![](images/cookingRNASeq/pca.scores.eda.gc.norm.png)

``` r
explo.plot(myexpPCA.norm.eda.gc, factor = "barcodes.condition", plottype = "loadings")
```

![](images/cookingRNASeq/pca.loadings.eda.gc.norm.png)

``` r
boxplot(log10(rna.eda.counts.gc.length[, 1:50])+1, outline = FALSE, las = 2)
```

![](images/cookingRNASeq/boxplot.eda.gc.norm.png)

``` r
myexpPCA.norm.eda.gc.length = dat(myexpdata.norm.eda.gc.length, type = "PCA", norm = TRUE, logtransf = FALSE)
par(cex = 0.75)
explo.plot(myexpPCA.norm.eda.gc.length, factor = "barcodes.condition", plottype = "scores")
```

![](images/cookingRNASeq/pca.scores.eda.gc.length.norm.png)

``` r
explo.plot(myexpPCA.norm.eda.gc.length, factor = "barcodes.condition", plottype = "loadings")
```

![](images/cookingRNASeq/pca.loadings.eda.gc.length.norm.png)

``` r
boxplot(log10(rna.eda.counts.length.gc[, 1:50])+1, outline = FALSE, las = 2)
```

![](images/cookingRNASeq/boxplot.eda.gc.length.norm.png)

``` r
myexpPCA.norm.eda.length.gc = dat(myexpdata.norm.eda.length.gc, type = "PCA", norm = TRUE, logtransf = FALSE)
par(cex = 0.75)
explo.plot(myexpPCA.norm.eda.length.gc, factor = "barcodes.condition", plottype = "scores")
```

![](images/cookingRNASeq/pca.scores.eda.length.gc.norm.png)

``` r
explo.plot(myexpPCA.norm.eda.length.gc, factor = "barcodes.condition", plottype = "loadings")
```

![](images/cookingRNASeq/pca.loadings.eda.length.gc.norm.png)

``` r
boxplot(log10(rna.eda.counts.length.gc[, 1:50])+1, outline = FALSE, las = 2)
```

![](images/cookingRNASeq/boxplot.eda.length.gc.norm.png)

Some of the samples are duplicated, so we need to either remove one of
them or take the average.

``` r
load("data/cooked/RNA-Seq/RNA.norm.rda")
rna.norm.expression <- rna.norm$y + rna.norm$offset

rna.patients <- substr(colnames(rna.norm.expression), 1, 15)
which(duplicated(rna.patients)) 
# 187 188 196 197 638 639 652 714 715
```

The duplicated samples are TCGA-A7-A26J-01 (3x), TCGA-A7-A26E-01 (3x),
TCGA-A7-A13E-01 (3x), TCGA-A7-A0DC-01 (2x) and TCGA-A7-A13D-01 (3x). To
help us decide what to do, we’ll first calculate the correlation between
these samples. Because they all end in “-01”, they all are primary solid
tumor samples.

``` r
plot(x = rna.norm.expression[,186],
y = rna.norm.expression[,187],
pch=19,col="grey",xlab=colnames(rna.norm.expression)[186],ylab=colnames(rna.norm.expression)[187],cex=0.5)

text(x=-5,y=10,
labels = paste(c(
"cor = ",
round(100*cor(rna.norm.expression[,186],
rna.norm.expression[,187]),
digits = 2),
"%"), collapse=""))
```

![](images/cookingRNASeq/correlation.186.187.png)

``` r
plot(x = rna.norm.expression[,186],
y = rna.norm.expression[,188],
pch=19,col="grey",xlab=colnames(rna.norm.expression)[186],ylab=colnames(rna.norm.expression)[188],cex=0.5)

text(x=-5,y=10,
labels = paste(c(
"cor = ",
round(100*cor(rna.norm.expression[,186],
rna.norm.expression[,188]),
digits = 2),
"%"), collapse=""))
```

![](images/cookingRNASeq/correlation.186.188.png)

``` r
plot(x = rna.norm.expression[,187],
y = rna.norm.expression[,188],
pch=19,col="grey",xlab=colnames(rna.norm.expression)[187],ylab=colnames(rna.norm.expression)[188],cex=0.5)

text(x=-5,y=10,
labels = paste(c(
"cor = ",
round(100*cor(rna.norm.expression[,187],
rna.norm.expression[,188]),
digits = 2),
"%"), collapse=""))
```

![](images/cookingRNASeq/correlation.187.188.png) [Broad Institute
recommends the
following:](https://gdac.broadinstitute.org/runs/stddata__2014_01_15/samples_report/READ_Replicate_Samples.html).

“In many instances there is more than one aliquot for a given
combination of individual, platform, and data type. However, only one
aliquot may be ingested into Firehose. Therefore, a set of precedence
rules are applied to select the most scientifically advantageous one
among them. Two filters are applied to achieve this aim: an Analyte
Replicate Filter and a Sort Replicate Filter.

-   Analyte Replicate Filter: the following precedence rules are applied
    when the aliquots have differing analytes. For RNA aliquots, T
    analytes are dropped in preference to H and R analytes, since T is
    the inferior extraction protocol. If H and R are encountered, H is
    the chosen analyte. This is somewhat arbitrary and subject to
    change, since it is not clear at present whether H or R is the
    better protocol. If there are multiple aliquots associated with the
    chosen RNA analyte, the aliquot with the later plate number is
    chosen. For DNA aliquots, D analytes (native DNA) are preferred over
    G, W, or X (whole-genome amplified) analytes, unless the G, W, or X
    analyte sample has a higher plate number.

-   Sort Replicate Filter: the following precedence rules are applied
    when the analyte filter still produces more than one sample. The
    sort filter chooses the aliquot with the highest lexicographical
    sort value, to ensure that the barcode with the highest portion
    and/or plate number is selected when all other barcode fields are
    identical.”

Since in the first case we’re dealing with different analytes
(TCGA-A7-A26J-01B and TCGA-A7-A26J-01A) and they both have the same
plate number (A277), which is the highest, we could discard
TCGA-A7-A26J-01A-11R-A169-0 and take the average of
TCGA-A7-A26J-01B-02R-A277-07 and TCGA-A7-A26J-01A-11R-A277-07, which
have a correlation of 86% between them (correlation between seemingly
unrelated samples sits at approximately 80%).

Another valid way to look at it is to discard FFPE samples, that is,
formalin-fixed paraffin-embedded samples, which are not suitable for
molecular analysis because the RNA and DNA are trapped in the nucleic
acid-protein cross linking from the fixation process.

``` r
rna.sample.info["TCGA-A7-A26J-01B-02R-A277-07", "is_ffpe"] # TRUE
rna.sample.info["TCGA-A7-A26J-01A-11R-A277-07", "is_ffpe"] # FALSE
rna.sample.info["TCGA-A7-A26J-01A-11R-A169-07", "is_ffpe"] # FALSE
```

We will use TCGA-A7-A26J-01A-11R-A277-07 and discard the other two.
Let’s move on to TCGA-A7-A26E-01.

``` r
plot(x = rna.norm.expression[,195],
y = rna.norm.expression[,196],
pch=19,col="grey",xlab=colnames(rna.norm.expression)[195],ylab=colnames(rna.norm.expression)[196],cex=0.5)

text(x=-5,y=10,
labels = paste(c(
"cor = ",
round(100*cor(rna.norm.expression[,195],
rna.norm.expression[,196]),
digits = 2),
"%"), collapse=""))
```

![](images/cookingRNASeq/correlation.195.196.png)

``` r
plot(x = rna.norm.expression[,196],
y = rna.norm.expression[,197],
pch=19,col="grey",xlab=colnames(rna.norm.expression)[196],ylab=colnames(rna.norm.expression)[197],cex=0.5)

text(x=-5,y=10,
labels = paste(c(
"cor = ",
round(100*cor(rna.norm.expression[,196],
rna.norm.expression[,197]),
digits = 2),
"%"), collapse=""))
```

![](images/cookingRNASeq/correlation.196.197.png)

``` r
plot(x = rna.norm.expression[,195],
y = rna.norm.expression[,197],
pch=19,col="grey",xlab=colnames(rna.norm.expression)[195],ylab=colnames(rna.norm.expression)[197],cex=0.5)

text(x=-5,y=10,
labels = paste(c(
"cor = ",
round(100*cor(rna.norm.expression[,195],
rna.norm.expression[,197]),
digits = 2),
"%"), collapse=""))
```

![](images/cookingRNASeq/correlation.195.197.png)

``` r
rna.sample.info["TCGA-A7-A26E-01A-11R-A277-07", "is_ffpe"] # FALSE
rna.sample.info["TCGA-A7-A26E-01A-11R-A169-07", "is_ffpe"] # FALSE
rna.sample.info["TCGA-A7-A26E-01B-06R-A277-07", "is_ffpe"] # TRUE
```

We will apply the same criteria and pick TCGA-A7-A26E-01A-11R-A277-07.
Now on to TCGA-A7-A13E-01.

``` r
plot(x = rna.norm.expression[,637],
y = rna.norm.expression[,638],
pch=19,col="grey",xlab=colnames(rna.norm.expression)[637],ylab=colnames(rna.norm.expression)[638],cex=0.5)

text(x=-5,y=10,
labels = paste(c(
"cor = ",
round(100*cor(rna.norm.expression[,637],
rna.norm.expression[,638]),
digits = 2),
"%"), collapse=""))
```

![](images/cookingRNASeq/correlation.637.638.png)

``` r
plot(x = rna.norm.expression[,637],
y = rna.norm.expression[,639],
pch=19,col="grey",xlab=colnames(rna.norm.expression)[637],ylab=colnames(rna.norm.expression)[639],cex=0.5)

text(x=-5,y=10,
labels = paste(c(
"cor = ",
round(100*cor(rna.norm.expression[,637],
rna.norm.expression[,639]),
digits = 2),
"%"), collapse=""))
```

![](images/cookingRNASeq/correlation.637.639.png)

``` r
plot(x = rna.norm.expression[,638],
y = rna.norm.expression[,639],
pch=19,col="grey",xlab=colnames(rna.norm.expression)[638],ylab=colnames(rna.norm.expression)[639],cex=0.5)

text(x=-5,y=10,
labels = paste(c(
"cor = ",
round(100*cor(rna.norm.expression[,638],
rna.norm.expression[,639]),
digits = 2),
"%"), collapse=""))
```

![](images/cookingRNASeq/correlation.638.639.png)

``` r
rna.sample.info["TCGA-A7-A13E-01A-11R-A12P-07", "is_ffpe"] # FALSE
rna.sample.info["TCGA-A7-A13E-01A-11R-A277-07", "is_ffpe"] # FALSE
rna.sample.info["TCGA-A7-A13E-01B-06R-A277-07", "is_ffpe"] # TRUE
```

We will apply the same criteria and pick TCGA-A7-A13E-01A-11R-A277-07.
Now on to TCGA-A7-A0DC-01.

``` r
plot(x = rna.norm.expression[,651],
y = rna.norm.expression[,652],
pch=19,col="grey",xlab=colnames(rna.norm.expression)[651],ylab=colnames(rna.norm.expression)[652],cex=0.5)

text(x=-5,y=10,
labels = paste(c(
"cor = ",
round(100*cor(rna.norm.expression[,651],
rna.norm.expression[,652]),
digits = 2),
"%"), collapse=""))
```

![](images/cookingRNASeq/correlation.651.652.png)

``` r
rna.sample.info["TCGA-A7-A0DC-01A-11R-A00Z-07", "is_ffpe"] # FALSE
rna.sample.info["TCGA-A7-A0DC-01B-04R-A22O-07", "is_ffpe"] # TRUE
```

Using the same logic we would use TCGA-A7-A0DC-01B-04R-A22O-07, but this
is a FFPE sample. On the other hand, the Broad Institute has blacklisted
TCGA-A7-A0DC-01A-11R-A00Z-07 as Katherine Hoadley identified it as
having a low read count or couldn’t confirm SNP comparison. Since there
isn’t a clear criterion, we’ll take the average of the two samples. Now
on to TCGA-A7-A13D-01.

``` r
plot(x = rna.norm.expression[,713],
y = rna.norm.expression[,714],
pch=19,col="grey",xlab=colnames(rna.norm.expression)[713],ylab=colnames(rna.norm.expression)[714],cex=0.5)

text(x=-5,y=10,
labels = paste(c(
"cor = ",
round(100*cor(rna.norm.expression[,713],
rna.norm.expression[,714]),
digits = 2),
"%"), collapse=""))
```

![](images/cookingRNASeq/correlation.713.714.png)

``` r
plot(x = rna.norm.expression[,713],
y = rna.norm.expression[,715],
pch=19,col="grey",xlab=colnames(rna.norm.expression)[713],ylab=colnames(rna.norm.expression)[715],cex=0.5)

text(x=-5,y=10,
labels = paste(c(
"cor = ",
round(100*cor(rna.norm.expression[,713],
rna.norm.expression[,715]),
digits = 2),
"%"), collapse=""))
```

![](images/cookingRNASeq/correlation.713.715.png)

``` r
plot(x = rna.norm.expression[,714],
y = rna.norm.expression[,715],
pch=19,col="grey",xlab=colnames(rna.norm.expression)[714],ylab=colnames(rna.norm.expression)[715],cex=0.5)

text(x=-5,y=10,
labels = paste(c(
"cor = ",
round(100*cor(rna.norm.expression[,714],
rna.norm.expression[,715]),
digits = 2),
"%"), collapse=""))
```

![](images/cookingRNASeq/correlation.714.715.png)

``` r
rna.sample.info["TCGA-A7-A13D-01A-13R-A277-07", "is_ffpe"] # FALSE
rna.sample.info["TCGA-A7-A13D-01B-04R-A277-07", "is_ffpe"] # TRUE
rna.sample.info["TCGA-A7-A13D-01A-13R-A12P-07", "is_ffpe"] # FALSE
```

Following the same criteria as before, we’ll go with
TCGA-A7-A13D-01A-13R-A277-07. All of these choices (except the second to
last, which doesn’t appear) can be confirmed to be recommended by the
Broad Institute
[here](https://gdac.broadinstitute.org/runs/tmp/sample_report__2018_01_17/Replicate_Samples.html).

Let’s make all these changes.

``` r
load("data/cooked/RNA-Seq/RNA.filt.rda")

rna.filt.counts$`TCGA-A7-A26J-01B-02R-A277-07` <- NULL
rna.filt.counts$`TCGA-A7-A26J-01A-11R-A169-07` <- NULL

rna.filt.counts$`TCGA-A7-A26E-01A-11R-A169-07` <- NULL
rna.filt.counts$`TCGA-A7-A26E-01B-06R-A277-07` <- NULL

rna.filt.counts$`TCGA-A7-A13E-01A-11R-A12P-07` <- NULL
rna.filt.counts$`TCGA-A7-A13E-01B-06R-A277-07` <- NULL

rna.filt.counts$`TCGA-A7-A0DC-01A-11R-A00Z-07` <- sum(rna.filt.counts$`TCGA-A7-A0DC-01A-11R-A00Z-07`, rna.filt.counts$`TCGA-A7-A0DC-01B-04R-A22O-07`)/2
rna.filt.counts$`TCGA-A7-A0DC-01B-04R-A22O-07` <- NULL

rna.filt.counts$`TCGA-A7-A13D-01A-13R-A12P-07` <- NULL
rna.filt.counts$`TCGA-A7-A13D-01B-04R-A277-07` <- NULL

save(rna.filt.counts, file = "data/cooked/RNA-Seq/RNA.filt.rda")
```

We will perform normalization again, just to be sure.

``` r
library(cqn)
load("results/preprocessing/cookingRNASeq/GC.length.RNA.rda")
gc.length.rna <- as.data.frame(gc.length.rna)
sizeFactors.rna <- colSums(rna.filt.counts)

rna.norm <- cqn(rna.filt.counts, lengths = gc.length.rna$length, x = gc.length.rna$gc, sizeFactors = sizeFactors.rna, verbose = TRUE)
# save(rna.norm, file = "data/cooked/RNA-Seq/RNA.norm.rda")

cqnOffset <- rna.norm$glm.offset
cqnNormFactors <- exp(cqnOffset)
# save(cqnNormFactors, file = "data/cooked/RNA-Seq/RNA.normFactors.rda")

rna.norm.expression <- rna.norm$y + rna.norm$offset
rna.norm.expression <- as.data.frame(rna.norm.expression)
```

# Analyzing differential expression

For all, I’ll select DEGs as those genes with a q.value/p.adj/FDR less
than 0.05; and with a logFoldChange bigger than 2. Standard is logFC of
1, but we are being more restrictive because of the large amount of data
we’ll be working with.

## DESeq2

First we create a `DESeqDataSet` object.

``` r
library(DESeq2)
library(TCGAbiolinks)

load("data/cooked/RNA-Seq/RNA.filt.rda")

barcodes <- get_IDs(rna.filt.counts)
myfactors <- data.frame(barcodes$tss, barcodes$portion, barcodes$plate, barcodes$condition)

dds <- DESeqDataSetFromMatrix(countData = rna.filt.counts,
                              colData = myfactors,
                              design = ~ barcodes.condition)

dds$barcodes.condition <- relevel(dds$barcodes.condition, ref = "normal")
```

We generate the normalization factor matrices.

``` r
load("data/cooked/RNA-Seq/RNA.normFactors.rda")
# Before inputing normalizationFactors into DESeq2, you should divide out the geometric mean
normFactors <- normFactors / exp(rowMeans(log(normFactors)))
normalizationFactors(dds) <- normFactors
# This is so that mean(counts(dds, normalized=TRUE)[gene,]) is roughly on the same scale as mean(counts(dds)[gene,])
```

And get the results. Calling results without any arguments will extract
the estimated log2 fold changes and p-values for the last variable in
the design formula. alpha: the significance cutoff used for optimizing
the independent filtering (by default 0.1). If the adjusted p-value
cutoff (FDR) will be a value other than 0.1, alpha should be set to that
value.

``` r
dds <- DESeq(dds)
resultsNames(dds) # lists the coefficients
# [1] "Intercept"                          
# [2] "barcodes.condition_cancer_vs_normal"
res.deseq2 <- results(dds, alpha = 0.05)
summary(res.deseq2)
```

``` r
out of 19318 with nonzero total read count
adjusted p-value < 0.05
LFC > 0 (up)       : 10328, 53%
LFC < 0 (down)     : 5224, 27%
outliers [1]       : 0, 0%
low counts [2]     : 0, 0%
(mean count < 4)
```

We’ll select as significant those genes with a p.adj \< 0.05 and a lFC
\> 2 or lFC \< -2. We get a total of 1203 upregulated DEGs (a 6% of the
filtered genes) in cancer samples, compared to normal ones; and 635
downregulated DEGs (a 3.3% of the filtered genes), from a total of
19,318 filtered genes (originally 60,660 genes in our raw data).

``` r
log.fold.change <- res.deseq2$log2FoldChange
q.value <- res.deseq2$padj
genes.ids <- rownames(rna.filt.counts)
names(log.fold.change) <- genes.ids
names(q.value) <- genes.ids
activated.genes.deseq2 <- genes.ids[log.fold.change > 2 & q.value < 0.05]
activated.genes.deseq2 <- activated.genes.deseq2[!is.na(activated.genes.deseq2)]
repressed.genes.deseq2 <- genes.ids[log.fold.change < - 2 & q.value < 0.05]
repressed.genes.deseq2 <- repressed.genes.deseq2[!is.na(repressed.genes.deseq2)]
length(activated.genes.deseq2) # 1203
length(repressed.genes.deseq2) # 635

log.q.val <- -log10(q.value)
plot(log.fold.change,log.q.val,pch=19,col="grey",cex=0.8,
xlim=c(-8,8),ylim = c(0,240),
xlab="log2(Fold-change)",ylab="-log10(q-value)",cex.lab=1.5)
points(x = log.fold.change[activated.genes.deseq2],
y = log.q.val[activated.genes.deseq2],col="red",cex=0.8,pch=19)
points(x = log.fold.change[repressed.genes.deseq2],
y = log.q.val[repressed.genes.deseq2],col="blue",cex=0.8,pch=19)
```

![](images/cookingRNASeq/volcano.plot.deseq2.png)

Our results table so far only contains the Ensembl gene IDs, but
alternative gene names may be more informative for interpretation.
Bioconductor’s annotation packages help with mapping various ID schemes
to each other. We load the `AnnotationDbi` package and the annotation
package `org.Hs.eg.db`. We can use the *mapIds* function to add
individual columns to our results table. We provide the row names of our
results table as a key, and specify that `keytype=ENSEMBL`. The `column`
argument tells the *mapIds* function which information we want, and the
`multiVals` argument tells the function what to do if there are multiple
possible values for a single input value. Here we ask to just give us
back the first one that occurs in the database. To add the gene symbol
and Entrez ID, we call *mapIds* twice.

``` r
# source(file = "scripts/preprocessing/DEGstoEntrez.R")
# DEGstoEntrez.result <- DEGstoEntrez(res = res.deseq2, activated.genes = activated.genes.deseq2, repressed.genes = repressed.genes.deseq2)
# res.deseq2 <- DEGstoEntrez.result[[1]]
# activated.genes.deseq2 <- DEGstoEntrez.result[[2]]
# repressed.genes.deseq2 <- DEGstoEntrez.result[[3]]
# write.table(activated.genes.deseq2$entrez, file = "results/preprocessing/cookingRNASeq/DESeq2.up.txt", row.names = FALSE, col.names = FALSE, quote = FALSE)
# write.table(repressed.genes.deseq2$entrez, file = "results/preprocessing/cookingRNASeq/DESeq2.down.txt", row.names = FALSE, col.names = FALSE, quote = FALSE)
```

Since some of the tested genes don’t have an Entrez ID associated with
them (they might have been discovered recently or be putative), after
converting we lose some of the DEGs; thus we’ll only use this list when
strictly necessary (if we were to perform GSEA with `ClusterProfiler`,
for example). The resulting DEGs lists will be in the form of ENSEMBL
IDs, as they are the same identifiers that our gene expression matrix
has. We then have 1203 activated genes and 635 supressed ones; if we
were to only keep DEGs with an Entrez ID, that would be only 1047
upregulated and 635 downregulated.

We’ll save both lists of DEGs (upregulated and downregulated) as well as
the DEA results for all of the tested genes, ordered by adjusted
p-value.

``` r
resOrdered <- res.deseq2[order(res.deseq2$padj),]
head(resOrdered)
resOrderedDF <- as.data.frame(resOrdered)
resOrderedDF <- na.omit(resOrderedDF)
write.table(resOrderedDF, file = "results/preprocessing/cookingRNASeq/DESeq2.ordered.csv", row.names=TRUE, col.names=TRUE, sep="\t", quote=FALSE)

write.table(activated.genes.deseq2, file = "results/preprocessing/cookingRNASeq/deseq2.up.txt", quote = FALSE, row.names = FALSE, col.names = FALSE)
write.table(repressed.genes.deseq2, file = "results/preprocessing/cookingRNASeq/deseq2.down.txt", quote = FALSE, row.names = FALSE, col.names = FALSE)
```

## limma

limma-voom transformation could not be applied, since the data is
normalized (and thus has negative values) and there is no way to provide
the function with the normalization factors like `DESeq2` and `edgeR`
have, so `limma` was directly applied to `cqn` normalized expression
data.

I left the design matrix as the intersection (without adding +1 or any
other constants) in order not to have to make the contrasts, thus
simplifying the process.

``` r
library(limma)
load("data/cooked/RNA-Seq/RNA.norm.rda")

rna.norm.expression <- rna.norm$y + rna.norm$offset

library(TCGAbiolinks)
barcodes <- get_IDs(rna.norm.expression)

barcodes$condition <- as.factor(barcodes$condition)
barcodes$condition <- relevel(barcodes$condition, ref = "normal")

design <- model.matrix(~ barcodes$condition)

fit1 <- lmFit(rna.norm.expression, design)

fit2 <- eBayes(fit1)

top.limma <- topTable(fit2, coef = 2, number = Inf)

log.fold.change <- top.limma$logFC
q.value <- top.limma$adj.P.Val
genes.ids <- rownames(top.limma)
names(log.fold.change) <- genes.ids
names(q.value) <- genes.ids

activated.genes.limma <- genes.ids[log.fold.change > 2 & q.value < 0.05]
repressed.genes.limma <- genes.ids[log.fold.change < -2 & q.value < 0.05]

length(activated.genes.limma) # 612
length(repressed.genes.limma) # 914

log.q.val <- -log10(q.value)
plot(log.fold.change,log.q.val,pch=19,col="grey",cex=0.8,
xlim=c(-8,8),ylim = c(0,160),
xlab="log2(Fold-change)",ylab="-log10(q-value)",cex.lab=1.5)
points(x = log.fold.change[activated.genes.limma],
y = log.q.val[activated.genes.limma],col="red",cex=0.8,pch=19)
points(x = log.fold.change[repressed.genes.limma],
y = log.q.val[repressed.genes.limma],col="blue",cex=0.8,pch=19)
```

![](images/cookingRNASeq/volcano.plot.limma.png)

`limma` gives us 612 upregulated genes and 914 downregulated genes. We
save the results.

``` r
# source(file = "scripts/preprocessing/DEGstoEntrez.R")
# DEGstoEntrez.result <- DEGstoEntrez(res = top.limma, activated.genes = activated.genes.limma, repressed.genes = repressed.genes.limma)
# top.limma <- DEGstoEntrez.result[[1]]
# activated.genes.limma <- DEGstoEntrez.result[[2]]
# repressed.genes.limma <- DEGstoEntrez.result[[3]]
# write.table(activated.genes.limma$entrez, file = "results/preprocessing/cookingRNASeq/limma.up.txt", row.names = FALSE, col.names = FALSE, quote = FALSE)
# write.table(repressed.genes.limma$entrez, file = "results/preprocessing/cookingRNASeq/limma.down.txt", row.names = FALSE, col.names = FALSE, quote = FALSE)

topOrdered <- top.limma[order(top.limma$adj.P.Val),]
topOrderedDF <- as.data.frame(topOrdered)
topOrderedDF <- na.omit(topOrderedDF)
write.table(topOrderedDF, file = "results/preprocessing/cookingRNASeq/limma.ordered.csv", row.names=TRUE, col.names=TRUE, sep="\t", quote=FALSE)

write.table(activated.genes.limma, file = "results/preprocessing/cookingRNASeq/limma.up.txt", quote = FALSE, row.names = FALSE, col.names = FALSE)
write.table(repressed.genes.limma, file = "results/preprocessing/cookingRNASeq/limma.down.txt", quote = FALSE, row.names = FALSE, col.names = FALSE)
```

## edgeR

First we construct a `DGEList`.

``` r
library(edgeR)
library(TCGAbiolinks)

load("data/cooked/RNA-Seq/RNA.normFactors.rda")
load("data/cooked/RNA-Seq/RNA.filt.rda")

barcodes <- get_IDs(rna.filt.counts)
barcodes$condition <- as.factor(barcodes$condition)
barcodes$condition <- relevel(barcodes$condition, ref = "normal")

y <- DGEList(counts = rna.filt.counts, lib.size = colSums(rna.filt.counts), group = barcodes$condition, genes = rownames(rna.filt.counts))

y$offset <- normFactors
```

Compute gene-wise exact tests for differences in the means between two
groups of negative-binomially distributed counts.

There are several options: classical (“Compute genewise exact tests for
differences in the means between two groups of negatively-binomially
distributed counts”) and GLM (“GLMs specify probability distributions
according to their mean-variance relationships”).

The classical one can only be used for single-factor designs, such as
ours; while GLMs are usually applied for more complicated experimental
designs, with more than 2 factors (although they can also be used if we
have only one factor, although in this case I have not found it
necessary).

We use classic `edgeR`, as opposed to GLM `edgeR`, as we only have one
factor. The first step is estimating the common and tagwise (gene after
gene) dispersion parameters. We also need to setup a design matrix.

These estimated dispersions can be plotted with a BCV plot, and thus it
can be checked whether the common dispersion really represents the
dispersion between genes.

``` r
design <- model.matrix(~ barcodes$condition)

y <- estimateCommonDisp(y, design = design)
y <- estimateTagwiseDisp(y, design = design)
plotBCV(y)
```

![](images/cookingRNASeq/plotBCV.edger.png)

In this case it is not correct to set the common dispersion for all
genes as there is a lot of difference between the two dispersions and
therefore this cannot be treated as a good representation of the
dispersion of each gene, so gene by gene dispersion was necessary.

After adjusting the dispersion parameters, we have to adjust the model
and perform a significance test. exactTest does the two-by-two
comparisons for the differential expression between the two groups and
topTags takes the output and adjusts the p-values using FDR correction,
and returns all the DEGs (because I have set n=Inf).

``` r
et <- exactTest(y) # performs pair-wise tests for differential expression between two groups
top.edger <- topTags(et, n = Inf) # takes the output from exactTest(), adjusts the raw p-values using the False Discovery Rate (FDR) correction, and returns the top differentially expressed genes

topSig <- top.edger[top.edger$table$FDR < 0.05, ] # we select DEGs with alpha=0.05
dim(topSig)
topSig <- topSig[abs(top.edger$table$logFC) >= 2, ] # we filter the output of dataDEGs by abs(LogFC) >=2
dim(topSig)

# this is equivalent to doing
de <- (decideTestsDGE(et, lfc = 2, p.value = 0.05))
summary(de)
```

``` r
       cancer-normal
Down             586
NotSig         17760
Up               972
```

``` r
detags <- rownames(y)[as.logical(de)]
plotSmear(et, de.tags=detags, main="plotSmear")
abline(h=c(-2,2), col="blue")

activated.genes.edger <- topSig$table$genes[topSig$table$logFC > 0]
length(activated.genes.edger) # 972
repressed.genes.edger <- topSig$table$genes[topSig$table$logFC < 0]
length(repressed.genes.edger) # 586
```

![](images/cookingRNASeq/plotsmear.edger.png)

``` r
top.edger <- top.edger$table
top.edger <- top.edger[order(top.edger$genes), ]
log.fold.change <- top.edger$logFC
q.value <- top.edger$FDR
genes.ids <- rownames(rna.filt.counts)
names(log.fold.change) <- genes.ids
names(q.value) <- genes.ids
activated.genes.edger <- genes.ids[log.fold.change > 2 & q.value < 0.05]
activated.genes.edger <- activated.genes.edger[!is.na(activated.genes.edger)]
repressed.genes.edger <- genes.ids[log.fold.change < - 2 & q.value < 0.05]
repressed.genes.edger <- repressed.genes.edger[!is.na(repressed.genes.edger)]

log.q.val <- -log10(q.value)
plot(log.fold.change,log.q.val,pch=19,col="grey",cex=0.8,
xlim=c(-8,8),ylim = c(0,240),
xlab="log2(Fold-change)",ylab="-log10(q-value)",cex.lab=1.5)
points(x = log.fold.change[activated.genes.edger],
y = log.q.val[activated.genes.edger],col="red",cex=0.8,pch=19)
points(x = log.fold.change[repressed.genes.edger],
y = log.q.val[repressed.genes.edger],col="blue",cex=0.8,pch=19)
```

![](images/cookingRNASeq/volcano.plot.edger.png)

In cancer samples (compared to normal samples) there would be 972
upregulated genes and 586 downregulated genes.

Again, let’s save the results.

``` r
top.edger <- as.data.frame(top.edger)
# source(file = "scripts/preprocessing/DEGstoEntrez.R")
# DEGstoEntrez.result <- DEGstoEntrez(res = top.edger, activated.genes = activated.genes.edger, repressed.genes = repressed.genes.edger)
# top.edger <- DEGstoEntrez.result[[1]]
# activated.genes.edger <- DEGstoEntrez.result[[2]]
# repressed.genes.edger <- DEGstoEntrez.result[[3]]
# write.table(activated.genes.edger$entrez, file = "results/preprocessing/cookingRNASeq/edgeR.up.txt", row.names = FALSE, col.names = FALSE, quote = FALSE)
# write.table(repressed.genes.edger$entrez, file = "results/preprocessing/cookingRNASeq/edgeR.down.txt", row.names = FALSE, col.names = FALSE, quote = FALSE)

topOrdered <- top.edger[order(top.edger$FDR),]
topOrderedDF <- as.data.frame(topOrdered)
topOrderedDF <- na.omit(topOrderedDF)
write.table(topOrderedDF, file = "results/preprocessing/cookingRNASeq/edgeR.ordered.csv", row.names=TRUE, col.names=TRUE, sep="\t", quote=FALSE)

write.table(activated.genes.edger, file = "results/preprocessing/cookingRNASeq/edger.up.txt", quote = FALSE, row.names = FALSE, col.names = FALSE)
write.table(repressed.genes.edger, file = "results/preprocessing/cookingRNASeq/edger.down.txt", quote = FALSE, row.names = FALSE, col.names = FALSE)
```

## Intersecting DEGs

``` r
# activated.genes.deseq2 <- read.table(file = "results/preprocessing/cookingRNASeq/DESeq2.up.txt")
# activated.genes.deseq2 <- as.vector(activated.genes.deseq2$V1)
# 
# repressed.genes.deseq2 <- read.table(file = "results/preprocessing/cookingRNASeq/DESeq2.down.txt")
# repressed.genes.deseq2 <- as.vector(repressed.genes.deseq2$V1)
# 
# activated.genes.limma <- read.table(file = "results/preprocessing/cookingRNASeq/limma.up.txt")
# activated.genes.limma <- as.vector(activated.genes.limma$V1)
# 
# repressed.genes.limma <- read.table(file = "results/preprocessing/cookingRNASeq/limma.down.txt")
# repressed.genes.limma <- as.vector(repressed.genes.limma$V1)
# 
# activated.genes.edger <- read.table(file = "results/preprocessing/cookingRNASeq/edgeR.up.txt")
# activated.genes.edger <- as.vector(activated.genes.edger$V1)
# 
# repressed.genes.edger <- read.table(file = "results/preprocessing/cookingRNASeq/edgeR.down.txt")
# repressed.genes.edger <- as.vector(repressed.genes.edger$V1)

common.activated <- intersect(intersect(activated.genes.deseq2, activated.genes.edger), activated.genes.limma) 
length(common.activated) # 535

common.repressed <- intersect(intersect(repressed.genes.deseq2, repressed.genes.edger), repressed.genes.limma) 
length(common.repressed) # 486

write.table(common.activated, file = "results/preprocessing/cookingRNASeq/common.up.txt", row.names = FALSE, col.names = FALSE, quote = FALSE)

write.table(common.repressed, file = "results/preprocessing/cookingRNASeq/common.down.txt", row.names = FALSE, col.names = FALSE, quote = FALSE)
```

|    DEGs     | DESeq2 | limma  | edgeR  | Common |
|:-----------:|:------:|:------:|:------:|:------:|
| *Activated* |  1203  |  612   |  972   |  535   |
| *Repressed* |  635   |  914   |  586   |  486   |
|   *Total*   | *1838* | *1526* | *1558* | *1021* |

We have 1021 DEGs, a 5.3% of the filtered genes (19,318) and a 1.7% of
the original genes (60,660).
