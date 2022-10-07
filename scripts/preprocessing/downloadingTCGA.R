#------------------------- loading required packages --------------------------#

library(SummarizedExperiment)
library(TCGAbiolinks)
library(dplyr)
library(DT)
library(sesame)

#-------------- querying for intersecting samples for all omics ---------------#

query.exp <- GDCquery(
  project = "TCGA-BRCA",
  data.category = "Transcriptome Profiling",
  data.type = "Gene Expression Quantification", 
  workflow.type = "STAR - Counts"
)

query.mirna <- GDCquery(
  project = "TCGA-BRCA",
  data.category = "Transcriptome Profiling",
  experimental.strategy = "miRNA-Seq",
  workflow.type = "BCGSC miRNA Profiling",
  data.type = "miRNA Expression Quantification",
  data.format = "TXT"
)

query.met <- GDCquery(
  project = "TCGA-BRCA",
  data.category = "DNA Methylation",
  legacy = FALSE,
  platform = c("Illumina Human Methylation 450")
)

query.prot <- GDCquery(
  project = "TCGA-BRCA",
  data.category = "Proteome Profiling",
  legacy = FALSE,
)

# we extract substrings of the barcode from 1 to 15 because the sample is 
# designated by TCGA-XX-XXXX-XX (15 characters in total)

common.samples <- intersect(
  substr(getResults(query.mirna, cols = "cases"), 1, 15),
  substr(getResults(query.exp, cols = "cases"), 1, 15))

common.samples <- intersect(common.samples, 
                            substr(getResults(query.met, cols = "cases"), 1, 15))

common.samples.prot <- intersect(common.samples, 
                                 substr(getResults(query.prot, cols = "cases"), 1, 15))

length(common.samples) # 853 samples
length(common.samples.prot) # 647 samples

#---------- downloading intersecting samples for all omics (- prot) -----------#

query.exp <- GDCquery(
  project = "TCGA-BRCA",
  data.category = "Transcriptome Profiling",
  data.type = "Gene Expression Quantification", 
  workflow.type = "STAR - Counts",
  barcode = common.samples,
  sample.type = c("Primary Tumor", "Solid Tissue Normal")
)

query.mirna <- GDCquery(
  project = "TCGA-BRCA",
  data.category = "Transcriptome Profiling",
  experimental.strategy = "miRNA-Seq",
  workflow.type = "BCGSC miRNA Profiling",
  data.type = "miRNA Expression Quantification",
  data.format = "TXT",
  barcode = common.samples,
  sample.type = c("Primary Tumor", "Solid Tissue Normal")
)

query.met <- GDCquery(
  project = "TCGA-BRCA",
  data.category = "DNA Methylation",
  legacy = FALSE,
  platform = c("Illumina Human Methylation 450"),
  barcode = common.samples,
  data.type = "Methylation Beta Value",
  sample.type = c("Primary Tumor", "Solid Tissue Normal")
)

GDCdownload(query = query.exp, method = "api")
dat <- GDCprepare(query = query.exp, save = TRUE, save.filename = "data/raw/RNA-Seq/RNA.rda")

GDCdownload(query = query.mirna, method = "api")
dat <- GDCprepare(query = query.mirna, save = TRUE, save.filename = "data/raw/miRNA-Seq/miRNA.rda")

GDCdownload(query = query.met, method = "api")
dat <- GDCprepare(query = query.met, save = TRUE, save.filename = "data/raw/met/met.rda")

#------------------------------- loading samples ------------------------------#

load("data/raw/RNA-Seq/RNA.rda")

load("data/raw/miRNA-Seq/miRNA.rda")

load("data/raw/met/met.rda")

#---------- downloading intersecting samples for all omics (+ prot) -----------#

query.exp <- GDCquery(
  project = "TCGA-BRCA",
  data.category = "Transcriptome Profiling",
  data.type = "Gene Expression Quantification", 
  workflow.type = "STAR - Counts",
  barcode = common.samples.prot,
  sample.type = c("Primary Tumor", "Solid Tissue Normal")
)

query.mirna <- GDCquery(
  project = "TCGA-BRCA",
  data.category = "Transcriptome Profiling",
  experimental.strategy = "miRNA-Seq",
  workflow.type = "BCGSC miRNA Profiling",
  data.type = "miRNA Expression Quantification",
  data.format = "TXT",
  barcode = common.samples.prot,
  sample.type = c("Primary Tumor", "Solid Tissue Normal")
)

query.met <- GDCquery(
  project = "TCGA-BRCA",
  data.category = "DNA Methylation",
  legacy = FALSE,
  platform = c("Illumina Human Methylation 450"),
  barcode = common.samples.prot,
  data.type = "Methylation Beta Value",
  sample.type = c("Primary Tumor", "Solid Tissue Normal")
)

query.prot <- GDCquery(
  project = "TCGA-BRCA",
  data.category = "Proteome Profiling",
  legacy = FALSE,
  data.type = "Protein Expression Quantification",
  experimental.strategy = "Reverse Phase Protein Array",
  barcode = common.samples.prot,
  sample.type = c("Primary Tumor", "Solid Tissue Normal")
)

GDCdownload(query = query.exp, method = "api")
rna.with.prot <- GDCprepare(query = query.exp, save = TRUE, save.filename = "data/raw/RNA-Seq/RNA.with.prot.rda")

GDCdownload(query = query.mirna, method = "api")
mirna.with.prot <- GDCprepare(query = query.mirna, save = TRUE, save.filename = "data/raw/miRNA-Seq/miRNA.with.prot.rda")

GDCdownload(query = query.met, method = "api")
met.with.prot <- GDCprepare(query = query.met, save = TRUE, save.filename = "data/raw/met/met.with.prot.rda")

GDCdownload(query = query.prot, method = "api")
prot <- GDCprepare(query = query.prot, save = TRUE, save.filename = "data/raw/prot/prot.rda")

#------------------------------- loading samples ------------------------------#

load("data/raw/RNA-Seq/RNA.with.prot.rda")

load("data/raw/miRNA-Seq/miRNA.with.prot.rda")

load("data/raw/met/met.with.prot.rda")

load("data/raw/prot/prot.rda")