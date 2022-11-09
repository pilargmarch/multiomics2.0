# multiomics2.0
![Charting a multi-omic universe. Image credits: Isabel Romero Calvo/EMBL.](reports/cover.jpg)
<sup>[Charting a multi-omic universe. Image credits: Isabel Romero Calvo/EMBL.](https://www.embl.org/news/science/charting-a-multi-omic-universe/)</sup>

This is the repository for my Master's Final Thesis. First, RNA-Seq, miRNA-Seq, methylation and proteomics data was selected and downloaded from TCGA. After the data was processed and tested for differential expression, regulatory associations between these omics were retrieved from various databases. Finally, relevant features, associations and gene expression were all put into different methods to try to come up with a regulatory model for gene expression in breast cancer.

## :dna: [Data](data/)
RNA-Seq, miRNA-Seq, methylation and proteomics data used in this project, both as raw and cooked.

---

## :computer: [Scripts](scripts/)
Final `.R` scripts used to generate the obtained results.

### 🍳 [Preprocessing](scripts/preprocessing/)
- [downloadingTCGA.R](scripts/preprocessing/downloadingTCGA.R)
- [cookingRNASeq.R](scripts/preprocessing/cookingRNASeq.R)
- [cookingmiRNASeq.R](scripts/preprocessing/cookingmiRNASeq.R)
- [cookingMet.R](scripts/preprocessing/cookingMet.R)
- [cookingProt.R](scripts/preprocessing/cookingProt.R)

### 👫 [Associations](scripts/associations/)
- [From transcription factors to genes](scripts/associations/TF2gene.R)
- [From miRNAs to genes](scripts/associations/miRNA2gene.R)
- [From methylation sites to genes](scripts/associations/met2gene.R)
- [From proteins to genes](scripts/associations/protein2gene.R)

### 📈 [Integration](scripts/integration/)
- [MORE]()
- [pls-Multiblock]()

---

## 🗺️ [Reports](reports/)

Exploratory `.Rmd` notebooks (knitted as `.md`) used to decide which pipelines to follow in the final scripts.

### 🍳 [Preprocessing](reports/preprocessing/)
- [downloadingTCGA.md](reports/preprocessing/downloadingTCGA.md)
- [cookingRNASeq.md](reports/preprocessing/cookingRNASeq.md)
- [cookingmiRNASeq.md](reports/preprocessing/cookingmiRNASeq.md)
- [cookingMet.md](reports/preprocessing/cookingMet.md)
- [cookingProt.md](reports/preprocessing/cookingProt.md)

### 👫 [Associations](reports/associations/)
- [From transcription factors to genes](reports/associations/TF-gene/TF2gene.md)
- [From miRNAs to genes](reports/associations/miRNA-gene/miRNA2gene.md)
- [From methylation sites to genes](reports/associations/met-gene/met2gene.md)
- [From proteins to genes](reports/associations/protein-gene/protein2gene.md)

### 📈 [Integration](reports/integration/)
- [PaintOmics 4](reports/integration/paintomics)
- [MORE]()
- [pls-Multiblock]()

---

## 📓 [Results](results/)

Results tables, relevant features lists, etc... along with a summary of the process.

### 🍳 [Preprocessing](results/preprocessing/)
- [Cooking RNA-Seq data](results/preprocessing/cookingRNASeq)
- [Cooking miRNA-Seq data](results/preprocessing/cookingmiRNASeq)
- [Cooking methylation data](results/preprocessing/cookingMet)
- [Cooking proteomics data](results/preprocessing/cookingProt)

### 👫 [Associations](results/associations/)
- [From transcription factors to genes](results/associations/TF-gene/)
- [From miRNAs to genes](results/associations/miRNA-gene/)
- [From methylation sites to genes](results/associations/met-gene/)
- [From proteins to genes](results/associations/protein-gene/)

### 📈 [Integration](results/integration/)
- [PaintOmics 4](results/integration/paintomics)
- [MORE]()
- [pls-Multiblock]()