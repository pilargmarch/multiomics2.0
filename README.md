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
- [Transcription factors to genes](scripts/associations/TF2gene.R)
- [miRNAs to genes](scripts/associations/miRNA2gene.R)
- [Methylation/CpG sites to genes](scripts/associations/met2gene.R)
- [Proteins to genes](scripts/associations/protein2gene.R)

### 📈 [Models](scripts/models/)
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
- [Transcription factors to genes](reports/associations/TF-gene/TF2gene.md)
- [miRNAs to genes](reports/associations/miRNA-gene/miRNA2gene.md)
- [Methylation/CpG sites to genes](reports/associations/met-gene/met2gene.md)
- [Proteins to genes](reports/associations/protein-gene/protein2gene.md)

### 📈 [Models](reports/models/)
- [PaintOmics 4](reports/models/paintomics)
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
- [Transcription factors to genes]()
- [miRNAs to genes]()
- [Methylation/CpG sites to genes]()
- [Proteins to genes]()

### 📈 [Models](results/models/)
- [PaintOmics 4](results/models/paintomics)
- [MORE]()
- [pls-Multiblock]()