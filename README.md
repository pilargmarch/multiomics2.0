# multiomics2.0
![Charting a multi-omic universe. Image credits: Isabel Romero Calvo/EMBL.](scripts/cover.jpg)
<sup>[Charting a multi-omic universe. Image credits: Isabel Romero Calvo/EMBL.](https://www.embl.org/news/science/charting-a-multi-omic-universe/)</sup>

This is the repository for my Master's Final Thesis. First, breast cancer RNA-Seq, miRNA-Seq, methylation and proteomics data was selected and downloaded from TCGA-BRCA. After the data was processed and tested for differential expression, regulatory associations between these omics were retrieved from various databases. Finally, relevant features, associations and gene expression were all put into different methods to try to come up with a regulatory model for gene expression in breast cancer.

## :dna: [Data](data/)
RNA-Seq, miRNA-Seq, methylation and proteomics data used in this project, both as raw and cooked.

---

## :computer: [Scripts](scripts/)
`.Rmd` notebooks (knitted as `.md`) used to generate the obtained results.

### 🍳 [Preprocessing](scripts/preprocessing/)
- [downloadingTCGA.md](scripts/preprocessing/downloadingTCGA.md)
- [cookingRNASeq.md](scripts/preprocessing/cookingRNASeq.md)
- [cookingmiRNASeq.md](scripts/preprocessing/cookingmiRNASeq.md)
- [cookingMet.md](scripts/preprocessing/cookingMet.md)
- [cookingProt.md](scripts/preprocessing/cookingProt.md)

### 👫 [Associations](scripts/associations/)
- [TF2gene.md](scripts/associations/TF-gene/TF2gene.md)
- [miRNA2gene.md](scripts/associations/miRNA-gene/miRNA2gene.md)
- [met2gene.md](scripts/associations/met-gene/met2gene.md)
- [protein2gene.md](scripts/associations/protein-gene/protein2gene.md)

### 📈 [Integration](scripts/integration/)
- [PaintOmics 4](scripts/integration/paintomics)

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