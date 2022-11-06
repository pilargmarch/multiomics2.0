# PaintOmics

The input files can be found in [this folder](/input). The data was input as seen in the following image.

![](/data.input.png)

Some features could not be mapped, as seen in the next image.

![/mapping.info.png]

- ![First attempt](https://www.paintomics.org/?jobID=4S57H2t4b6): no proteomics relevant features were deemed relevant.
- ![Second attempt](https://www.paintomics.org/?jobID=OYJ7cb0d0D): so, proteomics relevant features were selected as those with a q.value<0.05 (without applying a logFC filter, which was 0.5 before).

Note: even though Reactome database was also available for human, only KEGG was used, as KEGG + Reactome kept timing out due to the large amount of data.