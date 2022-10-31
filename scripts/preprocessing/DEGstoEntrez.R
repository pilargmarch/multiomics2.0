DEGstoEntrez <- function(res, activated.genes, repressed.genes){
  library(AnnotationDbi)
  library(org.Hs.eg.db)
  
  ens.str <- substr(rownames(res), 1, 15)
  res$symbol <- mapIds(org.Hs.eg.db,
                       keys=ens.str,
                       column="SYMBOL",
                       keytype="ENSEMBL",
                       multiVals="first")
  res$entrez <- mapIds(org.Hs.eg.db,
                       keys=ens.str,
                       column="ENTREZID",
                       keytype="ENSEMBL",
                       multiVals="first")
  
  activated.genes <- as.data.frame(activated.genes)
  ens.str <- substr(activated.genes$activated.genes, 1, 15)
  activated.genes$entrez <- mapIds(org.Hs.eg.db,
                                          keys=ens.str,
                                          column="ENTREZID",
                                          keytype="ENSEMBL",
                                          multiVals="first")
  activated.genes <- na.omit(activated.genes)

  repressed.genes <- as.data.frame(repressed.genes)
  ens.str <- substr(repressed.genes$repressed.genes, 1, 15)
  repressed.genes$entrez <- mapIds(org.Hs.eg.db,
                                          keys=ens.str,
                                          column="ENTREZID",
                                          keytype="ENSEMBL",
                                          multiVals="first")
  repressed.genes <- na.omit(repressed.genes)
  
  return(list(res, activated.genes, repressed.genes))
}