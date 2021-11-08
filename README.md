# Crump_zebrafish_scNCC

The essential codes and R scripts for building up the "Constellation" map in our publication:<br>
**"Lifelong single-cell profiling of cranial neural crest diversification"** <br>
preprint is available on [bioRxiv](https://www.biorxiv.org/content/10.1101/2021.08.19.456710v1)

Please refer to the manuscript for detailed pipeline of data collection, QC, and processing.

---

## Processed data loading

The processed data are available at [FaceBase](https://www.facebase.org/chaise/record/#1/isa:project/RID=3-KG12) as R objects (rds format). These rds objects are essentially Seurat objects with the following assays: <br>

- An `[['ATAC']]` assay of *peak-by-cell* matrix with peaks called by SnapATAC
- An `[['ACTIVITY']]` assay of *gene activity-by-cell* matrix with gene activities calculated from *bmat* of the corresponding Snap objects
- An `[['chromvar']]` assay of *motif-by-cell* matrix with motif accessibilities calculated by [chromVAR](https://greenleaflab.github.io/chromVAR/articles/Introduction.html) 

---


