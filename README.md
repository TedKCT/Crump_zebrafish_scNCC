# Crump_zebrafish_scNCC

The essential codes and R scripts for building up the "Constellation" map in our publication:<br>
**"Lifelong single-cell profiling of cranial neural crest diversification"** <br>
preprint is available on [bioRxiv](https://www.biorxiv.org/content/10.1101/2021.08.19.456710v1)

Please refer to the manuscript for detailed pipeline of data collection, QC, and processing.

---

## Step 0: Processed data loading

The processed data are available at [FaceBase](https://www.facebase.org/chaise/record/#1/isa:project/RID=3-KG12) as R objects (rds format). We are just going to use the snATAC_SnapATAC rds objects as we already impute the transcriptomes from scRNA to these snATAC_SnapATAC objects. These rds objects are essentially Seurat/Signac objects with the following assays: <br>

- An `[['ATAC']]` assay of *peak-by-cell* matrix with peaks called by [SnapATAC](https://github.com/r3fang/SnapATAC)
- An `[['ACTIVITY']]` assay of imputed *gene activity-by-cell* matrix from corresponding scRNA data anchored by *gmat* calculated from *bmat* of Snap objects
- A `[['chromvar']]` assay of *motif-by-cell* matrix with motif accessibilities calculated by [chromVAR](https://greenleaflab.github.io/chromVAR/articles/Introduction.html) 

---

## Step 1: Constructing Constellation Map

By running the script "*I. Constructing Constellation Map.R*", you will be able to use the imported rds objects to calculate the module scores and further use these scores to construct the Constellation Map as shown in the Fig. 5 in the manuscript. Here is the simple break down of the script:<br>

1. Import the required Seurat/Signac rds objects
2. Calculated the specific enrich peaks in each tissue (cluster) of 14 dpf data
3. Calculated the module scores of every library for each tissue based on their enriched peaks
4. Derived the distance matrix **between-time points matrix**
5. Derived the distance matrices **between-tissue matrices**
6. Calculated the skewness of each tissue module score of every library
7. Build up and output the skewed tissue-time table (will be used in *Step 3*)
8. Expand and combine **between-time points matrix** and **between-tissue matrices** and get the **final distance matrix**
9. Use the **final distance matrix** for dimensional reduction by UMAP (will be used in *Step 3*)
10. Output and plotting

---
