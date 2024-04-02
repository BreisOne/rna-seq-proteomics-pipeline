# Multiomic integration pipeline for the study of TGF-B pathway in Alström syndrome.

This multiomic integration project combines RNA-seq and LFQ proteomics data to explore gene and protein interactions in a cellular model context. The project focuses on studying TGF-B signaling in Alström syndrome.

## Objective

The main objective of this project is to identify differentially expressed genes and proteins between conditions and explore enriched biological pathways in Alström syndrome. Additionally, the project aims to identify common enriched pathways among the differentially expressed genes and proteins.

## Analysis Steps

1. **Gene/Protein Overlap**:
   - Perform an overlap analysis between the genes and proteins identified in the RNA-seq and LFQ proteomics data to determine the degree of coverage between the two omics.

2. **Differential Expression (DE) Analysis**:
   - Perform differential expression analysis for RNA-seq and LFQ proteomics data separately using the DESeq2 and DEP packages, respectively.

3. **Over-Representation Analysis (ORA) for RNA-seq and LFQ proteomics**:
   - Perform over-representation analysis (ORA) using gene ontologies (GO-terms) and different pathway databases (Bioplanet, KEEG, MSigDB and WikiPathway) to identify biological pathways enriched in differentially expressed genes and proteins. This step was performed through the web application [enrichR](https://maayanlab.cloud/Enrichr/).

4. **Over-Representation Analysis (ORA) for Common Genes**:
   - Repeat the over-representation analysis, this time using the common differentially expressed genes/proteins between RNA-seq and proteomics.

## File Structure

The project follows an R package structure, facilitating code organization and reusability. The directory structure is as follows:
```
project/
│
├── R/ # Directory for R code files
│
├── data/ # Directory for storing data
│
├── results/ # Directory for results
│ └── figures/ # Directory for figures of the analysis
│
├── scripts/ # Directory with other scripts used in the analysis
  └── rna-seq/ # Directory with the scripts used to 
               # preprocess fastq files and generate 
               # reads count matrix.
```
## Citation

If you find this repository useful for your research, please consider citing the corresponding article:

```
Bea-Mascato, B., Gómez-Castañeda, E., Sánchez-Corrales, Y.E. et al. Loss of the centrosomal protein ALMS1 
alters lipid metabolism and the regulation of extracellular matrix-related processes. Biol Direct 18, 84 (2023). 
https://doi.org/10.1186/s13062-023-00441-2
```