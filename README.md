## Introduction of snCUT&Tag bioinformatic pipeline

This is a pipeline for snCUT&Tag data analysis.


## Bioinformatics analysis procedures:

### Map the raw reads to reference genome
Use cellranger-atac count from 10x Genomics to map the reads and generate the possorted_bam.bam and fragments.tsv.gz files for downstream analysis.

### Dimension reduction and cell cluster identification
Use the Main_process.R scripts to conduct the general analysis, such as PCA analysis, dimension reduction, cell clustering, and marker peaker calling and so on.

### Generate bw file for each single cell
Use the Individual_scBW.sh scripts to generate bigwig files for each single cells.

### FRiP score
Use the FRiP_score_simulate_01.sh and FRiP_score_simulate_02.R scripts to calculate the FRiP score.

### ABC model
Predict chromatin interaction using the ABC_model.sh scripts.

## Note
The results obtained may slightly vary when different versions of softwares, such as Cell Ranger ATAC (In this pipeline, Cell Ranger ATAC v1.2 was used) and so on, are used.

For more information and experiment details, please refer to:
- Ouyang W, Luan S, Xiang X, Guo M, Zhang Y, Li G, Li X. Profiling plant histone modification at single-cell resolution using snCUT&Tag. Plant Biotechnol J. 2022, 20(3):420-422. doi: 10.1111/pbi.13768.
