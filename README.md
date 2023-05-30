# Introduction of snCUT&Tag bioinformatic pipeline

This is a pipeline for snCUT&Tag data analysis.


# Bioinformatics analysis procedures:

### Prepare a config file (RefGenome.config) for reference genome creation as below:

- {
  GENOME_FASTA_INPUT: "RefGenome.fasta",
  GENE_ANNOTATION_INPUT: "RefGenome.gff3",
  MOTIF_INPUT: "",
  ORGANISM: "Oryza Sativa",
  PRIMARY_CONTIGS: ["Chr01","Chr02","Chr03","Chr04","Chr05","Chr06","Chr07","Chr08","Chr09","Chr10","Chr11","Chr12"],
  NON_NUCLEAR_CONTIGS: []
  }

### Construct a custom reference genome using the cellranger-atac mkref function.

- $ cellranger-atac mkref RefGenomeName --config RefGenome.config

### Map the sequencing reads (fastq files in the DataDir directory) to reference genome and count barcodes using the cellranger-atac count function. The function also perform analysis including identification of Tn5 cut sites, detection of peaks, cell calling, and count matrix generation. As a result, files such as possorted_bam.bam and fragments.tsv.gz are generated for downstream analysis.

- $ cellranger-atac count --id=idName --reference=Path/RefGenomeName --fastqs=Path/DataDir

### General analysis using R packages. Input files generated from the cellreanger-atac count outputs in Step 45:

- $ fragments <− 'fragments.tsv.gz'
- $ all_barcodes_file <− 'all_barcodes.txt'
- $ peak_barcodes_file <− 'peaks_barcodes.txt'
- $ metadata_file <− 'metadata.csv'

### Create genome annotation files

- $ rice_gff <− import.gff(con = "RefGenome.gtf")
- $ seqlengths(rice_gff) <− c("Chr01" = 44512328,"Chr02"=36671280,"Chr03"=39351490,"Chr04"=36167251,"Chr05"=30881543,"Chr06"=31652805,"Chr07"=29891017,"Chr08"=29797537,"Chr09"=24332368,"Chr10"=25127214,"Chr11"=32883170,"Chr12"=26156356)
- $ genebody.coords <− rice_gff [rice_gff$type=="gene",]
- $ genebody.coords[is.na(genebody.coords$symbol),]$symbol <− genebody.coords[is.na(genebody.coords$symbol),]$Name
- $ genebody.coords.flat <− GenomicRanges::reduce(x = genebody.coords)
- $ genebodyandpromoter.coords.flat <− Signac::Extend(genebody.coords.flat,upstream = 2000)
- $ genebodyandpromoter.coords.flat$name <− genebody.coords[nearest(genebodyandpromoter.coords.flat,genebody.coords)]$symbol

### Generate a gene-activity matrix for identifying high-quality cells according to fragments.tsv.gz file using FeatureMatrix function of the R package Seurat.

- $ gene.matrix <− FeatureMatrix(fragments = fragments, features = genebodyandpromoter.coords.flat, cells = gsub(paste0(args$sample,"_"),"",colnames(seurat_object)))

### Create a Seurat object using the CreatSeuratObject function of the R package Seurat.

- $ seurat_object_gene <− CreateSeuratObject(counts = gene.matrix, assay = 'GA', min.features = min_features, min.cells = min_cells)

### Normalize the count data.

- $ seurat_object_gene <− NormalizeData(seurat_object_gene, normalization.method = 'LogNormalize', scale.factor=10000)

### Find high variable features between cells.

- $ seurat_object_gene <− FindVariableFeatures(seurat_object_gene, selection.method = "vst", nfeatures = 1000)

### Scale the data.

- $ seurat_object_gene <− ScaleData(seurat_object_gene, features = rownames(seurat_object_gene))

### PCA dimension reduction.

- $ seurat_object_gene <− RunPCA(seurat_object_gene, features = VariableFeatures(object = seurat_object_gene))
- $ VizDimLoadings(seurat_object_gene, dims = 1:2, reduction = "pca")
- $ DimPlot(seurat_object_gene, reduction = "pca")
- $ DimHeatmap(seurat_object_gene, dims = , cells = , balanced = TRUE)

### Run the Uniform Manifold Approximation and Projection (UMAP) dimensional reduction.

- $ seurat_object_gene <− FindNeighbors(seurat_object_gene, dims = 2:20)
- $ seurat_object_gene <− FindClusters(seurat_object_gene, resolution = 0.6)
- $ seurat_object_gene <− RunUMAP(seurat_object_gene, dims = 2:20)
- $ DimPlot(seurat_object_gene, reduction = "umap")

### Find marker peaks for each cluster and visualize the results either by bubble plot or heatmap.

- $ seurat_object_gene.markers <− FindAllMarkers(seurat_object_gene, only.pos = TRUE, min.pct = 0.1, logfc.threshold = 0.25)
- $ top5 <− seurat_object_gene.markers %>% dplyr::filter(p_val_adj < 0.05) %>% group_by(cluster) %>% top_n(n = 3, wt = avg_log2FC)
- $ DotPlot(seurat_object_gene, features = unique(top5$gene)) + RotatedAxis()
- $ DoHeatmap(seurat_object_gene, features = top5$gene) + NoLegend() + theme(legend.position = "none", axis.text.y = element_text(size = 6))

### Generate bw files for visualization of each single cells.

- $ bash Individual_scBW.sh

### Calculate FRiP cores and visualize using violin plots.

- $ bash FRiP_score_simulate_01.sh
- $ Rscript FRiP_score_simulate_02.R

### Predict chromatin interactions using the ABC model.
- $ bash ABC_model.sh

# Note
The results obtained may slightly vary when different versions of softwares, such as Cell Ranger ATAC (In this pipeline, Cell Ranger ATAC v1.2 was used) and so on, are used.

For more information and experiment details, please refer to:
- Ouyang W, Luan S, Xiang X, Guo M, Zhang Y, Li G, Li X. Profiling plant histone modification at single-cell resolution using snCUT&Tag. Plant Biotechnol J. 2022, 20(3):420-422. doi: 10.1111/pbi.13768.


