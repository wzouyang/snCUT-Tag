
cat("*** Loading libraries")
suppressMessages({
library(argparse);
library(gridExtra);
library(Seurat);
library(Signac);
library(viridis);
library(ggplot2);
library(dplyr);
library(rtracklayer);
library(yaml);
library(pheatmap);
library(GenomicRanges);
  }
)

set.seed(1)

#Note that this pipeline use cellranger-atac count(v1.2) results as input

########### Arguments parser
parser <- ArgumentParser()
parser$add_argument("-s", "--sample", type="character", default='snCUT&Tag', help="sample name [as in config file key]")
#parser$add_argument("-c", "--config", type="character", default='', help="maximum number of reads in cell")
parser$add_argument("-o", "--out_prefix", type="character", default="", help="folder for the output in clustering_snakemake folder")
parser$add_argument("-w", "--window", type="integer", default=1000, help="width of the genome window (if -t bins)")
#args <- parser$parse_args()
print(args)
config <- yaml::read_yaml(args$config)


########### Cell filter parameters
cutoff_reads_min            = 4000
cutoff_reads_max            = 300000
cutoff_peak_percentage_low  = 0.2
cutoff_peak_percentage_high = 0.98
ndim = 50
window  = args$window
assay = "peaksMB"


########### input files
fragments <- 'fragments.tsv.gz'
all_barcodes_file  <- 'all_barcodes.txt'
peak_barcodes_file <- 'peaks_barcodes.txt'
metadata_file      <- 'metadata.csv'


#### Create genome annotation
rice_gff <-import.gff(con = "MH63RS2.LNNK00000000.v2_4.gtf")
seqlengths(rice_gff) <- c("Chr01" = 44512328,"Chr02"=36671280,"Chr03"=39351490,"Chr04"=36167251,"Chr05"=30881543,"Chr06"=31652805,"Chr07"=29891017,"Chr08"=29797537,"Chr09"=24332368,"Chr10"=25127214,"Chr11"=32883170,"Chr12"=26156356)
rice_gff_chr <- renameSeqlevels(rice_gff, c(Chr01="chr1",Chr02="chr2",Chr03="chr3",Chr04="chr4",Chr05="chr5",Chr06="chr6",Chr07="chr7",Chr08="chr8",Chr09="chr9",Chr10="chr10",Chr11="chr11",Chr12="chr12"))
genebody.coords<-rice_gff_chr[rice_gff_chr$type=="gene",]
genebody.coords[is.na(genebody.coords$symbol),]$symbol<-genebody.coords[is.na(genebody.coords$symbol),]$Name
genebody.coords.flat <- GenomicRanges::reduce(x = genebody.coords)
genebodyandpromoter.coords.flat <- Signac::Extend(genebody.coords.flat,upstream = 2000)
genebodyandpromoter.coords.flat$name<- genebody.coords[nearest(genebodyandpromoter.coords.flat,genebody.coords)]$symbol


########## remould metadata file information
metadata = read.csv(metadata_file, header = 1)
metadata = metadata[2:nrow(metadata),]
metadata$logUMI = log10(metadata$passed_filters + 1)
metadata$promoter_ratio = (metadata$promoter_region_fragments+1) / (metadata$passed_filters + 1)
metadata$peak_region_ratio = (metadata$peak_region_fragments+1) / (metadata$passed_filters + 1)

all_barcodes <- read.table(file=all_barcodes_file)
peak_barcodes <- read.table(file=peak_barcodes_file)
bcd <- merge(all_barcodes,peak_barcodes,by="V2")  
colnames(bcd) <- c("barcode","all_unique_reads","peak_reads")
bcd$peak_ratio_reads <- bcd$peak_reads/bcd$all_unique_reads
bcd$sample <- args$sample

metadata <- merge(metadata,bcd,by='barcode')
metadata <- metadata[metadata$is__cell_barcode==1,]
write.table(metadata, file=paste0(args$out_prefix,'origine_cell_metadata.tsv'),row.names = FALSE,quote = FALSE, sep='\t')


########## Cell filter
metadata$is__cell_barcode <- as.factor(metadata$is__cell_barcode)
metadata[,"passed"] <- FALSE
metadata[metadata$all_unique_reads > 10^cutoff_reads_min &
         metadata$all_unique_reads < 10^cutoff_reads_max &
         metadata$peak_ratio_reads > cutoff_peak_percentage_low &
         metadata$peak_ratio_reads < cutoff_peak_percentage_high,"passed"] <- TRUE
metadata <- metadata[metadata$passed,]
rownames(metadata) <- metadata$barcode
write.table(metadata, file=paste0(args$out_prefix,'filtered_cell_metadata.tsv'),row.names = FALSE,quote = FALSE, sep='\t')
#cell_num <- nrow(metadata)


fragments_gr     <- rtracklayer::import(fragments,format = "bed")
fragments.pass   <- fragments_gr[fragments_gr$name %in% barcode_pass]
fragments.pass   <- as.data.frame(fragments.pass)
write.table(fragments.pass, file=paste0(args$out_prefix,'fragments_pass.tsv'), sep='\t', header=FALSE)


gene.matrix     <- FeatureMatrix(fragments = fragments,
                                 features = genebodyandpromoter.coords.flat,
                                 cells = gsub(paste0(args$sample,"_"),"",colnames(seurat_object)))


seurat_object_gene <- CreateSeuratObject(counts = gene.matrix,
                     project = args$sample,
                     assay = 'GA',
                     min.features = min_features,
                     min.cells = min_cells)


seurat_object_gene <- NormalizeData(seurat_object_gene,normalization.method = 'LogNormalize',scale.factor=10000)

###### Find and show high variable genes
seurat_object_gene <- FindVariableFeatures(seurat_object_gene, selection.method = "vst", nfeatures = 1000)
top15 <- head(VariableFeatures(seurat_object_gene),15)
plot1 <- VariableFeaturePlot(seurat_object_gene)
plot2 <- LabelPoints(plot = plot1, points = top15, repel = TRUE)
CombinePlots(plots = list(plot1, plot2))


###### Data Scale
all.genes <- rownames(seurat_object_gene)
seurat_object_gene <- ScaleData(seurat_object_gene, features = all.genes)

###### PCA
seurat_object_gene <- RunPCA(seurat_object_gene, features = VariableFeatures(object = seurat_object_gene))

###### Visualization
VizDimLoadings(seurat_object_gene, dims = 1:2, reduction = "pca")
DimPlot(seurat_object_gene, reduction = "pca")
DimHeatmap(seurat_object_gene, dims = 1:12, cells = 500, balanced = TRUE, ncol = 2)

###### Dimention decision
seurat_object_gene <- JackStraw(seurat_object_gene, num.replicate = 100)
seurat_object_gene <- ScoreJackStraw(seurat_object_gene, dims = 1:20)
JackStrawPlot(seurat_object_gene, dims = 1:15)
ElbowPlot(seurat_object_gene)


###### Cell clustering
seurat_object_gene <- FindNeighbors(seurat_object_gene, dims = 2:20)
seurat_object_gene <- FindClusters(seurat_object_gene, resolution = 0.6)
seurat_object_gene <- RunUMAP(seurat_object_gene, dims = 2:20)

DimPlot(seurat_object_gene, reduction = "umap")


###### Acquire DEX gene info
seurat_object_gene.markers <- FindAllMarkers(seurat_object_gene, only.pos = TRUE, min.pct = 0.1, logfc.threshold = 0.25)
top20_marker_gene <- seurat_object_gene.markers %>% dplyr::filter(p_val_adj < 0.1) %>% group_by(cluster) %>% top_n(n = 20, wt = avg_log2FC)
write.csv(top20_marker_gene, file = paste0(args$out_prefix, "top_markers_by_GA.csv"), sep="\t", quote=F)


###### Heatmap of top variable genes
DoHeatmap(seurat_object_gene, features = top5$gene) + NoLegend() + theme(legend.position = "none", axis.text.y = element_text(size = 6))
DotPlot(seurat_object_gene, features = unique(top2$gene)) + RotatedAxis()


###### Export bw file per cluster
fragments.path <- 'fragments.tsv.gz'
fragments <- rtracklayer::import(con = fragments.path,format = 'bed')
chrom.sizes <- read.table("Chromosize.txt",,sep="\t",stringsAsFactors = FALSE)
chrom.sizes <- chrom.sizes[1:12,]
exportBW <- function(object,cluster,fragments){
  if(class(object) == "Seurat"){
    cells <- rownames(object@meta.data[object@active.ident == cluster,])
  }
  
  new_read <- GRanges(seqnames = chrom.sizes[,1], 
        ranges =IRanges(start = as.numeric(chrom.sizes[,2]),
                        width=1),
        name = rep("in_silico_extra_read",dim(chrom.sizes)[1]),
        score = rep(0,dim(chrom.sizes)[1])
        )
  
  fragments.x <- fragments$name %in% cells
  fragments.x <- fragments[fragments.x]
  fragments.x <- c(fragments.x,new_read)
  
  coverage.x <- coverage(fragments.x)
  
  # Option A - normalize by number of reads per sample
  coverage.x <- coverage.x/(length(fragments.x)/1e6)
  # Option B - normalize by mean signal (~ enrichment of mean signal)
  # coverage.x <- coverage.x / mean(unlist(coverage.x))
  rtracklayer::export.bw(object = coverage.x,paste0(args$out_prefix,'bw_per_cluster',"/cluster_",cluster,".bw"))
}

lapply(levels(seurat_object_gene@active.ident),function(x){
  exportBW(seurat_object_gene,x,fragments)
})


