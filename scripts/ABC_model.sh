
module load MACS2/2.1.1

macs2 callpeak \
-t RMAS02.bam \
-n sample_name \
-f BAM \
-g 3.6e8 \
-p .1 \
--call-summits \
--outdir /PATH_TO_OUTDIR/ 


############### Step 1. Define candidate elemets
for i in *narrowPeak
do
sample=${i%_*}
bedtools sort -faidx MH63RS2_chromsize.txt -i ${sample}_peaks.narrowPeak > ${sample}_peaks.narrowPeak.sorted
done

## Call candidate regions
conda activate final-abc-env

python makeCandidateRegions.py \
--narrowPeak ${sample}_peaks.narrowPeak.sorted \
--bam ${sample}_sort.bam \
--outDir .../ABC/${sample}/ \
--chrom_sizes MH63RS2_chromsize.txt \
--peakExtendFromSummit 250 \
--nStrongestPeaks 150000


############### Step 2. Quantifying Enhancer Activity:
for i in *candidateRegions.bed
do
sample=${i%_*}

python run.neighborhoods.py \
--candidate_enhancer_regions ${sample}_peaks.narrowPeak.sorted.candidateRegions.bed \
--genes MH63RS2_gene.bed \
--H3K27ac H3K27ac_RMCS96.bam \
--ATAC ${sample}_sort.bam \
--chrom_sizes MH63RS2_chromsize \
--cellType ${sample} \
--outdir .../02_Quantifying_Enhancer_Activity/
done


############### Step 3. Computing the ABC Score:
python predict.py \
--enhancers EnhancerList.txt \
--genes GeneList.txt \
--HiCdir .../SRR6765292/bedpes/ \
--hic_type bedpe \
--hic_resolution 5000 \
--scale_hic_using_powerlaw \
--threshold .02 \
--cellType sample_cell \
--outdir .../03_Computing_ABC_Score/ \
--chrom_sizes MH63RS2_chromsize_file \
--make_all_putative



