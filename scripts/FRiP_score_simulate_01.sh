
cat turn_fragments_pass.tsv|awk '{sum+=$5} END {print "Average = ", sum/NR}' >> average_fragments_length.txt

############  Acquire fragments number per cell（passed_filtered）:
awk '{print $2"\t"$9"\t"$18}' origine_cell_metadata.tsv >> barcode_score.txt
awk '{print $1"\t"$3/$2}' barcode_score.txt >> barcode_real_score.txt
sed -i "1ibarcode\treal_FRiP" barcode_real_score.txt
#sed -i  's/-1//g' origine_cell_metadata.tsv
#sed -i '1d' origine_cell_metadata.tsv
awk '{print $2"_"$9}' origine_cell_metadata.tsv > nFragments_per_cell3679.txt


############  Calculate random FRiP file score:
mkdir random
cd random

cat ../nFragments_per_cell3679.txt |while read n;
do
id=${n%_*}
Nall=${n#*_}
bedtools random -g MH63RS2_chromsize.txt -l $num -n ${Nall} > radom_${id}.bed
bedtools intersect -a radom_${id}.bed -b peaks.narrowPeak  -wa -u | wc -l >> simulated_n_fragments_in_peaks.txt
Nsimu=`bedtools intersect -a radom_${id}.bed -b peaks.narrowPeak  -wa -u | wc -l`
echo ${id}    ${Nsimu}    ${Nall} >> for_FRiP.txt
awk '{print $1"\t"$2/$3}' for_FRiP.txt > barcode_simulated_score.txt
done

sed -i "1ibarcode\trandom_FRiP" barcode_simulated_score.txt


############  Prepare FRiP data for R：
barcode_simulated <- read.table('barcode_simulated_score.txt', header=TRUE,sep="\t")
barcode_real      <- read.table('barcode_real_score.txt', header=TRUE,sep="\t")
frip <- merge(barcode_simulated,barcode_real,by="barcode") 
write.table(frip, file=paste0(args$out_prefix,'merged_FRiP_data.tsv'), quote=FALSE, sep='\t', header=FALSE, row.names=FLASE)

awk '{print "real""\t"$3*100}' merged_FRiP_data.tsv >> all_input_scCT.txt
awk '{print "random""\t"$2*100}' merged_FRiP_data.tsv >> all_input_random.txt
cat all_input_scCT.txt all_input_random.txt >> FRiP_data_file.tsv
sed -i "1igroup\tvalue" FRiP_data_file.tsv

Rscript FRiP_score_simulate_02.R





