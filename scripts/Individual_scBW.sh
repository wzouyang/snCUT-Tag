
samtools view -h -o out.sam 10×out.bam

cat barcode_list.txt |while read barcode;
do
head -n 20  out.sam > sam/$barcode.sam
grep $barcode out.sam >> sam/$barcode.sam
done


module load picard/2.1.1-Java-1.8.0_92
module load deepTools/2.5.3

for i in *sam
do
sample=${i%.*}
samtools view -h -b -S  ${sample}.sam   |  samtools sort  -  |  samtools view -h -bq 30 -  >  ${sample}.bam
java -jar ${EBROOTPICARD}/picard.jar MarkDuplicates REMOVE_DUPLICATES=true I=./${sample}.bam O=./${sample}_rmdup_picard.bam M=./${sample}_rmdup_picard.txt > ./${sample}_picard.log 2>&1
num1=10000000
num2="$(samtools view -c  $i  2>&1 )"
res=$(printf "%.5f" `echo "scale=5;$num1/$num2"|bc`)
samtools index -b ${sample}_rmdup_picard.bam && bamCoverage --scaleFactor $res -b ${sample}_rmdup_picard.bam -o ./${sample}.10M.bw
done


