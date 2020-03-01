#============ VARIANT CALLING USING GATK BEST PRACTICES ============
# NIKOLAOS GIALITSIS MSC DATA SCIENCE AND INFORMATION TECHNOLOGIES
#==================================================================
#>>>>>>>>>>>>>>>>>> LOG FILE <<<<<<<<<<<<<<<<<<

#!/bin/env bash
FILES=/home/ubuntu/Desktop/workspace/*
echo '=================== PREPROCESSING AND MAPPING READS TO REF GENOME'


for file in $FILES;do

if [ ${file: -4} != ".bam" ]; then
echo $file
continue
else
temp=${file%%.*}
mv $file $temp.bam
file=$temp
echo $file
fi

echo 'sort bams...'
samtools sort -n $file.bam $file.sorted
echo '...done'

echo 'bam2fastq...'
bamToFastq -i $file.sorted.bam -fq $file.fq1 -fq2 $file.fq2
echo '...done'


sed '1,4d' $file.fq1
sed '1,4d' $file.fq2


echo 'bwa mem fq1 fq2 to sam (paired-end reads)...'
bwa mem -t 4 -M b37.fasta $file.fq1 $file.fq2 > $file.sam
echo '...done'

echo 'convert sam to bam...'
samtools view -bT b37.fasta $file.sam > $file.mapped.bam
echo '...done'

echo 'sort bam...'
samtools sort $file.mapped.bam $file.sorted.mapped
echo '...done'
done



#=======================================================
#=================== MARK DUPLICATES ===================
#=======================================================


echo -e '\n========= MARK DUPLICATES USING PICARD  =======\n\n'
for file in $FILES;do

if [ ${file: -10} != "mapped.bam" ]; then
echo -e '\t\t>>skip' $file
continue
else
file=${file%%.*}
echo -e '\t\t>>mark duplicates in file' $file.sorted.mapped.bam '...'
fi
#echo -e '\t\t>>...done\n'

java -Xmx2g  -jar MarkDuplicates.jar  INPUT=$file.sorted.mapped.bam  OUTPUT=$file.dedup_reads.bam  METRICS_file=$file.metrics.txt
echo -e '\t\t>>...done\n'
done





#========================================================
#================== FIX READ GROUPS =====================
#========================================================
#
for file in $FILES;do

if [ ${file: -22} != "dedup_reads.sorted.bam" ]; then
echo 'skip' $file
continue
else
file=${file%%.*}
echo 'Process Read Groups in file'$file.dedup_reads.sorted.bam
fi

java -jar AddOrReplaceReadGroups.jar \
      I=$file.dedup_reads.sorted.bam \
      O=$file.bam \
      RGID=1 \
      RGLB=lib1 \
      RGPL=illumina \
      RGPU=unit1 \
      RGSM=20
done


#========================================
#=========== Produce .dict file =========
#========================================

#java -jar NormalizeFasta.jar \
#      I=b37.fa \
#      O=normalized_b37.fa

java -jar CreateSequenceDictionary.jar R=b37.fasta O=b37.dict
samtools faidx b37.fasta



#======================================
#======  BASE RECALIBRATION ==========
#=====================================
#!/bin/env bash
echo -e '\n\tRetrieve known SNP sites...'
wget -c ftp://gsapubftp-anonymous@ftp.broadinstitute.org/bundle/b37/dbsnp_138.b37.vcf.gz
echo -e '\t>>loaded dbsnp_138.b37.vcf.gz'
wget -c ftp://gsapubftp-anonymous@ftp.broadinstitute.org/bundle/b37/Mills_and_1000G_gold_standard.indels.b37.vcf.gz
echo -e '\t>>loaded Mills_and_1000G_gold_standards_indels.b37.vcf.gz'

echo -e '\n Unpack Gzip files...'
gunzip -df dbsnp_138.b37.vcf.gz
gunzip -df Mills_and_1000G_gold_standard.indels.b37.vcf.gz
echo -e '\t>>ok\n'

echo -e '\n Index known SNP sites...\n'
gatk IndexFeatureFile -I dbsnp_138.b37.vcf
gatk IndexFeatureFile -I Mills_and_1000G_gold_standard.indels.b37.vcf
echo -e '\t>>ok\n'

FILES=/home/ubuntu/Desktop/workspace/*
echo -e '\n========= BASE RECALIBRATION  =======\n\n'
for file in $FILES;do

if [ ${file: -4} != ".bam" ]; then
echo -e '\t\t>>skip' $file
continue
else
file=${file%%.*}
#samtools index $file.bam
echo -e '\t\t>>Base Recalibration' $file.bam '...'
fi

if test -f $file.recal_data.table; then
    echo "$FILE.bam has already been recalibrated"
    continue
fi

samtools index /home/ubuntu/Desktop/workspace/$file.bam

gatk BaseRecalibrator \
-R b37.fasta \
-I $file.bam \
-L 11 \
--known-sites dbsnp_138.b37.vcf \
--known-sites Mills_and_1000G_gold_standard.indels.b37.vcf \
-O $file.recal_data.table
echo -e '\t\t>>...done\n'
done


#===============================================
#========  Analysis Ready Reads ===============
#===============================================
echo -e '\n========= PREPARE ANALYSIS  =======\n\n'
for file in $FILES;do

if [ ${file: -4} != ".bam" ]; then
echo -e '\t\t>>skip' $file
continue
else
file=${file%%.*}
fi

if test -f $file.recal_data.table; then
    echo -e '\n\tPrint Reads for file' $file.bam
    gatk PrintReads -R b37.fasta -I $file.bam -O $file.recal_reads.bam
fi


#========================================
#==========  HAPLOTYPE CALLER ===========
#========================================

echo -e '\n========= Haplotype Caller  =======\n\n'
for file in $FILES;do

if [ ${file: -4} != ".bam" ]; then
echo -e '\t\t>>skip' $file
continue
else
file=${file%%.*}
#samtools index $file.bam
echo -e '\t\t>>get .vcf for file' $file.bam '...'
fi

gatk HaplotypeCaller -R b37.fasta -I $file.recal_reads.bam -O $file.variants.g.vcf -ERC GVCF
done

