#!/bin/bash
set -e -o pipefail
##########################################
# load environment variables for LRSDAY
source ./../../env.sh

###########################################
# set project-specific variables
prefix="SK1" # file name prefix for the output files
raw_assembly="./../01.Long-read-based_Genome_Assembly/SK1.canu.fa" # path of the raw pacbio assembly before pilon correction
threads=1  # number of threads to use
mode="PE" # Illumina sequencing mode, "PE" for paired-end sequencing and "SE" for single-end sequencing
fixlist="snps,indels" # types of errors for pilon correction; see pilon's manual for more details
if [[ $mode == "PE" ]]
then
    reads_PE1="./../00.Illumina_Reads/SRR4074258_pass_1.fastq.gz" # please replace the PE reads file name for your own project
    reads_PE2="./../00.Illumina_Reads/SRR4074258_pass_2.fastq.gz" # please replace the PE reads file name for your own project
else
    reads_SE="./../00.Illumina_Reads/sample_pass_1.fastq.gz" # please replace the SE reads file name for your own project if you only have SE data
fi
debug="no" # use "yes" if prefer to keep intermediate files, otherwise use "no".

###########################################
# process the pipeline
if [[ $mode == "PE" ]]
then
    adapter="$trimmomatic_dir/adapters/TruSeq3-PE-2.fa" # adapter for PE reads
    ln -s $reads_PE1 raw.R1.fq.gz;
    ln -s $reads_PE2 raw.R2.fq.gz;
else
    adpater="$trimmomatic_dir/adapters/TruSeq3-SE.fa" # adapter for SE reads
    ln -s $reads_SE raw.fq.gz;
fi

ln -s $raw_assembly refseq.fa

cp $adapter adapter.fa

mkdir tmp
if [[ $mode == "PE" ]]
then
    java -Djava.io.tmpdir=./tmp -XX:ParallelGCThreads=$threads -jar $trimmomatic_dir/trimmomatic.jar PE -threads $threads -phred33  raw.R1.fq.gz  raw.R2.fq.gz  trimmed.R1.fq.gz trimmed.unpaired.R1.fq.gz trimmed.R2.fq.gz trimmed.unpaired.R2.fq.gz ILLUMINACLIP:adapter.fa:2:30:10  SLIDINGWINDOW:5:20 MINLEN:36
    rm trimmed.unpaired.R1.fq.gz
    rm trimmed.unpaired.R2.fq.gz 
else
    java -Djava.io.tmpdir=./tmp -XX:ParallelGCThreads=$threads -jar $trimmomatic_dir/trimmomatic.jar SE -threads $threads -phred33  raw.fq.gz  trimmed.fq.gz ILLUMINACLIP:adapter.fa:2:30:10 SLIDINGWINDOW:5:20 MINLEN:36
fi

# bwa mapping
$bwa_dir/bwa index refseq.fa

if [[ $mode == "PE" ]]
then
    $bwa_dir/bwa mem -t $threads -M refseq.fa  trimmed.R1.fq.gz trimmed.R2.fq.gz >$prefix.sam
else
    $bwa_dir/bwa mem -t $threads -M refseq.fa  trimmed.fq.gz >$prefix.sam
fi

# index reference sequence
$samtools_dir/samtools faidx refseq.fa
java -Djava.io.tmpdir=./tmp -XX:ParallelGCThreads=$threads -jar $picard_dir/picard.jar CreateSequenceDictionary  REFERENCE=refseq.fa OUTPUT=refseq.dict

# sort bam file by picard-tools SortSam
java -Djava.io.tmpdir=./tmp -XX:ParallelGCThreads=$threads -jar $picard_dir/picard.jar SortSam  INPUT=$prefix.sam OUTPUT=$prefix.sort.bam SORT_ORDER=coordinate

# fixmate
java -Djava.io.tmpdir=./tmp -XX:ParallelGCThreads=$threads -jar $picard_dir/picard.jar FixMateInformation INPUT=$prefix.sort.bam OUTPUT=$prefix.fixmate.bam

# add or replace read groups and sort
java -Djava.io.tmpdir=./tmp -XX:ParallelGCThreads=$threads -jar $picard_dir/picard.jar AddOrReplaceReadGroups INPUT=$prefix.fixmate.bam OUTPUT=$prefix.rdgrp.bam \
    SORT_ORDER=coordinate RGID=$prefix RGLB=$prefix RGPL="Illumina" RGPU=$prefix RGSM=$prefix RGCN="RGCN"

# remove duplicates
java -Djava.io.tmpdir=./tmp -XX:ParallelGCThreads=$threads -jar $picard_dir/picard.jar MarkDuplicates INPUT=$prefix.rdgrp.bam REMOVE_DUPLICATES=true  \
    METRICS_FILE=$prefix.dedup.matrics  OUTPUT=$prefix.dedup.bam 

# index the dedup.bam file
$samtools_dir/samtools index $prefix.dedup.bam

# GATK local realign
# find realigner targets
java -Djava.io.tmpdir=./tmp -XX:ParallelGCThreads=$threads -jar $gatk_dir/GenomeAnalysisTK.jar -R refseq.fa -T RealignerTargetCreator -I $prefix.dedup.bam  -o $prefix.realn.intervals
# run realigner
java -Djava.io.tmpdir=./tmp -XX:ParallelGCThreads=$threads -jar $gatk_dir/GenomeAnalysisTK.jar -R refseq.fa -T IndelRealigner  -I $prefix.dedup.bam -targetIntervals $prefix.realn.intervals  -o $prefix.realn.bam

# index final bam file
$samtools_dir/samtools index $prefix.realn.bam

# for PE sequencing
if [[ $mode == "PE" ]]
then
    java -Djava.io.tmpdir=./tmp -Xmx16G -XX:ParallelGCThreads=$threads -jar $pilon_dir/pilon.jar --genome $raw_assembly --frags $prefix.realn.bam  --fix $fixlist --vcf --changes --output $prefix.pilon  >$prefix.log
else
    java -Djava.io.tmpdir=./tmp -Xmx16G -XX:ParallelGCThreads=$threads -jar $pilon_dir/pilon.jar --genome $raw_assembly --unpaired $prefix.realn.bam  --fix $fixlist --vcf --changes --output $prefix.pilon  >$prefix.log
fi

perl $LRSDAY_HOME/scripts/summarize_pilon_correction.pl -i $prefix.pilon.changes

mv $prefix.pilon.fasta $prefix.pilon.fa
rm -r tmp

# clean up intermediate files
if [[ $debug == "no" ]]
then
    rm adapter.fa
    rm refseq.fa
    rm refseq.fa.fai
    rm refseq.dict
    rm refseq.fa.bwt
    rm refseq.fa.pac
    rm refseq.fa.ann
    rm refseq.fa.amb
    rm refseq.fa.sa
    rm $prefix.sam
    rm $prefix.sort.bam
    rm $prefix.fixmate.bam
    rm $prefix.rdgrp.bam
    rm $prefix.dedup.bam
    rm $prefix.dedup.matrics
    rm $prefix.dedup.bam.bai
    rm $prefix.realn.intervals
    rm *.fq.gz
    rm $prefix.realn.bam.bai
fi

   
############################
# checking bash exit status
if [[ $? -eq 0 ]]
then
    echo ""
    echo "LRSDAY message: This bash script has been successfully processed! :)"
    echo ""
    echo ""
    exit 0
fi
############################
