#!/bin/bash
set -e -o pipefail
##########################################
# load environment variables for LRSDAY
source ./../../env.sh

###########################################
# set project-specific variables

prefix="SK1" # The file name prefix for the processing sample. Default = "SK1" for the testing example.
input_assembly="./../02.Long-read-based_Assembly_Polishing/$prefix.assembly.long_read_polished.fa" # The file path of the input assembly before Illumina-based correction
trim_illumina_reads="yes" # Whether to trim the input Illumina reads. Use "yes" if prefer to perform trimming, otherwise use "no". Default = "yes". 
rounds_of_successive_polishing=1 # The number of total rounds of Illumina-read-based assembly polishing. Default = "1" for the testing example.
threads=1 # The number of threads to use. Default = "1".
mode="PE" # Illumina sequencing mode, "PE" for paired-end sequencing and "SE" for single-end sequencing. Default = "PE".
fixlist="snps,indels" # The types of errors for Illumina-read-based correction by Pilon; see Pilon's manual for more details. Default = "snps,indels".
if [[ $mode == "PE" ]]
then
    reads_PE1="./../00.Illumina_Reads/SRR4074258_pass_1.fastq.gz" # Please replace the PE reads file name for your own project
    reads_PE2="./../00.Illumina_Reads/SRR4074258_pass_2.fastq.gz" # Please replace the PE reads file name for your own project
else
    reads_SE="./../00.Illumina_Reads/sample_pass_1.fastq.gz" # Please replace the SE reads file name for your own project if you only have SE data
fi
debug="no" # Whether to keep intermediate files for debugging. Use "yes" if prefer to keep intermediate files, otherwise use "no". Default = "no".

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

cp  $input_assembly refseq.tmp.fa
cp $adapter adapter.fa

mkdir tmp

if [[ $trim_illumina_reads == "yes" ]]
then
    if [[ $mode == "PE" ]]
    then
	java -Djava.io.tmpdir=./tmp -XX:ParallelGCThreads=$threads -jar $trimmomatic_dir/trimmomatic.jar PE -threads $threads -phred33  raw.R1.fq.gz  raw.R2.fq.gz  trimmed.R1.fq.gz trimmed.unpaired.R1.fq.gz trimmed.R2.fq.gz trimmed.unpaired.R2.fq.gz ILLUMINACLIP:adapter.fa:2:30:10  SLIDINGWINDOW:5:20 MINLEN:36
	rm trimmed.unpaired.R1.fq.gz
	rm trimmed.unpaired.R2.fq.gz 
	mv trimmed.R1.fq.gz clean.R1.fq.gz
	mv trimmed.R2.fq.gz clean.R2.fq.gz
    else
	java -Djava.io.tmpdir=./tmp -XX:ParallelGCThreads=$threads -jar $trimmomatic_dir/trimmomatic.jar SE -threads $threads -phred33  raw.fq.gz  trimmed.fq.gz ILLUMINACLIP:adapter.fa:2:30:10 SLIDINGWINDOW:5:20 MINLEN:36
	mv trimmed.fq.gz clean.fq.gz
    fi
else
    if [[ $mode == "PE" ]]
    then
	cp raw.R1.fq.gz clean.R1.fq.gz
	cp raw.R2.fq.gz clean.R2.fq.gz
    else
	cp raw.fq.gz clean.fq.gz
    fi
fi

# bwa mapping
for i in $(seq 1 1 $rounds_of_successive_polishing)
do
    $bwa_dir/bwa index refseq.tmp.fa

    if [[ $mode == "PE" ]]
    then
	$bwa_dir/bwa mem -t $threads -M refseq.tmp.fa clean.R1.fq.gz clean.R2.fq.gz >$prefix.round_${i}.sam
    else
	$bwa_dir/bwa mem -t $threads -M refseq.tmp.fa clean.fq.gz >$prefix.round_${i}.sam
    fi

    # index reference sequence
    $samtools_dir/samtools faidx refseq.tmp.fa
    java -Djava.io.tmpdir=./tmp -Dpicard.useLegacyParser=false -XX:ParallelGCThreads=$threads -jar $picard_dir/picard.jar CreateSequenceDictionary \
	-REFERENCE refseq.tmp.fa \
	-OUTPUT refseq.tmp.dict

    # sort bam file by picard-tools SortSam
    java -Djava.io.tmpdir=./tmp -Dpicard.useLegacyParser=false -XX:ParallelGCThreads=$threads -jar $picard_dir/picard.jar SortSam \
	-INPUT $prefix.round_${i}.sam \
	-OUTPUT $prefix.round_${i}.sort.bam \
	-SORT_ORDER coordinate

    # fixmate
    java -Djava.io.tmpdir=./tmp -Dpicard.useLegacyParser=false -XX:ParallelGCThreads=$threads -jar $picard_dir/picard.jar FixMateInformation \
	-INPUT $prefix.round_${i}.sort.bam \
	-OUTPUT $prefix.round_${i}.fixmate.bam

    # add or replace read groups and sort
    java -Djava.io.tmpdir=./tmp -Dpicard.useLegacyParser=false -XX:ParallelGCThreads=$threads -jar $picard_dir/picard.jar AddOrReplaceReadGroups \
	-INPUT $prefix.round_${i}.fixmate.bam \
	-OUTPUT $prefix.round_${i}.rdgrp.bam \
	-SORT_ORDER coordinate \
	-RGID $prefix \
	-RGLB $prefix \
	-RGPL "Illumina" \
	-RGPU $prefix \
	-RGSM $prefix \
	-RGCN "RGCN"

    # remove duplicates
    java -Djava.io.tmpdir=./tmp -Dpicard.useLegacyParser=false -XX:ParallelGCThreads=$threads -jar $picard_dir/picard.jar MarkDuplicates \
	-INPUT $prefix.round_${i}.rdgrp.bam \
	-REMOVE_DUPLICATES true \
	-METRICS_FILE $prefix.round_${i}.dedup.matrics \
	-OUTPUT $prefix.round_${i}.dedup.bam 

    # index the dedup.bam file
    $samtools_dir/samtools index $prefix.round_${i}.dedup.bam

    # GATK local realign
    # find realigner targets
    java -Djava.io.tmpdir=./tmp -XX:ParallelGCThreads=$threads -jar $gatk3_dir/GenomeAnalysisTK.jar \
	-R refseq.tmp.fa \
	-T RealignerTargetCreator \
	-I $prefix.round_${i}.dedup.bam \
	-o $prefix.round_${i}.realn.intervals
    # run realigner
    java -Djava.io.tmpdir=./tmp -XX:ParallelGCThreads=$threads -jar $gatk3_dir/GenomeAnalysisTK.jar \
	-R refseq.tmp.fa \
	-T IndelRealigner \
	-I $prefix.round_${i}.dedup.bam \
	-targetIntervals $prefix.round_${i}.realn.intervals \
	-o $prefix.round_${i}.realn.bam

    # index final bam file
    $samtools_dir/samtools index $prefix.round_${i}.realn.bam

    # for PE sequencing
    if [[ $mode == "PE" ]]
    then
	java -Djava.io.tmpdir=./tmp -Xmx16G -XX:ParallelGCThreads=$threads -jar $pilon_dir/pilon.jar \
	    --genome refseq.tmp.fa \
	    --frags $prefix.round_${i}.realn.bam \
	    --fix $fixlist \
	    --vcf \
	    --changes \
	    --output $prefix.assembly.pilon.round_${i} \
	    >$prefix.pilon.round_${i}.log
    else
	java -Djava.io.tmpdir=./tmp -Xmx16G -XX:ParallelGCThreads=$threads -jar $pilon_dir/pilon.jar \
	    --genome refseq.tmp.fa \
	    --unpaired $prefix.round_${i}.realn.bam \
	    --fix $fixlist \
	    --vcf \
	    --changes \
	    --output $prefix.assembly.pilon.round_${i} \
	    >$prefix.pilon.round_${i}.log
    fi
    perl $LRSDAY_HOME/scripts/summarize_pilon_correction.pl -i $prefix.assembly.pilon.round_${i}.changes
    gzip $prefix.assembly.pilon.round_${i}.vcf
    rm refseq.tmp.*
    rm $prefix.round_${i}.sam
    rm $prefix.round_${i}.sort.bam
    rm $prefix.round_${i}.fixmate.bam
    rm $prefix.round_${i}.rdgrp.bam
    rm $prefix.round_${i}.dedup.bam
    rm $prefix.round_${i}.dedup.matrics
    rm $prefix.round_${i}.dedup.bam.bai
    rm $prefix.round_${i}.realn.intervals
    cp $prefix.assembly.pilon.round_${i}.fasta refseq.tmp.fa
done

ln -s $prefix.assembly.pilon.round_${rounds_of_successive_polishing}.fasta $prefix.assembly.illumina_read_polished.fa
rm refseq.tmp.fa
rm -r tmp

# clean up intermediate files
if [[ $debug == "no" ]]
then
    rm adapter.fa
    rm *.fq.gz
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
