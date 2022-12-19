#!/bin/bash
set -e -o pipefail
##########################################
# load environment variables for LRSDAY
source ./../../env.sh

###########################################
# set project-specific variables

prefix="CPG_1a"  # The file name prefix (only allowing strings of alphabetical letters, numbers, and underscores) for the processing sample. Default = "CPG_1a" for the testing example.         

input_assembly="./../02.Long-read-based_Assembly_Polishing/$prefix.assembly.long_read_polished.fa" # The file path of the input assembly before Illumina-based correction
trim_illumina_reads="yes" # Whether to trim the input Illumina reads. Use "yes" if prefer to perform trimming, otherwise use "no". Default = "yes". 
rounds_of_successive_polishing=3 # The number of total rounds of Illumina-read-based assembly polishing. Default = "3" for the testing example.
threads=8 # The number of threads to use. Default = "8".
mode="PE" # Illumina sequencing mode, "PE" for paired-end sequencing and "SE" for single-end sequencing. Default = "PE".
if [[ $mode == "PE" ]]
then
    reads_PE1="./../00.Short_Reads/$prefix.R1.fastq.gz" # Please replace the PE reads file name for your own project
    reads_PE2="./../00.Short_Reads/$prefix.R2.fastq.gz" # Please replace the PE reads file name for your own project
else
    reads_SE="./../00.Short_Reads/$prefix.fastq.gz" # Please replace the SE reads file name for your own project if you only have SE data
fi

fixlist="snps,indels" # The types of errors for Illumina-read-based correction by Pilon; see Pilon's manual for more details. Default = "snps,indels".
debug="no" # Whether to keep intermediate files for debugging. Use "yes" if prefer to keep intermediate files, otherwise use "no". Default = "no".

###########################################
# process the pipeline

if [[ $mode == "PE" ]]
then
    adapter="$trimmomatic_dir/adapters/TruSeq3-PE-2.fa" # adapter for PE reads
else
    adpater="$trimmomatic_dir/adapters/TruSeq3-SE.fa" # adapter for SE reads
fi

cp  $input_assembly refseq.tmp.fa
cp $adapter adapter.fa

mkdir tmp

if [[ $trim_illumina_reads == "yes" ]]
then
    if [[ $mode == "PE" ]]
    then
	$java_dir/java -Djava.io.tmpdir=./tmp -XX:ParallelGCThreads=$threads -jar $trimmomatic_dir/trimmomatic.jar PE -threads $threads -phred33  $reads_PE1 $reads_PE2  $prefix.trimmed.R1.fq.gz $prefix.trimmed.unpaired.R1.fq.gz $prefix.trimmed.R2.fq.gz $prefix.trimmed.unpaired.R2.fq.gz ILLUMINACLIP:adapter.fa:2:30:10  SLIDINGWINDOW:5:20 MINLEN:36
	rm $prefix.trimmed.unpaired.R1.fq.gz
	rm $prefix.trimmed.unpaired.R2.fq.gz 
	mv $prefix.trimmed.R1.fq.gz $prefix.clean.R1.fq.gz
	mv $prefix.trimmed.R2.fq.gz $prefix.clean.R2.fq.gz
    else
	$java_dir/java -Djava.io.tmpdir=./tmp -XX:ParallelGCThreads=$threads -jar $trimmomatic_dir/trimmomatic.jar SE -threads $threads -phred33 $reads_SE  $prefix.trimmed.fq.gz ILLUMINACLIP:adapter.fa:2:30:10 SLIDINGWINDOW:5:20 MINLEN:36
	mv $prefix.trimmed.fq.gz $prefix.clean.fq.gz
    fi
else
    if [[ $mode == "PE" ]]
    then
	cp $prefix.raw.R1.fq.gz $prefix.clean.R1.fq.gz
	cp $prefix.raw.R2.fq.gz $prefix.clean.R2.fq.gz
    else
	cp $prefix.raw.fq.gz $prefix.clean.fq.gz
    fi
fi

# bwa mapping
for i in $(seq 1 1 $rounds_of_successive_polishing)
do
    $bwa_dir/bwa index refseq.tmp.fa

    if [[ $mode == "PE" ]]
    then
	$bwa_dir/bwa mem -t $threads -M refseq.tmp.fa $prefix.clean.R1.fq.gz $prefix.clean.R2.fq.gz >$prefix.round_${i}.sam
    else
	$bwa_dir/bwa mem -t $threads -M refseq.tmp.fa $prefix.clean.fq.gz >$prefix.round_${i}.sam
    fi

    # index reference sequence
    $samtools_dir/samtools faidx refseq.tmp.fa
    $java_dir/java -Djava.io.tmpdir=./tmp -Dpicard.useLegacyParser=false -XX:ParallelGCThreads=$threads -jar $picard_dir/picard.jar CreateSequenceDictionary \
	-REFERENCE refseq.tmp.fa \
	-OUTPUT refseq.tmp.dict

    # sort bam file by picard-tools SortSam
    $java_dir/java -Djava.io.tmpdir=./tmp -Dpicard.useLegacyParser=false -XX:ParallelGCThreads=$threads -jar $picard_dir/picard.jar SortSam \
	-INPUT $prefix.round_${i}.sam \
	-OUTPUT $prefix.round_${i}.sort.bam \
	-SORT_ORDER coordinate

    # fixmate
    $java_dir/java -Djava.io.tmpdir=./tmp -Dpicard.useLegacyParser=false -XX:ParallelGCThreads=$threads -jar $picard_dir/picard.jar FixMateInformation \
	-INPUT $prefix.round_${i}.sort.bam \
	-OUTPUT $prefix.round_${i}.fixmate.bam

    # add or replace read groups and sort
    $java_dir/java -Djava.io.tmpdir=./tmp -Dpicard.useLegacyParser=false -XX:ParallelGCThreads=$threads -jar $picard_dir/picard.jar AddOrReplaceReadGroups \
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
    $java_dir/java -Djava.io.tmpdir=./tmp -Dpicard.useLegacyParser=false -XX:ParallelGCThreads=$threads -jar $picard_dir/picard.jar MarkDuplicates \
	-INPUT $prefix.round_${i}.rdgrp.bam \
	-REMOVE_DUPLICATES true \
	-METRICS_FILE $prefix.round_${i}.dedup.matrics \
	-OUTPUT $prefix.round_${i}.dedup.bam 

    # index the dedup.bam file
    $samtools_dir/samtools index $prefix.round_${i}.dedup.bam

    # GATK local realign
    # find realigner targets
    $java_dir/java -Djava.io.tmpdir=./tmp -XX:ParallelGCThreads=$threads -jar $gatk3_dir/GenomeAnalysisTK.jar \
	-R refseq.tmp.fa \
	-T RealignerTargetCreator \
	-I $prefix.round_${i}.dedup.bam \
	-o $prefix.round_${i}.realn.intervals
    # run realigner
    $java_dir/java -Djava.io.tmpdir=./tmp -XX:ParallelGCThreads=$threads -jar $gatk3_dir/GenomeAnalysisTK.jar \
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
	$java_dir/java -Djava.io.tmpdir=./tmp -Xmx16G -XX:ParallelGCThreads=$threads -jar $pilon_dir/pilon.jar \
	    --genome refseq.tmp.fa \
	    --frags $prefix.round_${i}.realn.bam \
	    --fix $fixlist \
	    --vcf \
	    --changes \
	    --output $prefix.assembly.pilon.round_${i} \
	    >$prefix.pilon.round_${i}.log
    else
	$java_dir/java -Djava.io.tmpdir=./tmp -Xmx16G -XX:ParallelGCThreads=$threads -jar $pilon_dir/pilon.jar \
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

ln -s $prefix.assembly.pilon.round_${rounds_of_successive_polishing}.fasta $prefix.assembly.short_read_polished.fa
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
    echo "#########################################################"
    echo ""
    echo "LRSDAY message: This bash script has been successfully processed! :)"
    echo ""
    echo "#########################################################"
    echo ""
    exit 0
fi
############################
