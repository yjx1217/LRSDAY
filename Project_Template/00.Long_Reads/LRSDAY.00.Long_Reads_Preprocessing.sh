#!/bin/bash
set -e -o pipefail
#######################################
# load environment variables
source ./../../env.sh

#######################################
# set project-specific variables

prefix="CPG_1a" # The file name prefix (only allowing strings of alphabetical letters, numbers, and underscores) for the processing sample. Default = "CPG_1a" for the testing example.
reads="./../00.Long_Reads/nanopore_basecalled_fastq_files/$prefix/$prefix.basecalled_reads.Q5.pass.fastq.gz" # The file path of the long reads file (in fastq or fastq.gz format).
reads_type="nanopore-raw" # The long reads data type: "pacbio-raw" or "pacbio-corrected" or "nanopore-raw" or "nanopore-corrected".
run_filtering="yes" # Whether to filter and downsample the reads: "yes" or "no". Default = "yes".
genome_size="12500000" # The haploid genome size (in bp) of sequenced organism. Default = "12500000" (i.e. 12.5 Mb for the budding yeast S. cereviaie genome). This is used to calculate targeted sequencing coverage after read filtering (see below). 
post_filtering_coverage="60" # Targeted sequencing coverage after read filtering and downsampling. Default = "60" (i.e. 60x coverage).
threads=4 # The number of threads to use. Default = "4".

#######################################
# process the pipeline

if [[ "$reads_type" == "nanopore-raw" || "$reads_type" == "nanopore-corrected" ]]
then
    $porechop_dir/porechop -i $reads -o $prefix.porechop.fastq.gz --discard_middle --threads $threads > $prefix.porechop.summary.txt
    if [[ "$run_filtering" == "yes" ]]
    then
	filtlong_target_bases=$(($genome_size * $post_filtering_coverage))
	echo ""
	echo "genome_size=$genome_size, post_filtering_coverage=$post_filtering_coverage, filtlong_target_bases=$filtlong_target_bases"
	echo ""
	$filtlong_dir/filtlong --min_length 1000 --mean_q_weight 10 --window_q_weight 5 --target_bases $filtlong_target_bases $prefix.porechop.fastq.gz | gzip > $prefix.filtlong.fastq.gz
    fi
else
    if [[ "$run_filtering" == "yes" ]]
    then
	filtlong_target_bases=$(($genome_size * $post_filtering_coverage))
	echo ""
	echo "genome_size=$genome_size, post_filtering_coverage=$post_filtering_coverage, filtlong_target_bases=$filtlong_target_bases"
	echo ""
	$filtlong_dir/filtlong --min_length 1000 --mean_q_weight 10 --window_q_weight 5 --target_bases $filtlong_target_bases $reads | gzip > $prefix.filtlong.fastq.gz
    fi
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
