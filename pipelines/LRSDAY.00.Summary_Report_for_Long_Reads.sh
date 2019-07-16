#!/bin/bash
set -e -o pipefail
#######################################
# load environment variables for LRSDAY
source ./../../env.sh

#######################################
# set project-specific variables
long_reads_in_fastq="SK1.filtered_subreads.fastq.gz" # The fastq file of long-reads (in fastq or fastq.gz format). Default = "SK1.filtered_subreads.fastq.gz"
prefix="SK1" # The file name prefix for output files of the testing example. 
threads=1 # The number of threads to use. Default = "1".

#######################################
# process the pipeline

#source $nanoplot_dir/activate
$nanoplot_dir/NanoPlot \
    --threads $threads \
    --fastq $long_reads_in_fastq \
    --minlength 0 \
    --drop_outliers \
    --N50 \
    -o "${prefix}_Long_Reads_Summary_Report_out"

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
