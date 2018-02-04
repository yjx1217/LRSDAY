#!/bin/bash
set -e -o pipefail
#######################################
# load environment variables for LRSDAY
source ./../../env.sh

#######################################
# set project-specific variables
file_name="SK1.filtered_subreads.bam" # the name of the ENA bam file
file_url="ftp://ftp.sra.ebi.ac.uk/vol1/ERZ448/ERZ448251/SK1.filtered_subreads.bam" # the URL of the ENA bam file
prefix="SK1" # file name prefix for output files

#######################################
# process the pipeline

echo "download the bam file from the ENA database"
wget $file_url
echo "bam2fastq ..."
$bedtools_dir/bedtools bamtofastq -i $file_name -fq $prefix.filtered_subreads.fastq 
echo "gzip fastq ..."
gzip $prefix.filtered_subreads.fastq
rm $file_name

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
