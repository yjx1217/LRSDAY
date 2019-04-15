#!/bin/bash
set -e -o pipefail
#######################################
# load environment variables for LRSDAY
source ./../../env.sh

#######################################
# set project-specific variables
project_name="Project_Example" # LRSDAY Project name. Default = "Project_Example".
run_basecalling="yes" # Whether to perform basecalling: "yes" or "no". Default = "yes". 
run_demultiplexing="yes" # Whether to perform demultiplexing: "yes" or "no". Default = "yes". 
run_nanoplotting="yes" # Whether to perform nanoplotting: "yes" or "no". Default = "yes". 

flowcell_id="FAKXXXXX" # The flowcell ID of the nanopore run. Default = "FAKXXXXX".
flowcell_version="FLO-MIN106" # The flowcell version of the nanopore run. Default = "FLO-MIN106".
sequencing_kit_version="SQK-LSK108" # The sequencing kit version of the nanopore run. Default = "SQK-LSK108".
barcode_kit_version="EXP-NBD103" # The barcode kit version of the nanopore run. Default = "EXP-NBD103".

raw_reads_directory="$LRSDAY_HOME/$project_name/00.Long_Reads/nanopore_raw_fast5_files" # The directory containing the raw nanopore reads before basecalling
basecalling_output_directory="$LRSDAY_HOME/$project_name/00.Long_Reads/nanopore_basecalled_fast5_files" # The directory containing the basecalled nanopore reads. This directory will be automatically generated when running basecalling.
threads=8 # The number of threads to use. Default = 8.

#############################
# normally no need to change the following
qual=5 # read quality filter for guppy basecalling
num_callers=$threads # num_callers for guppy
threads_per_caller=1 # threads_per_caller for guppy
demultiplexing_threads=$threads # threads to use for demultiplexing
demultiplexing_output_directory="$LRSDAY_HOME/$project_name/00.Long_Reads/nanopore_demultiplexed_fastq_files" # The directory containing the demultiplexed basecalled nanopore reads. This directory will be automatically generated when running demultiplexing. 

if [[ "$run_basecalling" == "yes" ]]
then
    echo "Check if $basecalling_output_directory is empty for running basecalling."
    if [[ "$(ls $basecalling_output_directory)" ]]
    then
	echo "Warning! The basecalling directory is not empty! Please empty its content if you want to run basecalling."
	echo "Exit!!!"
	exit
    else
	echo "Running basecalling."
	$guppy_dir/guppy_basecaller \
	    --flowcell $flowcell_version \
	    --kit $sequencing_kit_version \
	    --recursive \
	    --input_path $raw_reads_directory \
	    --save_path $basecalling_output_directory \
	    --fast5_out \
	    --qscore_filtering \
	    --min_qscore $qual \
	    --num_callers $num_callers \
	    --cpu_threads_per_caller $threads_per_caller 
	cd $basecalling_output_directory
	cat ./pass/*.fastq |gzip -c > $project_name.basecalled_reads.Q${qual}.pass.fastq.gz
	cat ./fail/*.fastq |gzip -c > $project_name.basecalled_reads.Q${qual}.fail.fastq.gz
    fi
fi

if [[ "$run_demultiplexing" == "yes" ]]
then
    echo "Check if $basecalling_output_directory/pass has basecalled reads for running demultiplexing."
    if [[ "$(ls $basecalling_output_directory/pass)" ]]
    then
	echo "Running demultiplexing."
	$guppy_dir/guppy_barcoder \
	    --barcode_kit $barcode_kit_version \
	    --recursive \
	    --input_path $basecalling_output_directory/pass \
	    --save_path $demultiplexing_output_directory \
	    --worker_threads $demultiplexing_threads
	
	cd $demultiplexing_output_directory
	for b in barcode*
	do 
	    echo "for demultiplexing: barcode=$b"
	    cat ./$b/*.fastq |gzip -c > $project_name.basecalled_reads.Q${qual}.pass.$b.fastq.gz
	done
	cat ./unclassified/*.fastq |gzip -c > $project_name.basecalled_reads.Q${qual}.pass.unclassified.fastq.gz
    else
	echo "There is no reads in $basecalling_output_directory/pass!"
	echo "Please put the basecalled reads in $basecalling_output_directory/pass for demultiplexing!"
	echo "Exit!!!"
	exit
    fi
fi

set +oe pipefail 

if [[ "$run_nanoplotting" == "yes" ]]
then
    echo "Check if $basecalling_output_directory/pass has basecalled reads for running nanoplotting."
    if [[ "$(ls $basecalling_output_directory/pass)" ]]
    then
	echo "Running nanoplotting."
	cd $basecalling_output_directory
	fastq_input="$project_name.basecalled_reads.Q${qual}.pass.fastq.gz"
	source $nanoplot_dir/activate
	$nanoplot_dir/NanoPlot \
	    --threads $threads \
	    --fastq $fastq_input \
	    --N50 \
	    -o "${project_name}_Q${qual}_pass_NanoPlot_out"
    fi
    if [[ "$run_demultiplexing" == "yes" ]]
    then
	cd $demultiplexing_output_directory
	for b in barcode*
	do
	    echo "for nanoplotting: barcode=$b"
	    fastq_input="$project_name.basecalled_reads.Q${qual}.pass.$b.fastq.gz"
	    source $nanoplot_dir/activate
	    $nanoplot_dir/NanoPlot \
		--threads $threads \
		--fastq $fastq_input \
		--N50 \
		-o "${project_name}_Q${qual}_pass_${b}_NanoPlot_out"
	done
	echo "for nanoplotting: unclassified"
	fastq_input="$project_name.basecalled_reads.Q${qual}.pass.unclassified.fastq.gz"
	$nanoplot_dir/NanoPlot \
            --threads $threads \
            --fastq $fastq_input \
            --N50 \
            -o "${project_name}_Q${qual}_pass_unclassified_NanoPlot_out" 
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
