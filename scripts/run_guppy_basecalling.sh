#!/bin/bash
# last update: 2019.03.07

source ./../env.sh
guppy_dir="$LRSDAY_HOME/build/ont-guppy-cpu/bin"

flowcell_version="FLO-MIN106" # Please check your own flowcell version
sequencing_kit_version="SQK-LSK108" # Please check your own flowcell_kit_version
barcode_kit_version="EXP-NBD104" # Please check your own barcode_kit_version

project_name="Project_YGL3210" # LRSDAY Project name
raw_reads_directory="$LRSDAY_HOME/$project_name/00.Long_Reads/RAW" # The directory containing the raw nanopore reads before basecalling
basecalling_output_directory="$LRSDAY_HOME/$project_name/00.Long_Reads/Guppy_Basecalling_out"
demultiplexing_output_directory="$LRSDAY_HOME/$project_name/00.Long_Reads/Guppy_Demultiplexing_out"

qual=5 # read quality filter for guppy basecalling
num_callers=1 # num_callers for guppy
threads_per_caller=1 # threads_per_caller for guppy
demultiplexing_threads=1 # threads to use for demultiplexing

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

$guppy_dir/guppy_barcoder \
    --barcode_kit $barcode_kit_version \
    --recursive \
    --input_path $basecalling_output_directory/pass \
    --save_path $demultiplexing_output_directory \
    --worker_threads $demultiplexing_threads


cd $demultiplexing_output_directory
for b in barcode*
do 
    echo "barcode=$b"
    cat ./$b/*.fastq |gzip -c $project_name.basecalled_reads.Q${qual}.pass.$b.fastq.gz
done

cat ./unclassified/*.fastq |gzip -c $project_name.basecalled_reads.Q${qual}.pass.$b.fastq.gz





