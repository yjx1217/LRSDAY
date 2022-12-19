#!/bin/bash
set -e -o pipefail
#######################################
# load environment variables
source ./../../env.sh

########################
guppy_run_mode="cpu" # The running mode of Guppy basecalling: "gpu" or "cpu". Default = "cpu".
if [[ $guppy_run_mode == "gpu" ]]
then
    # set CUDA environment for GPU
    gpu_bin_path="/public/software/cuda-11.4/bin" # The path of CUDA's bin directory. Default = "/public/software/cuda-11.4/bin".
    gpu_lib_path="/public/software/cuda-11.4/lib64" # The path of CUDA's lib directory. Default = "/public/software/cuda-11.4/lib64".
    gpu_include_path="/public/software/cuda-11.4/include" # The path of CUDA's include directory. Default = "/public/software/cuda-11.4/include".
fi
##########################3
export PATH=$gpu_bin_path:$PATH
export LD_LIBRARY_PATH=$gpu_lib_path:$LD_LIBRARY_PATH
export C_INCLUDE_PATH=$gpu_include_path:$C_INCLUDE_PATH
export CPLUS_INCLUDE_PATH=$C_INCLUDE_PATH

##k#####################################
# set project-specific variables

sample_id="CPG_1a" # The file name prefix (only allowing strings of alphabetical letters, numbers, and underscores) for the processing sample. Default = "CPG_1a" for the testing example.              
raw_fast5_dir="./nanopore_raw_fast5_files/$sample_id" # The directory containing the raw nanopore reads before basecalling. Default = "./raw_fast5/$sample_id".
flowcell_version="FLO-MIN106" # The flowcell version of the nanopore run. Default = "FLO-MIN106".
sequencing_kit_version="SQK-LSK108" # The sequencing kit version of the nanopore run. Default = "SQK-LSK109". 
barcode_kit_version="EXP-NBD104" # The barcode kit version of the nanopore run if barcoding has been introduced during library preparation. Default = "EXP-NBD104". For the testing example, we have not applied barcoding. So this information will not be used when we set run_demultiplexing="no" below.

run_basecalling="yes" # Whether to perform basecalling: "yes" or "no". Default = "yes". 
run_demultiplexing="no" # Whether to perform demultiplexing: "yes" or "no". Default = "no". For the testing example, we have not applied barcoding. So no need to run this step. 
run_nanoplotting="yes" # Whether to perform nanoplotting: "yes" or "no". Default = "yes". 

threads=16 # The number of CPU threads to use. Default = 16.
debug="no" # Whether to keep intermediate files for debugging: "yes" or "no". Default = "no".

#############################
# Normally no need to change the following parameter settings.
basecalled_fast5_dir="./nanopore_basecalled_fast5_files/$sample_id" # The directory containing the basecalled fast5 reads. This directory will be automatically generated when running basecalling.
basecalled_fastq_dir="./nanopore_basecalled_fastq_files/$sample_id" # The directory containing the basecalled fastq reads. This directory will be automatically generated when running basecalling.
basecalled_qc_dir="./nanopore_basecalled_qc_files/$sample_id" # The directory containing nanoplot QC outputs. This directory will be automatically generated when running nanoplot.
trim_strategy="dna" # Trimming strategy to apply: 'dna' or 'rna' or 'none' (to disable trimming). Default = "rna".
gpu_device="auto" # Which GPU device to use. Default = "auto". 
qual=5 # read quality filter for guppy basecalling. Default = 5.
num_callers_in_cpu_mode=$threads # The number of callers for guppy basecalling in the CPU mode. Default = "$threads".
num_callers_in_gpu_mode=$threads # The number of callers for guppy basecalling in the GPU mode. Default = "1".
gpu_runners_per_device=1 # The number of GPU runners per device for guppy basecalling in the GPU mode. Default = "1".
threads_per_caller=1 # The number of threads per caller for guppy basecalling. Default = "1".
demultiplexing_threads=$threads # The number of threads to use for guppy demultiplexing. Default = "$threads".
demultiplexed_fastq_dir="$basecalled_fastq_dir/demultiplexed_fastq" # The directory containing the demultiplexed basecalled nanopore reads. This directory will be automatically generated when running demultiplexing. 

##############################

wkdir=$(pwd)
if [[ "$run_demultiplexing" == "yes" && "barcode_kit_version" == "" ]]
then
    echo "The variable run_demultiplexing has been set to \"yes\" but the barcode_kit_version variable has not been specified!"
    echo "Please specified barcode_kit_version if you want to run demultiplexing."
    exit;
fi


if [[ "$run_basecalling" == "yes" ]]
then
    echo "Check if $basecalled_fast5_dir is empty for running basecalling."
    if [[ -d $basecalled_fast5_dir && "$(ls $basecalled_fast5_dir)" ]]
    then
	echo "Warning! The basecalled fast5 directory $basecalled_fast5_dir exists and it is not empty! Please empty its content if you want to run basecalling."
	echo "Exit!!!"
	exit
    elif [[ -d $basecalled_fastq_dir && "$(ls $basecalled_fastq_dir)" ]]
    then
        echo "Warning! The basecalled fastq directory $basecalled_fastq_dir exists and it is not empty! Please empty its content if you want to run basecalling."
        echo "Exit!!!"
        exit
    else
	echo "Check passed!"
	echo "Running basecalling .."
	if [[ ! -d $basecalled_fast5_dir ]]
	then
	    mkdir -p $basecalled_fast5_dir
	fi
	if [[ ! -d $basecalled_fastq_dir ]]
	then
	    mkdir -p $basecalled_fastq_dir
	fi

	if [[ "$guppy_run_mode" == "gpu" ]]
	then
	    $guppy_gpu_dir/guppy_basecaller \
		--flowcell $flowcell_version \
		--kit $sequencing_kit_version \
		--recursive \
		--trim_strategy $trim_strategy \
		--input_path $raw_fast5_dir \
		--save_path $basecalled_fast5_dir \
		--fast5_out \
		--min_qscore $qual \
		--device $gpu_device \
		--num_callers $num_callers_in_gpu_mode \
		--gpu_runners_per_device $gpu_runners_per_device \
		--compress_fastq	    
	else
	    $guppy_cpu_dir/guppy_basecaller \
		--flowcell $flowcell_version \
		--kit $sequencing_kit_version \
		--recursive \
		--trim_strategy $trim_strategy \
		--input_path $raw_fast5_dir \
		--save_path $basecalled_fast5_dir \
		--fast5_out \
		--min_qscore $qual \
		--num_callers $num_callers_in_cpu_mode \
		--cpu_threads_per_caller $threads_per_caller \
		--compress_fastq	    
	fi
	cat $basecalled_fast5_dir/pass/*.fastq.gz  > $basecalled_fastq_dir/$sample_id.basecalled_reads.Q${qual}.pass.fastq.gz
	# cat $basecalled_fast5_dir/fail/*.fastq.gz  > $basecalled_fastq_dir/$sample_id.basecalled_reads.Q${qual}.fail.fastq.gz
	if [[ $debug != "yes" ]]
	then
	    rm -r $basecalled_fast5_dir/pass
	    rm -r $basecalled_fast5_dir/fail
	fi
    fi
fi

cd $wkdir
if [[ "$run_demultiplexing" == "yes" ]]
then
    echo "Check if $basecalled_fastq_dir has basecalled reads for running demultiplexing."
    if [[ "$(ls $basecalled_fastq_dir)" ]]
    then
	echo "Running demultiplexing."
	if [[ "$guppy_run_mode" == "gpu" ]]
	then
	    $guppy_gpu_dir/guppy_barcoder \
		--barcode_kits $barcode_kit_version \
		--recursive \
		--input_path $basecalled_fastq_dir \
		--save_path $demultiplexed_fastq_dir \
		--worker_threads $demultiplexing_threads
	else
	    $guppy_cpu_dir/guppy_barcoder \
		--barcode_kits $barcode_kit_version \
		--recursive \
		--input_path $basecalled_fastq_dir \
		--save_path $demultiplexed_fastq_dir \
		--worker_threads $demultiplexing_threads
	fi

	cd $demultiplexed_fastq_dir
	for b in barcode*
	do 
	    echo "for demultiplexing: barcode=$b"
	    cat ./$b/*.fastq |gzip -c > $sample_id.basecalled_reads.Q${qual}.pass.$b.fastq.gz
	done
	cat ./unclassified/*.fastq |gzip -c > $sample_id.basecalled_reads.Q${qual}.pass.unclassified.fastq.gz
    else
	echo "There is no reads in $basecalled_fastq_dir!"
	echo "Please put the basecalled reads in $basecalled_fastq_dir for demultiplexing!"
	echo "Exit!!!"
	exit
    fi
fi

set +oe pipefail 

cd $wkdir
if [[ "$run_nanoplotting" == "yes" ]]
then
    echo "Check if $basecalled_fastq_dir has basecalled reads for running nanoplotting."
    if [[ "$(ls $basecalled_fastq_dir)" ]]
    then
	echo "Running nanoplotting."
	mkdir -p $basecalled_qc_dir
	fastq_input="$basecalled_fastq_dir/$sample_id.basecalled_reads.Q${qual}.pass.fastq.gz"
	$nanoplot_dir/NanoPlot \
	    --threads $threads \
	    --fastq $fastq_input \
	    --N50 \
	    -o "$basecalled_qc_dir/${sample_id}_basecalled_reads_Q${qual}_pass_NanoPlot_out"
    fi
    if [[ "$run_demultiplexing" == "yes" ]]
    then
	cd $demultiplexed_fastq_dir
	for b in barcode*
	do
	    echo "for nanoplotting: barcode=$b"
	    fastq_input="./$sample_id.basecalled_reads.Q${qual}.pass.$b.fastq.gz"
	    $nanoplot_dir/NanoPlot \
		--threads $threads \
		--fastq $fastq_input \
		--N50 \
		-o "./../../../$basecalled_qc_dir/${sample_id}_basecalled_reads_Q${qual}_pass_${b}_NanoPlot_out"
	done
	echo "for nanoplotting: unclassified"
	fastq_input="./$sample_id.basecalled_reads.Q${qual}.pass.unclassified.fastq.gz"
	$nanoplot_dir/NanoPlot \
            --threads $threads \
            --fastq $fastq_input \
            --N50 \
            -o "./../../../$basecalled_qc_dir/${sample_id}_basecalled_reads_Q${qual}_pass_unclassified_NanoPlot_out" 
    fi
fi



############################
# checking bash exit status
if [[ $? -eq 0 ]]
then
    echo "" 
    echo "##########################################################################"
    echo "" 
    echo "LRSDAY message: This bash script has been successfully processed! :)"
    echo ""
    echo "##########################################################################"
    echo ""
    exit 0
fi
############################
