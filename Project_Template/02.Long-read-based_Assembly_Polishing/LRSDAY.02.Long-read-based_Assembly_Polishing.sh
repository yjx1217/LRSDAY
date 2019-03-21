#!/bin/bash
set -e -o pipefail
##########################################
# load environment variables for LRSDAY
source ./../../env.sh

###########################################
# set project-specific variables

input_assembly="./../01.Long-read-based_Genome_Assembly/SK1.assembly.raw.fa" # The file path of the input raw long-read-based assembly for polishing.
long_reads_in_fastq="./../00.Long_Reads/SK1.filtered_subreads.fastq.gz" # The file path of the long-read fastq file. 
long_read_technology="pacbio" # The used long-read sequencing technology. Use "pacbio" or "nanopore". Default = "pacbio" for the testing example.

### When long_read_technology="pacbio" ####
pacbio_bam_fofn_file="./../00.Long_Reads/pacbio_fofn_files/SK1.merged.bam.fofn" # The file path to the fofn file containing the absolute path to the PacBio bam files. BAM file is the native output format for PacBio Sequel platform but this is not the case for the RSII platform. For RSII data, the bax2bam file conversion is needed. This can be done by running the LRSDAY.00.Retrieve_Sample_PacBio_Reads.sh script in the 00.Long_Reads directory.
pacbio_reads_type="RSII" # The sequencing machine used to generate the input PacBio reads . Use "RSII" or "Sequel". Default = "RSII" for the testing example. 

### When long_read_technology="nanopore" ###
nanopore_fast5_files="./../00.Long_Reads/nanopore_fast5_files" # The file path to the directory containing raw Oxford Nanopore FAST5 files.
nanopore_basecalling_sequencing_summary="./../00.Long_Reads/nanopore_fast5_files/sequencing_summary.txt" # The file path to the nanopore albacore/guppy basecaller sequencing summary output. This summary file is not necessary but it can help the polishing step to run much faster when available. When this file is unavailable, set nanopore_albacore_sequencing_summary="".

prefix="SK1" # The file name prefix for the output files. Default = "SK1" for the testing example.

threads=1  # The number of threads to use. Default = "1".
ploidy=1 # The ploidy status of the sequenced genome. use "1" for haploid genome and "2" for diploid genome. Default = "1" for the testing example.
rounds_of_successive_polishing=1 # The number of total rounds of long-read-based assembly polishing. Default = "1" for the testing example.
debug="no" # Use "yes" if prefer to keep intermediate files, otherwise use "no". Default = "no"

###########################################
# process the pipeline

ln -s $input_assembly $prefix.assembly.raw.fa
cp $prefix.assembly.raw.fa $prefix.assembly.tmp.fa

mkdir tmp

if [[ $long_read_technology == "pacbio" ]]
then
    # perform correction using PacBio's pbalign-arrow pipeline
    source $miniconda2_dir/activate $conda_pacbio_dir/../../conda_pacbio_env
    if [[ $pacbio_reads_type == "RSII" ]]
    then
	for i in $(seq 1 1 $rounds_of_successive_polishing)
	do
	    $samtools_dir/samtools faidx $prefix.assembly.tmp.fa
	    $conda_pacbio_dir/pbalign --nproc $threads --algorithm blasr --tmpDir ./tmp $pacbio_bam_fofn_file $prefix.assembly.tmp.fa $prefix.pbalign.round_${i}.bam
	    if [[ $ploidy == "1" ]]
	    then
		$conda_pacbio_dir/variantCaller --algorithm=quiver -x 5 -X 120 -q 20 -v -j $threads $prefix.pbalign.round_${i}.bam -r $prefix.assembly.tmp.fa -o $prefix.assembly.consensus.round_${i}.fa -o $prefix.assembly.consensus.round_${i}.fq -o $prefix.assembly.consensus.round_${i}.vcf
	    else
		$conda_pacbio_dir/variantCaller --algorithm=quiver -x 5 -X 120 -q 20 -v -j $threads $prefix.pbalign.round_${i}.bam -r $prefix.assembly.tmp.fa -o $prefix.assembly.consensus.round_${i}.fa -o $prefix.assembly.consensus.round_${i}.fq -o $prefix.assembly.consensus.round_${i}.vcf --diploid 
	    fi
	    rm $prefix.assembly.tmp.fa
            rm $prefix.assembly.tmp.fa.fai
	    cp $prefix.assembly.consensus.round_${i}.fa $prefix.assembly.tmp.fa
	    rm $prefix.assembly.consensus.round_${i}.fq
	    rm $prefix.assembly.consensus.round_${i}.vcf
	done
    else
	for i in $(seq 1 1 $rounds_of_successive_polishing)
        do
	    $samtools_dir/samtools faidx $prefix.assembly.tmp.fa
	    $conda_pacbio_dir/pbalign --nproc $threads --algorithm blasr --tmpDir ./tmp $pacbio_bam_fofn_file $prefix.assembly.tmp.fa $prefix.pbalign.round_${i}.bam
	    if [[ $ploidy == "1" ]]
	    then
		$conda_pacbio_dir/variantCaller --algorithm=arrow -x 5 -X 120 -q 20 -v -j $threads $prefix.pbalign.round_${i}.bam -r $prefix.assembly.tmp.fa -o $prefix.assembly.consensus.round_${i}.fa -o $prefix.assembly.consensus.round_${i}.fq -o $prefix.assembly.consensus.round_${i}.vcf
	    else
		$conda_pacbio_dir/variantCaller --algorithm=arrow -x 5 -X 120 -q 20 -v -j $threads $prefix.pbalign.round_${i}.bam -r $prefix.assembly.tmp.fa -o $prefix.assembly.consensus.round_${i}.fa -o $prefix.assembly.consensus.round_${i}.fq -o $prefix.assembly.consensus.round_${i}.vcf --diploid 
	    fi
	    rm $prefix.assembly.tmp.fa
            rm $prefix.assembly.tmp.fa.fai
	    cp $prefix.assembly.consensus.round_${i}.fa $prefix.assembly.tmp.fa
	    gzip $prefix.assembly.consensus.round_${i}.fq
	    gzip $prefix.assembly.consensus.round_${i}.vcf
	done
    fi
    ln -s $prefix.assembly.consensus.round_${rounds_of_successive_polishing}.fa $prefix.assembly.long_read_polished.fa
    rm $prefix.assembly.tmp.fa
    source $miniconda2_dir/deactivate
else
    # perform correction using the minimap2-nanopolish pipeline
    source $nanopolish_dir/py3_virtualenv_nanopolish/bin/activate
    if [[ -z "$nanopore_basecalling_sequencing_summary" ]]
    then
	$nanopolish_dir/nanopolish index -d $nanopore_fast5_files $long_reads_in_fastq
    else
	$nanopolish_dir/nanopolish index -d $nanopore_fast5_files -s $nanopore_basecalling_sequencing_summary $long_reads_in_fastq
    fi
    for i in $(seq 1 1 $rounds_of_successive_polishing)
    do
	java -Djava.io.tmpdir=./tmp -Dpicard.useLegacyParser=false -XX:ParallelGCThreads=$threads -jar $picard_dir/picard.jar CreateSequenceDictionary -REFERENCE $prefix.assembly.tmp.fa -OUTPUT $prefix.assembly.tmp.dict
	$minimap2_dir/minimap2 -ax map-ont $prefix.assembly.tmp.fa $long_reads_in_fastq > $prefix.minimap2.round_${i}.sam
	java -Djava.io.tmpdir=./tmp -Dpicard.useLegacyParser=false -XX:ParallelGCThreads=$threads -jar $picard_dir/picard.jar SortSam -INPUT $prefix.minimap2.round_${i}.sam -OUTPUT $prefix.minimap2.round_${i}.bam -SORT_ORDER coordinate
	$samtools_dir/samtools index $prefix.minimap2.round_${i}.bam
	rm $prefix.minimap2.round_${i}.sam
	python3 $nanopolish_dir/scripts/nanopolish_makerange.py $prefix.assembly.tmp.fa | $parallel_dir/parallel --results ${prefix}_nanopolish_round_${i}_results -P 1 \
     	$nanopolish_dir/nanopolish variants --consensus -o $prefix.polished.{1}.vcf -w {1} --ploidy $ploidy -r $long_reads_in_fastq -b $prefix.minimap2.round_${i}.bam -g $prefix.assembly.tmp.fa -t $threads --min-candidate-frequency 0.2  || true 
	$nanopolish_dir/nanopolish vcf2fasta -g $prefix.assembly.tmp.fa $prefix.polished.*.vcf > $prefix.assembly.nanopolish.round_${i}.fa
	rm $prefix.assembly.tmp.fa
	rm $prefix.assembly.tmp.dict
	cp $prefix.assembly.nanopolish.round_${i}.fa $prefix.assembly.tmp.fa
	mv $prefix.polished.*.vcf ${prefix}_nanopolish_round_${i}_results
    done
    ln -s $prefix.assembly.nanopolish.round_${rounds_of_successive_polishing}.fa $prefix.assembly.long_read_polished.fa
    rm $prefix.assembly.tmp.fa
fi

rm -r tmp

# clean up intermediate files
if [[ $debug == "no" ]]
then
    echo "clean up"
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
