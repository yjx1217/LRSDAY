#!/bin/bash
set -e -o pipefail
##########################################
# load environment variables for LRSDAY
source ./../../env.sh
PATH=$gnuplot_dir:$PATH

###########################################
# set project-specific variables
prefix="SK1" # file name prefix for the output files
reads="$LRSDAY_HOME/Project_Example/00.Long_Reads/SK1.filtered_subreads.fastq.gz" # full path of the long reads file (in fastq or fastq.gz format)
reads_type="pacbio-raw" # long reads data type: "pacbio-raw" or "pacbio-corrected" or "nanopore-raw" or "nanopore-corrected"
genome_size="12.5m" # estimated genome size with the format of <number>[g|m|k], e.g. 12.5m for 12.5 Mb
threads=1 # number of threads to use
vcf="yes" # use "yes" if prefer to have vcf file generated to show SNP and INDEL differences between the assembled genome and the reference genome.
dotplot="yes" # use "yes" if prefer to plot genome-wide dotplot based on the comparison with the reference genome below, otherwise use "no"
ref_genome_raw="./../00.Ref_Genome/S288C.ASM205763v1.fa" # path of the raw reference genome, only needed when dotplot="yes" or vcf="yes".
debug="no" # use "yes" if prefer to keep intermediate files, otherwise use "no".

###########################################
# process the pipeline

# check project-specific variables
if [[ $vcf == "yes" || $dotplot == "yes" ]]
then
    # check if ref_genome_raw is defined
    if [[ -z $ref_genome_raw ]]
    then
	echo "The vcf and doptlot outputs require the variable ref_genome_raw be defined!"
	echo "Please define this variable in the script; Please delete all the old output files and directories; and re-run this step!"
	echo "Script exit!"
	echo ""
	exit
    elif [[ ! -f $ref_genome_raw ]]
    then
	echo "The vcf and doptlot outputs require the $ref_genome_raw as defined with the variable \"ref_genome raw\" but this file cannot be found!"
	echo "Please make sure that this file truly exists; Please delete all the old output files and directories; and re-run this step!"
	echo "Script exit!"
	echo ""
	exit
    else
	ln -s $ref_genome_raw ref_genome.fa
    fi
fi

out_dir=${prefix}_canu_out

$canu_dir/canu -p $prefix -d $out_dir \
    useGrid=false \
    maxThreads=$threads \
    genomeSize=$genome_size \
    gnuplot=$gnuplot_dir/gnuplot \
    -${reads_type} $reads

# ln -s $out_dir/$prefix.correctedReads.fasta.gz .
perl $LRSDAY_HOME/scripts/simplify_seq_name.pl -i ./$out_dir/$prefix.contigs.fasta -o $prefix.canu.fa

# generate assembly statistics
perl $LRSDAY_HOME/scripts/cal_assembly_stats.pl -i $prefix.canu.fa -o $prefix.canu.stats.txt

# make the comparison between the assembled genome and the reference genome
$mummer_dir/nucmer -t $threads --maxmatch --nosimplify  -p $prefix.canu  $ref_genome_raw $prefix.canu.fa 
$mummer_dir/delta-filter -m  $prefix.canu.delta > $prefix.canu.delta_filter

# generate the vcf output
if [[ $vcf == "yes" ]]
then
    $mummer_dir/show-coords -b -T -r -c -l -d   $prefix.canu.delta_filter > $prefix.canu.filter.coords
    $mummer_dir/show-snps -C -T -l -r $prefix.canu.delta_filter > $prefix.canu.filter.snps
    perl $LRSDAY_HOME/scripts/mummer2vcf.pl -r ref_genome.fa -i $prefix.canu.filter.snps -t SNP -p $prefix.canu.filter
    perl $LRSDAY_HOME/scripts/mummer2vcf.pl -r ref_genome.fa -i $prefix.canu.filter.snps -t INDEL -p $prefix.canu.filter
    $samtools_dir/samtools faidx ref_genome.fa 
    awk '{printf("##contig=<ID=%s,length=%d>\n",$1,$2);}' ref_genome.fa.fai > $prefix.vcf_header.txt
    sed -i -e "/##reference/r $prefix.vcf_header.txt" $prefix.canu.filter.mummer2vcf.SNP.vcf
    sed -i -e "/##reference/r $prefix.vcf_header.txt" $prefix.canu.filter.mummer2vcf.INDEL.vcf
fi

# generate genome-wide dotplot
if [[ $dotplot == "yes" ]]
then
    $mummer_dir/mummerplot --large --postscript $prefix.canu.delta_filter -p $prefix.canu.filter
    perl $LRSDAY_HOME/scripts/fine_tune_gnuplot.pl -i $prefix.canu.filter.gp -o $prefix.canu.filter_adjust.gp -r ref_genome.fa -q ${prefix}.canu.fa
    $gnuplot_dir/gnuplot < $prefix.canu.filter_adjust.gp
fi

# clean up intermediate files
if [[ $debug == "no" ]]
then
    rm *.delta
    rm *.delta_filter
    # rm ref_genome.fa
    # rm ref_genome.fa.fai
    if [[ $vcf == "yes" ]] 
    then
	rm *.filter.coords
	rm $prefix.vcf_header.txt
	rm $prefix.canu.filter.snps
    fi
    if [[ $dotplot == "yes" ]]
    then
	rm *.filter.fplot
        rm *.filter.rplot
	rm *.filter.gp
        rm *.filter_adjust.gp
        rm *.filter.ps
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


