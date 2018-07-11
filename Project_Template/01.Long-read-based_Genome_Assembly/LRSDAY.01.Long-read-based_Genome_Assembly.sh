#!/bin/bash
set -e -o pipefail
##########################################
# load environment variables for LRSDAY
source ./../../env.sh
PATH=$gnuplot_dir:$PATH

###########################################
# set project-specific variables
prefix="SK1" # file name prefix for the output files
reads="./../00.Long_Reads/SK1.filtered_subreads.fastq.gz" # path of the long reads file (in fastq or fastq.gz format)
reads_type="pacbio-raw" # long reads data type: "pacbio-raw" or "pacbio-corrected" or "nanopore-raw" or "nanopore-corrected"
genome_size="12.5m" # estimated genome size with the format of <number>[g|m|k], e.g. 12.5m for 12.5 Mb
assembler="canu" # long-read assembler: "canu" (default) or "flye" or "smartdenovo". Based on our test, canu gives better results but takes longer to finish.
customized_canu_parameters="-correctedErrorRate=0.04" # users can set customized Canu assembly parameters here or simply leave it empty like "" to use Canu's default assembly parameter. For example you could set "-correctedErrorRate=0.04" for high coverage (>60X) PacBio data and "-correctedErrorRate=0.12 -overlapper=mhap -utgReAlign=true" for high coverage (>60X) Nanopore data to improve the assembly speed. More than one customized parameters can be set here as long as they are separeted by space (e.g. "-option1=XXX -option2=YYY -option3=ZZZ"). Please consult Canu's manual "http://canu.readthedocs.io/en/latest/faq.html#what-parameters-can-i-tweak" for advanced customization settings. 
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



out_dir=${prefix}_${assembler}_out

if [[ "$assembler" == "canu" ]]
then
    $canu_dir/canu -p $prefix -d $out_dir \
	useGrid=false \
	maxThreads=$threads \
	genomeSize=$genome_size \
	gnuplot=$gnuplot_dir/gnuplot \
	-${reads_type} $reads \
	$customized_canu_parameters
    
    perl $LRSDAY_HOME/scripts/simplify_seq_name.pl -i ./$out_dir/$prefix.contigs.fasta -o $prefix.$assembler.fa
elif [[ "$assembler" == "flye" ]]
then
    if [[ "$reads_type" == "pacbio-corrected" ]]
    then
	reads_type="pacbio-corr"
    elif [[ "$reads_type" == "nanopore-raw" ]]
    then
        reads_type="nano-raw"
    elif [[ "$reads_type" == "nanopore-corrected" ]]
    then
        reads_type="nano-corr"
    fi
    $flye_dir/flye -o $out_dir \
	-t $threads \
	-g $genome_size \
	--${reads_type} $reads \
	--min-overlap 7000 \
	-i 2
    perl $LRSDAY_HOME/scripts/simplify_seq_name.pl -i ./$out_dir/contigs.fasta -o $prefix.$assembler.fa
elif [[ "$assembler" == "smartdenovo" ]]
then
    mkdir $out_dir
    cd $out_dir
    $smartdenovo_dir/smartdenovo_customized.pl -p $prefix -t $threads -c 1 ./../$reads > $prefix.mak
    make -f $prefix.mak
    perl $LRSDAY_HOME/scripts/simplify_seq_name.pl -i $prefix.dmo.cns  -o ./../$prefix.$assembler.fa
    cd ..
fi

# generate assembly statistics
perl $LRSDAY_HOME/scripts/cal_assembly_stats.pl -i $prefix.$assembler.fa -o $prefix.$assembler.stats.txt

# make the comparison between the assembled genome and the reference genome
$mummer_dir/nucmer -t $threads --maxmatch --nosimplify  -p $prefix.$assembler  $ref_genome_raw $prefix.$assembler.fa 
$mummer_dir/delta-filter -m  $prefix.$assembler.delta > $prefix.$assembler.delta_filter

# generate the vcf output
if [[ $vcf == "yes" ]]
then
    $mummer_dir/show-coords -b -T -r -c -l -d   $prefix.$assembler.delta_filter > $prefix.$assembler.filter.coords
    $mummer_dir/show-snps -C -T -l -r $prefix.$assembler.delta_filter > $prefix.$assembler.filter.snps
    perl $LRSDAY_HOME/scripts/mummer2vcf.pl -r ref_genome.fa -i $prefix.$assembler.filter.snps -t SNP -p $prefix.$assembler.filter
    perl $LRSDAY_HOME/scripts/mummer2vcf.pl -r ref_genome.fa -i $prefix.$assembler.filter.snps -t INDEL -p $prefix.$assembler.filter
    $samtools_dir/samtools faidx ref_genome.fa 
    awk '{printf("##contig=<ID=%s,length=%d>\n",$1,$2);}' ref_genome.fa.fai > $prefix.vcf_header.txt
    sed -i -e "/##reference/r $prefix.vcf_header.txt" $prefix.$assembler.filter.mummer2vcf.SNP.vcf
    sed -i -e "/##reference/r $prefix.vcf_header.txt" $prefix.$assembler.filter.mummer2vcf.INDEL.vcf
fi

# generate genome-wide dotplot
if [[ $dotplot == "yes" ]]
then
    $mummer_dir/mummerplot --large --postscript $prefix.$assembler.delta_filter -p $prefix.$assembler.filter
    perl $LRSDAY_HOME/scripts/fine_tune_gnuplot.pl -i $prefix.$assembler.filter.gp -o $prefix.$assembler.filter_adjust.gp -r ref_genome.fa -q ${prefix}.$assembler.fa
    $gnuplot_dir/gnuplot < $prefix.$assembler.filter_adjust.gp
fi

# clean up intermediate files
if [[ $debug == "no" ]]
then
    rm ref_genome.fa
    rm ref_genome.fa.fai
    rm *.delta
    rm *.delta_filter
    # rm ref_genome.fa
    # rm ref_genome.fa.fai
    if [[ $vcf == "yes" ]] 
    then
	rm *.filter.coords
	rm $prefix.vcf_header.txt
	rm $prefix.$assembler.filter.snps
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


