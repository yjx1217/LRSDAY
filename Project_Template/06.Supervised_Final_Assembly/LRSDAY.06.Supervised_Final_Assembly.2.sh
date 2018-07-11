#!/bin/bash
set -e -o pipefail
#######################################
# load environment variables for LRSDAY
source ./../../env.sh
PATH=$gnuplot_dir:$PATH

#######################################
# set project-specific variables
genome="./../05.Mitochondrial_Genome_Assembly_Improvement/SK1.mt_improved.fa" # path of the input genome assembly
prefix="SK1" # file name prefix for the output files
vcf="yes" # use "yes" if prefer to have vcf file generated to show SNP and INDEL differences between the assembled genome and the reference genome.
dotplot="yes" # use "yes" if prefer to plot genome-wide dotplot based on the comparison with the reference genome below, otherwise use "no"
ref_genome_raw="./../00.Ref_Genome/S288C.ASM205763v1.fa" # path of the raw reference genome, only needed when dotplot="yes" or vcf="yes".
threads=1 # number of threads to use
debug="no" # use "yes" if prefer to keep intermediate files, otherwise use "no".

#######################################
# process the pipeline
# Please mark desiarable changes in the $prefix.modification.list file and comment the command lines for Step 1 before proceeding with Step 2.
# Step 2:
perl $LRSDAY_HOME/scripts/relabel_and_reorder_sequences.pl  -i $genome -m ${prefix}.modification.list -o ${prefix}.final.fa
# generate assembly statistics
perl $LRSDAY_HOME/scripts/cal_assembly_stats.pl -i $prefix.final.fa -o $prefix.final.stats.txt

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

# make the comparison between the assembled genome and the reference genome
$mummer_dir/nucmer -t $threads --maxmatch --nosimplify  -p $prefix.final  $ref_genome_raw $prefix.final.fa 
$mummer_dir/delta-filter -m  $prefix.final.delta > $prefix.final.delta_filter

# generate the vcf output
if [[ $vcf == "yes" ]]
then
    $mummer_dir/show-coords -b -T -r -c -l -d   $prefix.final.delta_filter > $prefix.final.filter.coords
    $mummer_dir/show-snps -C -T -l -r $prefix.final.delta_filter > $prefix.final.filter.snps
    perl $LRSDAY_HOME/scripts/mummer2vcf.pl -r ref_genome.fa -i $prefix.final.filter.snps -t SNP -p $prefix.final.filter
    perl $LRSDAY_HOME/scripts/mummer2vcf.pl -r ref_genome.fa -i $prefix.final.filter.snps -t INDEL -p $prefix.final.filter
    $samtools_dir/samtools faidx ref_genome.fa 
    awk '{printf("##contig=<ID=%s,length=%d>\n",$1,$2);}' ref_genome.fa.fai > $prefix.vcf_header.txt
    sed -i -e "/##reference/r $prefix.vcf_header.txt" $prefix.final.filter.mummer2vcf.SNP.vcf
    sed -i -e "/##reference/r $prefix.vcf_header.txt" $prefix.final.filter.mummer2vcf.INDEL.vcf
fi

# generate genome-wide dotplot
if [[ $dotplot == "yes" ]]
then
    $mummer_dir/mummerplot --large --postscript $prefix.final.delta_filter -p $prefix.final.filter
    perl $LRSDAY_HOME/scripts/fine_tune_gnuplot.pl -i $prefix.final.filter.gp -o $prefix.final.filter_adjust.gp -r ref_genome.fa -q $prefix.final.fa
    $gnuplot_dir/gnuplot < $prefix.final.filter_adjust.gp
fi

# clean up intermediate files
if [[ $debug == "no" ]]
then
    rm *.delta
    rm *.delta_filter
    rm ref_genome.fa
    rm ref_genome.fa.fai
    if [[ $vcf == "yes" ]] 
    then
        rm *.filter.coords
        rm $prefix.vcf_header.txt
        rm $prefix.final.filter.snps
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
