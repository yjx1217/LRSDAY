#!/bin/bash
set -e -o pipefail
#######################################
# load environment variables for LRSDAY
source ./../../env.sh
PATH=$gnuplot_dir:$PATH

#######################################
# set project-specific variables
prefix="SK1" # The file name prefix for the processing sample. Default = "SK1" for the testing example.
genome="./../06.Mitochondrial_Genome_Assembly_Improvement/$prefix.assembly.mt_improved.fa" # The file path of the input genome assembly.
vcf="yes" # Whether to generate a vcf file generated to show SNP and INDEL differences between the assembled genome and the reference genome for their uniquely alignable regions. Use "yes" if prefer to have vcf file generated to show SNP and INDEL differences between the assembled genome and the reference genome. Default = "yes".
dotplot="yes" # Whether to plot genome-wide dotplot based on the comparison with the reference genome below. Use "yes" if prefer to plot, otherwise use "no". Default = "yes".
ref_genome_raw="./../00.Ref_Genome/S288C.ASM205763v1.fa" # The path of the raw reference genome, only needed when dotplot="yes" or vcf="yes".
threads=1 # The number of threads to use. Default = 1.
debug="no" # Whether to keep intermediate files for debugging. Use "yes" if prefer to keep intermediate files, otherwise use "no". Default = "no".

#######################################
# process the pipeline
# Please mark desiarable changes in the $prefix.modification.list file and comment the command lines for Step 1 before proceeding with Step 2.
# Step 2:
perl $LRSDAY_HOME/scripts/relabel_and_reorder_sequences.pl  -i $genome -m $prefix.assembly.modification.list -o $prefix.assembly.final.fa
# generate assembly statistics
perl $LRSDAY_HOME/scripts/cal_assembly_stats.pl -i $prefix.assembly.final.fa -o $prefix.assembly.final.stats.txt

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
$mummer_dir/nucmer -t $threads --maxmatch --nosimplify  -p $prefix.assembly.final  $ref_genome_raw $prefix.assembly.final.fa 
$mummer_dir/delta-filter -m $prefix.assembly.final.delta > $prefix.assembly.final.delta_filter

# generate the vcf output
if [[ $vcf == "yes" ]]
then
    $mummer_dir/show-coords -b -T -r -c -l -d $prefix.assembly.final.delta_filter > $prefix.assembly.final.filter.coords
    $mummer_dir/show-snps -C -T -l -r $prefix.assembly.final.delta_filter > $prefix.assembly.final.filter.snps
    perl $LRSDAY_HOME/scripts/mummer2vcf.pl -r ref_genome.fa -i $prefix.assembly.final.filter.snps -t SNP -p $prefix.assembly.final.filter
    perl $LRSDAY_HOME/scripts/mummer2vcf.pl -r ref_genome.fa -i $prefix.assembly.final.filter.snps -t INDEL -p $prefix.assembly.final.filter
    $samtools_dir/samtools faidx ref_genome.fa 
    awk '{printf("##contig=<ID=%s,length=%d>\n",$1,$2);}' ref_genome.fa.fai > $prefix.vcf_header.txt
    sed -i -e "/##reference/r $prefix.vcf_header.txt" $prefix.assembly.final.filter.mummer2vcf.SNP.vcf
    sed -i -e "/##reference/r $prefix.vcf_header.txt" $prefix.assembly.final.filter.mummer2vcf.INDEL.vcf
fi

# generate genome-wide dotplot
if [[ $dotplot == "yes" ]]
then
    $mummer_dir/mummerplot --large --postscript $prefix.assembly.final.delta_filter -p $prefix.assembly.final.filter
    perl $LRSDAY_HOME/scripts/fine_tune_gnuplot.pl -i $prefix.assembly.final.filter.gp -o $prefix.assembly.final.filter_adjust.gp -r ref_genome.fa -q $prefix.assembly.final.fa
    $gnuplot_dir/gnuplot < $prefix.assembly.final.filter_adjust.gp
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
        rm $prefix.assembly.final.filter.snps
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
