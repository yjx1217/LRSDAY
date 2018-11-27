#!/bin/bash
set -e -o pipefail
##########################################
# load environment variables for LRSDAY
source ./../../env.sh
PATH=$gnuplot_dir:$PATH

###########################################
# set project-specific variables
prefix="SK1" # The file name prefix for the output files.
long_reads="./../00.Long_Reads/SK1.filtered_subreads.fastq.gz" # The file path of the long reads file (in fastq or fastq.gz format).
long_reads_type="pacbio-raw" # The long reads data type. Use "pacbio-raw" or "pacbio-corrected" or "nanopore-raw" or "nanopore-corrected". Default = "pacbio-raw" for the testing example
genome_size="12.5m" # The estimated genome size with the format of <number>[g|m|k], e.g. 12.5m for 12.5 Mb. Default = "12.5m".
assembler="canu" # The long-read assembler: Use "canu" or "flye" or "wtdbg2" or "smartdenovo" or "canu-flye" or "canu-wtdbg2" or "canu-smartdenovo". For "canu-flye", "canu-wtdbg2", and "canu-smartdenovo", the assembler canu is used first to generate error-corrected reads from the raw reads and then the assembler flye/wtdbg2/smartdenovo is used to assemble the genome. Based on our test, assembler="canu" generally gives the best result but will take substantially longer time than the other options.
customized_canu_parameters="correctedErrorRate=0.04" # For assembler="canu" only. Users can set customized Canu assembly parameters here or simply leave it empty like customized_canu_parameters="" to use Canu's default assembly parameter. For example you could set customized_canu_parameters="correctedErrorRate=0.04" for high coverage (>60X) PacBio data and customized_canu_parameters="overlapper=mhap;utgReAlign=true" for high coverage (>60X) Nanopore data to improve the assembly speed. When assembling genomes with high heterozygosity, you can could set customized_canu_parameters="corOutCoverage=200;batOptions=-dg 3 -db 3 -dr 1 -ca 500 -cp 50" to avoid collasping haplotypes. As shown in these examples, more than one customized parameters can be set here as long as they are separeted by a semicolon and contained in a pair of double quotes (e.g. customized_canu_parameters="option1=XXX;option2=YYY;option3=ZZZ"). Please consult Canu's manual "http://canu.readthedocs.io/en/latest/faq.html#what-parameters-can-i-tweak" for advanced customization settings. Default = "correctedErrorRate=0.04" for the testing example.
threads=2 # The number of threads to use. Default = 2.
vcf="yes" # Use "yes" if prefer to have vcf file generated to show SNP and INDEL differences between the assembled genome and the reference genome for their uniquely alignable regions. Otherwise use "no". Default = "yes".
dotplot="yes" # Use "yes" if prefer to plot genome-wide dotplot based on the comparison with the reference genome below. Otherwise use "no". Default = "yes".
ref_genome_raw="./../00.Ref_Genome/S288C.ASM205763v1.fa" # The file path of the raw reference genome. This is only needed when the option "dotplot=" or "vcf=" has been set as "yes".
debug="no" # Use "yes" if prefer to keep intermediate files. Otherwise use "no". Default = "no".

###########################################
# process the pipeline

# check project-specific variables
if [[ $vcf == "yes" || $dotplot == "yes" ]]
then
    # check if ref_genome_raw is defined
    if [[ -z "$ref_genome_raw" ]]
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
	cp $ref_genome_raw ref_genome.fa
    fi
fi

out_dir=${prefix}_${assembler}_out

if [[ "$assembler" == "canu" ]]
then
    OLDIFS=$IFS;
    IFS=";"
    customized_canu_parameters_array=($customized_canu_parameters)
    IFS=$OLDIFS;
    printf "%s\n" "${customized_canu_parameters_array[@]}" > $prefix.customized_canu_parameters.spec
    $canu_dir/canu -p $prefix -d $out_dir \
	-s $prefix.customized_canu_parameters.spec \
	useGrid=false \
	maxThreads=$threads \
	genomeSize=$genome_size \
	gnuplot=$gnuplot_dir/gnuplot \
	-${long_reads_type} $long_reads
    mv $prefix.customized_canu_parameters.spec ./$out_dir/
    perl $LRSDAY_HOME/scripts/simplify_seq_name.pl -i $out_dir/$prefix.contigs.fasta -o $prefix.assembly.$assembler.fa
elif [[ "$assembler" == "flye" ]]
then
    if [[ "$long_reads_type" == "pacbio-corrected" ]]
    then
	long_reads_type="pacbio-corr"
    elif [[ "$long_reads_type" == "nanopore-raw" ]]
    then
        long_reads_type="nano-raw"
    elif [[ "$long_reads_type" == "nanopore-corrected" ]]
    then
        long_reads_type="nano-corr"
    fi
    $flye_dir/flye -o $out_dir \
	-t $threads \
	-g $genome_size \
	--${long_reads_type} $long_reads \
	-i 2
    perl $LRSDAY_HOME/scripts/simplify_seq_name.pl -i $out_dir/scaffolds.fasta -o $prefix.assembly.$assembler.fa
elif [[ "$assembler" == "wtdbg2" ]]
then
    mkdir $out_dir
    cd $out_dir
    $wtdbg2_dir/wtdbg2 -t $threads -L 5000 -i ./../$long_reads -fo $prefix
    $wtdbg2_dir/wtpoa-cns -t $threads -i $prefix.ctg.lay.gz -fo $prefix.ctg.lay.fa
    cd ..
    perl $LRSDAY_HOME/scripts/simplify_seq_name.pl -i $out_dir/$prefix.ctg.lay.fa -o $prefix.assembly.$assembler.fa
elif [[ "$assembler" == "smartdenovo" ]]
then
    mkdir $out_dir
    cd $out_dir
    # $smartdenovo_dir/smartdenovo.pl -p $prefix -t $threads -c 1 ./../$long_reads > $prefix.mak
    $smartdenovo_dir/smartdenovo_customized.pl -p $prefix -t $threads -c 1 ./../$long_reads > $prefix.mak
    make -f $prefix.mak
    cd ..
    perl $LRSDAY_HOME/scripts/simplify_seq_name.pl -i $out_dir/$prefix.dmo.cns  -o $prefix.assembly.$assembler.fa

elif [[ "$assembler" == "canu-flye" ]]
then
    $canu_dir/canu -correct -p $prefix -d $out_dir/canu \
	useGrid=false \
	maxThreads=$threads \
	genomeSize=$genome_size \
	gnuplot=$gnuplot_dir/gnuplot \
	-${long_reads_type} $long_reads \
    
    if [[ "$long_reads_type" == "pacbio-raw" ]]
    then
	long_reads_type="pacbio-corr"
    elif [[ "$long_reads_type" == "pacbio-corrected" ]]
    then
	long_reads_type="pacbio-corr"
    elif [[ "$long_reads_type" == "nanopore-raw" ]]
    then
        long_reads_type="nano-corr"
    elif [[ "$long_reads_type" == "nanopore-corrected" ]]
    then
        long_reads_type="nano-corr"
    fi
    $flye_dir/flye -o $out_dir/flye \
	-t $threads \
	-g $genome_size \
	--${long_reads_type} $out_dir/canu/$prefix.correctedReads.fasta.gz \
	-i 2
    perl $LRSDAY_HOME/scripts/simplify_seq_name.pl -i $out_dir/flye/scaffolds.fasta -o $prefix.assembly.$assembler.fa
elif [[ "$assembler" == "canu-wtdbg2" ]]
then
    $canu_dir/canu -correct -p $prefix -d $out_dir/canu \
        useGrid=false \
        maxThreads=$threads \
        genomeSize=$genome_size \
        gnuplot=$gnuplot_dir/gnuplot \
        -${long_reads_type} $long_reads \

    mkdir -p $out_dir/wtdbg2
    cd $out_dir/wtdbg2
    $wtdbg2_dir/wtdbg2 -t $threads -L 5000 -i ./../canu/$prefix.correctedReads.fasta.gz -fo $prefix
    $wtdbg2_dir/wtpoa-cns -t $threads -i $prefix.ctg.lay.gz -fo $prefix.ctg.lay.fa
    cd ../..
    perl $LRSDAY_HOME/scripts/simplify_seq_name.pl -i $out_dir/wtdbg2/$prefix.ctg.lay.fa -o $prefix.assembly.$assembler.fa
elif [[ "$assembler" == "canu-smartdenovo" ]]
then
    $canu_dir/canu -correct -p $prefix -d $out_dir/canu \
        useGrid=false \
        maxThreads=$threads \
        genomeSize=$genome_size \
        gnuplot=$gnuplot_dir/gnuplot \
        -${long_reads_type} $long_reads \

    mkdir -p $out_dir/smartdenovo
    cd $out_dir/smartdenovo
    $smartdenovo_dir/smartdenovo.pl -p $prefix -t $threads -c 1 ./../canu/$prefix.correctedReads.fasta.gz  > $prefix.mak
    make -f $prefix.mak
    cd ../..
    perl $LRSDAY_HOME/scripts/simplify_seq_name.pl -i $out_dir/smartdenovo/$prefix.dmo.cns  -o $prefix.assembly.$assembler.fa
fi

ln -s $prefix.assembly.$assembler.fa $prefix.assembly.raw.fa

# generate assembly statistics
perl $LRSDAY_HOME/scripts/cal_assembly_stats.pl -i $prefix.assembly.raw.fa -o $prefix.assembly.raw.stats.txt

# make the comparison between the assembled genome and the reference genome
$mummer_dir/nucmer -t $threads --maxmatch --nosimplify  -p $prefix.assembly.raw  $ref_genome_raw $prefix.assembly.raw.fa 
$mummer_dir/delta-filter -m  $prefix.assembly.raw.delta > $prefix.assembly.raw.delta_filter

# generate the vcf output
if [[ $vcf == "yes" ]]
then
    $mummer_dir/show-coords -b -T -r -c -l -d   $prefix.assembly.raw.delta_filter > $prefix.assembly.raw.filter.coords
    $mummer_dir/show-snps -C -T -l -r $prefix.assembly.raw.delta_filter > $prefix.assembly.raw.filter.snps
    perl $LRSDAY_HOME/scripts/mummer2vcf.pl -r ref_genome.fa -i $prefix.assembly.raw.filter.snps -t SNP -p $prefix.assembly.raw.filter
    perl $LRSDAY_HOME/scripts/mummer2vcf.pl -r ref_genome.fa -i $prefix.assembly.raw.filter.snps -t INDEL -p $prefix.assembly.raw.filter
    $samtools_dir/samtools faidx ref_genome.fa 
    awk '{printf("##contig=<ID=%s,length=%d>\n",$1,$2);}' ref_genome.fa.fai > $prefix.vcf_header.txt
    sed -i -e "/##reference/r $prefix.vcf_header.txt" $prefix.assembly.raw.filter.mummer2vcf.SNP.vcf
    sed -i -e "/##reference/r $prefix.vcf_header.txt" $prefix.assembly.raw.filter.mummer2vcf.INDEL.vcf
fi

# generate genome-wide dotplot
if [[ $dotplot == "yes" ]]
then
    $mummer_dir/mummerplot --large --postscript $prefix.assembly.raw.delta_filter -p $prefix.assembly.raw.filter
    perl $LRSDAY_HOME/scripts/fine_tune_gnuplot.pl -i $prefix.assembly.raw.filter.gp -o $prefix.assembly.raw.filter_adjust.gp -r ref_genome.fa -q $prefix.assembly.raw.fa
    $gnuplot_dir/gnuplot < $prefix.assembly.raw.filter_adjust.gp
fi

# clean up intermediate files
if [[ $debug == "no" ]]
then
    rm *.delta
    rm *.delta_filter
    rm ref_genome.fa
    if [[ $vcf == "yes" ]] 
    then
	rm *.filter.coords
	rm $prefix.vcf_header.txt
	rm $prefix.assembly.raw.filter.snps
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


