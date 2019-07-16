#!/bin/bash
set -e -o pipefail
##########################################
# load environment variables for LRSDAY
source ./../../env.sh
PATH=$gnuplot_dir:$PATH

###########################################
# set project-specific variables
prefix="SK1" # The file name prefix for the processing sample. Default = "SK1" for the testing example.
long_reads="./../00.Long_Reads/$prefix.filtlong.fastq.gz" # The file path of the long reads file (in fastq or fastq.gz format).
long_reads_type="pacbio-raw" # The long reads data type. Use "pacbio-raw" or "pacbio-corrected" or "nanopore-raw" or "nanopore-corrected". Default = "pacbio-raw" for the testing example
genome_size="12.5m" # The estimated genome size with the format of <number>[g|m|k], e.g. 12.5m for 12.5 Mb. Default = "12.5m".
assembler="canu" # The long-read assembler to use. Use "canu" or "flye" or "wtdbg2" or "smartdenovo" or "ra" or "shasta" or "canu-flye" or "canu-wtdbg2" or "canu-smartdenovo" or "canu-ra" or "canu-shasta". For "canu-flye", "canu-wtdbg2", "canu-smartdenovo", "canu-ra", or "canu-shasta", the assembler canu is used first to generate error-corrected reads from the raw reads and then the assembler flye/wtdbg2/smartdenovo/ra/shasta is used to assemble the genome. Default = "canu".
customized_canu_parameters="" # When assembler="canu" or "canu-flye", "canu-wtdbg2", "canu-smartdenovo", "canu-ra", or "canu-shasta", users can set customized Canu assembly parameters here when needed. If users want to keep Canu's default parameters or when other assembler is used, simply leave it empty like customized_canu_parameters="". For example users could set customized_canu_parameters="correctedErrorRate=0.04" for high coverage (>60X) PacBio data and customized_canu_parameters="overlapper=mhap;utgReAlign=true" for high coverage (>60X) Nanopore data to improve the assembly speed. When assembling genomes with high heterozygosity, you can could set customized_canu_parameters="corOutCoverage=200;batOptions=-dg 3 -db 3 -dr 1 -ca 500 -cp 50" to avoid collasping haplotypes. As shown in these examples, more than one customized parameters can be set here as long as they are separeted by a semicolon and contained in a pair of double quotes (e.g. customized_canu_parameters="option1=XXX;option2=YYY;option3=ZZZ"). Please consult Canu's manual "http://canu.readthedocs.io/en/latest/faq.html#what-parameters-can-i-tweak" for advanced customization settings. Default = "" for the testing example.
threads=4 # The number of threads to use. Default = 4.
vcf="yes" # Use "yes" if prefer to have vcf file generated to show SNP and INDEL differences between the assembled genome and the reference genome for their uniquely alignable regions. Otherwise use "no". Default = "yes".
dotplot="yes" # Use "yes" if prefer to plot genome-wide dotplot based on the comparison with the reference genome below. Otherwise use "no". Default = "yes".

canu_triobinning_mode="no" # Whether to use Canu's TrioBinning pipeline when assembler=canu. Default = "no".
if [[ "$canu_triobinning_mode" == "no" ]]
then
    ref_genome_raw="./../00.Reference_Genome/S288C.ASM205763v1.fa" # The file path of the raw reference genome. This is only needed when the option "dotplot" or "vcf" has been set as "yes".
else
    parent1_tag="Parent1" # The name tag for parent1. Default = "Parent1".
    parent1_short_reads="./../00.Short_Reads/$parent1_tag.illumina.fastq.gz" # The relative path to Illumina reads of parent1 of the long-read sequenced hybrid. This is only needed when the "canu_triobinning_mode" option has been set as "yes". You can combine the R1 and R2 reads into a single fastq file here.
    parent1_ref_genome_raw="./../00.Reference_Genome/$parent1_tag.genome.raw.fa" # The relative path to the reference genome of parent1 of the hybrid. This is only needed when the options "canu_triobinning_mode" and "dotplot" (or "canu_triobinning_mode" and "vcf") have been set as "yes". When not available, it can be set to a general reference genome of the organism. 

    parent2_tag="Parent2" # The name tag for parent2. Default = "Parent2".
    parent2_short_reads="./../00.Short_Reads/$parent2_tag.illumina.fastq.gz" # The relative path to Illumina reads of parent2 of the long-read sequenced hybrid. This is only needed when the "canu_triobinning_mode" option has been set as "yes". You can combine the R1 and R2 reads into a single fastq file here.
    parent2_ref_genome_raw="./../00.Reference_Genome/$parent2_tag.genome.raw.fa" # The relative path to the reference genome of parent1 of the hybrid. This is only needed when the options "canu_triobinning_mode" and "dotplot" (or "canu_triobinning_mode" and "vcf") have been set as "yes". When not available, it can be set to a general reference genome of the organism. 
fi

debug="no" # Use "yes" if prefer to keep intermediate files. Otherwise use "no". Default = "no".

###########################################
# process the pipeline

# check project-specific variables
if [[ $vcf == "yes" || $dotplot == "yes" ]]
then
    if [[ $canu_triobinning_mode == "no" ]]
    then
	# check if ref_genome_raw has been defined
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
    else
	# check if parent1_ref_genome_raw and parent2_ref_genome_raw have been defined   
	if [[ -z "$parent1_ref_genome_raw" ]]
	then
	    echo "The vcf and doptlot outputs require the variable parent1_ref_genome_raw to be defined when canu_triobinning_mode=\"yes\"!"
	    echo "Please define this variable in the script; Please delete all the old output files and directories; and re-run this step!"
	    echo "Script exit!"
	    echo ""
	    exit
	elif [[ ! -f $parent1_ref_genome_raw ]]
	then
	    echo "The vcf and doptlot outputs require the $parent1_ref_genome_raw as defined with the variable \"parent1_ref_genome raw\" but this file cannot be found!"
	    echo "Please make sure that this file truly exists; Please delete all the old output files and directories; and re-run this step!"
	    echo "Script exit!"
	    echo ""
	    exit
	else
	    cp $parent1_ref_genome_raw parent1.ref_genome.raw.fa
	    sed -e "s/>/>${parent1_tag}_/g"  < parent1.ref_genome.raw.fa >  parent1.ref_genome.fa
	    if [[ "$debug" == "no" ]]
	    then
		rm parent1.ref_genome.raw.fa
	    fi
	    
	fi
	
	if [[ -z "$parent2_ref_genome_raw" ]]
	then
	    echo "The vcf and doptlot outputs require the variable parent2_ref_genome_raw to be defined when canu_triobinning_mode=\"yes\"!"
	    echo "Please define this variable in the script; Please delete all the old output files and directories; and re-run this step!"
	    echo "Script exit!"
	    echo ""
	    exit
	elif [[ ! -f $parent2_ref_genome_raw ]]
	then
	    echo "The vcf and doptlot outputs require the $parent2_ref_genome_raw as defined with the variable \"parent2_ref_genome raw\" but this file cannot be found!"
	    echo "Please make sure that this file truly exists; Please delete all the old output files and directories; and re-run this step!"
	    echo "Script exit!"
	    echo ""
	    exit
	else
	    cp $parent2_ref_genome_raw parent2.ref_genome.raw.fa
	    sed -e "s/>/>${parent2_tag}_/g"  < parent2.ref_genome.raw.fa >  parent2.ref_genome.fa 
	    if [[ "$debug" == "no" ]]
	    then
		rm parent2.ref_genome.raw.fa
	    fi
	fi
	cat parent1.ref_genome.fa parent2.ref_genome.fa > ref_genome.fa
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
    if [[ $canu_triobinning_mode == "no" ]]
    then
	$canu_dir/canu -p $prefix -d $out_dir \
	    -s $prefix.customized_canu_parameters.spec \
	    useGrid=false \
	    maxThreads=$threads \
	    genomeSize=$genome_size \
	    gnuplot=$gnuplot_dir/gnuplot \
	    -${long_reads_type} $long_reads 
	mv $prefix.customized_canu_parameters.spec ./$out_dir/
	perl $LRSDAY_HOME/scripts/simplify_seq_name.pl -i $out_dir/$prefix.contigs.fasta -o $prefix.assembly.$assembler.fa
    else
	$canu_dir/canu -p $prefix -d $out_dir \
            -s $prefix.customized_canu_parameters.spec \
            useGrid=false \
            maxThreads=$threads \
            genomeSize=$genome_size \
            gnuplot=$gnuplot_dir/gnuplot \
            -${long_reads_type} $long_reads \
	    -haplotype1 $parent1_short_reads \
	    -haplotype2 $parent2_short_reads 
	mv $prefix.customized_canu_parameters.spec ./$out_dir/
	perl $LRSDAY_HOME/scripts/simplify_seq_name.pl -i $out_dir/${prefix}-haplotype1/${prefix}-haplotype1.contigs.fasta -o $prefix.haplotype1.assembly.$assembler.raw.fa
	sed -e "s/>/>${parent1_tag}_/g" < $prefix.haplotype1.assembly.$assembler.raw.fa > $prefix.haplotype1.assembly.$assembler.fa
	perl $LRSDAY_HOME/scripts/simplify_seq_name.pl -i $out_dir/${prefix}-haplotype2/${prefix}-haplotype2.contigs.fasta -o $prefix.haplotype2.assembly.$assembler.raw.fa
	sed -e "s/>/>${parent2_tag}_/g" < $prefix.haplotype2.assembly.$assembler.raw.fa > $prefix.haplotype2.assembly.$assembler.fa
	cat $prefix.haplotype1.assembly.$assembler.fa $prefix.haplotype2.assembly.$assembler.fa > $prefix.assembly.$assembler.fa
	if [[ "$debug" == "no" ]]
	then
	    rm $prefix.haplotype1.assembly.$assembler.raw.fa
	    rm $prefix.haplotype2.assembly.$assembler.raw.fa
	fi
    fi
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
    perl $LRSDAY_HOME/scripts/simplify_seq_name.pl -i $out_dir/assembly.fasta -o $prefix.assembly.$assembler.fa
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
elif [[ "$assembler" == "ra" ]]
then
    mkdir $out_dir
    cd $out_dir
    if [[ "$long_reads_type" == "pacbio-raw" || "$long_reads_type" == "pacbio-corrected" ]]
    then
        long_reads_type="pb"
    elif [[ "$long_reads_type" == "nanopore-raw" || "$long_reads_type" == "nanopore-corrected" ]]
    then
        long_reads_type="ont"
    fi
    $ra_dir/ra -x $long_reads_type -t $threads ./../$long_reads > $prefix.assembly.$assembler.fa
    cd ..
    perl $LRSDAY_HOME/scripts/simplify_seq_name.pl -i $out_dir/$prefix.assembly.$assembler.fa -o $prefix.assembly.$assembler.fa
elif [[ "$assembler" == "shasta" ]]
then
    if [[ "$long_reads_type" == "pacbio-raw" || "$long_reads_type" == "pacbio-corrected" ]]
    then
        long_reads_type="pb"
    elif [[ "$long_reads_type" == "nanopore-raw" || "$long_reads_type" == "nanopore-corrected" ]]
    then
        long_reads_type="ont"
    fi
    perl $LRSDAY_HOME/scripts/fastq2fasta.pl -i $long_reads -o $prefix.long_reads.fasta
    $shasta_dir/shasta --input $prefix.long_reads.fasta  --output $out_dir
    perl $LRSDAY_HOME/scripts/simplify_seq_name.pl -i $out_dir/Assembly.fasta -o $prefix.assembly.$assembler.fa
    rm $prefix.long_reads.fasta
elif [[ "$assembler" == "canu-flye" ]]
then
    OLDIFS=$IFS;
    IFS=";"
    customized_canu_parameters_array=($customized_canu_parameters)
    IFS=$OLDIFS;
    printf "%s\n" "${customized_canu_parameters_array[@]}" > $prefix.customized_canu_parameters.spec
    $canu_dir/canu -correct -p $prefix -d $out_dir/canu \
	-s $prefix.customized_canu_parameters.spec \
	useGrid=false \
	maxThreads=$threads \
	genomeSize=$genome_size \
	gnuplot=$gnuplot_dir/gnuplot \
	-${long_reads_type} $long_reads
    mv $prefix.customized_canu_parameters.spec ./$out_dir/canu
    if [[ "$long_reads_type" == "pacbio-raw" || "$long_reads_type" == "pacbio-corrected" ]]
    then
	long_reads_type="pacbio-corr"
    elif [[ "$long_reads_type" == "nanopore-raw" || "$long_reads_type" == "nanopore-corrected" ]]
    then
        long_reads_type="nano-corr"
    fi
    $flye_dir/flye -o $out_dir/flye \
	-t $threads \
	-g $genome_size \
	--${long_reads_type} $out_dir/canu/$prefix.correctedReads.fasta.gz \
	-i 2
    perl $LRSDAY_HOME/scripts/simplify_seq_name.pl -i $out_dir/flye/assembly.fasta -o $prefix.assembly.$assembler.fa
elif [[ "$assembler" == "canu-wtdbg2" ]]
then
    OLDIFS=$IFS;
    IFS=";"
    customized_canu_parameters_array=($customized_canu_parameters)
    IFS=$OLDIFS;
    printf "%s\n" "${customized_canu_parameters_array[@]}" > $prefix.customized_canu_parameters.spec
    $canu_dir/canu -correct -p $prefix -d $out_dir/canu \
	-s $prefix.customized_canu_parameters.spec \
        useGrid=false \
        maxThreads=$threads \
        genomeSize=$genome_size \
        gnuplot=$gnuplot_dir/gnuplot \
        -${long_reads_type} $long_reads 
    mv $prefix.customized_canu_parameters.spec ./$out_dir/canu
    mkdir -p $out_dir/wtdbg2
    cd $out_dir/wtdbg2
    $wtdbg2_dir/wtdbg2 -t $threads -L 5000 -i ./../canu/$prefix.correctedReads.fasta.gz -fo $prefix
    $wtdbg2_dir/wtpoa-cns -t $threads -i $prefix.ctg.lay.gz -fo $prefix.ctg.lay.fa
    cd ../..
    perl $LRSDAY_HOME/scripts/simplify_seq_name.pl -i $out_dir/wtdbg2/$prefix.ctg.lay.fa -o $prefix.assembly.$assembler.fa
elif [[ "$assembler" == "canu-smartdenovo" ]]
then
    OLDIFS=$IFS;
    IFS=";"
    customized_canu_parameters_array=($customized_canu_parameters)
    IFS=$OLDIFS;
    printf "%s\n" "${customized_canu_parameters_array[@]}" > $prefix.customized_canu_parameters.spec

    $canu_dir/canu -correct -p $prefix -d $out_dir/canu \
	-s $prefix.customized_canu_parameters.spec \
        useGrid=false \
        maxThreads=$threads \
        genomeSize=$genome_size \
        gnuplot=$gnuplot_dir/gnuplot \
        -${long_reads_type} $long_reads
    mv $prefix.customized_canu_parameters.spec ./$out_dir/canu
    mkdir -p $out_dir/smartdenovo
    cd $out_dir/smartdenovo
    $smartdenovo_dir/smartdenovo.pl -p $prefix -t $threads -c 1 ./../canu/$prefix.correctedReads.fasta.gz  > $prefix.mak
    make -f $prefix.mak
    cd ../..
    perl $LRSDAY_HOME/scripts/simplify_seq_name.pl -i $out_dir/smartdenovo/$prefix.dmo.cns  -o $prefix.assembly.$assembler.fa
elif [[ "$assembler" == "canu-ra" ]]
then
    OLDIFS=$IFS;
    IFS=";"
    customized_canu_parameters_array=($customized_canu_parameters)
    IFS=$OLDIFS;
    printf "%s\n" "${customized_canu_parameters_array[@]}" > $prefix.customized_canu_parameters.spec

    $canu_dir/canu -correct -p $prefix -d $out_dir/canu \
	-s $prefix.customized_canu_parameters.spec \
        useGrid=false \
        maxThreads=$threads \
        genomeSize=$genome_size \
        gnuplot=$gnuplot_dir/gnuplot \
        -${long_reads_type} $long_reads 
    mv $prefix.customized_canu_parameters.spec ./$out_dir/canu
    if [[ "$long_reads_type" == "pacbio-raw" || "$long_reads_type" == "pacbio-corrected" ]]
    then
	long_reads_type="pb"
    elif [[ "$long_reads_type" == "nanopore-raw" || "$long_reads_type" == "nanopore-corrected" ]]
    then
        long_reads_type="ont"
    fi
    mkdir -p $out_dir/ra
    cd $out_dir/ra
    $ra_dir/ra -x $long_reads_type -t $threads ./../canu/$prefix.correctedReads.fasta.gz > $prefix.assembly.$assembler.fa
    cd ../..
    perl $LRSDAY_HOME/scripts/simplify_seq_name.pl -i $out_dir/ra/$prefix.assembly.$assembler.fa -o $prefix.assembly.$assembler.fa
elif [[ "$assembler" == "canu-shasta" ]]
then
    OLDIFS=$IFS;
    IFS=";"
    customized_canu_parameters_array=($customized_canu_parameters)
    IFS=$OLDIFS;
    printf "%s\n" "${customized_canu_parameters_array[@]}" > $prefix.customized_canu_parameters.spec

    $canu_dir/canu -correct -p $prefix -d $out_dir/canu \
	-s $prefix.customized_canu_parameters.spec \
        useGrid=false \
        maxThreads=$threads \
        genomeSize=$genome_size \
        gnuplot=$gnuplot_dir/gnuplot \
        -${long_reads_type} $long_reads 
    mv $prefix.customized_canu_parameters.spec ./$out_dir/canu
    if [[ "$long_reads_type" == "pacbio-raw" || "$long_reads_type" == "pacbio-corrected" ]]
    then
	long_reads_type="pb"
    elif [[ "$long_reads_type" == "nanopore-raw" || "$long_reads_type" == "nanopore-corrected" ]]
    then
        long_reads_type="ont"
    fi
    gunzip < $out_dir/canu/$prefix.correctedReads.fasta.gz > $out_dir/$prefix.correctedReads.fasta
    $shasta_dir/shasta --input $out_dir/$prefix.correctedReads.fasta  --output $out_dir/shasta
    rm $out_dir/$prefix.correctedReads.fasta
    perl $LRSDAY_HOME/scripts/simplify_seq_name.pl -i $out_dir/shasta/Assembly.fasta -o $prefix.assembly.$assembler.fa
fi

if [[ "$canu_triobinning_mode" == "no" ]]
then
     ln -s $prefix.assembly.$assembler.fa $prefix.assembly.raw.fa
    # generate assembly statistics
    perl $LRSDAY_HOME/scripts/cal_assembly_stats.pl -i $prefix.assembly.raw.fa -o $prefix.assembly.raw.stats.txt

    # make the comparison between the assembled genome and the reference genome
    $mummer4_dir/nucmer -t $threads --maxmatch --nosimplify -p $prefix.assembly.raw  ref_genome.fa $prefix.assembly.raw.fa 
    $mummer4_dir/delta-filter -m $prefix.assembly.raw.delta > $prefix.assembly.raw.delta_filter
    
    # generate the vcf output
    if [[ $vcf == "yes" ]]
    then
	$mummer4_dir/show-coords -b -T -r -c -l -d  $prefix.assembly.raw.delta_filter > $prefix.assembly.raw.filter.coords
	$mummer4_dir/show-snps -C -T -l -r $prefix.assembly.raw.delta_filter > $prefix.assembly.raw.filter.snps
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
	$mummer4_dir/mummerplot --large --postscript $prefix.assembly.raw.delta_filter -p $prefix.assembly.raw.filter
	perl $LRSDAY_HOME/scripts/fine_tune_gnuplot.pl -i $prefix.assembly.raw.filter.gp -o $prefix.assembly.raw.filter_adjust.gp -r ref_genome.fa -q $prefix.assembly.raw.fa
	$gnuplot_dir/gnuplot < $prefix.assembly.raw.filter_adjust.gp
    fi
else
    
    ln -s $prefix.assembly.$assembler.fa $prefix.assembly.raw.fa
    ln -s $prefix.haplotype1.assembly.$assembler.fa $prefix.haplotype1.assembly.raw.fa
    ln -s $prefix.haplotype2.assembly.$assembler.fa $prefix.haplotype2.assembly.raw.fa
    # generate assembly statistics
    perl $LRSDAY_HOME/scripts/cal_assembly_stats.pl -i $prefix.assembly.raw.fa -o $prefix.assembly.raw.stats.txt
    perl $LRSDAY_HOME/scripts/cal_assembly_stats.pl -i $prefix.haplotype1.assembly.raw.fa -o $prefix.haplotype1.assembly.raw.stats.txt
    perl $LRSDAY_HOME/scripts/cal_assembly_stats.pl -i $prefix.haplotype2.assembly.raw.fa -o $prefix.haplotype2.assembly.raw.stats.txt
    
    # make the comparison between the assembled genome and the reference genome
    $mummer4_dir/nucmer -t $threads --maxmatch --nosimplify  -p $prefix.assembly.raw ref_genome.fa $prefix.assembly.raw.fa 
    $mummer4_dir/delta-filter -m  $prefix.assembly.raw.delta > $prefix.assembly.raw.delta_filter
    
    # generate the vcf output
    if [[ $vcf == "yes" ]]
    then
	$mummer4_dir/show-coords -b -T -r -c -l -d   $prefix.assembly.raw.delta_filter > $prefix.assembly.raw.filter.coords
	$mummer4_dir/show-snps -C -T -l -r $prefix.assembly.raw.delta_filter > $prefix.assembly.raw.filter.snps
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
	$mummer4_dir/mummerplot --large --postscript $prefix.assembly.raw.delta_filter -p $prefix.assembly.raw.filter
	perl $LRSDAY_HOME/scripts/fine_tune_gnuplot.pl -i $prefix.assembly.raw.filter.gp -o $prefix.assembly.raw.filter_adjust.gp -r ref_genome.fa -q $prefix.assembly.raw.fa
	$gnuplot_dir/gnuplot < $prefix.assembly.raw.filter_adjust.gp
    fi
fi

# clean up intermediate files
if [[ $debug == "no" ]]
then
    rm *.delta
    rm *.delta_filter

    if [[ "$canu_triobinning_mode" == "yes" ]]
    then
	rm parent1.ref_genome.fa
        rm parent2.ref_genome.fa
    fi

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


