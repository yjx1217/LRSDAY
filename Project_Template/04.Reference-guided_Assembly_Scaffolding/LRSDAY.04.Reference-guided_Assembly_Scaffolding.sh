#!/bin/bash
set -e -o pipefail
#######################################
# load environment variables for LRSDAY
source ./../../env.sh
PATH=$gnuplot_dir:$hal_dir:$PATH

#######################################
# set project-specific variables
prefix="SK1" # The file name prefix for processing sample. Please avoid the character '.' in prefix. Default = "SK1" for the testing example.
input_assembly="./../03.Short-read-based_Assembly_Polishing/$prefix.assembly.short_read_polished.fa" # The file path of the input genome assembly.
ref_genome_raw="./../00.Reference_Genome/S288C.ASM205763v1.fa" # The file path of the raw reference genome.
ref_genome_noncore_masked="./../00.Reference_Genome/S288C.ASM205763v1.noncore_masked.fa" # The file path of the specially masked reference genome where subtelomeres and chromosome-ends were hard masked. When the subtelomere/chromosome-end information is unavailable for the organism that you are interested in, you can just put the path of the raw reference genome assembly here.
chrMT_tag="chrMT" # The sequence name for the mitochondrial genome in the raw reference genome file, if there are multiple reference mitochondrial genomes that you want to check, use a single ';' to separate them. e.g. "Sc_chrMT;Sp_chrMT". Default = "chrMT".
gap_size=5000 # The number of Ns to insert between adjacent contigs during scaffolding. Default = "5000".
scaffolder="ragout" # The reference-based assembly scaffolder to use: "ragout" or "ragoo". Default = "ragout". If the reference genome size is large (e.g. > 100 Mb), please use "ragoo" since extra dependency is needed for "ragout" to handle large genomes.
threads=1 # The number of threads to use. Default = "1".
debug="no" # Whether to keep intermediate files for debugging. Use "yes" if prefer to keep intermediate files, otherwise use "no". Default = "no".

######################################
# process the pipeline
#####################################
# run Ragout for scaffolding based on the reference genome

sed -e '/^[^>]/s/[^ATGCatgc]/N/g' $ref_genome_noncore_masked > ref_genome_noncore_masked.fa

if [[ $scaffolder == "ragout" ]]
then
    echo ".references = ref_genome" > ragout.recipe.txt
    echo ".target = $prefix" >> ragout.recipe.txt
    echo "ref_genome.fasta = ./ref_genome_noncore_masked.fa" >> ragout.recipe.txt
    echo "$prefix.fasta = $input_assembly" >> ragout.recipe.txt
    echo ".naming_ref = ref_genome" >> ragout.recipe.txt
    source $miniconda2_dir/activate $build_dir/conda_ragout_env
    python2 $ragout_dir/ragout -o ${prefix}_ragout_out --solid-scaffolds  -t $threads  ragout.recipe.txt
    source $miniconda2_dir/deactivate
    cat ./${prefix}_ragout_out/${prefix}_scaffolds.fasta | sed "s/^>chr_/>/g" > ./${prefix}_ragout_out/${prefix}_scaffolds.renamed.fasta 
    cat ./${prefix}_ragout_out/${prefix}_scaffolds.renamed.fasta ./${prefix}_ragout_out/${prefix}_unplaced.fasta > ./${prefix}_ragout_out/${prefix}.ragout.raw.fa    
    perl $LRSDAY_HOME/scripts/adjust_assembly_by_ragoutAGP.pl -i $input_assembly -p $prefix -a ./${prefix}_ragout_out/${prefix}_scaffolds.agp -g $gap_size
    ln -s ${prefix}.ragout.fa $prefix.assembly.ref_based_scaffolded.fa
elif [[ $scaffolder == "ragoo" ]]
then
    source $build_dir/py3_virtualenv_ragoo/bin/activate
    python3 $ragoo_dir/ragoo.py -t $threads -g $gap_size  -C -m $minimap2_dir/minimap2 $input_assembly ref_genome_noncore_masked.fa
    mv ragoo_output ${prefix}_ragoo_out
    cp ${prefix}_ragoo_out/ragoo.fasta ${prefix}.ragoo.fa
    ln -s ${prefix}.ragoo.fa $prefix.assembly.ref_based_scaffolded.fa
else
    echo "Unrecognized scaffolder: $scaffolder";
    echo "Please set scaffolder=\"ragout\" or \"ragoo\".";
    echo "LRSDAY Exit!";
    exit;
fi

# generate assembly statistics
perl $LRSDAY_HOME/scripts/cal_assembly_stats.pl -i $prefix.assembly.ref_based_scaffolded.fa -o $prefix.assembly.ref_based_scaffolded.stats.txt

# generate genome-wide dotplot
$mummer4_dir/nucmer -t $threads --maxmatch --nosimplify  -p $prefix.assembly.ref_based_scaffolded  $ref_genome_raw $prefix.assembly.ref_based_scaffolded.fa 
$mummer4_dir/delta-filter -m  $prefix.assembly.ref_based_scaffolded.delta > $prefix.assembly.ref_based_scaffolded.delta_filter
$mummer4_dir/show-coords -b -T -r -c -l -d   $prefix.assembly.ref_based_scaffolded.delta_filter > $prefix.assembly.ref_based_scaffolded.filter.coords
echo $chrMT_tag | sed -e "s/;/\n/g" > ref.chrMT.list
perl $LRSDAY_HOME/scripts/identify_contigs_for_RefChr_by_mummer.pl -i $prefix.assembly.ref_based_scaffolded.filter.coords -query_chr_list ref.chrMT.list -assembly_fasta $prefix.assembly.ref_based_scaffolded.fa -cov 75 -o $prefix.assembly.ref_based_scaffolded.mt_contig.list
$mummer4_dir/mummerplot --large --postscript $prefix.assembly.ref_based_scaffolded.delta_filter -p $prefix.assembly.ref_based_scaffolded.filter
perl $LRSDAY_HOME/scripts/fine_tune_gnuplot.pl -i $prefix.assembly.ref_based_scaffolded.filter.gp -o $prefix.assembly.ref_based_scaffolded.filter_adjust.gp -r $ref_genome_raw -q ${prefix}.assembly.ref_based_scaffolded.fa
$gnuplot_dir/gnuplot < $prefix.assembly.ref_based_scaffolded.filter_adjust.gp

# generate dotplot for the mitochondrial genome only
perl $LRSDAY_HOME/scripts/select_fasta_by_list.pl -i $ref_genome_raw -l ref.chrMT.list -m normal -o ref.chrMT.fa
perl $LRSDAY_HOME/scripts/select_fasta_by_list.pl -i $prefix.assembly.ref_based_scaffolded.fa -l $prefix.assembly.ref_based_scaffolded.mt_contig.list -m normal -o $prefix.assembly.ref_based_scaffolded.mt_contig.fa
$mummer4_dir/nucmer --maxmatch --nosimplify  -p $prefix.assembly.ref_based_scaffolded.chrMT ref.chrMT.fa $prefix.assembly.ref_based_scaffolded.mt_contig.fa
$mummer4_dir/delta-filter -m  $prefix.assembly.ref_based_scaffolded.chrMT.delta > $prefix.assembly.ref_based_scaffolded.chrMT.delta_filter
$mummer4_dir/mummerplot --large --postscript $prefix.assembly.ref_based_scaffolded.chrMT.delta_filter -p $prefix.assembly.ref_based_scaffolded.chrMT.filter
perl $LRSDAY_HOME/scripts/fine_tune_gnuplot.pl -i $prefix.assembly.ref_based_scaffolded.chrMT.filter.gp -o $prefix.assembly.ref_based_scaffolded.chrMT.filter_adjust.gp -r ref.chrMT.fa -q $prefix.assembly.ref_based_scaffolded.mt_contig.fa
$gnuplot_dir/gnuplot < $prefix.assembly.ref_based_scaffolded.chrMT.filter_adjust.gp

# clean up intermediate files
if [[ $debug == "no" ]]
then
    if [[ $scaffolder == "ragout" ]]
    then
	rm ragout.recipe.txt
    fi
    rm ref_genome_noncore_masked.fa
    rm *.filter.fplot
    rm *.filter.rplot
    rm *.delta
    rm *.delta_filter
    rm *.filter.gp
    rm *.filter_adjust.gp
    rm *.filter.ps
    rm ref.chrMT.list
    rm ref.chrMT.fa
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
