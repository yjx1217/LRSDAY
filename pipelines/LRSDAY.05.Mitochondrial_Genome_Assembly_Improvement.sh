#!/bin/bash
set -e
##########################################
# load environment variables for LRSDAY
source ./../../env.sh
PATH=$bwa_dir:$samtools_dir:$gnuplot_dir:$canu_dir:$mummer_legacy_dir:$spades_dir:$prodigal_dir:$PATH

###########################################
# set project-specific variables
genome="./../03.Reference-guided_Assembly_Scaffolding/SK1.ragout.fa" # path of the input genome assembly
prefix="SK1" # file name prefix for the output files
mt_contig_list="./../03.Reference-guided_Assembly_Scaffolding/SK1.mt_contig.list" # the mitochodnrial contig list generated in 03.ChromosomeAssignment_Provisional
gene_start="$LRSDAY_HOME/data/ATP6.cds.fa" # can be set to any gene as long as a fasta file containing the DNA sequence of the gene is provided.
ref_genome_raw="./../00.Ref_Genome/S288C.ASM205763v1.fa" # path of the raw reference genome 
chrMT_tag="chrMT" # sequence name for the mitochondrial genome in the raw reference genome file
threads=1 # number of threads to use
debug="no" # use "yes" if prefer to keep intermediate files, otherwise use "no".

##########################################
# process the pipeline
$LRSDAY_HOME/scripts/select_fasta_by_list.pl -i $genome -l $mt_contig_list -m reverse -o $prefix.non_mt_contig.fa
$LRSDAY_HOME/scripts/select_fasta_by_list.pl -i $genome -l $mt_contig_list -m normal -o $prefix.mt_contig.fa

# break the mitochondrial scaffold by the gaps introduced during reference-based scaffolding
perl $LRSDAY_HOME/scripts/break_scaffolds_by_N.pl -i $prefix.mt_contig.fa -o $prefix.mt_contig.descaffold.fa -g 5

if [[ $(egrep -c "^>" "$prefix.mt_contig.descaffold.fa") -eq 1 ]]
then
    echo "A single contig corresponding to chrMT is found, will run circularization directly ..."
    $circlator_dir/circlator fixstart --genes_fa $gene_start --min_id 90  $prefix.mt_contig.descaffold.fa  ${prefix}.mt_contig.circlator
else
    echo "Multiple contigs corresponding to chrMT are found, will run fragmentation, re-assembly and circularization ..."
    db="$prefix.mt_contig.descaffold.fa"
    db_tag="$prefix.mt_contig.descaffold"
    $blast_dir/makeblastdb -in $db  -dbtype nucl  -title $db_tag  -hash_index -out $db_tag
    $blast_dir/blastn -query $gene_start -db $db_tag -evalue 1E-6 \
	-outfmt '7 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore qlen qcovs' >$prefix.mt_contig.descaffold.tblastx.out

    perl $LRSDAY_HOME/scripts/filter_blast_result.pl -i $prefix.mt_contig.descaffold.tblastx.out -o $prefix.mt_contig.descaffold.tblastx.filtered.out -pct_identity_cutoff 90 -query_cov_cutoff 90
    perl $LRSDAY_HOME/scripts/break_contig_by_blast.pl -i $prefix.mt_contig.descaffold.fa -b $prefix.mt_contig.descaffold.tblastx.filtered.out -o $prefix.mt_contig.break_by_blast.fa

    # reassemble the mitochondrial genome using CAP3 to collapse smaller fragments
    $cap_dir/cap3 $prefix.mt_contig.break_by_blast.fa -k 0 
    cat $prefix.mt_contig.break_by_blast.fa.cap.contigs| sed "s/>/>chrMT_/" > $prefix.mt_contig.for_fixstart.fa

    db="$prefix.mt_contig.for_fixstart.fa"
    db_tag="$prefix.mt_contig.for_fixstart"
    $blast_dir/makeblastdb -in $db  -dbtype nucl  -title $db_tag  -hash_index -out $db_tag
    $blast_dir/blastn -query $gene_start -db $db_tag -evalue 1E-6 \
	-outfmt '7 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore qlen qcovs' >$prefix.mt_contig.for_fixstart.tblastx.out

    perl $LRSDAY_HOME/scripts/filter_blast_result.pl -i $prefix.mt_contig.for_fixstart.tblastx.out -o $prefix.mt_contig.for_fixstart.tblastx.filtered.out -pct_identity_cutoff 90 -query_cov_cutoff 90

    cat $prefix.mt_contig.for_fixstart.tblastx.filtered.out |egrep -v "^#" |cut -f 2 > $prefix.mt_contig.for_fixstart.list
    $LRSDAY_HOME/scripts/select_fasta_by_list.pl -i $prefix.mt_contig.for_fixstart.fa -l $prefix.mt_contig.for_fixstart.list -m reverse -o $prefix.for_fixstart.for_skip.fa
    cat $prefix.for_fixstart.for_skip.fa |egrep ">" |sed "s/>//" > $prefix.for_fixstart.for_skip.list

    $circlator_dir/circlator fixstart --genes_fa $gene_start --min_id 90 --ignore $prefix.for_fixstart.for_skip.list $prefix.mt_contig.for_fixstart.fa  ${prefix}.mt_contig.circlator
fi

mv ${prefix}.mt_contig.circlator.fasta ${prefix}.mt_contig.circlator.fa
cat $prefix.non_mt_contig.fa |egrep ">tig" |sed "s/>//" > non_primary_contig.list
$LRSDAY_HOME/scripts/select_fasta_by_list.pl -i $prefix.non_mt_contig.fa -l non_primary_contig.list -m reverse -o $prefix.primary_contig.fa
$LRSDAY_HOME/scripts/select_fasta_by_list.pl -i $prefix.non_mt_contig.fa -l non_primary_contig.list -m normal -o $prefix.non_primary_contig.fa
cat $prefix.primary_contig.fa $prefix.mt_contig.circlator.fa $prefix.non_primary_contig.fa > $prefix.mt_improved.fa

# generate dotplot for the mitochondrial genome only
echo $chrMT_tag > ref.chrMT.list
perl $LRSDAY_HOME/scripts/select_fasta_by_list.pl -i $ref_genome_raw -l ref.chrMT.list -m normal  -o ref.chrMT.fa
$mummer_dir/nucmer -t $threads --maxmatch --nosimplify  -p $prefix.mt_improved.chrMT ref.chrMT.fa $prefix.mt_contig.circlator.fa
$mummer_dir/delta-filter -m  $prefix.mt_improved.chrMT.delta > $prefix.mt_improved.chrMT.delta_filter
$mummer_dir/mummerplot --large --postscript $prefix.mt_improved.chrMT.delta_filter -p $prefix.mt_improved.chrMT.filter
perl $LRSDAY_HOME/scripts/fine_tune_gnuplot.pl -i $prefix.mt_improved.chrMT.filter.gp -o $prefix.mt_improved.chrMT.filter_adjust.gp -r ref.chrMT.fa -q $prefix.mt_contig.final.fa
$gnuplot_dir/gnuplot < $prefix.mt_improved.chrMT.filter_adjust.gp

# clean up intermediate files
if [[ $debug == "no" ]]
then
    rm *.filter.fplot
    rm *.filter.rplot
    rm *.delta
    rm *.delta_filter
    rm $prefix.mt_contig.*
    rm $prefix.non_mt_contig.*
    rm $prefix.mt_improved.chrMT.filter.ps
    rm $prefix.mt_improved.chrMT.filter.gp
    rm $prefix.mt_improved.chrMT.filter_adjust.gp
    rm ref.chrMT.list
    rm ref.chrMT.fa
    rm non_primary_contig.list
    rm $prefix.primary_contig.fa
    rm $prefix.non_primary_contig.fa
    if [ -e "$prefix.for_fixstart.for_skip.fa" ]
    then
	rm $prefix.for_fixstart.for_skip.fa
	rm $prefix.for_fixstart.for_skip.list
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
