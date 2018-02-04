#!/bin/bash
set -e -o pipefail
#######################################
# load environment variables for LRSDAY
source ./../../env.sh

#######################################
# set project-specific variables
genome_assembly="./../06.Supervised_Final_Assembly/SK1.final.fa" # path of the input genome assembly
genome_tag="SK1" # file name prefix for the output files
threads=1 # number of threads to use
maker_opts="$LRSDAY_HOME/misc/maker_opts.customized.ctl" # configuration file for MAKER. You can edit this file if you have native transciptome/EST data for the strain/species that you sequenced or if you want to adapt it to annotate other eukaryotic organisms. Otherwise, please keep it unchanged. Please note that if this file is in the same directory where this bash script is executed, the file name cannot be "maker_opts.ctl".
EVM_weights="$LRSDAY_HOME/misc/EVM_weights.customized.txt" # configuration file for EVM. A list of numeric weight values to be applied to each type of evidence.
debug="no" # use "yes" if prefer to keep intermediate files, otherwise use "no".

#######################################
# process the pipeline
echo "genome_assembly=$genome_assembly"
echo "genome_tag=$genome_tag"

# convert the genome assembly file to all uppercases
perl $LRSDAY_HOME/scripts/switch_letter_cases_in_fasta.pl -i $genome_assembly -o genome.uppercase.fa -c upper

cp $maker_opts maker_opts.ctl
cp $LRSDAY_HOME/misc/maker_exe.ctl .
cp $LRSDAY_HOME/misc/maker_bopts.ctl .
cp $LRSDAY_HOME/misc/maker_evm.ctl .

$maker_dir/maker -fix_nucleotides -genome genome.uppercase.fa -cpus $threads -base $genome_tag
# $maker_dir/fasta_merge -d $genome_tag.maker.output/${genome_tag}_master_datastore_index.log -o $genome_tag.maker.fasta
$maker_dir/gff3_merge -d $genome_tag.maker.output/${genome_tag}_master_datastore_index.log -n -g -o $genome_tag.maker.raw.gff3
perl $LRSDAY_HOME/scripts/collect_maker_evidences.pl -t $genome_tag

cat $genome_tag.maker.raw.gff3 |egrep "=trnascan" > $genome_tag.maker.raw.tRNA.gff3
# cat $genome_tag.maker.raw.gff3 |egrep "=snoscan" > $genome_tag.maker.raw.snoscan.gff3 # developmental feature
cat $genome_tag.maker.raw.gff3 |egrep -v "=trnascan"|egrep -v "=snoscan" > $genome_tag.maker.raw.protein_coding_gene.gff3

# use EVM to further polishing the annotation for multi-exon genes
perl $LRSDAY_HOME/scripts/filter_gff3_for_single_exon_genes.pl -i $genome_tag.maker.raw.protein_coding_gene.gff3 -o $genome_tag.maker.raw.protein_coding_gene.single_exon_gene.gff3

$bedtools_dir/bedtools intersect -v  -a $genome_tag.maker.raw.protein_coding_gene.gff3 \
		       -b $genome_tag.maker.raw.protein_coding_gene.single_exon_gene.gff3 \
		       > $genome_tag.maker.raw.protein_coding_gene.multiple_exon_gene.gff3
$bedtools_dir/bedtools intersect -v  -a $genome_tag.protein_evidence.gff3 \
		       -b $genome_tag.maker.raw.protein_coding_gene.single_exon_gene.gff3 \
		       > $genome_tag.protein_evidence.for_gene_model_refinement.gff3
$bedtools_dir/bedtools intersect -v  -a $genome_tag.est_evidence.gff3 \
		       -b $genome_tag.maker.raw.protein_coding_gene.single_exon_gene.gff3 \
		       > $genome_tag.est_evidence.for_gene_model_refinement.gff3

perl $LRSDAY_HOME/scripts/exonerate2gene_maker.pl -i $genome_tag.protein_evidence.for_gene_model_refinement.gff3 -o $genome_tag.protein_evidence.complementary_gene_model.gff3 
perl $LRSDAY_HOME/scripts/exonerate2gene_maker.pl -i $genome_tag.est_evidence.for_gene_model_refinement.gff3 -o $genome_tag.est_evidence.complementary_gene_model.gff3

cat $genome_tag.maker.raw.protein_coding_gene.gff3 $genome_tag.protein_evidence.complementary_gene_model.gff3 $genome_tag.est_evidence.complementary_gene_model.gff3 > $genome_tag.maker.combined.gff3

mkdir $genome_tag.EVM.output
cd $genome_tag.EVM.output
$EVM_HOME/EvmUtils/partition_EVM_inputs.pl \
    --genome ./../$genome_assembly \
    --gene_predictions ./../$genome_tag.maker.combined.gff3 \
    --protein_alignments ./../$genome_tag.protein_evidence.gff3 \
    --transcript_alignments ./../$genome_tag.est_evidence.gff3 \
    --segmentSize 100000 \
    --overlapSize 10000 \
    --partition_listing $genome_tag.partitions_list.out

$EVM_HOME/EvmUtils/write_EVM_commands.pl --genome ./../$genome_assembly  --weights $EVM_weights  \
					 --gene_predictions ./../$genome_tag.maker.combined.gff3 \
					 --output_file_name $genome_tag.evm.out \
					 --partitions $genome_tag.partitions_list.out > $genome_tag.commands.list

$EVM_HOME/EvmUtils/execute_EVM_commands.pl $genome_tag.commands.list | tee $genome_tag.EVM_run.log
$EVM_HOME/EvmUtils/recombine_EVM_partial_outputs.pl --partitions $genome_tag.partitions_list.out --output_file_name $genome_tag.evm.out
$EVM_HOME/EvmUtils/convert_EVM_outputs_to_GFF3.pl  --partitions $genome_tag.partitions_list.out --output $genome_tag.evm.out  --genome ./../$genome_assembly

perl $LRSDAY_HOME/scripts/collect_EVM_gff3.pl -p $genome_tag  -r ./../$genome_assembly

cat ${genome_tag}.EVM.raw.gff3 ./../$genome_tag.maker.raw.tRNA.gff3 > ${genome_tag}.EVM.raw.with_tRNA.gff3
perl $LRSDAY_HOME/scripts/tidy_maker_gff3.pl -i ${genome_tag}.EVM.raw.with_tRNA.gff3 -r ./../$genome_assembly -t $genome_tag -o ${genome_tag}.EVM.gff3
cp ${genome_tag}.EVM.gff3 ./../
cd ..
perl $LRSDAY_HOME/scripts/extract_cds_from_tidy_gff3.pl -r ${genome_assembly} -g ${genome_tag}.EVM.gff3 -o ${genome_tag}.EVM.cds.fa
perl $LRSDAY_HOME/scripts/cds2protein.pl -i ${genome_tag}.EVM.cds.fa -p ${genome_tag}.EVM

perl $LRSDAY_HOME/scripts/prepare_PoFFgff_simple.pl -i ${genome_tag}.EVM.gff3 -o ${genome_tag}.EVM.PoFF.gff 
perl $LRSDAY_HOME/scripts/prepare_PoFFfaa_simple.pl -i ${genome_tag}.EVM.trimmed_cds.fa -o ${genome_tag}.EVM.PoFF.ffn 
perl $LRSDAY_HOME/scripts/prepare_PoFFfaa_simple.pl -i ${genome_tag}.EVM.pep.fa -o ${genome_tag}.EVM.PoFF.faa

# filtered out snoRNA annotation since it is still an experimental features suffering from redundant annotations
cp ${genome_tag}.EVM.gff3 ${genome_tag}.EVM.gff3.tmp
cat ${genome_tag}.EVM.gff3.tmp | egrep -v "snoRNA" > ${genome_tag}.EVM.gff3

# clean up intermediate files
if [[ $debug == "no" ]]
then
    rm genome.uppercase.fa
    rm $genome_tag.maker.raw.tRNA.gff3
    rm $genome_tag.maker.raw.protein_coding_gene.gff3
    rm $genome_tag.maker.raw.protein_coding_gene.single_exon_gene.gff3
    rm $genome_tag.maker.raw.protein_coding_gene.multiple_exon_gene.gff3
    rm $genome_tag.protein_evidence.for_gene_model_refinement.gff3
    rm $genome_tag.est_evidence.for_gene_model_refinement.gff3
    rm $genome_tag.protein_evidence.complementary_gene_model.gff3
    rm $genome_tag.est_evidence.complementary_gene_model.gff3
    rm $genome_tag.EVM.gff3.tmp
    rm $genome_tag.maker.combined.gff3
    rm -rf _Inline
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
