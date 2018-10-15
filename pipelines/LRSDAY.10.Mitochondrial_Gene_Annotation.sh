#!/bin/bash
set -e -o pipefail

#######################################
# load environment variables for LRSDAY
source ./../../env.sh
PERL5LIB="$PERL5LIB:$pirobject_dir/lib"
export RNAFINDER_CFG_PATH="$rnafinder_dir"
export MF2SQN_LIB="$mf2sqn_dir/lib"
export MFANNOT_LIB_PATH="$mfannot_data_dir/protein_collections"
export MFANNOT_EXT_CFG_PATH="$mfannot_data_dir/config"
export MFANNOT_MOD_PATH="$mfannot_data_dir/models"
export BLASTMAT="$blast_matrices_dir"
export EGC="$mfannot_data_dir/EGC"
export ERPIN_MOD_PATH="$mfannot_data_dir/models/Erpin_models"
export PIR_DATAMODEL_PATH="$pirobject_dir/PirModels"
export PATH="$flip_dir:$blast_dir:$muscle_dir:$umac_dir:$hmmer_dir:$erpin_dir:$tbl2asn_dir:$pirobject_dir:$pirmodels_dir:$hmmsearchwc_dir:$exonerate_dir:$emboss_dir:$mf2sqn_dir:$mf2sqn_dir:$grab_fasta_dir:$rnafinder_dir:$mfannot_dir:$PATH"

#######################################
# set project-specific variables
genome="./../07.Supervised_Final_Assembly/SK1.assembly.final.fa" # The file path of the input genome assembly.
chrMT_tag="chrMT" # The sequence name for the mitochondrial genome in the input genome assembly, if there are multiple corresponding contigs/scaffolds, use a single ';' to separate them. e.g. "chrMT_1;chrMT_2". Default = "chrMT". 
genetic_code_table=3 # The NCBI genetic code table (https://www.ncbi.nlm.nih.gov/Taxonomy/Utils/wprintgc.cgi) for the annotated mitochondrial genome. Default = 3 (i.e. Yeast Mitochondria)

prefix="SK1" # The file name prefix for the output files.
debug="no" # Whehter to keep intermediate files for debugging. Use "yes" if prefer to keep intermediate files, otherwise use "no". Default = "no".

######################################
# process the pipeline

echo $chrMT_tag | sed -e "s/;/\n/g" > $prefix.assembly.chrMT.list

#$LRSDAY_HOME/scripts/select_fasta_by_list.pl -i $genome -l $prefix.assembly.chrMT.list -m reverse -o $prefix.assembly.nuclear_genome.fa
$LRSDAY_HOME/scripts/select_fasta_by_list.pl -i $genome -l $prefix.assembly.chrMT.list -m normal -o $prefix.assembly.mitochondrial_genome.fa

mkdir tmp

$mfannot_dir/mfannot \
    --genetic $genetic_code_table \
    --outputfile $prefix.mitochondrial_genome.mfannot.out \
    --logfile $prefix.mitochondrial_genome.mfannot.log \
    --T $(pwd)/tmp \
    $prefix.assembly.mitochondrial_genome.fa

perl $LRSDAY_HOME/scripts/mfannot2gff3.pl -i $prefix.mitochondrial_genome.mfannot.out -o $prefix.mitochondrial_genome.mfannot.gff3 -m lite
perl $LRSDAY_HOME/scripts/extract_cds_from_tidy_gff3.pl -r $prefix.assembly.mitochondrial_genome.fa -g $prefix.mitochondrial_genome.mfannot.gff3 -o $prefix.mitochondrial_genome.mfannot.cds.fa
perl $LRSDAY_HOME/scripts/cds2protein.pl -i $prefix.mitochondrial_genome.mfannot.cds.fa -t $genetic_code_table -p $prefix.mitochondrial_genome.mfannot

perl $LRSDAY_HOME/scripts/prepare_PoFFgff_simple.pl -i $prefix.mitochondrial_genome.mfannot.gff3 -o $prefix.mitochondrial_genome.mfannot.PoFF.gff 
perl $LRSDAY_HOME/scripts/prepare_PoFFfaa_simple.pl -i $prefix.mitochondrial_genome.mfannot.trimmed_cds.fa -o $prefix.mitochondrial_genome.mfannot.PoFF.ffn 
perl $LRSDAY_HOME/scripts/prepare_PoFFfaa_simple.pl -i $prefix.mitochondrial_genome.mfannot.pep.fa -o $prefix.mitochondrial_genome.mfannot.PoFF.faa

# clean up intermediate files
if [[ $debug == "no" ]]
then
    rm -r tmp
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
