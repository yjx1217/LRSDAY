#!/bin/bash
set -e -o pipefail
#######################################
# load environment variables for LRSDAY
source ./../../env.sh

#######################################
# set project-specific variables
prefix="CPG_1a" # The file name prefix (only allowing strings of alphabetical letters, numbers, and underscores) for the processing sample. Default = "CPG_1a" for the testing example.          
genome_assembly="./../07.Supervised_Final_Assembly/$prefix.assembly.final.fa" # The file path of the final input genome assembly.
chrMT_tag="chrMT" # The sequence name for the mitochondrial genome in the final assembly. If there are multiple sequences, use a single ';' to separate them. e.g. "chrMT_part1;chrMT_part2". Default = "chrMT".
debug="no" # Whether to keep intermediate files for debugging. Use "yes" if prefer to keep intermediate files, otherwise use "no". Default = "no".

#######################################
# process the pipeline
feature_type="X_element"
length_cutoff_for_completeness=300

echo $chrMT_tag | sed -e "s/;/\n/g" > $prefix.assembly.chrMT.list
perl $LRSDAY_HOME/scripts/select_fasta_by_list.pl -i $genome_assembly -l $prefix.assembly.chrMT.list -m reverse -o $prefix.assembly.nuclear_genome.fa
perl $LRSDAY_HOME/scripts/tidy_fasta.pl -i $prefix.assembly.nuclear_genome.fa -o $prefix.assembly.nuclear_genome.tidy.fa


$hmmer_dir/nhmmer -E 1 --tblout $prefix.$feature_type.nhmmer.out $LRSDAY_HOME/data/S288C.$feature_type.hmm $prefix.assembly.nuclear_genome.tidy.fa
perl $LRSDAY_HOME/scripts/nhmer2seq.pl -r $prefix.assembly.nuclear_genome.tidy.fa -i $prefix.$feature_type.nhmmer.out -e 0.0001 -p $prefix -ft $feature_type -l $length_cutoff_for_completeness
perl $LRSDAY_HOME/scripts/tidy_maker_gff3.pl -i $prefix.$feature_type.raw.gff3 -r $prefix.assembly.nuclear_genome.tidy.fa -o $prefix.nuclear_genome.$feature_type.gff3 -t $prefix
perl $LRSDAY_HOME/scripts/gff2seq_simple.pl -r $prefix.assembly.nuclear_genome.tidy.fa -g $prefix.nuclear_genome.$feature_type.gff3 -o $prefix.nuclear_genome.$feature_type.fa 
$muscle_dir/muscle -in $prefix.nuclear_genome.$feature_type.fa -out $prefix.nuclear_genome.$feature_type.aln.fa

# clean up intermediate files
if [[ $debug == "no" ]]
then
    rm $prefix.$feature_type.raw.gff3
    rm $prefix.$feature_type.raw.fa
    rm $prefix.assembly.nuclear_genome.fa
    rm $prefix.assembly.nuclear_genome.tidy.fa
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
