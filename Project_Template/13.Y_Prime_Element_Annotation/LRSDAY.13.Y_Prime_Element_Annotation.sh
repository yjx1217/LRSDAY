#!/bin/bash
set -e 
#######################################
# load environment variables for LRSDAY
source ./../../env.sh

#######################################
# set project-specific variables
prefix="SK1" # The file name prefix for the processing sample. Default = "SK1" for the testing example.
genome="./../07.Supervised_Final_Assembly/$prefix.assembly.final.fa" # The file path of the input genome assembly.
debug="no" # Whether to keep intermediate files for debugging. Use "yes" if prefer to keep intermediate files, otherwise use "no". Default = "no".

#######################################
# process the pipeline
feature_type="Y_prime_element"
query="$LRSDAY_HOME/data/query.Y_prime_element.long.fa"
length_cutoff_for_completeness=3500

$ucsc_dir/blat -maxIntron=1000 $genome $query  $prefix.$feature_type.blat.psl
$ucsc_dir/pslCDnaFilter  -minId=0.9 -minAlnSize=1000 -bestOverlap -filterWeirdOverlapped   $prefix.$feature_type.blat.psl  $prefix.$feature_type.blat.filtered.psl
$ucsc_dir/pslScore  $prefix.$feature_type.blat.filtered.psl | sort -nk5 -r > $prefix.$feature_type.blat.filtered.pslScore.out
perl $LRSDAY_HOME/scripts/psl2gff3.pl -i $prefix.$feature_type.blat.filtered.psl -o $prefix.$feature_type.raw.gff3 -ft $feature_type -t $prefix -r $genome -l $length_cutoff_for_completeness
perl $LRSDAY_HOME/scripts/tidy_maker_gff3.pl -i $prefix.$feature_type.raw.gff3 -r $genome -o $prefix.$feature_type.gff3 -t $prefix
perl $LRSDAY_HOME/scripts/gff2seq_simple.pl -r $genome -g $prefix.$feature_type.gff3 -o $prefix.$feature_type.fa
$muscle_dir/muscle -in $prefix.$feature_type.fa -out $prefix.$feature_type.aln.fa

# clean up intermediate files
if [[ $debug == "no" ]]
then
    rm $prefix.$feature_type.raw.gff3
    rm $prefix.Y_prime_element.blat.psl
    rm $prefix.Y_prime_element.blat.filtered.psl
    rm $prefix.Y_prime_element.blat.filtered.pslScore.out
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
