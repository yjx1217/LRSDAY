#!/bin/bash
set -e -o pipefail
#######################################
# load environment variables for LRSDAY
source ./../../env.sh

#######################################
# set project-specific variables
prefix="SK1" # file name prefix for the output files
threads=1 # number of threads to use
input_gff="./../08.Gene_Annotation/SK1.EVM.gff3" # path of the input maker gff3 file generated in the task 08.Gene_Annotation
query_PoFF_faa="./../08.Gene_Annotation/SK1.EVM.PoFF.faa" # path of the PoFF.faa file generated in the task 08.Gene_Annotation
query_PoFF_gff="./../08.Gene_Annotation/SK1.EVM.PoFF.gff" # path of the PoFF.gff file generated in the task 08.Gene_Annotation
ref_PoFF_faa="$LRSDAY_HOME/data/SGDref.PoFF.faa" # path of the reference proteome file in FASTA format: for S. cerevisiae and its close relatives, you can directly use the pre-shipped file: SGDref.PoFF.faa; if you work with other organisms, you can check ProteinOrtho's manual for details on how to prepare such file.
ref_PoFF_gff="$LRSDAY_HOME/data/SGDref.PoFF.gff" # path of the reference gene GFF file in GFF format: for S. cerevisiae and its close relatives, you can directly use the pre-shipped file: SGDref.PoFF.gff; if you work with other organisms, you can check ProteinOrtho's manual for details on how to prepare such file.
debug="no" # use "yes" if prefer to keep intermediate files, otherwise use "no".

#######################################
# process the pipeline
cp $query_PoFF_faa query.PoFF.faa
cp $query_PoFF_gff query.PoFF.gff
cp $ref_PoFF_faa ref.PoFF.faa
cp $ref_PoFF_gff ref.PoFF.gff
output_gff="$prefix.orthology_map.gff3"

test_file_existence () {
    filename=$1
    if [[ ! -f $filename ]]
    then
	echo "the file $filename does not exists! process terminated!"
	exit
    fi
}

test_file_existence $query_PoFF_faa
test_file_existence $query_PoFF_gff
test_file_existence $ref_PoFF_faa
test_file_existence $ref_PoFF_gff

proteinortho_input="query.PoFF.faa ref.PoFF.faa"
perl $proteinortho_dir/proteinortho5.pl -cpus=$threads -blastpath=$blast_dir -singles -synteny -project=$prefix  $proteinortho_input  
perl $LRSDAY_HOME/scripts/update_gff3_by_proteinortho.pl -i $input_gff -x $prefix.poff -r ref.PoFF.faa -q query.PoFF.faa -o $prefix.updated.gff3

# clean up intermediate files
if [[ $debug == "no" ]]
then
    rm query.*
    rm ref.*
    rm $prefix.ffadj-graph
    rm $prefix.blast-graph
    rm $prefix.proteinortho-graph
    rm $prefix.poff-graph
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
