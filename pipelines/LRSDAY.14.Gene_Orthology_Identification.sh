#!/bin/bash
set -e -o pipefail
#######################################
# load environment variables for LRSDAY
source ./../../env.sh

#######################################
# set project-specific variables
prefix="CPG_1a" # The file name prefix (only allowing strings of alphabetical letters, numbers, and underscores) for the processing sample. Default = "CPG_1a" for the testing example.          

threads=8 # The number of threads to use. Default = "8".

input_nuclear_genome_gff="./../09.Nuclear_Gene_Annotation/$prefix.nuclear_genome.gff3" # The file path of the input nuclear gene gff3 file generated in the task 09.Nuclear_Gene_Annotation. Set this variable as well as the following two variables to "" if you want to skip this step for the nuclear gene annotation. 
query_nuclear_genome_PoFF_faa="./../09.Nuclear_Gene_Annotation/$prefix.nuclear_genome.PoFF.faa" # The file path of the PoFF.faa file generated in the task 09.Nuclear_Gene_Annotation.Po
query_nuclear_genome_PoFF_gff="./../09.Nuclear_Gene_Annotation/$prefix.nuclear_genome.PoFF.gff" # The file path of the PoFF.gff file generated in the task 09.Nuclear_Gene_Annotation.

input_mitochondrial_genome_gff="./../10.Mitochondrial_Gene_Annotation/$prefix.mitochondrial_genome.gff3" # The file path of the input mitochondrial gene gff3 file generated in the task 10.Mitochondrial_Gene_Annotation. Set this variable as well as the following two variables to "" if you want to skip this step for the mitochondrial gene annotation. 
query_mitochondrial_genome_PoFF_faa="./../10.Mitochondrial_Gene_Annotation/$prefix.mitochondrial_genome.PoFF.faa" # The file path of the PoFF.faa file generated in the task 10.Mitochondrial_Gene_Annotation.
query_mitochondrial_genome_PoFF_gff="./../10.Mitochondrial_Gene_Annotation/$prefix.mitochondrial_genome.PoFF.gff" # The file path of the PoFF.gff file generated in the task 10.Mitochondrial_Gene_Annotation.

ref_PoFF_faa="$LRSDAY_HOME/data/SGDref.PoFF.faa" # The file path of the reference proteome file in FASTA format: for S. cerevisiae and its close relatives, you can directly use the pre-shipped file: SGDref.PoFF.faa; if you work with other organisms, you can check ProteinOrtho's manual for details on how to prepare such file.
ref_PoFF_gff="$LRSDAY_HOME/data/SGDref.PoFF.gff" # The path of the reference gene GFF file in GFF format: for S. cerevisiae and its close relatives, you can directly use the pre-shipped file: SGDref.PoFF.gff; if you work with other organisms, you can check ProteinOrtho's manual for details on how to prepare such file.

debug="no" # Whether to keep intermediate files for debugging. Use "yes" if prefer to keep intermediate files, otherwise use "no". Default = "no".

#######################################
# process the pipeline

test_file_existence () {
    filename=$1
    if [[ ! -f $filename ]]
    then
	echo "The file $filename does not exists! Process terminated!"
	exit
    fi
}

test_file_existence $ref_PoFF_faa
test_file_existence $ref_PoFF_gff

cp $ref_PoFF_faa ref.PoFF.faa
cp $ref_PoFF_gff ref.PoFF.gff

# process nuclear gene annotation
if [[ ! -z $input_nuclear_genome_gff ]]
then 
    echo "The variable input_nuclear_genome_gff has been specified, so expecting input files specified for query_nuclear_genome_PoFF_faa and query_nuclear_genome_PoFF_gff"
    test_file_existence $query_nuclear_genome_PoFF_faa
    test_file_existence $query_nuclear_genome_PoFF_gff
    cp $query_nuclear_genome_PoFF_faa query.nuclear_genome.PoFF.faa
    cp $query_nuclear_genome_PoFF_gff query.nuclear_genome.PoFF.gff
    proteinortho_nuclear_genome_input="query.nuclear_genome.PoFF.faa ref.PoFF.faa"
    $proteinortho_dir/proteinortho -cpus=$threads -p=blastp -binpath=$blast_dir -singles -synteny -project="$prefix.nuclear_genome" $proteinortho_nuclear_genome_input
    perl $LRSDAY_HOME/scripts/update_gff3_by_proteinortho.pl -i $input_nuclear_genome_gff -x $prefix.nuclear_genome.proteinortho.tsv -r ref.PoFF.faa -q query.nuclear_genome.PoFF.faa -o $prefix.nuclear_genome.SGD_orthology_mapped.gff3


    if [[ $debug == "no" ]]
    then
	rm $prefix.nuclear_genome.ffadj-graph
	rm $prefix.nuclear_genome.blast-graph
	rm $prefix.nuclear_genome.proteinortho-graph
	rm $prefix.nuclear_genome.poff-graph
    fi
fi

# process mitochondrial gene annotation
if [[ ! -z $input_mitochondrial_genome_gff ]]
then 
    echo "The variable input_mitochondrial_genome_gff has been specified, so expecting input files specified for query_mitochondrial_genome_PoFF_faa and query_mitochondrial_genome_PoFF_gff"
    test_file_existence $query_mitochondrial_genome_PoFF_faa
    test_file_existence $query_mitochondrial_genome_PoFF_gff
    cp $query_mitochondrial_genome_PoFF_faa query.mitochondrial_genome.PoFF.faa
    cp $query_mitochondrial_genome_PoFF_gff query.mitochondrial_genome.PoFF.gff
    proteinortho_mitochondrial_genome_input="query.mitochondrial_genome.PoFF.faa ref.PoFF.faa"
    perl $proteinortho_dir/proteinortho -cpus=$threads -p=blastp -binpath=$blast_dir -singles -synteny -project="$prefix.mitochondrial_genome"  $proteinortho_mitochondrial_genome_input  
    perl $LRSDAY_HOME/scripts/update_gff3_by_proteinortho.pl -i $input_mitochondrial_genome_gff -x $prefix.mitochondrial_genome.proteinortho.tsv -r ref.PoFF.faa -q query.mitochondrial_genome.PoFF.faa -o $prefix.mitochondrial_genome.SGD_orthology_mapped.gff3

    if [[ $debug == "no" ]]
    then
	rm $prefix.mitochondrial_genome.ffadj-graph
	rm $prefix.mitochondrial_genome.blast-graph
	rm $prefix.mitochondrial_genome.proteinortho-graph
	rm $prefix.mitochondrial_genome.poff-graph
    fi
fi


# clean up intermediate files
if [[ $debug == "no" ]]
then
    rm query.*
    rm ref.*
    rm proteinortho_cache_${prefix}.log
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
