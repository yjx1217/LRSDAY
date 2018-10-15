#!/bin/bash
set -e -o pipefail
#######################################
# load environment variables for LRSDAY
source ./../../env.sh

#######################################
# set project-specific variables
prefix="SK1" # The file name prefix for the output files.
threads=1 # The number of threads to use. Default = "1".

input_nuclear_gene_gff="./../09.Nuclear_Gene_Annotation/SK1.nuclear_genome.EVM.gff3" # The file path of the input nuclear gene gff3 file generated in the task 09.Nuclear_Gene_Annotation. Set this variable as well as the following two variables to "" if you want to skip this step for the nuclear gene annotation. 
query_nuclear_gene_PoFF_faa="./../09.Nuclear_Gene_Annotation/SK1.nuclear_genome.EVM.PoFF.faa" # The file path of the PoFF.faa file generated in the task 09.Nuclear_Gene_Annotation.
query_nuclear_gene_PoFF_gff="./../09.Nuclear_Gene_Annotation/SK1.nuclear_genome.EVM.PoFF.gff" # The file path of the PoFF.gff file generated in the task 09.Nuclear_Gene_Annotation.

input_mitochondrial_gene_gff="./../10.Mitochondrial_Gene_Annotation/SK1.mitochondrial_genome.mfannot.gff3" # The file path of the input mitochondrial gene gff3 file generated in the task 10.Mitochondrial_Gene_Annotation. Set this variable as well as the following two variables to "" if you want to skip this step for the mitochondrial gene annotation. 
query_mitochondrial_gene_PoFF_faa="./../10.Mitochondrial_Gene_Annotation/SK1.mitochondrial_genome.mfannot.PoFF.faa" # The file path of the PoFF.faa file generated in the task 10.Mitochondrial_Gene_Annotation.
query_mitochondrial_gene_PoFF_gff="./../10.Mitochondrial_Gene_Annotation/SK1.mitochondrial_genome.mfannot.PoFF.gff" # The file path of the PoFF.gff file generated in the task 10.Mitochondrial_Gene_Annotation.

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
if [[ ! -z $input_nuclear_gene_gff ]]
then 
    echo "The variable input_nuclear_gene_gff has been specified, so expecting input files specified for query_nuclear_gene_PoFF_faa and query_nuclear_gene_PoFF_gff"
    test_file_existence $query_nuclear_gene_PoFF_faa
    test_file_existence $query_nuclear_gene_PoFF_gff
    cp $query_nuclear_gene_PoFF_faa query.nuclear_gene.PoFF.faa
    cp $query_nuclear_gene_PoFF_gff query.nuclear_gene.PoFF.gff
    output_nuclear_gene_gff="$prefix.nuclear_gene.orthology_map.gff3"
    proteinortho_nuclear_gene_input="query.nuclear_gene.PoFF.faa ref.PoFF.faa"
    perl $proteinortho_dir/proteinortho5.pl -cpus=$threads -blastpath=$blast_dir -singles -synteny -project="$prefix.nuclear_gene" $proteinortho_nuclear_gene_input  
    perl $LRSDAY_HOME/scripts/update_gff3_by_proteinortho.pl -i $input_nuclear_gene_gff -x $prefix.nuclear_gene.poff -r ref.PoFF.faa -q query.nuclear_gene.PoFF.faa -o $prefix.nuclear_gene.updated.gff3

    if [[ $debug == "no" ]]
    then
	rm $prefix.nuclear_gene.ffadj-graph
	rm $prefix.nuclear_gene.blast-graph
	rm $prefix.nuclear_gene.proteinortho-graph
	rm $prefix.nuclear_gene.poff-graph
    fi
fi

# process mitochondrial gene annotation
if [[ ! -z $input_mitochondrial_gene_gff ]]
then 
    echo "The variable input_mitochondrial_gene_gff has been specified, so expecting input files specified for query_mitochondrial_gene_PoFF_faa and query_mitochondrial_gene_PoFF_gff"
    test_file_existence $query_mitochondrial_gene_PoFF_faa
    test_file_existence $query_mitochondrial_gene_PoFF_gff
    cp $query_mitochondrial_gene_PoFF_faa query.mitochondrial_gene.PoFF.faa
    cp $query_mitochondrial_gene_PoFF_gff query.mitochondrial_gene.PoFF.gff
    output_mitochondrial_gene_gff="$prefix.mitochondrial_gene.orthology_map.gff3"
    proteinortho_mitochondrial_gene_input="query.mitochondrial_gene.PoFF.faa ref.PoFF.faa"
    perl $proteinortho_dir/proteinortho5.pl -cpus=$threads -blastpath=$blast_dir -singles -synteny -project="$prefix.mitochondrial_gene"  $proteinortho_mitochondrial_gene_input  
    perl $LRSDAY_HOME/scripts/update_gff3_by_proteinortho.pl -i $input_mitochondrial_gene_gff -x $prefix.mitochondrial_gene.poff -r ref.PoFF.faa -q query.mitochondrial_gene.PoFF.faa -o $prefix.mitochondrial_gene.updated.gff3

    if [[ $debug == "no" ]]
    then
	rm $prefix.mitochondrial_gene.ffadj-graph
	rm $prefix.mitochondrial_gene.blast-graph
	rm $prefix.mitochondrial_gene.proteinortho-graph
	rm $prefix.mitochondrial_gene.poff-graph
    fi
fi


# clean up intermediate files
if [[ $debug == "no" ]]
then
    rm query.*
    rm ref.*
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
