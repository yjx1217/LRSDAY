#!/bin/bash
set -e 
#######################################
# load environment variables for LRSDAY
source ./../../env.sh

#######################################
# set project-specific variables
genome="./../06.Supervised_Final_Assembly/SK1.final.fa" # path of the input genome assembly
prefix="SK1" # file name prefix for the output files
threads=1 # number of threads to use
debug="no" # use "yes" if prefer to keep intermediate files, otherwise use "no".

#######################################
# process the pipeline

ln -s $genome genome.fa
$repeatmasker_dir/RepeatMasker -pa $threads -lib $LRSDAY_HOME/data/TY_lib.Yue_et_al_2017_NG.fa  genome.fa  -xsmall  -gff
$reannotate_dir/REannotate_longname  -g -f $LRSDAY_HOME/data/fuzzy_defragmentation.txt  -k $clustalw_dir/clustalw2 -r 0.48764  -d 10000  -t genome.fa.out  genome.fa
mv REannotate_output ${prefix}_REannotate_out

cd ${prefix}_REannotate_out 

perl $LRSDAY_HOME/scripts/parse_REannotate_gff.pl -o $prefix.REannotate.raw.gff
mv $prefix.REannotate.raw.gff $prefix.REannotate.gff
perl $LRSDAY_HOME/scripts/parse_REannotate_out.pl -i genome.fa.REannotation -p $prefix.TY_REannotate -r ./../genome.fa

cat $prefix.TY_REannotate.complete.raw.fa  |egrep ">" |egrep  '(TY1|TY2)' |sed "s/>//g" >$prefix.TY_REannotate.complete.TY1TY2.raw.list
perl $LRSDAY_HOME/scripts/select_fasta_by_list.pl -i $prefix.TY_REannotate.complete.raw.fa -l $prefix.TY_REannotate.complete.TY1TY2.raw.list  -o $prefix.TY_REannotate.complete.TY1TY2.raw.fa

cat $prefix.TY_REannotate.truncated.raw.fa  |egrep ">" |egrep  '(TY1|TY2)' |sed "s/>//g" >$prefix.TY_REannotate.truncated.TY1TY2.raw.list
perl $LRSDAY_HOME/scripts/select_fasta_by_list.pl -i $prefix.TY_REannotate.truncated.raw.fa -l $prefix.TY_REannotate.truncated.TY1TY2.raw.list  -o $prefix.TY_REannotate.truncated.TY1TY2.raw.fa

for i in {3..5}
do
    cat $prefix.TY_REannotate.complete.raw.fa  |egrep ">" |egrep "TY${i}" |sed "s/>//g"  >$prefix.TY_REannotate.complete.TY${i}.final.list
    cat $prefix.TY_REannotate.truncated.raw.fa  |egrep ">" |egrep "TY${i}" |sed "s/>//g"  >$prefix.TY_REannotate.truncated.TY${i}.final.list
done

TY2_query="$LRSDAY_HOME/data/TY2_specific_region.fa"

db="$prefix.TY_REannotate.complete.TY1TY2.raw.fa"
db_tag="complete_TY1TY2_db"
$blast_dir/makeblastdb  -in $db  -dbtype nucl -title $db_tag -hash_index -out $db_tag
$blast_dir/blastn -task blastn -query $TY2_query -num_threads $threads -db $db_tag -outfmt 7 >$prefix.$db_tag.blastn.fmt7.out
perl $LRSDAY_HOME/scripts/filter_blast_result.pl -i $prefix.$db_tag.blastn.fmt7.out -pct_identity_cutoff 95 -aln_length_cutoff 500 -o $prefix.$db_tag.blastn.fmt7.I95L500.out
cat $prefix.$db_tag.blastn.fmt7.I95L500.out |egrep -v "^#" |cut -f 2  > $prefix.TY_REannotate.complete.TY2.list
comm -23 <(sort $prefix.TY_REannotate.complete.TY1TY2.raw.list) <(sort $prefix.TY_REannotate.complete.TY2.list) > $prefix.TY_REannotate.complete.TY1.list
cat $prefix.TY_REannotate.complete.TY1.list | sed "s/TY2/TY1/g" >$prefix.TY_REannotate.complete.TY1.final.list
cat $prefix.TY_REannotate.complete.TY2.list | sed "s/TY1/TY2/g" >$prefix.TY_REannotate.complete.TY2.final.list

db="$prefix.TY_REannotate.truncated.TY1TY2.raw.fa"
db_tag="truncated_TY1TY2_db"
$blast_dir/makeblastdb  -in $db  -dbtype nucl -title $db_tag -hash_index -out $db_tag
$blast_dir/blastn -task blastn -query $TY2_query -num_threads $threads -db $db_tag -outfmt 7 >$prefix.$db_tag.blastn.fmt7.out
perl $LRSDAY_HOME/scripts/filter_blast_result.pl -i $prefix.$db_tag.blastn.fmt7.out -pct_identity_cutoff 95 -aln_length_cutoff 100 -o $prefix.$db_tag.blastn.fmt7.I95L100.out
cat $prefix.$db_tag.blastn.fmt7.I95L100.out |egrep -v "^#" |cut -f 2 > $prefix.TY_REannotate.truncated.TY2.list
comm -23 <(sort $prefix.TY_REannotate.truncated.TY1TY2.raw.list) <(sort $prefix.TY_REannotate.truncated.TY2.list) > $prefix.TY_REannotate.truncated.TY1.list
cat $prefix.TY_REannotate.truncated.TY1.list |sed "s/TY2/TY1/g" >$prefix.TY_REannotate.truncated.TY1.final.list
cat $prefix.TY_REannotate.truncated.TY2.list |sed "s/TY1/TY2/g" >$prefix.TY_REannotate.truncated.TY2.final.list

db="$LRSDAY_HOME/data/TY_lib.Yue_et_al_2017_NG.LTRonly.fa"
db_tag="soloLTR_db";
LTR_query="$prefix.TY_REannotate.soloLTR.raw.fa"
$blast_dir/makeblastdb  -in $db  -dbtype nucl -title $db_tag -hash_index -out $db_tag
$blast_dir/blastn  -task blastn -query $LTR_query -num_threads $threads -db $db_tag -outfmt 7 >$prefix.$db_tag.soloLTR.blastn.fmt7.out
perl $LRSDAY_HOME/scripts/trim_soloLTR_by_blast.pl -q $LTR_query -b  $prefix.$db_tag.soloLTR.blastn.fmt7.out -p $prefix  -i 70 -l 100
$bedtools_dir/bedtools sort  -i  $prefix.TY_REannotate.soloLTR.refined.gff > $prefix.TY_REannotate.soloLTR.refined.sorted.gff
perl $LRSDAY_HOME/scripts/rm_overlap_features_from_gff_simple.pl  -r ./../genome.fa  -i $prefix.TY_REannotate.soloLTR.refined.sorted.gff -o $prefix.TY_soloLTR.refined.nr.gff 

for i in {1..5}
do
    perl $LRSDAY_HOME/scripts/TY_list2gff3.pl -i $prefix.TY_REannotate.complete.TY${i}.final.list -o $prefix.TY_REannotate.complete.TY${i}.final.gff -t $prefix
    perl $LRSDAY_HOME/scripts/TY_list2gff3.pl -i $prefix.TY_REannotate.truncated.TY${i}.final.list -o $prefix.TY_REannotate.truncated.TY${i}.final.gff -t $prefix
done

cat $prefix.TY_REannotate.*.final.gff > $prefix.TY.complete_plus_truncated.final.gff

$bedtools_dir/bedtools sort -i $prefix.TY.complete_plus_truncated.final.gff > $prefix.TY.complete_plus_truncated.final.sorted.gff 
$bedtools_dir/bedtools intersect -v -a $prefix.TY_soloLTR.refined.nr.gff -b $prefix.TY.complete_plus_truncated.final.sorted.gff  >$prefix.TY.soloLTR.final.gff

cat $prefix.TY.complete_plus_truncated.final.gff $prefix.TY.soloLTR.final.gff > $prefix.TY.all.final.gff

perl $LRSDAY_HOME/scripts/tidy_maker_gff3.pl -r ./../genome.fa -i $prefix.TY.all.final.gff -o $prefix.TE.gff3 -t $prefix 

cp $prefix.TE.gff3 ./../
cp $prefix.REannotate.gff ./../
rm $prefix.*.final.gff
cd ..

# clean up intermediate files
if [[ $debug == "no" ]]
then
    rm genome.fa*
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
