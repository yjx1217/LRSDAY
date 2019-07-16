#!/bin/bash
set -e -o pipefail
##########################################
# load environment variables for LRSDAY
source ./../../env.sh
PATH=$gnuplot_dir:$PATH

###########################################
query_genome="SK1.assembly.raw.fa"
ref_genome="S288C.ASM205763v1.fa"
prefix="SK1.assembly.raw"
color_scheme="by_match_direction" # The color scheme used for the dotplots: "by_match_direction" or "by_match_identity". Default = "by_match_direction".
threads=4
debug="no" # Whether to keep intermediate files for debugging. Default = "no".

###########################################
# process the pipeline

# make the comparison between the assembled genome and the reference genome
$mummer4_dir/nucmer -t $threads --maxmatch --nosimplify  -p $prefix.mummerplot  $ref_genome $query_genome
$mummer4_dir/delta-filter -m  $prefix.mummerplot.delta > $prefix.mummerplot.delta_filter
if [[ $color_scheme == "by_match_direction" ]]
then
    $mummer4_dir/mummerplot --large --postscript $prefix.mummerplot.delta_filter -p $prefix.mummerplot.filter
else
    $mummer4_dir/mummerplot --large --postscript --color $prefix.mummerplot.delta_filter -p $prefix.mummerplot.filter
fi
perl $LRSDAY_HOME/scripts/fine_tune_gnuplot.pl -i $prefix.mummerplot.filter.gp -o $prefix.mummerplot.filter_adjust.gp -r $ref_genome -q $query_genome
$gnuplot_dir/gnuplot < $prefix.mummerplot.filter_adjust.gp

if [[ $debug == "no" ]]
then
    rm *.delta
    rm *.delta_filter
    rm *.filter.fplot
    rm *.filter.rplot
    rm *.filter.gp
    rm *.filter_adjust.gp
    rm *.filter.ps
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


