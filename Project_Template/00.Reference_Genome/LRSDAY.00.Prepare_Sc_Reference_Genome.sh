#!/bin/bash
set -e -o pipefail

#######################################
# load environment variables for LRSDAY
source ./../../env.sh

#######################################
# set project-specific variables

#######################################
# process the pipeline

cp $LRSDAY_HOME/data/S288C.ASM205763v1.fa.gz .
cp $LRSDAY_HOME/data/S288C.ASM205763v1.noncore_masked.fa.gz .
gunzip *.gz

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
