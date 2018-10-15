#!/bin/bash
set -e -o pipefail
#######################################
# load environment variables for LRSDAY
source ./../../env.sh

#######################################
# set project-specific variables
file_name="SK1.filtered_subreads.bam" # The name of the ENA bam file for the testing example. 
file_url="ftp://ftp.sra.ebi.ac.uk/vol1/ERZ448/ERZ448251/SK1.filtered_subreads.bam" # The URL of the ENA bam file for the testing example.
prefix="SK1" # The file name prefix for output files of the testing example. 

#######################################
# process the pipeline

echo "download the bam file from the ENA database ..."
wget $file_url
echo "bam2fastq ..."
$bedtools_dir/bedtools bamtofastq -i $file_name -fq $prefix.filtered_subreads.fastq 
echo "gzip fastq ..."
gzip $prefix.filtered_subreads.fastq
rm $file_name

cd pacbio_fofn_files
echo "download the metadata and raw PacBio reads in .h5 format ..."

wget ftp://ftp.sra.ebi.ac.uk/vol1/ERA525/ERA525888/pacbio_hdf5/m150811_092723_00127_c100844062550000001823187612311514_s1_p0.metadata.xml
wget ftp://ftp.sra.ebi.ac.uk/vol1/ERA525/ERA525888/pacbio_hdf5/m150814_233250_00127_c100823152550000001823177111031542_s1_p0.metadata.xml
wget ftp://ftp.sra.ebi.ac.uk/vol1/ERA535/ERA535258/pacbio_hdf5/m150911_220012_00127_c100861772550000001823190702121671_s1_p0.metadata.xml
wget ftp://ftp.sra.ebi.ac.uk/vol1/ERA525/ERA525888/pacbio_hdf5/m150813_110541_00127_c100823112550000001823177111031581_s1_p0.metadata.xml

if [[ ! -d Analysis_Results ]]
then
    mkdir Analysis_Results
fi

cd Analysis_Results
wget ftp://ftp.sra.ebi.ac.uk/vol1/ERA525/ERA525888/pacbio_hdf5/m150811_092723_00127_c100844062550000001823187612311514_s1_p0.1.bax.h5
wget ftp://ftp.sra.ebi.ac.uk/vol1/ERA525/ERA525888/pacbio_hdf5/m150811_092723_00127_c100844062550000001823187612311514_s1_p0.2.bax.h5
wget ftp://ftp.sra.ebi.ac.uk/vol1/ERA525/ERA525888/pacbio_hdf5/m150811_092723_00127_c100844062550000001823187612311514_s1_p0.3.bax.h5
wget ftp://ftp.sra.ebi.ac.uk/vol1/ERA525/ERA525888/pacbio_hdf5/m150811_092723_00127_c100844062550000001823187612311514_s1_p0.bas.h5

wget ftp://ftp.sra.ebi.ac.uk/vol1/ERA525/ERA525888/pacbio_hdf5/m150814_233250_00127_c100823152550000001823177111031542_s1_p0.1.bax.h5
wget ftp://ftp.sra.ebi.ac.uk/vol1/ERA525/ERA525888/pacbio_hdf5/m150814_233250_00127_c100823152550000001823177111031542_s1_p0.2.bax.h5
wget ftp://ftp.sra.ebi.ac.uk/vol1/ERA525/ERA525888/pacbio_hdf5/m150814_233250_00127_c100823152550000001823177111031542_s1_p0.3.bax.h5
wget ftp://ftp.sra.ebi.ac.uk/vol1/ERA525/ERA525888/pacbio_hdf5/m150814_233250_00127_c100823152550000001823177111031542_s1_p0.bas.h5

wget ftp://ftp.sra.ebi.ac.uk/vol1/ERA535/ERA535258/pacbio_hdf5/m150911_220012_00127_c100861772550000001823190702121671_s1_p0.1.bax.h5
wget ftp://ftp.sra.ebi.ac.uk/vol1/ERA535/ERA535258/pacbio_hdf5/m150911_220012_00127_c100861772550000001823190702121671_s1_p0.2.bax.h5
wget ftp://ftp.sra.ebi.ac.uk/vol1/ERA535/ERA535258/pacbio_hdf5/m150911_220012_00127_c100861772550000001823190702121671_s1_p0.3.bax.h5
wget ftp://ftp.sra.ebi.ac.uk/vol1/ERA535/ERA535258/pacbio_hdf5/m150911_220012_00127_c100861772550000001823190702121671_s1_p0.bas.h5

wget ftp://ftp.sra.ebi.ac.uk/vol1/ERA525/ERA525888/pacbio_hdf5/m150813_110541_00127_c100823112550000001823177111031581_s1_p0.1.bax.h5
wget ftp://ftp.sra.ebi.ac.uk/vol1/ERA525/ERA525888/pacbio_hdf5/m150813_110541_00127_c100823112550000001823177111031581_s1_p0.2.bax.h5
wget ftp://ftp.sra.ebi.ac.uk/vol1/ERA525/ERA525888/pacbio_hdf5/m150813_110541_00127_c100823112550000001823177111031581_s1_p0.3.bax.h5
wget ftp://ftp.sra.ebi.ac.uk/vol1/ERA525/ERA525888/pacbio_hdf5/m150813_110541_00127_c100823112550000001823177111031581_s1_p0.bas.h5

echo $(pwd)/m150811_092723_00127_c100844062550000001823187612311514_s1_p0.1.bax.h5 >> ./../$prefix.SMRTCell.1.RSII_bax.fofn
echo $(pwd)/m150811_092723_00127_c100844062550000001823187612311514_s1_p0.2.bax.h5 >> ./../$prefix.SMRTCell.1.RSII_bax.fofn
echo $(pwd)/m150811_092723_00127_c100844062550000001823187612311514_s1_p0.3.bax.h5 >> ./../$prefix.SMRTCell.1.RSII_bax.fofn
echo $(pwd)/m150814_233250_00127_c100823152550000001823177111031542_s1_p0.1.bax.h5 >> ./../$prefix.SMRTCell.2.RSII_bax.fofn
echo $(pwd)/m150814_233250_00127_c100823152550000001823177111031542_s1_p0.2.bax.h5 >> ./../$prefix.SMRTCell.2.RSII_bax.fofn
echo $(pwd)/m150814_233250_00127_c100823152550000001823177111031542_s1_p0.3.bax.h5 >> ./../$prefix.SMRTCell.2.RSII_bax.fofn
echo $(pwd)/m150911_220012_00127_c100861772550000001823190702121671_s1_p0.1.bax.h5 >> ./../$prefix.SMRTCell.3.RSII_bax.fofn
echo $(pwd)/m150911_220012_00127_c100861772550000001823190702121671_s1_p0.2.bax.h5 >> ./../$prefix.SMRTCell.3.RSII_bax.fofn
echo $(pwd)/m150911_220012_00127_c100861772550000001823190702121671_s1_p0.3.bax.h5 >> ./../$prefix.SMRTCell.3.RSII_bax.fofn
echo $(pwd)/m150813_110541_00127_c100823112550000001823177111031581_s1_p0.1.bax.h5 >> ./../$prefix.SMRTCell.4.RSII_bax.fofn
echo $(pwd)/m150813_110541_00127_c100823112550000001823177111031581_s1_p0.2.bax.h5 >> ./../$prefix.SMRTCell.4.RSII_bax.fofn
echo $(pwd)/m150813_110541_00127_c100823112550000001823177111031581_s1_p0.3.bax.h5 >> ./../$prefix.SMRTCell.4.RSII_bax.fofn

cd ..
cd ..

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
