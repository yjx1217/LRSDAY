#!/bin/bash
# last update: 2019/01/18

set -e -o pipefail

LRSDAY_HOME=$(pwd)
BUILD="build"

SRA_VERSION="2.9.2" # released on 2018.09.26
PORECHOP_VERSION="0.2.4" # 
PORECHOP_GITHUB_COMMIT_VERSION="109e437" # committed on 2018.10.19
FILTLONG_VERSION="0.2.0" #
FILTLONG_GITHUB_COMMIT_VERSION="d1bb46d" # committed on 2018.05.11
MINIMAP2_VERSION="2.13" # released on 2018.10.11
CANU_VERSION="1.8" # released on 2018.10.23
FLYE_VERSION="2.4" # released on 2018.01.15
WTDBG2_VERSION="2.3" # 
WTDBG2_GITHUB_COMMIT_VERSION="4a39621" # committed on 2019.01.15
SMARTDENOVO_VERSION="" # 
SMARTDENOVO_GITHUB_COMMIT_VERSION="5cc1356" # committed on 2018.02.19
RAGOUT_VERSION="2.1.1" # released on 2018.07.30
#QUAST_VERSION="5.0.1" # one of its dependency needs "csh" to be pre-installed
HDF_VERSION="1.10.1" # 
SONLIB_VERSION="" # 
SONLIB_GITHUB_COMMIT_VERSION="1afbd97" # committed on 2017.08.09
HAL_VERSION="" # not available, so we use the github comit hash below for version control
HAL_GITHUB_COMMIT_VERSION="a2ad656" # committed on 2017.09.09
MUMMER_VERSION="4.0.0beta2" # released on 2017.10.14
GNUPLOT_VERSION="4.6.6" # released on 2015.02.18
BEDTOOLS_VERSION="2.27.1" # released on 2017.12.14
SPADES_VERSION="3.13.0" # released on 2018.10.16
PRODIGAL_VERSION="2.6.3" # released on 2016.02.12
CAP_VERSION="" # see http://seq.cs.iastate.edu/cap3.html
BWA_VERSION="0.7.17" # released on 2017.10.23
SAMTOOLS_VERSION="1.9" # released on 2018.07.18
CIRCLATOR_VERSION="1.5.5" # released on 2018.01.31
TRIMMOMATIC_VERSION="0.36" # 
GATK_VERSION="3.6-6" #
PICARD_VERSION="2.18.15" # released on 2018.10.22
PILON_VERSION="1.22" # released on 2017.03.15
EXONERATE_VERSION="2.2.0" # 
BLAST_VERSION="2.2.31" # 
RMBLAST_VERSION="2.2.28" # 
SNAP_VERSION="" # 
SNAP_GITHUB_COMMIT_VERSION="a89d68e" # committed on 2017.05.18
RAPSEARCH_VERSION="2.24" #
TRNASCAN_VERSION="1.3.1" #
SNOSCAN_VERSION="0.9.1" #
REPEATMASKER_VERSION="open-4-0-7" #
TRF_VERSION="409" #
REANNOTATE_VERSION="17.03.2015-LongQueryName"
CLUSTALW_VERSION="2.1" #
MUSCLE_VERSION="3.8.31" #
HMMER_VERSION="3.2.1" # released on 2018.06.13
BAMTOOLS_VERSION="2.4.2" # released on 2017.11.02
AUGUSTUS_VERSION="3.2.3" # 
#AUGUSTUS_GITHUB_COMMIT_VERSION="79960c5"
EVM_VERSION="1.1.1" # released on 2015.07.03
PROTEINORTHO_VERSION="5.16b" # released on 2017.09.22
MAKER_VERSION="3.00.0-beta" #
MINICONDA2_VERSION="4.5.11" #
PB_ASSEMBLY_VERSION="0.0.2" #
BAX2BAM_VERSION="0.0.9" #
NANOPOLISH_VERSION="0.11.0" #
NANOPOLISH_GITHUB_COMMIT_VERSION="3180474" # commited on 2019.01.18 
PARALLEL_VERSION="20180722" # released on 2018.07.22
# for MFannot
EMBOSS_VERSION="6.5.7" # released on 2012.07.25
ERPIN_VERSION="5.5.4" # 
TBL2ASN_VERSION="" #
PIROBJECT_VERSION="1.19" #
PIRMODELS_GITHUB_COMMIT_VERSION="6b223ec" # committed on 2016.08.30
FLIP_GITHUB_COMMIT_VERSION="00a57cb" # committed on 2016.04.07
UMAC_GITHUB_COMMIT_VERSION="cae618e" # committed on 2016.08.30
HMMSEARCHWC_GITHUB_COMMIT_VERSION="9e3b461" # committed on 2016.11.05
RNAFINDER_GITHUB_COMMIT_VERSION="579dc58" # committed on 2016.12.07
MF2SQN_GITHUB_COMMIT_VERSION="6faf9f4" # committed on 2016.12.07
GRAB_FASTA_GITHUB_COMMIT_VERSION="accd32d" # committed on 2017.02.14
MFANNOT_DATA_GITHUB_COMMIT_VERSION="b039ac5" # committed on 2016.12.07
MFANNOT_VERSION="1.35" #
MFANNOT_GITHUB_COMMIT_VERSION="6472b97" # committed on 2018.10.31

# downloading URLs for dependencies
SRA_DOWNLOAD_URL="https://ftp-trace.ncbi.nlm.nih.gov/sra/sdk/${SRA_VERSION}/sratoolkit.${SRA_VERSION}-centos_linux64.tar.gz"
PORECHOP_DOWNLOAD_URL="https://github.com/rrwick/Porechop.git"
FILTLONG_DOWNLOAD_URL="https://github.com/rrwick/Filtlong.git"
MINIMAP2_DOWNLOAD_URL="https://github.com/lh3/minimap2/releases/download/v${MINIMAP2_VERSION}/minimap2-${MINIMAP2_VERSION}_x64-linux.tar.bz2"
CANU_DOWNLOAD_URL="https://github.com/marbl/canu/archive/v${CANU_VERSION}.tar.gz"
FLYE_DOWNLOAD_URL="https://github.com/fenderglass/Flye/archive/${FLYE_VERSION}.tar.gz"
WTDBG2_DOWNLOAD_URL="https://github.com/ruanjue/wtdbg2.git"
SMARTDENOVO_DOWNLOAD_URL="https://github.com/ruanjue/smartdenovo"
#QUAST_DOWNLOAD_URL="https://downloads.sourceforge.net/project/quast/quast-${QUAST_VERSION}.tar.gz"
RAGOUT_DOWNLOAD_URL="https://github.com/fenderglass/Ragout/archive/${RAGOUT_VERSION}.tar.gz"
HDF_VERSION_prefix=${HDF_VERSION%.*}
HDF_DOWNLOAD_URL="https://support.hdfgroup.org/ftp/HDF5/releases/hdf5-${HDF_VERSION_prefix}/hdf5-${HDF_VERSION}/src/hdf5-${HDF_VERSION}.tar.gz"
SONLIB_DOWNLOAD_URL="https://github.com/benedictpaten/sonLib.git"
HAL_DOWNLOAD_URL="https://github.com/glennhickey/hal.git"
MUMMER_DOWNLOAD_URL="https://github.com/gmarcais/mummer/releases/download/v${MUMMER_VERSION}/mummer-${MUMMER_VERSION}.tar.gz"
GNUPLOT_DOWNLOAD_URL="https://sourceforge.net/projects/gnuplot/files/gnuplot/${GNUPLOT_VERSION}/gnuplot-${GNUPLOT_VERSION}.tar.gz"
BEDTOOLS_DOWNLOAD_URL="https://github.com/arq5x/bedtools2/releases/download/v${BEDTOOLS_VERSION}/bedtools-${BEDTOOLS_VERSION}.tar.gz"
SPADES_DOWNLOAD_URL="http://cab.spbu.ru/files/release${SPADES_VERSION}/SPAdes-${SPADES_VERSION}-Linux.tar.gz"
PRODIGAL_DOWNLOAD_URL="https://github.com/hyattpd/Prodigal/archive/v${PRODIGAL_VERSION}.tar.gz"
CAP_DOWNLOAD_URL="http://seq.cs.iastate.edu/CAP3/cap3.linux.x86_64.tar"
CIRCLATOR_DOWNLOAD_URL="https://github.com/sanger-pathogens/circlator/archive/v${CIRCLATOR_VERSION}.tar.gz"
TRIMMOMATIC_DOWNLOAD_URL="http://www.usadellab.org/cms/uploads/supplementary/Trimmomatic/Trimmomatic-${TRIMMOMATIC_VERSION}.zip"
BWA_DOWNLOAD_URL="http://downloads.sourceforge.net/project/bio-bwa/bwa-${BWA_VERSION}.tar.bz2"
SAMTOOLS_DOWNLOAD_URL="https://github.com/samtools/samtools/releases/download/${SAMTOOLS_VERSION}/samtools-${SAMTOOLS_VERSION}.tar.bz2"
#GATK_DOWNLOAD_URL="https://github.com/broadgsa/gatk/archive/${GATK_VERSION}.tar.gz"
MAKER_DOWNLOAD_URL="http://topaz.genetics.utah.edu/maker_downloads/static/maker-${MAKER_VERSION}.tgz"
PICARD_DOWNLOAD_URL="https://github.com/broadinstitute/picard/releases/download/${PICARD_VERSION}/picard.jar"
PILON_DOWNLOAD_URL="https://github.com/broadinstitute/pilon/releases/download/v${PILON_VERSION}/pilon-${PILON_VERSION}.jar"
EXONERATE_DOWNLOAD_URL="http://ftp.ebi.ac.uk/pub/software/vertebrategenomics/exonerate/exonerate-${EXONERATE_VERSION}-x86_64.tar.gz"
BLAST_DOWNLOAD_URL="ftp://ftp.ncbi.nlm.nih.gov/blast/executables/blast+/${BLAST_VERSION}/ncbi-blast-${BLAST_VERSION}+-x64-linux.tar.gz"
RMBLAST_DOWNLOAD_URL="ftp://ftp.ncbi.nlm.nih.gov/blast/executables/rmblast/${RMBLAST_VERSION}/ncbi-rmblastn-${RMBLAST_VERSION}-x64-linux.tar.gz"
SNAP_DOWNLOAD_URL="https://github.com/KorfLab/SNAP"
RAPSEARCH_DOWNLOAD_URL="https://sourceforge.net/projects/rapsearch2/files/RAPSearch${RAPSEARCH_VERSION}_64bits.tar.gz"
TRNASCAN_DOWNLOAD_URL="http://eddylab.org/software/tRNAscan-SE/tRNAscan-SE.tar.gz"
SNOSCAN_DOWNLOAD_URL="http://eddylab.org/software/snoscan/snoscan.tar.gz"
REPEATMASKER_DOWNLOAD_URL="http://repeatmasker.org/RepeatMasker-${REPEATMASKER_VERSION}.tar.gz"
TRF_DOWNLOAD_URL="http://tandem.bu.edu/trf/downloads/trf${TRF_VERSION}.linux64"
REANNOTATE_DOWNLOAD_URL="https://github.com/yjx1217/REannotate_LongQueryName/archive/version_${REANNOTATE_VERSION}.tar.gz"
CLUSTALW_DOWNLOAD_URL="http://www.clustal.org/download/${CLUSTALW_VERSION}/clustalw-${CLUSTALW_VERSION}.tar.gz"
MUSCLE_DOWNLOAD_URL="http://www.drive5.com/muscle/downloads${MUSCLE_VERSION}/muscle${MUSCLE_VERSION}_i86linux64.tar.gz"
#HMMER_DOWNLOAD_URL="http://eddylab.org/software/hmmer3/${HMMER_VERSION}/hmmer-${HMMER_VERSION}-linux-intel-x86_64.tar.gz"
HMMER_DOWNLOAD_URL="http://eddylab.org/software/hmmer/hmmer-${HMMER_VERSION}.tar.gz"
BAMTOOLS_DOWNLOAD_URL="https://github.com/pezmaster31/bamtools/archive/v${BAMTOOLS_VERSION}.tar.gz"
AUGUSTUS_DOWNLOAD_URL="http://bioinf.uni-greifswald.de/augustus/binaries/old/augustus-${AUGUSTUS_VERSION}.tar.gz"
#AUGUSTUS_DOWNLOAD_URL="https://github.com/Gaius-Augustus/Augustus.git"
EVM_DOWNLOAD_URL="https://github.com/EVidenceModeler/EVidenceModeler/archive/v${EVM_VERSION}.tar.gz"
PROTEINORTHO_DOWNLOAD_URL="https://www.bioinf.uni-leipzig.de/Software/proteinortho/proteinortho_v${PROTEINORTHO_VERSION}.tar.gz"
MINICONDA2_DOWNLOAD_URL="https://repo.continuum.io/miniconda/Miniconda2-${MINICONDA2_VERSION}-Linux-x86_64.sh"
NANOPOLISH_DOWNLOAD_URL="https://github.com/jts/nanopolish.git"
PARALLEL_DOWNLOAD_URL="http://ftp.gnu.org/gnu/parallel/parallel-${PARALLEL_VERSION}.tar.bz2"

# UCSC Utilities
BLAT_DOWNLOAD_URL="http://hgdownload.soe.ucsc.edu/admin/exe/linux.x86_64/blat/blat"
FASPLIT_DOWNLOAD_URL="http://hgdownload.soe.ucsc.edu/admin/exe/linux.x86_64/faSplit"
PSLSORT_DOWNLOAD_URL="http://hgdownload.soe.ucsc.edu/admin/exe/linux.x86_64/pslSort"
PSLSCORE_DOWNLOAD_URL="http://hgdownload.soe.ucsc.edu/admin/exe/linux.x86_64/pslScore"
PSLCDNAFILTER_DOWNLOAD_URL="http://hgdownload.soe.ucsc.edu/admin/exe/linux.x86_64/pslCDnaFilter"

# for MFannot
EMBOSS_VERSION_prefix="${EMBOSS_VERSION:0:3}"
EMBOSS_DOWNLOAD_URL="ftp://emboss.open-bio.org/pub/EMBOSS/old/${EMBOSS_VERSION_prefix}.0/EMBOSS-${EMBOSS_VERSION}.tar.gz"
ERPIN_DOWNLOAD_URL="http://rna.igmors.u-psud.fr/download/Erpin/erpin${ERPIN_VERSION}.serv.tar.gz"
TBL2ASN_DOWNLOAD_URL="ftp://ftp.ncbi.nih.gov/toolbox/ncbi_tools/converters/by_program/tbl2asn/linux64.tbl2asn.gz"
PIROBJECT_DOWNLOAD_URL="https://github.com/prioux/PirObject/archive/v${PIROBJECT_VERSION}.tar.gz"
PIRMODELS_DOWNLOAD_URL="https://github.com/BFL-lab/PirModels.git"
FLIP_DOWNLOAD_URL="https://github.com/BFL-lab/Flip.git"
UMAC_DOWNLOAD_URL="https://github.com/BFL-lab/Umac.git"
HMMSEARCHWC_DOWNLOAD_URL="https://github.com/BFL-lab/HMMsearchWC.git"
RNAFINDER_DOWNLOAD_URL="https://github.com/BFL-lab/RNAfinder.git"
MF2SQN_DOWNLOAD_URL="https://github.com/BFL-lab/Mf2sqn.git"
GRAB_FASTA_DOWNLOAD_URL="https://github.com/BFL-lab/grab-fasta.git"
MFANNOT_DATA_DOWNLOAD_URL="https://github.com/BFL-lab/MFannot_data.git"
MFANNOT_DOWNLOAD_URL="https://github.com/BFL-lab/MFannot.git"


# Create the $BUILD directory for dependency installation
if [[ -d $BUILD ]]
then
    echo "The old \"$BUILD\" directory need to deleted for this installation to continue!"
    while true; do
	read -p $'Do you confirm to delete this directory? Please answer yes or no.\n' yn
	case $yn in
            [Yy]* ) echo "The answer \"yes\" is received. The old directory deleted and the new installation start now!"; rm -rf $BUILD; break;;
            [Nn]* ) echo "The answer \"no\" is received. Quit the installation!"; exit;;
            * ) echo "Please answer yes or no.";;
	esac
    done
fi

echo ""
echo "Create the new $BUILD directory"
mkdir $BUILD
cd $BUILD
build_dir=$(pwd)

echo ""
echo "Download and install all the dependencies"
echo ""

# Downloading all the dependencies
download () {
  url=$1
  download_location=$2
  echo "Downloading $url to $download_location"
  wget -nv --no-check-certificate $url -O $download_location
}

# ---------- set Perl & Python environment variables -------------
PYTHONPATH="$build_dir"
PERL5LIB="$build_dir:$PERL5LIB"
PERL5LIB="$build_dir/cpanm/perlmods/lib/perl5:$PERL5LIB"

mkdir -p  $build_dir/cpanm
cpanm_dir=$build_dir/cpanm
cd $cpanm_dir
wget -nv --no-check-certificate -O - https://cpanmin.us/ > cpanm
chmod +x cpanm
mkdir perlmods
$cpanm_dir/cpanm -l $cpanm_dir/perlmods --skip-installed Test::More@1.302086
$cpanm_dir/cpanm -l $cpanm_dir/perlmods --skip-installed Text::Soundex@3.05
$cpanm_dir/cpanm -l $cpanm_dir/perlmods --skip-installed Env@1.04
$cpanm_dir/cpanm -l $cpanm_dir/perlmods --skip-installed File::Which@1.21
$cpanm_dir/cpanm -l $cpanm_dir/perlmods --skip-installed Term::ReadKey@2.37
$cpanm_dir/cpanm -l $cpanm_dir/perlmods --skip-installed Carp@1.38
$cpanm_dir/cpanm -l $cpanm_dir/perlmods --skip-installed Perl::Unsafe::Signals@0.03
$cpanm_dir/cpanm -l $cpanm_dir/perlmods --skip-installed Carp::Clan@6.06 # dependency for Bit::Vector@7.4
$cpanm_dir/cpanm -l $cpanm_dir/perlmods --skip-installed Bit::Vector@7.4
$cpanm_dir/cpanm -l $cpanm_dir/perlmods --skip-installed Inline@0.80
$cpanm_dir/cpanm -l $cpanm_dir/perlmods --skip-installed Inline::C@0.78
$cpanm_dir/cpanm -l $cpanm_dir/perlmods --skip-installed List::MoreUtils@0.419 # dependency for forks@0.36
$cpanm_dir/cpanm -l $cpanm_dir/perlmods --skip-installed Acme::Damn@0.08 # dependency for forks@0.36
$cpanm_dir/cpanm -l $cpanm_dir/perlmods --skip-installed Sys::SigAction@0.23 # dependency for forks@0.36
$cpanm_dir/cpanm -l $cpanm_dir/perlmods --skip-installed forks@0.36
$cpanm_dir/cpanm -l $cpanm_dir/perlmods --skip-installed forks::shared@0.36
$cpanm_dir/cpanm -l $cpanm_dir/perlmods --skip-installed Want@0.29 # dependency for IO::Prompt@0.997004
$cpanm_dir/cpanm -l $cpanm_dir/perlmods --skip-installed IO::All@0.86
$cpanm_dir/cpanm -l $cpanm_dir/perlmods --skip-installed IO::Prompt@0.997004
$cpanm_dir/cpanm -l $cpanm_dir/perlmods --skip-installed DBI@1.636
$cpanm_dir/cpanm -l $cpanm_dir/perlmods --skip-installed DBD::SQLite@1.54
# $cpanm_dir/cpanm -l $cpanm_dir/perlmods --skip-installed DBD::Pg  # need $POSTGRES_HOME pre-defined
$cpanm_dir/cpanm -l $cpanm_dir/perlmods --skip-installed Proc::ProcessTable@0.53
$cpanm_dir/cpanm -l $cpanm_dir/perlmods --skip-installed threads@2.16
$cpanm_dir/cpanm -l $cpanm_dir/perlmods --skip-installed PerlIO::gzip@0.20
$cpanm_dir/cpanm -l $cpanm_dir/perlmods --skip-installed ExtUtils::CBuilder@0.280224 
$cpanm_dir/cpanm -l $cpanm_dir/perlmods --skip-installed LWP::Simple@6.26
$cpanm_dir/cpanm -l $cpanm_dir/perlmods --skip-installed LWP::UserAgent
$cpanm_dir/cpanm -l $cpanm_dir/perlmods --skip-installed Bio::Perl@1.007001


# ------------- SRA Toolkit -------------------
cd $build_dir
echo "Download SRAtoolkit-v${SRA_VERSION}"
download $SRA_DOWNLOAD_URL sratoolkit.${SRA_VERSION}-centos_linux64.tar.gz
tar -zxf sratoolkit.${SRA_VERSION}-centos_linux64.tar.gz
sra_dir="$build_dir/sratoolkit.${SRA_VERSION}-centos_linux64/bin"
rm sratoolkit.${SRA_VERSION}-centos_linux64.tar.gz

# --------------- Porechop ------------------
cd $build_dir
echo "Download Porechop-v${PORECHOP_VERSION}"
git clone $PORECHOP_DOWNLOAD_URL
cd Porechop
git checkout -f -q $PORECHOP_GITHUB_COMMIT_VERSION
virtualenv -p $(which python3) py3_virtualenv_porechop
source py3_virtualenv_porechop/bin/activate
py3_virtualenv_porechop/bin/python3 ./setup.py install
deactivate
porechop_dir="$build_dir/Porechop/py3_virtualenv_porechop/bin"

# --------------- Porechop ------------------
cd $build_dir
echo "Download Filtlong-v${FILTLONG_VERSION}"
git clone $FILTLONG_DOWNLOAD_URL
cd Filtlong
git checkout -f -q $FILTLONG_GITHUB_COMMIT_VERSION
make -j
filtlong_dir="$build_dir/Filtlong/bin"

# --------------- minimap2 ------------------
cd $build_dir
echo "Download minimap2-v${MINIMAP2_VERSION}"
download $MINIMAP2_DOWNLOAD_URL "minimap2-${MINIMAP2_VERSION}.tar.bz2"
tar -xjf minimap2-${MINIMAP2_VERSION}.tar.bz2
minimap2_dir="$build_dir/minimap2-${MINIMAP2_VERSION}_x64-linux"
rm minimap2-${MINIMAP2_VERSION}.tar.bz2

# ------------- Canu -------------------
cd $build_dir
echo "Download Canu-v${CANU_VERSION}"
download $CANU_DOWNLOAD_URL "canu-${CANU_VERSION}.tar.gz"
tar -xzf canu-${CANU_VERSION}.tar.gz
canu_dir="$build_dir/canu-${CANU_VERSION}"
cd $canu_dir
cd src
make -j 8
canu_dir="$build_dir/canu-${CANU_VERSION}/Linux-amd64/bin"
cd $canu_dir
ln -s $minimap2_dir/minimap2 .
cd $build_dir
rm canu-${CANU_VERSION}.tar.gz

# ------------- Flye -------------------
cd $build_dir
echo "Download Flye-v${FLYE_VERSION}"
download $FLYE_DOWNLOAD_URL "Flye-${FLYE_VERSION}.tar.gz"
tar -xzf Flye-${FLYE_VERSION}.tar.gz
cd Flye-${FLYE_VERSION}
python2 setup.py build
cd ..
flye_dir="$build_dir/Flye-${FLYE_VERSION}/bin"
rm Flye-${FLYE_VERSION}.tar.gz

# --------------- wtdbg2 ------------------
cd $build_dir
echo "Download wtdbg2-v${WTDBG2_VERSION}"
git clone $WTDBG2_DOWNLOAD_URL
cd wtdbg2
git checkout -f -q $WTDBG2_GITHUB_COMMIT_VERSION
C_INCLUDE_PATH="" 
make
wtdbg2_dir="$build_dir/wtdbg2"

# --------------- smartdenovo ------------------
cd $build_dir
echo "Download smartdenovo-v${SMARTDENOVO_VERSION}"
git clone $SMARTDENOVO_DOWNLOAD_URL
cd smartdenovo
git checkout -f -q $SMARTDENOVO_GITHUB_COMMIT_VERSION
cp wtlay.h wtlay.h.bk
cat wtlay.h.bk |sed s/inline//g > wtlay.h
C_INCLUDE_PATH="" 
make
cp $LRSDAY_HOME/misc/smartdenovo_customized.pl .
cd ..
smartdenovo_dir="$build_dir/smartdenovo"

# --------------- Ragout ------------------
cd $build_dir
echo "Download Ragout-v${RAGOUT_VERSION}"
download $RAGOUT_DOWNLOAD_URL Ragout-${RAGOUT_VERSION}.tar.gz
tar -zxf Ragout-${RAGOUT_VERSION}.tar.gz
cd Ragout-${RAGOUT_VERSION}
python2 setup.py build
python2 ./scripts/install-sibelia.py
ragout_dir="$build_dir/Ragout-${RAGOUT_VERSION}/bin"
cd ..
rm Ragout-${RAGOUT_VERSION}.tar.gz

# --------------- HDF ------------------
cd $build_dir
echo "Download HDF-v${HDF_VERSION}"
download $HDF_DOWNLOAD_URL "hdf5-${HDF_VERSION}.tar.gz"
tar -zxf hdf5-${HDF_VERSION}.tar.gz
mkdir hdf5
cd hdf5-${HDF_VERSION}
./configure --enable-cxx  --prefix $build_dir/hdf5
make
make install
hdf_dir="$build_dir/hdf5/bin"
PATH="$hdf_dir:${PATH}"
h5prefix="-prefix=$build_dir/hdf5"
cd $build_dir
rm hdf5-${HDF_VERSION}.tar.gz

# --------------- sonLib ------------------
cd $build_dir
git clone $SONLIB_DOWNLOAD_URL
cd sonLib
git checkout -f -q $SONLIB_GITHUB_COMMIT_VERSION
make

# ---------------- HAL -------------------
cd $build_dir
echo "Download HAL-v${HAL_VERSION}"
git clone $HAL_DOWNLOAD_URL
cd hal
git checkout -f -q $HAL_GITHUB_COMMIT_VERSION
make
hal_dir="$build_dir/hal/bin"

# --------------- samtools -----------------
cd $build_dir
echo "Download samtools-v${SAMTOOLS_VERSION}"
download $SAMTOOLS_DOWNLOAD_URL "samtools-${SAMTOOLS_VERSION}.tar.bz2"
tar -xjf samtools-${SAMTOOLS_VERSION}.tar.bz2
samtools_dir="$build_dir/samtools-${SAMTOOLS_VERSION}"
cd $samtools_dir
C_INCLUDE_PATH=""
./configure --without-curses;
make
cd $build_dir
rm samtools-${SAMTOOLS_VERSION}.tar.bz2

# --------------- gnuplot ------------------
cd $build_dir
echo "Download gnuplot-v${GNUPLOT_VERSION}"
download $GNUPLOT_DOWNLOAD_URL "gnuplot-${GNUPLOT_VERSION}.tar.gz"
tar -zxf gnuplot-${GNUPLOT_VERSION}.tar.gz
cd "$build_dir/gnuplot-${GNUPLOT_VERSION}"
./configure --prefix="$build_dir/gnuplot-${GNUPLOT_VERSION}" --disable-wxwidgets
make
make install
gnuplot_dir="$build_dir/gnuplot-${GNUPLOT_VERSION}/bin"
cd $build_dir
rm gnuplot-${GNUPLOT_VERSION}.tar.gz
PATH="$gnuplot_dir:${PATH}"

# # ------------- QUAST --------------------
# cd $build_dir
# echo "Download QUAST-v${QUAST_VERSION}"
# download $QUAST_DOWNLOAD_URL "QUAST-${QUAST_VERSION}.tar.gz"
# tar -xzf QUAST-${QUAST_VERSION}.tar.gz
# quast_dir="$build_dir/quast-${QUAST_VERSION}"
# cd $quast_dir
# virtualenv -p $(which python3) py3_virtualenv_quast
# source py3_virtualenv_quast/bin/activate
# py3_virtualenv_quast/bin/pip install joblib
# py3_virtualenv_quast/bin/pip install simplejson
# py3_virtualenv_quast/bin/python3 -mpip install -U matplotlib
# py3_virtualenv_quast/bin/python3 ./setup.py install
# deactivate
# cd ..
# rm QUAST-${QUAST_VERSION}.tar.gz

# --------------- mummer ------------------
cd $build_dir
echo "Download mummer-v${MUMMER_VERSION}"
download $MUMMER_DOWNLOAD_URL "mummer-${MUMMER_VERSION}.tar.gz"
tar -zxf mummer-${MUMMER_VERSION}.tar.gz
mummer_dir="$build_dir/mummer-${MUMMER_VERSION}"
echo "$mummer_dir"
cd $mummer_dir
./configure
make
PATH="$mummer_dir:${PATH}"
cd $build_dir
rm mummer-${MUMMER_VERSION}.tar.gz

# --------------- bedtools ------------------
cd $build_dir
echo "Download bedtools-v${BEDTOOLS_VERSION}"
download $BEDTOOLS_DOWNLOAD_URL "bedtools-${BEDTOOLS_VERSION}.tar.gz"
tar -zxf bedtools-${BEDTOOLS_VERSION}.tar.gz
cd "$build_dir/bedtools2"
make
bedtools_dir="$build_dir/bedtools2/bin"
cd $build_dir
rm bedtools-${BEDTOOLS_VERSION}.tar.gz

# --------------- SPAdes ------------------
cd $build_dir
echo "Download SPAdes-v${SPADES_VERSION}"
download $SPADES_DOWNLOAD_URL "SPAdes-${SPADES_VERSION}-Linux.tar.gz"
tar -zxf SPAdes-${SPADES_VERSION}-Linux.tar.gz
spades_dir="$build_dir/SPAdes-${SPADES_VERSION}-Linux/bin"
rm SPAdes-${SPADES_VERSION}-Linux.tar.gz

# --------------- Prodigal ------------------
cd $build_dir
echo "Download Prodigal-v${PRODIGAL_VERSION}"
download $PRODIGAL_DOWNLOAD_URL "v${PRODIGAL_VERSION}.tar.gz"
tar -zxf v${PRODIGAL_VERSION}.tar.gz
prodigal_dir="$build_dir/Prodigal-${PRODIGAL_VERSION}"
cd $prodigal_dir
make
cd $build_dir
rm v${PRODIGAL_VERSION}.tar.gz

# --------------- CAP3 ------------------
cd $build_dir
echo "Download CAP3-v${CAP3_VERSION}"
download $CAP_DOWNLOAD_URL "cap3.linux.x86_64.tar"
tar -xf cap3.linux.x86_64.tar
cap_dir="$build_dir/CAP3"
rm cap3.linux.x86_64.tar

# ------------- BWA -------------------
cd $build_dir
echo "Download BWA-v${BWA_VERSION}"
download $BWA_DOWNLOAD_URL "bwa-${BWA_VERSION}.tar.bz2"
tar -xjf bwa-${BWA_VERSION}.tar.bz2
bwa_dir="$build_dir/bwa-${BWA_VERSION}"
cd $bwa_dir
make
cd $build_dir
rm bwa-${BWA_VERSION}.tar.bz2

# --------------- Circlator ------------------
cd $build_dir
echo "Creating local virtual python3 environment and install Circlator-v${CIRCLATOR_VERSION}"
virtualenv -p $(which python3) py3_virtualenv_circlator
source py3_virtualenv_circlator/bin/activate
py3_virtualenv_circlator/bin/pip3 install "circlator==${CIRCLATOR_VERSION}"
circlator_dir="$build_dir/py3_virtualenv_circlator/bin"
deactivate

# --------------- Trimmomatic -----------------
cd $build_dir
echo "Download Trimmomatic-v${TRIMMOMATIC_VERSION}"
download $TRIMMOMATIC_DOWNLOAD_URL "Trimmomatic-${TRIMMOMATIC_VERSION}.zip"
unzip Trimmomatic-${TRIMMOMATIC_VERSION}.zip
trimmomatic_dir="$build_dir/Trimmomatic-${TRIMMOMATIC_VERSION}"
cd $trimmomatic_dir
chmod 755 trimmomatic-${TRIMMOMATIC_VERSION}.jar
ln -s trimmomatic-${TRIMMOMATIC_VERSION}.jar trimmomatic.jar 
cd $build_dir
rm Trimmomatic-${TRIMMOMATIC_VERSION}.zip

# --------------- Picard -----------------
cd $build_dir
echo "Download Picard-v${PICARD_VERSION}"
download $PICARD_DOWNLOAD_URL "picard.jar"
mkdir Picard-v${PICARD_VERSION}
picard_dir="$build_dir/Picard-v${PICARD_VERSION}"
mv picard.jar $picard_dir
cd $picard_dir
chmod 755 picard.jar

# --------------- Pilon -----------------
cd $build_dir
echo "Download Pilon-v${PILON_VERSION}"
download $PILON_DOWNLOAD_URL "pilon-${PILON_VERSION}.jar"
mkdir Pilon-v${PILON_VERSION}
pilon_dir="$build_dir/Pilon-v${PILON_VERSION}"
mv pilon-${PILON_VERSION}.jar $pilon_dir/pilon.jar
cd $pilon_dir
chmod 755 pilon.jar

# --------------- exonerate ------------------
cd $build_dir
echo "Download exonerate-v${EXONERATE_VERSION}"
download $EXONERATE_DOWNLOAD_URL "exonerate-${EXONERATE_VERSION}-x86_64.tar.gz"
tar -zxf exonerate-${EXONERATE_VERSION}-x86_64.tar.gz
exonerate_dir="$build_dir/exonerate-${EXONERATE_VERSION}-x86_64/bin"
rm exonerate-${EXONERATE_VERSION}-x86_64.tar.gz

# --------------- ncbi-blast+ ------------------
cd $build_dir
echo "Download ncbi-blast-v${BLAST_VERSION}"
download $BLAST_DOWNLOAD_URL "ncbi-blast-${BLAST_VERSION}+-x64-linux.tar.gz"
tar -zxf ncbi-blast-${BLAST_VERSION}+-x64-linux.tar.gz
rm ncbi-blast-${BLAST_VERSION}+-x64-linux.tar.gz
cd "ncbi-blast-${BLAST_VERSION}+"
mkdir matrices
cd matrices
wget ftp://ftp.ncbi.nlm.nih.gov/blast/matrices/*
blast_dir="$build_dir/ncbi-blast-${BLAST_VERSION}+/bin"
blast_matrices_dir="$build_dir/ncbi-blast-${BLAST_VERSION}+/matrices"

# --------------- ncbi-rmblast ------------------
cd $build_dir
echo "Download ncbi-rmblastn-v${BLAST_VERSION}"
download $RMBLAST_DOWNLOAD_URL "ncbi-rmblastn-${RMBLAST_VERSION}-x64-linux.tar.gz"
tar -zxf ncbi-rmblastn-${RMBLAST_VERSION}-x64-linux.tar.gz
rmblast_dir="$build_dir/ncbi-rmblastn-${RMBLAST_VERSION}/bin"
# copy rmblastn binary file to ncbi-blast+ directory for easy RepeatMasker configuration
cp $rmblast_dir/rmblastn $blast_dir
rm ncbi-rmblastn-${RMBLAST_VERSION}-x64-linux.tar.gz

# --------------- snap ------------------
cd $build_dir
echo "Download snap-v${SNAP_VERSION}"
git clone $SNAP_DOWNLOAD_URL
snap_dir="$build_dir/SNAP"
cd $snap_dir
git checkout -f -q $SNAP_GITHUB_COMMIT_VERSION
ZOE="$snap_dir/Zoe"
cp $LRSDAY_HOME/misc/snap.c .  # temporary fix for snap with gcc-8
make

# --------------- RAPSearch2 ------------------
cd $build_dir
echo "Download RAPsearch-v${RAPSEARCH_VERSION}"
download $RAPSEARCH_DOWNLOAD_URL "RAPSearch${RAPSEARCH_VERSION}_64bits.tar.gz"
tar -zxf RAPSearch${RAPSEARCH_VERSION}_64bits.tar.gz
rapsearch_dir="$build_dir/RAPSearch${RAPSEARCH_VERSION}_64bits/bin"
rm RAPSearch${RAPSEARCH_VERSION}_64bits.tar.gz

# --------------- tRNAscan-SE ------------------
cd $build_dir
echo "Download tRNAscan-SE-v${TRNASCAN_VERSION}"
download $TRNASCAN_DOWNLOAD_URL "tRNAscan-SE-${TRNASCAN_VERSION}.tar.gz"
tar -zxf tRNAscan-SE-${TRNASCAN_VERSION}.tar.gz
trnascan_dir="$build_dir/tRNAscan-SE-${TRNASCAN_VERSION}"
cd $trnascan_dir
mkdir bin
mkdir -p lib/tRNAscan-SE
mkdir -p man
cp $LRSDAY_HOME/misc/tRNAscan-SE.Makefile Makefile
PERL5LIB=$trnascan_dir/bin:$PERL5LIB
make BINDIR="$trnascan_dir/bin" LIBDIR="$trnascan_dir/lib/tRNAscan-SE" MANDIR="$trnascan_dir"
make install BINDIR="$trnascan_dir/bin" LIBDIR="$trnascan_dir/lib/tRNAscan-SE" MANDIR="$trnascan_dir"
trnascan_dir="$build_dir/tRNAscan-SE-${TRNASCAN_VERSION}/bin"
cd $build_dir
rm tRNAscan-SE-${TRNASCAN_VERSION}.tar.gz

# --------------- snoscan ------------------
cd $build_dir
echo "Download snoscan-v${SNOSCAN_VERSION}"
download $SNOSCAN_DOWNLOAD_URL "snoscan-${SNOSCAN_VERSION}.tar.gz"
tar -zxf snoscan-${SNOSCAN_VERSION}.tar.gz
snoscan_dir="$build_dir/snoscan-${SNOSCAN_VERSION}"
cd $snoscan_dir
cd squid-1.5.11
rm *.o
make
cd ..
cp $LRSDAY_HOME/misc/snoscan.Makefile Makefile
rm *.o
make
cd $build_dir
rm snoscan-${SNOSCAN_VERSION}.tar.gz

# --------------- RepeatMasker ------------------
cd $build_dir
echo "Download Repeatmasker-v${REPEATMASKER_VERSION}"
download $REPEATMASKER_DOWNLOAD_URL "RepeatMasker-${REPEATMASKER_VERSION}.tar.gz"
tar -zxf RepeatMasker-${REPEATMASKER_VERSION}.tar.gz
repeatmasker_dir="$build_dir/RepeatMasker"
cd $repeatmasker_dir
echo "Download and setup RepBase library"
REPBASE_VERSION="20170127"
wget https://github.com/yjx1217/RMRB/raw/master/RepBaseRepeatMaskerEdition-${REPBASE_VERSION}.tar.gz
tar xzf RepBaseRepeatMaskerEdition-${REPBASE_VERSION}.tar.gz
rm RepBaseRepeatMaskerEdition-${REPBASE_VERSION}.tar.gz
cd .. 
rm RepeatMasker-${REPEATMASKER_VERSION}.tar.gz

# --------------- TRF ------------------
cd $repeatmasker_dir
echo "Download TRF-v${TRF_VERSION}"
download $TRF_DOWNLOAD_URL "trf${TRF_VERSION}.linux64"
mv trf${TRF_VERSION}.linux64 trf
chmod 755 trf
trf_dir=$repeatmasker_dir

# --------------- REannotate ------------------
cd $build_dir
echo "Download REannotate-v${REANNOTATE_VERSION}"
download $REANNOTATE_DOWNLOAD_URL "version_${REANNOTATE_VERSION}.tar.gz"
tar -zxf version_${REANNOTATE_VERSION}.tar.gz
reannotate_dir="$build_dir/REannotate_LongQueryName-version_${REANNOTATE_VERSION}"
cd $reannotate_dir
chmod 755 REannotate_longname
ln -s REannotate_longname REannotate
cd $build_dir
rm version_${REANNOTATE_VERSION}.tar.gz

# --------------- ClustalW ------------------
cd $build_dir
echo "Download ClustalW-v${CLUSTALW_VERSION}"
download $CLUSTALW_DOWNLOAD_URL "clustalw-${CLUSTALW_VERSION}.tar.gz"
tar -zxf clustalw-${CLUSTALW_VERSION}.tar.gz
cd clustalw-${CLUSTALW_VERSION}
./configure --prefix="$build_dir/clustalw-${CLUSTALW_VERSION}"
make
make install
clustalw_dir="$build_dir/clustalw-${CLUSTALW_VERSION}/bin"
cd $build_dir
rm clustalw-${CLUSTALW_VERSION}.tar.gz

# --------------- MUSCLE ------------------
cd $build_dir
echo "Download MUSCLE-v${MUSCLE_VERSION}"
download $MUSCLE_DOWNLOAD_URL "muscle-${MUSCLE_VERSION}_i86linux64.tar.gz"
tar -zxf muscle-${MUSCLE_VERSION}_i86linux64.tar.gz
mkdir muscle-${MUSCLE_VERSION}
mv muscle${MUSCLE_VERSION}_i86linux64 ./muscle-${MUSCLE_VERSION}/muscle
muscle_dir="$build_dir/muscle-${MUSCLE_VERSION}"
rm muscle-${MUSCLE_VERSION}_i86linux64.tar.gz

# --------------- HMMER ------------------
cd $build_dir
echo "Download hmmer-v${HMMER_VERSION}"
download $HMMER_DOWNLOAD_URL "hmmer-${HMMER_VERSION}-linux-intel-x86_64.tar.gz"
tar -zxf hmmer-${HMMER_VERSION}-linux-intel-x86_64.tar.gz
hmmer_dir="$build_dir/hmmer-${HMMER_VERSION}"
cd $hmmer_dir
./configure --prefix=$hmmer_dir
make
make install
cd easel
make install
cd ..
cd ..
hmmer_dir="$build_dir/hmmer-${HMMER_VERSION}/bin"
rm hmmer-${HMMER_VERSION}-linux-intel-x86_64.tar.gz

# --------------- bamtools ------------------
cd $build_dir
echo "Download bamtools-v${BAMTOOLS_VERSION}"
download $BAMTOOLS_DOWNLOAD_URL "v${BAMTOOLS_VERSION}.tar.gz"
tar -zxf v${BAMTOOLS_VERSION}.tar.gz
cd $build_dir/bamtools-${BAMTOOLS_VERSION}
mkdir build
cd build
cmake -DCMAKE_INSTALL_PREFIX="$build_dir/bamtools-${BAMTOOLS_VERSION}" ..
make
make install
bamtools_dir="$build_dir/bamtools-${BAMTOOLS_VERSION}/bin"
cd $build_dir
rm v${BAMTOOLS_VERSION}.tar.gz

# --------------- Augustus ------------------
cd $build_dir
echo "Download Augustus-v${AUGUSTUS_VERSION}"
download $AUGUSTUS_DOWNLOAD_URL "augustus-${AUGUSTUS_VERSION}.tar.gz"
tar -zxf augustus-${AUGUSTUS_VERSION}.tar.gz
cd $build_dir/augustus-${AUGUSTUS_VERSION}/auxprogs/bam2hints/
cp $LRSDAY_HOME/misc/bam2hints.Makefile Makefile
cd $build_dir/augustus-${AUGUSTUS_VERSION}/auxprogs/filterBam/src/
cp $LRSDAY_HOME/misc/filterBam.Makefile Makefile
cd $build_dir/augustus-${AUGUSTUS_VERSION}
make BAMTOOLS="$build_dir/bamtools-${BAMTOOLS_VERSION}"
augustus_dir="$build_dir/augustus-${AUGUSTUS_VERSION}/bin"
export AUGUSTUS_CONFIG_PATH="$build_dir/augustus-${AUGUSTUS_VERSION}/config"
cd $build_dir
rm augustus-${AUGUSTUS_VERSION}.tar.gz

# cd $build_dir
# echo "Download Augustus-v${AUGUSTUS_VERSION}"
# git clone $AUGUSTUS_DOWNLOAD_URL
# cd Augustus
# git submodule update --init
# git checkout -f -q $AUGUSTUS_GITHUB_COMMIT_VERSION
# cd $build_dir/augustus-${AUGUSTUS_VERSION}/auxprogs/bam2hints/
# cp $LRSDAY_HOME/misc/bam2hints.Makefile Makefile
# cd $build_dir/augustus-${AUGUSTUS_VERSION}/auxprogs/filterBam/src/
# cp $LRSDAY_HOME/misc/filterBam.Makefile Makefile
# cd $build_dir/augustus-${AUGUSTUS_VERSION}
# make BAMTOOLS="$build_dir/bamtools-${BAMTOOLS_VERSION}"
# augustus_dir="$build_dir/augustus-${AUGUSTUS_VERSION}/bin"
# export AUGUSTUS_CONFIG_PATH="$build_dir/augustus-${AUGUSTUS_VERSION}/config"
# cd $build_dir
# rm augustus-${AUGUSTUS_VERSION}.tar.gz


# --------------- EVidenceModeler ------------------
cd $build_dir
echo "Download EvidenceModeler-v${EVM_VERSION}"
download $EVM_DOWNLOAD_URL "v${EVM_VERSION}.tar.gz"
tar -zxf v${EVM_VERSION}.tar.gz
evm_dir="$build_dir/EVidenceModeler-${EVM_VERSION}"
rm v${EVM_VERSION}.tar.gz

# --------------- Proteinortho ------------------
cd $build_dir
echo "Download Proteinortho-v${PROTEINORTHO_VERSION}"
download $PROTEINORTHO_DOWNLOAD_URL "proteinortho_v${PROTEINORTHO_VERSION}.tar.gz"
tar -zxf proteinortho_v${PROTEINORTHO_VERSION}.tar.gz
proteinortho_dir="$build_dir/proteinortho_v${PROTEINORTHO_VERSION}"
cp $LRSDAY_HOME/misc/proteinortho5_better_robustness.pl $proteinortho_dir/proteinortho5.pl
rm proteinortho_v${PROTEINORTHO_VERSION}.tar.gz

# --------------- GATK ------------------
cd $build_dir
echo "Create GATK3 folder for users' manual installation"
mkdir GATK3
cd GATK3
wget https://ftp-trace.ncbi.nlm.nih.gov/sra/sdk/${SRA_VERSION}/GenomeAnalysisTK.jar
chmod 755 GenomeAnalysisTK.jar
gatk_dir="$build_dir/GATK3"

# --------------- MAKER -----------------
cd $build_dir
echo "Download MAKER"
download $MAKER_DOWNLOAD_URL "maker-v${MAKER_VERSION}.tgz"
tar -zxf maker-v${MAKER_VERSION}.tgz
rm maker-v${MAKER_VERSION}.tgz
cd $build_dir/maker/src/
cp $LRSDAY_HOME/misc/maker_Build.PL .
echo "no"|perl maker_Build.PL
./Build install
maker_dir="$build_dir/maker/bin"

# --------------- UCSC Utilities -----------------
cd $build_dir
mkdir UCSC_Utilities
ucsc_dir="$build_dir/UCSC_Utilities"
cd $ucsc_dir
download $BLAT_DOWNLOAD_URL "blat"
download $FASPLIT_DOWNLOAD_URL "faSplit"
download $PSLSORT_DOWNLOAD_URL "pslSort"
download $PSLSCORE_DOWNLOAD_URL "pslScore"
download $PSLCDNAFILTER_DOWNLOAD_URL "pslCDnaFilter"
chmod 755 $ucsc_dir/*

# ------------- Conda-PacBio --------------------
cd $build_dir
download $MINICONDA2_DOWNLOAD_URL "Miniconda2-${MINICONDA2_VERSION}-Linux-x86_64.sh"
bash Miniconda2-${MINICONDA2_VERSION}-Linux-x86_64.sh -b -p $build_dir/miniconda2
#export PATH="$build_dir/miniconda2/bin:$PATH"
miniconda2_dir="$build_dir/miniconda2/bin"
$miniconda2_dir/conda config --add channels defaults
$miniconda2_dir/conda config --add channels bioconda
$miniconda2_dir/conda config --add channels conda-forge
$miniconda2_dir/conda create -y -p $build_dir/conda_pacbio_env
source $miniconda2_dir/activate $build_dir/conda_pacbio_env
$miniconda2_dir/conda install -y -c bioconda pb-assembly=${PB_ASSEMBLY_VERSION}
$miniconda2_dir/conda install -y -c bioconda bax2bam=${BAX2BAM_VERSION}
source $miniconda2_dir/deactivate
rm Miniconda2-${MINICONDA2_VERSION}-Linux-x86_64.sh 
conda_pacbio_dir=$build_dir/conda_pacbio_env/bin

# ------------- NANOPOLISH --------------------
cd $build_dir
echo "Download nanopolish-v${NANOPOLISH_VERSION}"
git clone --recursive $NANOPOLISH_DOWNLOAD_URL
nanopolish_dir="$build_dir/nanopolish"
cd $nanopolish_dir
git checkout -f -q $NANOPOLISH_GITHUB_COMMIT_VERSION
make
virtualenv -p $(which python3) py3_virtualenv_nanopolish
source py3_virtualenv_nanopolish/bin/activate
py3_virtualenv_nanopolish/bin/pip install biopython
py3_virtualenv_nanopolish/bin/pip install pysam
deactivate
rm *.tar.gz
rm *.tar.bz2

# --------------- parallel ------------------
cd $build_dir
echo "Download parallel"
download $PARALLEL_DOWNLOAD_URL "parallel_v${PARALLEL_VERSION}.tar.bz2"
tar -jxf parallel_v${PARALLEL_VERSION}.tar.bz2
cd parallel-${PARALLEL_VERSION}
./configure --prefix="$build_dir/parallel-${PARALLEL_VERSION}"
make
make install
parallel_dir="$build_dir/parallel-${PARALLEL_VERSION}/bin"
cd ..
rm parallel_v${PARALLEL_VERSION}.tar.bz2

# --------------- EMBOSS ------------------
cd $build_dir
echo "Download EMBOSS"
download $EMBOSS_DOWNLOAD_URL "emboss_v${EMBOSS_VERSION}.tar.gz"
tar -zxf emboss_v${EMBOSS_VERSION}.tar.gz
cd EMBOSS-${EMBOSS_VERSION}
./configure
make
emboss_dir="$build_dir/EMBOSS-${EMBOSS_VERSION}/emboss"
cd ..
rm emboss_v${EMBOSS_VERSION}.tar.gz

# --------------- ERPIN ------------------
cd $build_dir
echo "Download ERPIN"
download $ERPIN_DOWNLOAD_URL "erpin_v${ERPIN_VERSION}.tar.gz"
tar -zxf erpin_v${ERPIN_VERSION}.tar.gz
cd erpin${ERPIN_VERSION}.serv
make
erpin_dir="$build_dir/erpin${ERPIN_VERSION}.serv/bin"
cd $build_dir
rm erpin_v${ERPIN_VERSION}.tar.gz

# --------------- tbl2asn ------------------
cd $build_dir
echo "Download tbl2asn"
download $TBL2ASN_DOWNLOAD_URL "tbl2asn.gz"
mkdir tbl2asn_dir
gunzip tbl2asn.gz
chmod 755 tbl2asn
mv tbl2asn ./tbl2asn_dir/
tbl2asn_dir="$build_dir/tbl2asn_dir"
cd $build_dir

# --------------- PirObject ----------------
cd $build_dir
echo "Download PirObject"
download $PIROBJECT_DOWNLOAD_URL "pirobject_v${PIROBJECT_VERSION}.tar.gz"
tar -zxf pirobject_v${PIROBJECT_VERSION}.tar.gz
cd PirObject-${PIROBJECT_VERSION}
pirobject_dir="$build_dir/PirObject-${PIROBJECT_VERSION}"
ln -s ./../lib/PirObject.pm .
cd $build_dir
rm pirobject_v${PIROBJECT_VERSION}.tar.gz

# --------------- PirModels ------------------
cd $build_dir
echo "Download PirModels"
git clone $PIRMODELS_DOWNLOAD_URL
cd PirModels
git checkout -f -q $PIRMODELS_GITHUB_COMMIT_VERSION
cd ..
cp -r PirModels $pirobject_dir
pirmodels_dir="$perlobject_dir/PirModels"

# --------------- Flip ------------------
cd $build_dir
echo "Download Flip"
git clone $FLIP_DOWNLOAD_URL
cd Flip
git checkout -f -q $FLIP_GITHUB_COMMIT_VERSION
cd src
make
cp flip ./../
flip_dir="$build_dir/Flip"

# --------------- Umac ------------------
cd $build_dir
echo "Download Umac"
git clone $UMAC_DOWNLOAD_URL
cd Umac
git checkout -f -q $UMAC_GITHUB_COMMIT_VERSION
umac_dir="$build_dir/Umac"

# --------------- HMMsearchWC ------------------
cd $build_dir
echo "Download HMMsearchWC"
git clone $HMMSEARCHWC_DOWNLOAD_URL
cd HMMsearchWC
git checkout -f -q $HMMSEARCHWC_GITHUB_COMMIT_VERSION
hmmsearchwc_dir="$build_dir/HMMsearchWC"

# --------------- RNAfinder ------------------
cd $build_dir
echo "Download RNAfinder"
git clone $RNAFINDER_DOWNLOAD_URL
cd RNAfinder
git checkout -f -q $RNAFINDER_GITHUB_COMMIT_VERSION
rnafinder_dir="$build_dir/RNAfinder"

# --------------- Mf2sqn ------------------
cd $build_dir
echo "Download Mf2sqn"
git clone $MF2SQN_DOWNLOAD_URL
cd Mf2sqn
git checkout -f -q $MF2SQN_GITHUB_COMMIT_VERSION
mf2sqn_dir="$build_dir/Mf2sqn"
cp qualifs.pl $build_dir/cpanm/perlmods/lib/perl5

# --------------- grab-fasta ------------------
cd $build_dir
echo "Download grab-fasta"
git clone $GRAB_FASTA_DOWNLOAD_URL
cd grab-fasta
git checkout -f -q $GRAB_FASTA_GITHUB_COMMIT_VERSION
grab_fasta_dir="$build_dir/grab-fasta"

# --------------- MFannot_data ------------------
cd $build_dir
echo "Download MFannot_data"
git clone $MFANNOT_DATA_DOWNLOAD_URL
cd MFannot_data
git checkout -f -q $MFANNOT_DATA_GITHUB_COMMIT_VERSION
mfannot_data_dir="$build_dir/MFannot_data"

# --------------- MFannot ------------------
cd $build_dir
echo "Download MFannot"
git clone $MFANNOT_DOWNLOAD_URL
cd MFannot
git checkout -f -q $MFANNOT_GITHUB_COMMIT_VERSION
mfannot_dir="$build_dir/MFannot"


# Configure executable paths
cd $LRSDAY_HOME
echo "Configuring executable paths ..."
echo "export LRSDAY_HOME=${LRSDAY_HOME}" > env.sh
echo "export PYTHONPATH=${PYTHONPATH}" >> env.sh
echo "export PERL5LIB=${PERL5LIB}" >> env.sh 
echo "export cpanm_dir=${cpanm_dir}" >> env.sh
echo "export sra_dir=${sra_dir}" >> env.sh
echo "export porechop_dir=${porechop_dir}" >> env.sh
echo "export filtlong_dir=${filtlong_dir}" >> env.sh
echo "export minimap2_dir=${minimap2_dir}" >> env.sh
echo "export canu_dir=${canu_dir}" >> env.sh
echo "export flye_dir=${flye_dir}" >> env.sh
echo "export wtdbg2_dir=${wtdbg2_dir}" >> env.sh
echo "export smartdenovo_dir=${smartdenovo_dir}" >> env.sh
#echo "export quast_dir=${quast_dir}" >> env.sh
echo "export ragout_dir=${ragout_dir}" >> env.sh
echo "export hdf_dir=${hdf_dir}" >> env.sh
echo "export h5prefix=${h5prefix}" >> env.sh
echo "export hal_dir=${hal_dir}" >> env.sh
echo "export mummer_dir=${mummer_dir}" >> env.sh
echo "export gnuplot_dir=${gnuplot_dir}" >> env.sh
echo "export bedtools_dir=${bedtools_dir}" >> env.sh
echo "export spades_dir=${spades_dir}" >> env.sh
echo "export prodigal_dir=${prodigal_dir}" >> env.sh
echo "export cap_dir=${cap_dir}" >> env.sh
echo "export circlator_dir=${circlator_dir}" >> env.sh
echo "export trimmomatic_dir=${trimmomatic_dir}" >> env.sh
echo "export bwa_dir=${bwa_dir}" >> env.sh
echo "export samtools_dir=${samtools_dir}" >> env.sh
echo "export picard_dir=${picard_dir}" >> env.sh
echo "export pilon_dir=${pilon_dir}" >> env.sh
echo "export exonerate_dir=${exonerate_dir}" >> env.sh
echo "export blast_dir=${blast_dir}" >> env.sh
echo "export blast_matrices_dir=${blast_matrices_dir}" >> env.sh
echo "export rmblast_dir=${blast_dir}" >> env.sh
echo "export snap_dir=${snap_dir}" >> env.sh
echo "export ZOE=${snap_dir}/Zoe" >> env.sh
echo "export rapsearch_dir=${rapsearch_dir}" >> env.sh
echo "export trnascan_dir=${trnascan_dir}" >> env.sh
echo "export snoscan_dir=${snoscan_dir}" >> env.sh
echo "export repeatmasker_dir=${repeatmasker_dir}" >> env.sh
echo "export trf_dir=${trf_dir}" >> env.sh
echo "export reannotate_dir=${reannotate_dir}" >> env.sh
echo "export clustalw_dir=${clustalw_dir}" >> env.sh
echo "export muscle_dir=${muscle_dir}" >> env.sh
echo "export hmmer_dir=${hmmer_dir}" >> env.sh
echo "export bamtools_dir=${bamtools_dir}" >> env.sh
echo "export augustus_dir=${augustus_dir}" >> env.sh
echo "export AUGUSTUS_CONFIG_PATH=${build_dir}/augustus-${AUGUSTUS_VERSION}/config" >> env.sh
echo "export evm_dir=${evm_dir}" >> env.sh
echo "export EVM_HOME=${evm_dir}" >> env.sh
echo "export maker_dir=${maker_dir}" >> env.sh
echo "export proteinortho_dir=${proteinortho_dir}" >> env.sh
echo "export gatk_dir=${gatk_dir}" >> env.sh
echo "export ucsc_dir=${ucsc_dir}" >> env.sh
echo "export miniconda2_dir=${miniconda2_dir}" >> env.sh
echo "export conda_pacbio_dir=${conda_pacbio_dir}" >> env.sh
echo "export nanopolish_dir=${nanopolish_dir}" >> env.sh
echo "export parallel_dir=${parallel_dir}" >> env.sh

######### for Mfannot  ###########
echo "export emboss_dir=${emboss_dir}" >> env.sh
echo "export erpin_dir=${erpin_dir}" >> env.sh
echo "export tbl2asn_dir=${tbl2asn_dir}" >> env.sh
echo "export pirobject_dir=${pirobject_dir}" >> env.sh
echo "export pirmodels_dir=${pirmodels_dir}" >> env.sh
echo "export flip_dir=${flip_dir}" >> env.sh
echo "export umac_dir=${umac_dir}" >> env.sh
echo "export hmmsearchwc_dir=${hmmsearchwc_dir}" >> env.sh
echo "export rnafinder_dir=${rnafinder_dir}" >> env.sh
echo "export mf2sqn_dir=${mf2sqn_dir}" >> env.sh
echo "export grab_fasta_dir=${grab_fasta_dir}" >> env.sh
echo "export mfannot_data_dir=${mfannot_data_dir}" >> env.sh
echo "export mfannot_dir=${mfannot_dir}" >> env.sh


echo ""
echo "uncompress large supporting files ..."
gunzip $LRSDAY_HOME/data/Proteome_DB_for_annotation.CDhit_I95.fa.gz
gunzip $LRSDAY_HOME/data/SGDref.PoFF.ffn.gz
gunzip $LRSDAY_HOME/data/SGDref.PoFF.faa.gz
gunzip $LRSDAY_HOME/data/te_proteins.fasta.gz
gunzip $LRSDAY_HOME/Example_Outputs/SK1.assembly.final.fa.gz
gunzip $LRSDAY_HOME/Example_Outputs/SK1.final.gff3.gz
gunzip $LRSDAY_HOME/Example_Outputs/SK1.final.trimmed_cds.fa.gz
gunzip $LRSDAY_HOME/Example_Outputs/SK1.final.cds.fa.gz
gunzip $LRSDAY_HOME/Example_Outputs/SK1.final.pep.fa.gz
gunzip $LRSDAY_HOME/Example_Outputs/SK1.assembly.final.filter.mummer2vcf.SNP.vcf.gz
gunzip $LRSDAY_HOME/Example_Outputs/SK1.assembly.final.filter.mummer2vcf.INDEL.vcf.gz
echo "done!"
echo ""
echo ""
echo "#################### IMPORTANT !!! #######################"
echo ""
echo "1) Automatic dependencies installation finished! Please run \"source env.sh\" to load installation paths."
echo ""
echo "2) The installed UCSC Uitilites (i.e. blat, faSplit, pslSort, pslScore, and pslCDnaFilter) are free for nonprofit use."
echo "A license is required for commercial use."
echo "For information about commercial licensing of the Genome Browser software, see http://genome.ucsc.edu/license/."
echo ""
echo "3) The installed CAP3 is free for nonprofit use."
echo "CAP3 is not available for commercial use without a license. "
echo "If wishing to license CAP3 for commercial use, please contact Robin Kolehmainen at Michigan Tech by email at rakolehm@mtu.edu."
echo "Michigan Tech handles licensing agreements on PCAP for Iowa State."
echo ""
echo "4) The installed MAKER is free for nonprofit use under either the Artistic License 2.0 developed by the Perl Foundation"
echo "or the GNU General Public License developed by the Free Software Foundation."
echo "MAKER is not available for commercial use without a license."
echo "If wishing to license MAKER for commercial use, please contact Aaron Duffy at University of Utah TVC by email at Aaron.Duffy@tvc.utah.edu."
echo ""
echo "5) RepeatMasker need manual configuration."
echo "Please refer to our latest LRSDAY Manual for proper configuration of RepeatMasker."
echo "When making configuration for RepeatMasker, please use $trf_dir for configuring TRF and use $rmblast_dir for configuring rmblast."
echo ""
echo "6) The RepeatMasker dependency library Repbase (https://www.girinst.org/repbase/) will be switched to subscription-based model soon. Please consider subscribing this library for keeping up with the latest version of their product."
echo ""

echo "#########################################################"


############################
# checking Bash exit status
if [[ $? -eq 0 ]]
then
    echo ""
    echo "LRSDAY message: This bash script has been successfully processed! :)"
    echo ""
    echo ""
    exit 0
fi
############################
