#!/bin/bash
# last update: 2019/03/07

set -e -o pipefail

LRSDAY_HOME=$(pwd)
BUILD="build"

if [ -z "$MAKE_JOBS" ]
then
    echo "Defaulting to 2 concurrent jobs when executing make. Override with MAKE_JOBS=<NUM>"
    MAKE_JOBS=2
fi

if [ ! -z "$INSTALL_DEPS" ]; then
    echo "Installing LRSDAY build dependencies for Debian/Ubuntu."
    echo "sudo privileges are required and you will be prompted to enter your password"
    sudo apt-get update
    xargs -a debiandeps sudo apt-get install -y
fi

SRA_VERSION="2.9.2" # released on 2018.09.26
SRA_DOWNLOAD_URL="https://ftp-trace.ncbi.nlm.nih.gov/sra/sdk/${SRA_VERSION}/sratoolkit.${SRA_VERSION}-centos_linux64.tar.gz"

PORECHOP_VERSION="0.2.4" # 
PORECHOP_GITHUB_COMMIT_VERSION="109e437" # committed on 2018.10.19
PORECHOP_DOWNLOAD_URL="https://github.com/rrwick/Porechop.git"

FILTLONG_VERSION="0.2.0" #
FILTLONG_GITHUB_COMMIT_VERSION="d1bb46d" # committed on 2018.05.11
FILTLONG_DOWNLOAD_URL="https://github.com/rrwick/Filtlong.git"

MINIMAP2_VERSION="2.16" # released on 2019.2.28
MINIMAP2_DOWNLOAD_URL="https://github.com/lh3/minimap2/releases/download/v${MINIMAP2_VERSION}/minimap2-${MINIMAP2_VERSION}_x64-linux.tar.bz2"

CANU_VERSION="1.8" # released on 2018.10.23
CANU_DOWNLOAD_URL="https://github.com/marbl/canu/releases/download/v${CANU_VERSION}/canu-${CANU_VERSION}.Linux-amd64.tar.xz"

FLYE_VERSION="2.4.1" # released on 2019.03.07
FLYE_DOWNLOAD_URL="https://github.com/fenderglass/Flye/archive/${FLYE_VERSION}.tar.gz"

WTDBG2_VERSION="2.3" # 
WTDBG2_GITHUB_COMMIT_VERSION="59a39a6" # committed on 2019.03.06
WTDBG2_DOWNLOAD_URL="https://github.com/ruanjue/wtdbg2.git"

SMARTDENOVO_VERSION="" # 
SMARTDENOVO_GITHUB_COMMIT_VERSION="5cc1356" # committed on 2018.02.19
SMARTDENOVO_DOWNLOAD_URL="https://github.com/ruanjue/smartdenovo"

RAGOUT_VERSION="2.1.1" # released on 2018.07.30
RAGOUT_DOWNLOAD_URL="https://github.com/fenderglass/Ragout/archive/${RAGOUT_VERSION}.tar.gz"
# GUPPY_VERSION="2.3.5" # released on 2019.02.26
# QUAST_VERSION="5.0.1" # one of its dependency needs "csh" to be pre-installed

HDF_VERSION="1.10.1" # 
HDF_VERSION_prefix=${HDF_VERSION%.*}
HDF_DOWNLOAD_URL="https://support.hdfgroup.org/ftp/HDF5/releases/hdf5-${HDF_VERSION_prefix}/hdf5-${HDF_VERSION}/src/hdf5-${HDF_VERSION}.tar.gz"

SONLIB_VERSION="" # 
SONLIB_GITHUB_COMMIT_VERSION="1afbd97" # committed on 2017.08.09
SONLIB_DOWNLOAD_URL="https://github.com/benedictpaten/sonLib.git"

HAL_VERSION="" # not available, so we use the github comit hash below for version control
HAL_GITHUB_COMMIT_VERSION="a2ad656" # committed on 2017.09.09
HAL_DOWNLOAD_URL="https://github.com/glennhickey/hal.git"

MUMMER_VERSION="4.0.0beta2" # released on 2017.10.14
MUMMER_DOWNLOAD_URL="https://github.com/gmarcais/mummer/releases/download/v${MUMMER_VERSION}/mummer-${MUMMER_VERSION}.tar.gz"

GNUPLOT_VERSION="4.6.6" # released on 2015.02.18
GNUPLOT_DOWNLOAD_URL="https://sourceforge.net/projects/gnuplot/files/gnuplot/${GNUPLOT_VERSION}/gnuplot-${GNUPLOT_VERSION}.tar.gz"

BEDTOOLS_VERSION="2.27.1" # released on 2017.12.14
BEDTOOLS_DOWNLOAD_URL="https://github.com/arq5x/bedtools2/releases/download/v${BEDTOOLS_VERSION}/bedtools-${BEDTOOLS_VERSION}.tar.gz"

SPADES_VERSION="3.13.0" # released on 2018.10.16
SPADES_DOWNLOAD_URL="http://cab.spbu.ru/files/release${SPADES_VERSION}/SPAdes-${SPADES_VERSION}-Linux.tar.gz"

PRODIGAL_VERSION="2.6.3" # released on 2016.02.12
PRODIGAL_DOWNLOAD_URL="https://github.com/hyattpd/Prodigal/archive/v${PRODIGAL_VERSION}.tar.gz"

CAP_VERSION="" # see http://seq.cs.iastate.edu/cap3.html
CAP_DOWNLOAD_URL="http://seq.cs.iastate.edu/CAP3/cap3.linux.x86_64.tar"

BWA_VERSION="0.7.17" # released on 2017.10.23
BWA_DOWNLOAD_URL="http://downloads.sourceforge.net/project/bio-bwa/bwa-${BWA_VERSION}.tar.bz2"

SAMTOOLS_VERSION="1.9" # released on 2018.07.18
SAMTOOLS_DOWNLOAD_URL="https://github.com/samtools/samtools/releases/download/${SAMTOOLS_VERSION}/samtools-${SAMTOOLS_VERSION}.tar.bz2"

CIRCLATOR_VERSION="1.5.5" # released on 2018.01.31
CIRCLATOR_DOWNLOAD_URL="https://github.com/sanger-pathogens/circlator/archive/v${CIRCLATOR_VERSION}.tar.gz"

TRIMMOMATIC_VERSION="0.38" # 
TRIMMOMATIC_DOWNLOAD_URL="http://www.usadellab.org/cms/uploads/supplementary/Trimmomatic/Trimmomatic-${TRIMMOMATIC_VERSION}.zip"

#GATK_VERSION="3.6-6" #
#GATK_DOWNLOAD_URL="https://github.com/broadgsa/gatk/archive/${GATK_VERSION}.tar.gz"

PICARD_VERSION="2.18.23" # released on 2019.02.25
PICARD_DOWNLOAD_URL="https://github.com/broadinstitute/picard/releases/download/${PICARD_VERSION}/picard.jar"

PILON_VERSION="1.23" # released on 2018.11.27
PILON_DOWNLOAD_URL="https://github.com/broadinstitute/pilon/releases/download/v${PILON_VERSION}/pilon-${PILON_VERSION}.jar"

EXONERATE_VERSION="2.2.0" # 
EXONERATE_DOWNLOAD_URL="http://ftp.ebi.ac.uk/pub/software/vertebrategenomics/exonerate/exonerate-${EXONERATE_VERSION}-x86_64.tar.gz"

BLAST_VERSION="2.2.31" # 
RMBLAST_VERSION="2.2.28" # 
BLAST_DOWNLOAD_URL="ftp://ftp.ncbi.nlm.nih.gov/blast/executables/blast+/${BLAST_VERSION}/ncbi-blast-${BLAST_VERSION}+-x64-linux.tar.gz"
RMBLAST_DOWNLOAD_URL="ftp://ftp.ncbi.nlm.nih.gov/blast/executables/rmblast/${RMBLAST_VERSION}/ncbi-rmblastn-${RMBLAST_VERSION}-x64-linux.tar.gz"

SNAP_VERSION="" # 
SNAP_GITHUB_COMMIT_VERSION="a89d68e" # committed on 2017.05.18
SNAP_DOWNLOAD_URL="https://github.com/KorfLab/SNAP"

RAPSEARCH_VERSION="2.24" #
RAPSEARCH_DOWNLOAD_URL="https://sourceforge.net/projects/rapsearch2/files/RAPSearch${RAPSEARCH_VERSION}_64bits.tar.gz"

TRNASCAN_VERSION="1.3.1" #
TRNASCAN_DOWNLOAD_URL="http://eddylab.org/software/tRNAscan-SE/tRNAscan-SE.tar.gz"

SNOSCAN_VERSION="0.9.1" #
SNOSCAN_DOWNLOAD_URL="http://eddylab.org/software/snoscan/snoscan.tar.gz"

REPEATMASKER_VERSION="open-4-0-7" #
REPEATMASKER_DOWNLOAD_URL="http://repeatmasker.org/RepeatMasker-${REPEATMASKER_VERSION}.tar.gz"

TRF_VERSION="409" #
TRF_DOWNLOAD_URL="http://tandem.bu.edu/trf/downloads/trf${TRF_VERSION}.linux64"

REANNOTATE_VERSION="17.03.2015-LongQueryName"
REANNOTATE_DOWNLOAD_URL="https://github.com/yjx1217/REannotate_LongQueryName/archive/version_${REANNOTATE_VERSION}.tar.gz"

CLUSTALW_VERSION="2.1" #
CLUSTALW_DOWNLOAD_URL="http://www.clustal.org/download/${CLUSTALW_VERSION}/clustalw-${CLUSTALW_VERSION}.tar.gz"

MUSCLE_VERSION="3.8.31" #
MUSCLE_DOWNLOAD_URL="http://www.drive5.com/muscle/downloads${MUSCLE_VERSION}/muscle${MUSCLE_VERSION}_i86linux64.tar.gz"

HMMER_VERSION="3.2.1" # released on 2018.06.13
#HMMER_DOWNLOAD_URL="http://eddylab.org/software/hmmer3/${HMMER_VERSION}/hmmer-${HMMER_VERSION}-linux-intel-x86_64.tar.gz"
HMMER_DOWNLOAD_URL="http://eddylab.org/software/hmmer/hmmer-${HMMER_VERSION}.tar.gz"

BAMTOOLS_VERSION="2.4.2" # released on 2017.11.02
BAMTOOLS_DOWNLOAD_URL="https://github.com/pezmaster31/bamtools/archive/v${BAMTOOLS_VERSION}.tar.gz"

AUGUSTUS_VERSION="3.2.3" # 
#AUGUSTUS_GITHUB_COMMIT_VERSION="79960c5"
AUGUSTUS_DOWNLOAD_URL="http://bioinf.uni-greifswald.de/augustus/binaries/old/augustus-${AUGUSTUS_VERSION}.tar.gz"
#AUGUSTUS_DOWNLOAD_URL="https://github.com/Gaius-Augustus/Augustus.git"

EVM_VERSION="1.1.1" # released on 2015.07.03
EVM_DOWNLOAD_URL="https://github.com/EVidenceModeler/EVidenceModeler/archive/v${EVM_VERSION}.tar.gz"

PROTEINORTHO_VERSION="5.16b" # released on 2017.09.22
PROTEINORTHO_DOWNLOAD_URL="https://www.bioinf.uni-leipzig.de/Software/proteinortho/proteinortho_v${PROTEINORTHO_VERSION}.tar.gz"

MAKER_VERSION="3.00.0-beta" #
MAKER_DOWNLOAD_URL="http://topaz.genetics.utah.edu/maker_downloads/static/maker-${MAKER_VERSION}.tgz"

MINICONDA2_VERSION="4.5.11" #
MINICONDA2_DOWNLOAD_URL="https://repo.continuum.io/miniconda/Miniconda2-${MINICONDA2_VERSION}-Linux-x86_64.sh"
PB_ASSEMBLY_VERSION="0.0.2" #
BAX2BAM_VERSION="0.0.9" #

NANOPOLISH_VERSION="0.11.0" #
NANOPOLISH_GITHUB_COMMIT_VERSION="8b5e605" # commited on 2019.02.22 
NANOPOLISH_DOWNLOAD_URL="https://github.com/jts/nanopolish.git"

PARALLEL_VERSION="20180722" # released on 2018.07.22
PARALLEL_DOWNLOAD_URL="http://ftp.gnu.org/gnu/parallel/parallel-${PARALLEL_VERSION}.tar.bz2"
# for MFannot
EMBOSS_VERSION="6.5.7" # released on 2012.07.25
EMBOSS_VERSION_prefix="${EMBOSS_VERSION:0:3}"
EMBOSS_DOWNLOAD_URL="ftp://emboss.open-bio.org/pub/EMBOSS/old/${EMBOSS_VERSION_prefix}.0/EMBOSS-${EMBOSS_VERSION}.tar.gz"

ERPIN_VERSION="5.5.4" # 
ERPIN_DOWNLOAD_URL="http://rna.igmors.u-psud.fr/download/Erpin/erpin${ERPIN_VERSION}.serv.tar.gz"

TBL2ASN_VERSION="" #
TBL2ASN_DOWNLOAD_URL="ftp://ftp.ncbi.nih.gov/toolbox/ncbi_tools/converters/by_program/tbl2asn/linux64.tbl2asn.gz"

PIROBJECT_VERSION="1.19" #
PIROBJECT_DOWNLOAD_URL="https://github.com/prioux/PirObject/archive/v${PIROBJECT_VERSION}.tar.gz"
PIRMODELS_GITHUB_COMMIT_VERSION="6b223ec" # committed on 2016.08.30
PIRMODELS_DOWNLOAD_URL="https://github.com/BFL-lab/PirModels.git"

FLIP_GITHUB_COMMIT_VERSION="00a57cb" # committed on 2016.04.07
FLIP_DOWNLOAD_URL="https://github.com/BFL-lab/Flip.git"

UMAC_GITHUB_COMMIT_VERSION="cae618e" # committed on 2016.08.30
UMAC_DOWNLOAD_URL="https://github.com/BFL-lab/Umac.git"

HMMSEARCHWC_GITHUB_COMMIT_VERSION="9e3b461" # committed on 2016.11.05
HMMSEARCHWC_DOWNLOAD_URL="https://github.com/BFL-lab/HMMsearchWC.git"

RNAFINDER_GITHUB_COMMIT_VERSION="579dc58" # committed on 2016.12.07
RNAFINDER_DOWNLOAD_URL="https://github.com/BFL-lab/RNAfinder.git"

MF2SQN_GITHUB_COMMIT_VERSION="6faf9f4" # committed on 2016.12.07
MF2SQN_DOWNLOAD_URL="https://github.com/BFL-lab/Mf2sqn.git"

GRAB_FASTA_GITHUB_COMMIT_VERSION="accd32d" # committed on 2017.02.14
GRAB_FASTA_DOWNLOAD_URL="https://github.com/BFL-lab/grab-fasta.git"

# for MFannot
MFANNOT_DATA_GITHUB_COMMIT_VERSION="b039ac5" # committed on 2016.12.07
MFANNOT_VERSION="1.35" #
MFANNOT_GITHUB_COMMIT_VERSION="6472b97" # committed on 2018.10.31
MFANNOT_DATA_DOWNLOAD_URL="https://github.com/BFL-lab/MFannot_data.git"
MFANNOT_DOWNLOAD_URL="https://github.com/BFL-lab/MFannot.git"

# downloading URLs for dependencies
# GUPPY_DOWNLOAD_URL="https://mirror.oxfordnanoportal.com/software/analysis/ont-guppy-cpu_${GUPPY_VERSION}_linux64.tar.gz"
# QUAST_DOWNLOAD_URL="https://downloads.sourceforge.net/project/quast/quast-${QUAST_VERSION}.tar.gz"


# UCSC Utilities
BLAT_DOWNLOAD_URL="http://hgdownload.soe.ucsc.edu/admin/exe/linux.x86_64/blat/blat"
FASPLIT_DOWNLOAD_URL="http://hgdownload.soe.ucsc.edu/admin/exe/linux.x86_64/faSplit"
PSLSORT_DOWNLOAD_URL="http://hgdownload.soe.ucsc.edu/admin/exe/linux.x86_64/pslSort"
PSLSCORE_DOWNLOAD_URL="http://hgdownload.soe.ucsc.edu/admin/exe/linux.x86_64/pslScore"
PSLCDNAFILTER_DOWNLOAD_URL="http://hgdownload.soe.ucsc.edu/admin/exe/linux.x86_64/pslCDnaFilter"

# Create the $BUILD directory for dependency installation
echo ""
echo "Create $BUILD directory if it does not already exist"
mkdir -p $BUILD
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

download_and_extract() {
    url=$1
    download_location=$2
    echo "Downloading $url to $download_location"
    wget -nv --no-check-certificate $url -O $download_location
    download $url $download_location
    if [[ $download_location =~ \.bz2$ ]]; then
        extract_command="tar -xjf"
    elif [[ $download_location =~ \.xz$ || $download_location =~ \.tar$ ]]; then
        extract_command="tar -xf"
    else
        extract_command="tar -zxf"
    fi
    $($extract_command $download_location)
    rm $download_location

}

check_installed() {
    if [ -e "$1/installed" ]; then
        echo "installed"
    else
        echo ""
    fi
}

note_installed() {
    touch "$1/installed"
}

# ---------- set Perl & Python environment variables -------------
PYTHONPATH="$build_dir"
PERL5LIB="$build_dir:$PERL5LIB"
PERL5LIB="$build_dir/cpanm/perlmods/lib/perl5:$PERL5LIB"
cpanm_dir=$build_dir/cpanm

if [ ! -e "$build_dir/cpanm" ]; then
    mkdir -p  $build_dir/cpanm
    cd $cpanm_dir
    wget -nv --no-check-certificate -O - https://cpanmin.us/ > cpanm
    chmod +x cpanm
    mkdir -p perlmods
fi

# Testing the packages when they are getting installed slows them down. Is this absolutely required?
xargs -a "$LRSDAY_HOME/perldeps" $cpanm_dir/cpanm -l $cpanm_dir/perlmods --notest --skip-installed || (echo "Could not install some Perl modules. See logs for failures" && exit 1)
if [ ! -z "$USE_POSTGRES" ]; then
    # need $POSTGRES_HOME pre-defined
    $cpanm_dir/cpanm -l $cpanm_dir/perlmods --notest --skip-installed DBD::Pg
else
    $cpanm_dir/cpanm -l $cpanm_dir/perlmods --notest --skip-installed DBD::SQLite@1.54
fi

# ------------- SRA Toolkit -------------------
sra_dir="$build_dir/sratoolkit.${SRA_VERSION}-centos_linux64/bin"
if [ -z $(check_installed $sra_dir) ]; then
    cd $build_dir
    echo "Download SRAtoolkit-v${SRA_VERSION}"
    download_and_extract $SRA_DOWNLOAD_URL sratoolkit.${SRA_VERSION}-centos_linux64.tar.gz
    note_installed $sra_dir
fi

# --------------- Porechop ------------------
porechop_dir="$build_dir/Porechop/py3_virtualenv_porechop/bin"
if [ -z $(check_installed $porechop_dir) ]; then
    cd $build_dir
    echo "Download Porechop-v${PORECHOP_VERSION}"
    git clone $PORECHOP_DOWNLOAD_URL
    cd Porechop
    git checkout -f -q $PORECHOP_GITHUB_COMMIT_VERSION
    virtualenv -p $(which python3) py3_virtualenv_porechop
    source py3_virtualenv_porechop/bin/activate
    py3_virtualenv_porechop/bin/python3 ./setup.py install
    deactivate
    note_installed $porechop_dir
fi


# --------------- Filtlong ------------------
filtlong_dir="$build_dir/Filtlong/bin"
if [ -z $(check_installed $filtlong_dir) ]; then
    cd $build_dir
    echo "Download Filtlong-v${FILTLONG_VERSION}"
    git clone $FILTLONG_DOWNLOAD_URL
    cd Filtlong
    git checkout -f -q $FILTLONG_GITHUB_COMMIT_VERSION
    make -j $MAKE_JOBS
    note_installed $filtlong_dir
fi

# --------------- minimap2 ------------------
minimap2_dir="$build_dir/minimap2-${MINIMAP2_VERSION}_x64-linux"
if [ -z $(check_installed $minimap2_dir) ]; then
    cd $build_dir
    echo "Download minimap2-v${MINIMAP2_VERSION}"
    download_and_extract $MINIMAP2_DOWNLOAD_URL "minimap2-${MINIMAP2_VERSION}.tar.bz2"
    note_installed $minimap2_dir
fi

# ------------- Canu -------------------
canu_dir="$build_dir/canu-${CANU_VERSION}/Linux-amd64/bin"
if [ -z $(check_installed $canu_dir) ]; then
    cd $build_dir
    echo "Download Canu-v${CANU_VERSION}"
    download_and_extract $CANU_DOWNLOAD_URL "canu-${CANU_VERSION}.tar.xz"
    canu_dir="$build_dir/canu-${CANU_VERSION}"
    cd $canu_dir
    canu_dir="$build_dir/canu-${CANU_VERSION}/Linux-amd64/bin"
    cd $canu_dir
    ln -s $minimap2_dir/minimap2 .
    note_installed $canu_dir
fi

# ------------- Flye -------------------
flye_dir="$build_dir/Flye-${FLYE_VERSION}/bin"
if [ -z $(check_installed $flye_dir) ]; then
    cd $build_dir
    echo "Download Flye-v${FLYE_VERSION}"
    download_and_extract $FLYE_DOWNLOAD_URL "Flye-${FLYE_VERSION}.tar.gz"
    cd Flye-${FLYE_VERSION}
    python2 setup.py build
    note_installed $flye_dir
fi

# --------------- wtdbg2 ------------------
wtdbg2_dir="$build_dir/wtdbg2"
if [ -z $(check_installed $wtdbg2_dir) ]; then
    cd $build_dir
    echo "Download wtdbg2-v${WTDBG2_VERSION}"
    git clone $WTDBG2_DOWNLOAD_URL
    cd wtdbg2
    git checkout -f -q $WTDBG2_GITHUB_COMMIT_VERSION
    C_INCLUDE_PATH="" 
    make -j $MAKE_JOBS
    note_installed $wtdbg2_dir
fi

# --------------- smartdenovo ------------------
smartdenovo_dir="$build_dir/smartdenovo"
if [ -z $(check_installed $smartdenovo_dir) ]; then
    cd $build_dir
    echo "Download smartdenovo-v${SMARTDENOVO_VERSION}"
    git clone $SMARTDENOVO_DOWNLOAD_URL
    cd smartdenovo
    git checkout -f -q $SMARTDENOVO_GITHUB_COMMIT_VERSION
    cp wtlay.h wtlay.h.bk
    cat wtlay.h.bk |sed s/inline//g > wtlay.h
    C_INCLUDE_PATH="" 
    make -j $MAKE_JOBS
    cp $LRSDAY_HOME/misc/smartdenovo_customized.pl .
    note_installed $smartdenovo_dir
fi

# --------------- Ragout ------------------
ragout_dir="$build_dir/Ragout-${RAGOUT_VERSION}/bin"
if [ -z $(check_installed $ragout_dir) ]; then
    cd $build_dir
    echo "Download Ragout-v${RAGOUT_VERSION}"
    download_and_extract $RAGOUT_DOWNLOAD_URL Ragout-${RAGOUT_VERSION}.tar.gz
    cd Ragout-${RAGOUT_VERSION}
    python2 setup.py build
    python2 ./scripts/install-sibelia.py
    note_installed $ragout_dir
fi

# --------------- HDF ------------------
hdf_dir="$build_dir/hdf5/bin"
h5prefix="-prefix=$build_dir/hdf5"
if [ -z $(check_installed $hdf_dir) ]; then
    cd $build_dir
    echo "Download HDF-v${HDF_VERSION}"
    download_and_extract $HDF_DOWNLOAD_URL "hdf5-${HDF_VERSION}.tar.gz"
    mkdir hdf5
    cd hdf5-${HDF_VERSION}
    ./configure --enable-cxx  --prefix $build_dir/hdf5
    make -j $MAKE_JOBS
    make install
    note_installed $hdf_dir
fi
PATH="$hdf_dir:${PATH}"

# --------------- sonLib ------------------
if [ -z $(check_installed "sonLib") ]; then
    cd $build_dir
    git clone $SONLIB_DOWNLOAD_URL
    cd sonLib
    git checkout -f -q $SONLIB_GITHUB_COMMIT_VERSION
    make -j $MAKE_JOBS
    note_installed "$build_dir/sonLib"
fi


# ---------------- HAL -------------------
hal_dir="$build_dir/hal/bin"
if [ -z $(check_installed $hal_dir) ]; then
    cd $build_dir
    echo "Download HAL-v${HAL_VERSION}"
    git clone $HAL_DOWNLOAD_URL
    cd hal
    git checkout -f -q $HAL_GITHUB_COMMIT_VERSION
    # Compiling with multiple workers results in build failures
    # Possible race conditions / build order issues within the project?
    make
    note_installed $hal_dir
fi

# --------------- samtools -----------------
samtools_dir="$build_dir/samtools-${SAMTOOLS_VERSION}"
if [ -z $(check_installed $samtools_dir) ]; then
    cd $build_dir
    echo "Download samtools-v${SAMTOOLS_VERSION}"
    download_and_extract $SAMTOOLS_DOWNLOAD_URL "samtools-${SAMTOOLS_VERSION}.tar.bz2"
    cd $samtools_dir
    C_INCLUDE_PATH=""
    ./configure --without-curses;
    make -j $MAKE_JOBS
    note_installed $samtools_dir
fi

# --------------- gnuplot ------------------
gnuplot_dir="$build_dir/gnuplot-${GNUPLOT_VERSION}/bin"
if [ -z $(check_installed $gnuplot_dir) ]; then
    cd $build_dir
    echo "Download gnuplot-v${GNUPLOT_VERSION}"
    download_and_extract $GNUPLOT_DOWNLOAD_URL "gnuplot-${GNUPLOT_VERSION}.tar.gz"
    cd "$build_dir/gnuplot-${GNUPLOT_VERSION}"
    ./configure --prefix="$build_dir/gnuplot-${GNUPLOT_VERSION}" --disable-wxwidgets
    make
    make install
    note_installed $gnuplot_dir
fi
PATH="$gnuplot_dir:${PATH}"

# # --------------- Guppy --------------------
# cd $build_dir
# echo "Download Guppy-v${GUPPY_VERSION}"
# download $GUPPY_DOWNLOAD_URL "ont-guppy-cpu_${GUPPY_VERSION}_linux64.tar.gz"
# tar -xzf ont-guppy-cpu_${GUPPY_VERSION}_linux64.tar.gz
# guppy_dir="$build_dir/ont-guppy-cpu/bin"
# rm "ont-guppy-cpu_${GUPPY_VERSION}_linux64.tar.gz"

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
mummer_dir="$build_dir/mummer-${MUMMER_VERSION}"
if [ -z $(check_installed $mummer_dir) ]; then
    cd $build_dir
    echo "Download mummer-v${MUMMER_VERSION}"
    download_and_extract $MUMMER_DOWNLOAD_URL "mummer-${MUMMER_VERSION}.tar.gz"
    echo "$mummer_dir"
    cd $mummer_dir
    ./configure
    make -j $MAKE_JOBS
    note_installed $mummer_dir
fi
PATH="$mummer_dir:${PATH}"

# --------------- bedtools ------------------
bedtools_dir="$build_dir/bedtools2/bin"
if [ -z $(check_installed $bedtools_dir) ]; then
    cd $build_dir
    echo "Download bedtools-v${BEDTOOLS_VERSION}"
    download_and_extract $BEDTOOLS_DOWNLOAD_URL "bedtools-${BEDTOOLS_VERSION}.tar.gz"
    cd "$build_dir/bedtools2"
    make -j $MAKE_JOBS
    note_installed $bedtools_dir
fi

# --------------- SPAdes ------------------
spades_dir="$build_dir/SPAdes-${SPADES_VERSION}-Linux/bin"
if [ -z $(check_installed $spades_dir) ]; then
    cd $build_dir
    echo "Download SPAdes-v${SPADES_VERSION}"
    # This sometimes times out, could this be throttling from the server?
    # Allow administrator to manually upload a copy of the tarball for installation
    if [ ! -e "SPAdes-${SPADES_VERSION}-Linux.tar.gz" ]; then
        download_and_extract $SPADES_DOWNLOAD_URL "SPAdes-${SPADES_VERSION}-Linux.tar.gz"
    else
        tar -zxf "SPAdes-${SPADES_VERSION}-Linux.tar.gz"
    fi
    note_installed $spades_dir
fi

# --------------- Prodigal ------------------
prodigal_dir="$build_dir/Prodigal-${PRODIGAL_VERSION}"
if [$(check_installed $prodigal_dir) ]; then
    cd $build_dir
    echo "Download Prodigal-v${PRODIGAL_VERSION}"
    download_and_extract $PRODIGAL_DOWNLOAD_URL "v${PRODIGAL_VERSION}.tar.gz"
    cd $prodigal_dir
    make -j $MAKE_JOBS
    note_installed $prodigal_dir
fi

# --------------- CAP3 ------------------
cap_dir="$build_dir/CAP3"
if [ -z $(check_installed $cap_dir) ]; then
    cd $build_dir
    echo "Download CAP3-v${CAP3_VERSION}"
    download_and_extract $CAP_DOWNLOAD_URL "cap3.linux.x86_64.tar"
    note_installed $cap_dir
fi

# ------------- BWA -------------------
bwa_dir="$build_dir/bwa-${BWA_VERSION}"
if [$(check_installed $bwa_dir) ]; then
    cd $build_dir
    echo "Download BWA-v${BWA_VERSION}"
    download_and_extract $BWA_DOWNLOAD_URL "bwa-${BWA_VERSION}.tar.bz2"
    cd $bwa_dir
    make -j $MAKE_JOBS
    note_installed $bwa_dir
fi

# --------------- Circlator ------------------
circlator_dir="$build_dir/py3_virtualenv_circlator/bin"
if [ -z $(check_installed $circlator_dir) ]; then
    cd $build_dir
    echo "Creating local virtual python3 environment and install Circlator-v${CIRCLATOR_VERSION}"
    virtualenv -p $(which python3) py3_virtualenv_circlator
    source py3_virtualenv_circlator/bin/activate
    py3_virtualenv_circlator/bin/pip3 install "circlator==${CIRCLATOR_VERSION}"
    deactivate
    note_installed $circlator_dir
fi

# --------------- Trimmomatic -----------------
trimmomatic_dir="$build_dir/Trimmomatic-${TRIMMOMATIC_VERSION}"
if [ -z $(check_installed $trimmomatic_dir) ]; then
    cd $build_dir
    echo "Download Trimmomatic-v${TRIMMOMATIC_VERSION}"
    download $TRIMMOMATIC_DOWNLOAD_URL "Trimmomatic-${TRIMMOMATIC_VERSION}.zip"
    unzip Trimmomatic-${TRIMMOMATIC_VERSION}.zip
    rm Trimmomatic-${TRIMMOMATIC_VERSION}.zip

    cd $trimmomatic_dir
    chmod 755 trimmomatic-${TRIMMOMATIC_VERSION}.jar
    ln -s trimmomatic-${TRIMMOMATIC_VERSION}.jar trimmomatic.jar 
    note_installed $trimmomatic_dir
fi

# --------------- Picard -----------------
picard_dir="$build_dir/Picard-v${PICARD_VERSION}"
if [ -z $(check_installed $picard_dir) ]; then
    cd $build_dir
    echo "Download Picard-v${PICARD_VERSION}"
    download $PICARD_DOWNLOAD_URL "picard.jar"
    mkdir Picard-v${PICARD_VERSION}

    mv picard.jar $picard_dir
    cd $picard_dir
    chmod 755 picard.jar
    note_installed $picard_dir
fi

# --------------- Pilon -----------------
pilon_dir="$build_dir/Pilon-v${PILON_VERSION}"
if [ -z $(check_installed $pilon_dir) ]; then
    cd $build_dir
    echo "Download Pilon-v${PILON_VERSION}"
    download $PILON_DOWNLOAD_URL "pilon-${PILON_VERSION}.jar"
    mkdir Pilon-v${PILON_VERSION}
    mv pilon-${PILON_VERSION}.jar $pilon_dir/pilon.jar
    cd $pilon_dir
    chmod 755 pilon.jar
    note_installed $pilon_dir
fi

# --------------- exonerate ------------------
exonerate_dir="$build_dir/exonerate-${EXONERATE_VERSION}-x86_64/bin"
if [ -z $(check_installed $exonerate_dir) ]; then
    cd $build_dir
    echo "Download exonerate-v${EXONERATE_VERSION}"
    download_and_extract $EXONERATE_DOWNLOAD_URL "exonerate-${EXONERATE_VERSION}-x86_64.tar.gz"
    note_installed $exonerate_dir
fi

# --------------- ncbi-blast+ ------------------
blast_dir="$build_dir/ncbi-blast-${BLAST_VERSION}+/bin"
blast_matrices_dir="$build_dir/ncbi-blast-${BLAST_VERSION}+/matrices"
if [ -z $(check_installed $blast_dir) ]; then
    cd $build_dir
    echo "Download ncbi-blast-v${BLAST_VERSION}"
    download_and_extract $BLAST_DOWNLOAD_URL "ncbi-blast-${BLAST_VERSION}+-x64-linux.tar.gz"
    cd "ncbi-blast-${BLAST_VERSION}+"
    mkdir matrices
    cd matrices
    wget -nv --no-check-certificate ftp://ftp.ncbi.nlm.nih.gov/blast/matrices/*
    note_installed $blast_dir
fi

# --------------- ncbi-rmblast ------------------
rmblast_dir="$build_dir/ncbi-rmblastn-${RMBLAST_VERSION}/bin"
if [ -z $(check_installed $rmblast_dir) ]; then
    cd $build_dir
    echo "Download ncbi-rmblastn-v${BLAST_VERSION}"
    download_and_extract $RMBLAST_DOWNLOAD_URL "ncbi-rmblastn-${RMBLAST_VERSION}-x64-linux.tar.gz"
    # copy rmblastn binary file to ncbi-blast+ directory for easy RepeatMasker configuration
    cp $rmblast_dir/rmblastn $blast_dir
    note_installed $rmblast_dir
fi

# --------------- snap ------------------
snap_dir="$build_dir/SNAP"
if [ -z $(check_installed $snap_dir) ]; then
    cd $build_dir
    echo "Download snap-v${SNAP_VERSION}"
    git clone $SNAP_DOWNLOAD_URL
    cd $snap_dir
    git checkout -f -q $SNAP_GITHUB_COMMIT_VERSION
    ZOE="$snap_dir/Zoe"
    cp $LRSDAY_HOME/misc/snap.c .  # temporary fix for snap with gcc-8
    make -j $MAKE_JOBS
    note_installed $snap_dir
fi

# --------------- RAPSearch2 ------------------
rapsearch_dir="$build_dir/RAPSearch${RAPSEARCH_VERSION}_64bits/bin"
if [ -z $(check_installed $rapsearch_dir) ]; then
    cd $build_dir
    echo "Download RAPsearch-v${RAPSEARCH_VERSION}"
    download_and_extract $RAPSEARCH_DOWNLOAD_URL "RAPSearch${RAPSEARCH_VERSION}_64bits.tar.gz"
    note_installed $rapsearch_dir
fi

# --------------- tRNAscan-SE ------------------
trnascan_root="$build_dir/tRNAscan-SE-${TRNASCAN_VERSION}"
trnascan_dir="$trnascan_root/bin"
PERL5LIB=$trnascan_dir:$PERL5LIB
if [ -z $(check_installed $trnascan_dir) ]; then
    cd $build_dir
    echo "Download tRNAscan-SE-v${TRNASCAN_VERSION}"
    download_and_extract $TRNASCAN_DOWNLOAD_URL "tRNAscan-SE-${TRNASCAN_VERSION}.tar.gz"
    cd $trnascan_root
    mkdir bin
    mkdir -p lib/tRNAscan-SE
    mkdir -p man
    cp $LRSDAY_HOME/misc/tRNAscan-SE.Makefile Makefile
    make -j $MAKE_JOBS BINDIR="$trnascan_root/bin" LIBDIR="$trnascan_root/lib/tRNAscan-SE" MANDIR="$trnascan_root"
    make install BINDIR="$trnascan_root/bin" LIBDIR="$trnascan_root/lib/tRNAscan-SE" MANDIR="$trnascan_root"
    note_installed $trnascan_dir
fi

# --------------- snoscan ------------------
snoscan_dir="$build_dir/snoscan-${SNOSCAN_VERSION}"
if [ -z $(check_installed $snoscan_dir) ]; then
    cd $build_dir
    echo "Download snoscan-v${SNOSCAN_VERSION}"
    download_and_extract $SNOSCAN_DOWNLOAD_URL "snoscan-${SNOSCAN_VERSION}.tar.gz"
    cd $snoscan_dir
    cd squid-1.5.11
    rm *.o
    make -j $MAKE_JOBS
    cd ..
    cp $LRSDAY_HOME/misc/snoscan.Makefile Makefile
    rm *.o
    make -j $MAKE_JOBS
    note_installed $snoscan_dir
fi

# --------------- RepeatMasker ------------------
repeatmasker_dir="$build_dir/RepeatMasker"
if [ -z $(check_installed $repeatmasker_dir) ]; then
    cd $build_dir
    echo "Download Repeatmasker-v${REPEATMASKER_VERSION}"
    download_and_extract $REPEATMASKER_DOWNLOAD_URL "RepeatMasker-${REPEATMASKER_VERSION}.tar.gz"
    cd $repeatmasker_dir
    echo "Download and setup RepBase library"
    REPBASE_VERSION="20170127"
    wget -nv --no-check-certificate https://github.com/yjx1217/RMRB/raw/master/RepBaseRepeatMaskerEdition-${REPBASE_VERSION}.tar.gz
    tar xzf RepBaseRepeatMaskerEdition-${REPBASE_VERSION}.tar.gz
    rm RepBaseRepeatMaskerEdition-${REPBASE_VERSION}.tar.gz
    cd .. 
    note_installed $repeatmasker_dir
fi

# --------------- TRF ------------------
trf_dir=$repeatmasker_dir
if [! -e "$repeatmasker_dir/trf"]; then
    cd $repeatmasker_dir
    echo "Download TRF-v${TRF_VERSION}"
    download $TRF_DOWNLOAD_URL "trf${TRF_VERSION}.linux64"
    mv trf${TRF_VERSION}.linux64 trf
    chmod 755 trf
fi


# --------------- REannotate ------------------
reannotate_dir="$build_dir/REannotate_LongQueryName-version_${REANNOTATE_VERSION}"
if [ -z $(check_installed $reannotate_dir) ]; then
    cd $build_dir
    echo "Download REannotate-v${REANNOTATE_VERSION}"
    download_and_extract $REANNOTATE_DOWNLOAD_URL "version_${REANNOTATE_VERSION}.tar.gz"
    cd $reannotate_dir
    chmod 755 REannotate_longname
    ln -s REannotate_longname REannotate
    note_installed $reannotate_dir
fi

# --------------- ClustalW ------------------
clustalw_dir="$build_dir/clustalw-${CLUSTALW_VERSION}/bin"
if [ -z $(check_installed $clustalw_dir) ]; then
    cd $build_dir
    echo "Download ClustalW-v${CLUSTALW_VERSION}"
    download_and_extract $CLUSTALW_DOWNLOAD_URL "clustalw-${CLUSTALW_VERSION}.tar.gz"
    cd clustalw-${CLUSTALW_VERSION}
    ./configure --prefix="$build_dir/clustalw-${CLUSTALW_VERSION}"
    make -j $MAKE_JOBS
    make install
    note_installed $clustalw_dir
fi

# --------------- MUSCLE ------------------
muscle_dir="$build_dir/muscle-${MUSCLE_VERSION}"
if [ -z $(check_installed $muscle_dir) ]; then
    cd $build_dir
    echo "Download MUSCLE-v${MUSCLE_VERSION}"
    download_and_extract $MUSCLE_DOWNLOAD_URL "muscle-${MUSCLE_VERSION}_i86linux64.tar.gz"
    mkdir muscle-${MUSCLE_VERSION}
    mv muscle${MUSCLE_VERSION}_i86linux64 ./muscle-${MUSCLE_VERSION}/muscle
    note_installed $muscle_dir
fi

# --------------- HMMER ------------------
hmmer_root="$build_dir/hmmer-${HMMER_VERSION}"
hmmer_dir="$build_dir/hmmer-${HMMER_VERSION}/bin"
if [ -z $(check_installed $hmmer_dir) ]; then
    cd $build_dir
    echo "Download hmmer-v${HMMER_VERSION}"
    download_and_extract $HMMER_DOWNLOAD_URL "hmmer-${HMMER_VERSION}-linux-intel-x86_64.tar.gz"
    cd $hmmer_root
    ./configure --prefix=$hmmer_root
    make -j $MAKE_JOBS
    make install
    cd easel
    make install
    note_installed $hmmer_dir
fi

# --------------- bamtools ------------------
bamtools_dir="$build_dir/bamtools-${BAMTOOLS_VERSION}/bin"
if [ -z $(check_installed $bamtools_dir) ]; then
    cd $build_dir
    echo "Download bamtools-v${BAMTOOLS_VERSION}"
    download_and_extract $BAMTOOLS_DOWNLOAD_URL "v${BAMTOOLS_VERSION}.tar.gz"
    cd $build_dir/bamtools-${BAMTOOLS_VERSION}
    mkdir build
    cd build
    cmake -DCMAKE_INSTALL_PREFIX="$build_dir/bamtools-${BAMTOOLS_VERSION}" ..
    make -j $MAKE_JOBS
    make install
    cd $build_dir/bamtools-${BAMTOOLS_VERSION}
    ln -sf lib lib64
    note_installed $bamtools_dir
fi

# --------------- Augustus ------------------
augustus_dir="$build_dir/augustus-${AUGUSTUS_VERSION}/bin"
if [ -z $(check_installed $augustus_dir) ]; then
    cd $build_dir
    echo "Download Augustus-v${AUGUSTUS_VERSION}"
    download_and_extract $AUGUSTUS_DOWNLOAD_URL "augustus-${AUGUSTUS_VERSION}.tar.gz"
    cd $build_dir/augustus-${AUGUSTUS_VERSION}/auxprogs/bam2hints/
    cp $LRSDAY_HOME/misc/bam2hints.Makefile Makefile
    cd $build_dir/augustus-${AUGUSTUS_VERSION}/auxprogs/filterBam/src/
    cp $LRSDAY_HOME/misc/filterBam.Makefile Makefile
    cd $build_dir/augustus-${AUGUSTUS_VERSION}
    make -j $MAKE_JOBS BAMTOOLS="$build_dir/bamtools-${BAMTOOLS_VERSION}"
    note_installed $augustus_dir
fi
export AUGUSTUS_CONFIG_PATH="$build_dir/augustus-${AUGUSTUS_VERSION}/config"

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
evm_dir="$build_dir/EVidenceModeler-${EVM_VERSION}"
if [ -z $(check_installed $evm_dir) ]; then
    cd $build_dir
    echo "Download EvidenceModeler-v${EVM_VERSION}"
    download_and_extract $EVM_DOWNLOAD_URL "v${EVM_VERSION}.tar.gz"
    note_installed $evm_dir
fi

# --------------- Proteinortho ------------------
proteinortho_dir="$build_dir/proteinortho_v${PROTEINORTHO_VERSION}"
if [ -z $(check_installed $proteinortho_dir) ]; then
    cd $build_dir
    echo "Download Proteinortho-v${PROTEINORTHO_VERSION}"
    download_and_extract $PROTEINORTHO_DOWNLOAD_URL "proteinortho_v${PROTEINORTHO_VERSION}.tar.gz"
    cp $LRSDAY_HOME/misc/proteinortho5_better_robustness.pl $proteinortho_dir/proteinortho5.pl
    note_installed $proteinortho_dir
fi

# --------------- GATK ------------------
gatk_dir="$build_dir/GATK3"
if [ -z $(check_installed $gatk_dir) ]; then
    cd $build_dir
    echo "Create GATK3 folder for users' manual installation"
    mkdir GATK3
    cd GATK3
    wget -nv --no-check-certificate https://ftp-trace.ncbi.nlm.nih.gov/sra/sdk/${SRA_VERSION}/GenomeAnalysisTK.jar
    chmod 755 GenomeAnalysisTK.jar
    note_installed $gatk_dir
fi


# --------------- MAKER -----------------
maker_dir="$build_dir/maker/bin"
if [ -z $(check_installed $maker_dir) ]; then
    cd $build_dir
    echo "Download MAKER"
    download_and_extract $MAKER_DOWNLOAD_URL "maker-v${MAKER_VERSION}.tgz"
    cd $build_dir/maker/src/
    cp $LRSDAY_HOME/misc/maker_Build.PL .
    echo "no"|perl maker_Build.PL
    ./Build install
    note_installed $maker_dir
fi

# --------------- UCSC Utilities -----------------
ucsc_dir="$build_dir/UCSC_Utilities"
if [ -z $(check_installed $ucsc_dir) ]; then
    cd $build_dir
    mkdir UCSC_Utilities
    cd $ucsc_dir
    download $BLAT_DOWNLOAD_URL "blat"
    download $FASPLIT_DOWNLOAD_URL "faSplit"
    download $PSLSORT_DOWNLOAD_URL "pslSort"
    download $PSLSCORE_DOWNLOAD_URL "pslScore"
    download $PSLCDNAFILTER_DOWNLOAD_URL "pslCDnaFilter"
    chmod 755 $ucsc_dir/*
    note_installed $ucsc_dir
fi

# ------------- Conda-PacBio --------------------
miniconda2_dir="$build_dir/miniconda2/bin"
conda_pacbio_dir=$build_dir/conda_pacbio_env/bin
if [ -z $(check_installed $miniconda2_dir) ]; then
    cd $build_dir
    download $MINICONDA2_DOWNLOAD_URL "Miniconda2-${MINICONDA2_VERSION}-Linux-x86_64.sh"
    bash Miniconda2-${MINICONDA2_VERSION}-Linux-x86_64.sh -b -p $build_dir/miniconda2
    #export PATH="$build_dir/miniconda2/bin:$PATH"
    $miniconda2_dir/conda config --add channels defaults
    $miniconda2_dir/conda config --add channels bioconda
    $miniconda2_dir/conda config --add channels conda-forge
    $miniconda2_dir/conda create -y -p $build_dir/conda_pacbio_env
    source $miniconda2_dir/activate $build_dir/conda_pacbio_env
    $miniconda2_dir/conda install -y hdf5=${HDF_VERSION}
    $miniconda2_dir/conda install -y -c bioconda samtools=${SAMTOOLS_VERSION} openssl=1.0
    $miniconda2_dir/conda install -y -c bioconda pb-assembly=${PB_ASSEMBLY_VERSION}
    $miniconda2_dir/conda install -y -c bioconda bax2bam=${BAX2BAM_VERSION}
    source $miniconda2_dir/deactivate
    rm Miniconda2-${MINICONDA2_VERSION}-Linux-x86_64.sh 
    note_installed $miniconda2_dir
fi

# ------------- NANOPOLISH --------------------
nanopolish_dir="$build_dir/nanopolish"
if [ -z $(check_installed $nanopolish_dir) ]; then
    cd $build_dir
    echo "Download nanopolish-v${NANOPOLISH_VERSION}"
    git clone --recursive $NANOPOLISH_DOWNLOAD_URL
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
    note_installed $nanopolish_dir
fi

# --------------- parallel ------------------
parallel_dir="$build_dir/parallel-${PARALLEL_VERSION}/bin"
if [ -z $(check_installed $parallel_dir) ]; then
    cd $build_dir
    echo "Download parallel"
    download_and_extract $PARALLEL_DOWNLOAD_URL "parallel_v${PARALLEL_VERSION}.tar.bz2"
    cd parallel-${PARALLEL_VERSION}
    ./configure --prefix="$build_dir/parallel-${PARALLEL_VERSION}"
    make -j $MAKE_JOBS
    make install
    parallel_dir="$build_dir/parallel-${PARALLEL_VERSION}/bin"
    note_installed $parallel_dir
fi

# --------------- EMBOSS ------------------
emboss_dir="$build_dir/EMBOSS-${EMBOSS_VERSION}/emboss"
if [ -z $(check_installed $emboss_dir) ]; then
    cd $build_dir
    echo "Download EMBOSS"
    download_and_extract $EMBOSS_DOWNLOAD_URL "emboss_v${EMBOSS_VERSION}.tar.gz"
    cd EMBOSS-${EMBOSS_VERSION}
    ./configure --without-x
    make -j $MAKE_JOBS
    note_installed $emboss_dir
fi

# --------------- ERPIN ------------------
erpin_dir="$build_dir/erpin${ERPIN_VERSION}.serv/bin"
if [ -z $(check_installed $erpin_dir) ]; then
    cd $build_dir
    echo "Download ERPIN"
    download_and_extract $ERPIN_DOWNLOAD_URL "erpin_v${ERPIN_VERSION}.tar.gz"
    cd erpin${ERPIN_VERSION}.serv
    make -j $MAKE_JOBS
    note_installed $erpin_dir
fi

# --------------- tbl2asn ------------------
tbl2asn_dir="$build_dir/tbl2asn_dir"
if [ -z $(check_installed $tbl2asn_dir) ]; then
    cd $build_dir
    echo "Download tbl2asn"
    download $TBL2ASN_DOWNLOAD_URL "tbl2asn.gz"
    mkdir tbl2asn_dir
    gunzip tbl2asn.gz
    chmod 755 tbl2asn
    mv tbl2asn ./tbl2asn_dir/
    note_installed $tbl2asn_dir
fi

# --------------- PirObject ----------------
pirobject_dir="$build_dir/PirObject-${PIROBJECT_VERSION}"
if [ -z $(check_installed $pirobject_dir) ]; then
    cd $build_dir
    echo "Download PirObject"
    download_and_extract $PIROBJECT_DOWNLOAD_URL "pirobject_v${PIROBJECT_VERSION}.tar.gz"
    cd PirObject-${PIROBJECT_VERSION}
    ln -s ./lib/PirObject.pm .
    note_installed $pirobject_dir
fi

# --------------- PirModels ------------------
pirmodels_dir="$pirobject_dir/PirModels"
if [ -z $(check_installed $pirmodels_dir) ]; then
    cd $build_dir
    echo "Download PirModels"
    git clone $PIRMODELS_DOWNLOAD_URL
    cd PirModels
    git checkout -f -q $PIRMODELS_GITHUB_COMMIT_VERSION
    cd ..
    cp -r PirModels $pirobject_dir
    note_installed $pirmodels_dir
fi


# --------------- Flip ------------------
flip_dir="$build_dir/Flip"
if [ -z $(check_installed $flip_dir) ]; then
    cd $build_dir
    echo "Download Flip"
    git clone $FLIP_DOWNLOAD_URL
    cd Flip
    git checkout -f -q $FLIP_GITHUB_COMMIT_VERSION
    cd src
    make -j $MAKE_JOBS
    cp flip ./../
    note_installed $flip_dir
fi

# --------------- Umac ------------------
umac_dir="$build_dir/Umac"
if [ -z $(check_installed $umac_dir) ]; then
    cd $build_dir
    echo "Download Umac"
    git clone $UMAC_DOWNLOAD_URL
    cd Umac
    git checkout -f -q $UMAC_GITHUB_COMMIT_VERSION
    note_installed $umac_dir
fi


# --------------- HMMsearchWC ------------------
hmmsearchwc_dir="$build_dir/HMMsearchWC"
if [ -z $(check_installed $hmmsearchwc_dir) ]; then
    cd $build_dir
    echo "Download HMMsearchWC"
    git clone $HMMSEARCHWC_DOWNLOAD_URL
    cd HMMsearchWC
    git checkout -f -q $HMMSEARCHWC_GITHUB_COMMIT_VERSION
    note_installed $hmmsearchwc_dir
fi


# --------------- RNAfinder ------------------
rnafinder_dir="$build_dir/RNAfinder"
if [ -z $(check_installed $rnafinder_dir) ]; then
    cd $build_dir
    echo "Download RNAfinder"
    git clone $RNAFINDER_DOWNLOAD_URL
    cd RNAfinder
    git checkout -f -q $RNAFINDER_GITHUB_COMMIT_VERSION
    note_installed $rnafinder_dir
fi


# --------------- Mf2sqn ------------------
mf2sqn_dir="$build_dir/Mf2sqn"
if [ -z $(check_installed $mf2sqn_dir) ]; then
    cd $build_dir
    echo "Download Mf2sqn"
    git clone $MF2SQN_DOWNLOAD_URL
    cd Mf2sqn
    git checkout -f -q $MF2SQN_GITHUB_COMMIT_VERSION
    cp qualifs.pl $build_dir/cpanm/perlmods/lib/perl5
    note_installed $mf2sqn_dir
fi

# --------------- grab-fasta ------------------
grab_fasta_dir="$build_dir/grab-fasta"
if [ -z $(check_installed $grab_fasta_dir) ]; then
    cd $build_dir
    echo "Download grab-fasta"
    git clone $GRAB_FASTA_DOWNLOAD_URL
    cd grab-fasta
    git checkout -f -q $GRAB_FASTA_GITHUB_COMMIT_VERSION
    note_installed $grab_fasta_dir
fi

# --------------- MFannot_data ------------------
mfannot_data_dir="$build_dir/MFannot_data"
if [ -z $(check_installed $mfannot_data_dir) ]; then
    cd $build_dir
    echo "Download MFannot_data"
    git clone $MFANNOT_DATA_DOWNLOAD_URL
    cd MFannot_data
    git checkout -f -q $MFANNOT_DATA_GITHUB_COMMIT_VERSION
    note_installed $mfannot_data_dir
fi


# --------------- MFannot ------------------
mfannot_dir="$build_dir/MFannot"
if [ -z $(check_installed $mfannot_dir) ]; then
    cd $build_dir
    echo "Download MFannot"
    git clone $MFANNOT_DOWNLOAD_URL
    cd MFannot
    git checkout -f -q $MFANNOT_GITHUB_COMMIT_VERSION
    note_installed $mfannot_dir
fi


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
# echo "export guppy_dir=${guppy_dir}" >> env.sh
# echo "export quast_dir=${quast_dir}" >> env.sh
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
