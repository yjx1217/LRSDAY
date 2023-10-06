#!/bin/bash
# last update: 2022/12/21
set -e -o pipefail

LRSDAY_HOME=$(pwd)
BUILD="build"
lite_installation="no";
mainland_china_installation="no";
#########################                                                                                                       

timestamp () {
  date +"%F %T"
}

clean () {
    dir=$1
    if [ -d $dir ] 
    then
	echo "remove previously failed installation in $BUILD/$dir"
	rm -rf $dir
    fi
}

download () {
  url=$1
  download_location=$2
  echo "Downloading $url to $download_location"
  wget -c --no-check-certificate $url -O $download_location
}

clone () {
  url=$1
  dir=$(basename $url)
  echo "run clone for \"git clone $url\""
  git clone $url --depth 1
  cd $dir
  git fetch --unshallow
}

tidy_version () { 
    echo "$1" | awk -F. '{ printf("%d%03d%03d%03d\n", $1,$2,$3,$4); }';
}

check_installed () {
    if [ -e "$1/installed" ]; then
        echo "installed"
    else
        echo ""
    fi
}

note_installed () {
    touch "$1/installed"
}

echo ""
echo ""
echo "##################################################################"
echo "###                                                            ###"
echo "###                  Welcome to LRSDAY                         ###"
echo "###                                                            ###"
echo "##################################################################"
echo ""
echo ""
echo "[$(timestamp)] Installation starts ..."
echo ""


while getopts ":hlc" opt
do
    case "${opt}" in
        h)
            echo "Usage:"
            echo "bash install_dependencies.sh"
            echo "When installing within mainland China, please run this script with the '-c' option >"
            echo "bash install_dependencies.sh -c";;
        l)
	    echo "Detected the '-l' option >"
	    echo "Set installation mode as 'lite' by skipping the installation of a few non-essential dependencies"
	    lite_installation="yes";;
        c)
	    echo "Detected the '-c' option >"
	    echo "Set installation location as 'mainland_china'"
	    mainland_china_installation="yes";;
    esac
done
echo "";

# echo "mainland_china_installation=$mainland_china_installation"

if [ -z "$MAKE_JOBS" ]
then
    echo "[$(timestamp)] Defaulting to 2 concurrent jobs when executing make. Override with MAKE_JOBS=<NUM>"
    MAKE_JOBS=2
fi

if [ ! -z "$INSTALL_DEPS" ]; then
    echo "Installing LRSDAY build dependencies for Debian/Ubuntu."
    echo "sudo privileges are required and you will be prompted to enter your password"
    sudo apt-get update
    xargs -a debiandeps sudo apt-get install -y
fi

MINICONDA2_VERSION="py27_4.8.3" # released on 2020.06.06
if [[ "$mainland_china_installation" == "no" ]]
then
    MINICONDA2_DOWNLOAD_URL="https://repo.anaconda.com/miniconda/Miniconda2-${MINICONDA2_VERSION}-Linux-x86_64.sh"
else
    MINICONDA2_DOWNLOAD_URL="https://mirrors.tuna.tsinghua.edu.cn/anaconda/miniconda/Miniconda2-${MINICONDA2_VERSION}-Linux-x86_64.sh"
fi

MINICONDA3_VERSION="py37_4.9.2" # released on 2020.11.23                                                                                                                
if [[ "$mainland_china_installation" == "yes" ]]
then
    MINICONDA3_DOWNLOAD_URL="https://repo.anaconda.com/miniconda/Miniconda3-${MINICONDA3_VERSION}-Linux-x86_64.sh"
else
    MINICONDA3_DOWNLOAD_URL="https://mirrors.tuna.tsinghua.edu.cn/anaconda/miniconda/Miniconda3-${MINICONDA3_VERSION}-Linux-x86_64.sh"
fi

# for reads preparation and preprocessing
GUPPY_VERSION="6.2.1" # released on 2022.06.07
GUPPY_GPU_DOWNLOAD_URL="https://mirror.oxfordnanoportal.com/software/analysis/ont-guppy_${GUPPY_VERSION}_linux64.tar.gz"
GUPPY_CPU_DOWNLOAD_URL="https://mirror.oxfordnanoportal.com/software/analysis/ont-guppy-cpu_${GUPPY_VERSION}_linux64.tar.gz"

PORECHOP_VERSION="0.2.4" # 
PORECHOP_GITHUB_COMMIT_VERSION="109e437" # committed on 2018.10.19
PORECHOP_DOWNLOAD_URL="https://github.com/rrwick/Porechop.git"

FILTLONG_VERSION="0.2.1" #
FILTLONG_GITHUB_COMMIT_VERSION="c56f02b" # committed on 2021.07.01
FILTLONG_DOWNLOAD_URL="https://github.com/rrwick/Filtlong.git"

NANOPLOT_VERSION="1.40.0" # released on 2022.04.08

# for genome assembly
MINIMAP2_VERSION="2.22" # released on 2021.08.07
MINIMAP2_DOWNLOAD_URL="https://github.com/lh3/minimap2/releases/download/v${MINIMAP2_VERSION}/minimap2-${MINIMAP2_VERSION}_x64-linux.tar.bz2"

CANU_VERSION="2.2" # released on 2021.08.27
CANU_DOWNLOAD_URL="https://github.com/marbl/canu/releases/download/v${CANU_VERSION}/canu-${CANU_VERSION}.Linux-amd64.tar.xz"

# WHATSHAP_VERSION="1.1"

FLYE_VERSION="2.9.1" # released on 2021.08.22
FLYE_DOWNLOAD_URL="https://github.com/fenderglass/Flye/archive/${FLYE_VERSION}.tar.gz"

WTDBG2_VERSION="2.5" # 
WTDBG2_GITHUB_COMMIT_VERSION="b77c565" # committed on 2020.12.11
WTDBG2_DOWNLOAD_URL="https://github.com/ruanjue/wtdbg2.git"


SMARTDENOVO_VERSION="" # 
SMARTDENOVO_GITHUB_COMMIT_VERSION="8488de9" # committed on 2021.02.24
SMARTDENOVO_DOWNLOAD_URL="https://github.com/ruanjue/smartdenovo"

RAVEN_VERSION="1.8.1" # released on 2022.08.04
RAVEN_GITHUB_COMMIT_VERSION="0ea10b3"
#RAVEN_DOWNLOAD_URL="https://github.com/lbcb-sci/raven/releases/download/${RAVEN_VERSION}/raven-v${RAVEN_VERSION}.tar.gz"

SHASTA_VERSION="0.11.1" # 
SHASTA_GITHUB_COMMIT_VERSION="b739b1e"
SHASTA_DOWNLOAD_URL="https://github.com/paoloshasta/shasta/releases/download/${SHASTA_VERSION}/shasta-Linux-${SHASTA_VERSION}"

# for assembly polishing

PB_ASSEMBLY_VERSION="0.0.8" #
BAX2BAM_VERSION="0.0.9" #
PBMM2_VERSION="1.9.0" #

NANOPOLISH_VERSION="0.14.0" # released on 2021.04.06

PARALLEL_VERSION="20180722" # released on 2018.07.22
PARALLEL_DOWNLOAD_URL="http://ftp.gnu.org/gnu/parallel/parallel-${PARALLEL_VERSION}.tar.bz2"

RACON_VERSION="1.5.0" # released on 2021.03.26
#RACON_GITHUB_COMMIT_VERSION="9119181"
#RACON_DOWNLOAD_URL="https://github.com/lbcb-sci/racon/releases/download/${RACON_VERSION}/racon-v${RACON_VERSION}.tar.gz"
#RACON_DOWNLOAD_URL="https://github.com/lbcb-sci/racon.git"

MEDAKA_VERSION="1.7.2" # released on 2022.09.21
MEDAKA_DOWNLOAD_URL="https://github.com/nanoporetech/medaka/archive/v${MEDAKA_VERSION}.tar.gz"

# MARGINPOLISH_VERSION="1.3.0" # released on 2020.03.04
# MARGINPOLISH_GITHUB_COMMIT_VERSION="5492204" # commited on 2020.03.25 
# # MARGINPOLISH_DOWNLOAD_URL="https://github.com/UCSC-nanopore-cgl/marginPolish.git"
# MARGINPOLISH_DOWNLOAD_URL="https://github.com/UCSC-nanopore-cgl/MarginPolish/archive/refs/tags/v${MARGINPOLISH_VERSION}.tar.gz"

# HELEN_VERSION="0.0.1" # released on 2020.12.10
# HELEN_GITHUB_COMMIT_VERSION="a075e9f" # commited on 2020.10.28 
# HELEN_DOWNLOAD_URL="https://github.com/kishwarshafin/helen.git"

# HOMOPOLISH_VERSION="0.2.1" # released on 2020.12.10
# HOMOPOLISH_DOWNLOAD_URL="https://github.com/ythuang0522/homopolish"

# MARGIN_VERSION="2.3.1" # released on 2022.03.25
# MARGIN_DOWNLOAD_URL="https://github.com/UCSC-nanopore-cgl/margin"
# MARGIN_GITHUB_COMMIT_VERSION="9bf845f" # commited on 2022.03.25 

# PEPPER_VERSION="0.4" # released on 2020.12.10
# PEPPER_GITHUB_COMMIT_VERSION="734d226" # commited on 2021.03.08 
# PEPPER_DOWNLOAD_URL="https://github.com/kishwarshafin/pepper.git"

# for assembly scaffolding
RAGOUT_VERSION="2.3" # released on 2020.03.18
RAGOUT_GITHUB_COMMIT_VERSION="7b92fe7" # commited on 2020.12.09"
RAGOUT_DOWNLOAD_URL="https://github.com/fenderglass/Ragout.git"

RAGTAG_VERSION="2.1.0" # released on 2021.11.01

# QUAST_VERSION="5.0.2" # one of its dependency needs "csh" to be pre-installed

HDF_VERSION="1.10.6" # 
HDF_VERSION_prefix=${HDF_VERSION%.*}
HDF_DOWNLOAD_URL="https://support.hdfgroup.org/ftp/HDF5/releases/hdf5-${HDF_VERSION_prefix}/hdf5-${HDF_VERSION}/src/hdf5-${HDF_VERSION}.tar.gz"

# SONLIB_VERSION="" # 
# SONLIB_GITHUB_COMMIT_VERSION="1afbd97" # committed on 2017.08.09
# SONLIB_DOWNLOAD_URL="https://github.com/benedictpaten/sonLib.git"

# HAL_VERSION="" # not available, so we use the github comit hash below for version control
# HAL_GITHUB_COMMIT_VERSION="a2ad656" # committed on 2017.09.09
# HAL_DOWNLOAD_URL="https://github.com/glennhickey/hal.git"

MUMMER4_VERSION="4.0.0rc1" # released on 2020.10.03
MUMMER4_DOWNLOAD_URL="https://github.com/gmarcais/mummer/releases/download/v${MUMMER4_VERSION}/mummer-${MUMMER4_VERSION}.tar.gz"

GNUPLOT_VERSION="4.6.6" # released on 2015.02.18
GNUPLOT_DOWNLOAD_URL="https://sourceforge.net/projects/gnuplot/files/gnuplot/${GNUPLOT_VERSION}/gnuplot-${GNUPLOT_VERSION}.tar.gz"

BEDTOOLS_VERSION="2.30.0" # released on 2021.01.24
BEDTOOLS_DOWNLOAD_URL="https://github.com/arq5x/bedtools2/releases/download/v${BEDTOOLS_VERSION}/bedtools-${BEDTOOLS_VERSION}.tar.gz"

SPADES_VERSION="3.14.1" # released on 2020.05.02
SPADES_DOWNLOAD_URL="http://cab.spbu.ru/files/release${SPADES_VERSION}/SPAdes-${SPADES_VERSION}-Linux.tar.gz"

PRODIGAL_VERSION="2.6.3" # released on 2016.02.12
PRODIGAL_DOWNLOAD_URL="https://github.com/hyattpd/Prodigal/archive/v${PRODIGAL_VERSION}.tar.gz"

CAP3_VERSION="" # see http://seq.cs.iastate.edu/cap3.html
CAP3_DOWNLOAD_URL="https://github.com/yjx1217/CAP3.git" # "http://seq.cs.iastate.edu/CAP3/cap3.linux.x86_64.tar"

BWA_VERSION="0.7.17" # released on 2017.10.23
BWA_DOWNLOAD_URL="http://downloads.sourceforge.net/project/bio-bwa/bwa-${BWA_VERSION}.tar.bz2"

SAMTOOLS_VERSION="1.16.1" # released on 2022.09.02
HTSLIB_VERSION="1.16" # released on 2022.09.02
SAMTOOLS_DOWNLOAD_URL="https://github.com/samtools/samtools/releases/download/${SAMTOOLS_VERSION}/samtools-${SAMTOOLS_VERSION}.tar.bz2"

CIRCLATOR_VERSION="1.5.5" # released on 2020.10.29
CIRCLATOR_DOWNLOAD_URL="https://github.com/sanger-pathogens/circlator/archive/v${CIRCLATOR_VERSION}.tar.gz"

TRIMMOMATIC_VERSION="0.38" # 
TRIMMOMATIC_DOWNLOAD_URL="http://www.usadellab.org/cms/uploads/supplementary/Trimmomatic/Trimmomatic-${TRIMMOMATIC_VERSION}.zip"

GATK3_VERSION="3.6-6" #
GATK3_DOWNLOAD_URL="https://github.com/yjx1217/GATK3_Archive/raw/main/gatk3.jar.gz"

PICARD_VERSION="2.23.4" # released on 2020.09.03
PICARD_DOWNLOAD_URL="https://github.com/broadinstitute/picard/releases/download/${PICARD_VERSION}/picard.jar"

# HAPOG_VERSION="1.2" # released on 2021.10.01
 
PILON_VERSION="1.24" # released on 2021.01.29
PILON_DOWNLOAD_URL="https://github.com/broadinstitute/pilon/releases/download/v${PILON_VERSION}/pilon-${PILON_VERSION}.jar"

EXONERATE_VERSION="2.2.0" # 
EXONERATE_DOWNLOAD_URL="http://ftp.ebi.ac.uk/pub/software/vertebrategenomics/exonerate/exonerate-${EXONERATE_VERSION}-x86_64.tar.gz"

BLAST_VERSION="2.2.31" # 
RMBLAST_VERSION="2.2.28" # 
BLAST_DOWNLOAD_URL="http://ftp.ncbi.nlm.nih.gov/blast/executables/blast+/${BLAST_VERSION}/ncbi-blast-${BLAST_VERSION}+-x64-linux.tar.gz"
RMBLAST_DOWNLOAD_URL="http://ftp.ncbi.nlm.nih.gov/blast/executables/rmblast/${RMBLAST_VERSION}/ncbi-rmblastn-${RMBLAST_VERSION}-x64-linux.tar.gz"

SNAP_VERSION="" # 
SNAP_GITHUB_COMMIT_VERSION="a89d68e" # committed on 2017.05.18
# SNAP_GITHUB_COMMIT_VERSION="bbf2050" # commited on 2022.04.12 
SNAP_DOWNLOAD_URL="https://github.com/KorfLab/SNAP.git"

RAPSEARCH_VERSION="2.24" #
RAPSEARCH_DOWNLOAD_URL="https://sourceforge.net/projects/rapsearch2/files/RAPSearch${RAPSEARCH_VERSION}_64bits.tar.gz"

TRNASCAN_VERSION="1.3.1" #
TRNASCAN_DOWNLOAD_URL="http://eddylab.org/software/tRNAscan-SE/tRNAscan-SE.tar.gz"

SNOSCAN_VERSION="0.9.1" #
SNOSCAN_DOWNLOAD_URL="http://eddylab.org/software/snoscan/snoscan.tar.gz"

REPEATMASKER_VERSION="4.1.4" #
# REPEATMASKER_VERSION="4.1.0" #
REPEATMASKER_DOWNLOAD_URL="http://www.repeatmasker.org/RepeatMasker/RepeatMasker-${REPEATMASKER_VERSION}.tar.gz"

#REPBASE_VERSION="20181026"
REPBASE_VERSION="20170127"
REPBASE_DOWNLOAD_URL="https://github.com/yjx1217/RMRB.git"

TRF_VERSION="409" #
TRF_DOWNLOAD_URL="https://github.com/Benson-Genomics-Lab/TRF/releases/download/v4.09.1/trf409.linux64"
REANNOTATE_VERSION="17.03.2015-LongQueryName"
REANNOTATE_GITHUB_COMMIT_VERSION="7fdb339"
REANNOTATE_DOWNLOAD_URL="https://github.com/yjx1217/REannotate_LongQueryName"

CLUSTALW_VERSION="2.1" #
CLUSTALW_DOWNLOAD_URL="http://www.clustal.org/download/${CLUSTALW_VERSION}/clustalw-${CLUSTALW_VERSION}.tar.gz"

MUSCLE_VERSION="3.8.31" #
MUSCLE_DOWNLOAD_URL="http://www.drive5.com/muscle/downloads${MUSCLE_VERSION}/muscle${MUSCLE_VERSION}_i86linux64.tar.gz"

HMMER_VERSION="3.2.1" # released on 2018.06.13
HMMER_DOWNLOAD_URL="http://eddylab.org/software/hmmer/hmmer-${HMMER_VERSION}.tar.gz"

BAMTOOLS_VERSION="2.4.2" # released on 2017.11.02
BAMTOOLS_DOWNLOAD_URL="https://github.com/pezmaster31/bamtools/archive/v${BAMTOOLS_VERSION}.tar.gz"

AUGUSTUS_VERSION="3.2.3" # 
#AUGUSTUS_GITHUB_COMMIT_VERSION="79960c5"
AUGUSTUS_DOWNLOAD_URL="http://bioinf.uni-greifswald.de/augustus/binaries/old/augustus-${AUGUSTUS_VERSION}.tar.gz"
#AUGUSTUS_DOWNLOAD_URL="https://github.com/Gaius-Augustus/Augustus.git"

EVM_VERSION="1.1.1" # released on 2015.07.03
EVM_DOWNLOAD_URL="https://github.com/EVidenceModeler/EVidenceModeler/archive/v${EVM_VERSION}.tar.gz"

PROTEINORTHO_VERSION="6.0.35" # released on 2022.09.15
DIAMOND_VERSION="2.0.6"

#MAKER_VERSION="3.01.03" #
MAKER_VERSION="3.00.0-beta" #
MAKER_DOWNLOAD_URL="https://github.com/yjx1217/Maker3beta_archive.git"

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

RNAFINDER_GITHUB_COMMIT_VERSION="ee5b7de" # committed on 2019.11.26
RNAFINDER_DOWNLOAD_URL="https://github.com/BFL-lab/RNAfinder.git"

MF2SQN_GITHUB_COMMIT_VERSION="6faf9f4" # committed on 2016.12.07
MF2SQN_DOWNLOAD_URL="https://github.com/BFL-lab/Mf2sqn.git"

GRAB_FASTA_GITHUB_COMMIT_VERSION="accd32d" # committed on 2017.02.14
GRAB_FASTA_DOWNLOAD_URL="https://github.com/BFL-lab/grab-fasta.git"

# for MFannot
MFANNOT_DATA_GITHUB_COMMIT_VERSION="b039ac5" # committed on 2016.12.07
MFANNOT_VERSION="1.35" #
MFANNOT_GITHUB_COMMIT_VERSION="6472b97" # committed on 2018.10.31 # "104e1c5" # committed on 2021.02.10
MFANNOT_DATA_DOWNLOAD_URL="https://github.com/BFL-lab/MFannot_data.git"
MFANNOT_DOWNLOAD_URL="https://github.com/BFL-lab/MFannot.git"

# # GFF3
# GFF3TOOLKIT_VERSION="2.1.0" # released on 2021.12.09
# GFF3TOOLKIT_DOWNLOAD_URL="https://github.com/NAL-i5K/GFF3toolkit"

# UCSC Utilities
BLAT_DOWNLOAD_URL="http://hgdownload.soe.ucsc.edu/admin/exe/linux.x86_64/blat/blat"
FASPLIT_DOWNLOAD_URL="http://hgdownload.soe.ucsc.edu/admin/exe/linux.x86_64/faSplit"
PSLSORT_DOWNLOAD_URL="http://hgdownload.soe.ucsc.edu/admin/exe/linux.x86_64/pslSort"
PSLSCORE_DOWNLOAD_URL="http://hgdownload.soe.ucsc.edu/admin/exe/linux.x86_64/pslScore"
PSLCDNAFILTER_DOWNLOAD_URL="http://hgdownload.soe.ucsc.edu/admin/exe/linux.x86_64/pslCDnaFilter"


# Create the $BUILD directory for dependency installation

if [ -d $BUILD ]
then
    echo ""
    echo "[$(timestamp)] Detected previously generated $BUILD directory."
else
    echo "[$(timestamp)] Create the new $BUILD directory."
    mkdir $BUILD
    echo ""
fi

cd $BUILD
build_dir=$(pwd)

echo ""
echo "[$(timestamp)] Download and install all the dependencies"
echo ""

# Downloading all the dependencies

# ---------- set Perl & Python environment variables -------------
PYTHONPATH="$build_dir"
PERL5LIB="$build_dir:$PERL5LIB"
PERL5LIB="$build_dir/cpanm/perlmods/lib/perl5:$PERL5LIB"

echo ""
echo "[$(timestamp)] Installing Perl modules ..."
cpanm_dir=$build_dir/cpanm
if [ -z $(check_installed $cpanm_dir) ]; then
    clean "$build_dir/cpanm"
    mkdir -p  $cpanm_dir
    cd $cpanm_dir
    # wget -c --no-check-certificate -O - https://cpanmin.us/ > cpanm
    # work around for the unstable downloading issue
    cp $LRSDAY_HOME/misc/cpanm .

    chmod +x cpanm
    mkdir -p perlmods
    $cpanm_dir/cpanm -l $cpanm_dir/perlmods --notest --skip-installed Test::More@1.302086
    $cpanm_dir/cpanm -l $cpanm_dir/perlmods --notest --skip-installed Text::Soundex@3.05
    $cpanm_dir/cpanm -l $cpanm_dir/perlmods --notest --skip-installed Env@1.04
    $cpanm_dir/cpanm -l $cpanm_dir/perlmods --notest --skip-installed File::Which@1.21
    $cpanm_dir/cpanm -l $cpanm_dir/perlmods --notest --skip-installed Term::ReadKey@2.37
    $cpanm_dir/cpanm -l $cpanm_dir/perlmods --notest --skip-installed Carp@1.38
    $cpanm_dir/cpanm -l $cpanm_dir/perlmods --notest --skip-installed Perl::Unsafe::Signals@0.03
    $cpanm_dir/cpanm -l $cpanm_dir/perlmods --notest --skip-installed Carp::Clan@6.06 
    $cpanm_dir/cpanm -l $cpanm_dir/perlmods --notest --skip-installed Bit::Vector@7.4
    $cpanm_dir/cpanm -l $cpanm_dir/perlmods --notest --skip-installed Inline@0.80
    $cpanm_dir/cpanm -l $cpanm_dir/perlmods --notest --skip-installed Inline::C@0.78
    $cpanm_dir/cpanm -l $cpanm_dir/perlmods --notest --skip-installed List::MoreUtils@0.419
    $cpanm_dir/cpanm -l $cpanm_dir/perlmods --notest --skip-installed Acme::Damn@0.08
    $cpanm_dir/cpanm -l $cpanm_dir/perlmods --notest --skip-installed Sys::SigAction@0.23
    $cpanm_dir/cpanm -l $cpanm_dir/perlmods --notest --skip-installed forks@0.36
    $cpanm_dir/cpanm -l $cpanm_dir/perlmods --notest --skip-installed forks::shared@0.36
    $cpanm_dir/cpanm -l $cpanm_dir/perlmods --notest --skip-installed Want@0.29
    $cpanm_dir/cpanm -l $cpanm_dir/perlmods --notest --skip-installed IO::All@0.86
    $cpanm_dir/cpanm -l $cpanm_dir/perlmods --notest --skip-installed IO::Prompt@0.997004
    $cpanm_dir/cpanm -l $cpanm_dir/perlmods --notest --skip-installed DBI@1.636
    $cpanm_dir/cpanm -l $cpanm_dir/perlmods --notest --skip-installed Proc::ProcessTable@0.53
    $cpanm_dir/cpanm -l $cpanm_dir/perlmods --notest --skip-installed threads@2.16
    $cpanm_dir/cpanm -l $cpanm_dir/perlmods --notest --skip-installed PerlIO::gzip@0.20
    $cpanm_dir/cpanm -l $cpanm_dir/perlmods --notest --skip-installed ExtUtils::CBuilder@0.280224 
    $cpanm_dir/cpanm -l $cpanm_dir/perlmods --notest --skip-installed LWP::Simple@6.26
    $cpanm_dir/cpanm -l $cpanm_dir/perlmods --notest --skip-installed LWP::UserAgent
    $cpanm_dir/cpanm -l $cpanm_dir/perlmods --notest --skip-installed Bio::Perl@1.007001
    if [ ! -z "$USE_POSTGRES" ]; then
	# need $POSTGRES_HOME pre-defined
	$cpanm_dir/cpanm -l $cpanm_dir/perlmods --notest --skip-installed DBD::Pg
    else
	$cpanm_dir/cpanm -l $cpanm_dir/perlmods --notest --skip-installed DBD::SQLite@1.54
    fi
    note_installed $cpanm_dir
fi

# ------------- Miniconda2 --------------------
if [[ $lite_installation == "no" ]]
then
    echo ""
    echo "[$(timestamp)] Installing miniconda2 ..."
    miniconda2_dir="$build_dir/miniconda2/bin"
    if [ -z $(check_installed $miniconda2_dir) ]; then
	cd $build_dir
	clean "$build_dir/miniconda2"
	download $MINICONDA2_DOWNLOAD_URL "Miniconda2-${MINICONDA2_VERSION}-Linux-x86_64.sh"
	bash Miniconda2-${MINICONDA2_VERSION}-Linux-x86_64.sh -b -p $build_dir/miniconda2
	if [[ "$mainland_china_installation" == "yes" ]]
	then
	    $miniconda2_dir/conda config --add channels https://mirrors.bfsu.edu.cn/anaconda/pkgs/main
            $miniconda2_dir/conda config --add channels https://mirrors.bfsu.edu.cn/anaconda/pkgs/free
            $miniconda2_dir/conda config --add channels https://mirrors.bfsu.edu.cn/anaconda/pkgs/pro
            $miniconda2_dir/conda config --add channels https://mirrors.bfsu.edu.cn/anaconda/pkgs/msys2
            $miniconda2_dir/conda config --add channels https://mirrors.bfsu.edu.cn/anaconda/cloud/bioconda
            $miniconda2_dir/conda config --add channels https://mirrors.bfsu.edu.cn/anaconda/cloud/conda-forge
	    
            $miniconda2_dir/conda config --add channels https://anaconda.mirrors.sjtug.sjtu.edu.cn/pkgs/main
            $miniconda2_dir/conda config --add channels https://anaconda.mirrors.sjtug.sjtu.edu.cn/pkgs/free
            $miniconda2_dir/conda config --add channels https://anaconda.mirrors.sjtug.sjtu.edu.cn/pkgs/pro
            $miniconda2_dir/conda config --add channels https://anaconda.mirrors.sjtug.sjtu.edu.cn/pkgs/msys2
            $miniconda2_dir/conda config --add channels https://anaconda.mirrors.sjtug.sjtu.edu.cn/cloud/bioconda
            $miniconda2_dir/conda config --add channels https://anaconda.mirrors.sjtug.sjtu.edu.cn/cloud/conda-forge
	    
	    # $miniconda2_dir/conda config --add channels https://mirrors.tuna.tsinghua.edu.cn/anaconda/pkgs/free/
	    # $miniconda2_dir/conda config --add channels https://mirrors.tuna.tsinghua.edu.cn/anaconda/pkgs/main/
	    # $miniconda2_dir/conda config --add channels https://mirrors.tuna.tsinghua.edu.cn/anaconda/cloud/pytorch/
	    # $miniconda2_dir/conda config --add channels https://mirrors.tuna.tsinghua.edu.cn/anaconda/cloud/bioconda/
	    # $miniconda2_dir/conda config --add channels https://mirrors.tuna.tsinghua.edu.cn/anaconda/cloud/conda-forge/
	else 
	    $miniconda2_dir/conda config --add channels defaults
	    $miniconda2_dir/conda config --add channels bioconda
	    $miniconda2_dir/conda config --add channels conda-forge
	fi
	$miniconda2_dir/conda config --set show_channel_urls yes
	$miniconda2_dir/conda config --set ssl_verify no
	$miniconda2_dir/conda config --set channel_priority flexible
	$miniconda2_dir/pip install --timeout 1000 cython==0.29.14 numpy==1.13.1
	rm Miniconda2-${MINICONDA2_VERSION}-Linux-x86_64.sh 
	note_installed $miniconda2_dir
    fi
fi

# ------------- Miniconda3 --------------------
echo ""
echo "[$(timestamp)] Installing miniconda3 ..."
miniconda3_dir="$build_dir/miniconda3/bin"
if [ -z $(check_installed $miniconda3_dir) ]; then
    cd $build_dir
    clean "$build_dir/miniconda3"
    download $MINICONDA3_DOWNLOAD_URL "Miniconda3-${MINICONDA3_VERSION}-Linux-x86_64.sh"
    bash Miniconda3-${MINICONDA3_VERSION}-Linux-x86_64.sh -b -p $build_dir/miniconda3
    if [[ "$mainland_china_installation" == "yes" ]]
    then
	$miniconda3_dir/conda config --add channels https://mirrors.bfsu.edu.cn/anaconda/pkgs/main
        $miniconda3_dir/conda config --add channels https://mirrors.bfsu.edu.cn/anaconda/pkgs/free
        $miniconda3_dir/conda config --add channels https://mirrors.bfsu.edu.cn/anaconda/pkgs/pro
        $miniconda3_dir/conda config --add channels https://mirrors.bfsu.edu.cn/anaconda/pkgs/msys2
        $miniconda3_dir/conda config --add channels https://mirrors.bfsu.edu.cn/anaconda/cloud/bioconda
        $miniconda3_dir/conda config --add channels https://mirrors.bfsu.edu.cn/anaconda/cloud/conda-forge

        $miniconda3_dir/conda config --add channels https://anaconda.mirrors.sjtug.sjtu.edu.cn/pkgs/main
        $miniconda3_dir/conda config --add channels https://anaconda.mirrors.sjtug.sjtu.edu.cn/pkgs/free
        $miniconda3_dir/conda config --add channels https://anaconda.mirrors.sjtug.sjtu.edu.cn/pkgs/pro
        $miniconda3_dir/conda config --add channels https://anaconda.mirrors.sjtug.sjtu.edu.cn/pkgs/msys2
        $miniconda3_dir/conda config --add channels https://anaconda.mirrors.sjtug.sjtu.edu.cn/cloud/bioconda
        $miniconda3_dir/conda config --add channels https://anaconda.mirrors.sjtug.sjtu.edu.cn/cloud/conda-forge

	# $miniconda3_dir/conda config --add channels https://mirrors.tuna.tsinghua.edu.cn/anaconda/pkgs/free/
	# $miniconda3_dir/conda config --add channels https://mirrors.tuna.tsinghua.edu.cn/anaconda/pkgs/main/
	# $miniconda3_dir/conda config --add channels https://mirrors.tuna.tsinghua.edu.cn/anaconda/cloud/pytorch/
	# $miniconda3_dir/conda config --add channels https://mirrors.tuna.tsinghua.edu.cn/anaconda/cloud/bioconda/
	# $miniconda3_dir/conda config --add channels https://mirrors.tuna.tsinghua.edu.cn/anaconda/cloud/conda-forge/
    else 
	$miniconda3_dir/conda config --add channels defaults
	$miniconda3_dir/conda config --add channels bioconda
	$miniconda3_dir/conda config --add channels conda-forge
    fi
    $miniconda3_dir/conda config --set show_channel_urls yes
    $miniconda3_dir/conda config --set ssl_verify no
    $miniconda3_dir/conda config --set channel_priority flexible
    cd $build_dir
    rm Miniconda3-${MINICONDA3_VERSION}-Linux-x86_64.sh 
    note_installed $miniconda3_dir
fi


# ----------------- PB-ASSEMBLY ----------------------
if [[ $lite_installation == "no" ]]
then
    echo ""
    echo "[$(timestamp)] Installing pb-assembly ..."
    pbassembly_dir=$build_dir/pbassembly_conda_env/bin
    if [ -z $(check_installed $pbassembly_dir) ]; then
	cd $build_dir
	clean "$build_dir/pbassembly_conda_env"
	$miniconda3_dir/conda create -y -p $build_dir/pbassembly_conda_env 
	source $miniconda3_dir/activate $build_dir/pbassembly_conda_env
	$miniconda3_dir/conda install -y -c bioconda pb-assembly=${PB_ASSEMBLY_VERSION}
	source $miniconda3_dir/deactivate
	note_installed $pbassembly_dir
    fi
fi

# ----------------- BAX2BAM  ----------------------
if [[ $lite_installation == "no" ]]
then
    echo ""
    echo "[$(timestamp)] Installing bax2bam ..."
    bax2bam_dir="$build_dir/bax2bam_conda_env/bin"
    if [ -z $(check_installed $bax2bam_dir) ]; then
	cd $build_dir
	clean "$build_dir/bax2bam_conda_env"
	$miniconda2_dir/conda create -y -p $build_dir/bax2bam_conda_env 
	source $miniconda2_dir/activate $build_dir/bax2bam_conda_env
	#$miniconda2_dir/conda install -y -c bioconda bax2bam=${BAX2BAM_VERSION}
	$miniconda2_dir/conda install -y -c "bioconda/label/cf201901" bax2bam 
	source $miniconda2_dir/deactivate
	note_installed $bax2bam_dir
    fi
fi

# --------------- Porechop ------------------
echo ""
echo "[$(timestamp)] Installing Porechop ..."
porechop_dir="$build_dir/porechop_conda_env/bin"
if [ -z $(check_installed $porechop_dir) ]; then
    cd $build_dir
    clean "$build_dir/porechop_conda_env"
    echo "Download Porechop-v${PORECHOP_VERSION}"
    $miniconda3_dir/conda create -y -p $build_dir/porechop_conda_env python=3.7
    source $miniconda3_dir/activate $build_dir/porechop_conda_env
    $miniconda3_dir/conda install -y -c bioconda porechop=${PORECHOP_VERSION} 
    source $miniconda3_dir/deactivate
    note_installed $porechop_dir
fi

# --------------- Filtlong ------------------
echo ""
echo "[$(timestamp)] Installing Filtlong ..."
filtlong_dir="$build_dir/Filtlong/bin"
if [ -z $(check_installed $filtlong_dir) ]; then
    cd $build_dir
    clean "$build_dir/Filtlong"
    echo "Download Filtlong-v${FILTLONG_VERSION}"
    git clone $FILTLONG_DOWNLOAD_URL
    cd Filtlong
    git checkout -f -q $FILTLONG_GITHUB_COMMIT_VERSION
    make -j $MAKE_JOBS
    note_installed $filtlong_dir
fi

# --------------- minimap2 ------------------
echo ""
echo "[$(timestamp)] Installing minimap2 ..."
minimap2_dir="$build_dir/minimap2-${MINIMAP2_VERSION}_x64-linux"
if [ -z $(check_installed $minimap2_dir) ]; then
    cd $build_dir
    clean "$build_dir/minimap2-${MINIMAP2_VERSION}_x64-linux"
    echo "Download minimap2-v${MINIMAP2_VERSION}"
    download $MINIMAP2_DOWNLOAD_URL "minimap2-${MINIMAP2_VERSION}.tar.bz2"
    tar xvjf minimap2-${MINIMAP2_VERSION}.tar.bz2
    rm minimap2-${MINIMAP2_VERSION}.tar.bz2
    note_installed $minimap2_dir
fi
PATH=$minimap2_dir:${PATH}

# --------------- samtools -----------------
echo ""
echo "[$(timestamp)] Installing samtools ..."
samtools_dir="$build_dir/samtools-${SAMTOOLS_VERSION}"
htslib_dir="$samtools_dir/htslib-${HTSLIB_VERSION}"
tabix_dir="$samtools_dir/htslib-${HTSLIB_VERSION}"

if [ -z $(check_installed $samtools_dir) ]; then
    cd $build_dir
    clean "$build_dir/samtools-${SAMTOOLS_VERSION}"
    echo "Download samtools-v${SAMTOOLS_VERSION}"
    download $SAMTOOLS_DOWNLOAD_URL "samtools-${SAMTOOLS_VERSION}.tar.bz2"
    tar xvjf samtools-${SAMTOOLS_VERSION}.tar.bz2
    cd $samtools_dir
    C_INCLUDE_PATH=""
    ./configure --without-curses;
    make -j $MAKE_JOBS
    cd $htslib_dir
    ./configure
    make -j $MAKE_JOBS
    cd $build_dir
    rm samtools-${SAMTOOLS_VERSION}.tar.bz2
    note_installed $samtools_dir
fi
PATH="$samtools_dir:$htslib_dir:$tabix_dir:${PATH}"

# ------------- Canu -------------------
echo ""
echo "[$(timestamp)] Installing canu ..."
canu_dir="$build_dir/canu-${CANU_VERSION}/bin"
if [ -z $(check_installed $canu_dir) ]; then
    cd $build_dir
    clean "$build_dir/canu-${CANU_VERSION}"
    echo "Download Canu-v${CANU_VERSION}"
    download $CANU_DOWNLOAD_URL "canu-${CANU_VERSION}.tar.xz"
    tar xvf canu-${CANU_VERSION}.tar.xz
    cd $canu_dir
    ln -s $minimap2_dir/minimap2 .
    cd $build_dir
    rm canu-${CANU_VERSION}.tar.xz
    note_installed $canu_dir
fi

# # --------------- WhatsHap -----------------
# echo ""
# echo "[$(timestamp)] Installing WhatsHap ..."
# whatshap_dir="$build_dir/whatshap_conda_env/bin"
# if [ -z $(check_installed $whatshap_dir) ]; then
#     cd $build_dir
#     # clean "$build_dir/whatshap_conda_env"
#     echo "Download whatshap-v${WHATSHAP_VERSION}"
#     # $miniconda3_dir/conda create -y -p $build_dir/whatshap_conda_env python=3.9
#     source $miniconda3_dir/activate $build_dir/whatshap_conda_env
#     $miniconda3_dir/conda install -y -c bioconda whatshap=${WHATSHAP_VERSION} nomkl 
#     source $miniconda3_dir/deactivate
#     note_installed $whatshap_dir
# fi


# ------------- Flye -------------------
echo ""
echo "[$(timestamp)] Installing flye ..."
flye_dir="$build_dir/flye_conda_env/bin"
if [ -z $(check_installed $flye_dir) ]; then
    cd $build_dir
    clean "$build_dir/flye_conda_env"
    $miniconda3_dir/conda create -y -p $build_dir/flye_conda_env python=3.7
    source $miniconda3_dir/activate $build_dir/flye_conda_env
    $miniconda3_dir/conda install -y -c bioconda flye=${FLYE_VERSION} 
    source $miniconda3_dir/deactivate
    note_installed $flye_dir
fi

# --------------- wtdbg2 ------------------
echo ""
echo "[$(timestamp)] Installing wtdbg2 ..."
wtdbg2_dir="$build_dir/wtdbg2"
if [ -z $(check_installed $wtdbg2_dir) ]; then
    cd $build_dir
    clean "$build_dir/wtdbg2"
    echo "Download wtdbg2-v${WTDBG2_VERSION}"
    git clone $WTDBG2_DOWNLOAD_URL
    cd wtdbg2
    git checkout -f -q $WTDBG2_GITHUB_COMMIT_VERSION
    C_INCLUDE_PATH="" 
    make -j $MAKE_JOBS
    note_installed $wtdbg2_dir
fi

# --------------- smartdenovo ------------------
echo ""
echo "[$(timestamp)] Installing smartdenovo ..."
smartdenovo_dir="$build_dir/smartdenovo"
if [ -z $(check_installed $smartdenovo_dir) ]; then
    cd $build_dir
    clean "$build_dir/smartdenovo"
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

# # --------------- Raven ------------------
# echo ""
# echo "[$(timestamp)] Installing raven ..."
# raven_dir="$build_dir/raven/build/bin"
# if [ -z $(check_installed $raven_dir) ]; then
#     cd $build_dir
#     clean "$build_dir/raven"
#     echo "Download Raven-v${RAVEN_VERSION}"
#     git clone --recursive https://github.com/lbcb-sci/raven.git raven
#     cd raven
#     git checkout -f -q $RAVEN_GITHUB_COMMIT_VERSION
#     mkdir build
#     cd build
#     cmake -DCMAKE_BUILD_TYPE=Release ..   #requires cmake >=3.11
#     make -j $MAKE_JOBS
#     note_installed $raven_dir
# fi

# --------------- Shasta ------------------
echo ""
echo "[$(timestamp)] Installing shasta ..."
shasta_dir="$build_dir/shasta-${SHASTA_VERSION}"
if [ -z $(check_installed $shasta_dir) ]; then
    cd $build_dir
    clean "$build_dir/shasta-${SHASTA_VERSION}"
    echo "Download Shasta-v${SHASTA_VERSION}"
    mkdir shasta-${SHASTA_VERSION}
    cd shasta-${SHASTA_VERSION}
    wget $SHASTA_DOWNLOAD_URL
    chmod ugo+x shasta-Linux-${SHASTA_VERSION}
    ln -s shasta-Linux-${SHASTA_VERSION} shasta
    note_installed $shasta_dir
fi

# --------------- Ragout ------------------
echo ""
echo "[$(timestamp)] Installing ragout ..."
ragout_dir="$build_dir/ragout_conda_env/bin"
#ragout_dir="$build_dir/ragout_conda_env/Ragout/bin"
if [ -z $(check_installed $ragout_dir) ]; then
    cd $build_dir
    clean "$build_dir/ragout_conda_env"
    $miniconda3_dir/conda create -y -p $build_dir/ragout_conda_env python=3.8
    source $miniconda3_dir/activate $build_dir/ragout_conda_env
    $miniconda3_dir/conda install -y -c bioconda ragout=${RAGOUT_VERSION}

    # cd $build_dir/ragout_conda_env
    # git clone $RAGOUT_DOWNLOAD_URL
    # cd Ragout
    # git checkout -f -q $RAGOUT_GITHUB_COMMIT_VERSION 
    # python setup.py build
    # pip install -r requirements.txt
    # python scripts/install-sibelia.py

    source $miniconda3_dir/deactivate
    note_installed $ragout_dir
fi

# --------------- Ragtag ------------------
echo ""
echo "[$(timestamp)] Installing ragtag ..."
ragtag_dir="$build_dir/ragtag_conda_env/bin"
if [ -z $(check_installed $ragtag_dir) ]; then
    cd $build_dir
    clean "$build_dir/ragout_ragtag_env"
    $miniconda3_dir/conda create -y -p $build_dir/ragtag_conda_env python=3.7 pip
    source $miniconda3_dir/activate $build_dir/ragtag_conda_env
    #$miniconda3_dir/conda install -y -c bioconda ragtag=${RAGTAG_VERSION} 
    $ragtag_dir/pip3 install --timeout 1000 ragtag==${RAGTAG_VERSION}
    cd $ragtag_dir
    ln -s $minimap2_dir/minimap2 .
    source $miniconda3_dir/deactivate
    note_installed $ragtag_dir
fi

# --------------- gnuplot ------------------
echo ""
echo "[$(timestamp)] Installing gnuplot ..."
gnuplot_dir="$build_dir/gnuplot-${GNUPLOT_VERSION}/bin"
if [ -z $(check_installed $gnuplot_dir) ]; then
    cd $build_dir
    clean "$build_dir/gnuplot-${GNUPLOT_VERSION}"
    echo "Download gnuplot-v${GNUPLOT_VERSION}"
    download $GNUPLOT_DOWNLOAD_URL "gnuplot-${GNUPLOT_VERSION}.tar.gz"
    tar xvzf gnuplot-${GNUPLOT_VERSION}.tar.gz
    cd "$build_dir/gnuplot-${GNUPLOT_VERSION}"
    ./configure --prefix="$build_dir/gnuplot-${GNUPLOT_VERSION}" --exec-prefix="$build_dir/gnuplot-${GNUPLOT_VERSION}" --disable-plugins --disable-wxwidgets --with-gd=no --without-latex --with-x=no
    make -j $MAKE_JOBS
    make install
    cd $build_dir
    rm gnuplot-${GNUPLOT_VERSION}.tar.gz
    note_installed $gnuplot_dir
fi
PATH="$gnuplot_dir:${PATH}"

# --------------- Guppy-GPU --------------------
echo ""
echo "[$(timestamp)] Installing guppy-GPU ..."
guppy_gpu_dir="$build_dir/ont-guppy-gpu/bin"
if [ -z $(check_installed $guppy_gpu_dir) ]; then
    cd $build_dir
    clean "$build_dir/ont-guppy-gpu"
    echo "Download Guppy-v${GUPPY_VERSION}"
    download $GUPPY_GPU_DOWNLOAD_URL "ont-guppy_${GUPPY_VERSION}_linux64.tar.gz"
    tar xvzf ont-guppy_${GUPPY_VERSION}_linux64.tar.gz
    mv ont-guppy ont-guppy-gpu
    rm ont-guppy_${GUPPY_VERSION}_linux64.tar.gz
    note_installed $guppy_gpu_dir
fi

# --------------- Guppy-CPU --------------------
echo ""
echo "[$(timestamp)] Installing guppy-CPU ..."
guppy_cpu_dir="$build_dir/ont-guppy-cpu/bin"
if [ -z $(check_installed $guppy_cpu_dir) ]; then
    cd $build_dir
    clean "$build_dir/ont-guppy-cpu"
    echo "Download Guppy-v${GUPPY_VERSION}"
    download $GUPPY_CPU_DOWNLOAD_URL "ont-guppy_${GUPPY_VERSION}_linux64.tar.gz"
    tar xvzf ont-guppy_${GUPPY_VERSION}_linux64.tar.gz
    rm ont-guppy_${GUPPY_VERSION}_linux64.tar.gz
    note_installed $guppy_cpu_dir
fi


# # --------------- Nanofilt --------------------
# echo ""
# echo "[$(timestamp)] Installing NanoFilt ..."
# nanofilt_dir="$build_dir/nanofilt_conda_env/bin"
# if [ -z $(check_installed $nanofilt_dir) ]; then
#     cd $build_dir
#     $miniconda3_dir/conda create -y -p $build_dir/nanofilt_conda_env 
#     source $miniconda3_dir/activate $build_dir/nanofilt_conda_env
#     $miniconda3_dir/conda install -y -c bioconda nanofilt=${NANOFILT_VERSION}
#     source $miniconda3_dir/deactivate
#     note_installed $nanofilt_dir
# fi

# --------------- Nanoplot --------------------
echo ""
echo "[$(timestamp)] Installing NanoPlot ..."
nanoplot_dir="$build_dir/nanoplot_conda_env/bin"
if [ -z $(check_installed $nanoplot_dir) ]; then
    cd $build_dir
    $miniconda3_dir/conda create -y -p $build_dir/nanoplot_conda_env python=3.9
    source $miniconda3_dir/activate $build_dir/nanoplot_conda_env
    $miniconda3_dir/conda install -y -c bioconda nanoplot=${NANOPLOT_VERSION}
    source $miniconda3_dir/deactivate
    note_installed $nanoplot_dir
fi

# --------------- Nanopolish --------------------
echo ""
echo "[$(timestamp)] Installing nanopolish ..."
nanopolish_dir="$build_dir/nanopolish_conda_env/bin"
if [ -z $(check_installed $nanopolish_dir) ]; then
    cd $build_dir
    clean "$build_dir/nanopolish_conda_env"
    $miniconda3_dir/conda create -y -p $build_dir/nanopolish_conda_env python=3.9
    source $miniconda3_dir/activate $build_dir/nanopolish_conda_env
    $miniconda3_dir/conda install -y -c bioconda nanopolish=${NANOPOLISH_VERSION}
    source $miniconda3_dir/deactivate
    note_installed $nanopolish_dir
fi

# --------------- Parallel ------------------
echo ""
echo "[$(timestamp)] Installing parallele ..."
parallel_dir="$build_dir/parallel-${PARALLEL_VERSION}/bin"
if [ -z $(check_installed $parallel_dir) ]; then
    cd $build_dir
    echo "Download parallel"
    clean "$build_dir/parallele-${PARALLEL_VERSION}"
    download $PARALLEL_DOWNLOAD_URL "parallel_v${PARALLEL_VERSION}.tar.bz2"
    tar xvjf parallel_v${PARALLEL_VERSION}.tar.bz2
    cd parallel-${PARALLEL_VERSION}
    ./configure --prefix="$build_dir/parallel-${PARALLEL_VERSION}"
    make -j $MAKE_JOBS
    make install
    cd $build_dir
    rm parallel_v${PARALLEL_VERSION}.tar.bz2
    note_installed $parallel_dir
fi

# --------------- Racon -----------------
echo ""
echo "[$(timestamp)] Installing racon ..."
racon_dir="$build_dir/racon_conda_env/bin"
if [ -z $(check_installed $racon_dir) ]; then
    cd $build_dir
    clean "$build_dir/racon_conda_env"
    $miniconda3_dir/conda create -y -p $build_dir/racon_conda_env python=3.7
    source $miniconda3_dir/activate $build_dir/racon_conda_env
    $miniconda3_dir/conda install -y -c bioconda racon=${RACON_VERSION}
    source $miniconda3_dir/deactivate
    note_installed $racon_dir
fi

# --------------- Medaka -----------------
echo ""
echo "[$(timestamp)] Installing medaka ..."
medaka_dir="$build_dir/medaka_conda_env/bin"
if [ -z $(check_installed $medaka_dir) ]; then
    cd $build_dir
    clean "$build_dir/medaka_conda_env"
    $miniconda3_dir/conda create -y -p $build_dir/medaka_conda_env 
    source $miniconda3_dir/activate $build_dir/medaka_conda_env
    $miniconda3_dir/conda install -y -c conda-forge -c bioconda medaka=${MEDAKA_VERSION}
    #$miniconda3_dir/pip3 install medaka==${MEDAKA_VERSION}
    source $miniconda3_dir/deactivate
    note_installed $medaka_dir
fi

# # --------------- MarginPolish -----------------
# echo ""
# echo "[$(timestamp)] Installing marginpolish ..."
# marginpolish_dir="$build_dir/marginPolish/build"
# if [ -z $(check_installed $marginpolish_dir) ]; then
#     cd $build_dir
#     echo "Download MarginPolish-v${MARGINPOLISH_VERSION}"
#     clean "$build_dir/marginPolish"
#     download $MARGINPOLISH_DOWNLOAD_URL "marginPolish.tar.gz"
#     tar xzf marginPolish.tar.gz
#     cd MarginPolish-${MARGINPOLISH_VERSION}/externalTools
#     git clone --recursive https://github.com/samtools/htslib.git
#     cd htslib
#     git checkout -f -q 1832d3a
#     cd ..
#     cd git clone --resursive https://github.com/benedictpaten/sonLib.git
#     git checkout -f -q 8e403f6
#     cd ..
#     cd ..

#     # git clone --recursive $MARGINPOLISH_DOWNLOAD_URL
#     # cd marginPolish
#     # git checkout -f -q $MARGINPOLISH_GITHUB_COMMIT_VERSION
#     # # git submodule update --init
#     mkdir build
#     cd build
#     cmake ..
#     make -j $MAKE_JOBS
#     note_installed $marginpolish_dir
# fi

# # --------------- Helen -----------------
# echo ""
# echo "[$(timestamp)] Installing helen ..."
# helen_dir="$build_dir/helen/build"
# if [ -z $(check_installed $helen_dir) ]; then
#     cd $build_dir
#     echo "Download Helen-v${HELEN_VERSION}"
#     clean "$build_dir/helen"
#     git clone $HELEN_DOWNLOAD_URL
#     cd helen
#     git checkout -f -q $HELEN_GITHUB_COMMIT_VERSION
#     make install
#     . ./venv/bin/activate
#     cd build
#     note_installed $helen_dir
# fi

# # --------------- PEPPER -----------------
# pepper_dir="$build_dir/pepper_conda_env/bin"
# if [ -z $(check_installed $pepper_dir) ]; then
#     cd $build_dir
#     $miniconda3_dir/conda create -y -p $build_dir/pepper_conda_env python=3.7
#     source $miniconda3_dir/activate $build_dir/pepper_conda_env
#     $miniconda3_dir/pip3 install pepper-polish
#     source $miniconda3_dir/deactivate
#     note_installed $pepper_dir
# fi

# # --------------- HOMOPOLISH -----------------
# homopolish_dir="$build_dir/homopolish_conda_env/bin"
# if [ -z $(check_installed $homopolish_dir) ]; then
#     cd $build_dir
#     clean "$build_dir/homopolish_conda_env"
#     $miniconda3_dir/conda create -y -p $build_dir/homopolish_conda_env python=3.7
#     source $miniconda3_dir/activate $build_dir/homopolish_conda_env
#     $miniconda3_dir/conda install -y -c bioconda homopolish=${HOMOPOLISH_VERSION}
#     source $miniconda3_dir/deactivate
#     note_installed $homopolish_dir
# fi

# # ------------- QUAST --------------------
# quast_dir="$build_dir/quast_conda_env/bin"
# if [ -z $(check_installed $quast_dir) ]; then
#     cd $build_dir
#     $miniconda3_dir/conda create -y -p $build_dir/quast_conda_env python=3.7
#     source $miniconda3_dir/activate $build_dir/quast_conda_env
#     $miniconda3_dir/conda install -y -c bioconda quast=${QUAST_VERSION}
#     source $miniconda3_dir/deactivate
#     note_installed $quast_dir
# fi

# --------------- mummer4 ------------------
echo ""
echo "[$(timestamp)] Installing mummer4 ..."
mummer4_dir="$build_dir/mummer-${MUMMER4_VERSION}"
if [ -z $(check_installed $mummer4_dir) ]; then
    cd $build_dir
    clean "$build_dir/mummer-${MUMMER4_VERSION}"
    echo "Download mummer-v${MUMMER4_VERSION}"
    download $MUMMER4_DOWNLOAD_URL "mummer-${MUMMER4_VERSION}.tar.gz"
    tar xvzf mummer-${MUMMER4_VERSION}.tar.gz
    echo "$mummer4_dir"
    cd $mummer4_dir
    ./configure
    make -j $MAKE_JOBS
    cd $build_dir
    rm mummer-${MUMMER4_VERSION}.tar.gz
    note_installed $mummer4_dir
fi
PATH="$mummer4_dir:${PATH}"

# --------------- bedtools ------------------
echo ""
echo "[$(timestamp)] Installing bedtools ..."
bedtools_dir="$build_dir/bedtools2/bin"
if [ -z $(check_installed $bedtools_dir) ]; then
    cd $build_dir
    clean "$build_dir/bedtools2"
    echo "Download bedtools-v${BEDTOOLS_VERSION}"
    download $BEDTOOLS_DOWNLOAD_URL "bedtools-${BEDTOOLS_VERSION}.tar.gz"
    tar xvzf bedtools-${BEDTOOLS_VERSION}.tar.gz
    cd "$build_dir/bedtools2"
    make -j $MAKE_JOBS
    cd $build_dir
    rm bedtools-${BEDTOOLS_VERSION}.tar.gz
    note_installed $bedtools_dir
fi

# --------------- SPAdes ------------------
echo ""
echo "[$(timestamp)] Installing spades ..."                                                                                       
spades_dir="$build_dir/SPAdes-${SPADES_VERSION}-Linux/bin"
if [ -z $(check_installed $spades_dir) ]; then
    cd $build_dir
    clean "$build_dir/SPAdes-${SPADES_VERSION}-Linux"
    echo "Download SPAdes-v${SPADES_VERSION}"
    download $SPADES_DOWNLOAD_URL "SPAdes-${SPADES_VERSION}-Linux.tar.gz"
    tar xvzf SPAdes-${SPADES_VERSION}-Linux.tar.gz
    rm SPAdes-${SPADES_VERSION}-Linux.tar.gz
    note_installed $spades_dir
fi

# --------------- Prodigal ------------------
prodigal_dir="$build_dir/Prodigal-${PRODIGAL_VERSION}"
echo ""
echo "[$(timestamp)] Installing Prodigal ..."
if [ -z $(check_installed $prodigal_dir) ]; then
    cd $build_dir
    clean "$build_dir/Prodigal-${PRODIGAL_VERSION}"
    echo "Download Prodigal-v${PRODIGAL_VERSION}"
    download $PRODIGAL_DOWNLOAD_URL "v${PRODIGAL_VERSION}.tar.gz"
    tar xvzf v${PRODIGAL_VERSION}.tar.gz
    cd $prodigal_dir
    make -j $MAKE_JOBS
    cd $build_dir
    rm v${PRODIGAL_VERSION}.tar.gz
    note_installed $prodigal_dir
fi

# --------------- CAP3 ------------------
cap3_dir="$build_dir/CAP3/CAP3"
echo ""
echo "[$(timestamp)] Installing CAP3 ..."
if [ -z $(check_installed $cap3_dir) ]; then
    cd $build_dir
    clean "$build_dir/CAP3"
    echo "Download CAP3-v${CAP3_VERSION}"
    # download $CAP3_DOWNLOAD_URL "cap3.linux.x86_64.tar"
    git clone $CAP3_DOWNLOAD_URL
    cd CAP3
    tar xvf cap3.linux.x86_64.tar
    rm cap3.linux.x86_64.tar
    note_installed $cap3_dir
fi

# ------------- BWA -------------------
bwa_dir="$build_dir/bwa-${BWA_VERSION}"
echo ""
echo "[$(timestamp)] Installing BWA ..."
if [ -z $(check_installed $bwa_dir) ]; then
    cd $build_dir
    clean "$build_dir/bwa-${BWA_VERSION}"
    echo "Download BWA-v${BWA_VERSION}"
    download $BWA_DOWNLOAD_URL "bwa-${BWA_VERSION}.tar.bz2"
    tar xvjf bwa-${BWA_VERSION}.tar.bz2
    cd $bwa_dir
    make -j $MAKE_JOBS
    cd $build_dir
    rm bwa-${BWA_VERSION}.tar.bz2
    note_installed $bwa_dir
fi

# --------------- Circlator ------------------
circlator_dir="$build_dir/circlator_conda_env/bin"
echo ""
echo "[$(timestamp)] Installing Circlator ..."
if [ -z $(check_installed $circlator_dir) ]; then
    cd $build_dir
    clean "$build_dir/circlator_conda_env"
    $miniconda3_dir/conda create -y -p $build_dir/circlator_conda_env
    source $miniconda3_dir/activate $build_dir/circlator_conda_env
    $miniconda3_dir/conda install -y -c bioconda circlator=${CIRCLATOR_VERSION}
    source $miniconda3_dir/deactivate
    note_installed $circlator_dir
fi

# --------------- Trimmomatic -----------------
trimmomatic_dir="$build_dir/Trimmomatic-${TRIMMOMATIC_VERSION}"
echo ""
echo "[$(timestamp)] Installing Trimmomatic ..."
if [ -z $(check_installed $trimmomatic_dir) ]; then
    cd $build_dir
    clean "$build_dir/Trimmomatic-${TRIMMOMATIC_VERSION}"
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
echo ""
echo "[$(timestamp)] Installing Picard ..."
if [ -z $(check_installed $picard_dir) ]; then
    cd $build_dir
    clean "$build_dir/Picard-v${PICARD_VERSION}"
    echo "Download Picard-v${PICARD_VERSION}"
    download $PICARD_DOWNLOAD_URL "picard.jar"
    mkdir Picard-v${PICARD_VERSION}

    mv picard.jar $picard_dir
    cd $picard_dir
    chmod 755 picard.jar
    note_installed $picard_dir
fi

# # --------------- hapo-G ------------------
# hapog_dir="$build_dir/hapog_conda_env/bin"
# echo ""
# echo "[$(timestamp)] Installing hapo-G ..."
# if [ -z $(check_installed $hapog_dir) ]; then
#     cd $build_dir
#     $miniconda3_dir/conda create -y -p $build_dir/hapog_conda_env
#     source $miniconda3_dir/activate $build_dir/hapog_conda_env
#     $miniconda3_dir/conda install -y -c lbgb_cea hapog=${HAPOG_VERSION}
#     source $miniconda3_dir/deactivate
#     note_installed $hapog_dir
# fi

# --------------- Pilon -----------------
pilon_dir="$build_dir/Pilon-v${PILON_VERSION}"
echo ""
echo "[$(timestamp)] Installing Pilon ..."
if [ -z $(check_installed $pilon_dir) ]; then
    cd $build_dir
    clean "$build_dir/Pilon-v${PILON_VERSION}"
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
echo ""
echo "[$(timestamp)] Installing exonerate ..."
if [ -z $(check_installed $exonerate_dir) ]; then
    cd $build_dir
    clean "$build_dir/exonerate-${EXONERATE_VERSION}-x86_64"
    echo "Download exonerate-v${EXONERATE_VERSION}"
    download $EXONERATE_DOWNLOAD_URL "exonerate-${EXONERATE_VERSION}-x86_64.tar.gz"
    tar xvzf exonerate-${EXONERATE_VERSION}-x86_64.tar.gz
    rm exonerate-${EXONERATE_VERSION}-x86_64.tar.gz
    note_installed $exonerate_dir
fi

# --------------- ncbi-blast+ ------------------
blast_dir="$build_dir/ncbi-blast-${BLAST_VERSION}+/bin"
blast_matrices_dir="$build_dir/ncbi-blast-${BLAST_VERSION}+/matrices"
echo ""
echo "[$(timestamp)] Installing ncbi-blast+ ..."
if [ -z $(check_installed $blast_dir) ]; then
    cd $build_dir
    clean "$build_dir/ncbi-blast-${BLAST_VERSION}+"
    echo "Download ncbi-blast-v${BLAST_VERSION}"
    download $BLAST_DOWNLOAD_URL "ncbi-blast-${BLAST_VERSION}+-x64-linux.tar.gz"
    tar xvzf ncbi-blast-${BLAST_VERSION}+-x64-linux.tar.gz
    cd "ncbi-blast-${BLAST_VERSION}+"
    mkdir matrices
    cd matrices
    wget -c --no-check-certificate ftp://ftp.ncbi.nlm.nih.gov/blast/matrices/*
    cd $build_dir
    rm ncbi-blast-${BLAST_VERSION}+-x64-linux.tar.gz
    note_installed $blast_dir
fi

# --------------- ncbi-rmblast ------------------
rmblast_dir="$build_dir/ncbi-rmblastn-${RMBLAST_VERSION}/bin"
echo ""
echo "[$(timestamp)] Installing rmblast ..."
if [ -z $(check_installed $rmblast_dir) ]; then
    cd $build_dir
    clean "$build_dir/ncbi-rmblastn-${RMBLAST_VERSION}"
    echo "Download ncbi-rmblastn-v${BLAST_VERSION}"
    download $RMBLAST_DOWNLOAD_URL "ncbi-rmblastn-${RMBLAST_VERSION}-x64-linux.tar.gz"
    tar xvzf ncbi-rmblastn-${RMBLAST_VERSION}-x64-linux.tar.gz
    # copy rmblastn binary file to ncbi-blast+ directory for easy RepeatMasker configuration
    cp $rmblast_dir/rmblastn $blast_dir
    rm ncbi-rmblastn-${RMBLAST_VERSION}-x64-linux.tar.gz
    note_installed $rmblast_dir
fi

# --------------- snap ------------------
snap_dir="$build_dir/SNAP"
echo ""
echo "[$(timestamp)] Installing snap ..."
if [ -z $(check_installed $snap_dir) ]; then
    cd $build_dir
    clean "$build_dir/SNAP"
    echo "Download snap-v${SNAP_VERSION}"
    git clone $SNAP_DOWNLOAD_URL
    cd $snap_dir
    git checkout -f -q $SNAP_GITHUB_COMMIT_VERSION
    ZOE="$snap_dir/Zoe"
    cp $LRSDAY_HOME/misc/snap.c .  # temporary fix for snap with gcc-8
    make -j $MAKE_JOBS
    note_installed $snap_dir
fi

# # --------------- RAPSearch2 ------------------
# rapsearch_dir="$build_dir/RAPSearch${RAPSEARCH_VERSION}_64bits/bin"
# echo ""
# echo "[$(timestamp)] Installing RAPSearch2 ..."
# if [ -z $(check_installed $rapsearch_dir) ]; then
#     cd $build_dir
#     echo "Download RAPsearch-v${RAPSEARCH_VERSION}"
#     download $RAPSEARCH_DOWNLOAD_URL "RAPSearch${RAPSEARCH_VERSION}_64bits.tar.gz"
#     tar xvzf RAPSearch${RAPSEARCH_VERSION}_64bits.tar.gz
#     rm RAPSearch${RAPSEARCH_VERSION}_64bits.tar.gz
#     note_installed $rapsearch_dir
# fi

# --------------- tRNAscan-SE ------------------
trnascan_root="$build_dir/tRNAscan-SE-${TRNASCAN_VERSION}"
trnascan_dir="$trnascan_root/bin"
PERL5LIB=$trnascan_dir:$PERL5LIB
echo ""
echo "[$(timestamp)] Installing tRNAscan-SE ..."
if [ -z $(check_installed $trnascan_dir) ]; then
    cd $build_dir
    clean "$build_dir/tRNAscan-SE-${TRNASCAN_VERSION}"
    echo "Download tRNAscan-SE-v${TRNASCAN_VERSION}"
    download $TRNASCAN_DOWNLOAD_URL "tRNAscan-SE-${TRNASCAN_VERSION}.tar.gz"
    tar xvzf tRNAscan-SE-${TRNASCAN_VERSION}.tar.gz
    cd $trnascan_root
    mkdir bin
    mkdir -p lib/tRNAscan-SE
    mkdir -p man
    cp $LRSDAY_HOME/misc/tRNAscan-SE.Makefile Makefile
    make -j $MAKE_JOBS BINDIR="$trnascan_root/bin" LIBDIR="$trnascan_root/lib/tRNAscan-SE" MANDIR="$trnascan_root"
    make install BINDIR="$trnascan_root/bin" LIBDIR="$trnascan_root/lib/tRNAscan-SE" MANDIR="$trnascan_root"
    cd $build_dir
    rm tRNAscan-SE-${TRNASCAN_VERSION}.tar.gz
    note_installed $trnascan_dir
fi

# --------------- snoscan ------------------
snoscan_dir="$build_dir/snoscan-${SNOSCAN_VERSION}"
echo ""
echo "[$(timestamp)] Installing snoscan ..."
if [ -z $(check_installed $snoscan_dir) ]; then
    cd $build_dir
    clean "$build_dir/snoscan-${SNOSCAN_VERSION}"
    echo "Download snoscan-v${SNOSCAN_VERSION}"
    download $SNOSCAN_DOWNLOAD_URL "snoscan-${SNOSCAN_VERSION}.tar.gz"
    tar xvzf snoscan-${SNOSCAN_VERSION}.tar.gz
    cd $snoscan_dir
    cd squid-1.5.11
    rm *.o
    make -j $MAKE_JOBS
    cd ..
    cp $LRSDAY_HOME/misc/snoscan.Makefile Makefile
    rm *.o
    make -j $MAKE_JOBS
    cd $build_dir
    rm snoscan-${SNOSCAN_VERSION}.tar.gz
    note_installed $snoscan_dir
fi

# --------------- RepeatMasker ------------------
repeatmasker_dir="$build_dir/RepeatMasker"
echo ""
echo "[$(timestamp)] Installing RepeatMasker ..."
if [ -z $(check_installed $repeatmasker_dir) ]; then
    cd $build_dir
    clean "$build_dir/RepeatMasker"
    echo "Download Repeatmasker-v${REPEATMASKER_VERSION}"
    download $REPEATMASKER_DOWNLOAD_URL "RepeatMasker-${REPEATMASKER_VERSION}.tar.gz"
    tar xvzf RepeatMasker-${REPEATMASKER_VERSION}.tar.gz
    cd $repeatmasker_dir
    echo "Download and setup RepBase library"
    # download $REPBASE_DOWNLOAD_URL "RepBaseRepeatMaskerEdition-${REPBASE_VERSION}.tar.gz"
    git clone $REPBASE_DOWNLOAD_URL
    mv ./RMRB/RepBaseRepeatMaskerEdition-${REPBASE_VERSION}.tar.gz .
    tar xzf RepBaseRepeatMaskerEdition-${REPBASE_VERSION}.tar.gz
    rm RepBaseRepeatMaskerEdition-${REPBASE_VERSION}.tar.gz
    cd $build_dir
    rm RepeatMasker-${REPEATMASKER_VERSION}.tar.gz
    note_installed $repeatmasker_dir
fi

# --------------- TRF ------------------
trf_dir=$repeatmasker_dir
echo ""
echo "[$(timestamp)] Installing trf ..."
if [ ! -e "$repeatmasker_dir/trf" ]; then
    cd $repeatmasker_dir
    echo "Download TRF-v${TRF_VERSION}"
    download $TRF_DOWNLOAD_URL "trf${TRF_VERSION}.linux64"
    mv trf${TRF_VERSION}.linux64 trf
    chmod 755 trf
fi

# --------------- REannotate ------------------
reannotate_dir="$build_dir/REannotate_LongQueryName"
echo ""
echo "[$(timestamp)] Installing REannotate ..."
if [ -z $(check_installed $reannotate_dir) ]; then
    cd $build_dir
    clean "$build_dir/REannotate_LongQueryName"
    echo "Download REannotate-v${REANNOTATE_VERSION}"
    git clone $REANNOTATE_DOWNLOAD_URL
    cd REannotate_LongQueryName 
    git checkout -f -q $REANNOTATE_GITHUB_COMMIT_VERSION
    chmod 755 REannotate_longname
    ln -s REannotate_longname REannotate
    cd $build_dir
    note_installed $reannotate_dir
fi

# --------------- ClustalW ------------------
clustalw_dir="$build_dir/clustalw-${CLUSTALW_VERSION}/bin"
echo ""
echo "[$(timestamp)] Installing ClustalW ..."
if [ -z $(check_installed $clustalw_dir) ]; then
    cd $build_dir
    clean "$build_dir/clustalw-${CLUSTALW_VERSION}"
    echo "Download ClustalW-v${CLUSTALW_VERSION}"
    download $CLUSTALW_DOWNLOAD_URL "clustalw-${CLUSTALW_VERSION}.tar.gz"
    tar xvzf clustalw-${CLUSTALW_VERSION}.tar.gz
    cd clustalw-${CLUSTALW_VERSION}
    ./configure --prefix="$build_dir/clustalw-${CLUSTALW_VERSION}"
    make -j $MAKE_JOBS
    make install
    cd $build_dir
    rm clustalw-${CLUSTALW_VERSION}.tar.gz
    note_installed $clustalw_dir
fi

# --------------- MUSCLE ------------------
muscle_dir="$build_dir/muscle-${MUSCLE_VERSION}"
echo ""
echo "[$(timestamp)] Installing muscle ..."
if [ -z $(check_installed $muscle_dir) ]; then
    cd $build_dir
    clean "$build_dir/muscle-${MUSCLE_VERSION}"
    echo "Download MUSCLE-v${MUSCLE_VERSION}"
    download $MUSCLE_DOWNLOAD_URL "muscle-${MUSCLE_VERSION}_i86linux64.tar.gz"
    tar xvzf muscle-${MUSCLE_VERSION}_i86linux64.tar.gz
    mkdir muscle-${MUSCLE_VERSION}
    mv muscle${MUSCLE_VERSION}_i86linux64 ./muscle-${MUSCLE_VERSION}/muscle
    cd $build_dir
    rm muscle-${MUSCLE_VERSION}_i86linux64.tar.gz
    note_installed $muscle_dir
fi

# --------------- HMMER ------------------
hmmer_dir="$build_dir/hmmer-${HMMER_VERSION}/bin"
echo ""
echo "[$(timestamp)] Installing HMMER ..."
if [ -z $(check_installed $hmmer_dir) ]; then
    cd $build_dir
    clean "$build_dir/hmmer-${HMMER_VERSION}"
    echo "Download hmmer-v${HMMER_VERSION}"
    download $HMMER_DOWNLOAD_URL "hmmer-${HMMER_VERSION}-linux-intel-x86_64.tar.gz"
    tar xvzf hmmer-${HMMER_VERSION}-linux-intel-x86_64.tar.gz
    cd hmmer-${HMMER_VERSION}
    ./configure --prefix="$build_dir/hmmer-${HMMER_VERSION}"
    make -j $MAKE_JOBS
    make install
    cd easel
    make install
    cd $build_dir
    rm hmmer-${HMMER_VERSION}-linux-intel-x86_64.tar.gz
    note_installed $hmmer_dir
fi

# --------------- bamtools ------------------
bamtools_dir="$build_dir/bamtools-${BAMTOOLS_VERSION}/bin"
echo ""
echo "[$(timestamp)] Installing bamtools ..."
if [ -z $(check_installed $bamtools_dir) ]; then
    cd $build_dir
    clean "$build_dir/bamtools-${BAMTOOLS_VERSION}"
    echo "Download bamtools-v${BAMTOOLS_VERSION}"
    download $BAMTOOLS_DOWNLOAD_URL "v${BAMTOOLS_VERSION}.tar.gz"
    tar xvzf v${BAMTOOLS_VERSION}.tar.gz
    cd bamtools-${BAMTOOLS_VERSION}
    mkdir build
    cd build
    cmake -DCMAKE_INSTALL_PREFIX="$build_dir/bamtools-${BAMTOOLS_VERSION}" ..
    make -j $MAKE_JOBS
    make install
    cd $build_dir/bamtools-${BAMTOOLS_VERSION}
    ln -sf lib lib64
    cd $build_dir
    rm v${BAMTOOLS_VERSION}.tar.gz
    note_installed $bamtools_dir
fi

# --------------- Augustus ------------------
augustus_dir="$build_dir/augustus-${AUGUSTUS_VERSION}/bin"
echo ""
echo "[$(timestamp)] Installing augustus ..."
if [ -z $(check_installed $augustus_dir) ]; then
    cd $build_dir
    clean "$build_dir/augustus-${AUGUSTUS_VERSION}"
    echo "Download Augustus-v${AUGUSTUS_VERSION}"
    download $AUGUSTUS_DOWNLOAD_URL "augustus-${AUGUSTUS_VERSION}.tar.gz"
    tar xvzf augustus-${AUGUSTUS_VERSION}.tar.gz
    cd $build_dir/augustus-${AUGUSTUS_VERSION}/auxprogs/bam2hints/
    cp $LRSDAY_HOME/misc/bam2hints.Makefile Makefile
    cd $build_dir/augustus-${AUGUSTUS_VERSION}/auxprogs/filterBam/src/
    cp $LRSDAY_HOME/misc/filterBam.Makefile Makefile
    cd $build_dir/augustus-${AUGUSTUS_VERSION}
    make -j $MAKE_JOBS BAMTOOLS="$build_dir/bamtools-${BAMTOOLS_VERSION}"
    cd $build_dir
    rm augustus-${AUGUSTUS_VERSION}.tar.gz
    note_installed $augustus_dir
fi
export AUGUSTUS_CONFIG_PATH="$build_dir/augustus-${AUGUSTUS_VERSION}/config"

# --------------- EVidenceModeler ------------------
evm_dir="$build_dir/EVidenceModeler-${EVM_VERSION}"
echo ""
echo "[$(timestamp)] Installing EVM ..."
if [ -z $(check_installed $evm_dir) ]; then
    cd $build_dir
    clean "$build_dir/EVidenceModeler-${EVM_VERSION}"
    echo "Download EvidenceModeler-v${EVM_VERSION}"
    download $EVM_DOWNLOAD_URL "v${EVM_VERSION}.tar.gz"
    tar xvzf v${EVM_VERSION}.tar.gz
    rm v${EVM_VERSION}.tar.gz
    note_installed $evm_dir
fi

# --------------- Proteinortho ------------------
proteinortho_dir="$build_dir/proteinortho_conda_env/bin"
diamond_dir="$build_dir/proteinortho_conda_env/bin"
echo ""
echo "[$(timestamp)] Installing proteinortho ..."
if [ -z $(check_installed $proteinortho_dir) ]; then
    cd $build_dir
    clean "$build_dir/proteinortho_conda_env"
    echo "Download Porteinortho-v${PROTEINORTHO_VERSION}"
    $miniconda3_dir/conda create -y -p $build_dir/proteinortho_conda_env 
    source $miniconda3_dir/activate $build_dir/proteinortho_conda_env
    $miniconda3_dir/conda install -y proteinortho=${PROTEINORTHO_VERSION}
    # $miniconda3_dir/conda install -y -c bioconda diamond=${DIAMOND_VERSION}
    # $miniconda3_dir/conda install -y -c bioconda blast=2.10.1
    # $miniconda3_dir/conda install -y -c bioconda last=${LAST_VERSION}
    source $miniconda3_dir/deactivate
    note_installed $proteinortho_dir
fi

# --------------- GATK3 ------------------
gatk3_dir="$build_dir/GATK3"
echo ""
echo "[$(timestamp)] Installing gatk3 ..."
if [ -z $(check_installed $gatk3_dir) ]; then
    cd $build_dir
    clean "$build_dir/GATK3"
    echo "Create the GATK3 folder for installation"
    mkdir GATK3
    cd GATK3
    #git clone $GATK3_DOWNLOAD_URL 
    download $GATK3_DOWNLOAD_URL gatk3.jar.gz
    gunzip gatk3.jar.gz 
    mv gatk3.jar GenomeAnalysisTK.jar 
    chmod 755 GenomeAnalysisTK.jar
    note_installed $gatk3_dir
fi

# # --------------- MAKER -----------------
# maker_dir="$build_dir/maker_conda_env/bin"
# echo ""
# echo "[$(timestamp)] Installing maker ..."
# if [ -z $(check_installed $maker_dir) ]; then
#     cd $build_dir
#     echo "Download MAKER-${MAKER_VERSION}"
#     $miniconda3_dir/conda create -y -p $build_dir/maker_conda_env python=3.8 
#     source $miniconda3_dir/activate $build_dir/maker_conda_env
#     $miniconda3_dir/conda install -y -c bioconda maker==${MAKER_VERSION}    
#     source $miniconda3_dir/deactivate
#     note_installed $maker_dir
# fi

# --------------- MAKER -----------------
maker_dir="$build_dir/maker/bin"
if [ -z $(check_installed $maker_dir) ]; then
    cd $build_dir
    clean "$build_dir/maker"
    echo "Download MAKER"
    cp -r $LRSDAY_HOME/misc/maker-${MAKER_VERSION}.tar.gz .
    tar xvzf maker-${MAKER_VERSION}.tar.gz
    mv maker-${MAKER_VERSION} maker
    cd $build_dir/maker/src/
    cp $LRSDAY_HOME/misc/maker_Build.PL .
    echo "no"|perl maker_Build.PL
    ./Build install
    cd $build_dir
    rm maker-${MAKER_VERSION}.tar.gz
    note_installed $maker_dir
fi


# # --------------- Gff3toolkit ------------------
# gff3toolkit_dir="$build_dir/gff3toolkit_conda_env/bin"
# echo ""
# echo "[$(timestamp)] Installing gff3toolkit ..."
# if [ -z $(check_installed $gff3toolkit_dir) ]; then
#     cd $build_dir
#     clean "$build_dir/gff3toolkit_conda_env"
#     echo "Download GFF3toolkit-v${GFF3TOOLKIT_VERSION}"
#     $miniconda3_dir/conda create -y -p $build_dir/gff3toolkit_conda_env python=3.7
#     source $miniconda3_dir/activate $build_dir/gff3toolkit_conda_env
#     cd gff3toolkit_conda_env
#     git clone $GFF3TOOLKIT_DOWNLOAD_URL
#     cd  GFF3toolkit
#     python setup.py install
#     #$miniconda3_dir/conda install -y pip
#     #pip3 install gff3tool==${GFF3TOOLKIT_VERSION}
#     source $miniconda3_dir/deactivate
#     note_installed $gff3toolkit_dir
# fi

# --------------- UCSC Utilities -----------------
ucsc_dir="$build_dir/UCSC_Utilities"
echo ""
echo "[$(timestamp)] Installing UCSC utilities ..."
if [ -z $(check_installed $ucsc_dir) ]; then
    cd $build_dir
    clean "$build_dir/UCSC_Utilities"
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

# --------------- EMBOSS ------------------
emboss_dir="$build_dir/EMBOSS-${EMBOSS_VERSION}/emboss"
echo ""
echo "[$(timestamp)] Installing emboss ..."
if [ -z $(check_installed $emboss_dir) ]; then
    cd $build_dir
    clean "$build_dir/EMBOSS-${EMBOSS_VERSION}"
    echo "Download EMBOSS"
    download $EMBOSS_DOWNLOAD_URL "emboss_v${EMBOSS_VERSION}.tar.gz"
    tar xvzf emboss_v${EMBOSS_VERSION}.tar.gz
    cd EMBOSS-${EMBOSS_VERSION}
    ./configure --without-x
    make -j $MAKE_JOBS
    cd $build_dir
    rm emboss_v${EMBOSS_VERSION}.tar.gz
    note_installed $emboss_dir
fi

# --------------- ERPIN ------------------
erpin_dir="$build_dir/erpin${ERPIN_VERSION}.serv/bin"
echo ""
echo "[$(timestamp)] Installing erpin ..."
if [ -z $(check_installed $erpin_dir) ]; then
    cd $build_dir
    clean "$build_dir/erpin${ERPIN_VERSION}.serv"
    echo "Download ERPIN"
    download $ERPIN_DOWNLOAD_URL "erpin_v${ERPIN_VERSION}.tar.gz"
    tar xvzf erpin_v${ERPIN_VERSION}.tar.gz
    cd erpin${ERPIN_VERSION}.serv
    make -j $MAKE_JOBS
    cd $build_dir
    rm erpin_v${ERPIN_VERSION}.tar.gz
    note_installed $erpin_dir
fi

# --------------- tbl2asn ------------------
tbl2asn_dir="$build_dir/tbl2asn_dir"
echo ""
echo "[$(timestamp)] Installing tbl2asn ..."
if [ -z $(check_installed $tbl2asn_dir) ]; then
    cd $build_dir
    clean "$build_dir/tbl2asn_dir"
    echo "Download tbl2asn"
    mkdir tbl2asn_dir
    cd tbl2asn_dir
    wget -c  $TBL2ASN_DOWNLOAD_URL # linux64.tbl2asn.gz
    mv linux64.tbl2asn.gz tbl2asn.gz
    gunzip tbl2asn.gz
    chmod 755 tbl2asn
    note_installed $tbl2asn_dir
fi

# --------------- PirObject ----------------
pirobject_dir="$build_dir/PirObject-${PIROBJECT_VERSION}"
echo ""
echo "[$(timestamp)] Installing PirObject ..."
if [ -z $(check_installed $pirobject_dir) ]; then
    cd $build_dir
    clean "$build_dir/PirObject-${PIROBJECT_VERSION}"
    echo "Download PirObject"
    download $PIROBJECT_DOWNLOAD_URL "pirobject_v${PIROBJECT_VERSION}.tar.gz"
    tar xvzf pirobject_v${PIROBJECT_VERSION}.tar.gz
    cd PirObject-${PIROBJECT_VERSION}
    ln -s ./lib/PirObject.pm .
    cd $build_dir
    rm pirobject_v${PIROBJECT_VERSION}.tar.gz
    note_installed $pirobject_dir
fi

# --------------- PirModels ------------------
pirmodels_dir="$pirobject_dir/PirModels"
echo ""
echo "[$(timestamp)] Installing PirModels ..."
if [ -z $(check_installed $pirmodels_dir) ]; then
    cd $build_dir
    clean "$pirobject_dir/PirModels"
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
echo ""
echo "[$(timestamp)] Installing Flip ..."
if [ -z $(check_installed $flip_dir) ]; then
    cd $build_dir
    clean "$build_dir/Flip"
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
echo ""
echo "[$(timestamp)] Installing Umac ..."
if [ -z $(check_installed $umac_dir) ]; then
    cd $build_dir
    clean "$build_dir/Umac"
    echo "Download Umac"
    git clone $UMAC_DOWNLOAD_URL
    cd Umac
    git checkout -f -q $UMAC_GITHUB_COMMIT_VERSION
    note_installed $umac_dir
fi


# --------------- HMMsearchWC ------------------
hmmsearchwc_dir="$build_dir/HMMsearchWC"
echo ""
echo "[$(timestamp)] Installing HMMsearchWC ..."
if [ -z $(check_installed $hmmsearchwc_dir) ]; then
    cd $build_dir
    clean "$build_dir/HMMsearchWC"
    echo "Download HMMsearchWC"
    git clone $HMMSEARCHWC_DOWNLOAD_URL
    cd HMMsearchWC
    git checkout -f -q $HMMSEARCHWC_GITHUB_COMMIT_VERSION
    note_installed $hmmsearchwc_dir
fi


# --------------- RNAfinder ------------------
rnafinder_dir="$build_dir/RNAfinder"
echo ""
echo "[$(timestamp)] Installing RNAfinder ..."
if [ -z $(check_installed $rnafinder_dir) ]; then
    cd $build_dir
    clean "$build_dir/RNAfinder"
    echo "Download RNAfinder"
    git clone $RNAFINDER_DOWNLOAD_URL
    cd RNAfinder
    git checkout -f -q $RNAFINDER_GITHUB_COMMIT_VERSION
    cp DOT_RNAfinder.cfg .RNAfinder.cfg
    note_installed $rnafinder_dir
fi


# --------------- Mf2sqn ------------------
mf2sqn_dir="$build_dir/Mf2sqn"
echo ""
echo "[$(timestamp)] Installing Mf2sqn ..."
if [ -z $(check_installed $mf2sqn_dir) ]; then
    cd $build_dir
    clean "$build_dir/Mf2sqn"
    echo "Download Mf2sqn"
    git clone $MF2SQN_DOWNLOAD_URL
    cd Mf2sqn
    git checkout -f -q $MF2SQN_GITHUB_COMMIT_VERSION
    cp qualifs.pl $build_dir/cpanm/perlmods/lib/perl5
    note_installed $mf2sqn_dir
fi

# --------------- grab-fasta ------------------
grab_fasta_dir="$build_dir/grab-fasta"
echo ""
echo "[$(timestamp)] Installing grab-fasta ..."
if [ -z $(check_installed $grab_fasta_dir) ]; then
    cd $build_dir
    clean "$build_dir/grab-fasta"
    echo "Download grab-fasta"
    git clone $GRAB_FASTA_DOWNLOAD_URL
    cd grab-fasta
    git checkout -f -q $GRAB_FASTA_GITHUB_COMMIT_VERSION
    note_installed $grab_fasta_dir
fi

# --------------- MFannot_data ------------------
mfannot_data_dir="$build_dir/MFannot_data"
echo ""
echo "[$(timestamp)] Installing MFannot_data ..."
if [ -z $(check_installed $mfannot_data_dir) ]; then
    cd $build_dir
    clean "$build_dir/MFannot_data"
    echo "Download MFannot_data"
    git clone $MFANNOT_DATA_DOWNLOAD_URL
    cd MFannot_data
    git checkout -f -q $MFANNOT_DATA_GITHUB_COMMIT_VERSION
    note_installed $mfannot_data_dir
fi

# --------------- MFannot ------------------
mfannot_dir="$build_dir/MFannot"
echo ""
echo "[$(timestamp)] Installing MFannot ..."
if [ -z $(check_installed $mfannot_dir) ]; then
    cd $build_dir
    clean "$build_dir/MFannot"
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
echo "export build_dir=${build_dir}" >> env.sh
echo "export PYTHONPATH=${PYTHONPATH}" >> env.sh
echo "export PERL5LIB=${PERL5LIB}" >> env.sh 
echo "export cpanm_dir=${cpanm_dir}" >> env.sh
if [[ $lite_installation == "no" ]]
then
    echo "export miniconda2_dir=${miniconda2_dir}" >> env.sh
    echo "export pbassembly_dir=${pbassembly_dir}" >> env.sh
    echo "export bax2bam_dir=${bax2bam_dir}" >> env.sh
fi
echo "export miniconda3_dir=${miniconda3_dir}" >> env.sh
# echo "export sra_dir=${sra_dir}" >> env.sh
echo "export porechop_dir=${porechop_dir}" >> env.sh
echo "export filtlong_dir=${filtlong_dir}" >> env.sh
echo "export minimap2_dir=${minimap2_dir}" >> env.sh
echo "export canu_dir=${canu_dir}" >> env.sh
# echo "export whatshap_dir=${whatshap_dir}" >> env.sh
echo "export flye_dir=${flye_dir}" >> env.sh
echo "export wtdbg2_dir=${wtdbg2_dir}" >> env.sh
echo "export smartdenovo_dir=${smartdenovo_dir}" >> env.sh
# echo "export raven_dir=${raven_dir}" >> env.sh
echo "export shasta_dir=${shasta_dir}" >> env.sh
echo "export guppy_cpu_dir=${guppy_cpu_dir}" >> env.sh
echo "export guppy_gpu_dir=${guppy_gpu_dir}" >> env.sh
# echo "export bonito_dir=${bonito_dir}" >> env.sh
# echo "export nanofilt_dir=${nanofilt_dir}" >> env.sh
echo "export nanoplot_dir=${nanoplot_dir}" >> env.sh
echo "export nanopolish_dir=${nanopolish_dir}" >> env.sh
echo "export parallel_dir=${parallel_dir}" >> env.sh
echo "export medaka_dir=${medaka_dir}" >> env.sh
echo "export racon_dir=${racon_dir}" >> env.sh
# echo "export marginpolish_dir=${marginpolish_dir}" >> env.sh
# echo "export helen_dir=${helen_dir}" >> env.sh
# echo "export pepper_dir=${pepper_dir}" >> env.sh
# echo "export homopolish_dir=${homopolish_dir}" >> env.sh
# echo "export quast_dir=${quast_dir}" >> env.sh
echo "export ragout_dir=${ragout_dir}" >> env.sh
echo "export ragtag_dir=${ragtag_dir}" >> env.sh
# echo "export hdf_dir=${hdf_dir}" >> env.sh
# echo "export h5prefix=${h5prefix}" >> env.sh
# echo "export hal_dir=${hal_dir}" >> env.sh
echo "export mummer4_dir=${mummer4_dir}" >> env.sh
echo "export gnuplot_dir=${gnuplot_dir}" >> env.sh
echo "export bedtools_dir=${bedtools_dir}" >> env.sh
echo "export spades_dir=${spades_dir}" >> env.sh
echo "export prodigal_dir=${prodigal_dir}" >> env.sh
echo "export cap3_dir=${cap3_dir}" >> env.sh
echo "export circlator_dir=${circlator_dir}" >> env.sh
echo "export trimmomatic_dir=${trimmomatic_dir}" >> env.sh
echo "export bwa_dir=${bwa_dir}" >> env.sh
echo "export samtools_dir=${samtools_dir}" >> env.sh
echo "export picard_dir=${picard_dir}" >> env.sh
#echo "export hapog_dir=${hapog_dir}" >> env.sh
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
echo "export diamond_dir=${diamond_dir}" >> env.sh
echo "export gatk3_dir=${gatk3_dir}" >> env.sh
#echo "export gff3toolkit_dir=${gff3toolkit_dir}" >> env.sh
echo "export ucsc_dir=${ucsc_dir}" >> env.sh

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


# test java configuration: requireds java 1.8 
echo ""
echo "##########################################"
echo "Testing java configuration ..."
echo ""
java_bin=""
if type -p java
then 
    java_bin=$(which java)
    echo "found java executable in PATH: $java_bin"
elif [[ -n "$JAVA_HOME" ]] && [[ -x "$JAVA_HOME/bin/java" ]]
then 
    java_bin="$JAVA_HOME/bin/java"
    echo "found java executable in JAVA_HOME: $java_bin" 
else 
    echo "";
    echo "!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!";
    echo "Failed to detect Java installation in the system!"
    echo "Please install java 1.8, which is a dependency of RecombineX!\n";
    echo "After the java installation, please manually set the directory path to java 1.8 executable on the last line of the env.sh file generated by this installation script!"
    echo "export java_dir=" >> env.sh
fi  

if [[ -n "$java_bin" ]]
then
    java_version=$("$java_bin" -version 2>&1 | awk -F '"' '/version/ {print $2}')
    echo "detected java_version: $java_version"
    if [ $(tidy_version "$java_version") -eq $(tidy_version "1.8") ]
    then
        java_dir=$(dirname $java_bin)
        echo "export java_dir=${java_dir}" >> env.sh
        echo "You have the correct java version for LRSDAY! LRSDAY will take care of the configuration."
    else
        echo "";
        echo "!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!";
        echo "Your java version is not the version required by LRSDAY (java v1.8)!"
        echo "Please manually set the directory path to java 1.8 executable on the last line of the env.sh file generated by this installation script!"
        echo "!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!";
        echo "export java_dir=" >> env.sh
    fi
fi



echo ""
echo "uncompress large supporting files ..."
if [[ -e $LRSDAY_HOME/data/Proteome_DB_for_annotation.CDhit_I95.fa.gz ]]
then
    gunzip $LRSDAY_HOME/data/Proteome_DB_for_annotation.CDhit_I95.fa.gz
fi

if [[ -e $LRSDAY_HOME/data/SGDref.PoFF.ffn.gz ]]
then
    gunzip $LRSDAY_HOME/data/SGDref.PoFF.ffn.gz
fi

if [[ -e $LRSDAY_HOME/data/SGDref.PoFF.faa.gz ]]
then
    gunzip $LRSDAY_HOME/data/SGDref.PoFF.faa.gz
fi

if [[ -e $LRSDAY_HOME/data/SGDref.PoFF.gff.gz ]]
then
    gunzip $LRSDAY_HOME/data/SGDref.PoFF.gff.gz
fi

if [[ -e $LRSDAY_HOME/data/te_proteins.fasta.gz ]]
then
    gunzip $LRSDAY_HOME/data/te_proteins.fasta.gz
fi

if [[ -e $LRSDAY_HOME/Example_Outputs/CPG_1a.nuclear_genome.tidy.fa.gz ]]
then
    gunzip $LRSDAY_HOME/Example_Outputs/CPG_1a.nuclear_genome.tidy.fa.gz
fi

if [[ -e $LRSDAY_HOME/Example_Outputs/CPG_1a.nuclear_genome.tidy.gff3.gz ]]
then
    gunzip $LRSDAY_HOME/Example_Outputs/CPG_1a.nuclear_genome.tidy.gff3.gz
fi

if [[ -e $LRSDAY_HOME/Example_Outputs/CPG_1a.nuclear_genome.tidy.cds.fa.gz ]]
then
    gunzip $LRSDAY_HOME/Example_Outputs/CPG_1a.nuclear_genome.tidy.cds.fa.gz
fi

if [[ -e $LRSDAY_HOME/Example_Outputs/CPG_1a.nuclear_genome.tidy.pep.fa.gz ]]
then
    gunzip $LRSDAY_HOME/Example_Outputs/CPG_1a.nuclear_genome.tidy.pep.fa.gz
fi

if [[ -e $LRSDAY_HOME/Example_Outputs/CPG_1a.mitochondrial_genome.tidy.fa.gz ]]
then
    gunzip $LRSDAY_HOME/Example_Outputs/CPG_1a.mitochondrial_genome.tidy.fa.gz
fi

if [[ -e $LRSDAY_HOME/Example_Outputs/CPG_1a.mitochondrial_genome.tidy.gff3.gz ]]
then
    gunzip $LRSDAY_HOME/Example_Outputs/CPG_1a.mitochondrial_genome.tidy.gff3.gz
fi

if [[ -e $LRSDAY_HOME/Example_Outputs/CPG_1a.mitochondrial_genome.tidy.cds.fa.gz ]]
then
    gunzip $LRSDAY_HOME/Example_Outputs/CPG_1a.mitochondrial_genome.tidy.cds.fa.gz
fi

if [[ -e $LRSDAY_HOME/Example_Outputs/CPG_1a.mitochondrial_genome.tidy.pep.fa.gz ]]
then
    gunzip $LRSDAY_HOME/Example_Outputs/CPG_1a.mitochondrial_genome.tidy.pep.fa.gz
fi

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
