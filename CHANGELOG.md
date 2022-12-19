# Changelog

All notable changes to this project will be documented in this file.

The format is based on [Keep a Changelog](http://keepachangelog.com/en/1.0.0/)
and this project adheres to [Semantic Versioning](http://semver.org/spec/v2.0.0.html).

## [1.7.0] - 2022-12-31
### Changed
- Better robunsteness for the dependency installation script.
- Functional enhancement for the basecalling, demultiplexing, and quality summary of ONT raw fast5 reads.
- Version updates for multiple genome assemblers.
- Functional enhancement for the mitochondrial genome assembly.
- Functional enhancement for the mitochondrial gene annotation.
- Functional enhancement for the nuclear gene annotation.
- Functional enhancement for the TE annotation.
- Functional enhancement for SGD_based orthology identification and tagging.
- Removal of support for the long-read-assembler: Ra.
- Removal of support for the long-read-polisher: MarginPolish.

## [1.6.0] - 2019-10-03
### Added
- Supports for canu-triobinning assembly mode for phased hybrid genome assembly.
- Supports for new alternative long-read-assemblers: Ra and Shasta.
- Supports for new alternative long-read-based polishers: Racon, Medaka, and MarginPolish.
- Supports for new alternative reference-genome-based scaffolder: RaGOO.

### Changed
- Better robustness for the dependency installation script.
- Better treatment for identifying mitochondrial contigs when they contain assembly gaps.
- Minor adjustment for a few task-specific directory names.
- Software version updates for a number of dependencies.
### Fixed
- A bug that will disrupt long-read-based assembly run when assembler="canu-flye", "canu-wtdbg2" or "canu-smartdenovo" because of deleting a non-existing file.

## [1.5.0] - 2019-05-13
### Added
- Supports for native ONT nanopore basecalling, demultiplexing, and read profile plotting.
- An lite-weight bash script for generating dotplot for any pairs of fasta files.
### Changed
- Better robustness for the dependency installation script.
- Reducing the requirements of setting output prefix by assuming the same prefix was used across different modules.
- Setting adjustment for long-read filtering and downsampling.
- Applying read filtering/downsampling for the testing example.
- Software version updates for a number of dependencies.
### Fixed
- A bug that might lead to missing gene annotation in certain genomic region.
- Typos in the installation script.
- Typos in the manual.

## [1.4.0] - 2019-03-21
### Changed
- Supports for multi-round assembly polishing using both long and short reads.
- Setting adjustment for long-read filtering and downsampling.
- Software version updates for a number of dependencies. 
### Fixed
- Compatibility issues due to recent version updates of conda an bioconda.
- Typos in the installation script.

## [1.3.1] - 2019-01-22
### Added
- A script for generated demultiplexed fastq reads based on nanopore's guppy demultiplexing summary file.
- A script for extracting all read IDs from the input fastq file.
- A script for calculating the sequence length of each enclosed sequences from the input fasta file.
### Changed
- Software version updates for a number of dependencies. Importantly, nanopolish was bumped up to v0.11.0 to support the new multi-fast5 nanopore reads.
- Downloading URL updates for testing data due to changes on the EBI data server. 
- Adopting picard tools' new commandline syntax.
### Fixed
- A bug that prevents correct parameter parsing when multiple customized canu parameters are specified.
- Minor bugs related with deleting intermediate files in special cases.

## [1.3.0] - 2018-11-13
### Added
- Support for one more alternative assembler: wtdbg2.
### Changed
- Substantially more automated installation/setup process.
- Software version or downloading URL updates for a number of dependencies.
### Fixed
- Bugs introduced due to changes made for file/parameter names in the LRSDAY.01.Long-read-based_Genome_Assembly.sh script when using some alternative assemblers.
- Mismatched step numbers and file names in the manual due to previous version changes.
- Typos in the manual.

## [1.2.0] - 2018-10-15
### Added
- Support for adapter trimming for Nanopore reads (via Porechop).
- Support for long-read filtering based on both quality and length (via Filtlong).
- Support for long-read-based polishing for PacBio and Nanopore reads (via Quiver/Arrow for PacBio reads and nanopolish for Nanopore reads).
- Support for the bax2bam format conversion for the PacBio RSII reads to make it compatible with PacBio's current SMRT pipeline.
- Support for dedicated mitochondrial gene annotation (via Mfannot).
### Changed
- Treat nuclear genome and mitochondrial genome separately in the annotation phase.
- Better robustness for various processing scripts.
- Software version or downloading URL updates for a number of dependencies.
### Fixed
- Typos in the manual.

## [1.1.0] - 2018-07-11
### Added
- This change log file: Changelog.md
- Support for customized parameter settings for Canu.
- Support for alternative assemblers: Flye and smartdenovo.
- The script subsampling_seqeunces.pl for downsampling input reads.
- Support for the report of exact numbers of SNP and INDEL corrections made based on the Illumina data.
### Changed
- Better resource management for the Illumina-based polishing step.
- Better robustness for the reference-based scaffolding step.
- Better performance and robustness for the mitochondrial genome assembly improvement step.
- Better robustness for the TE annotation step.
- Better robustness for the reference-based gene orthology identification step.
- Software version or downloading URL updates for a number of dependencies.
- Preparation for future incorporation of long-read-based polishing step.
### Fixed
- A bug in the script identify_contigs_for_RefChr_by_mummer.pl that ignores user-specified options.
- Typos in the manual.

## [1.0.0] - 2018-02-04
### Added
- Initial version of LRSDAY.
