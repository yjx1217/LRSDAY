# Changelog

All notable changes to this project will be documented in this file.

The format is based on [Keep a Changelog](http://keepachangelog.com/en/1.0.0/)
and this project adheres to [Semantic Versioning](http://semver.org/spec/v2.0.0.html).

## [Unreleased]

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
