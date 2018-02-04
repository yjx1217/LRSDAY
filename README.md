# LRSDAY
**LRSDAY: Long-read Sequencing Data Analysis for Yeasts**

A highly transparent, automated and powerful computational framework for high-quality genome assembly and annotation.


## Description
Long-read sequencing technologies have become increasingly popular in genome projects due to their strengths in resolving complex genomic regions. As a leading model organism with small genome size and great biotechnological importance, the budding yeast, *Saccharomyces cerevisiae*, has many isolates currently being sequenced with long reads. However, analyzing long-read sequencing data to produce high-quality genome assembly and annotation remains challenging. Here we present LRSDAY, the first one-stop solution to streamline this process. LRSDAY can produce chromosome-level end-to-end genome assembly and comprehensive annotations for various genomic features (including centromeres, protein-coding genes, tRNAs, transposable elements and telomere-associated elements) that are ready for downstream analysis. Although tailored for *S. cerevisiae*, we designed LRSDAY to be highly modular and customizable, making it adaptable for virtually any eukaryotic organisms. Applying LRSDAY to a *S. cerevisiae* strain takes ~41 hrs to generate a complete and well-annotated genome from ~100X Pacific Biosciences (PacBio) reads using 4 threads.


## Citations
Jia-Xing Yue & Gianni Liti. (2018) Long-read sequencing data analysis for yeasts (LRSDAY). *Nature Protocols* (accepted pending minor revision; preprint available at http://www.biorxiv.org/content/early/2017/09/05/184572)

Jia-Xing Yue, Jing Li, Louise Aigrain, Johan Hallin, Karl Persson, Karen Oliver, Anders Bergström, Paul Coupland, Jonas Warringer, Marco Cosentino Lagomarsino, Gilles Fischer, Richard Durbin, Gianni Liti. (2017) Contrasting evolutionary genome dynamics between domesticated and wild yeasts. *Nature Genetics*, 49:913-924.

## Current version
* v1.0.0 Released on 2018/02/04

## License
LRSDAY itself is distributed under the MIT license. A number of LRSDAY's dependencies (e.g. CAP3, MAKER, GATK, blat, RepBase, etc) are under more restricted licenses, for which commerical use of the software needs to be discussed with the corresponding developers.


## Requirements
### Hardware, operating system and network
This protocol is designed for a desktop or computing server running an x86-64-bit Linux operating system. Multithreaded processors are preferred to speed up the process since many steps can be configured to use multiple threads in parallel. For assembling and analyzing the budding yeast genomes (genome size = ~12 Mb), at least 16 Gb of RAM and 100 Gb of free disk space are reccomended. When adapted for other eukaryotic organisms with larger genome sizes, the RAM and disk space consumption will scale up, majorly during de novo genome assembly (performed by [Canu](https://github.com/marbl/canu). Plese refer to [Canu’s manual](http://canu.readthedocs.io/en/latest/) for suggested RAM and disk space consumption for assembling large genomes. Stable Internet connection is required for the installation and configuration of LRSDAY as well as for retrieving the test data.


### Software or library requirements
* Bash (https://www.gnu.org/software/bash/)
* Bzip2 (http://www.bzip.org/)
* Cmake (https://cmake.org/)
* GCC and G++ v4.7 or newer (https://gcc.gnu.org/)
* Ghostscript (https://www.ghostscript.com)
* Git (https://git-scm.com/)
* GNU make (https://www.gnu.org/software/make/)
* Gzip (https://www.gnu.org/software/gzip/)
* Java runtime environment (JRE) v1.8.0 or newer (https://www.java.com)
* Perl v5.12 or newer (https://www.perl.org/)
* Python v2.7.9 or newer (https://www.python.org/)
* Python v3.6 or newer (https://www.python.org/)
* Tar (https://www.gnu.org/software/tar/)
* Unzip (http://www.info-zip.org/UnZip.html)
* Virtualenv v15.1.0 or newer (https://virtualenv.pypa.io)
* Wget (https://www.gnu.org/software/wget/)
* Zlib (https://zlib.net/)


