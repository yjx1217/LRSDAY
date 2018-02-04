#-----Location of Executables Used by MAKER/EVALUATOR
makeblastdb=$blast_dir/makeblastdb #location of NCBI+ makeblastdb executable
blastn=$blast_dir/blastn #location of NCBI+ blastn executable
blastx=$blast_dir/blastx #location of NCBI+ blastx executable
tblastx=$blast_dir/tblastx #location of NCBI+ tblastx executable
formatdb= #location of NCBI formatdb executable
blastall= #location of NCBI blastall executable
xdformat= #location of WUBLAST xdformat executable
blasta= #location of WUBLAST blasta executable
prerapsearch=$rapsearch_dir/prerapsearch #location of prerapsearch executable
rapsearch=$rapsearch_dir/rapsearch #location of rapsearch executable
RepeatMasker=$repeatmasker_dir/RepeatMasker #location of RepeatMasker executable
exonerate=$exonerate_dir/exonerate #location of exonerate executable

#-----Ab-initio Gene Prediction Algorithms
snap=$snap_dir/snap #location of snap executable
gmhmme3= #location of eukaryotic genemark executable
gmhmmp= #location of prokaryotic genemark executable
augustus=$augustus_dir/augustus #location of augustus executable
fgenesh= #location of fgenesh executable
evm=$evm_dir/evidence_modeler.pl #location of EvidenceModeler executable
tRNAscan-SE=$trnascan_dir/tRNAscan-SE #location of trnascan executable
snoscan=$snoscan_dir/snoscan #location of snoscan executable

#-----Other Algorithms
probuild= #location of probuild executable (required for genemark)
