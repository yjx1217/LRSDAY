#!/usr/bin/perl
use warnings;
use strict;
use Getopt::Long;

##############################################################
#  script: mfannot2gff3.pl
#  author: Jia-Xing Yue (GitHub ID: yjx1217)
#  last edited: 2022.12.09
#  description: converting the mfannot and trnascan output files into a final GFF3 output 
#  example: perl mfannot2gff3.pl -tag strain_id -mfannot_out genome.mfannot.out.txt -trnascan_out genome.trnascan.out.txt -o genome.mfannot.gff3 -m lite 
##############################################################

my ($mfannot_out, $trnascan_out, $tag, $output, $mode);
$mode = "lite"; # use lite to filter out gene models annotated within introns
GetOptions('mfannot_out|mfannot_out:s' => \$mfannot_out,
	   'trnascan_out|trnascan_out:s' => \$trnascan_out,
	   'tag|t:s' => \$tag,
           'output|o:s' => \$output,
           'mode|m:s' => \$mode); # "full" or "lite"

my $trnascan_out_fh = read_file($trnascan_out);
my %trnascan_out = parse_trnascan_out_file($trnascan_out_fh);

my $mfannot_out_fh = read_file($mfannot_out);
my %chr = ();
my @chr = ();
my %features = ();
my @features = ();
my $feature = "";
my %gene2exon_intron = ();

while (<$mfannot_out_fh>) {
    chomp;
    /^\s*$/ and next;
    /^#/ and next;
    /^;;/ and next;
    if (/^>(.*) gc=(\d+)/) {
	my $chr = $1;
	$chr{$chr}{'genetic_code'} = $2;
	print "chr=$chr, genetic_code=$chr{$chr}{'genetic_code'}\n";
	push @chr, $chr;
    } elsif (/^;\s+G-orf/) {
	next;
    } elsif (/^;\s+G-(\S+)\s+(\S+)\s+(start|end)\s*(\S?\.*)/) {
	$feature = $1;
	my $feature_type;
	my $feature_strand = $2;
	my $feature_mark = $3;
	my $feature_note = $4;
	if ($feature =~ /(\S+)-E(\d+)$/) {
	    $gene2exon_intron{$1}{'exon'}{$2}{'feature'} = $feature;
	    $feature_type = "exon";
	} elsif ($feature =~ /(\S+)-I(\d+)$/) {
	    $gene2exon_intron{$1}{'intron'}{$2}{'feature'} = $feature;
            $feature_type = "intron";
	} else {
	    $feature_type = "gene";
	}
	if ($feature_strand eq "==>") {
	    $feature_strand = "+";
	} else {
	    $feature_strand = "-";
	}

	print "feature=$feature, feature_type=$feature_type, strand=$feature_strand, mark=$feature_mark, note=$feature_note\n";
	if (not exists $features{$feature}{'type'}) {
	    $features{$feature}{'type'} = $feature_type;
	    push @features, $feature;
	    # print "add new feature: $feature\n";
	}
	if (not exists $features{$feature}{'strand'}) {
	    $features{$feature}{'strand'} = $feature_strand;
	}
	if (not exists $features{$feature}{'note'}) {
	    @{$features{$feature}{'note'}} = ($feature_note);
	} else {
	    push @{$features{$feature}{'note'}}, $feature_note;
	}
	if (not exists $features{$feature}{'seq'}) {
            $features{$feature}{'seq'} = "";
	}
	if (not exists $features{$feature}{'start'}) {
            $features{$feature}{'start'} = "";
        }
	if (not exists $features{$feature}{'end'}) {
            $features{$feature}{'end'} = "";
        }	
	$features{$feature}{'mark'} = $feature_mark;
	my $flag = 0;
	if (($features{$feature}{'mark'} eq "end" ) and ($features{$feature}{'strand'} eq "+")) {
	    $flag = 1;
	} 
	if (($features{$feature}{'mark'} eq "start" ) and ($features{$feature}{'strand'} eq "-")) {
	    $flag = 1;
	} 
	if ($flag == 1) {
	    $feature = "";
	}
    } elsif (/^\s+(\d+)\s+([ATCGatcgNn]+)/) {
	my ($pos_start, $seq) = ($1, $2);
	my $pos_end = length($seq) + $pos_start - 1;
	if ($feature ne "") {
	    if (($features{$feature}{'mark'} eq "start") and ($features{$feature}{'strand'} eq "+")) {
		if ($features{$feature}{'start'} eq "") {
		    $features{$feature}{'start'} = $pos_start;
		}
		$features{$feature}{'seq'} = $features{$feature}{'seq'} . $seq;
		$features{$feature}{'end'} = $pos_end;
	    } 
	    if (($features{$feature}{'mark'} eq "end") and ($features{$feature}{'strand'} eq "-")) {
		if ($features{$feature}{'end'} eq "") {
		    $features{$feature}{'end'} = $pos_start;
		}
		$features{$feature}{'seq'} = revcom($seq) . $features{$feature}{'seq'};
		$features{$feature}{'start'} = $pos_end;
	    }
	}
    }
}   


foreach my $f (sort keys %trnascan_out) {
    print "trn-$f\n";
    my ($chr, $start_end, $strand) = split /:/, $f;
    my ($start, $end) = split /-/, $start_end;
    $features{"trn-$f"}{'type'} = "gene";
    $features{"trn-$f"}{'start'} = $start;
    $features{"trn-$f"}{'end'} = $end;
    $features{"trn-$f"}{'strand'} = $strand;
}




my $output_fh = write_file($output);
print $output_fh "##gff-version 3\n";

foreach my $chr (@chr) {
    print "output chr=$chr\n";
    print $output_fh "##chr=$chr, genetic_code=$chr{$chr}{'genetic_code'}\n";
    # adjustment before output
    foreach my $feature (@features) {
	my $feature_note = join ";", @{$features{$feature}{'note'}};
	if ($features{$feature}{'strand'} eq "-") {
	    ($features{$feature}{'start'}, $features{$feature}{'end'}) = ($features{$feature}{'end'}, $features{$feature}{'start'});
	}
    }
    # begin to output
    my %gene2cds_phase = ();
    foreach my $feature (@features) {
        print "output feature=$feature\n";
	if ($features{$feature}{'type'} eq "gene") {
	    if ($feature =~ /^(rnl|rns|rps)/) {
		print $output_fh "$chr\t$tag\trRNA\t$features{$feature}{'start'}\t$features{$feature}{'end'}\t.\t$features{$feature}{'strand'}\t.\tID=${tag}_${feature};Name=${tag}_${feature}\n";
	    } elsif ($feature =~ /^rnpB/) {
		print $output_fh "$chr\t$tag\trnpB\t$features{$feature}{'start'}\t$features{$feature}{'end'}\t.\t$features{$feature}{'strand'}\t.\tID=${tag}_RPM1;Name=${tag}_RPM1\n";
	    } elsif ($feature =~ /^trn/) {
		if (exists $trnascan_out{"$chr:$features{$feature}{'start'}-$features{$feature}{'end'}:$features{$feature}{'strand'}"}) {
		    my $id = "${tag}_tRNA:$chr:$features{$feature}{'start'}-$features{$feature}{'end'}:$features{$feature}{'strand'}";
		    my $name = $trnascan_out{"$chr:$features{$feature}{'start'}-$features{$feature}{'end'}:$features{$feature}{'strand'}"}{"name"};
		    $trnascan_out{"$chr:$features{$feature}{'start'}-$features{$feature}{'end'}:$features{$feature}{'strand'}"}{"flag"} = 1;
		    print $output_fh "$chr\t$tag\ttRNA\t$features{$feature}{'start'}\t$features{$feature}{'end'}\t.\t$features{$feature}{'strand'}\t.\tID=$id;Name=$name\n";
		} else {
		    print "skip feature: $feature!\n";
		}
	    } else {
		if ($mode eq "lite") {
		    if ($feature =~ /-I(\d+)-/) {
			next;
		    }
		    if ($feature =~ /^orf/) {
			next;
		    }
		}
		my $gene_name = uc $feature;
		$gene_name = $tag . "_" . $gene_name;
		if (not exists $gene2exon_intron{$feature}) {
		    
		    # single exon gene
		    print $output_fh "$chr\t$tag\tgene\t$features{$feature}{'start'}\t$features{$feature}{'end'}\t.\t$features{$feature}{'strand'}\t.\tID=$gene_name;Name=$gene_name\n";
		    print $output_fh "$chr\t$tag\tmRNA\t$features{$feature}{'start'}\t$features{$feature}{'end'}\t.\t$features{$feature}{'strand'}\t.\tID=${gene_name}.mRNA.1;Name=${gene_name}.mRNA.1;Parent=$gene_name\n";
		    # add back exon/CDS annotation
		    my $gene_id = $feature;
		    my $exon_index = 1;
 		    $gene2cds_phase{$feature}{'phase'} = 0;
		    $gene2cds_phase{$feature}{'cds_walk'} = $features{$feature}{'end'} - $features{$feature}{'start'} + 1;

		    print $output_fh "$chr\t$tag\texon\t$features{$feature}{'start'}\t$features{$feature}{'end'}\t.\t$features{$feature}{'strand'}\t.\tID=${gene_name}.exon.$exon_index;Name=${gene_name}.exon.$exon_index;Parent=${gene_name}.mRNA.1\n";
		    print $output_fh "$chr\t$tag\tCDS\t$features{$feature}{'start'}\t$features{$feature}{'end'}\t.\t$features{$feature}{'strand'}\t0\tID=${gene_name}.CDS.$exon_index;Name=${gene_name}.CDS.$exon_index;Parent=${gene_name}.mRNA.1\n";
		} else {
		    # multi-exon gene
		    # add back gene annotation
		    my $gene_start;
		    my $gene_end;
		    my @exon_index = sort {$a <=> $b} keys %{$gene2exon_intron{$feature}{'exon'}};
		    my $last_exon_index = $exon_index[-1];
		    $gene_start = $features{"${feature}-E1"}{'start'};
		    $gene_end = $features{"${feature}-E${last_exon_index}"}{'end'};
		    if ($features{$feature}{'strand'} eq "+") {
			$features{$feature}{'start'} = $gene_start;
			$features{$feature}{'end'} = $gene_end;
		    } else {
			$features{$feature}{'start'} = $gene_end;
                        $features{$feature}{'end'} = $gene_start;
		    }
		    $gene2cds_phase{$feature}{'phase'} = 0;
		    $gene2cds_phase{$feature}{'cds_walk'} = 0;

		    print $output_fh "$chr\t$tag\tgene\t$features{$feature}{'start'}\t$features{$feature}{'end'}\t.\t$features{$feature}{'strand'}\t.\tID=$gene_name;Name=$gene_name\n";
		    print $output_fh "$chr\t$tag\tmRNA\t$features{$feature}{'start'}\t$features{$feature}{'end'}\t.\t$features{$feature}{'strand'}\t.\tID=${gene_name}.mRNA.1;Name=${gene_name}.mRNA.1;Parent=$gene_name\n";
		}
	    }
	} elsif ($features{$feature}{'type'} eq "exon") {
	    my ($gene_id, $exon_index) = ($feature =~ /(\S+)-E(\d+)/);
	    my $gene_name = uc $gene_id;
	    $gene_name = $tag . "_" . $gene_name;
	    print $output_fh "$chr\t$tag\texon\t$features{$feature}{'start'}\t$features{$feature}{'end'}\t.\t$features{$feature}{'strand'}\t.\tID=${gene_name}.exon.$exon_index;Name=${gene_name}.exon.$exon_index;Parent=${gene_name}.mRNA.1\n";
	    if ($exon_index == 1) {
		$gene2cds_phase{$gene_id}{'phase'} = 0;
		$gene2cds_phase{$gene_id}{'cds_walk'} = $features{$feature}{'end'} - $features{$feature}{'start'} + 1;
		print $output_fh "$chr\t$tag\tCDS\t$features{$feature}{'start'}\t$features{$feature}{'end'}\t.\t$features{$feature}{'strand'}\t0\tID=${gene_name}.CDS.$exon_index;Name=${gene_name}.CDS.$exon_index;Parent=${gene_name}.mRNA.1\n";
	    } else {
		if ($gene2cds_phase{$gene_id}{'cds_walk'} % 3 != 0) {
		    $gene2cds_phase{$gene_id}{'phase'} = 3 - $gene2cds_phase{$gene_id}{'cds_walk'} % 3;
		} else {
		    $gene2cds_phase{$gene_id}{'phase'} = 0;
		}
		$gene2cds_phase{$gene_id}{'cds_walk'} = $gene2cds_phase{$gene_id}{'cds_walk'} + $features{$feature}{'end'} - $features{$feature}{'start'} + 1;
		print $output_fh "$chr\t$tag\tCDS\t$features{$feature}{'start'}\t$features{$feature}{'end'}\t.\t$features{$feature}{'strand'}\t$gene2cds_phase{$gene_id}{'phase'}\tID=${gene_name}.CDS.$exon_index;Name=${gene_name}.CDS.$exon_index;Parent=${gene_name}.mRNA.1\n";
	    }
	}
    }
}
foreach my $f (sort keys %trnascan_out) {
    if ($trnascan_out{$f}{"flag"} == 0) {
	my ($chr, $start_end, $strand) = split /:/, $f;
	my ($start, $end) = split /-/, $start_end;
	my $name = $trnascan_out{$f}{"name"};
	my $id = "${tag}_tRNA:$chr:$start_end:$strand";
	print $output_fh "$chr\t$tag\ttRNA\t$start\t$end\t.\t$strand\t.\tID=$id;Name=$name\n";
    }
}




sub read_file {
    my $file = shift @_;
    my $fh;
    if ($file =~ /\.gz$/) {
        open($fh, "gunzip -c $file |") or die "can't open pipe to $file";
    } else {
        open($fh, $file) or die "can't open $file";
    }
    return $fh;
}

sub write_file {
    my $file = shift @_;
    my $fh;
    if ($file =~ /.gz$/) {
        open($fh, "| gzip -c >$file") or die "can't open $file\n";
    } else {
        open($fh, ">$file") or die "can't open $file\n";
    }
    return $fh;
}  

sub revcom {
    my $seq = shift @_;
    my $seq_revcom = reverse $seq;
    $seq_revcom =~ tr/ATGCNatgcn/TACGNtacgn/;
    return $seq_revcom;
}


sub parse_trnascan_out_file {
    my $fh = shift @_;
    my %trnascan = ();
    while (<$fh>) {
        chomp;
        /^Sequence/ and next;
        /^Name/ and next;
        /^------/ and next;
        my ($chr, $id, $start, $end, $trna_type, $anti_codon, $intron_start, $intron_end, $cov_score) = split /\s+/, $_;
        my $strand;
        if ($start < $end) {
            $strand = "+";
        } else {
            ($start, $end) = ($end, $start);
            $strand = "-";
        }
        $trnascan{"$chr:$start-$end:$strand"}{'name'} = "tRNA_${trna_type}(${anti_codon})";
	$trnascan{"$chr:$start-$end:$strand"}{'flag'} = 0;
    }
    return %trnascan;
}
