#!/usr/bin/perl
use warnings;
use strict;
use Getopt::Long;

##############################################################
#  script: mfannot2gff3.pl
#  author: Jia-Xing Yue (GitHub ID: yjx1217)
#  last edited: 2018.09.01
#  description: converting the mfannot output to the GFF3 format 
#  example: perl mfannot2gff3.pl -i genome.mfannot.out -g genome.mfannot.gff3 -m lite 
##############################################################

my ($input, $output, $mode);
$mode = "lite"; # use lite to filter out gene models annotated within introns
GetOptions('input|i:s' => \$input,
           'output|o:s' => \$output,
           'mode|m:s' => \$mode); # "full" or "lite"

my $input_fh = read_file($input);

my %chr = ();
my @chr = ();
my %features = ();
my @features = ();
my $feature = "";
my %gene2exon_intron = ();

while (<$input_fh>) {
    chomp;
    /^\s*$/ and next;
    /^#/ and next;
    /^;;/ and next;
    if (/^>(.*) gc=(\d+)/) {
	my $chr = $1;
	$chr{$chr}{'genetic_code'} = $2;
	print "chr=$chr, genetic_code=$chr{$chr}{'genetic_code'}\n";
	push @chr, $chr;
    } elsif (/^;\s+G-(\S+)\s+(\S+)\s+(start|end)\s*(\S?\.*)/) {
	$feature = $1;
	my $feature_type;
	my $feature_strand = $2;
	my $feature_mark = $3;
	my $feature_note = $4;
	if ($feature =~ /(\S+)-E(\d+)$/) {
	    $gene2exon_intron{$1}{'exon'}{$2} = $feature;
	    $feature_type = "exon";
	} elsif ($feature =~ /(\S+)-I(\d+)$/) {
	    $gene2exon_intron{$1}{'intron'}{$2} = $feature;
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
    foreach my $feature (@features) {
        print "output feature=$feature\n";
	if ($features{$feature}{'type'} eq "gene") {
	    if ($feature =~ /(rnl|rns|rps)/) {
		print $output_fh "$chr\tmfannot\trRNA\t$features{$feature}{'start'}\t$features{$feature}{'end'}\t.\t$features{$feature}{'strand'}\t.\tID=${feature};Name=${feature}\n";
	    } else {
		if ($mode eq "lite") {
		    if ($feature =~ /-I(\d+)-/) {
			next;
		    }
		}
		if (not exists $gene2exon_intron{$feature}) {
		    # single exon gene
		    print $output_fh "$chr\tmfannot\tgene\t$features{$feature}{'start'}\t$features{$feature}{'end'}\t.\t$features{$feature}{'strand'}\t.\tID=$feature;Name=$feature\n";
		    print $output_fh "$chr\tmfannot\tmRNA\t$features{$feature}{'start'}\t$features{$feature}{'end'}\t.\t$features{$feature}{'strand'}\t.\tID=${feature}.mRNA.1;Name=${feature}.mRNA.1;Parent=$feature\n";
		    # add back exon/CDS annotation
		    my $gene_id = $feature;
		    my $exon_index = 1;
		    print $output_fh "$chr\tmfannot\texon\t$features{$feature}{'start'}\t$features{$feature}{'end'}\t.\t$features{$feature}{'strand'}\t.\tID=${gene_id}.exon.$exon_index;Name=${gene_id}.exon.$exon_index;Parent=${gene_id}.mRNA.1\n";
		    print $output_fh "$chr\tmfannot\tCDS\t$features{$feature}{'start'}\t$features{$feature}{'end'}\t.\t$features{$feature}{'strand'}\t.\tID=${gene_id}.CDS.$exon_index;Name=${gene_id}.CDS.$exon_index;Parent=${gene_id}.mRNA.1\n";
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
		    print $output_fh "$chr\tmfannot\tgene\t$features{$feature}{'start'}\t$features{$feature}{'end'}\t.\t$features{$feature}{'strand'}\t.\tID=$feature;Name=$feature\n";
		    print $output_fh "$chr\tmfannot\tmRNA\t$features{$feature}{'start'}\t$features{$feature}{'end'}\t.\t$features{$feature}{'strand'}\t.\tID=${feature}.mRNA.1;Name=${feature}.mRNA.1;Parent=$feature\n";
		}
	    }
	} elsif ($features{$feature}{'type'} eq "exon") {
	    my ($gene_id, $exon_index) = ($feature =~ /(\S+)-E(\d+)/);
	    print $output_fh "$chr\tmfannot\texon\t$features{$feature}{'start'}\t$features{$feature}{'end'}\t.\t$features{$feature}{'strand'}\t.\tID=${gene_id}.exon.$exon_index;Name=${gene_id}.exon.$exon_index;Parent=${gene_id}.mRNA.1\n";
	    print $output_fh "$chr\tmfannot\tCDS\t$features{$feature}{'start'}\t$features{$feature}{'end'}\t.\t$features{$feature}{'strand'}\t.\tID=${gene_id}.CDS.$exon_index;Name=${gene_id}.CDS.$exon_index;Parent=${gene_id}.mRNA.1\n";
	}	
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
