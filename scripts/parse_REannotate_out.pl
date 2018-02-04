#!/usr/bin/perl
use warnings;
use strict;
use Getopt::Long;

##############################################################
#  script: parse_REannotate_out.pl
#  author: Jia-Xing Yue (GitHub ID: yjx1217)
#  last edited: 2017.06.17
#  description: parse REannotate tabular output and extract fasta sequences for annotated TEs
#  example: perl parse_REannotate_out.pl -r genome.fa(.gz) -i prefix.REannotate.out -p prefix
##############################################################

my ($refseq, $input, $prefix);

GetOptions('r|refseq:s' => \$refseq, # reference genome in fasta format, could be compressed by gzip (e.g. *.fa, *.fasta, *.fa.gz, *.fasta.gz)
	   'i|input:s' => \$input, # REannotate tabular output file
	   'p|prefix:s' => \$prefix); # file name prefix for the output file


my $refseq_fh = read_file($refseq);
my %refseq = parse_fasta_file($refseq_fh);

my $input_fh = read_file($input);
#my $output_complete_gff = "$prefix.complete.raw.gff";
#my $output_complete_gff_fh = write_file($output_complete_gff);
my $output_complete_seq = "$prefix.complete.raw.fa";
my $output_complete_seq_fh = write_file($output_complete_seq);
#my $output_truncated_gff = "$prefix.truncated.raw.gff";
#my $output_truncated_gff_fh = write_file($output_truncated_gff);
my $output_truncated_seq = "$prefix.truncated.raw.fa";
my $output_truncated_seq_fh = write_file($output_truncated_seq);
#my $output_soloLTR_gff = "$prefix.soloLTR.raw.gff";
#my $output_soloLTR_gff_fh = write_file($output_soloLTR_gff);
my $output_soloLTR_seq = "$prefix.soloLTR.raw.fa";
my $output_soloLTR_seq_fh = write_file($output_soloLTR_seq);


while (<$input_fh>) {
    chomp;
    /^#/ and next;
    /^\s*$/ and next;
    /^id/ and next;
    my @line = split /\s+/, $_;
    my $id = $line[0];
    my $chr = $line[1];
    my ($type) = ($id =~ /([a-z]+)\d+/);
    my $family = $line[2];
    my $start = $line[6];
    my $end = $line[9];
    my $strand = $line[20];
    # print "id=$id, chr=$chr, type=$type, family=$family, start=$start, end=$end, strand=$strand\n";
    my $seq = substr $refseq{$chr}, $start - 1, $end-$start+1;
    if ($strand eq 'C') {
	$strand = '-';
	$seq = revcom($seq);
    }
    
    if (($type eq 'u') or ($type eq 'i')) {
	# print $output_complete_gff_fh "$chr\tcomplete\t$family\t$start\t$end\t.\t$strand\t.\t$id:$family:$chr:$start-$end:$strand\n";
	print $output_complete_seq_fh ">$id:$family:$chr:$start-$end:$strand\n$seq\n";
    } elsif ($type eq 't') {
	# print $output_truncated_gff_fh "$chr\ttruncated\t$family\t$start\t$end\t.\t$strand\t.\t$id:$family:$chr:$start-$end:$strand\n";
	print $output_truncated_seq_fh ">$id:$family:$chr:$start-$end:$strand\n$seq\n";
    } else {
	# print $output_soloLTR_gff_fh "$chr\tsoloLTR\t$family\t$start\t$end\t.\t$strand\t.\t$id:$family:$chr:$start-$end:$strand\n";
	print $output_soloLTR_seq_fh ">$id:$family:$chr:$start-$end:$strand\n$seq\n";

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


sub parse_fasta_file {
    my $fh = shift @_;
    my %seq;
    my $seq_name = "";
    while (<$fh>) {
	chomp;
	if (/^\s*$/) {
	    next;
	} elsif (/^\s*\#/) {
	    next;
	} elsif (/^>(\S+)/) {
	    $seq_name = $1;
	    $seq{$seq_name} = "";
	} else {
	    $seq{$seq_name} .= $_;
	}
    }
    return %seq;
}

sub revcom {
    my $seq = shift @_;
    my $seq_revcom = reverse $seq;
    $seq_revcom =~ tr/ATGCNatgcn/TACGNtacgn/;
    return $seq_revcom;
}
