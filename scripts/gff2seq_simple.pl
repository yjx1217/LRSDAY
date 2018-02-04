#!/usr/bin/perl
use warnings;
use strict;
use Getopt::Long;

##############################################################
#  script: gff2seq_simple.pl
#  author: Jia-Xing Yue (GitHub ID: yjx1217)
#  last edited: 2017.06.17
#  description: extract feature sequences from refrence genome based GFF3 annotation 
#  example: perl gff2seq_simple.pl -r genome.fa(.gz) -g feature.gff3 -o feature.fa(.gz)
##############################################################


my ($refseq_file, $gff_file, $output);
GetOptions('refseq|r:s' => \$refseq_file, # reference genome file in fasta format with or without compression (e.g. *.fa or *.fasta or *.fa.gz ...)
	   'gff|g:s' => \$gff_file, # feature annotation file in gff3 format
	   'output|o:s' => \$output); # feature sequence output file in fasta format 

my $refseq_fh = read_file($refseq_file);
my %refseq = parse_fasta_file($refseq_fh);
my $gff_fh = read_file($gff_file);
my %gff = parse_simple_gff_file($gff_fh);
my %output = ();
my $output_fh = write_file($output);

foreach my $id (sort keys %gff) {
    # print "gene_id = $gene_id\n";
    my ($chr, $start, $end, $strand) = ($gff{$id} =~ /(\S+):(\d+)-(\d+):(\+|\-)/);
    if (not exists $refseq{$chr}) {
	print "error! unrecognized chromosome: $chr\n";
	exit;
    } else {
	my $seq = substr $refseq{$chr}, $start - 1, $end - $start + 1;
	if ($strand eq '-') {
	    $seq = revcom($seq);
	}
	$output{$id} = $seq;
	print $output_fh ">$id\n$output{$id}\n"; 
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
    if ($file =~ /\.gz$/) {
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
    my $flag = 0;
    while (<$fh>) {
	chomp;
	if (/^\s*$/) {
	    next;
	} elsif (/^\s*#/) {
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


sub parse_simple_gff_file {
    my $fh = shift @_;
    my %gff = ();
    while (<$fh>) {
	chomp;
	(/^\s*\#/) and next;
	(/^\#/) and next;
	my @line = split /\t/, $_;
	my $chr = $line[0];
	my $type = $line[2];
	my $start = $line[3];
	my $end = $line[4];
	my $strand = $line[6];
	my $phase = $line[7];
	my $annotation = $line[8];
	# print "annotation = $annotation\n";
	my ($id, $name);
	($id, $name) = ($annotation =~ /ID\=([^;]+);Name\=([^;]+)/);
	$gff{$id} = "$chr:$start-$end:$strand";
	}
    return %gff;
}


sub revcom {
    my $seq = shift @_;
    my $seq_revcom = reverse $seq;
    $seq_revcom =~ tr/ATGCatgc/TACGtacg/;
    return $seq_revcom;
}
