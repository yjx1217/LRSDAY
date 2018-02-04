#!/usr/bin/perl
use warnings;
use strict;
use Getopt::Long;

##############################################################
#  script: filter_gff_by_feature_type.pl
#  author: Jia-Xing Yue (GitHub ID: yjx1217)
#  last edited: 2018.01.30
#  description: filter gff by selected feature type
#  example: perl filter_gff_by_feature_type.pl  -i input.gff3(.gz) -o output.gff3(.gz) -f gene 
##############################################################


my ($input, $output, $feature);
GetOptions('input|i:s' => \$input,
	   'output|o:s' => \$output,
	   'feature|f:s' => \$feature); 

my $input_fh = read_file($input);
my $output_fh = write_file($output);

while (<$input_fh>) {
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
    if ($type eq $feature) {
	print $output_fh "$_\n";
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



sub revcom {
    my $seq = shift @_;
    my $seq_revcom = reverse $seq;
    $seq_revcom =~ tr/ATGCatgc/TACGtacg/;
    return $seq_revcom;
}
