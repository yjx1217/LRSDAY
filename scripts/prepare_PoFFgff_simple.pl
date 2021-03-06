#!/usr/bin/perl
use warnings;
use strict;
use Getopt::Long;

##############################################################
#  script: prepare_PoFFgff_simple.pl
#  author: Jia-Xing Yue (GitHub ID: yjx1217)
#  last edited: 2017.06.17
#  description: prepare the gff3 file for the proteinortho run with PoFF option support
#  example: perl prepare_PoFFgff_simple.pl -i prefix.raw.gff3 -o prefix.PoFF.gff3 
##############################################################

my ($input, $output);

GetOptions('input|i:s' => \$input, # input raw gff3
           'output|o:s' => \$output); # output gff3 for proteinortho

my $input_fh = read_file($input);
my $output_fh = write_file($output);

my %seen = ();
while (<$input_fh>) {
    chomp;
    /^##FASTA/ and exit;
    /^#/ and next;
    /^\s*$/ and next;
    my ($chr, $source, $type, $start, $end, $score, $strand, $frame, $attributes) = split /\t/, $_;
    if ($type eq 'gene') {
	my ($id, $name) = ($attributes =~ /ID=([^;]+);\S*Name=([^;]+)/);
	if (not exists $seen{$id}) {
	    print $output_fh "$chr\t$source\tCDS\t$start\t$end\t$score\t$strand\t$frame\tID=$id;Name=$name\n";
	    $seen{$id} = 1;
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
    if ($file =~ /\.gz$/) {
        open($fh, "| gzip -c >$file") or die "can't open $file\n";
    } else {
        open($fh, ">$file") or die "can't open $file\n";
    }
    return $fh;
}  
