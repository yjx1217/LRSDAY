#!/usr/bin/perl
use warnings;
use strict;
use Getopt::Long;

##############################################################
#  script: fastq2id.pl
#  author: Jia-Xing Yue (GitHub ID: yjx1217)
#  last edited: 2018.12.28
#  description: convert uncompressed or compressed fastq files into the fasta format.
#  example: perl fastq2id.pl -i input.fastq(.gz) -o output.id.list.txt
##############################################################

my ($input, $output);

GetOptions('input|i:s' => \$input,
	   'output|o:s' => \$output);


my $input_fh = read_file($input);
my $count = 0;
my $output_fh = write_file($output);

while (<$input_fh>) {
    chomp;
    if ($. % 4 == 1) {
	my ($id) = ($_ =~ /^\@(\S+)/);
	print $output_fh "$id\n";
    } elsif ($. % 4 == 2) {
	# print $output_fh "$_\n";
	$count++;
    }
}

print "\nDone! Processed $count reads in total!\n\n";

sub read_file {
    my $file = shift @_;
    my $fh;
    if ($file =~ /\.gz$/) {
        open($fh, "gunzip -c $file |") or die "can't open pipe to $file";
    }
    else {
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


