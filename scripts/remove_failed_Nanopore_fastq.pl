#!/usr/bin/perl
use warnings;
use strict;
use Getopt::Long;

##############################################################
#  script: filter_failed_Nanopore_fastq.pl
#  author: Jia-Xing Yue (GitHub ID: yjx1217)
#  last edited: 2017.10.12
#  description: convert uncompressed or compressed fastq files into the fasta format.
#  example: perl filter_failed_ONT_fastq.pl -i input.fastq(.gz) -o output.fastq(.gz)
##############################################################

my ($input, $output);

GetOptions('input|i:s' => \$input,
	   'output|o:s' => \$output);


my $input_fh = read_file($input);
my $output_fh = write_file($output);
my $flag = 0;
my $total_count = 0;
my $pass_count = 0;

while (<$input_fh>) {
    chomp;
    if ($. % 4 == 1) {
	$flag = 0;
	$total_count++;
	my ($id) = ($_ =~ /^\@(.*)/);
	if ($id =~ /_fail_/) {
	    $flag = 1;
	} else {
	    $pass_count++;
	    print $output_fh "\@$id\n";
	}
    } else {
	if ($flag == 0) {
	    print $output_fh "$_\n";
	}
    }
}

print "\nDone! Processed $total_count reads in total, $pass_count of them passed!\n\n";

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


