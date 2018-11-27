#!/usr/bin/perl
use warnings;
use strict;
use Getopt::Long;

##############################################################
#  script: guppy_reads_classifier.pl
#  author: Jia-Xing Yue (GitHub ID: yjx1217)
#  last edited: 2018.12.28
#  description: classify guppy basecalled fastq reads based on the guppy debarcoding summary file 
#  example: perl guppy_reads_classifier.pl -fq basecalled_reads.fastq(.gz) -s barcoding_summary.txt 
##############################################################

my ($fq, $summary);
$fq = "./../Basecalling_Guppy_out/basecalled_reads.fastq.gz";
$summary = "barcoding_summary.txt";

GetOptions('fastq|fq:s' => \$fq,
	   'summary|s:s' => \$summary);

my $summary_fh = read_file($summary);
my %summary = parse_summary_file($summary_fh);

my $fq_fh = read_file($fq);

my %barcodes = ();
my %output = ();
my $count = 0;
my $b;
my $file;
my $fh;
while (<$fq_fh>) {
    chomp;
    if ($. % 4 == 1) {
	my ($id) = ($_ =~ /^\@(\S+)/);
	# print "id = $id\n";
	$count++;
	if (exists $summary{$id}) {
	    $b = $summary{$id};
	    # print "barcode=$b\n";
	    $file = "$b.basecalled_reads.fastq";
	    open($fh, ">>$file");
	} else {
	    die "unexpected barcode: $b\n";
	}
    }
    print $fh "$_\n";
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


sub parse_summary_file {
    my $fh = shift @_;
    my %summary = ();
    while (<$fh>) {
	chomp;
	/^\s*$/ and next;
	/^#/ and next;
	/^read_id\tbarcode_arrangement/ and next;
	my @line = split /\t/, $_;
	my $read_id = $line[0];
	my $barcode_arrangement = $line[1];
	$summary{$read_id} = $barcode_arrangement;
    }
    return %summary;
}





