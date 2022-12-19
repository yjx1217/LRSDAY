#!/usr/bin/perl
use warnings;
use strict;
use Getopt::Long;

##############################################################
#  script: parse_REannotate_gff.pl
#  author: Jia-Xing Yue (GitHub ID: yjx1217)
#  last edited: 2017.07.06
#  description: reformat REannotate GFF output files
#  example: perl parse_REannotate_gff.pl  -o output
##############################################################

my $output;

GetOptions('output|o:s' => \$output); # file name for the output file


my @inputs = glob "*.REannotate.gff";
my $output_fh = write_file($output);


foreach my $input (@inputs) {
    my ($chr) = ($input =~ /(\S+)\.REannotate.gff/);
    if ($chr !~ /^\+/) {
	print "chr=$chr\n";
	my $input_fh = read_file($input);
	my $index = 0;
	while (<$input_fh>) {
	    chomp;
	    /^#/ and next;
	    /^\s*$/ and next;
	    my @line = split /\t/, $_;
	    if ($line[0] =~ /(_TY|_TSU4)/) {
		$index++;
		$line[0] = $chr;
		($line[1], $line[2]) = ($line[1] =~ /(REannotate)_(\S+)/);
		$line[8] = "ID=$line[8]"."-"."$line[2]"."$index";
		my $line_out = join "\t", @line;
		print $output_fh "$line_out\n";
	    }
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

