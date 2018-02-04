#!/usr/bin/perl
use warnings;
use strict;
use Getopt::Long;

##############################################################
#  script: filter_blast_result.pl
#  author: Jia-Xing Yue (GitHub ID: yjx1217)
#  last edited: 2017.06.17
#  description: filter blast tabular output by percent of identity, alignment length, e-value, bit score, and query coverage 
#  example: perl filter_blast_result.pl -i input.blast.out -o output.blast.filtered.out -pct_identity_cutoff 70 -aln_length_cutoff 100 -e_value_cutoff 0.01 -bit_score_cutoff 50
##############################################################


my ($input, $output, $pct_identity_cutoff, $aln_length_cutoff, $e_value_cutoff, $bit_score_cutoff, $query_cov_cutoff);
$pct_identity_cutoff = 0;
$aln_length_cutoff = 0;
$e_value_cutoff = 100;
$bit_score_cutoff = 0;
$query_cov_cutoff = 0;

GetOptions('input|i:s' => \$input, # input blast tabular output
	   'output|o:s' => \$output, # filtered blast tabular output
	   'pct_identity_cutoff|p:s' => \$pct_identity_cutoff, # percent of identity cutoff [0-100]
	   'aln_length_cutoff|l:s' => \$aln_length_cutoff, # alignment length cutoff
	   'e_value_cutoff|e:s' => \$e_value_cutoff, # e-value cutoff
	   'bit_score_cutoff|b:s' => \$bit_score_cutoff, # bit score cutoff
	   'query_cov_cutoff|q:s' =>\$query_cov_cutoff); # query coverage cutoff if this output is specified when doing blast

my $input_fh = read_file($input);
my $output_fh = write_file($output);

while (<$input_fh>) {
    chomp;
    if (/^\#/) {
	print $output_fh "$_\n";
    } elsif (/^\s*$/) {
	print $output_fh "$_\n";
    } else {
	my @line = split /\t/, $_;
	my $pct_identity = $line[2];
	my $aln_length = $line[3];
	my $e_value = $line[10];
	my $bit_score = $line[11];
	my $query_cov = "";
	if (scalar @line == 13) {
	    # when specified extra blast output columns: -outfmt '7 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore qlen qcovs'"
	    $query_cov = $line[13];
	}
	my $flag = 0;
	if ($pct_identity < $pct_identity_cutoff) {
	    $flag = 1;
	} elsif ($aln_length < $aln_length_cutoff) {
	    $flag = 1;
	} elsif ($e_value > $e_value_cutoff) {
	    $flag = 1;
	} elsif ($bit_score < $bit_score_cutoff) {
	    $flag = 1;
	} elsif (($query_cov ne "") and ($query_cov < $query_cov_cutoff)) {
	    $flag = 1;
	}
	if ($flag == 0) {
	    print $output_fh "$_\n";
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




