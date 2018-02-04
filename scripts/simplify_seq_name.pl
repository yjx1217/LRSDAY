#!/usr/bin/env perl
use warnings;
use strict;
use Getopt::Long;

##############################################################
#  script: simplify_seq_name.pl
#  author: Jia-Xing Yue (GitHub ID: yjx1217)
#  last edited: 2017.06.17
#  description: simplify sequence ids in fasta file
#  example: perl simplify_seq_name.pl -i input.fa(.gz) -o output.fa(.gz)
##############################################################

my ($input, $output);

GetOptions('i|input:s' => \$input, # input fasta file
	   'o|output:s' => \$output); # output fasta file


my $input_fh = read_file($input);
my @seq_names = ();
my %input = parse_fasta_file($input_fh, \@seq_names);

my $output_fh = write_file($output);
foreach my $seq_name (@seq_names) {
    print $output_fh ">$seq_name\n$input{$seq_name}\n";
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
    my ($fh, $seq_names_arrayref) = @_;
    my %seq;
    my $seq_name = "";
    while (<$fh>) {
        chomp;
        if (/^\s*$/) {
            next;
        } elsif(/^\s*\#/) {
            next;
        } elsif(/^>(\S+)/) {
            $seq_name = $1;
	    push @$seq_names_arrayref, $seq_name;
            $seq{$seq_name} = "";
        } else {
            $seq{$seq_name} .= $_;
        }
    }
    return %seq;
}
