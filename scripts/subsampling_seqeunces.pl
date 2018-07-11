#!/usr/bin/perl
use warnings;
use strict;
use Getopt::Long;
use List::Util qw/shuffle/;

##############################################################
#  script: subsampling_sequences.pl
#  author: Jia-Xing Yue (GitHub ID: yjx1217)
#  last edited: 2018.05.01
#  description: subsampling input fasta or fastq sequences
#  example: perl subsampling_sequences.pl -i input.fq(.gz) -f fastq -s 0.1 -m longest -p output # sampling the 10% longest sequences; 
##############################################################

my ($input, $prefix, $format, $mode, $sampling_portion);
$prefix = "output";
$format = "fasta";
$mode = "random";
$sampling_portion = 1;

GetOptions('i|input:s' => \$input,
	   'p|prefix:s' => \$prefix,
	   'f|format:s' => \$format, # fasta or fastq
	   'm|mode:s' => \$mode, # longest or random
	   's|sampling_portion:s' => \$sampling_portion); # 0.05 for top 5%

my @input = ();
my %input = ();
my $input_fh = read_file($input);

if ($format eq "fasta") {
    parse_fasta_file($input_fh, \%input, \@input);
} elsif ($format eq "fastq") {
    parse_fastq_file($input_fh, \%input, \@input);
} else {
    die "Error! Unknown format: $format! Please set the '-f' option as either 'fasta' or 'fastq'. Exit!\n";
}
    

my $total_length = 0;
my %length = ();
foreach my $id (@input) {
    if ($format eq "fasta") {
	$length{$id} = length $input{$id};
    } else {
	$length{$id} = length $input{$id}{'seq'};
    }
    $total_length += $length{$id};
}

my $length_cutoff = $total_length * $sampling_portion;

my $output_sampled;
my $output_sampled_fh;
my $output_remained;
my $output_remained_fh;
if ($format eq "fasta") {
    $output_sampled = "$prefix.sampled.fa.gz";
    $output_remained = "$prefix.remained.fa.gz";
} else {
    $output_sampled = "$prefix.sampled.fq.gz";
    $output_remained = "$prefix.remained.fq.gz";
}    

$output_sampled_fh = write_file($output_sampled);
$output_remained_fh = write_file($output_remained);

my $cummulative_length = 0;
my @output = ();
if ($mode eq "longest") {
    @output = sort {$length{$b} <=> $length{$a}} keys %length;
} elsif ($mode eq "random") {
    @output = shuffle keys %length;
} else {
    die "Error! Unknown mode: $mode! Please set the '-m' option as either 'longest' or 'random'. Exit!\n";
}

foreach my $id (@output) {
    $cummulative_length += $length{$id};
    if ($cummulative_length < $length_cutoff) {
	if ($format eq "fasta") {
	    print $output_sampled_fh ">$id\n$input{$id}\n";
	} else {
	    print $output_sampled_fh "\@${id}\n$input{$id}{'seq'}\n\+\n$input{$id}{'qual'}\n";
	}
    } else {
	if ($format eq "fasta") {    
	    print $output_remained_fh ">$id\n$input{$id}\n";
	} else {
	    print $output_remained_fh "\@${id}\n$input{$id}{'seq'}\n\+\n$input{$id}{'qual'}\n";
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


sub parse_fasta_file {
    my ($fh, $seq_hashref, $seq_arraryref) = @_;
    my $id = "";
    while (<$fh>) {
	chomp;
	if (/^\s*$/) {
	    next;
	} elsif (/^\s*\#/) {
	    next;
	} elsif (/^>(.*)/) {
	    $id = $1;
	    $$seq_hashref{$id} = "";
	    push @$seq_arraryref, $id;
	} else {
	    $$seq_hashref{$id} .= $_;
	}
    }
}


sub parse_fastq_file {
    my ($fh, $seq_hashref, $seq_arrayref) = @_;
    my $id = "";
    while (<$fh>) {
	chomp;
	if ($. % 4 == 1) {
	    ($id) = ($_ =~ /^\@(.*)/);
	    push @$seq_arrayref, $id;
	} elsif ($. % 4 == 2) {
	    $$seq_hashref{$id}{'seq'} = $_;
	} elsif ($. % 4 == 0) {
	    $$seq_hashref{$id}{'qual'} = $_;
	}
    }
}

