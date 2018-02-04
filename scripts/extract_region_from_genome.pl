#!/usr/bin/perl
use warnings;
use strict;
use Getopt::Long;

##############################################################
#  script: extract_region_from_genome.pl
#  author: Jia-Xing Yue (GitHub ID: yjx1217)
#  last edited: 2017.06.19
#  description: extract any given region from a multi-fasta file by the query string: "chr:start-end:strand"
#  example: perl extract_region_from_genome.pl -i genome.fa(.gz) -o query.fa(.gz) -q chrI:1000-4000:+ -f 100
##############################################################


my ($input, $output, $query, $flanking);
$flanking = 0;
GetOptions('input|i:s' => \$input, # input multi fasta file
	   'output|o:s' => \$output, # output query fasta file
	   'query|q:s' => \$query, # query string, format: chr:start-end:strand
	   'flanking|f:s' => \$flanking); # size of flanking regions for the output

my $input_fh = read_file($input);
my %input = ();
my @input = ();
parse_fasta_file($input_fh, \%input, \@input);
close $input_fh;

my $output_fh = write_file($output);

my ($chr, $start, $end, $strand) = ($query =~ /(\S+):(\d+)\-(\d+):(\+|\-)/);
print "chr=$chr, start=$start, end=$end, strand=$strand, flanking=$flanking\n";
my $chr_len = length $input{$chr};
my $chr_left = $start - $flanking;
my $chr_right = $chr_len - ($end + $flanking);
my $subject = substr $input{$chr}, $start - 1, $end - $start + 1;
my $left_flank;
my $right_flank;
if ($chr_left >=0) {
    $left_flank = substr $input{$chr}, $start - 1 -$flanking, $flanking;
} else {
    $left_flank = substr $input{$chr}, 0, $start -1;
    my $left_padding = 'N' x (-$chr_left);
    $left_flank = $left_padding.$left_flank;
}

if ($chr_right >= 0) {
    $right_flank = substr $input{$chr}, $end, $flanking;
} else {
    $right_flank = substr $input{$chr}, $end, $chr_len - $end;
    my $right_padding = 'N' x (-$chr_right);
    $right_flank = $right_flank . $right_padding;
}
$subject = $left_flank.$subject.$right_flank;


if ($strand eq '-') {
    $subject = revcom($subject);
}

if ($flanking > 0) {
    my $new_start = $start - $flanking;
    my $new_end = $end + $flanking;
    $query = "$chr:${new_start}-${new_end}:+";
}
print $output_fh  ">$query\n$subject\n";




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
    my $seq_name = "";
    while (<$fh>) {
        chomp;
        if (/^\s*$/) {
            next;
        } elsif (/^\s*\#/) {
            next;
        } elsif (/^>(.*)/) {
            $seq_name = $1;
            $$seq_hashref{$seq_name} = "";
            push @$seq_arraryref, $seq_name;
        } else {
            $$seq_hashref{$seq_name} .= $_;
        }
    }
}


sub revcom {
    my $seq = shift @_;
    my $seq_revcom = reverse $seq;
    $seq_revcom =~ tr/ATGCNatgcn/TACGNtacgn/;
    return $seq_revcom;
}
