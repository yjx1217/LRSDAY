#!/usr/bin/perl
use warnings;
use strict;
use Getopt::Long;

##############################################################
#  script: join_contigs_by_N.pl
#  author: Jia-Xing Yue (GitHub ID: yjx1217)
#  last edited: 2017.06.19
#  description: join a multi-fasta file by order into a single sequence with the specified gaps
#  example: perl join_contigs_by_N.pl -i input.fa(.gz) -o output.fa(.gz) -g 5000 -t output_tag
##############################################################

my ($input, $output, $gap_size, $tag);
$gap_size = 100;
GetOptions('input|i:s' => \$input, # input multi-fasta file
           'output|o:s' => \$output, # file name for the output
	   'gap|g:i' => \$gap_size, # gap size
	   'tag|t:s' => \$tag); # sequence name for the output

my $input_fh = read_file($input);
my @input = ();
my %input = ();
parse_fasta_file($input_fh, \%input, \@input);

my $outseq = "";
my $seq_left = scalar @input;

my $spacer = 'N' x $gap_size;

foreach my $id (@input) {
    $outseq .= $input{$id};
    $seq_left--;
    if ($seq_left > 0) {
	$outseq .= $spacer;
    }
}

my $output_fh = write_file($output);
print $output_fh ">$tag\n$outseq\n";



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
