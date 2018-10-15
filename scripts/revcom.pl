#!/usr/bin/perl
use warnings;
use strict;
use Getopt::Long;

##############################################################
#  script: revcom.pl
#  author: Jia-Xing Yue (GitHub ID: yjx1217)
#  last edited: 2018.09.01
#  description: generate reverse complementary sequences for the input FASTA sequences.
#  example: perl revcom.pl -i input.fa(.gz) -o output.fa(.gz)
##############################################################


my ($input, $output);
GetOptions('input|i:s' => \$input,
           'output|o:s' => \$output);

my $input_fh = read_file($input);
my %input = ();
my @input = ();
parse_fasta_file($input_fh, \%input, \@input);

my $output_fh = write_file($output);

foreach my $id (@input) {
    my $original_seq = $input{$id};
    my $revcom_seq = revcom($original_seq);
    print $output_fh ">$id\n$revcom_seq\n";
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
    my ($fh, $input_hashref, $input_arrayref) = @_;
    my $seq_name = "";
    while (<$fh>) {
        chomp;
        if (/^\s*$/) {
            next;
        } elsif (/^\s*#/) {
            next;
        } elsif (/^>(.*)/) {
            $seq_name = $1;
            push @$input_arrayref, $seq_name;
            $$input_hashref{$seq_name} = "";
        } else {
            $$input_hashref{$seq_name} .= $_;
        }
    }
}


sub revcom {
    my $seq = shift @_;
    my $seq_revcom = reverse $seq;
    $seq_revcom =~ tr/ATGCNatgcn/TACGNtacgn/;
    return $seq_revcom;
}
