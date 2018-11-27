#!/usr/bin/perl
use warnings;
use strict;
use Getopt::Long;

##############################################################
#  script: fasta_length.pl
#  author: Jia-Xing Yue (GitHub ID: yjx1217)
#  last edited: 2018.01.14
#  description: calculate the length of each sequences in the input fasta file
#  example: perl fasta_length.pl -i input.fa(.gz) -o output.length.txt(.gz) -sort_by_length yes
##############################################################

my ($input, $output, $sort_by_length);
$sort_by_length = "no"; # default: no resorting, keep the input sequence order

GetOptions('input|i:s' => \$input,
	   'output|o:s' => \$output,
	   'sort_by_length|s:s' => \$sort_by_length);


my $input_fh = read_file($input);
my @input = ();
my %input = ();
parse_fasta_file($input_fh, \%input, \@input);

my $output_fh = write_file($output);

my %length = ();    
foreach my $id (@input) {
    my $length = length $input{$id};
    $length{$id} = $length;
}

if ($sort_by_length eq "no") {
    foreach my $id (@input) {
	print $output_fh "$id\t$length{$id}\n";
    }
} else {
    foreach my $id (sort {$length{$b} <=> $length{$a}} keys %length) {
	print $output_fh "$id\t$length{$id}\n";
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


sub parse_list_file {
    my $fh = shift @_;
    my %list = ();
    while (<$fh>) {
	chomp;
	$_ =~ s/\s+$//g;
	$list{$_}++;
    }
    return %list;
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


