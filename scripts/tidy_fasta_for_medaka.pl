#! /usr/bin/perl
use strict;
use warnings;
use Getopt::Long;

##############################################################
#  script: tidy_fasta.pl
#  author: Jia-Xing Yue (GitHub ID: yjx1217)
#  last edited: 2019.08.16
#  description: tidy up multi-fasta file by only keeping sequence id (but not the descriptive associated information) and wrapping the sequence with a width of 70 character
#  example: perl tidy_fasta.pl -i input.fa(.gz) -o output.fa(.gz)
##############################################################

my ($input, $output, $width);
$width = 60;
GetOptions('input|i:s' => \$input, # input genome fasta file
	   'output|o:s' => \$output,
	   'width|w:i' => \$width); # output genome fasta file

my $input_fh = read_file($input);
my @input = ();
my %input = ();
parse_fasta_file($input_fh, \%input, \@input);

my $output_fh = write_file($output);

foreach my $id (@input) {
    print $output_fh ">$id\n";
    my $total_length = length $input{$id};
    my $offset = 0;
    while ($total_length - $offset > $width) {
	my $subseq = substr $input{$id}, $offset, $width;
	print $output_fh "$subseq\n";
	$offset += $width;
    }
    my $remainder = substr $input{$id}, $offset;
    print $output_fh "$remainder\n";
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
    my $seq_id = "";
    while (<$fh>) {
        chomp;
        if (/^\s*$/) {
            next;
        } elsif (/^\s*\#/) {
            next;
        } elsif (/^>(\S+)/) {
	    $seq_id = $1;
	    $seq_id =~ s/\:/_/gi;
	    $seq_id =~ s/\./_/gi;
            $$seq_hashref{$seq_id} = "";
            push @$seq_arraryref, $seq_id;
        } else {
            $$seq_hashref{$seq_id} .= (uc $_);
        }
    }
}

