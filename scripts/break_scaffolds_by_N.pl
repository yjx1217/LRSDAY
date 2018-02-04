#!/usr/bin/perl
use warnings;
use strict;
use Getopt::Long;

##############################################################
#  script: break_scaffolds_by_N.pl
#  author: Jia-Xing Yue (GitHub ID: yjx1217)
#  last edited: 2017.06.19
#  description: break scaffolds containing assembly gaps into multiple contigs
#  example: perl break_scaffolds_by_N.pl -i input.fa(.gz) -o output.fa(.gz) -g 5000 
##############################################################


my ($input, $output, $gap_size);
$gap_size = 100;
GetOptions('input|i:s' => \$input, # input genome assembly name
	   'gap_size|g:i' => \$gap_size, # the minimal size of assembly gaps that you want to break
           'output|o:s' => \$output); # output genome assembly name after the break

my $input_fh = read_file($input);
my @input = ();
my %input = ();
parse_fasta_file($input_fh, \%input, \@input);

my $output_fh = write_file($output);
my $count = 0;
foreach my $id (@input) {
    if ($input{$id} =~ /N{$gap_size,}/i) {
	$count = 0;
	my @subseq = split /N{$gap_size,}/i, $input{$id};
	foreach my $subseq (@subseq) {
	    $count++;
	    print $output_fh ">${id}_${count}\n$subseq\n";
	}
    } else {
	print $output_fh ">$id\n$input{$id}\n";
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
