#!/usr/bin/perl
use warnings;
use strict;
use Getopt::Long;

##############################################################
#  script: break_contig_by_blast.pl
#  author: Jia-Xing Yue (GitHub ID: yjx1217)
#  last edited: 2018.01.24
#  description: break contigs containing blast hit into multiple contigs
#  example: perl break_contig_by_blast.pl -i input.fa(.gz) -o output.fa(.gz) -b blast.result.outfmt7.txt 
##############################################################


my ($input, $output, $blast);

GetOptions('input|i:s' => \$input, # input genome assembly name
	   'blast|b:s' => \$blast,
           'output|o:s' => \$output); # output genome assembly name after the break

my $input_fh = read_file($input);
my @input = ();
my %input = ();
parse_fasta_file($input_fh, \%input, \@input);

my $output_fh = write_file($output);

my $blast_fh = read_file($blast);
my %blast = parse_blast_file($blast_fh);

foreach my $id (@input) {
    if (exists $blast{$id}) {
	print $output_fh ">$id.1\n$blast{$id}{'part1'}\n";
	print $output_fh ">$id.2\n$blast{$id}{'part2'}\n";
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


sub parse_blast_file {
    my $fh = shift @_;
    my %blast = ();
    while (<$fh>) {
	chomp;
	(/^\#/) and next;
	(/^\s*$/) and next;
        my @line = split /\t/, $_;
	my $subject = $line[1];
	my $subject_start = $line[8];
	my $subject_end = $line[9];
        my $pct_identity = $line[2];
        my $aln_length = $line[3];
        my $e_value = $line[10];
        my $bit_score = $line[11];
	if ($subject_start <= $subject_end) {
	    # hit on the positive strand
	    $blast{$subject}{'part1'} = substr $input{$subject}, 0, $subject_start - 1;
	    $blast{$subject}{'part2'} = substr $input{$subject}, $subject_start - 1;
	} else {
	    # hit on the negative strand
	    $blast{$subject}{'part1'} = substr $input{$subject}, 0, $subject_end;
	    $blast{$subject}{'part2'} = substr $input{$subject}, $subject_end;
	}
    }
    return %blast;
}

