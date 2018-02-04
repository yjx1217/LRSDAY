#!/usr/bin/perl
use warnings;
use strict;
use Getopt::Long;

##############################################################
#  script: switch_letter_cases_in_fasta.pl
#  author: Jia-Xing Yue (GitHub ID: yjx1217)
#  last edited: 2018.01.11
#  description: convert sequences in a fasta file to all uppercases or to all lowercases
#  example: perl switch_letter_cases_in_fasta.pl -i input.fa(.gz) -o output.fa(.gz) -c upper # to transform all lowercases to uppercases 
#  example: perl switch_letter_cases_in_fasta.pl -i input.fa(.gz) -o output.fa(.gz) -c lower # to transform all uppercases to lowercases 
##############################################################

my ($input, $output, $case);
$case = "upper"; # normal mode, output sequences found in the list. 

GetOptions('input|i:s' => \$input,
	   'output|o:s' => \$output,
	   'case|c:s' => \$case); # case can be 'lower'(default) or 'upper'. 

my $input_fh = read_file($input);
my @input = ();
my %input = ();
parse_fasta_file($input_fh, \%input, \@input);


my $output_fh = write_file($output);
    
foreach my $id (@input) {
    if ($case =~ /^upper/) {
	my $seq = uc $input{$id};
	print $output_fh ">$id\n$seq\n";
    } elsif ($case =~ /^lower/) {
	my $seq = lc $input{$id};
	print $output_fh ">$id\n$seq\n";
    } else {
	die "unrecognized case (-c) option: $case! Please only use 'upper' or 'lower'!\n";
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


