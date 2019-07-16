#!/usr/bin/perl
use warnings;
use strict;
use Getopt::Long;

##############################################################
#  script: relabel_and_reorder_sequences.pl
#  author: Jia-Xing Yue (GitHub ID: yjx1217)
#  last edited: 2017.06.17
#  description: relabel and reorder sequences
#  example: perl relabel_and_reorder_sequences.pl -i input.fa -o output.fa -m modification_instruction_list
##############################################################


my ($input, $output, $modification);

GetOptions('i|input:s' => \$input, # input sequences in fasta format 
	   'o|output:s' => \$output, # output sequences in fasta format
           'm|modification:s' => \$modification); # manually edited modification instruction list (each row goes like: old_id,orientation,new_id(if renaming is needed))


my $input_fh = read_file($input);
my %input = parse_fasta_file($input_fh);

my $modification_fh = read_file($modification);

my $output_fh = write_file($output);

while (<$modification_fh>) {
    chomp;
    /^#/ and next;
    /^\s*$/ and next;
    if ($_ =~ /,/) {
	my ($old_name, $orientation, $new_name) = split /,/;
	my $seq = $input{$old_name};
	if ($orientation eq '-') {
	    $seq  = revcom($seq);
	}
	print $output_fh ">$new_name\n$seq\n";
    } else {
	my $name = $_;
	my $seq = $input{$name};
	print $output_fh ">$name\n$seq\n";
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
    my $fh = shift @_;
    my %seq;
    my $seq_name = "";
    my $flag = 0;
    while (<$fh>) {
        chomp;
        if(/^\s*$/) {
            next;
        } elsif (/^\s*\#/) {
            next;
        } elsif (/^>(\S+)/) {
            $seq_name = $1;
            $seq{$seq_name} = "";
        } else {
            $seq{$seq_name} .= $_;
        }
    }
    return %seq;
}

sub revcom {
    my $seq = shift @_;
    my $seq_revcom = reverse $seq;
    $seq_revcom =~ tr/ATGCNatgcn/TACGNtacgn/;
    return $seq_revcom;
}
