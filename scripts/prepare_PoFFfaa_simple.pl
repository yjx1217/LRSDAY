#!/usr/bin/perl
use warnings;
use strict;
use Getopt::Long;

##############################################################
#  script: prepare_PoFFfaa_simple.pl
#  author: Jia-Xing Yue (GitHub ID: yjx1217)
#  last edited: 2017.06.17
#  description: prepare the fasta-formated proteome (.faa) or CDSs (.ffn) sequence files for proteinortho run with PoFF option support 
#  example: perl prepare_PoFFffn_simple.pl -i prefix.pep.fa -o prefix.PoFF.faa  
##############################################################

my ($input, $output);

GetOptions('input|i:s' => \$input, # input proteome or CDSs sequence files 
           'output|o:s' => \$output); # output proteome or CDSs sequence files for proteinortho

my $input_fh = read_file($input);
my $output_fh = write_file($output);

my %seq = ();
my $gene_id;
my $mRNA_id;



while (<$input_fh>) {
    chomp;
    /^#/ and next;
    /^\s*$/ and next;
    if (/^>/) {
	($gene_id, $mRNA_id) = ($_ =~ />(\S+?)\|(\S+?)\|?/);
	$seq{$gene_id}{$mRNA_id} = "";
    } else {
	$seq{$gene_id}{$mRNA_id} .= $_;
    }
}


# select longest transcript to represent the gene
foreach my $gene_id (sort keys %seq) {
    my @mRNA_id = sort keys %{$seq{$gene_id}};
    my %length = ();
    foreach my $mRNA_id (@mRNA_id) {
	$length{$mRNA_id} = length $seq{$gene_id}{$mRNA_id};
    }
    foreach my $mRNA_id (sort {$length{$b}<=>$length{$a}} keys %length) {
	print $output_fh  ">${gene_id}\n$seq{$gene_id}{$mRNA_id}\n";
	last;
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
    if ($file =~ /\.gz$/) {
        open($fh, "| gzip -c >$file") or die "can't open $file\n";
    } else {
        open($fh, ">$file") or die "can't open $file\n";
    }
    return $fh;
}  
