#!/usr/bin/perl
use warnings;
use strict;
use Getopt::Long;

##############################################################
#  script: exonerate2gene_maker.pl
#  author: Jia-Xing Yue (GitHub ID: yjx1217)
#  last edited: 2018.01.30
#  description: turn exonerate hits in maker evidence to gene models
#  example: perl exonerate2gene_maker.pl -i input.exonerate_hit.gff3 -o output.exonerate_hit.gene_model.gff3
##############################################################

my ($input, $output);
GetOptions('input|i:s' => \$input, # input est2genome or protein2genome GFF3 file from maker
	   'output|o:s' => \$output); # output est2genome or protein2genome gene model GFF3 file


my $input_fh = read_file($input);
my $output_fh = write_file($output);
print $output_fh "##gff-version 3\n";

my $gene_id;
my $gene_name;
my $mRNA_id;
my $exon_id;
my $cds_id;
my $exon_index = 0;
while (<$input_fh>) {
    chomp;
    /^##FASTA/ and last;
    /^#/ and next;
    my ($chr, $source, $type, $start, $end, $score, $strand, $phase, $attributes) = split /\t/, $_;
    if ($type =~ /(protein_match|expressed_sequence_match)/) {
	($gene_id, $gene_name) = ($attributes =~ /ID=([^;]+);\S*Name=([^;]+)/);
	print $output_fh "$chr\tmaker_complementary\tgene\t$start\t$end\t$score\t$strand\t$phase\tID=$gene_id;Name=$gene_name\n";
	$mRNA_id = "$gene_id.mRNA";
	print $output_fh "$chr\tmaker_complementary\tmRNA\t$start\t$end\t$score\t$strand\t$phase\tID=$mRNA_id;Parent=$gene_id\n";
	$exon_index = 0;
    } elsif ($type =~ /match_part/) {
	$exon_index++;
	$exon_id = "$mRNA_id.exon.$exon_index";
	$cds_id = "$mRNA_id.cds.$exon_index";
	print $output_fh "$chr\tmaker_complementary\texon\t$start\t$end\t$score\t$strand\t$phase\tID=$exon_id;Parent=$mRNA_id\n";
	print $output_fh "$chr\tmaker_complementary\tCDS\t$start\t$end\t$score\t$strand\t$phase\tID=$cds_id;Parent=$mRNA_id\n";
    } else {
	die "unrecognized feature type: $type\n";
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

