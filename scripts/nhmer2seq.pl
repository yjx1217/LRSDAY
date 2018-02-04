#!/usr/bin/perl
use warnings;
use strict;
use Getopt::Long;

##############################################################
#  script: nhmer2seq.pl
#  author: Jia-Xing Yue (GitHub ID: yjx1217)
#  last edited: 2017.06.17
#  description: extract nhmer hit sequences 
#  example: perl nhmer2seq.pl -i S288C.X_element.nhmer.out -p S288C -r genome.fa(.gz) -e 0.001 -f X_element
##############################################################


my ($input, $prefix, $refseq, $evalue_cutoff, $feature_type, $length_cutoff_for_completeness);
$evalue_cutoff = 0.001;
$length_cutoff_for_completeness = 0;
GetOptions('input|i:s' => \$input, # nhmmer output
	   'prefix|p:s' => \$prefix, # file name prefix for output files
	   'refseq|r:s' => \$refseq, # reference genome sequences
	   'evalue|e:s' => \$evalue_cutoff, # e-value cutoff
	   'feature_type|ft:s' => \$feature_type, # feature type tag for naming the output sequences
	   'length_cutoff_for_completeness|l:s' => \$length_cutoff_for_completeness); # length cutoff for labeling whether the feature is complete or partial

my $input_fh = read_file($input);
my $output_seq = "$prefix.$feature_type.raw.fa";
my $output_seq_fh = write_file($output_seq);

my $output_gff = "$prefix.$feature_type.raw.gff3";
my $output_gff_fh = write_file($output_gff);

my @candidates = ();

while (<$input_fh>) {
    chomp;
    /^#/ and next;
    /^\s*$/ and next;
    my @line = split /\s+/, $_;
    my $chr_t = $line[0];
    my $start_t = $line[8];
    my $end_t = $line[9];
    my $strand_t = $line[11];
    my $eVal = $line[12];
    if ($eVal < $evalue_cutoff) {
	my $score = $line[13];
	if ($start_t > $end_t) {
	    ($start_t, $end_t) = ($end_t, $start_t);
	}
	my $new_feature_type = $feature_type;
	my $feature_length = $end_t - $start_t + 1;
	if ($feature_length < $length_cutoff_for_completeness) {
	    $new_feature_type = "${feature_type}_partial";
	}
	my $id = "$new_feature_type:$chr_t:${start_t}-${end_t}:$strand_t";
	print $output_gff_fh "$chr_t\t$prefix\t$new_feature_type\t${start_t}\t${end_t}\t.\t${strand_t}\t.\tID=$id;Name=$id\n";
	push @candidates, "$id:$chr_t:${start_t}-${end_t}:$strand_t:$eVal:$score";
    }
}

my $refseq_fh = read_file($refseq);
my %refseq = parse_fasta_file($refseq_fh);

foreach my $candidate (@candidates) {
    my ($id, $chr, $pos, $strand, $eVal, $score) = split ':', $candidate;
    my ($start, $end) = split '-', $pos;
    ####
    my $seq = substr $refseq{$chr}, $start-1, $end-$start+1;
    $seq = uc $seq;
    if ($strand eq '-') {
	$seq = revcom($seq);
    }
    print $output_seq_fh ">$id\n$seq\n";
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
    while (<$fh>) {
        chomp;
        if (/^\s*$/) {
            next;
        } elsif (/^\s*#/) {
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
