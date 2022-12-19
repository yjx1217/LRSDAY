#!/usr/bin/perl
use warnings;
#use warnings FATAL => 'all';
use strict;
use Getopt::Long;

##############################################################
#  script: tidy_TE_gff3.pl
#  author: Jia-Xing Yue (GitHub ID: yjx1217)
#  last edited: 2022.12.11
#  description: tidy TE gff3 by re-sorting and re-naming feature IDs
#  example: perl tidy_TE_gff3.pl -i raw.gff3 -t genome_tag -o tidy.gff3 -r genome.fa(.gz)
##############################################################

my ($refseq, $input, $output, $tag);
GetOptions('refseq|r:s' => \$refseq, # reference genome file in fasta format; can be compressed by gzip
	   'input|i:s' => \$input, # input raw gff3
	   'output|o:s' => \$output, # output tidy gff3
	   'tag|t:s' => \$tag); # genome tag

my $refseq_fh = read_file($refseq);
my %refseq = ();
my @refseq = ();
parse_fasta_file($refseq_fh, \%refseq, \@refseq);

my $input_fh = read_file($input);
my %gff = parse_gff_file($input_fh);
close $input_fh;

my $output_fh = write_file($output);

print $output_fh "##gff-version 3\n";

foreach my $chr (@refseq) {
    my $chr_length = length $refseq{$chr};
    print $output_fh "##sequence-region $chr 1 $chr_length\n";
}

my $feature_index = 0;
foreach my $chr (@refseq) { 
    my @features_on_chr = ();
    foreach my $feature_id (sort keys %gff) {
	if ($gff{$feature_id}{'chr'} eq $chr) {
	    push @features_on_chr, $feature_id;
	}
    }
    foreach my $feature_id (sort {$gff{$a}{'start'} <=> $gff{$b}{'start'} or $gff{$a}{'end'} <=> $gff{$b}{'end'}} @features_on_chr) {
	my $feature_name = $gff{$feature_id}{'feature_name'};
	my $feature_note = $gff{$feature_id}{'feature_note'};
	my $feature_type = $gff{$feature_id}{'feature_type'};
	my $feature_start = $gff{$feature_id}{'start'};
	my $feature_end = $gff{$feature_id}{'end'};
	my $feature_score = $gff{$feature_id}{'score'};
	my $feature_strand = $gff{$feature_id}{'strand'};
	my $feature_phase = $gff{$feature_id}{'phase'};
	print $output_fh "$chr\t$tag\t$feature_type\t$feature_start\t$feature_end\t$feature_score\t$feature_strand\t$feature_phase\tID=$feature_id;Name=$feature_name;$feature_note\n";
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

sub parse_gff_file {
    my $fh = shift @_;
    my %gff = ();
    my $feature_id;
    my $feature_name;
    my $feature_type;
    while (<$fh>) {
	chomp;
	/^##FASTA/ and last;
	/^#/ and next;
	my ($chr, $source, $type, $start, $end, $score, $strand, $phase, $attributes) = split /\t/, $_;
	# print "$_\n";
	($feature_id, $feature_name) = ($attributes =~ /ID=([^;]+);\S*Name=([^;]+)/);
	$gff{$feature_id}{'feature_type'} = "mobile_element";
	$gff{$feature_id}{'feature_name'} = $feature_name;
	$gff{$feature_id}{'feature_note'} = "mobile_element_type=$type";
	$gff{$feature_id}{'chr'} = $chr;
	$gff{$feature_id}{'start'} = $start;
	$gff{$feature_id}{'end'} = $end;
	$gff{$feature_id}{'strand'} = $strand;
	$gff{$feature_id}{'source'} = $source;
	$gff{$feature_id}{'score'} = $score;
	$gff{$feature_id}{'phase'} = $phase;
    }
    return %gff;
}

