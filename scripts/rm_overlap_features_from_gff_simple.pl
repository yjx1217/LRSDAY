#!/usr/bin/perl
use warnings;
use strict;
use Getopt::Long;

##############################################################
#  script: rm_overlap_features_from_gff_simple.pl
#  author: Jia-Xing Yue (GitHub ID: yjx1217)
#  last edited: 2017.06.17
#  description: remove overlap feature records from gff file regardless of feature types
#  example: perl rm_overlap_features_from_gff_simple.pl  -i input.raw.gff -o output.non_overlap.gff -r genome.fa(.gz)
##############################################################

my ($input, $output, $refseq);

GetOptions('i|input:s' => \$input, # input raw gff file that may contain overlap features
	   'o|output:s' => \$output, # output gff file with only non_overlap features
	   'r|refseq:s' => \$refseq); # reference genome for sorting the gff file

my %rep = ();
my %current = ();
$current{"chr"} = "";
$current{"start"} = "";
$current{"end"} = "";
$current{"record"} = "";
my $input_fh = read_file($input);

my @output = ();
while (<$input_fh>) {
    chomp;
    /^#/ and next;
    /^\s*$/ and next;
    my ($chr, $source, $type, $start, $end, $score, $strand, $frame, $attributes) = split /\t/, $_;
    if ($current{"chr"} eq "") {
	$current{"chr"} = $chr;
	$current{"end"} = $end;
	$current{"start"} = $start;
	$current{"record"} = $_;
    } elsif ($current{"chr"} ne $chr) {
	push @output, $current{"record"};
	$current{"chr"} = $chr;
	$current{"end"} = $end;
	$current{"start"} = $start;
	$current{"record"} = $_;
    } elsif (($current{"end"} < $start) or ($current{"start"} > $end)) {
	push @output, $current{"record"};
	$current{"chr"} = $chr;
	$current{"end"} = $end;
	$current{"start"} = $start;
	$current{"record"} = $_;
    } else {
	my $current_length = $current{"end"} - $current{"start"} + 1;
	my $length = $end - $start + 1;
	if ($current_length < $length) {
	    $current{"chr"} = $chr;
	    $current{"end"} = $end;
	    $current{"start"} = $start;
	    $current{"record"} = $_;
	}
    }
}

my %gff = ();

my $refseq_fh = read_file($refseq);
my @chr = get_chr_list($refseq_fh);


foreach my $record (@output) {
    my ($chr, $source, $type, $start, $end, $score, $strand, $frame, $attributes) = split /\t/, $record;
    $gff{$type}{$chr}{$start} = $record;

}

my $output_fh = write_file($output);	
foreach my $t (sort keys %gff) {
    my %gff_t = %{$gff{$t}};
    foreach my $chr (@chr) {
        my %gff_t_c = ();
        foreach my $c (sort keys %gff_t) {
            if ($chr eq $c) {
                my %gff_t_c_s = %{$gff_t{$c}};
                foreach my $s (sort {$a<=>$b} keys %gff_t_c_s) {
                    print $output_fh "$gff_t_c_s{$s}\n";
                }
            }
        }
    }
}
    



sub read_file {
    my $filename = shift @_;
    my $fh;
    open($fh, $filename) or die "cannot open the input file $filename\n";
    return $fh;
}

sub write_file {
    my $filename = shift @_;
    my $fh;
    open($fh, ">$filename") or die "cannot open the output file $filename\n";
    return $fh;
}  

sub get_chr_list {
    my $fh = shift @_;
    my @chr;
    my $seq_name = "";
    my $flag = 0;
    while (<$fh>) {
        chomp;
        if (/^\s*$/) {
            next;
        } elsif (/^\s*\#/) {
            next;
        } elsif (/^>(.*)/) {
            $seq_name = $1;
            push @chr, $seq_name;
        } 
    }
    return @chr;
}
