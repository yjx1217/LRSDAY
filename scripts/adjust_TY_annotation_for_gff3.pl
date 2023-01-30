#!/usr/bin/perl
use warnings;
use strict;
use Getopt::Long;

##############################################################
#  script: adjust_TY_annotation_for_gff3.pl
#  author: Jia-Xing Yue (GitHub ID: yjx1217)
#  last edited: 2023.01.30
#  description: change TY4_soloLTR, TSU4_soloLTR, TY4_truncated, TSU4_truncated to TY4/TSU4_soloLTR, TY4/TSU4_soloLTR, TY4/TSU4_truncated, and TY4/TSU4_truncated respectively; remove "p" suffix after TY.
#  example: perl adjust_TY4TSU4_annotation_for_gff3.pl -i input.gff3(.gz) -o output.gff3(.gz)
##############################################################

my ($input, $output);
GetOptions('input|i:s' => \$input, # input TY list file
	   'output|o:s' => \$output); # output TY GFF3 file

my $input_fh = read_file($input);
my $output_fh = write_file($output);
while (<$input_fh>) {
    chomp;
    if ((/^\s*\#/) or (/^\#/)) {
	print $output_fh "$_\n";
    } else {
	my $flag = 0;
	my @line = split /\t/, $_;
	my $chr = $line[0];
	my $source = $line[1];
	my $feature = $line[2];
	my $start = $line[3];
	my $end = $line[4];
	my $score = $line[5];
	my $strand = $line[6];
	my $phase = $line[7];
	my $attributes = $line[8];
	if ($feature eq "mobile_element") {
	    my ($mobile_element_type) = ($attributes =~ /mobile_element_type=([^;]+)/);
	    $mobile_element_type =~ s/TY1p/TY1/gi;
	    $mobile_element_type =~ s/TY2p/TY2/gi;
	    $mobile_element_type =~ s/TY3p/TY3/gi;
	    $mobile_element_type =~ s/TY4p/TY4/gi;
	    $mobile_element_type =~ s/TY5p/TY5/gi;
	    $mobile_element_type =~ s/TSU4p/TSU4/gi;
	    $mobile_element_type =~ s/TSU4b/TSU4/gi;
	    if (($mobile_element_type =~ /(TSU4|TY4|TY4\/TSU4|TSU4\/TY4)/) and ($mobile_element_type =~ /truncated/)) {
		$mobile_element_type = "TY4/TSU4_truncated";
	    } elsif (($mobile_element_type =~ /(TSU4|TY4|TY4\/TSU4|TSU4\/TY4)/) and ($mobile_element_type =~ /soloLTR/)) {
		$mobile_element_type = "TY4/TSU4_soloLTR";
	    }
	    print $output_fh "$chr\t$source\tmobile_element\t$start\t$end\t$score\t$strand\t$phase\tID=$mobile_element_type:$chr:$start-$end:$strand;Name=$mobile_element_type:$chr:$start-$end:$strand;mobile_element_type=$mobile_element_type\n";
	} else {
	    print $output_fh "$_\n";
	}
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


sub parse_fasta_file {
    my $fh = shift @_;
    my %seq;
    my $seq_name = "";
    my $flag = 0;
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
    $seq_revcom =~ tr/ATGCatgc/TACGtacg/;
    return $seq_revcom;
}
