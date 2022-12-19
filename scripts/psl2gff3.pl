#!/usr/bin/perl
use warnings;
use strict;
use Getopt::Long;

##############################################################
#  script: psl2gff3.pl
#  author: Jia-Xing Yue (GitHub ID: yjx1217)
#  last edited: 2017.06.17
#  description: convert psl format to gff3 format
#  example: perl psl2gff3.pl -i input.psl -o output.gff3 -r genome.fa(.gz) -f feature -t genome_tag -l 100
##############################################################

my ($input, $output, $tag, $refseq, $feature_type, $length_cutoff_for_completeness);
$length_cutoff_for_completeness = 0;

GetOptions('input|i:s' => \$input, # input psl file
	   'output|o:s' => \$output, # output gff3 file
	   'refseq|r:s' => \$refseq, # reference genome file in fasta format, could be compressed with gzip
	   'feature_type|ft:s' => \$feature_type, # feature type label for the output
	   'tag|t:s' => \$tag, # genome_tag label for the output
	   'length_cutoff_for_completeness|l:s' => \$length_cutoff_for_completeness); # length cutoff for labeling whether the feature is complete or partial


my $refseq_fh = read_file($refseq);
my %refseq = parse_fasta_file($refseq_fh);
my $input_fh = read_file($input);
my $output_fh = write_file($output);

while (<$input_fh>) {
    chomp;
    if (/^\d+/) {
	my @line = split /\t/, $_;
	my $strand = $line[8];
	my $chr = $line[13];
	my $start = $line[15] + 1; 
	my $end = $line[16];
	if ($start > $end) {
	    ($start, $end) = ($end, $start);
	}
	my $feature_length = $end - $start + 1;
	if ($feature_length >= $length_cutoff_for_completeness) {
	    my $new_feature_type = $feature_type;
	    if ($feature_length < $length_cutoff_for_completeness) {
		$new_feature_type = "${feature_type}_partial";
	    }
	    my $id = "$new_feature_type:$chr:$start-$end:$strand";
	    print $output_fh "$chr\t$tag\t$new_feature_type\t$start\t$end\t.\t$strand\t.\tID=$id;Name=$id\n";
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
        } elsif (/^\s*\#/) {
            next;
        } elsif (/^>(.*)/) {
            $seq_name = $1;
            $seq{$seq_name} = "";
        } else {
            $seq{$seq_name} .= $_;
        }
    }
    return %seq;
}
