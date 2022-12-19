#!/usr/bin/perl
use warnings;
use strict;
use Getopt::Long;

##############################################################
#  script: tidy_truncated_TY_gff.pl
#  author: Jia-Xing Yue (GitHub ID: yjx1217)
#  last edited: 2017.06.17
#  description: tidy truncated TY gff annotated by REannotate
#  example: perl tidy_truncated_TY_gff.pl -i raw.truncated_TY.gff -o tidy.TY_truncated.gff3 -t genome_tag
##############################################################

my ($input, $output, $tag);
GetOptions('input|i:s' => \$input, # raw input gff
	   'output|o:s' => \$output, # tidy output gff3
	   'tag|t:s' => \$tag); # genome_tag

my $input_fh = read_file($input);
my $output_fh = write_file($output);

while (<$input_fh>) {
    chomp;
    /^#/ and next;
    /^\s*$/ and next;
    my ($chr, $source, $type, $start, $end, $score, $strand, $phase, $attributes) = split /\t/, $_;
    my ($new_type) = ($type =~ /(TY\d+|TSU4)/);
    $new_type .= "_truncated";
    my $new_id = "$new_type:$chr:$start-$end:$strand";
    print $output_fh "$chr\t$tag\t$new_type\t$start\t$end\t.\t$strand\t.\tID=${new_id};Name=${new_id}\n";
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


