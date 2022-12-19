#!/usr/bin/perl
use warnings;
use strict;
use Getopt::Long;

##############################################################
#  script: TY_list2gff3.pl
#  author: Jia-Xing Yue (GitHub ID: yjx1217)
#  last edited: 2022.12.11
#  description: covert TY_list to GFF3
#  example: perl TY_list2gff3.pl -i ty.list -o ty.gff3 -t genome_tag
##############################################################

my ($input, $output, $tag);
GetOptions('input|i:s' => \$input, # input TY list file
	   'output|o:s' => \$output, # output TY GFF3 file
	   'tag|t:s' => \$tag); # genome_tag

my $input_fh = read_file($input);
my $output_fh = write_file($output);

while (<$input_fh>) {
    chomp;
    /^#/ and next;
    /^\s*$/ and next;
    $_ =~ s/\s+$//g;
    my ($id, $type, $chr, $start, $end, $strand) = ($_ =~ /^(\S+?):(\S+?):(\S+?):(\d+)\-(\d+):(\+|\-)/);
    if ($id =~ /^t/) {
	$type .= "_truncated";
    } elsif ($id =~ /^s/) {
	$type .= "_soloLTR";
    } elsif ($id !~ /(u|i)/) {
	print "unrecognized TE id: $id\n";
	next;
    }
    my $new_id = "$type:$chr:$start-$end:$strand";
    print $output_fh "$chr\t$tag\t$type\t$start\t$end\t.\t$strand\t.\tID=${new_id};Name=${new_id}\n";
    # print $output_fh "$chr\t$tag\tmobile_element\t$start\t$end\t.\t$strand\t.\tID=${new_id};Name=${new_id};mobile_element_type=$type\n";
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


