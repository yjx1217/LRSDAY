#!/usr/bin/perl
use warnings;
use strict;
use Getopt::Long;

##############################################################
#  script: summarize_pilon_correction.pl
#  author: Jia-Xing Yue (GitHub ID: yjx1217)
#  last edited: 2018.07.11
#  description: summarize the total number of SNP, insertion, and deletion corrections made by Pilon
#  example: perl summarize_pilon_correction.pl -i input.pilon.changes
##############################################################

my ($input);

GetOptions('input|i:s' => \$input); # input pilon changes file

my $input_fh = read_file($input);
my $snp_count = 0;
my $insertion_count = 0;
my $deletion_count = 0;

while (<$input_fh>) {
    chomp;
    /^#/ and next;
    /^\s*$/ and next;
    my ($original_pos, $new_pos, $original_base, $new_base) = split /\s+/, $_;
    if (($original_base ne '.') and ($new_base ne '.')) {
	$snp_count += (length $original_base);
    } elsif ($original_base eq '.') {
	$insertion_count++;
    } else {
	$deletion_count++;
    }
}
my $indel_count = $insertion_count + $deletion_count;

print "\n";
print "### Pilon correction summary ###\n";
print "total number of SNP correction: $snp_count\n";
print "total number of INDEL correction: $indel_count\n";
print "  including $insertion_count insertions and $deletion_count deletions relative to the original assembly\n\n";
print "\n";

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
