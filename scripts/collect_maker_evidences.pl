#!/usr/bin/perl
use warnings;
use strict;
use Getopt::Long;

##############################################################
#  script: collect_maker_evidences.pl
#  author: Jia-Xing Yue (GitHub ID: yjx1217)
#  last edited: 2017.06.17
#  description: collect MAKER annotation evidences (EST and protein alignment in GFF3 format).
#  example: perl collect_maker_evidences.pl -t genome_tag
##############################################################


my ($tag);
GetOptions('tag|t:s' => \$tag); # genome_tag 

my $index = "\.\/$tag.maker.output\/${tag}_master_datastore_index\.log";
my $index_fh = read_file($index);

my %data = ();
my $output1 = "$tag.est_evidence.gff3";
my $output1_fh = write_file($output1);
my $output2 = "$tag.protein_evidence.gff3";
my $output2_fh = write_file($output2);


while (<$index_fh>) {
    chomp;
    /^#/ and next;
    /^\s*$/ and next;
    my ($chr, $dir, $stat) = split '\t', $_;
    if (not exists $data{$chr}) {
	$data{$chr}{'dir'} = $dir;
	my $evidence_dir = "\.\/$tag.maker.output\/${dir}theVoid\.$chr";
	my @evidence = glob "$evidence_dir\/evidence_*.gff";
	foreach my $evidence (@evidence) {
	    my $fh = read_file($evidence);
	    while (<$fh>) {
		chomp;
		/^#/ and next;
		my @line = split '\t', $_;
		if ($line[1] eq 'est2genome') {
		#if (($line[1] eq 'est2genome') and ($line[2] eq 'match_part')) {
		    print $output1_fh "$_\n";
		}
		if ($line[1] eq 'protein2genome') {
		#if (($line[1] eq 'protein2genome') and ($line[2] eq 'match_part')) {
		    print $output2_fh "$_\n";
		}
	    }
	}
    }
}



sub read_file {
    my $filename = shift @_;
    my $fh;
    open ($fh, $filename) or die "cannot open the file $filename\n";
    return $fh;
}


sub write_file {
    my $filename = shift @_;
    my $fh;
    open ($fh, ">$filename") or die "cannot open the file $filename\n";
    return $fh;
}  
