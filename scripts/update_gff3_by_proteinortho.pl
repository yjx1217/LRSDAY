#!/usr/bin/perl
use warnings;
use strict;
use Getopt::Long;

##############################################################
#  script: update_gff3_by_proteinortho.pl
#  author: Jia-Xing Yue (GitHub ID: yjx1217)
#  last edited: 2017.06.17
#  description: update gff3 file based on proteinortho result: attach gene name based on gene orthology 
#  example: perl update_gff3_by_proteinortho.pl -i query.maker.gff3 -x query.poff -r ref.PoFF.faa -q query.maker.PoFF.faa -o query.final.gff
##############################################################

my ($input, $orthology, $output, $ref, $query);

GetOptions('input|i:s' => \$input, # input annotation in GFF3 format
	   'orthology|x:s' => \$orthology, # proteinortho output
	   'ref|r:s' => \$ref, # reference genome file name showed in proteinortho output e.g. ref.PoFF.faa 
	   'query|q:s' => \$query, # query genome file name showed in proteinortho output e.g. SK1.PoFF.faa 
	   'output|o:s' => \$output);  # output prefix

my $orthology_fh = read_file($orthology);
my %orthology = parse_orthology_file($orthology_fh);
my %ortho_map = ();
my $input_fh = read_file($input);
my $output_fh = write_file($output);

foreach my $index (sort {$a<=>$b} keys %orthology) {
    # print "index=$index\n";
    if (($orthology{$index}{$query} ne '*') and ($orthology{$index}{$ref} ne '*')) {
	my @query_gene = ();
	my @ref_gene = ();
	# print "query_gene_array: $orthology{$index}{$query}\n";
	# print "ref_gene_array: $orthology{$index}{$ref}\n";

	if ($orthology{$index}{$query} =~ /,/) {
	    @query_gene = split /,/, $orthology{$index}{$query};
	} else {
	    @query_gene = ($orthology{$index}{$query});
	}

	if ($orthology{$index}{$ref} =~ /,/) {
            @ref_gene = split /,/, $orthology{$index}{$ref};
        } else {
            @ref_gene = ($orthology{$index}{$ref});
        }

	my $ref_gene;
	if ((scalar @ref_gene) > 1) {
	    $ref_gene = join '/', @ref_gene;
	} else {
	    $ref_gene = $ref_gene[0];
	}

	foreach my $query_gene (@query_gene) {
	    $ortho_map{$query_gene} = $ref_gene;
	}
    }
}

while (<$input_fh>) {
    chomp;
    /^#/ and next;
    /^\s*$/ and next;
    my ($chr, $source, $type, $start, $end, $score, $strand, $phase, $attributes) = split '\t', $_;
    my ($id, $name) = ($attributes =~ /ID=([^;]+);\S*Name=([^;]+)/);
    if (($type eq "gene") and (exists $ortho_map{$id})) {
	$name = $ortho_map{$id};
        print $output_fh "$chr\t$source\t$type\t$start\t$end\t$score\t$strand\t$phase\tID=$id;Name=$name\n";
    } else {
        print $output_fh "$_\n";
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

sub parse_orthology_file {
    my $fh = shift @_;
    my %orthology = ();
    my @sp_list = ();
    my $index = 0;
    while (<$fh>) {
	chomp;
	/^\s*$/ and next;
	my @line = split /\t/, $_;
	if (/^\#/) {
	    @sp_list = splice @line, 3;
	} else {
	    $index++;
	    my $sp_count = $line[0];
	    my $gene_count = $line[1];
	    my $connectivity = $line[2];
	    my @genes = splice @line, 3;
	    for (my $i = 0; $i < (scalar @sp_list); $i++) {
		$orthology{$index}{$sp_list[$i]} = $genes[$i];
		# print "index =$index, i=$i, sp=$sp_list[$i], gene=$genes[$i]\n";
	    }
	}
    }
    return %orthology;
}


