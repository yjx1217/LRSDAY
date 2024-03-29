#!/usr/bin/perl
use warnings;
use strict;
use Getopt::Long;

##############################################################
#  script: label_pseudogene_in_gff3.pl
#  author: Jia-Xing Yue (GitHub ID: yjx1217)
#  last edited: 2022.12.12
#  description: label pseudogenes in gff3 by manual_check_list generated by cds2protein.pl
#  example: perl label_pseudogene_in_gff3.pl -i input.gff3 -o output.gff3 -l check_list
##############################################################

my ($input, $output, $check_list);

GetOptions('input|i:s' => \$input,
	   'output|o:s' => \$output,
	   'check_list|l:s' => \$check_list);

my $list_fh = read_file($check_list);
my %list = parse_list_file($list_fh);
close $list_fh;

my $input_fh = read_file($input);
my $output_fh = write_file($output);

my $gene_id;
my $gene_name;
my $new_gene_id;
my $new_gene_name;
my $new_gene_type;
my $new_gene_attributes;

while (<$input_fh>) {
    chomp;
    if (/^#/) {
	print $output_fh "$_\n";
    } else {
	my ($chr, $source, $type, $start, $end, $score, $strand, $phase, $attributes) = split /\t/, $_;
        if ($type eq "gene") {
            ($gene_id, $gene_name) = ($attributes =~ /ID=([^;]+);\S*Name=([^;]+)/);
	    if (exists $list{'gene_id'}{$gene_id}) {
		$new_gene_type = "pseudogene";
		$new_gene_id = "${new_gene_type}:${chr}:${start}-${end}:$strand";
		$new_gene_name = $new_gene_id;
		$new_gene_attributes = "ID=$new_gene_id;Name=$new_gene_name";
		print $output_fh "$chr\t$source\t$new_gene_type\t$start\t$end\t$score\t$strand\t$phase\t$new_gene_attributes\n";
	    } else {
		print $output_fh "$_\n";
	    }
	} elsif ($type eq "mRNA") {
	    my ($mrna_id) = ($attributes =~ /ID=([^;]+)/);
            if (exists $list{'mrna_id'}{$mrna_id}) {
		my $new_mrna_type = "pseudogenic_transcript";
		my $new_mrna_id = "${new_mrna_type}:${chr}:${start}-${end}:$strand";
		my $new_mrna_name = $new_mrna_id;
		my $new_mrna_attributes = "ID=$new_mrna_id;Name=$new_mrna_name;Parent=$new_gene_id";
		print $output_fh "$chr\t$source\t$new_mrna_type\t$start\t$end\t$score\t$strand\t$phase\t$new_mrna_attributes\n";
            } else {
		print $output_fh "$_\n";
            }
	} elsif ($type =~ /(exon|CDS)/)  { 
	    my ($parent_id) = ($attributes =~ /Parent=([^;]+)/);
	    if (not exists $list{'mrna_id'}{$parent_id}) {
		print $output_fh "$_\n";
	    }
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

sub parse_list_file {
    my $fh = shift @_;
    my %list = ();
    while (<$fh>) {
	chomp;
	/^##FASTA/ and last;
	/^#/ and next;
	/^\s*$/ and next;
	my ($transcript, $cause) = split /\t/, $_;
	my ($gene_id, $transcript_id) = ($transcript =~ /([^\|]+)\|([^\|]+)/);
	print "gene_id=$gene_id, transcript_id = $transcript_id\n";
	$list{'gene_id'}{$gene_id} = $cause;
	$list{'mrna_id'}{$transcript_id} = $cause;

    }
    return %list;
}
