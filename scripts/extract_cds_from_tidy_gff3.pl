#!/usr/bin/perl
use warnings;
use strict;
use Getopt::Long;

##############################################################
#  script: extract_cds_from_tidy_gff3.pl
#  author: Jia-Xing Yue (GitHub ID: yjx1217)
#  last edited: 2017.06.17
#  description: extract CDS sequences from GFF3 file
#  example: perl extract_cds_from_tidy_gff3.pl -r genome.fa.gz -g genome.gff3 -o genome.cds.fa
##############################################################

my ($refseq, $gff, $output);
GetOptions('refseq|r:s' => \$refseq, # reference genome
	   'gff|g:s' => \$gff, # input GFF3 file
	   'output|o:s' => \$output); # output file

my $refseq_fh = read_file($refseq);
my %refseq = ();
my @refseq = ();
parse_fasta_file($refseq_fh, \%refseq, \@refseq);

my $gff_fh = read_file($gff);
my %gff = parse_gff_file($gff_fh);
close $gff_fh;

my $output_fh = write_file($output);


foreach my $chr (@refseq) { 
    my @genes_on_chr = ();
    foreach my $gene_id (sort keys %gff) {
	if (($gff{$gene_id}{'gene_chr'} eq $chr) and ($gff{$gene_id}{'gene_type'} eq 'gene')) {
	    push @genes_on_chr, $gene_id;
	}
    }
    foreach my $gene_id (sort {$gff{$a}{'gene_start'} <=> $gff{$b}{'gene_start'} or $gff{$a}{'gene_end'} <=> $gff{$b}{'gene_end'}} @genes_on_chr) {
	foreach my $mRNA_id (sort {$gff{$gene_id}{'mRNA'}{$a}{'mRNA_index'} <=> $gff{$gene_id}{'mRNA'}{$b}{'mRNA_index'}} keys %{$gff{$gene_id}{'mRNA'}}) {
	    my $mRNA_index = $gff{$gene_id}{'mRNA'}{$mRNA_id}{'mRNA_index'};
	    my $mRNA_strand = $gff{$gene_id}{'mRNA'}{$mRNA_id}{'mRNA_strand'};
	    my $mRNA_seq = "";
	    my @cds_indices = sort {$a <=> $b} keys %{$gff{$gene_id}{'mRNA'}{$mRNA_id}{'cds'}};
	    foreach my $cds_index (@cds_indices) {
		my $cds_start = $gff{$gene_id}{'mRNA'}{$mRNA_id}{'cds'}{$cds_index}{'cds_start'};
		my $cds_end = $gff{$gene_id}{'mRNA'}{$mRNA_id}{'cds'}{$cds_index}{'cds_end'};
		my $cds_seq = substr $refseq{$chr}, $cds_start - 1, $cds_end - $cds_start + 1;
		if ($mRNA_strand eq '+') {
		    $mRNA_seq = $mRNA_seq . $cds_seq;
		} else {
		    $mRNA_seq = $cds_seq . $mRNA_seq;
		}
	    }
	    
	    if ($mRNA_strand eq '-') {
		$mRNA_seq = revcom($mRNA_seq);
	    }
	    print $output_fh ">$gene_id|$mRNA_id\n$mRNA_seq\n";
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

sub parse_gff_file {
    my $fh = shift @_;
    my %gff = ();
    my $gene_id;
    my $gene_name;
    my $gene_type;
    my $mRNA_index;
    while (<$fh>) {
	chomp;
	/^##FASTA/ and last;
	/^#/ and next;
	my ($chr, $source, $type, $start, $end, $score, $strand, $phase, $attributes) = split /\t/, $_;
	if ($type eq "gene") {
	    ($gene_id, $gene_name) = ($attributes =~ /ID=([^;]+);\S*Name=([^;]+)/);
	    if ($gene_id =~ /^trnascan/) {
		$gene_type = "tRNA";
		next;
	    } elsif ($gene_id =~ /^snoscan/) {
		$gene_type = "snoRNA";
		next;
	    } else {
		$gene_type = "gene";
	    }
	    $gff{$gene_id}{'gene_type'} = $gene_type;
	    $gff{$gene_id}{'gene_chr'} = $chr;
	    $gff{$gene_id}{'gene_start'} = $start;
	    $gff{$gene_id}{'gene_end'} = $end;
	    $gff{$gene_id}{'gene_strand'} = $strand;
	    $gff{$gene_id}{'gene_source'} = $source;
	    $gff{$gene_id}{'gene_score'} = $score;
	    $gff{$gene_id}{'gene_phase'} = $phase;
	    $mRNA_index = 0;
	} elsif ($type !~ /(exon|CDS|mRNA|UTR)/) {
	    # e.g. type = centromere, TY, X-element, Y_prime_element ...
	    $gene_type = $type;
	    $gene_id = "$gene_type:$chr:${start}-${end}:$strand";
	    $gff{$gene_id}{'gene_type'} = $gene_type;
	    $gff{$gene_id}{'gene_type'} = $type;
	    $gff{$gene_id}{'gene_chr'} = $chr;
	    $gff{$gene_id}{'gene_start'} = $start;
	    $gff{$gene_id}{'gene_end'} = $end;
	    $gff{$gene_id}{'gene_strand'} = $strand;
	    $gff{$gene_id}{'gene_source'} = $source;
	    $gff{$gene_id}{'gene_score'} = $score;
	    $gff{$gene_id}{'gene_phase'} = $phase;
	} elsif ($type eq 'mRNA') {
	    my ($mRNA_id, $gene_id) = ($attributes =~ /ID=([^;]+);\S*Parent=([^;]+)/);
	    if (exists $gff{$gene_id}) {
		$mRNA_index++;
		$gff{$gene_id}{'mRNA'}{$mRNA_id}{'mRNA_index'} = $mRNA_index;
		$gff{$gene_id}{'mRNA'}{$mRNA_id}{'mRNA_chr'} = $chr;
		$gff{$gene_id}{'mRNA'}{$mRNA_id}{'mRNA_start'} = $start;
		$gff{$gene_id}{'mRNA'}{$mRNA_id}{'mRNA_end'} = $end;
		$gff{$gene_id}{'mRNA'}{$mRNA_id}{'mRNA_strand'} = $strand;
		$gff{$gene_id}{'mRNA'}{$mRNA_id}{'mRNA_source'} = $source;
		$gff{$gene_id}{'mRNA'}{$mRNA_id}{'mRNA_score'} = $score;
		$gff{$gene_id}{'mRNA'}{$mRNA_id}{'mRNA_phase'} = $phase;
	    } else {
		die "cannot find matching gene record for the mRNA $mRNA_id derived from the gene $gene_id\n"
	    }
	} elsif ($type eq 'exon') {
	    my ($exon_id, $mRNA_id) = ($attributes =~ /ID=([^;]+);\S*Parent=([^;]+)/);
	    if ($exon_id =~ /^trnascan/) {
		next;
	    } elsif ($exon_id =~ /^snoscan/) {
		next;
	    } elsif (exists $gff{$gene_id}{'mRNA'}{$mRNA_id}) {
		my $exon_index = $start;
		$gff{$gene_id}{'mRNA'}{$mRNA_id}{'exon'}{$exon_index}{'exon_chr'} = $chr;
		$gff{$gene_id}{'mRNA'}{$mRNA_id}{'exon'}{$exon_index}{'exon_start'} = $start;
		$gff{$gene_id}{'mRNA'}{$mRNA_id}{'exon'}{$exon_index}{'exon_end'} = $end;
		$gff{$gene_id}{'mRNA'}{$mRNA_id}{'exon'}{$exon_index}{'exon_strand'} = $strand;
		$gff{$gene_id}{'mRNA'}{$mRNA_id}{'exon'}{$exon_index}{'exon_source'} = $source;
		$gff{$gene_id}{'mRNA'}{$mRNA_id}{'exon'}{$exon_index}{'exon_score'} = $score;
		$gff{$gene_id}{'mRNA'}{$mRNA_id}{'exon'}{$exon_index}{'exon_phase'} = $phase;
	    } else {
		die "cannot find matching mRNA record for the exon $exon_id derived from the mRNA $mRNA_id\n"
	    }
	} elsif ($type eq 'CDS') {
            my ($cds_id, $mRNA_id) = ($attributes =~ /ID=([^;]+);\S*Parent=([^;]+)/);
            if (exists $gff{$gene_id}{'mRNA'}{$mRNA_id}) {
                my $cds_index = $start;
                $gff{$gene_id}{'mRNA'}{$mRNA_id}{'cds'}{$cds_index}{'cds_chr'} = $chr;
                $gff{$gene_id}{'mRNA'}{$mRNA_id}{'cds'}{$cds_index}{'cds_start'} = $start;
                $gff{$gene_id}{'mRNA'}{$mRNA_id}{'cds'}{$cds_index}{'cds_end'} = $end;
                $gff{$gene_id}{'mRNA'}{$mRNA_id}{'cds'}{$cds_index}{'cds_strand'} = $strand;
                $gff{$gene_id}{'mRNA'}{$mRNA_id}{'cds'}{$cds_index}{'cds_source'} = $source;
                $gff{$gene_id}{'mRNA'}{$mRNA_id}{'cds'}{$cds_index}{'cds_score'} = $score;
                $gff{$gene_id}{'mRNA'}{$mRNA_id}{'cds'}{$cds_index}{'cds_phase'} = $phase;
            } else {
                die "cannot find matching mRNA record for the CDS $cds_id derived from the mRNA $mRNA_id\n"
            }
	}
    }
    return %gff;
}


sub revcom {
    my $seq = shift @_;
    my $seq_revcom = reverse $seq;
    $seq_revcom =~ tr/ATGCatgc/TACGtacg/;
    return $seq_revcom;
}
