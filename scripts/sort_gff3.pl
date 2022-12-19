#!/usr/bin/perl
use warnings;
use strict;
use Getopt::Long;

##############################################################
#  script: sort_gff3.pl
#  author: Jia-Xing Yue (GitHub ID: yjx1217)
#  last edited: 2022.12.11
#  description: sort gff3 file by the input multi-fasta genome file
#  example: perl sort_gff3.pl -i raw.gff3 -t genome_tag -o sorted.gff3 -r genome.fa(.gz)
##############################################################

my ($refseq, $input, $output, $tag);
GetOptions('refseq|r:s' => \$refseq,
	   'input|i:s' => \$input,
	   'output|o:s' => \$output,
	   'tag|t:s' => \$tag);

my $refseq_fh = read_file($refseq);
my %refseq = ();
my @refseq = ();
parse_fasta_file($refseq_fh, \%refseq, \@refseq);

my $input_fh = read_file($input);
my %gff = parse_gff_file($input_fh);
close $input_fh;

my $output_fh = write_file($output);

print $output_fh "##gff-version 3\n";

foreach my $chr (@refseq) {
    my $chr_length = length $refseq{$chr};
    print $output_fh "##sequence-region $chr 1 $chr_length\n";
}

my $gene_index = 0;
foreach my $chr (@refseq) { 
    my @genes_on_chr = ();
    foreach my $gene_id (sort keys %gff) {
	if ($gff{$gene_id}{'gene_chr'} eq $chr) {
	    push @genes_on_chr, $gene_id;
	}
    }
    foreach my $gene_id (sort {$gff{$a}{'gene_start'} <=> $gff{$b}{'gene_start'} or $gff{$a}{'gene_end'} <=> $gff{$b}{'gene_end'}} @genes_on_chr) {
	my $gene_name = $gff{$gene_id}{'gene_name'};
	my $gene_type = $gff{$gene_id}{'gene_type'};
	my $gene_start = $gff{$gene_id}{'gene_start'};
	my $gene_end = $gff{$gene_id}{'gene_end'};
	my $gene_score = $gff{$gene_id}{'gene_score'};
	my $gene_strand = $gff{$gene_id}{'gene_strand'};
	my $gene_phase = $gff{$gene_id}{'gene_phase'};
	print $output_fh "$chr\t$tag\t$gene_type\t$gene_start\t$gene_end\t$gene_score\t$gene_strand\t$gene_phase\tID=$gene_id;Name=$gene_name\n";
	if ($gene_type eq "gene") {
	    foreach my $mRNA_id (sort {$gff{$gene_id}{'mRNA'}{$a}{'mRNA_index'} <=> $gff{$gene_id}{'mRNA'}{$b}{'mRNA_index'}} keys %{$gff{$gene_id}{'mRNA'}}) {
		my $mRNA_index = $gff{$gene_id}{'mRNA'}{$mRNA_id}{'mRNA_index'};
		my $mRNA_start = $gff{$gene_id}{'mRNA'}{$mRNA_id}{'mRNA_start'};
		my $mRNA_end = $gff{$gene_id}{'mRNA'}{$mRNA_id}{'mRNA_end'};
		my $mRNA_score = $gff{$gene_id}{'mRNA'}{$mRNA_id}{'mRNA_score'};
		my $mRNA_strand = $gff{$gene_id}{'mRNA'}{$mRNA_id}{'mRNA_strand'};
		my $mRNA_phase = $gff{$gene_id}{'mRNA'}{$mRNA_id}{'mRNA_phase'};
		my $new_mRNA_id = "$gene_id.mRNA.$mRNA_index";
		print $output_fh "$chr\t$tag\tmRNA\t$mRNA_start\t$mRNA_end\t$mRNA_score\t$mRNA_strand\t$mRNA_phase\tID=$new_mRNA_id;Name=$new_mRNA_id;Parent=$gene_id\n";
		my @exon_indices = sort {$a <=> $b} keys %{$gff{$gene_id}{'mRNA'}{$mRNA_id}{'exon'}};
		my $exon_num = scalar @exon_indices;
		my $new_exon_index;
		if ($mRNA_strand eq '+') {
		    $new_exon_index = 0;
		} else {
		    $new_exon_index = $exon_num + 1;
		}
		foreach my $exon_index (@exon_indices) {
		    my $exon_start = $gff{$gene_id}{'mRNA'}{$mRNA_id}{'exon'}{$exon_index}{'exon_start'};
		    my $exon_end = $gff{$gene_id}{'mRNA'}{$mRNA_id}{'exon'}{$exon_index}{'exon_end'};
		    my $exon_score = $gff{$gene_id}{'mRNA'}{$mRNA_id}{'exon'}{$exon_index}{'exon_score'};
		    my $exon_strand = $gff{$gene_id}{'mRNA'}{$mRNA_id}{'exon'}{$exon_index}{'exon_strand'};
		    my $exon_phase = $gff{$gene_id}{'mRNA'}{$mRNA_id}{'exon'}{$exon_index}{'exon_phase'};
		    if ($mRNA_strand eq '+') {
			$new_exon_index++;
		    } else {
			$new_exon_index--;
		    }
		    my $new_exon_id = "$new_mRNA_id.exon.$new_exon_index";
		    print $output_fh "$chr\t$tag\texon\t$exon_start\t$exon_end\t$exon_score\t$exon_strand\t$exon_phase\tID=$new_exon_id;Name=$new_exon_id;Parent=$new_mRNA_id\n";
		}
		my @cds_indices = sort {$a <=> $b} keys %{$gff{$gene_id}{'mRNA'}{$mRNA_id}{'cds'}};
		my $cds_num = scalar @cds_indices;
		my $new_cds_index;
		if ($mRNA_strand eq '+') {
		    $new_cds_index = 0;
		} else {
		    $new_cds_index = $cds_num + 1;
		}
		foreach my $cds_index (@cds_indices) {
		    my $cds_start = $gff{$gene_id}{'mRNA'}{$mRNA_id}{'cds'}{$cds_index}{'cds_start'};
		    my $cds_end = $gff{$gene_id}{'mRNA'}{$mRNA_id}{'cds'}{$cds_index}{'cds_end'};
		    my $cds_score = $gff{$gene_id}{'mRNA'}{$mRNA_id}{'cds'}{$cds_index}{'cds_score'};
		    my $cds_strand = $gff{$gene_id}{'mRNA'}{$mRNA_id}{'cds'}{$cds_index}{'cds_strand'};
		    my $cds_phase = $gff{$gene_id}{'mRNA'}{$mRNA_id}{'cds'}{$cds_index}{'cds_phase'};
		    if ($mRNA_strand eq '+') {
			$new_cds_index++;
		    } else {
			$new_cds_index--;
		    }
		    my $new_cds_id = "$new_mRNA_id.CDS.$new_cds_index";
		    print $output_fh "$chr\t$tag\tCDS\t$cds_start\t$cds_end\t$cds_score\t$cds_strand\t$cds_phase\tID=$new_cds_id;Name=$new_cds_id;Parent=$new_mRNA_id\n";
		}
	    }
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
	    $gff{$gene_id}{'gene_name'} = $gene_name;
	    $gff{$gene_id}{'gene_type'} = $gene_type;
	    $gff{$gene_id}{'gene_chr'} = $chr;
	    $gff{$gene_id}{'gene_start'} = $start;
	    $gff{$gene_id}{'gene_end'} = $end;
	    $gff{$gene_id}{'gene_strand'} = $strand;
	    $gff{$gene_id}{'gene_source'} = $source;
	    $gff{$gene_id}{'gene_score'} = $score;
	    $gff{$gene_id}{'gene_phase'} = $phase;
	    $mRNA_index = 0;
	} elsif ($type =~ /(centromere|rRNA|tRNA)/) {
	    $gene_type = $type;
	    ($gene_id, $gene_name) = ($attributes =~ /ID=([^;]+);\S*Name=([^;]+)/);
	    $gff{$gene_id}{'gene_name'} = $gene_name;
	    $gff{$gene_id}{'gene_type'} = $gene_type;
            $gff{$gene_id}{'gene_chr'} = $chr;
            $gff{$gene_id}{'gene_start'} = $start;
            $gff{$gene_id}{'gene_end'} = $end;
            $gff{$gene_id}{'gene_strand'} = $strand;
            $gff{$gene_id}{'gene_source'} = $source;
            $gff{$gene_id}{'gene_score'} = $score;
            $gff{$gene_id}{'gene_phase'} = $phase;
	} elsif ($type !~ /(exon|CDS|mRNA|UTR)/) {
	    # e.g. type = TY, X_element, Y_prime_element ...
	    $gene_type = $type;
	    $gene_id = "$gene_type:$chr:${start}-${end}:$strand";
	    $gff{$gene_id}{'gene_name'} = $gene_id;
	    $gff{$gene_id}{'gene_type'} = $gene_type;
	    $gff{$gene_id}{'gene_chr'} = $chr;
	    $gff{$gene_id}{'gene_start'} = $start;
	    $gff{$gene_id}{'gene_end'} = $end;
	    $gff{$gene_id}{'gene_strand'} = $strand;
	    $gff{$gene_id}{'gene_source'} = $source;
	    $gff{$gene_id}{'gene_score'} = $score;
	    $gff{$gene_id}{'gene_phase'} = $phase;
	} elsif ($type eq "mRNA") {
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
	} elsif ($type eq "exon") {
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
	} elsif ($type eq "CDS") {
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

