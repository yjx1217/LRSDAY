#!/usr/bin/perl
use warnings;
#use warnings FATAL => 'all';
use strict;
use Getopt::Long;

##############################################################
#  script: tidy_maker_gff3.pl
#  author: Jia-Xing Yue (GitHub ID: yjx1217)
#  last edited: 2022.12.09
#  description: tidy maker gff3 by re-sorting and re-naming feature IDs
#  example: perl tidy_maker_gff3.pl -i raw.gff3 -t genome_tag -o tidy.gff3 -r genome.fa(.gz)
##############################################################

my ($refseq, $input, $output, $tag);
GetOptions('refseq|r:s' => \$refseq, # reference genome file in fasta format; can be compressed by gzip
	   'input|i:s' => \$input, # input raw gff3
	   'output|o:s' => \$output, # output tidy gff3
	   'tag|t:s' => \$tag); # genome tag

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

my $feature_index = 0;
foreach my $chr (@refseq) { 
    my @features_on_chr = ();
    foreach my $feature_id (sort keys %gff) {
	if ($gff{$feature_id}{'chr'} eq $chr) {
	    push @features_on_chr, $feature_id;
	}
    }
    foreach my $feature_id (sort {$gff{$a}{'start'} <=> $gff{$b}{'start'} or $gff{$a}{'end'} <=> $gff{$b}{'end'}} @features_on_chr) {
	my $feature_name = $gff{$feature_id}{'feature_name'};
	my $feature_type = $gff{$feature_id}{'feature_type'};
	my $feature_start = $gff{$feature_id}{'start'};
	my $feature_end = $gff{$feature_id}{'end'};
	my $feature_score = $gff{$feature_id}{'score'};
	my $feature_strand = $gff{$feature_id}{'strand'};
	my $feature_phase = $gff{$feature_id}{'phase'};
	if ($feature_type ne "gene") {
	    print $output_fh "$chr\t$tag\t$feature_type\t$feature_start\t$feature_end\t$feature_score\t$feature_strand\t$feature_phase\tID=$feature_id;Name=$feature_name\n";
	} else {
	    $feature_index += 10;
	    my $new_feature_id = sprintf("%07d", $feature_index);
	    $new_feature_id = "${tag}_" . "G" . $new_feature_id;
	    print $output_fh "$chr\t$tag\t$feature_type\t$feature_start\t$feature_end\t$feature_score\t$feature_strand\t$feature_phase\tID=$new_feature_id;Name=$new_feature_id\n";
	    foreach my $mRNA_id (sort {$gff{$feature_id}{'mRNA'}{$a}{'mRNA_index'} <=> $gff{$feature_id}{'mRNA'}{$b}{'mRNA_index'}} keys %{$gff{$feature_id}{'mRNA'}}) {
		my $mRNA_index = $gff{$feature_id}{'mRNA'}{$mRNA_id}{'mRNA_index'};
		my $mRNA_start = $gff{$feature_id}{'mRNA'}{$mRNA_id}{'start'};
		my $mRNA_end = $gff{$feature_id}{'mRNA'}{$mRNA_id}{'end'};
		my $mRNA_score = $gff{$feature_id}{'mRNA'}{$mRNA_id}{'score'};
		my $mRNA_strand = $gff{$feature_id}{'mRNA'}{$mRNA_id}{'strand'};
		my $mRNA_phase = $gff{$feature_id}{'mRNA'}{$mRNA_id}{'phase'};
		my $new_mRNA_id = "$new_feature_id.mRNA.$mRNA_index";
		print $output_fh "$chr\t$tag\tmRNA\t$mRNA_start\t$mRNA_end\t$mRNA_score\t$mRNA_strand\t$mRNA_phase\tID=$new_mRNA_id;Name=$new_mRNA_id;Parent=$new_feature_id\n";
		my @exon_indices = sort {$a <=> $b} keys %{$gff{$feature_id}{'mRNA'}{$mRNA_id}{'exon'}};
		my $exon_num = scalar @exon_indices;
		my $new_exon_index;
		if ($mRNA_strand eq '+') {
		    $new_exon_index = 0;
		} else {
		    $new_exon_index = $exon_num + 1;
		}
		foreach my $exon_index (@exon_indices) {
		    my $exon_start = $gff{$feature_id}{'mRNA'}{$mRNA_id}{'exon'}{$exon_index}{'start'};
		    my $exon_end = $gff{$feature_id}{'mRNA'}{$mRNA_id}{'exon'}{$exon_index}{'end'};
		    my $exon_score = $gff{$feature_id}{'mRNA'}{$mRNA_id}{'exon'}{$exon_index}{'score'};
		    my $exon_strand = $gff{$feature_id}{'mRNA'}{$mRNA_id}{'exon'}{$exon_index}{'strand'};
		    my $exon_phase = $gff{$feature_id}{'mRNA'}{$mRNA_id}{'exon'}{$exon_index}{'phase'};
		    if ($mRNA_strand eq '+') {
			$new_exon_index++;
		    } else {
			$new_exon_index--;
		    }
		    my $new_exon_id = "$new_mRNA_id.exon.$new_exon_index";
		    print $output_fh "$chr\t$tag\texon\t$exon_start\t$exon_end\t$exon_score\t$exon_strand\t$exon_phase\tID=$new_exon_id;Name=$new_exon_id;Parent=$new_mRNA_id\n";
		}
		my @cds_indices = sort {$a <=> $b} keys %{$gff{$feature_id}{'mRNA'}{$mRNA_id}{'cds'}};
		my $cds_num = scalar @cds_indices;
		my $new_cds_index;
		if ($mRNA_strand eq '+') {
		    $new_cds_index = 0;
		} else {
		    $new_cds_index = $cds_num + 1;
		}
		foreach my $cds_index (@cds_indices) {
		    my $cds_start = $gff{$feature_id}{'mRNA'}{$mRNA_id}{'cds'}{$cds_index}{'start'};
		    my $cds_end = $gff{$feature_id}{'mRNA'}{$mRNA_id}{'cds'}{$cds_index}{'end'};
		    my $cds_score = $gff{$feature_id}{'mRNA'}{$mRNA_id}{'cds'}{$cds_index}{'score'};
		    my $cds_strand = $gff{$feature_id}{'mRNA'}{$mRNA_id}{'cds'}{$cds_index}{'strand'};
		    my $cds_phase = $gff{$feature_id}{'mRNA'}{$mRNA_id}{'cds'}{$cds_index}{'phase'};
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
    my $feature_id;
    my $feature_name;
    my $feature_type;
    my $mRNA_index;
    while (<$fh>) {
	chomp;
	/^##FASTA/ and last;
	/^#/ and next;
	my ($chr, $source, $type, $start, $end, $score, $strand, $phase, $attributes) = split /\t/, $_;
	# print "$_\n";
	if ($type eq "gene") {
	    ($feature_id, $feature_name) = ($attributes =~ /ID=([^;]+);\S*Name=([^;]+)/);
	    if ($feature_id =~ /^trnascan/) {
		$feature_type = "tRNA";
		next;
	    } elsif ($feature_id =~ /^snoscan/) {
		$feature_type = "snoRNA";
		next;
	    } else {
		$feature_type = "gene";
	    }
	    $gff{$feature_id}{'feature_type'} = $feature_type;
	    $gff{$feature_id}{'feature_name'} = $feature_name;
	    $gff{$feature_id}{'chr'} = $chr;
	    $gff{$feature_id}{'start'} = $start;
	    $gff{$feature_id}{'end'} = $end;
	    $gff{$feature_id}{'strand'} = $strand;
	    $gff{$feature_id}{'source'} = $source;
	    $gff{$feature_id}{'score'} = $score;
	    $gff{$feature_id}{'phase'} = $phase;
	    $mRNA_index = 0;
	} elsif ($type eq "centromere") {
	    $feature_type = $type;
	    ($feature_id, $feature_name) = ($attributes =~ /ID=([^;]+);\S*Name=([^;]+)/);
	    $gff{$feature_id}{'feature_type'} = $feature_type;
	    $gff{$feature_id}{'feature_name'} = $feature_name;
            $gff{$feature_id}{'chr'} = $chr;
            $gff{$feature_id}{'start'} = $start;
            $gff{$feature_id}{'end'} = $end;
            $gff{$feature_id}{'strand'} = $strand;
            $gff{$feature_id}{'source'} = $source;
            $gff{$feature_id}{'score'} = $score;
            $gff{$feature_id}{'phase'} = $phase;
	# } elsif ($type eq "mobile_element") {
	#     $feature_type = $type;
	#     ($feature_id, $feature_name) = ($attributes =~ /ID=([^;]+);\S*Name=([^;]+\S*mobile_element_type=[^;]+)/);
	#     $gff{$feature_id}{'feature_type'} = $feature_type;
	#     $gff{$feature_id}{'feature_name'} = $feature_name;
        #     $gff{$feature_id}{'chr'} = $chr;
        #     $gff{$feature_id}{'start'} = $start;
        #     $gff{$feature_id}{'end'} = $end;
        #     $gff{$feature_id}{'strand'} = $strand;
        #     $gff{$feature_id}{'source'} = $source;
        #     $gff{$feature_id}{'score'} = $score;
        #     $gff{$feature_id}{'phase'} = $phase;
	} elsif ($type eq "tRNA") {
	    $feature_type = $type;
	    ($feature_id, $feature_name) = ($attributes =~ /ID=([^;]+);\S*Name=([^;]+)/);
	    # print "feature_id=$feature_id\n";
	    # print "feature_name=$feature_name\n";
	    my ($aa, $anticodon) = ($feature_id =~ /noncoding\-([^\_]+)\_([^\_]+)\-gene/);
	    # print "aa=$aa\n";
	    # print "anticodon=$anticodon\n";
	    $feature_id = "$feature_type:$chr:${start}-${end}:$strand";
	    $feature_name = "tRNA_${aa}($anticodon)";
	    # print "feature_id=$feature_id, feature_name=$feature_name";
	    $gff{$feature_id}{'feature_type'} = $feature_type;
	    $gff{$feature_id}{'feature_name'} = $feature_name;
            $gff{$feature_id}{'chr'} = $chr;
            $gff{$feature_id}{'start'} = $start;
            $gff{$feature_id}{'end'} = $end;
            $gff{$feature_id}{'strand'} = $strand;
            $gff{$feature_id}{'source'} = $source;
            $gff{$feature_id}{'score'} = $score;
            $gff{$feature_id}{'phase'} = $phase;
	} elsif ($type !~ /(exon|CDS|mRNA|UTR)/) {
	    # e.g. type = TY, X-element, Y_prime_element ...
	    $feature_type = $type;
	    $feature_id = "$feature_type:$chr:${start}-${end}:$strand";
	    $feature_name = "$feature_type:$chr:${start}-${end}:$strand";
	    $gff{$feature_id}{'feature_type'} = $feature_type;
	    $gff{$feature_id}{'feature_name'} = $feature_name;
	    $gff{$feature_id}{'chr'} = $chr;
	    $gff{$feature_id}{'start'} = $start;
	    $gff{$feature_id}{'end'} = $end;
	    $gff{$feature_id}{'strand'} = $strand;
	    $gff{$feature_id}{'source'} = $source;
	    $gff{$feature_id}{'score'} = $score;
	    $gff{$feature_id}{'phase'} = $phase;
	} elsif ($type eq "mRNA") {
	    my ($mRNA_id, $feature_id) = ($attributes =~ /ID=([^;]+);\S*Parent=([^;]+)/);
	    if (exists $gff{$feature_id}) {
		$mRNA_index++;
		$gff{$feature_id}{'mRNA'}{$mRNA_id}{'mRNA_index'} = $mRNA_index;
		$gff{$feature_id}{'mRNA'}{$mRNA_id}{'chr'} = $chr;
		$gff{$feature_id}{'mRNA'}{$mRNA_id}{'start'} = $start;
		$gff{$feature_id}{'mRNA'}{$mRNA_id}{'end'} = $end;
		$gff{$feature_id}{'mRNA'}{$mRNA_id}{'strand'} = $strand;
		$gff{$feature_id}{'mRNA'}{$mRNA_id}{'source'} = $source;
		$gff{$feature_id}{'mRNA'}{$mRNA_id}{'score'} = $score;
		$gff{$feature_id}{'mRNA'}{$mRNA_id}{'phase'} = $phase;
	    } else {
		die "cannot find matching gene record for the mRNA $mRNA_id derived from the gene $feature_id\n"
	    }
	} elsif ($type eq "exon") {
	    my ($exon_id, $mRNA_id) = ($attributes =~ /ID=([^;]+);\S*Parent=([^;]+)/);
	    if ($exon_id =~ /^trnascan/) {
		next;
	    } elsif ($exon_id =~ /^snoscan/) {
		next;
	    } elsif (exists $gff{$feature_id}{'mRNA'}{$mRNA_id}) {
		my $exon_index = $start;
		$gff{$feature_id}{'mRNA'}{$mRNA_id}{'exon'}{$exon_index}{'chr'} = $chr;
		$gff{$feature_id}{'mRNA'}{$mRNA_id}{'exon'}{$exon_index}{'start'} = $start;
		$gff{$feature_id}{'mRNA'}{$mRNA_id}{'exon'}{$exon_index}{'end'} = $end;
		$gff{$feature_id}{'mRNA'}{$mRNA_id}{'exon'}{$exon_index}{'strand'} = $strand;
		$gff{$feature_id}{'mRNA'}{$mRNA_id}{'exon'}{$exon_index}{'source'} = $source;
		$gff{$feature_id}{'mRNA'}{$mRNA_id}{'exon'}{$exon_index}{'score'} = $score;
		$gff{$feature_id}{'mRNA'}{$mRNA_id}{'exon'}{$exon_index}{'phase'} = $phase;
	    } else {
		die "cannot find matching mRNA record for the exon $exon_id derived from the mRNA $mRNA_id\n"
	    }
	} elsif ($type eq "CDS") {
            my ($cds_id, $mRNA_id) = ($attributes =~ /ID=([^;]+);\S*Parent=([^;]+)/);
            if (exists $gff{$feature_id}{'mRNA'}{$mRNA_id}) {
                my $cds_index = $start;
                $gff{$feature_id}{'mRNA'}{$mRNA_id}{'cds'}{$cds_index}{'chr'} = $chr;
                $gff{$feature_id}{'mRNA'}{$mRNA_id}{'cds'}{$cds_index}{'start'} = $start;
                $gff{$feature_id}{'mRNA'}{$mRNA_id}{'cds'}{$cds_index}{'end'} = $end;
                $gff{$feature_id}{'mRNA'}{$mRNA_id}{'cds'}{$cds_index}{'strand'} = $strand;
                $gff{$feature_id}{'mRNA'}{$mRNA_id}{'cds'}{$cds_index}{'source'} = $source;
                $gff{$feature_id}{'mRNA'}{$mRNA_id}{'cds'}{$cds_index}{'score'} = $score;
                $gff{$feature_id}{'mRNA'}{$mRNA_id}{'cds'}{$cds_index}{'phase'} = $phase;
            } else {
                die "cannot find matching mRNA record for the CDS $cds_id derived from the mRNA $mRNA_id\n"
            }
	}
    }
    return %gff;
}

