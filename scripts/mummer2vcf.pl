#!/usr/bin/perl
use warnings;
use strict;
use Getopt::Long;

##############################################################
#  script: mummer2vcf.pl
#  author: Jia-Xing Yue (GitHub ID: yjx1217)
#  last edited: 2017.09.11
#  description: convert mummer snps output to vcf format.
#  example: perl mummer2vcf.pl -i mummer.snps(.gz) -r genome.fa(.gz) -t SNP -p prefix
#  example: perl mummer2vcf.pl -i mummer.snps(.gz) -r genome.fa(.gz) -t INDEL -p prefix
##############################################################

my ($refseq, $input, $type, $prefix);
$type = "BOTH";

GetOptions('r|refseq:s' => \$refseq,
	   'i|input:s' => \$input,
	   't|type:s' => \$type, # case insensitive: SNP or INDEL or BOTH
	   'p|output:s' => \$prefix);


my $refseq_fh = read_file($refseq);
my %refseq = ();
my @refseq = ();
parse_fasta_file($refseq_fh, \%refseq, \@refseq);

my $input_fh = read_file($input);
my $output = "$prefix.mummer2vcf.$type.vcf";
my $output_fh = write_file($output);

print $output_fh "##fileformat=VCFv4.2\n";
print $output_fh "##source=mummer2vcf.pl\n";
print $output_fh "##reference=$refseq\n";
print $output_fh "##FILTER=<ID=PASS,Description=\"All filters passed\">\n";
print $output_fh "##FILTER=<ID=1to1_aln,Description=\"1to1_genome_alignment\">\n";
print $output_fh "##INFO=<ID=SAMPLE,Number=1 Type=String,Description=\"Samples id\">\n";
print $output_fh "##INFO=<ID=type,Number=1,Type=String,Description=\"variant type: SNP or INDEL\">\n";
print $output_fh "##INFO=<ID=ref_chr,Number=1,Type=String,Description=\"chr in the refrence\">\n";
print $output_fh "##INFO=<ID=ref_start,Number=1,Type=String,Description=\"start position in the refrence\">\n";
print $output_fh "##INFO=<ID=ref_end,Number=1,Type=String,Description=\"end position in the refrence\">\n";
print $output_fh "##INFO=<ID=query_chr,Number=1,Type=String,Description=\"chr in the query\">\n";
print $output_fh "##INFO=<ID=query_start,Number=1,Type=String,Description=\"start position in the query\">\n";
print $output_fh "##INFO=<ID=query_end,Number=1,Type=String,Description=\"end position in the query\">\n";
print $output_fh "##INFO=<ID=query_orientation,Number=1,Type=String,Description=\"matching direction in the query\">\n";
print $output_fh "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n";

my $indel_ref_chr;
my $indel_ref_start;
my $indel_ref_end;
my $indel_query_chr;
my $indel_query_start;
my $indel_query_end;
my $indel_type;
my $indel_ref_base;
my $indel_query_base;
my $indel_query_orientation;

my $snp_count = 0;
my $indel_count = 0;

while (<$input_fh>) {
    chomp;
    (1..4) and next;
    /^#/ and next;
    /^\s*$/ and next;
    my ($ref_pos, $ref_base, $query_base, $query_pos, $buff, $dist, $ref_length, $query_length, $ref_orientation, $query_orientation, $ref_chr, $query_chr) = split /\t/, $_;
    if (($ref_base ne ".") and ($query_base ne ".")) {
	# this is a SNP
	$snp_count++;
	if ($type =~ /(SNP|BOTH)/i) {
	    record_SNP($output_fh, $ref_chr, $ref_pos, $ref_base, $query_base, $query_chr, $query_pos, $query_orientation);
	} 
    } else {
	if (not defined $indel_type) {
	  # this is a new INDEL
	  NEW_INDEL:
	    if ($ref_base eq ".") {
		$indel_type = "INSERTION";
	    } else {
		$indel_type = "DELETION";
	    }
	    $indel_ref_chr = $ref_chr;
	    $indel_ref_start = $ref_pos;
	    $indel_ref_end = $ref_pos;
	    $indel_query_chr = $query_chr;
	    $indel_query_start = $query_pos;
	    $indel_query_end = $query_pos;
	    $indel_ref_base = $ref_base;
	    $indel_query_base = $query_base;
	    $indel_query_orientation = $query_orientation;
	    $indel_count++;
	} elsif (($indel_type eq "INSERTION") and ($ref_base eq ".")) {
	    if ((($indel_ref_chr eq $ref_chr) and ($indel_query_chr eq $query_chr)) and (($indel_ref_start == $ref_pos) and ($buff == 0))) {
		# this is the continuation of the previously encountered INSERTION
		$indel_query_end = $query_pos;
		if ($indel_query_orientation == 1) {
		    $indel_query_base = $indel_query_base . $query_base;
		} else {
		    $indel_query_base = $query_base . $indel_query_base;
		}
	    } else {
		# this is the start of a new INSERTION
		# record previous INSERTION
		if ($type =~ /(INDEL|BOTH)/i) {
		    record_INDEL($output_fh, $indel_type, $indel_ref_chr, $indel_ref_start, $indel_ref_end, $indel_query_chr, $indel_query_start, $indel_query_end, $indel_ref_base, $indel_query_base, $indel_query_orientation);
		}
		undef $indel_type;
		    goto NEW_INDEL;
	    }
	} elsif (($indel_type eq "DELETION") and ($query_base eq ".")) {
	    # this is the continuation of the previously encountered DELETION
	    if ((($indel_ref_chr eq $ref_chr) and ($indel_query_chr eq $query_chr)) and (($indel_query_start == $query_pos) and ($buff == 1))) {
		$indel_ref_base .= $ref_base;
		$indel_ref_end = $ref_pos;
	    } else {
		# this is the start of a new DELETION
		# record previous DELETION
		if ($type =~ /(INDEL|BOTH)/i) {
		    record_INDEL($output_fh, $indel_type, $indel_ref_chr, $indel_ref_start, $indel_ref_end, $indel_query_chr, $indel_query_start, $indel_query_end, $indel_ref_base, $indel_query_base, $indel_query_orientation);
		}
                    undef $indel_type;
		goto NEW_INDEL;
	    }
	} else {
	    # this is the start of a new DELETION                                                                                                               
	    # record previous DELETION                                                                                                                          
	    if ($type =~ /(INDEL|BOTH)/i) {
		record_INDEL($output_fh, $indel_type, $indel_ref_chr, $indel_ref_start, $indel_ref_end, $indel_query_chr, $indel_query_start, $indel_query_end, $indel_ref_base, $indel_query_base, $indel_query_orientation);
	    }
	    undef $indel_type;
	    goto NEW_INDEL;
	}
    }
}

# record the last INDEL                                                                                                               
if ($type =~ /(INDEL|BOTH)/i) {
    record_INDEL($output_fh, $indel_type, $indel_ref_chr, $indel_ref_start, $indel_ref_end, $indel_query_chr, $indel_query_start, $indel_query_end, $indel_ref_base, $indel_query_base, $indel_query_orientation);
}
undef $indel_type;


print "SNP count: $snp_count\n";
print "INDEL count: $indel_count\n";
		    
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

sub record_SNP {
    my ($output_fh, $ref_chr, $ref_pos, $ref_base, $query_base, $query_chr, $query_pos, $query_orientation) = @_;
    print $output_fh "$ref_chr\t$ref_pos\t.\t$ref_base\t$query_base\t.\tPASS\ttype=SNP;ref_chr=$ref_chr;ref_start=$ref_pos;ref_end=$ref_pos;query_chr=$query_chr;query_start=$query_pos;query_end=$query_pos;query_orientation=$query_orientation\n";
}

sub record_INDEL {
    my ($output_fh, $indel_type, $indel_ref_chr, $indel_ref_start, $indel_ref_end, $indel_query_chr, $indel_query_start, $indel_query_end, $indel_ref_base, $indel_query_base, $indel_query_orientation) = @_;
    my $ref_representation;
    my $query_representation;
    if ($indel_type eq "INSERTION") {
	$ref_representation = substr $refseq{$indel_ref_chr}, $indel_ref_start - 1, 1;
	    $query_representation = $ref_representation . $indel_query_base;
	if ($indel_query_orientation == 1) {
	    $indel_query_start -= 1; 
	} else {
	    $indel_query_end += 1; 
	}
    } else {
	$indel_ref_start -= 1; # tested when query orientation = 1
	$query_representation = substr $refseq{$indel_ref_chr}, $indel_ref_start - 1, 1;
	$ref_representation = $query_representation . $indel_ref_base;
	if ($indel_query_orientation == -1) {
	    $indel_query_start += 1;
	    $indel_query_end += 1;
	}
    }
    print $output_fh "$indel_ref_chr\t$indel_ref_start\t.\t$ref_representation\t$query_representation\t.\tPASS\ttype=$indel_type;ref_chr=$indel_ref_chr;ref_start=$indel_ref_start;ref_end=$indel_ref_end;query_chr=$indel_query_chr;query_start=$indel_query_start;query_end=$indel_query_end;query_orientation=$indel_query_orientation\n";
}
