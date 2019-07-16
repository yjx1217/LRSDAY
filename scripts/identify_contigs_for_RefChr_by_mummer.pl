#!/usr/bin/perl
use warnings;
use strict;
use Getopt::Long;

##############################################################
#  script: identify_contigs_for_RefChr_by_mummer.pl
#  author: Jia-Xing Yue (GitHub ID: yjx1217)
#  last edited: 2019.08.14
#  description: identify contigs corresponding to the specified reference chromosome based on mummer show-coords output
#  example: perl identify_contigs_for_RefChr_by_mummer.pl -i mummer.coords -query_chr_list query_chr.list -assembly_fasta assembly.fasta -cov 90 -o chrMT.match.list
##############################################################


my ($input, $output, $query_chr_list, $assembly_fasta, $cov);
$cov = 75;
GetOptions('input|i:s' => \$input, # input blast tabular output
	   'query_chr_list|query_chr_list:s' => \$query_chr_list, # a simple list file containing query chr ids from the reference genome
	   'assembly_fasta|assembly_fasta:s' => \$assembly_fasta, # a simple list file containing query chr sequences from the reference genome
	   'coverage|cov:s' => \$cov, # cummulative query contig coverage
	   'output|o:s' => \$output); # filtered blast tabular output


my $input_fh = read_file($input);
my %match = ();
my $output_fh = write_file($output);

my $query_chr_list_fh = read_file($query_chr_list);
my %query_chr_list = ();
my @query_chr_list = ();
parse_list_file($query_chr_list_fh, \%query_chr_list, \@query_chr_list);

my $assembly_fasta_fh = read_file($assembly_fasta);
my %assembly_fasta = ();
my @assembly_fasta = ();
parse_fasta_file($assembly_fasta_fh, \%assembly_fasta, \@assembly_fasta);

while (<$input_fh>) {
    chomp;
    (1..4) and next;
    /^\#/ and next;
    /^\s*$/ and next;
    my ($ref_start, $ref_end, $query_start, $query_end, $ref_match_length, $query_match_length, $ref_length, $query_length, $ref_cov, $query_cov, $ref_id, $query_id) = split /\t/, $_;
    foreach my $query_chr (@query_chr_list) {
	if ($ref_id =~ /$query_chr/) {
	    if (exists $match{$query_id}) {
		$match{$query_id} += $query_cov;
	    } else {
		$match{$query_id} = $query_cov;
	    }
	}
    }
}

foreach my $query_id (sort {$match{$b} <=> $match{$a}} keys %match) {
    my $query_seq = $assembly_fasta{$query_id};
    my $query_seq_length = length $query_seq;
    my $query_seq_N_length = () = $query_seq =~ /N/gi;
    my $raw_query_cov = $match{$query_id};
    my $adjusted_query_cov = ($query_seq_length * $raw_query_cov)/($query_seq_length - $query_seq_N_length);
    print "query_id=$query_id, query_seq_length=$query_seq_length, query_seq_N_length=$query_seq_N_length, raw_query_cov=$raw_query_cov, adjusted_query_cov=$adjusted_query_cov\n";
    if ($adjusted_query_cov >= $cov) {
	print $output_fh "$query_id\n";
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

sub parse_list_file {
    my ($list_fh, $list_hashref, $list_arrayref) = @_;
    while (<$list_fh>) {
        chomp;
	/^#/ and next;
	/^\s*$/ and next;
        $_ =~ s/\s+$//g;
        if (not exists $$list_hashref{$_}) {
            $$list_hashref{$_} = 1;
            push @$list_arrayref, $_;
        } else {
            $$list_hashref{$_}++;
        }
    }
}


sub parse_fasta_file {
    my ($fasta_fh, $fasta_hashref, $fasta_arrayref) = @_;
    my $fasta_id = "";
    while (<$fasta_fh>) {
	chomp;
        if (/^\s*$/) {
            next;
        } elsif (/^\s*#/) {
		 next;
	    } elsif (/^>(\S+)/) {
		$fasta_id = $1;
		push @$fasta_arrayref, $fasta_id;
		$$fasta_hashref{$fasta_id} = "";
	    } else {
		$$fasta_hashref{$fasta_id} .= $_;
	    }
    }
}

