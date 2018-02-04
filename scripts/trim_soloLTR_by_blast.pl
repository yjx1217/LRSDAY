#!/usr/bin/perl
use warnings;
use strict;
use Getopt::Long;

##############################################################
#  script: trim_soloLTR_by_blast.pl
#  author: Jia-Xing Yue (GitHub ID: yjx1217)
#  last edited: 2017.06.17
#  description: trim soloLTR based on BLAST search result 
#  example: perl trim_soloLTR_by_blast.pl -q raw.soloLTR.fa  -i 70 -l 100 -p refined
##############################################################

my ($query, $blast, $pct_identity_cutoff, $prefix, $aln_length_cutoff);
$pct_identity_cutoff = 70;
$aln_length_cutoff = 100;

GetOptions('query|q:s' => \$query, # query sequences in fasta format used for BLAST search
	   'blast|b:s' => \$blast, # BLAST result in tabular format
	   'pct_identity_cutoff|i:s' => \$pct_identity_cutoff, # percent of identity cutoff, e.g. 70
	   'aln_length|l:s' => \$aln_length_cutoff, # alignment length cutoff, e.g. 100
	   'prefix|p:s' => \$prefix); # file name prefix for the output files

my $query_fh = read_file($query);
my @query = ();
my %query = ();
parse_fasta_file($query_fh, \%query, \@query);

my $blast_fh = read_file($blast);
my %blast = parse_blast_file($blast_fh);

my $output = "$prefix.TY_REannotate.soloLTR.refined.gff";
my $output_fh = write_file($output);
foreach my $query_id (@query) {
    my ($id, $type, $chr, $region, $strand) = split /:/, $query_id;
    my ($start, $end) = split /-/, $region;
    if (exists $blast{$query_id}) {
	my @blast = @{$blast{$query_id}};
	foreach $blast (@blast) {
	    my @hit = split /\t/, $blast;
	    my $subject = $hit[1];
	    ($type) = ($subject =~ /^(\S+?)\-/);
	    $type = $type . "_soloLTR";
	    my $q_start = $hit[6];
	    my $q_end = $hit[7];
	    my $s_start = $hit[8];
	    my $s_end = $hit[9];
	    my $match_direction;
	    if ($s_start <= $s_end) {
		$match_direction = "+";
	    } else {
		$match_direction = "-";
	    }
	    if (($strand eq "+") and ($match_direction eq "+")) {
		$strand = "+";
	    } elsif (($strand eq "-") and ($match_direction eq "+")) {
		$strand = "-";
	    } elsif (($strand eq "+") and ($match_direction eq "-")) {
		$strand = "-";
	    } elsif (($strand eq "-") and ($match_direction eq "-")) {
		$strand = "+";
	    }
	    my $new_start = $start + $q_start - 1;
	    my $new_end = $start + $q_end - 1;
	    my $new_id = "$type:$chr:$new_start-$new_end:$strand";
	    print $output_fh "$chr\t$prefix\t$type\t$new_start\t$new_end\t.\t$strand\t.\tID=$new_id;Name=$new_id\n";
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
    if ($file =~ /\.gz$/) {
        open($fh, "| gzip -c >$file") or die "can't open $file\n";
    } else {
        open($fh, ">$file") or die "can't open $file\n";
    }
    return $fh;
}  



sub parse_fasta_file {
    my ($fh, $seq_hashref, $seq_arrayref) = @_;
    my $seq_name = "";
    while (<$fh>) {
        chomp;
        if (/^\s*$/) {
            next;
        } elsif(/^\#/) {
            next;
        } elsif(/^>(\S+)/) {
            $seq_name = $1;
            push @$seq_arrayref, $seq_name;
            $$seq_hashref{$seq_name} = "";
        } else {
            $$seq_hashref{$seq_name} .= $_;
        }
    }
}



sub parse_blast_file {
    my $fh = shift @_;
    my %blast = ();
    while (<$fh>) {
	chomp;
	/^#/ and next;
	/^\s*$/ and next;
	my @line = split /\t/, $_;
	my $query_id = $line[0];
	my $pct_identity = $line[2];
	my $aln_length = $line[3];
	my $e_value = $line[10];
	my $bit_score = $line[11];
	$pct_identity =~ s/\s//gi;
	if (($pct_identity >= $pct_identity_cutoff) and ($aln_length >= $aln_length_cutoff)){
	    if (not exists $blast{$query_id}) {
		@{$blast{$query_id}} = ("$_");
	    } else {
		push @{$blast{$query_id}}, $_;
	    }
	}
    }
    return %blast;
}
