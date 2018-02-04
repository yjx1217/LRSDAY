#!/usr/bin/perl
use warnings;
use strict;
use Getopt::Long;

##############################################################
#  script: cds2protein.pl
#  author: Jia-Xing Yue (GitHub ID: yjx1217)
#  last edited: 2018.01.24
#  description: translate cds sequences into protein sequences
#  example: perl cds2protein.pl -i input.cds.fa(.gz) -p prefix
##############################################################

my ($input, $prefix);
GetOptions('input|i:s' => \$input, # input file
	   'prefix|p:s' => \$prefix); # file name prefix for output files

my $input_fh = read_file($input);
my @input = ();
my %input = ();
parse_fasta_file($input_fh, \%input, \@input);
my $trimmed_cds_output = "$prefix.trimmed_cds.fa";
my $trimmed_cds_output_fh = write_file($trimmed_cds_output);
my $pep_output = "$prefix.pep.fa";
my $pep_output_fh = write_file($pep_output);
my $trimmed_cds_log = "$prefix.trimmed_cds.log";
my $trimmed_cds_log_fh = write_file($trimmed_cds_log);

my $pseudogene_output = "$prefix.manual_check.list";
my $pseudogene_output_fh = write_file($pseudogene_output);

foreach my $id (@input) {
    # print "id = $id\n";
    my @pseudogene_test = ();
    my $cds = uc $input{$id};
    if (($cds !~ /^ATG/) or ($cds !~ /(TAA|TAG|TGA)$/)) {
	push @pseudogene_test, "unexpected start/end codons"
    }
    my $cds_length = length $cds;
    my $length_test = $cds_length % 3;
    my $cds_trim = $cds;

    if ($length_test != 0) {
	push @pseudogene_test, "incorrect CDS length";
	my $cds_phase0 = substr $cds, 0, $cds_length - $length_test;
	my $pep_phase0 = translate($cds_phase0);
	my $pep_phase0_stop_count = () = $pep_phase0 =~ /\*/g;
	my $cds_phase1 = substr $cds, 1, $cds_length - $length_test + 1;
	my $pep_phase1 = translate($cds_phase1);
	my $pep_phase1_stop_count = () = $pep_phase1 =~ /\*/g;
	my $cds_phase2 = substr $cds, 2, $cds_length - $length_test + 2;
	my $pep_phase2 = translate($cds_phase2);
	my $pep_phase2_stop_count = () = $pep_phase2 =~ /\*/g;
	if (($pep_phase0_stop_count <= $pep_phase1_stop_count) and ($pep_phase0_stop_count <= $pep_phase2_stop_count)) {
	    # print "my guess is phase 0, resulting a total of $pep_phase0_stop_count stop codons\n";
	    $cds_trim = $cds_phase0;
	    print $trimmed_cds_log_fh "incorrect_length: $id, phase_guess: 0\n";
	} elsif (($pep_phase1_stop_count <= $pep_phase0_stop_count) and ($pep_phase1_stop_count <= $pep_phase2_stop_count)) {
            # print "my guess is phase 1, resulting a total of $pep_phase1_stop_count stop codons\n";
            $cds_trim = $cds_phase1;
	    print $trimmed_cds_log_fh "incorrect_length: $id, phase_guess: 1\n";
	} elsif (($pep_phase2_stop_count <= $pep_phase0_stop_count) and ($pep_phase2_stop_count <= $pep_phase1_stop_count)) {
            # print "my guess is phase 2, resulting a total of $pep_phase2_stop_count stop codons\n";
            $cds_trim = $cds_phase2;
	    print $trimmed_cds_log_fh "incorrect_length: $id, phase_guess: 2\n";
	}
    }
    print $trimmed_cds_output_fh ">$id\n$cds_trim\n";
    my $protein = translate($cds_trim);
    # get rid of the ending stop codon for translated protein sequences
    $protein =~ s/\*$//;
    print $pep_output_fh ">$id\n$protein\n";
    if (($protein =~ /\*/) and ($length_test == 0)) {
	push @pseudogene_test, "internal stop codon(s)";
    }
    if ((scalar @pseudogene_test) > 0) {
	my $pseudogene_test = join ";", @pseudogene_test;
	print $pseudogene_output_fh "$id\t$pseudogene_test\n";
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


sub translate {
    my $dna = shift @_;
    my $protein = "";
    for (my $i = 0; $i < (length($dna) - 2); $i += 3) {
	$protein .= codon2aa(substr($dna, $i, 3));
    }
    return $protein;
}


sub codon2aa {
    my $codon  = shift @_;
    $codon = uc $codon;
    my %genetic_code = (
	'TCA' => 'S',    # Serine
	'TCC' => 'S',    # Serine
	'TCG' => 'S',    # Serine
	'TCT' => 'S',    # Serine
	'TCN' => 'S',    # Serine
	'TTC' => 'F',    # Phenylalanine
	'TTT' => 'F',    # Phenylalanine
	'TTA' => 'L',    # Leucine
	'TTG' => 'L',    # Leucine
	'TAC' => 'Y',    # Tyrosine
	'TAT' => 'Y',    # Tyrosine
	'TAA' => '*',    # Stop
	'TAG' => '*',    # Stop
	'TGC' => 'C',    # Cysteine
	'TGT' => 'C',    # Cysteine
	'TGA' => '*',    # Stop
	'TGG' => 'W',    # Tryptophan
	'CTA' => 'L',    # Leucine
	'CTC' => 'L',    # Leucine
	'CTG' => 'L',    # Leucine
	'CTT' => 'L',    # Leucine
	'CTN' => 'L',    # Leucine
	'CCA' => 'P',    # Proline
	'CCC' => 'P',    # Proline
	'CCG' => 'P',    # Proline
	'CCT' => 'P',    # Proline
	'CCN' => 'P',    # Proline
	'CAC' => 'H',    # Histidine
	'CAT' => 'H',    # Histidine
	'CAA' => 'Q',    # Glutamine
	'CAG' => 'Q',    # Glutamine
	'CGA' => 'R',    # Arginine
	'CGC' => 'R',    # Arginine
	'CGG' => 'R',    # Arginine
	'CGT' => 'R',    # Arginine
	'CGN' => 'R',    # Arginine
	'ATA' => 'I',    # Isoleucine
	'ATC' => 'I',    # Isoleucine
	'ATT' => 'I',    # Isoleucine
	'ATG' => 'M',    # Methionine
	'ACA' => 'T',    # Threonine
	'ACC' => 'T',    # Threonine
	'ACG' => 'T',    # Threonine
	'ACT' => 'T',    # Threonine
	'ACN' => 'T',    # Threonine
	'AAC' => 'N',    # Asparagine
	'AAT' => 'N',    # Asparagine
	'AAA' => 'K',    # Lysine
	'AAG' => 'K',    # Lysine
	'AGC' => 'S',    # Serine
	'AGT' => 'S',    # Serine
	'AGA' => 'R',    # Arginine
	'AGG' => 'R',    # Arginine
	'GTA' => 'V',    # Valine
	'GTC' => 'V',    # Valine
	'GTG' => 'V',    # Valine
	'GTT' => 'V',    # Valine
	'GTN' => 'V',    # Valine
	'GCA' => 'A',    # Alanine
	'GCC' => 'A',    # Alanine
	'GCG' => 'A',    # Alanine
	'GCT' => 'A',    # Alanine
	'GCN' => 'A',    # Alanine
	'GAC' => 'D',    # Aspartic Acid
	'GAT' => 'D',    # Aspartic Acid
	'GAA' => 'E',    # Glutamic Acid
	'GAG' => 'E',    # Glutamic Acid
	'GGA' => 'G',    # Glycine
	'GGC' => 'G',    # Glycine
	'GGG' => 'G',    # Glycine
	'GGT' => 'G',    # Glycine
	'GGN' => 'G',    # Glycine
	'MGR' => 'R',    # Arginine
	'AAY' => 'N',    # Asparagine
	'GAY' => 'D',    # Aspartic Acid
	'TGY' => 'C',    # Cysteine
	'CAR' => 'Q',    # Glutamine
	'GAR' => 'E',    # Glutamic Acid
	'CAY' => 'H',    # Histidine
	'ATH' => 'I',    # Isoleucine
	'YTR' => 'L',    # Leucine
	'AAR' => 'K',    # Lysine
	'TTY' => 'F',    # Phenylalanine
	'TCN' => 'S',    # Serine
	'AGY' => 'S',    # Serine
	'ACN' => 'T',    # Threonine
	'TAY' => 'Y',    # Tyrosine
	'TAR' => '*',    # Stop
	'TRA' => '*',    # Stop
	'---' => '-'     # alignment gap
	);
    
    if (exists $genetic_code{$codon}) {
        return $genetic_code{$codon};
    } else {
	#print STDERR "Bad codon \"$codon\"!!\n";
	return 'X';
    }
}

