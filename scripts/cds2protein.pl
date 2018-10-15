#!/usr/bin/perl
use warnings;
use strict;
use Getopt::Long;

##############################################################
#  script: cds2protein.pl
#  author: Jia-Xing Yue (GitHub ID: yjx1217)
#  last edited: 2018.10.03
#  description: translate cds sequences into protein sequences
#  example: perl cds2protein.pl -i input.cds.fa(.gz) -t 1 -p prefix
#  See this link (https://www.ncbi.nlm.nih.gov/Taxonomy/Utils/wprintgc.cgi) for the genetic code table options. Only "1" and "3" are supported so far. 
##############################################################

my ($input, $genetic_code_table, $prefix);
$genetic_code_table = 1;
GetOptions('input|i:s' => \$input, # input file
	   'genetic_code_table|t:s' => \$genetic_code_table, # "1" for the standard code and "3" for the yeast mitochondrial code.
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
    $cds =~ tr/U/T/;
    if (($cds !~ /^ATG/) and ($cds =~ /(TAA|TAG|TGA)$/)) {
	push @pseudogene_test, "unexpected start codon based on standard genentic code;your selected code table is $genetic_code_table";
    } elsif (($cds =~ /^ATG/) and ($cds !~ /(TAA|TAG|TGA)$/)) {
	push @pseudogene_test, "unexpected stop codon based on standard genentic code;your selected code table is $genetic_code_table";
    } elsif (($cds !~ /^ATG/) and ($cds !~ /(TAA|TAG|TGA)$/)) {
	push @pseudogene_test, "unexpected start & end codons based on standard genentic code;your selected code table is $genetic_code_table";
    }
    my $cds_length = length $cds;
    my $length_test = $cds_length % 3;
    my $cds_trim = $cds;

    if ($length_test != 0) {
	push @pseudogene_test, "incorrect CDS length";
	my $cds_phase0 = substr $cds, 0, $cds_length - $length_test;
	my $pep_phase0 = translate($genetic_code_table, $cds_phase0);
	my $pep_phase0_stop_count = () = $pep_phase0 =~ /\*/g;
	my $cds_phase1 = substr $cds, 1, $cds_length - $length_test + 1;
	my $pep_phase1 = translate($genetic_code_table, $cds_phase1);
	my $pep_phase1_stop_count = () = $pep_phase1 =~ /\*/g;
	my $cds_phase2 = substr $cds, 2, $cds_length - $length_test + 2;
	my $pep_phase2 = translate($genetic_code_table, $cds_phase2);
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
    my $protein = translate($genetic_code_table, $cds_trim);
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
    my ($genetic_code_table, $dna) = @_;
    my $protein = "";
    for (my $i = 0; $i < (length($dna) - 2); $i += 3) {
	$protein .= codon2aa($genetic_code_table, substr($dna, $i, 3));
    }
    return $protein;
}


sub codon2aa {
    my ($genetic_code_table, $codon) = @_;
    $codon = uc $codon;
    my %genetic_code = ();
    %{$genetic_code{"1"}} = (
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
    for (my $i = 2; $i <= 31; $i++) {
	%{$genetic_code{$i}} =  %{$genetic_code{"1"}};
	if ($i eq "2") {
	    $genetic_code{$i}{'AGA'} = '*'; # Stop
	    $genetic_code{$i}{'AGG'} = '*'; # Stop
	    $genetic_code{$i}{'AGR'} = '*'; # Stop
	    $genetic_code{$i}{'ATA'} = 'M'; # Methionine
	    $genetic_code{$i}{'TGA'} = 'W'; # Tryptophan
	} elsif ($i eq "3") {
	    $genetic_code{$i}{'ATA'} = 'M'; # Methionine
	    $genetic_code{$i}{'CTT'} = 'T'; # Threonine
	    $genetic_code{$i}{'CTC'} = 'T'; # Threonine
	    $genetic_code{$i}{'CTA'} = 'T'; # Threonine
	    $genetic_code{$i}{'CTG'} = 'T'; # Threonine
	    $genetic_code{$i}{'CTR'} = 'T'; # Threonine
	    $genetic_code{$i}{'CTY'} = 'T'; # Threonine
	    $genetic_code{$i}{'CTS'} = 'T'; # Threonine
	    $genetic_code{$i}{'CTW'} = 'T'; # Threonine
	    $genetic_code{$i}{'CTK'} = 'T'; # Threonine
	    $genetic_code{$i}{'CTM'} = 'T'; # Threonine
	    $genetic_code{$i}{'CTB'} = 'T'; # Threonine
	    $genetic_code{$i}{'CTD'} = 'T'; # Threonine
	    $genetic_code{$i}{'CTH'} = 'T'; # Threonine
	    $genetic_code{$i}{'CTN'} = 'T'; # Threonine
	    $genetic_code{$i}{'TGA'} = 'W'; # Trptophan	    
	} elsif ($i eq "4") {
	    $genetic_code{$i}{'TGA'} = 'W'; # Tryptophan
	} elsif ($i eq "5") {
	    $genetic_code{$i}{'AGA'} = 'S'; # Serine
	    $genetic_code{$i}{'AGG'} = 'S'; # Serine
	    $genetic_code{$i}{'AGR'} = 'S'; # Serine
	    $genetic_code{$i}{'ATA'} = 'M'; # Methionine
	    $genetic_code{$i}{'TGA'} = 'W'; # Tryptophan
	} elsif ($i eq "6") {
	    $genetic_code{$i}{'TAA'} = 'Q'; # Glutamine
	    $genetic_code{$i}{'TAG'} = 'Q'; # Glutamine
	    $genetic_code{$i}{'TAR'} = 'Q'; # Glutamine
	} elsif ($i eq "9") {
	    $genetic_code{$i}{'AAA'} = 'N'; # Asparagine
	    $genetic_code{$i}{'AGA'} = 'S'; # Serine
	    $genetic_code{$i}{'AGG'} = 'S'; # Serine
	    $genetic_code{$i}{'AGR'} = 'S'; # Serine
	    $genetic_code{$i}{'TGA'} = 'W'; # Tryptophan
	} elsif ($i eq "10") {
	    $genetic_code{$i}{'TGA'} = 'C'; # Cysteine
	} elsif ($i eq "12") {
	    $genetic_code{$i}{'CTG'} = 'S'; # Serine
	} elsif ($i eq "13") {
	    $genetic_code{$i}{'AGA'} = 'G'; # Glycine
	    $genetic_code{$i}{'AGG'} = 'G'; # Glycine
	    $genetic_code{$i}{'AGR'} = 'G'; # Glycine
	    $genetic_code{$i}{'ATA'} = 'M'; # Methionine
	    $genetic_code{$i}{'TGA'} = 'W'; # Tryptophan
	} elsif ($i eq "14") {
	    $genetic_code{$i}{'AAA'} = 'N'; # Asparagine
	    $genetic_code{$i}{'AGA'} = 'S'; # Serine
	    $genetic_code{$i}{'AGG'} = 'S'; # Serine
	    $genetic_code{$i}{'AGR'} = 'S'; # Serine
	    $genetic_code{$i}{'TAA'} = 'Y'; # Tyrosine
	    $genetic_code{$i}{'TGA'} = 'W'; # Tryptophan
	} elsif ($i eq "16") {
	    $genetic_code{$i}{'TAG'} = 'L'; # Leucine
	} elsif ($i eq "21") {
	    $genetic_code{$i}{'TGA'} = 'W'; # Tryptophan
	    $genetic_code{$i}{'ATA'} = 'M'; # Methionine
	    $genetic_code{$i}{'AGA'} = 'S'; # Serine
	    $genetic_code{$i}{'AGG'} = 'S'; # Serine
	    $genetic_code{$i}{'AGR'} = 'S'; # Serine
	    $genetic_code{$i}{'AAA'} = 'N'; # Asparagine
	} elsif ($i eq "22") {
	    $genetic_code{$i}{'TAG'} = 'L'; # Leucine
	} elsif ($i eq "24") {
	    $genetic_code{$i}{'AGA'} = 'S'; # Serine
	    $genetic_code{$i}{'AGG'} = 'K'; # Lysine
	    $genetic_code{$i}{'TGA'} = 'W'; # Tryptophan
	} elsif ($i eq "25") {
	    $genetic_code{$i}{'TGA'} = 'G'; # Glycine
	} elsif ($i eq "26") {
	    $genetic_code{$i}{'ATG'} = 'A'; # Alanine
	} elsif ($i eq "27") {
	    $genetic_code{$i}{'TAG'} = 'Q'; # Glutamine
	    $genetic_code{$i}{'TAA'} = 'Q'; # Glutamine
	    $genetic_code{$i}{'TAR'} = 'Q'; # Glutamine
	} elsif ($i eq "29") {
	    $genetic_code{$i}{'TAG'} = 'Y'; # Tyrosine
            $genetic_code{$i}{'TAA'} = 'Y'; # Tyrosine
            $genetic_code{$i}{'TAR'} = 'Y'; # Tyrosine
	} elsif ($i eq "30") {
	    $genetic_code{$i}{'TAG'} = 'E'; # Glutamic Acid
	    $genetic_code{$i}{'TAA'} = 'E'; # Glutamic Acid
	    $genetic_code{$i}{'TAR'} = 'E'; # Glutamic Acid
	} elsif ($i eq "31") {
	    $genetic_code{$i}{'TGA'} = 'W'; # Tryptophan
	}
    }
    
    if (exists $genetic_code{$genetic_code_table}{$codon}) {
        return $genetic_code{$genetic_code_table}{$codon};
    } else {
	#print STDERR "Bad codon \"$codon\"!!\n";
	return 'X';
    }
}

