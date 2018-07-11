#!/usr/bin/perl
use warnings;
use strict;
use Getopt::Long;

##############################################################
#  script: select_sequences_by_list.pl
#  author: Jia-Xing Yue (GitHub ID: yjx1217)
#  last edited: 2018.05.11
#  description: select (-m normal) or reversely select (-m reverse) fasta or fastq sequences based on the sequence id list
#  example: perl select_sequences_by_list.pl -i input.fa(.gz) -l keep.list -o keep.fa(.gz) -m normal -f fasta 
##############################################################

my ($list, $input, $output, $mode, $format);
$mode = "normal"; # normal mode, output sequences found in the list. 
$format = "fasta";
GetOptions('list|l:s' => \$list,
	   'input|i:s' => \$input,
	   'output|o:s' => \$output,
	   'format|f:s' => \$format,
	   'mode|m:s' => \$mode); # mode = normal or mode = reverse, if mode = reverse, will output sequences not found in the list;

my $list_fh = read_file($list);
my %list = parse_list_file($list_fh);

my $input_fh = read_file($input);
my @input = ();
my %input = ();
if ($format eq "fasta") {
    parse_fasta_file($input_fh, \%input, \@input);
} elsif ($format eq "fastq") {
    parse_fastq_file($input_fh, \%input, \@input);
} else {
    die "Error! Unknown format: $format! Please set the '-f' option as either 'fasta' or 'fastq'. Exit!\n";
}

my $output_fh = write_file($output);

    
foreach my $id (@input) {
    if (exists $list{$id}) {
	if ($mode eq "normal") {
	    if ($format eq "fasta") {
		print $output_fh ">$id\n$input{$id}\n";
	    } else {
		print $output_fh "\@$id\n$input{$id}{'seq'}\n\+\n$input{$id}{'qual'}\n";
	    }
	}
    } else {
	if ($mode eq "reverse") {
	    if ($format eq "fasta") {
		print $output_fh ">$id\n$input{$id}\n";
	    } else {
		print $output_fh "\@$id\n$input{$id}{'seq'}\n\+\n$input{$id}{'qual'}\n";
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

sub parse_list_file {
    my $fh = shift @_;
    my %list = ();
    while (<$fh>) {
	chomp;
	$_ =~ s/\s+$//g;
	$list{$_}++;
    }
    return %list;
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


sub parse_fastq_file {
    my ($fh, $input_hashref, $input_arrayref) = @_;
    my $id = "";
    while (<$fh>) {
	chomp;
	if ($. % 4 == 1) {
	    ($id) = ($_ =~ /^\@(.*)/);
	    push @$input_arrayref, $id;
	} elsif ($. % 4 == 2) {
	    $$input_hashref{$id}{'seq'} = $_;
	} elsif ($. % 4 == 0) {
	    $$input_hashref{$id}{'qual'} = $_;
	}
    }
}

