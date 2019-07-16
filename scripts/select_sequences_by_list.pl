#!/usr/bin/perl
use warnings;
use strict;
use Getopt::Long;

##############################################################
#  script: select_sequences_by_list.pl
#  author: Jia-Xing Yue (GitHub ID: yjx1217)
#  last edited: 2019.08.22
#  description: select (-m normal) or reversely select (-m reverse) fasta or fastq sequences based on the sequence id list
#  example: perl select_sequences_by_list.pl -i input.fa(.gz) -l keep.list -o keep.fa(.gz) -m normal -f fasta -r by_sequence
##############################################################

my ($list, $input, $output, $selection_mode, $format, $ranking_order);
$selection_mode = "normal"; # normal mode, output sequences found in the list. 
$ranking_order = "by_sequence";
$format = "fasta";
GetOptions('list|l:s' => \$list,
	   'input|i:s' => \$input,
	   'output|o:s' => \$output,
	   'format|f:s' => \$format,
	   'selection_mode|m:s' => \$selection_mode,
	   'ranking_order|r:s' => \$ranking_order); # mode = normal or mode = reverse, if mode = reverse, will output sequences not found in the list;

print "Current option settings:\n";
print "input: $input\n";
print "output: $output\n";
print "list: $list\n";
print "format: $format\n";
print "selection_mode: $selection_mode\n";
print "ranking_order: $ranking_order\n";

if ($format !~ /(fasta|fastq)/) {
    die "Unrecognized format: $format. Please reset it to either \"fasta\" (default) or \"fastq\" via the option -f";
}
if ($selection_mode !~ /(normal|reverse)/) {
    die "Unrecognized selection_mode: $selection_mode. Please reset it to either \"normal\" (default) or \"reverse\" via the option -m";
}
if ($ranking_order !~ /(by_sequence|by_list)/) {
    die "Unrecognized ranking_order: $ranking_order. Please reset it to either \"by_sequence\" (default) or \"by_list\" via the option -r";
}

my $list_fh = read_file($list);
my %list = ();
my @list = ();
parse_list_file($list_fh, \%list, \@list);

my $input_fh = read_file($input);
my %input = ();
my @input = ();

if ($format eq "fasta") {
    parse_fasta_file($input_fh, \%input, \@input);
} else {
    parse_fastq_file($input_fh, \%input, \@input);
}

my $output_fh = write_file($output);

if ($ranking_order eq "by_sequence") {    
    foreach my $id (@input) {
	if (exists $list{$id}) {
	    if ($selection_mode eq "normal") {
		if ($format eq "fasta") {
		    print $output_fh ">$id\n$input{$id}\n";
		} else {
		    print $output_fh "\@$id\n$input{$id}{'seq'}\n\+\n$input{$id}{'qual'}\n";
		}
	    }
	} else {
	    if ($selection_mode eq "reverse") {
		if ($format eq "fasta") {
		    print $output_fh ">$id\n$input{$id}\n";
		} else {
		    print $output_fh "\@$id\n$input{$id}{'seq'}\n\+\n$input{$id}{'qual'}\n";
		}
	    }
	}
    }
} else {
    foreach my $id (@list) {
	if (exists $input{$id}) {
	    if ($selection_mode eq "normal") {
		if ($format eq "fasta") {
		    print $output_fh ">$id\n$input{$id}\n";
		} else {
		    print $output_fh "\@$id\n$input{$id}{'seq'}\n\+\n$input{$id}{'qual'}\n";
		}
	    }
	} else {
	    if ($selection_mode eq "reverse") {
		if ($format eq "fasta") {
		    print $output_fh ">$id\n$input{$id}\n";
		} else {
		    print $output_fh "\@$id\n$input{$id}{'seq'}\n\+\n$input{$id}{'qual'}\n";
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

sub parse_list_file {
    my ($fh, $input_hashref, $input_arrayref) = @_;
    while (<$fh>) {
        chomp;
        $_ =~ s/\s+$//g;
        if (not exists $$input_hashref{$_}) {
            $$input_hashref{$_} = 1;
            push @$input_arrayref, $_;
        } else {
            $$input_hashref{$_}++;
        }
    }
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

