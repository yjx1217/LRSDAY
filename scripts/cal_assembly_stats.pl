#! /usr/bin/perl
use strict;
use warnings;
use Getopt::Long;
use List::Util qw(sum min max);

##############################################################
#  script: cal_assembly_stats.pl
#  author: Jia-Xing Yue (GitHub ID: yjx1217)
#  last edited: 2022.07.22
#  description: calculate basice statistics for genome assemblies, such as GC%, N50, N90, etc ...
#  example: perl cal_assembly_stats.pl -i input.fa(.gz) -o output.txt
##############################################################

my ($input, $output);

GetOptions('input|i:s' => \$input, # input genome fasta file
	   'output|o:s' => \$output); # output summary statistics file

my $input_fh = read_file($input);
my @input = ();
my %input = ();
parse_fasta_file($input_fh, \%input, \@input);

my $output_fh = write_file($output);

my %base_count = ();
my @base = qw(A T G C N);
my %length = ();

foreach my $id (@input) {
    my $seq = uc $input{$id};
    $length{$id} = length $seq;
    foreach my $base (@base) {
	my $bc  = () = $seq =~ /$base/g;
	# print "id=$id, base=$base, base_count=$bc\n";
	if (exists $base_count{$base}) {
	    $base_count{$base} += $bc
	} else {
	    $base_count{$base} = $bc;
	}
    }
}

my @length = values %length;
my $total_sequence_count = scalar @length;
my $total_sequence_length = sum(@length);
my $min_length = min(@length);
my $max_length = max(@length);
my $mean_length = cal_mean(\@length);
$mean_length = sprintf("%.2f", $mean_length); 
my $median_length = cal_median(\@length);
$median_length = sprintf("%.2f", $median_length); 
my ($N50, $L50) = cal_N50(\@length, 50);
my ($N90, $L90) = cal_N50(\@length, 90);

my %base_pct = ();
foreach my $base (@base) {
    my $base_pct = $base_count{$base} * 100 / $total_sequence_length;
    $base_pct{$base} =  sprintf("%.2f", $base_pct);
}

$base_count{'AT'} = $base_count{'A'} + $base_count{'T'};
$base_count{'GC'} = $base_count{'G'} + $base_count{'C'};
$base_pct{'AT'} =  ($base_count{'AT'} * 100) / $total_sequence_length;
$base_pct{'GC'} =  ($base_count{'GC'} * 100) / $total_sequence_length;
$base_pct{'AT'} =  sprintf("%.2f", $base_pct{'AT'});
$base_pct{'GC'} =  sprintf("%.2f", $base_pct{'GC'});
$base_pct{'AT_adjusted'} = ($base_count{'AT'} * 100)/($base_count{'AT'} + $base_count{'GC'});
$base_pct{'GC_adjusted'} = ($base_count{'GC'} * 100)/($base_count{'AT'} + $base_count{'GC'});
$base_pct{'AT_adjusted'} =  sprintf("%.2f", $base_pct{'AT_adjusted'});
$base_pct{'GC_adjusted'} =  sprintf("%.2f", $base_pct{'GC_adjusted'});

print $output_fh "total sequence count: $total_sequence_count\n";
print $output_fh "total sequence length: $total_sequence_length\n";
print $output_fh "min sequence length: $min_length\n";
print $output_fh "max sequence length: $max_length\n";
print $output_fh "mean sequence length: $mean_length\n";
print $output_fh "median sequence length: $median_length\n";
print $output_fh "###################\n";
print $output_fh "N50: $N50\n";
print $output_fh "L50: $L50\n";
print $output_fh "N90: $N90\n";
print $output_fh "L90: $L90\n";
print $output_fh "###################\n";
print $output_fh "A count: $base_count{'A'}\n";
print $output_fh "T count: $base_count{'T'}\n";
print $output_fh "G count: $base_count{'G'}\n";
print $output_fh "C count: $base_count{'C'}\n";
print $output_fh "(A+T) count: $base_count{'AT'}\n";
print $output_fh "(G+C) count: $base_count{'GC'}\n";
print $output_fh "N count: $base_count{'N'}\n";
print $output_fh "###################\n";
print $output_fh "A%: $base_pct{'A'}\n";
print $output_fh "T%: $base_pct{'T'}\n";
print $output_fh "G%: $base_pct{'G'}\n";
print $output_fh "C%: $base_pct{'C'}\n";
print $output_fh "(A+T)%: $base_pct{'AT'}\n";
print $output_fh "(G+C)%: $base_pct{'GC'}\n";
print $output_fh "N%: $base_pct{'N'}\n";
print $output_fh "((A+T)/(A+T+G+C))%: $base_pct{'AT_adjusted'}\n";
print $output_fh "((G+C)/(A+T+G+C))%: $base_pct{'GC_adjusted'}\n";


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
    my ($fh, $seq_hashref, $seq_arraryref) = @_;
    my $seq_name = "";
    while (<$fh>) {
        chomp;
        if (/^\s*$/) {
            next;
        } elsif (/^\s*\#/) {
            next;
        } elsif (/^>(.*)/) {
            $seq_name = $1;
            $$seq_hashref{$seq_name} = "";
            push @$seq_arraryref, $seq_name;
        } else {
            $$seq_hashref{$seq_name} .= $_;
        }
    }
}

sub cal_mean {
	my $data_arrayref = shift @_;
	my $data_num = scalar @$data_arrayref;
	my $data_sum = sum @$data_arrayref;
	my $mean;
	if ($data_num == 0) {
	    $mean = "NA";
	} else {
	    $mean = $data_sum/$data_num;
	}
	return $mean;
}


sub cal_median {
	my $data_arrayref = shift @_;
	my @data_sorted = sort{$a<=>$b} @$data_arrayref;
	my $median;
	my $data_num = scalar @data_sorted;
	if ($data_num == 0) {
	    $median = "NA";
	} elsif ($data_num % 2 == 0) {
	    $median = ($data_sorted[$data_num/2 - 1] + $data_sorted[$data_num/2])/2;
	} else {
	    $median = $data_sorted[$data_num/2];
	}
	return $median;
}

sub cal_stdev {
    my $data_arrayref = shift @_;
    my $data_num = scalar @$data_arrayref;
    my $stdev;
    if ($data_num == 1){
        $stdev = 0;
    } else {
        my $mean = cal_mean($data_arrayref);
        my $sq_sum = 0;
        foreach my $i (@$data_arrayref) {
            $sq_sum += ($mean -$i) ** 2;
        }
        $stdev = ($sq_sum / ($data_num - 1)) ** 0.5;
    }
    return $stdev;
}




sub cal_N50 {
	my ($data_arrayref, $k) = @_;
	my @data_sorted = sort{$b<=>$a} @$data_arrayref;
        my $data_num = scalar @data_sorted;
	my $data_sum = sum @data_sorted;
	my $cumulative_sum = 0;
	my $N50;
	my $L50;
	for (my $i = 0; $i < $data_num; $i++){
	    $cumulative_sum += $data_sorted[$i];
	    if ($cumulative_sum >= ($data_sum * $k / 100)){
		$N50 = $data_sorted[$i];
		$L50 = $i + 1;
		last;
	    }
	}
	return ($N50, $L50);
}

