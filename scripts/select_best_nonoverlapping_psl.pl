#!/usr/bin/env perl

use strict;
use warnings;
use Getopt::Long;

my ($input, $output);
GetOptions(
    'input|i=s'  => \$input,
    'output|o=s' => \$output,
) or die "Usage: $0 -i input.psl -o output.psl\n";

die "Usage: $0 -i input.psl -o output.psl\n"
    unless defined $input && defined $output;

open my $input_fh, '<', $input or die "Cannot open $input: $!\n";

my @header;
my @hits;
my $input_index = 0;

while (my $line = <$input_fh>) {
    if ($line !~ /^\d+\t/) {
        push @header, $line;
        next;
    }

    chomp $line;
    my @field = split /\t/, $line, -1;
    die "Invalid PSL record with fewer than 21 columns in $input:\n$line\n"
        if @field < 21;

    my ($matches, $mismatches, $repeat_matches) = @field[0, 1, 2];
    my ($query_inserts, $target_inserts) = @field[4, 6];
    my ($query_name, $query_size, $query_start, $query_end) = @field[9 .. 12];
    my ($target_name, $target_start, $target_end) = @field[13, 15, 16];

    for my $value (
        $matches, $mismatches, $repeat_matches, $query_inserts,
        $target_inserts, $query_size, $query_start, $query_end,
        $target_start, $target_end,
    ) {
        die "Invalid numeric field in PSL record from $input:\n$line\n"
            unless defined $value && $value =~ /^\d+$/;
    }
    die "Invalid target interval in PSL record from $input:\n$line\n"
        unless $target_start < $target_end;
    die "Invalid query interval in PSL record from $input:\n$line\n"
        unless $query_size > 0 && $query_start < $query_end;

    my $identity_denominator = $matches + $repeat_matches + $mismatches;
    my $identity = $identity_denominator > 0
        ? ($matches + $repeat_matches) / $identity_denominator
        : 0;
    my $aligned_bases = $query_end - $query_start;

    push @hits, {
        line          => "$line\n",
        input_index   => $input_index++,
        target_name   => $target_name,
        target_start  => $target_start + 0,
        target_end    => $target_end + 0,
        query_name    => $query_name,
        score         => $matches + int($repeat_matches / 2)
                         - $mismatches - $query_inserts - $target_inserts,
        identity      => $identity,
        query_coverage => $aligned_bases / $query_size,
        aligned_bases => $aligned_bases,
    };
}
close $input_fh;

# Greedy non-maximum suppression across all query names. PSL target intervals are
# zero-based, half-open, so adjacent intervals (end == start) do not overlap.
my @ranked_hits = sort {
       $b->{score}          <=> $a->{score}
    || $b->{identity}       <=> $a->{identity}
    || $b->{query_coverage} <=> $a->{query_coverage}
    || $b->{aligned_bases}  <=> $a->{aligned_bases}
    || $a->{target_name}    cmp $b->{target_name}
    || $a->{target_start}   <=> $b->{target_start}
    || $a->{target_end}     <=> $b->{target_end}
    || $a->{query_name}     cmp $b->{query_name}
    || $a->{input_index}    <=> $b->{input_index}
} @hits;

my %selected_by_target;
my @selected;
HIT:
for my $hit (@ranked_hits) {
    for my $kept (@{ $selected_by_target{$hit->{target_name}} // [] }) {
        if ($hit->{target_start} < $kept->{target_end}
            && $kept->{target_start} < $hit->{target_end}) {
            next HIT;
        }
    }
    push @{ $selected_by_target{$hit->{target_name}} }, $hit;
    push @selected, $hit;
}

@selected = sort { $a->{input_index} <=> $b->{input_index} } @selected;

open my $output_fh, '>', $output or die "Cannot write $output: $!\n";
print {$output_fh} @header;
print {$output_fh} $_->{line} for @selected;
close $output_fh or die "Cannot close $output: $!\n";

my $removed = @hits - @selected;
print STDERR "Selected " . scalar(@selected) . " non-overlapping PSL hits; "
    . "removed $removed overlapping lower-ranked hits.\n";
