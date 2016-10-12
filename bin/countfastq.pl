#!/usr/bin/perl -w

use strict;
use lib '/home/redwards/bioinformatics/Modules';
use Rob;
my $rob = new Rob;

my $f=shift || die "fastq file?";
my $count = 0;
my ($longest, $shortest, $total)=(0, 1e12, 0);
my @lengths;
while (my $tple = $rob->stream_fastq($f)) {
	$count += scalar(keys %$tple);
	map {
		my $l = length($tple->{$_}->[0]);
		($l > $longest) ? ($longest=$l):1;
		($l < $shortest) ? ($shortest=$l):1;
		$total += $l;
		push @lengths, $l;
	} keys %$tple;

}
print "$count sequences. Longest: $longest Shortest: $shortest: Total bp: $total Median: ", $rob->median(\@lengths), " Mean: ", (int(($total/$count)*100)/100), "\n";
