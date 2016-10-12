#!/usr/bin/perl -w

=pod 

=head1 split_by_lines

Split a file into smaller files based on how many lines you specifiy. This will just break at line breaks.

=cut

use strict;

my ($file, $lines)=@ARGV;
die "$0 <file> <lines>" unless ($file && $lines);

my $fc=1;
my $out=$file;
$out =~ s/.txt$//;

open(IN, $file) || die "Can't open $file";
while (-e "$out.split.$fc.txt") {$fc++}
open(OUT, ">$out.split.$fc.txt") || die "Can't open out.split.$fc.txt";
my $count=0;
while (<IN>)
{
	$count++;
	unless ($count % $lines)
	{
		close OUT;
		while (-e "$out.split.$fc.txt") {$fc++}
		open(OUT, ">$out.split.$fc.txt") || die "Can't open out.split.$fc.txt";
	}
	print OUT;
}

