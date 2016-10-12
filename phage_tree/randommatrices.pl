#!/usr/bin/perl -w

#    Copyright 2001, 20002 Rob Edwards
#    For updates, more information, or to discuss the scripts
#    please contact Rob Edwards at redwards@utmem.edu or via http://www.salmonella.org/
#
#    This file is part of The Phage Proteome Scripts developed by Rob Edwards.
#
#    Tnese scripts are free software; you can redistribute and/or modify
#    them under the terms of the GNU General Public License as published by
#    the Free Software Foundation; either version 2 of the License, or
#    (at your option) any later version.
#
#    They are distributed in the hope that they will be useful,
#    but WITHOUT ANY WARRANTY; without even the implied warranty of
#    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#    GNU General Public License for more details.
#
#    You should have received a copy of the GNU General Public License
#    in the file (COPYING) along with these scripts; if not, write to the Free Software
#    Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA
#
#



# randomize the matrix file

use strict;

my $file = shift || die "randommatrix.pl <file> <Number of times to randomize> <-s for subreplicates>\n";
my $randtimes = shift || die "randommatrix.pl <file> <Number of times to randomize> <-s for subreplicates>\n";
my $args;
if (@ARGV) {$args = join (" ", @ARGV)}
else {print STDERR "No subreplicates in file\n"}


my @matrix; my %genomenumber;
open (IN, $file) || die "Can't open $file\n";
while (<IN>) {
	next unless (/^genome/);
	chomp;
	s/genome//;
	my @line = split /\s+/;
	my $genome = shift @line;
	if ($args =~ /-s/) {
		my @templine;
		for (my $i=0; $i <= $#line; $i +=2) {
			 my $linetemp = join (" ", $line[$i], $line[$i+1]);
			 push (@templine, $linetemp);
			 }
		@line = @templine;
	}
	foreach my $i (0 .. $#line) {$matrix[$genome][$i+1] = $line[$i]}
	$genomenumber{$genome}=1;
	}
close IN;


open (OUT, ">random.matrices") || die "can't open random.matrices\n";

my @genomes = sort {$a <=> $b} keys %genomenumber;
foreach my $i (0 .. $randtimes) {
	my $array;
	if ($i) {$array = &randomize(\@genomes)} else {$array = \@genomes}
	print OUT $#genomes+1, "\n";
	foreach my $x (@$array) {
		my $genomestring =  "genome".$x;
		my $addspaces = " " x (10 - length($genomestring));
		print OUT $genomestring,$addspaces;
		foreach my $y (@$array) {
			if ($matrix[$x][$y]) {print OUT $matrix[$x][$y], " "}
			else {print OUT "0 "}
			}
		print OUT "\n";
		}
	}
close OUT;
exit (0);


sub randomize {
	my $array = shift;
	for (my $i = @$array; --$i; ) {
		my $j = int rand ($i+1);
		next if ($i == $j);
		@$array[$i, $j] = @$array[$j, $i];
		}
	return $array;
	}
