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



use strict;

my $dir = shift || die "Need a dir of prot dists\n";

my %max; my $max =1;
my $min = 10000;

opendir(DIR, $dir);
while (my $file = readdir(DIR)) {
	next if ($file =~ /^\./);
	open (IN, "$dir/$file") || die "Can't open $dir/$file\n";
	while (<IN>) {
		next unless (/^\d+_\d/);
		my @a = split /\s+/;
		foreach my $x (1 .. $#a) {
			if ($a[$x] > $max) {$max{$file} = $a[$x]; $max=$a[$x]}
			#if ($a[$x] < $min) {$min = $a[$x]}
			}
		}
	}


my @max = sort {$max{$a} <=> $max{$b}} keys %max;
print "Max: $max{$max[$#max]}\n";
foreach my $max (@max) {
	print "$max\t$max{$max}\n";
	}
print "\n";
