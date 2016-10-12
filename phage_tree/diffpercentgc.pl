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



#generate a tree based on the difference in percent GC

use strict;
use DBI;

my $dbh=DBI->connect("DBI:mysql:phage", "SQLUSER", "SQLPASSWORD") or die "Can't connect to database\n";

#calculate the percentgc
my %percent;
{
my $exc = $dbh->prepare("SELECT count, sequence from phage") or croak $dbh->errstr;
$exc->execute or die $dbh->errstr;
while (my @retrieved = $exc->fetchrow_array) {
	my $seq = $retrieved[1];
	my $gc; my $at;
	while ($seq) {
		my $base = chop($seq);
		if ((uc($base) eq "G") || (uc($base) eq "C")) {$gc++} else {$at++}
	}
	$percent{$retrieved[0]} = ($gc/($gc+$at)) *100;
	}
}
my $max;
foreach my $x (keys %percent) {
	foreach my $y (keys %percent) {
		if (abs($percent{$x} - $percent{$y}) > $max) {$max = abs($percent{$x} - $percent{$y})};
		}
	}


my @diff;
foreach my $x (keys %percent) {
	$diff[$x][0]='';
	foreach my $y (keys %percent) {
		$diff[$x][$y] = (abs($percent{$x} - $percent{$y})/$max)*100;
		}
	}

print "   ", $#diff, "\n";
foreach my $x (sort {$a <=> $b} keys %percent) {
	my $name = "genome".$x;
	my $space = " " x (9-(length($name)));
	$name .= $space;
	print "$name", join ("  ", @{$diff[$x]}), "\n";
	}

$dbh->disconnect;
