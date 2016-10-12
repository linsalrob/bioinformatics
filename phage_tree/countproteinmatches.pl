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



use DBI;
use strict;

my $dbh=DBI->connect("DBI:mysql:phage", "SQLUSER", "SQLPASSWORD") or die "Can't connect to database\n";
my $dir= shift || &niceexit("countproteinmatches.pl <dir of prot dists>\n");

my @count; my %max;
# read each file one at a time, and add the data to an array
opendir(DIR, $dir) || &niceexit("Can't open $dir");
while (my $file = readdir(DIR)) {
	next if ($file =~ /^\./);
	open (IN, "$dir/$file") || &niceexit("Can't open $dir/$file\n");
	my $get;
	while (my $line = <IN>) {
		if ($line =~ /Sequences producing/) {$get = 1; next}
		next unless ($get);
		last if ($line =~ /^>/);
		next unless ($line =~ /\d_\d/);
		my @a = split (/\s+/, $line);
		my ($gene, $source) =split (/_/, $a[0]);
		unless ($gene && $source) {die "problem parsing $line\n"}
		$count[$source][$gene]++;
		unless ($max{$source}) {$max{$source}=$count[$source][$gene]}
		if ($count[$source][$gene] > $max{$source}) {$max{$source}=$count[$source][$gene]}
		
	}
	close IN;
}
closedir(DIR);

foreach my $source (0 .. $#count) {
	foreach my $gene (0 .. $#{$count[$source]}) {
		if ($count[$source][$gene]) {
			if ($count[$source][$gene] > 1) {print "$source\t$gene\t",$count[$source][$gene] - 1,"\n"}}
		}
	
	}
&niceexit();

#my $exc = $dbh->prepare("SELECT translation from protein where count = $gene" ) or croak $dbh->errstr;
#		$exc->execute or die $dbh->errstr;
#		while (my @retrieved = $exc->fetchrow_array) {print OUT ">$a[0]\n$retrieved[0]\n"}



sub niceexit {
	my $reason = shift;
	$dbh->disconnect;
	if ($reason) {print STDERR $reason; exit(-1)}
	else {exit(0)}
	}
