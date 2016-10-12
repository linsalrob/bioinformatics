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

my $treefile = shift || die "Need a tree file to work with\n";

my %color;

$color{'Tectiviridae'} = "#000FF";
$color{'Leviviridae'} = "#B22222";
$color{'Plasmaviridae'} = "#9932CC";
$color{'Inoviridae'} = "#32CD32";
$color{'Fuselloviridae'} = "#FF7F50";
$color{'Corticoviridae'} = "#006400";
$color{'Myoviridae'} = "#778899";
$color{'Podoviridae'} = "#FF0000";
$color{'Microviridae'} = "#B8860B";
$color{'Siphoviridae'} = "#20B2AA";


# get all the phage genome names and their appropriate number from the database
my $dbh=DBI->connect("DBI:mysql:phage", "SQLUSER", "SQLPASSWORD") or die "Can't connect to database\n";

my %genomename;

my $exc = $dbh->prepare("select count, organism, family from phage" ) or croak $dbh->errstr;
$exc->execute or die $dbh->errstr;
while (my @retrieved = $exc->fetchrow_array) {
	my $temp = "genome".$retrieved[0];
	$genomename{$temp} = $retrieved[1]." (".$retrieved[2].")";
	}

open (IN, $treefile) || &niceexit("Can't open $treefile\n");
open (OUT, ">$treefile.html")  || &niceexit("Can't open $treefile.corrected for writing\n");
my $table;
print OUT "<html><title>Phage Trees</title></head><body bgcolor=\"#FFFFFF\">\n<pre>\n"; 

while (<IN>) {
	s/(genome\d+)/$genomename{$1}/g;
	if (/remember/i) {print OUT "</pre>\n$_"; next}
	if (/between/i) {$table=1; print OUT "<table>\n"; s/^/     /}
	if ($table) {
		next if (/\-\-\-/);
		s/\s+/<tr><td>/;
		s/\s*$/<\/td><\/tr>\n/;
		s/\s\s+/<\/td><td>/g;
	}
	foreach my $key (keys %color) {if (/$key/) {s/^/<font color=\"$color{$key}\">/; s/$/<\/font>/}}
	print OUT;
	}

print OUT "</table><p>\n";
print OUT "\n\n<p><p><hr><p><p>\n\nCOLOR CODES:\n<p>\n";
foreach my $key (sort {$a cmp $b} keys %color) {print OUT "<font color=\"$color{$key}\">$key</font><br>\n"}

print OUT "</body></html>\n";




&niceexit;



sub niceexit {
	my $reason = shift;
	$dbh->disconnect;
	if ($reason) {print STDERR $reason; exit(-1)}
	else {exit(0)}
	}

