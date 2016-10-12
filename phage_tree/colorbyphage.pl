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

my %family;
my %phage;

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
$color{'unclassified'} = "#000000";



print "<html><head><title>Colors and Phage</title></head></body>\n";

my $exc = $dbh->prepare("select organism, family from phage");
$exc->execute or die $dbh->errstr;

while (my @ret = $exc->fetchrow_array) {
	$family{$ret[0]} = $ret[1];
	push (@{$phage{$ret[1]}}, $ret[0]);
	}

foreach my $family (sort {uc($a) cmp uc($b)} keys %phage) {
	print "<font color=\"$color{$family}\"><b>$family</b><br><ul>\n<li>", join ("\n<li>", @{$phage{$family}}), "\n</ul><p>\n";
	}
print "</font>\n</ul><p><hr><p>\n";

print "<table border=1><p>";
foreach my $phage (sort {uc($a) cmp uc($b)} keys %family) {
	print "<tr><td><font color=\"", $color{$family{$phage}}, "\">$phage</td><td>", 
	  "<font color=\"", $color{$family{$phage}}, "\">$family{$phage}<font></td></tr>\n";
	}

print "</body></html>\n";

$dbh->disconnect;

