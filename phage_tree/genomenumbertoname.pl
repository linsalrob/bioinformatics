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


# get all the phage genome names and their appropriate number from the database
my $dbh=DBI->connect("DBI:mysql:phage", "SQLUSER", "SQLPASSWORD") or die "Can't connect to database\n";

my %genomename;

my $exc = $dbh->prepare("select count, organism, family from phage" ) or croak $dbh->errstr;
$exc->execute or die $dbh->errstr;
while (my @retrieved = $exc->fetchrow_array) {
	my $temp = "genome".$retrieved[0];
	$genomename{$temp} = $retrieved[1]." ".$retrieved[2];
	}

open (IN, $treefile) || &niceexit("Can't open $treefile\n");
open (OUT, ">$treefile.corrected")  || &niceexit("Can't open $treefile.corrected for writing\n");
while (<IN>) {
	unless (/genome\d+/) {print OUT; next}
	while (/(genome\d+)/) {s/(genome\d+)/$genomename{$1}/g}
	print OUT;
	}
&niceexit;



sub niceexit {
	my $reason = shift;
	$dbh->disconnect;
	if ($reason) {print STDERR $reason; exit(-1)}
	else {exit(0)}
	}

