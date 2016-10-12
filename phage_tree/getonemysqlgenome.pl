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
use DBI;

my $dbh=DBI->connect("DBI:mysql:phage", "SQLUSER", "SQLPASSWORD") or die "Can't connect to database\n";
my $genome = shift or &niceexit("which number genome to get?\n");
my $count;

open OUT, ">genome.$genome";
my $exc = $dbh->prepare("SELECT count, translation from protein where organism = $genome" ) or croak $dbh->errstr;
$exc->execute or die $dbh->errstr;
while (my @retrieved = $exc->fetchrow_array) {
	print OUT ">$retrieved[0]\n$retrieved[1]\n";
	$count++;
	}
close OUT;


$exc = $dbh->prepare("SELECT organism from phage where count = $genome" ) or croak $dbh->errstr;
$exc->execute or die $dbh->errstr;

my @retrieved;
while (@retrieved = $exc->fetchrow_array) {
print "Genome written to genome.$genome ($retrieved[0], $count proteins)\n";
}

&niceexit(0);

sub niceexit {
	my $reason = shift;
	$dbh->disconnect;
	if ($reason) {print STDERR $reason; exit(-1)}
	else {exit(0)}
	}

