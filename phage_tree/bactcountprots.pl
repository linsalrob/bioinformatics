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

# get all the phage genome names and their appropriate number from the database
my $dbh=DBI->connect("DBI:mysql:bacteria", "SQLUSER", "SQLPASSWORD") or die "Can't connect to database\n";

my %genomename;

my $exc = $dbh->prepare("select count, organism from bacteria" ) or croak $dbh->errstr;
$exc->execute or die $dbh->errstr;
while (my @retrieved = $exc->fetchrow_array) {
	my $exce = $dbh->prepare("select count from protein where organism=$retrieved[0]" ) or croak $dbh->errstr;
	$exce->execute or die $dbh->errstr;
	my @prots;
	while (my @ret = $exce->fetchrow_array) {push (@prots, $ret[0])}
	print "$retrieved[1]: ", $#prots+1, "\n";
}

$dbh->disconnect;
