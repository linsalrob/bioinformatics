#!/usr/bin/perl -w


#    Copyright 2001, 20002 Rob Edwards
#    For updates, more information, or to discuss the scripts
#    please contact Rob Edwards at redwards@utmem.edu or via
#    http://www.salmonella.org/
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
#    in the file (COPYING) along with these scripts; if not, write to the Free
#    Software
#    Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA
#

# make fasta databases.

# extract the data from mysql and write it to a fasta file

use strict;
use DBI;
use Term::ReadKey;


print "\n\nWe need a user authorized to access data\nPlease enter the mySQL user name:  ";
my $sqluser = ReadLine 0;
chomp $sqluser;

print "Please enter the mySQL user password:  ";
ReadMode 2;
my $sqlpassword = ReadLine 0;
chomp $sqlpassword;
ReadMode 1;
print "\n";

my $dbh = DBI->connect('DBI:mysql:phage', "$sqluser", "$sqlpassword");

my $exc = $dbh->prepare("SELECT count,sequence FROM phage") or die $dbh->errstr;
$exc->execute or die $dbh->errstr;

my $count; my $sequence; my $org;

open (OUT, ">/seqs/databases/phage.nt.dbs") or die "Can't open /seqs/databases/phage.dna.dbs for writing\n";
while (($count, $sequence)= $exc->fetchrow_array) {
	$sequence =~ s/(.{60})/$1\n/g;
	print OUT ">$count\n$sequence\n"}
close (OUT);

$exc = $dbh->prepare("SELECT count,organism,translation FROM protein") or die $dbh->errstr;
$exc->execute or die $dbh->errstr;

open (OUT, ">/seqs/databases/phage.aa.dbs") or die "Can't open /seqs/databases/phage.aa.dbs for writing\n";
while (($count, $org, $sequence)= $exc->fetchrow_array) {print OUT ">", $count, "_", $org, "\n", $sequence, "\n"}

$exc = $dbh->prepare("SELECT count,ac,seq from swiss") or die $dbh->errstr;
$exc->execute or die $dbh->errstr;
while (($count, $org, $sequence)= $exc->fetchrow_array) {
	$org =~ s/\W//g;
	print OUT ">", $count, "__", $org, "\n", $sequence, "\n"}


close (OUT);


system ('pressdb /seqs/databases/phage.nt.dbs');
system ('setdb /seqs/databases/phage.aa.dbs');

$dbh->disconnect;
