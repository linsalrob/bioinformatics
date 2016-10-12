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

# get functions from mysql and then tree things with similar functions.

if (-e "products") {die "products dir already exists\n"} else {mkdir "products", 0755}

my $dbh=DBI->connect("DBI:mysql:phage", "SQLUSER", "SQLPASSWORD") or die "Can't connect to database\n";

my %product;

{
my $exc = $dbh->prepare("select product from protein" ) or croak $dbh->errstr;
$exc->execute or die $dbh->errstr;
while (my @retrieved = $exc->fetchrow_array) {
	next unless ($retrieved[0]);
	next if (($retrieved[0] =~ /unknown/i) || ($retrieved[0] eq "NULL") || ($retrieved[0] =~ /orf/i));
	$product{$retrieved[0]} ++;
	}
}

open (PRODS, ">products.txt") || die "can't open products.txt\n";
foreach my $product (sort {$product{$b} <=> $product{$a}} keys %product) {
	next unless ($product{$product}>1);
	print PRODS "$product\t$product{$product}\n";
	my $exc = $dbh->prepare("select count, organism, translation from protein where product = '$product'" ) or croak $dbh->errstr;
	$exc->execute or die $dbh->errstr;
	open (OUT, ">products/$product") || die "Can't open products/$product for writing\n";
	while (my @ret = $exc->fetchrow_array) {
		print OUT ">",$ret[0],"_",$ret[1],"\n",$ret[2],"\n";
		}
	}


unless (-e "products.clustal") {mkdir "products.clustal", 0755}
unless (-e "products.protdist") {mkdir "products.protdist", 0755}
unless (-e "yes") {open YES, ">yes"; print YES "y\n"; close YES}



opendir(DIR, "products") || &niceexit("Can't open products for reading\n");
while (my $file = readdir(DIR)) {
	next if ($file =~ /^\./);
	next if ($file =~ /\.dnd$/);
	system "/usr/local/genome/bin/clustalw -INFILE=products/$file -OUTFILE=products.clustal/infile -OUTPUT=PHYLIP";
	chdir "products.clustal";
	system "/usr/local/genome/bin/protdist < ../yes";
	system "mv infile $file";
	system "mv outfile ../products.protdist/$file";
	chdir "..";
	}

&niceexit();





sub niceexit {
	my $reason = shift;
	$dbh->disconnect;
	print STDERR "Done in ", time-$^T, " seconds\n";
	if ($reason) {print STDERR $reason; exit(-1)}
	else {exit(0)}
	}



	
