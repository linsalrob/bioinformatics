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

my $file =shift || die "Need a file\n";

open (IN, $file) || die "Can't open $file\n";

my $dbh = DBI->connect('DBI:mysql:phage', "SQLUSER", "SQLPASSWORD");

my %function;

my %two;
$two{"NP_049691"} = 1;
$two{"NP_049656"} = 1;
$two{"NP_049799"} = 1;
$two{"NP_049849"} = 1;
$two{"NP_059623"} = 1;
$two{"AAA98578"} = 1;
$two{"AAA98599"} = 1;


while (my $line = <IN>) {
	chomp($line);
	if ($line =~ /^\d+\:\s+(\S+)$/) {
		my $id = $1;
		if (exists $two{$id}) {$id = $id.".2"} else {$id=$id.".1"}
		unless ($id)  {&niceexit("Can't get an id at $line\n")}
		my $exc = $dbh->prepare("select count from protein where proteinid = '$id'") or &niceexit("$dbh->errstr");
		$exc->execute or &niceexit("$dbh->errstr");
		my @retrieved = $exc->fetchrow_array;
		unless ($retrieved[0]) {&niceexit("Can't find $id\n")}
		
		if ($function{$id}) {&niceexit("Function $function{$id} exists for $id at $line\n")}
		$line = <IN>;
		if ($line =~ /(.*)\[/) {
			$function{$id} = $1;
			$function{$id} =~ s/\s+$//;
			$function{$id} =~ s/^\s+//;
			$function{$id} =~ s/\~/ /g;
			}
		else {&niceexit("Can't get a function from line: $line\n")}
		}
	elsif ($line =~ /^\d+/) {die "Can't get an id for $line\n"}
	}

my @ids = keys %function;
print "Found ", $#ids+1, " functions\n";
print "Adding data to database\n";
my $funcadded;
foreach my $id (@ids) {
	my $prepare = "update protein set function='".$function{$id}."' where proteinid = '".$id."'";
	my $exc = $dbh->prepare($prepare) or &niceexit("$dbh->errstr");
	$exc->execute or &niceexit("$dbh->errstr");
	$funcadded++;
	}

print "$funcadded functions added\n";


&niceexit(0);


sub niceexit {
	my $reason = shift;
	$dbh->disconnect;
	if ($reason) {print STDERR $reason; exit(-1)}
	else {exit(0)}
	}
