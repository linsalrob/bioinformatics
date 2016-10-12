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

my $file = shift || die "Need a fasta file from artemis to parse\n";
my $org = shift || die "What number organism is this?\n";
my $do = shift; unless ($do) {$do = '-1'}

my $dbh = DBI->connect('DBI:mysql:phage', "SQLUSER", "SQLPASSWORD");


my %translation;
{
open (IN, $file) || die "Can't open $file\n";
my $tag; my $seq;
while (<IN>) {
	chomp;
	if (/^>/) {
		if ($seq) {
			$translation{$tag} = $seq;
			undef $seq;
			}
		s/^>//;
		$tag = $_;
		}
	else {$seq .= $_}
	}
$translation{$tag} = $seq;
}

my @keys = keys %translation;

my $count;
foreach my $key (keys %translation) {
print "PROCESSING $key\n";
	#orf.10, undefined product 6768:7079 forward MW:10883 
	unless ($key =~ /^orf/) {print STDERR "Error with $key, not added\n"; next}
	$key =~ /(orf\.\d+).*?,\s(.*?) (\d+)\:(\d+)\s+(\w+)\s+MW/;
	my ($gene, $function, $start, $stop, $dir) = ($1, $2, $3, $4, $5);
	my $complement;
	unless ($dir eq "forward") {$complement = 1} else {$complement=0}
	print "Adding:\n\tgene: $gene\n\tfunc: $function\n\tstart: $start\n\tstop: $stop\n\tdir: $complement\n";
	print "Translation: $translation{$key}\n";
	my $putdata = "INSERT INTO protein (count, organism, start, stop, complement, gene, function, translation)";
	$putdata .= " VALUES ('NULL', '".$org."', '".$start."', '".$stop."', '".$complement."', '".$gene;
	$putdata .= "', '".$function."', '".$translation{$key}."')";
	
	print "\n$putdata\n\n\n";
	$count++;
	if ($do == 2) {
		my $exc = $dbh->prepare($putdata) or die $dbh->errstr;
		$exc->execute or die $dbh->errstr;
		}
}


print $#keys+1, " entries and $count added to database\n";
$dbh->disconnect;
