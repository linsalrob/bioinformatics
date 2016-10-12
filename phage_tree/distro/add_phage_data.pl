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

# extract all the data from Genbank entries for mySQL.
# The only data we want from the start is the GI number.
# we will then cut down to the features, and get a bunch of information:

# from source we will get start and stop (should be the length), organism, and db_xref

# for each CDS, we want: 
#	start and stop position, and complement (?)
#	gene name
#	standard_name
#	note
#	product
#	protein_id
#	db_xref
#	translation
#	function

# until we get to origin. Then we want all the sequence until //

# note that this will run straight into mysql. You need to make sure that you allow
# all privileges for the user on the appropriate database.

use strict;
use DBI;
use Term::ReadKey;

my $usage = "add_phage <genbank file> <file of corrected names> <swiss file>\n";
my $genomefile = shift || die $usage;
my $goodnamesfile = shift || die $usage;
my $swissfile = shift || die $usage;


print "\n\nNow we need a user authorized to add data\nPlease enter the mySQL user name:  ";
my $sqluser = ReadLine 0;
chomp $sqluser;

print "Please enter the mySQL user password:  ";
ReadMode 2;
my $sqlpassword = ReadLine 0;
chomp $sqlpassword;
ReadMode 1;
print "\n";


my @line; my %locus; my %accession; my %from; my %to; my %counts;
	my %org; my %gene; my %function; my %product; my %proteinid;
	my %dbxref; my %translation; my %sequence; my %complement;
	my %standardname; my %note; my %start; my %end; my %organism;
	my %seenlocus; my %seenacc;


my $dbh = DBI->connect('DBI:mysql:phage', "$sqluser", "$sqlpassword") || die "Can't connect to the database\n";

print "\n\nREADING GENOME DATA\n";
open (IN, $genomefile) || die "Can't open $genomefile\n";
my $count=1; my $orgcount=0;
while (my $line=<IN>) {
	next until ($line =~ /VERSION/);
	$orgcount++;
	@line = split (/\s+/, $line);
	unless ($line[1]) {$line[1]=$orgcount} unless ($line[2]) {$line[2]=$orgcount}
	if ((exists $seenlocus{$line[1]}) || (exists $seenacc{$line[2]})) {die "$line[1], $line[2] already seen\n"}
	$seenlocus{$line[1]}=$seenacc{$line[2]}=1;
	$locus{$orgcount} = $line[1]; $accession{$orgcount}=$line[2];

	until ($line =~ /FEATURES/) {$line=<IN>}
	$line = <IN>;
	$line =~ s/\.\./  /;
	@line=split (/\s+/, $line);
	$from{$orgcount} = $line[$#line-1]; $to{$orgcount} = $line[$#line];
	if (($from{$orgcount} =~ /\D/) || ($to{$orgcount} =~ /\D/)) {print STDERR "WARNING ERROR WITH $from{$orgcount} and $to{$orgcount} from $line\n"}
	$line = <IN>;
	$line =~ /\"(.*)\"/; $org{$orgcount} = $1; $organism{$count} = $org{$orgcount};
	until (($line =~ /\s+CDS\s+/) || ($line =~ /BASE COUNT/)) {$line = <IN>}
	if ($line =~ /BASE COUNT/) {print STDERR "No CDS found for $org{$orgcount}\n"}
	else {
		if ($line =~ /complement/i) {$complement{$count}=1}
		$line =~ /(\d+)\.\.(\d+)/; $start{$count}=$1; $end{$count}=$2;
		push (@{$counts{$orgcount}}, $count);
		}
	until ($line =~ /\/\//) {
		$line = <IN>;
		if ($line=~ /\/gene=\"(\w+)\"/) {unless ($gene{$count}) {$gene{$count}=$1}}
		if ($line=~ /\/function=\"(\w+)\"/) {$function{$count}=$1}
		if ($line=~ /\/product=\"(\w+)\"/) {$product{$count}=$1}
		if ($line=~ /\/protein_id=\"(.*)\"/) {$proteinid{$count}=$1}
		if ($line=~ /\/db_xref=\"(.*)\"/) {$dbxref{$count}=$1}
		if ($line=~ /\/translation/) {
			$line =~ s/\/translation=\"//;
			$line =~ s/\s+//g;
			chomp($line);
			$translation{$count} = $line;
			until ($line =~ /\"/) {
				$line = <IN>;
				$line =~ s/\s+//g;
				chomp($line);
				$translation{$count} .= $line;
				}
			$translation{$count} =~ s/\"//g;
			}
		if ($line=~ /\/note=\"(.*)\"/) {$note{$count}=$1; $note{$count}=~ s/\'//g;}
		if ($line=~ /\/standard_name=\"(.*)\"/) {$standardname{$count}=$1}
		if ($line=~ /\s+CDS\s+/) {
			if ($line =~ />/) {$line =~ s/>//}
			if ($line =~ /</) {$line =~ s/<//}
			$count++;
			$organism{$count}=$organism{$count-1};
			push (@{$counts{$orgcount}}, $count);
			if ($line =~ /complement/i) {$complement{$count}=1}
			$line =~ /(\d+)\.\.(\d+)/; $start{$count}=$1; $end{$count}=$2;
			}
		if ($line =~ /ORIGIN/) {
			until ($line =~ /\/\//) {
				$line=<IN>;
				$line =~ s/\s+//g;
				$line =~ s/\d+//g;
				$sequence{$orgcount} .= $line;
				}
			$sequence{$orgcount} =~ s/\/\///g;
			}
		}
	$count++;
	}
close IN;

print STDERR "Reading good names\n";
my %goodname; my %genbankname;
open (IN, $goodnamesfile) || die "Can't open $goodnamesfile\n";
while (<IN>) {
	my @line=split /\t/;
	$goodname{$line[0]}=$line[1];
	$genbankname{$line[0]}=$line[2];
}
close IN;

print STDERR "Reading phylogeny\n";
my %phylo; my %phylofamily;
open (PHYLO, "phylogeny.txt") || die "can't open phylogeny.txt\n";
while (my $line = <PHYLO>) {
	chomp($line);
	my @parts = split (/\s+\-\-\>\s+/, $line);
	$phylo{$parts[$#parts]} = $line;
	if ($line =~ /\b(\w+viridae)\b/) {$phylofamily{$parts[$#parts]} = $1}
	elsif ($line =~ /unclassified/) {$phylofamily{$parts[$#parts]} = 'unclassified'}
	elsif ($line =~ /\b(\w+virus)\b/) {$phylofamily{$parts[$#parts]} = $1}
}
close PHYLO;




foreach my $org (sort {uc($org{$a}) cmp uc($org{$b})} keys %org) {
	unless (exists $phylo{$org{$org}}) {print STDERR "No phylogeny found for $org{$org}\n"}
	unless (exists $phylofamily{$org{$org}}) {print STDERR "No phylogeny family found for $org{$org}\n"}
	my $goodname; my $genbankname;
	if (exists $goodname{$accession{$org}}) {
		$goodname = $goodname{$accession{$org}};
		$genbankname = $genbankname{$accession{$org}};
	}
	else {
		$goodname = $genbankname = $org{$org};
		print STDERR "No Good name found for $accession{$org}, using $org{$org}\n";
	}

	my $newput ="INSERT INTO phage (count, organism, genbankname, accession, locus, beginning, end, sequence, phylogeny, family) ";
	$newput .= "VALUES ('NULL', '".$goodname."', '".$genbankname."', '".$accession{$org}."', '".$locus{$org}."', '".$from{$org};
	$newput .= "', '".$to{$org}."', '".$sequence{$org}."', '".$phylo{$org{$org}}."', '".$phylofamily{$org{$org}}."')";
	my $exc = $dbh->prepare("$newput") or die $dbh->errstr;
	$exc->execute or die $dbh->errstr;
	
	}

print STDERR "Nullifying data\n";
foreach my $org (keys %org) {
foreach my $key (@{$counts{$org}}) {
	unless (exists $gene{$key}) {$gene{$key}=""}
	unless (exists $function{$key}) {$function{$key}=""}
	unless (exists $note{$key}) {$note{$key}=""}
	unless (exists $product{$key}) {$product{$key}=""}
	unless (exists $proteinid{$key}) {$proteinid{$key}=""}
	unless (exists $dbxref{$key}) {$dbxref{$key}=""}
}
}

{
my $count=1;
print STDERR "\nADDING DATA\nCount\tORG\t# CDS\n";

foreach my $org (sort keys %org) {
	my $orgtotaladded;
	foreach my $key (@{$counts{$org}}) {
		my $exc = $dbh->prepare("SELECT count FROM phage WHERE genbankname LIKE '%$organism{$key}%'") or die $dbh->errstr;
		$exc->execute or die $dbh->errstr;

		my @retrieved = $exc->fetchrow_array;
		my $retrieved = $retrieved[0];

	#	print "Got $retrieved for $organism{$key} using $org\n";

		my $putdata = "INSERT INTO protein (count, organism, start, stop, complement, gene,";
		$putdata .= "function, note, product, proteinid, dbxref, translation)";

		$putdata .= " VALUES ('NULL', '".$retrieved."', '".$start{$key}."', '".$end{$key};
		if (exists $complement{$key}) {$putdata .= "', '1', '"} else {$putdata .= "', 'NULL', '"} 
		$putdata .= $gene{$key}."','".$function{$key}."','".$note{$key}."','".$product{$key}."','".$proteinid{$key}."','".$dbxref{$key};
		$putdata .= "','".$translation{$key}."')";

		$exc = $dbh->prepare($putdata) or die $dbh->errstr;
		$exc->execute or die $dbh->errstr;
		$orgtotaladded++;
		}
	print STDERR "$count\t$org{$org}\t$orgtotaladded\n";
	$count++;
}
}

print "\n\nADDING SWISS DATA\n\n";

open (SWISS, $swissfile) || die "Can't open $swissfile\n";

my $cols = <SWISS>;
$cols =~ s/\t/\,/g;

while (my $line = <SWISS>) {
	$line =~ s/\'//g;
	$line =~ s/^\d+\t/NULL\t/;
	$line =~ s/\t/\'\,\'/g;
	$line = "'".$line."'";
	my $exc = $dbh->prepare("INSERT INTO swiss ($cols) VALUES ($line)") or die $dbh->errstr;
	$exc->execute or die $dbh->errstr;
	}




$dbh->disconnect;

