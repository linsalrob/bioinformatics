#!/usr/bin/perl 

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



# extract all the data from Genbank entries for mySQL.
# The only data we want from the start is the GI number.
# we will then cut down to the features, and get a bunch of information:

# This will add all genbank data from a single file, separated by //

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


my @line; my %locus; my %accession; my %from; my %to; my %counts;
	my %org; my %gene; my %function; my %product; my %proteinid;
	my %dbxref; my %translation; my %sequence; my %complement;
	my %standardname; my %note; my %start; my %end; my %organism;
	

my $genomefile = shift || die "Usage: gbk2sql.pl <single gbk file with genomes> 1 (to do it)\n";
my $do = shift;



my $dbh = DBI->connect('DBI:mysql:phage', "SQLUSER", "SQLPASSWORD");

#my $dbh = DBI->connect('DBI:mysql:phage', 'apache');


print "\n\nADDING GENOME DATA\n";
open (IN, $genomefile) || die "Can't open $genomefile\n";
my $count=1; my $orgcount=0;
while (my $line=<IN>) {
	next until ($line =~ /VERSION/);
	$orgcount++;
	@line = split (/\s+/, $line);
	$locus{$orgcount} = $line[1]; $accession{$orgcount}=$line[2];
	until ($line =~ /FEATURES/) {$line=<IN>}
	$line = <IN>;
	$line =~ s/\.\./  /;
	@line=split (/\s+/, $line);
	$from{$orgcount} = $line[$#line-1]; $to{$orgcount} = $line[$#line];
	if (($from{$orgcount} =~ /\D/) || ($to{$orgcount} =~ /\D/)) {print STDERR "WARNING ERROR WITH $from{$orgcount} and $to{$orgcount} from $line\n"}
	$line = <IN>;
	until ($line=~ /\".*\"/) {chomp($line); $line .= <IN>}
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
	my $newput ="INSERT INTO phage (count, organism, accession, locus, beginning, end, sequence, phylogeny, family) ";
	$newput .= "VALUES ('NULL', '".$org{$org}."', '".$accession{$org}."', '".$locus{$org}."', '".$from{$org};
	$newput .= "', '".$to{$org}."', '".$sequence{$org}."', '".$phylo{$org{$org}}."', '".$phylofamily{$org{$org}}."')";
	if ($do) {
		my $exc = $dbh->prepare("$newput") or die $dbh->errstr;
		$exc->execute or die $dbh->errstr;
	}
	else {print "$newput\n"}
	
	}

print STDERR "\n\nORG\t# CDS\n";

foreach my $org (sort keys %org) {
	my $orgtotaladded;
	foreach my $key (@{$counts{$org}}) {
		my $retrieved;
		if ($do) {
			my $exc = $dbh->prepare("SELECT count FROM phage WHERE organism LIKE '%$organism{$key}%'") or die $dbh->errstr;
			$exc->execute or die $dbh->errstr;
		
			my @retrieved = $exc->fetchrow_array;
			$retrieved = $retrieved[0];

		#print "Got $retrieved for $organism{$key} using $org\n";
		}
		else {$retrieved="xxx"}
		
		my $putdata = "INSERT INTO protein (count, organism, start, stop, complement, gene,";
		$putdata .= "function, note, product, proteinid, dbxref, translation)";

		$putdata .= " VALUES ('NULL', '".$retrieved."', '".$start{$key}."', '".$end{$key};
		if (exists $complement{$key}) {$putdata .= "', '1', '"} else {$putdata .= "', 'NULL', '"} 
		$putdata .= $gene{$key}."','".$function{$key}."','".$note{$key}."','".$product{$key}."','".$proteinid{$key}."','".$dbxref{$key};
		$putdata .= "','".$translation{$key}."')";
		
		if ($do) {
			my $exc = $dbh->prepare($putdata) or die $dbh->errstr;
			$exc->execute or die $dbh->errstr;
		}
		else {print "$putdata\n"}
		
		$orgtotaladded++;
		}
	print STDERR "$org{$org}\t$orgtotaladded\n";
}


$dbh->disconnect;

