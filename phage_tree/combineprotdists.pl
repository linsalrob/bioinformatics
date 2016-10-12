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



# combinprotdists.pl

# new version. This will assign a distance of 100 to any sequence that does not match, and to all empty
# spaces. There is a good rationale for this. The The Dayhoff PAM matrix scoring system returns a percent of
# the amino acids that are likely to have changed. Therefore a 100% score means that they have all changed.

# We will make an average, and include the number of scores used to calculate the average. Then we can fitch
# it with the subreplicas option.

# the subreplicate number will be added if the protein appears, but is not similar. (i.e. an average score of 100 2
# means that two proteins were found to be similar but were to distant for a score. But an average score of
# 100 0 means that no proteins were found!

# I may change this, and add 100 to all genomes that don't match, not sure yet. This may be unnecessarily
# penalizing genomes that have something very distantly related.

use DBI;
use strict;

my $dbh=DBI->connect("DBI:mysql:phage", "SQLUSER", "SQLPASSWORD") or die "Can't connect to database\n";
my $dir= shift || &niceexit("combineprotdists.pl <dir of prot dists> <number of genomes used>\n");
my $nogenomes = shift || &niceexit("combineprotdists.pl <dir of prot dists> <number of genomes used>\n");

my @matrix;
my @count;
my %linegenomecount;
# read each file one at a time, and add the data to an array
opendir(DIR, $dir) || &niceexit("Can't open $dir");
while (my $file=readdir(DIR)) {
	next if ($file =~ /^\./);
	open (IN, "$dir/$file") || &niceexit("Can't open $dir/$file");
	my @genomes = ('0');
	my @dists;
	my %genomeshash;
	my %checkdupgenomes;
	my %genegenome;
	# we need to know all the genomes before we can store the data. Therefore
	# read and store each line in @dists
	# then get all the genome numbers and store them in an array
	while (my $line = <IN>) {
		chomp($line);
		my @line = split (/\s+/, $line);
		next unless ($#line);
		next if ($line =~ /^\s+/);
		# added loop to get the genome number from the database
		unless ($line[0] =~ /_/) {
			unless ($genegenome{$line[0]}) {
				my $exc = $dbh->prepare("select organism from protein where count like '$line[0]'");
				$exc->execute or die $dbh->errstr;
				my @retrieved = $exc->fetchrow_array;
				$genegenome{$line[0]} = $retrieved[0];
				}
			$line[0] .= "_".$genegenome{$line[0]};
			$line = join ("  ", @line); # this just corrects $line if we change it
		}
		push (@dists, $line);
		my ($gene, $genome) = split (/_/, $line[0]);
		unless ($gene && $genome) {&niceexit("Can't parse $line in $file\n")}
		
		push (@genomes, $genome);
		$checkdupgenomes{$genome}++;
		}


	# now we loop through all the lines, and split them on white space. 
	# then we add each value to the pre-existing value in the matrix
	# note that because the genomes are represented as numbers we can just
	# use these numbers for the position in the matrix.
	# we are going to also count the number of times that we save each data
	# point for the final average.
	# Finally, we only do this in one direction because the input
	# matrices are complete (and identical) on both halves.
	
	# note that column zero of the matrix is empty (there is no genome 0)
	foreach my $dist (@dists) {
		my @line = split (/\s+/, $dist);
		unless ($#line == $#genomes) {
			my $x; foreach my $y (0 .. $#dists) {if ($dists[$y] eq $dist) {$x = $y}}
			&niceexit("PROBLEM WITH \n@line AND \n@genomes\nIN $file BECAUSE $#line AND $#genomes AT LINE $x\n");
			}
		my ($gene, $linegenome) = split (/_/, $line[0]);
		unless ($gene && $linegenome) {&niceexit("CAN'T PARSE @line SECOND TIME AROUND\n")}
$linegenomecount{$linegenome}++;
		foreach my $x (1 .. $#genomes) {
			if ($line[$x] == -1) {$line[$x] = 100}
			if ($genomes[$x] == $linegenome) {$line[$x] = '0.000'}
			if ($matrix[$linegenome][$genomes[$x]]) {
				$matrix[$linegenome][$genomes[$x]] += $line[$x];
				$count[$linegenome][$genomes[$x]] ++;
				}
			else {
				$matrix[$linegenome][$genomes[$x]] = $line[$x];
				$count[$linegenome][$genomes[$x]] ++;
				}
			}
		# now we need to pad out all the missing genomes with 100's
		foreach my $x (1 .. $nogenomes) {
			next if ($checkdupgenomes{$x});
			if ($matrix[$linegenome][$x]) {
				$matrix[$linegenome][$x] += 100;
				$count[$linegenome][$x] ++;
				}
			else {
				$matrix[$linegenome][$x] = 100;
				$count[$linegenome][$x] ++;
				}
			}
			
		}
	}


# now we will average the matrix based on the count.
foreach my $y (0 .. $#matrix) {
	next unless ($matrix[$y]);
	foreach my $x (1 .. $#{$matrix[$y]}) {
		next unless ($count[$y][$x] && $matrix[$y][$x]);
		$matrix[$y][$x] = $matrix[$y][$x]/$count[$y][$x];
	}
}

# now we have all the data, lets just print out the matrix
print $#matrix, "\n";
#foreach my $y (1 .. $#matrix) {print STDERR "\t$y"}
#print STDERR "\n";
foreach my $y (1 .. $#matrix) {
	my $tempstring = "genome".$y;
	if (length($tempstring) > 10) {print STDERR "$tempstring is too long\n"}
	my $spacestoadd = " " x (10 - length($tempstring));
	print $tempstring,$spacestoadd;
	foreach my $x (1 .. $#matrix) {
		unless (defined $matrix[$y][$x]) {print "100 0  "; next}
		unless ($matrix[$y][$x]) {
			print "0 ";
			if ($count[$y][$x]) {print "$count[$y][$x]  "}
			else {print "0  "}
			next;
		}
		print $matrix[$y][$x], " ", $count[$y][$x], "  ";
		}
	print "\n";
	}
			


&niceexit(0);






sub niceexit {
	my $reason = shift;
	$dbh->disconnect;
	if ($reason) {print STDERR $reason; exit(-1)}
	else {exit(0)}
	}
