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



# compareprotdists.pl

# compare whether the protein distances are always the same across the alignments.
use DBI;
use strict;

my $usage = "compareprotdists.pl <dir of prot dists> <number of genomes used>\n";


my $dbh=DBI->connect("DBI:mysql:phage", "apache") or die "Can't connect to database\n";
my $dir= shift || &niceexit($usage);
my $nogenomes = shift || &niceexit($usage);


my @matrix; my @files;
my @count;
my @proteinmatches;
my %linegenomecount;
# read each file one at a time, and add the data to an array
opendir(DIR, $dir) || &niceexit("Can't open $dir");
while (my $file=readdir(DIR)) {
	next if ($file =~ /^\./);
	open (IN, "$dir/$file") || &niceexit("Can't open $dir/$file");
	my @genomes = ('0');
	my @genes = ('0');
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
		
		push (@genes, $gene);
		push (@genomes, $genome);
		$checkdupgenomes{$genome}++;
		}


	# now we loop through all the lines, and split them on white space. 
	# then we make an array with each value to the pre-existing value in the matrix
	# note that because the genomes are represented as numbers we can just
	# use these numbers for the position in the matrix.
	
	# this version is fixed to count one half of the matrix and ignore the zero diagonal
	
	# note that column zero of the matrix is empty (there is no genome 0)
	foreach my $z (0 .. $#dists) {
		my @line = split (/\s+/, $dists[$z]);
		unless ($#line == $#genomes) {
			my $x; foreach my $y (0 .. $#dists) {if ($dists[$y] eq $dists[$z]) {$x = $y}}
			&niceexit("PROBLEM WITH \n@line AND \n@genomes\nIN $file BECAUSE $#line AND $#genomes AT LINE $x\n");
			}
		my ($gene, $linegenome) = split (/_/, $line[0]);
		unless ($gene && $linegenome) {&niceexit("CAN'T PARSE @line SECOND TIME AROUND\n")}
		foreach my $x (1 .. $#genomes) {
			next if ($x <= $z+1);
			my ($low, $high) = ($gene, $genes[$x]);
			if ($low > $high) {($low, $high) = ($high, $low)}
			push (@{$matrix[$low][$high]}, $line[$x]);
			push (@{$files[$low][$high]}, $file);
			}
		}
	}



my $totaldiff; my $diffcount; my $totalsame; my $maxdiff=1; my $mindiff=100; my $maxprot; my $minprot;
# now we will loop through the matrix and look for differences in each array.
foreach my $y (1 .. $#matrix) {
  next unless ($matrix[$y]);
  foreach my $x (1 .. $#{$matrix[$y]}) {
    next if ($y == $x);
    next unless ($matrix[$y][$x]);
    my @matches = @{$matrix[$y][$x]};
    print "Checking $y and $x (", $#matches+1," records, ", join (" ", @{$files[$y][$x]}), " files): ";
    my ($same, $diff) = (0, 0);
    for my $c (0 .. $#matches) {
      for my $d ($c .. $#matches) {
	next if ($c == $d);
	unless ($matches[$c] == $matches[$d]) {
	my ($trash, $file1) = split /_/, ${$files[$y][$x]}[$c];
	my ($trash, $file2) = split /_/, ${$files[$y][$x]}[$d];
	next if ($file1 == $file2);
	  my $diffvalue =  abs($matches[$c] - $matches[$d]);
	  $totaldiff+= $diffvalue; $diffcount++;
	  if ($diffvalue > $maxdiff) {$maxdiff = $diffvalue; $maxprot = ${$files[$y][$x]}[$c]." and ".${$files[$y][$x]}[$d]}
	  if ($diffvalue < $mindiff) {$mindiff = $diffvalue; $minprot =  ${$files[$y][$x]}[$c]." and ".${$files[$y][$x]}[$d]}
	  print "\n\tDifferent in ", ${$files[$y][$x]}[$c]," and ", ${$files[$y][$x]}[$d], " because $matches[$c] and $matches[$d]: difference $diffvalue";
	  $diff++;
	  }
	 else {$same++; $totalsame++}
	}
       }
     if ($diff) {print "\n"}
     print "SAME: $same DIFF: $diff\n";
   }
  }
 
print STDERR "Total difference: $totaldiff, number different: $diffcount, number same: $totalsame, average difference: ", $totaldiff/$diffcount, "\n";
print "Total difference: $totaldiff, number different: $diffcount, number same: $totalsame, average difference: ", $totaldiff/$diffcount, "\n";
print "Maximum difference in scores for a single protein: $maxdiff (for $maxprot) and Minimum difference in scores for a single protein: $mindiff (for $minprot).\n";
print STDERR "Maximum difference in scores for a single protein: $maxdiff (for $maxprot) and Minimum difference in scores for a single protein: $mindiff (for $minprot).\n";





&niceexit(0);






sub niceexit {
	my $reason = shift;
	$dbh->disconnect;
	if ($reason) {print STDERR $reason; exit(-1)}
	else {exit(0)}
	}
