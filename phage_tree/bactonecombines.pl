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



# a new version of combine protein distances. This one will only use the target gene
# line for the matrix. That way, each protein is counted once against it's peers, and we
# should avoid having 1,000,000's of comparisons. Esp. important with bacterial genomes

# We will make an average, and include the number of scores used to calculate the average. Then we can fitch
# it with the subreplicas option.

# the subreplicate number will be added if the protein appears, but is not similar. (i.e. an average score of 100 2
# means that two proteins were found to be similar but were to distant for a score. But an average score of
# 100 0 means that no proteins were found!

# in this version you can select -n for no padding of missing proteins. Padding runs through the genome and if
# no match to another genome is found it increments the score by 100 for each protein in the query (line) genome

use DBI;
use strict;

my $usage = "combineprotdists.pl <dir of prot dists> <number of genomes used> <options>\nOPTIONS\n";
$usage .= "\t-n DON'T pad missing proteins with worst score (supplied with the -p option)\n\t-m print out all protein matches";
$usage .= "\n\t-p # penalty for being bad Rob. Default is 100\n";


my $dbh=DBI->connect("DBI:mysql:bacteria", "SQLUSER", "SQLPASSWORD") or die "Can't connect to database\n";
my $dir= shift || &niceexit($usage);
my $nogenomes = shift || &niceexit($usage);

my $args= join (" ", @ARGV);
my $pad=1; my $penalty =100;
if ($args =~ /-n/) {$pad=0}
if ($args =~ /-p\s+(\d+)/) {$penalty=$1}
print STDERR "Using PENALTY of $penalty and PAD of $pad\n";

my %noorfs;
&getnoorfs;

my @matrix; my @count; my @protein;


{
my $filecount; my $oldtime=$^T;
# read each file one at a time, and add the data to an array
opendir(DIR, $dir) || &niceexit("Can't open $dir");
while (my $file=readdir(DIR)) {
	my %hit; my @genomes= '0'; my $hitgenome; my %proteinhits;
	next if ($file =~ /^\./);
	open (IN, "$dir/$file") || &niceexit("Can't open $dir/$file");
	$filecount++; unless ($filecount % 1000) {print STDERR "$filecount\t", time-$^T, "\t", time-$oldtime, "\n"; $oldtime = time}
	# we need to know all the genomes before we can store the data. Therefore
	# read and store each line in @dists
	# then get all the genome numbers and store them in an array
	while (<IN>) {
		next if (/^\s+/);
		chomp;
		my @line = split;
		next unless ($#line);
		# store a reference to the array in a hash. Now we can retrieve the array that matches the filename!
		my ($gene, $genome) = split (/_/, $line[0]);
		unless ($gene && $genome) {&niceexit("Can't parse $_ in $file\n")}
		$hit{$line[0]}=\@line;
		push (@genomes, $genome);
		push (@{$proteinhits{$genome}}, $gene);
		if ($line[0] eq $file) {$hitgenome = $genome}
		}
	close IN;
	
	&niceexit("Can't find the genome match for $file\n") unless ($hitgenome);
	
	# now we have all the data, and all the genomes. We just need to get the line that 
	# matches the file name, and set up the matrix.
	
	my @match = @{$hit{$file}};

	&niceexit("Can't find the line for $file\n") unless (@match);
	my %seen;
	$seen{$hitgenome}=1;
	foreach my $x (1 .. $#match) {
		next if ($genomes[$x] == $hitgenome);
		$matrix[$hitgenome][$genomes[$x]]+=$match[$x];
		$count[$hitgenome][$genomes[$x]]++;
		$seen{$genomes[$x]}=1;
		# if we just push everything into one array and then sort it later, the array gets
		# too big for the memory. Lets sort everything here.
		{
			# remove any duplicate protein hits
			my %seenproteins;
			foreach (@{$protein[$hitgenome][$genomes[$x]]}) {$seenproteins{$_}=1}
			foreach (@{$proteinhits{$genomes[$x]}}) {$seenproteins{$_}=1}
			@{$protein[$hitgenome][$genomes[$x]]} = keys %seenproteins;
		}
	
	}
	if ($pad) {
		foreach my $x (1 .. $nogenomes) {
			next if (exists $seen{$x});
			$matrix[$hitgenome][$x]+=$penalty;
			$count[$hitgenome][$x]++;
			}
		}
	}

}

# now add a penalty for each genome for missing proteins
{
	foreach my $y (0 .. $#matrix) {
		next unless ($matrix[$y]);
		foreach my $x (1 .. $#{$matrix[$y]}) {
			next if ($y == $x);
			my $proteinsseen = $#{$protein[$y][$x]}+1;
print STDERR "For $y and $x: $proteinsseen proteins seen\n";
			my $diff = $noorfs{$y} - $proteinsseen;
			$matrix[$y][$x]+=($diff * $penalty);
			$count[$y][$x]+=$diff;
		}
	}
}



{

my %seen;
# now we will average the matrix based on the count.
foreach my $y (0 .. $#matrix) {
	next unless ($matrix[$y]);
	foreach my $x (1 .. $#{$matrix[$y]}) {
		next unless ($count[$y][$x] && $matrix[$y][$x]);
		my $temp = $x."+".$y; my $temp1 = $y."+".$x;
		next if ($seen{$temp} || $seen{$temp1});
		$seen{$temp} = $seen{$temp1} =1;

		if ($matrix[$y][$x] && $matrix[$x][$y]) {
			unless ($matrix[$y][$x] == $matrix[$x][$y]) {
			  $matrix[$y][$x] = $matrix[$x][$y] = $matrix[$y][$x] + $matrix[$x][$y];
			  $count[$y][$x] = $count[$x][$y] = $count[$y][$x] + $count[$x][$y];
			  }
		}
		elsif ($matrix[$y][$x]) {
			unless ($matrix[$x][$y]) {
				$matrix[$x][$y] = $matrix[$y][$x];
				$count[$x][$y] = $count[$y][$x];
				}
			}
		elsif ($matrix[$x][$y]) {
			unless ($matrix[$y][$x]) {
				$matrix[$y][$x] = $matrix[$x][$y];
				$count[$y][$x] = $count[$x][$y];
				}
			}
		else {&niceexit("Can't figure out matrix $y and $x ($matrix[$x][$y] and $matrix[$y][$x] when merging 3\n")}		
		# finally, take the average!
		$matrix[$x][$y] = $matrix[$y][$x] = $matrix[$y][$x]/$count[$y][$x];
		}
	}
}		

{
# now we have all the data, lets just print out the matrix
print $#matrix, "\n";
foreach my $y (1 .. $#matrix) {
	my $tempstring = "genome".$y;
	if (length($tempstring) > 10) {print STDERR "$tempstring is too long\n"}
	my $spacestoadd = " " x (10 - length($tempstring));
	print $tempstring,$spacestoadd;
	foreach my $x (1 .. $#matrix) {
		if ($y == $x)  {print "0 $noorfs{$x}  "; next}
		unless (defined $matrix[$y][$x]) {print "$penalty 0  "; next}
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
}



&niceexit(0);

sub getnoorfs {
	my $exc = $dbh->prepare("select organism from protein");
	$exc->execute or die $dbh->errstr;
	while (my @retrieved = $exc->fetchrow_array) {$noorfs{$retrieved[0]}++}
}




sub niceexit {
	my $reason = shift;
	$dbh->disconnect;
	if ($reason) {print STDERR $reason; exit(-1)}
	else {exit(0)}
	}
		

