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

# modified to work with bacterial assemblies!!

# new version. This will assign a distance of $penalty to any sequence that does not match, and to all empty
# spaces. There is a good rationale for this. The The Dayhoff PAM matrix scoring system returns a percent of
# the amino acids  that are likely to have changed. Therefore a 100% score means that they have all changed.

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

my $filecount; my $oldtime=$^T;
my @matrix;
my @count;
my @proteinmatches; my @genematches;
my %linegenomecount;
# read each file one at a time, and add the data to an array
opendir(DIR, $dir) || &niceexit("Can't open $dir");
while (my $file=readdir(DIR)) {
	next if ($file =~ /^\./);
	open (IN, "$dir/$file") || &niceexit("Can't open $dir/$file");
	$filecount++; unless ($filecount % 1000) {print STDERR "$filecount\t", time-$^T, "\t", time-$oldtime, "\n"; $oldtime = time}
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
	# then we add each value to the pre-existing value in the matrix
	# note that because the genomes are represented as numbers we can just
	# use these numbers for the position in the matrix.
	# we are going to also count the number of times that we save each data
	# point for the final average.
	# Finally, we only do this in one direction because the input
	# matrices are complete (and identical) on both halves.
	
	# note that column zero of the matrix is empty (there is no genome 0)
	foreach my $z (0 .. $#dists) {
#print STDERR "\nSKIPPING: ";
		my @line = split (/\s+/, $dists[$z]);
		unless ($#line == $#genomes) {
			my $x; foreach my $y (0 .. $#dists) {if ($dists[$y] eq $dists[$z]) {$x = $y}}
			&niceexit("PROBLEM WITH \n@line AND \n@genomes\n\nIN FILE: $file\n\nBECAUSE $#line AND $#genomes AT LINE $x\n");
			}
		my ($gene, $linegenome) = split (/_/, $line[0]);
		unless ($gene && $linegenome) {&niceexit("CAN'T PARSE @line SECOND TIME AROUND\n")}
		$linegenomecount{$linegenome}++;
		my @seengenome;
		foreach my $x (1 .. $#genomes) {
#if ($x <= $z+1) {print STDERR " $line[$x]"; next}
			next if ($x <= $z+1);
			# If we are padding the table with 100s where there is no match, we
			# need to convert the -1's to 100. Otherwise we will ignore it.
			if ($line[$x] == -1) {if ($pad) {$line[$x] = $penalty} else {next}}
			
			#if it is itself, we want to make it zero. Otherwise, we'll save the protein numbers that match
			# note that we can store all the pairwise protein matches, but we need the gene matches
			# so that we can pad out missing genes correctly
			if ($genomes[$x] == $linegenome) {$line[$x] = '0.000'}
#			else {
#				my $genematch;
				# save the protein matches, but I only want to save them one way around
				# to make it easier
#				if ($gene <$genes[$x]) {$genematch = $gene.",".$genes[$x].";".$line[$x]}
#					else {$genematch = $genes[$x].",".$gene.";".$line[$x]}
				# protein match is a two dimensional array where each element is an array.
				# but it is called with an array! 4 dimensions?
#				${$proteinmatches[$linegenome][$genomes[$x]]}{$genematch} =1;
				# gene matches is all the genes from $linegenome that match genome. This will
#				# be used to calculate the penalty for ORFs that are missed.
#				}
			$matrix[$linegenome][$genomes[$x]] += $line[$x];
			$count[$linegenome][$genomes[$x]] ++;
			$seengenome[$linegenome][$genomes[$x]] ++;
			${$genematches[$linegenome][$genomes[$x]]}{$gene} =1;
			}
		# now we need to pad out all the missing genomes with 100's
		if ($pad) {
			foreach my $x (1 .. $nogenomes) {
				next if ($checkdupgenomes{$x});
				next if ($seengenome[$linegenome][$x]);
				$matrix[$linegenome][$x] += $penalty;
				$count[$linegenome][$x] ++;
				}
			}
			
		}
	}


{

#now we need to penalize genomes that have only a few macthes.
# we will go through gene matches for each pair in the matrix, and
# add a penalty based on the number of missing orfs.
# note, we are adding this above, each time around

unless (1) {
if ($pad) {
	open (MISS, ">missing.seqs.txt") || &niceexit("Can't open missing.seqs.txt\n");
	print MISS "Original\t#ORFs\tCompared to\t#ORFS\t# similar\t# different\tcurr. score\tcurr. count\tpenalty\n";
	foreach my $y (0 .. $#genematches) {
	  next unless (exists $noorfs{$y}); # this just checks we have orfs for genome $y
	  foreach my $x (1 .. $#{$matrix[$y]}) {
	    next unless (exists $noorfs{$x});
	    next if ($y == $x);
	    my @similar = keys %{$genematches[$y][$x]};
	    my $difference = $noorfs{$y} - ($#similar+1);
	    print MISS "$y\t$noorfs{$y}\t$x\t$noorfs{$x}\t",$#similar+1, "\t$difference\t$matrix[$y][$x]\t$count[$y][$x]\t",($penalty * $difference),"\n";
	    next unless ($difference);
	    $matrix[$y][$x] += ($penalty * $difference);
	    $count[$y][$x] += $difference;
	    }
	  }
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
		
		# because we are only looking at one half of the matrix (see above)
		# we need to be sure that both halves are the same.
		# this loop will take care of that.
		
		if ($matrix[$y][$x] && $matrix[$x][$y]) {
			unless ($matrix[$y][$x] == $matrix[$x][$y]) {
#print STDERR "$matrix[$y][$x] and $matrix[$x][$y] (counts:  $count[$y][$x] AND $count[$x][$y]) ";
			  $matrix[$y][$x] = $matrix[$x][$y] = $matrix[$y][$x] + $matrix[$x][$y];
			  $count[$y][$x] = $count[$x][$y] = $count[$y][$x] + $count[$x][$y];
print STDERR " $matrix[$y][$x] and $matrix[$x][$y] (counts:  $count[$y][$x] AND $count[$x][$y]) Fixed at 1\n";
			  }
		}
		elsif ($matrix[$y][$x]) {
			unless ($matrix[$x][$y]) {
				$matrix[$x][$y] = $matrix[$y][$x];
				$count[$x][$y] = $count[$y][$x];
#print STDERR "Fixed at 2  ";
				}
			else {&niceexit("Can't figure out matrix $y and $x ($matrix[$x][$y] and $matrix[$y][$x] when merging 1\n")}
			}
		elsif ($matrix[$x][$y]) {
			unless ($matrix[$y][$x]) {
				$matrix[$y][$x] = $matrix[$x][$y];
				$count[$y][$x] = $count[$x][$y];
#print STDERR "Fixed at 3  ";
				}
			else {&niceexit("Can't figure out matrix $y and $x ($matrix[$x][$y] and $matrix[$y][$x] when merging 2\n")}
			}
		else {&niceexit("Can't figure out matrix $y and $x ($matrix[$x][$y] and $matrix[$y][$x] when merging 3\n")}
		
		# finally, take the average!
		
		$matrix[$x][$y] = $matrix[$y][$x] = $matrix[$y][$x]/$count[$y][$x];
#print STDERR "AVERAGES : $matrix[$x][$y] AND $matrix[$y][$x]\n";
	}
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
			

if ($args=~ /-m/) {
	open (PROT, ">$dir.protein.matches") || &niceexit("Can't open $dir.protein.matches for writing\n");
	
	#print out all the protein matches
	foreach my $y (1 .. $nogenomes) {
	my $tempstring = "genome".$y;
	if (length($tempstring) > 10) {print STDERR "$tempstring is too long\n"}
	my $spacestoadd = " " x (10 - length($tempstring));
	print PROT $tempstring,$spacestoadd, "\t";
	foreach my $x (1 .. $nogenomes) {
		unless (defined $proteinmatches[$y][$x]) {print PROT "\t"; next}
		unless ($proteinmatches[$y][$x]) {print PROT "\t"; next}
		my @allmatches = (keys %{$proteinmatches[$y][$x]}, keys %{$proteinmatches[$x][$y]});
		my %allmatches;
		@allmatches{@allmatches}=1;
		@allmatches = sort keys %allmatches;
		print PROT join (" ", sort @allmatches), "\t";
		}
	print PROT "\n";
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
