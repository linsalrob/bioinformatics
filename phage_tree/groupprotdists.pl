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



# groupprotdists.pl

# based on combineprotdists.pl

use DBI;
use strict;
$|=1;

my $usage = "groupprotdists.pl <dir of prot dists> <matrix filename> <number of genomes used> <options>\nOPTIONS\n";
$usage .= "\t-n DON'T pad missing proteins with worst score (100)\n\t-m print out all protein matches\n";
$usage .= "\t-gr <filename> read groups to file\n\t-gf begin by grouping by family\n\t-rg # randomize genomes into # groups\n";
$usage .= "\t-ro randomize proteins into the same number of genomes as in the database, giving each group a similar # proteins\n";
$usage .= "\t-rp # randomize proteins into # groups\n\t-w <filename> write groups to file\n";



my $dbh=DBI->connect("DBI:mysql:phage", "apache") or die "Can't connect to database\n";
my $dir= shift || &niceexit($usage);
my $matrixfile= shift || &niceexit($usage);
my $nogenomes = shift || &niceexit($usage);

my $args= join (" ", @ARGV);
my $pad=1;
if ($args =~ /-n/) {$pad=0}



my %noorfs; my %genomename; my %genfam;
{
my ($noorfs, $genomename, $genfam)  = &getnoorfs;
%noorfs = %{$noorfs}; %genomename=%{$genomename}; %genfam=%{$genfam};
}

my %group; my %proteingroups; my $doprots;


if ($args =~ /-gr\s+(\S+)/) {
	my $filename = $1;
	open (GRP, $filename) || &niceexit("Can't open $filename\n");
	while (<GRP>) {
		next if (/^Org/);
		my ($genome, $group, @trash) = split /\t/;
		$group{$genome}=$group;
		}
	close GRP;
	}


if ($args =~ /-ro/) {
	my $group = &randomizebygenome();
	%proteingroups = %{$group};
	$doprots=1;
	}

elsif ($args =~ /-rp\s+(\d+)/) {
	my $groups=$1;
	my $group = &randomizeproteins($groups);
	%proteingroups = %{$group};
	$doprots=1;
	}
elsif ($args =~ /-rg\s+(\d+)/) {
	my $groupno = $1; my @genomes =  keys %genomename;
	my $groups = &randomizegenomes($groupno, \@genomes);
	%group = %$groups;
	}
else {
	unless ($args =~ /-gr/) { # only do this if we haven't read it in from a file.
	   if ($args =~ /-gf/) {
		   my $groups = &groupbyfamily();
		   %group = %$groups;
		   }
	    else {
		print "For each of the following genomes, please enter the number of the group to which it belongs\n";
		foreach my $genome (sort {$a <=> $b} keys %genomename) {
 		       print "$genome. $genomename{$genome} [".$genfam{$genome}."] : ";
		       my $group = <STDIN>;
		       chomp($group);
		       $group{$genome} = $group;
		       }
	   }
	}
	 my $yn; my $test=1;
	while ($test) {
		my %allgroups; my %genbyname;
		foreach my $key (keys %group) {
			$allgroups{$group{$key}}++;
			push (@{$genbyname{$group{$key}}}, $key);
			}
		print "SUMMARY OF DATA ENTERED\n";
		my @keys = sort {$a <=> $b} keys %allgroups;
		foreach my $key (@keys) {
 		  print  "\nGroup$key. $allgroups{$key} genomes:\n\t";
		  my @printgenomes = sort {uc($genomename{$a}) cmp uc($genomename{$b})} @{$genbyname{$key}};
		  foreach my $printgenome (@printgenomes) {print "\n\t$printgenome: $genomename{$printgenome} [".$genfam{$printgenome}."]"}
		}
		 print "\n\n\nIS THIS CORRECT? (y/n) ";
		 $yn = <STDIN>;
		 chomp($yn);
		 if (uc($yn) eq "N") {
		    undef $yn;
		    print "Genome number to change: ";
		    my $change = <STDIN>;
		    chomp($change);
		    print "$change. $genomename{$change}. Current group: $group{$change}. New group: ";
		    my $group = <STDIN>;
		    chomp($group);
		    $group{$change} = $group;
		    }
		  if (uc($yn) eq "Y") {undef $test}

	 }
	if ($args =~ /-w\s+(\S+)/) { #print out the data and correct noorfs simultaneously
		my $filename=$1; my %neworfs;
		open (GRP, ">$filename") || &niceexit("can't open $filename\n");
		print GRP "Org\tGroup\tName\tFamily\n";
		foreach my $key (keys %group) {
			print GRP "$key\t$group{$key}\t$genomename{$key}\t$genfam{$key}\n";
			}
		close GRP;
		}
	#correct noorfs
	my %neworfs;
	foreach my $key (keys %group) {$neworfs{$group{$key}} += $noorfs{$key}}
	%noorfs = %neworfs;
	
 }


# get all the genomes that make up a group

my %genomebygroup;
foreach my $genome (keys %group) {
	push (@{$genomebygroup{$group{$genome}}}, $genome)
	}



my @allscores;
my @matrix;
my @count;
my @proteinmatches; my @genematches;
my %linegenomecount; my %savematchesbyfile; my %savematchesbygroup;
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
	my %genegenome; my %genomesinfile;
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
		$genomesinfile{$genome}=1;
		
		unless ($gene && $genome) {&niceexit("Can't parse $line in $file\n")}
		
		if ($doprots) {$genome = $proteingroups{$gene}}
		else {
			$genome = $group{$genome};
			unless ($group{$genome}) {&niceexit("Can't get a group for genome |$genome|\n")}
		}
		
		push (@genes, $gene);
		push (@genomes, $genome);
		$checkdupgenomes{$genome}++;
		}


	# check to see whether all the genomes from any one group are in this matrix
	my %savematchesthistime;
	{
	foreach my $group (keys %genomebygroup) {
		my @genomesingroup = @{$genomebygroup{$group}};
		my $allthere = 1;
		foreach my $genomesingroup (@genomesingroup) {
			unless (exists $genomesinfile{$genomesingroup}) {undef $allthere}
			}
		if ($allthere) {$savematchesthistime{$group}=1}

		}
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
			&niceexit("PROBLEM WITH \n@line AND \n@genomes\nIN $file BECAUSE $#line AND $#genomes AT LINE $x\n");
			}
		my ($gene, $genome) = split (/_/, $line[0]);
		
		
		unless ($gene && $genome) {&niceexit("CAN'T PARSE @line SECOND TIME AROUND\n")}
		my $linegenome;
		if ($doprots) {$linegenome = $proteingroups{$gene}}
		else {$linegenome = $group{$genome}}

		#save all the matches if we want them. Note here $linegenome is now the group #
		if ($savematchesthistime{$linegenome}) {
			push (@{$savematchesbyfile{$file}{$linegenome}}, $gene);
			push (@{$savematchesbygroup{$linegenome}{$file}}, $gene);
			}

		$linegenomecount{$linegenome}++;
		my @seengenome;
		foreach my $x (1 .. $#genomes) {
#if ($x <= $z+1) {print STDERR " $line[$x]"; next}
			next if ($x <= $z+1);
			# If we are padding the table with 100s where there is no match, we
			# need to convert the -1's to 100. Otherwise we will ignore it.
			if ($line[$x] == -1) {if ($pad) {$line[$x] = 100} else {next}}
			
			#if it is itself, we want to make it zero. Otherwise, we'll save the protein numbers that match
			if ($genomes[$x] == $linegenome) {$line[$x] = '0.000'}
			else {
				my $genematch;
				# save the protein matches, but I only want to save them one way around
				# to make it easier
				if ($gene <$genes[$x]) {$genematch = $gene.",".$genes[$x].";".$line[$x]}
					else {$genematch = $genes[$x].",".$gene.";".$line[$x]}
				# protein match is a two dimensional array where each element is an array.
				# but it is called with an array! 4 dimensions?
				${$proteinmatches[$linegenome][$genomes[$x]]}{$genematch} =1;
				
				# save all the scores for later
				unless ($line[$x] == 100) {push (@{$allscores[$linegenome][$genomes[$x]]}, $line[$x])}
				# gene matches is all the genes from $linegenome that match genome. This will
				# be used to calculate the penalty for ORFs that are missed.
				${$genematches[$linegenome][$genomes[$x]]}{$gene} =1;
				
				}
			if ($matrix[$linegenome][$genomes[$x]]) {
				$matrix[$linegenome][$genomes[$x]] += $line[$x];
				$count[$linegenome][$genomes[$x]] ++;
				$seengenome[$linegenome][$genomes[$x]] ++;
				}
			else {
				$matrix[$linegenome][$genomes[$x]] = $line[$x];
				$count[$linegenome][$genomes[$x]] ++;
				$seengenome[$linegenome][$genomes[$x]] ++;
				}
			}
		# now we need to pad out all the missing genomes with 100's

		if ($pad) {
			foreach my $f (1 .. $nogenomes) {
				my $g;
				if ($doprots) {$g = $proteingroups{$gene}}
				else {$g = $group{$genome}}
				
				next if ($checkdupgenomes{$g});
				next if ($seengenome[$linegenome][$g]);
				if ($matrix[$linegenome][$g]) {
					$matrix[$linegenome][$g] += 100;
					$count[$linegenome][$g] ++;
					}
				else {
					$matrix[$linegenome][$g] = 100;
					$count[$linegenome][$g] ++;
					}
				}
			}
			
		}
	}


{

#now we need to penalize genomes that have only a few macthes.
# we will go through gene matches for each pair in the matrix, and
# add a penalty based on the number of missing orfs.

if ($pad) {
	open (MISS, ">missing.seqs.txt") || &niceexit("Can't open missing.seqs.txt\n");
	print MISS "Original\t#ORFs\tCompared to\t# similar\t# different\n";
	foreach my $y (0 .. $#genematches) {
	  next unless (exists $noorfs{$y}); # this just checks we have orfs for genome $y
	  foreach my $x (1 .. $#{$matrix[$y]}) {
	    next unless (exists $noorfs{$x});
	    next if ($y == $x);
	    my @similar = keys %{$genematches[$y][$x]};
	    my $difference = $noorfs{$y} - ($#similar+1);
	    print MISS "$y\t$noorfs{$y}\t$x\t",$#similar+1, "\t$difference\n";
	    next unless ($difference);
	    $matrix[$y][$x] += (100 * $difference);
	    $count[$y][$x] += $difference;
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
			  $matrix[$y][$x] = $matrix[$x][$y] = $matrix[$y][$x] + $matrix[$x][$y];
			  $count[$y][$x] = $count[$x][$y] = $count[$y][$x] + $count[$x][$y];
			  }
		}
		elsif ($matrix[$y][$x]) {
			unless ($matrix[$x][$y]) {
				$matrix[$x][$y] = $matrix[$y][$x];
				$count[$x][$y] = $count[$y][$x];
				}
			else {&niceexit("Can't figure out matrix $y and $x ($matrix[$x][$y] and $matrix[$y][$x] when merging 1\n")}
			}
		elsif ($matrix[$x][$y]) {
			unless ($matrix[$y][$x]) {
				$matrix[$y][$x] = $matrix[$x][$y];
				$count[$y][$x] = $count[$x][$y];
				}
			else {&niceexit("Can't figure out matrix $y and $x ($matrix[$x][$y] and $matrix[$y][$x] when merging 2\n")}
			}
		else {&niceexit("Can't figure out matrix $y and $x ($matrix[$x][$y] and $matrix[$y][$x] when merging 3\n")}
		
		# finally, take the average!
		
		$matrix[$x][$y] = $matrix[$y][$x] = $matrix[$y][$x]/$count[$y][$x];
	}
}

}

open (OUT, ">$matrixfile") || &niceexit("Can't open $matrixfile for writing\n");
# now we have all the data, lets just print out the matrix
print OUT $#matrix, "\n";
#foreach my $y (1 .. $#matrix) {print STDERR "\t$y"}
#print STDERR "\n";
foreach my $y (1 .. $#matrix) {
	my $tempstring = "group".$y;
	if (length($tempstring) > 10) {print STDERR "$tempstring is too long\n"}
	my $spacestoadd = " " x (10 - length($tempstring));
	print OUT $tempstring,$spacestoadd;
	foreach my $x (1 .. $#matrix) {
		if ($y == $x)  {print OUT "0 $noorfs{$x}  "; next}
		unless (defined $matrix[$y][$x]) {print OUT "100 0  "; next}
		unless ($matrix[$y][$x]) {
			print OUT "0 ";
			if ($count[$y][$x]) {print OUT "$count[$y][$x]  "}
			else {print OUT "0  "}
			next;
		}
		print OUT $matrix[$y][$x], " ", $count[$y][$x], "  ";
		}
	print OUT "\n";
	}
close OUT;

if ($args=~ /-m/) {
	open (PROT, ">$dir.protein.matches") || &niceexit("Can't open $dir.protein.matches for writing\n");
	
	#print out all the protein matches
	foreach my $y (1 .. $nogenomes) {
	my $tempstring = "group".$y;
	if (length($tempstring) > 10) {print STDERR "$tempstring is too long\n"}
	my $spacestoadd = " " x (10 - length($tempstring));
	print PROT $tempstring,$spacestoadd, "\t";
	foreach my $x (1 .. $nogenomes) {
		unless (defined $proteinmatches[$y][$x]) {print PROT "\t"; next}
		unless ($proteinmatches[$y][$x]) {print PROT "\t"; next}
		print PROT join (" ", sort keys %{$proteinmatches[$y][$x]}), "\t";
		}
	print PROT "\n";
	}
}

{
# now we will figure out all the scores. We want max, min, and average, to begin with.
  open (MIN, ">min.scores") || &niceexit("Can't open min.scores for writing\n");
  open (MAX, ">max.scores") || &niceexit("Can't open max.scores for writing\n");
  open (AV, ">average.scores") || &niceexit("Can't open average.scores for writing\n");
  open (ALL, ">all.scores") || &niceexit("Can't open all.scores for writing\n");
  
  foreach my $y (1 .. $#allscores) {
  print MIN "group$y\t"; print MAX "group$y\t"; print AV "group$y\t"; print ALL "group$y\t";
	foreach my $x (1 .. $#allscores) {
		unless ($allscores[$y][$x]) {print MIN "\t"; print MAX "\t"; print AV "\t"; print ALL "\t"; next}
		else {
			my @scores = sort {$a <=> $b} @{$allscores[$y][$x]};
			print MIN "$scores[0]\t"; print MAX "$scores[$#scores]\t"; 
			my $average; my $scorecount;
			foreach my $score (@scores) {$average += $score; $scorecount++} 
			$average = $average/$scorecount;
			print AV "$average\t"; print ALL join (' ', @scores), "\t";
		}
	}
  print MIN "\n"; print MAX "\n"; print AV "\n"; print ALL "\n";
  }
  close MIN; close MAX; close AV; close ALL;
}

{
# print out all the saved matches (if there ever are any)
open (MATCHES, ">saved.matches.by.file") || &niceexit("Can't open saved.matches.by.file for writing\n");
foreach my $file (sort {$a cmp $b} keys %savematchesbyfile) {
	print MATCHES "FILE: ", $file, "\n";
  foreach my $group (sort {$a <=> $b} keys %{$savematchesbyfile{$file}}) {
	print MATCHES $group, ": ", join (" ",sort {$a <=> $b}  @{$savematchesbyfile{$file}{$group}}), "\n";
	}
  }
 close MATCHES;
 #print matches by group
 open (MATCHES, ">saved.matches.by.group") || &niceexit("Can't open saved.matches.by.group for writing\n");
foreach my $group (sort {$a <=> $b} keys %savematchesbygroup) {
	print MATCHES "GROUP: ", $group, "\n";
  foreach my $file (sort {$a cmp $b} keys %{$savematchesbygroup{$group}}) {
	print MATCHES $file, ": ", join (" ",sort {$a <=> $b} @{$savematchesbygroup{$group}{$file}}), "\n";
	}
  }
 close MATCHES;

}
&niceexit(0);

sub getnoorfs {
	my %noorfs; my %genomename; my %genfam;
	my $exc = $dbh->prepare("select organism from protein");
	$exc->execute or die $dbh->errstr;
	while (my @retrieved = $exc->fetchrow_array) {$noorfs{$retrieved[0]}++}

	$exc = $dbh->prepare("select count, organism, family from phage");
	$exc->execute or die $dbh->errstr;
	while (my @ret = $exc->fetchrow_array) {
		$genomename{$ret[0]}=$ret[1];
		$genfam{$ret[0]}=$ret[2];
		}
	return \%noorfs, \%genomename, \%genfam;
}




sub niceexit {
	my $reason = shift;
	$dbh->disconnect;
	if ($reason) {print STDERR $reason; exit(-1)}
	else {exit(0)}
	}


sub randomizegenomes {
	# randomly assign genomes to groups
	my $groups = shift;
	my $genomes = shift;
	for (my $i = @$genomes; --$i; ) {
		my $j = int rand ($i+1);
		next if ($i == $j);
		@$genomes[$i, $j] = @$genomes[$j, $i];
		}
	my $split = int(($#{$genomes}+1)/$groups);
	my $gencount=1; my %group; my %gencount; my %groupsbygen;
	while ($#{$genomes} > $split) {
		foreach my $x (1 .. $split) {
			my $gen = shift (@$genomes);
			$group{$gen} = $gencount;
			$gencount{$gencount}++;
			push (@{$groupsbygen{$gencount}}, $gen);
			}
		$gencount++;
		}
	foreach my $gen (@{$genomes}) {$group{$gen} = $gencount; $gencount{$gencount}++; push (@{$groupsbygen{$gencount}}, $gen);}
	foreach my $x (1 .. $gencount) {print "Group $x : $gencount{$x} genomes: ", join (" ", @{$groupsbygen{$x}}), "\n"}

	return \%group;

}

sub randomizeproteins {
	# randomly assign proteins to groups
	my $groups = shift;
	print "MAKING $groups GROUPS\n";
	my @proteins;
	my $exc = $dbh->prepare("select count from protein");
	$exc->execute or die $dbh->errstr;
	while (my @retrieved = $exc->fetchrow_array) {push (@proteins, $retrieved[0])}
	
	print $#proteins+1, " proteins in total\n";
	
	for (my $i = @proteins; --$i; ) {
		my $j = int rand ($i+1);
		next if ($i == $j);
		@proteins[$i, $j] = @proteins[$j, $i];
		}
	
	
	my $split = int(($#proteins+1)/$groups);
	print "Putting $split proteins into each group\n";
	my $count; my %group; my %count;
	for ($count=1; $count <= $groups; $count++) {
		foreach my $x (1 .. $split) {
			my $prot = shift (@proteins);
			$group{$prot} = $count;
			$count{$count}++;
			}
		}
	$count--;
	foreach my $prot (@proteins) {$group{$prot} = $count; $count{$count}++}
	foreach my $x (1 .. $count) {print "Group $x : $count{$x} proteins.\n"}

	return \%group;

}
	


sub randomizebygenome {
	# randomly assign proteins to groups
	# but make each group have the same number of proteins as the genome with the same number.

	my $groups; # number of groups to make. Get this from database
	my $exc = $dbh->prepare("select count from phage");
	$exc->execute or die $dbh->errstr;
	while (my @retrieved = $exc->fetchrow_array) {$groups++}

	print "MAKING $groups GROUPS\n";
	my @proteins; my %protspergroup;
	$exc = $dbh->prepare("select count, organism from protein");
	$exc->execute or die $dbh->errstr;
	while (my @retrieved = $exc->fetchrow_array) {
		push (@proteins, $retrieved[0]);
		$protspergroup{$retrieved[1]}++;
		}
	
	print $#proteins+1, " proteins in total\n";
	
	for (my $i = @proteins; --$i; ) {
		my $j = int rand ($i+1);
		next if ($i == $j);
		@proteins[$i, $j] = @proteins[$j, $i];
		}
	
	my %group; my %count; my $count; my $max=1;
	foreach $count (keys %protspergroup) {
		foreach my $x (1 .. $protspergroup{$count}) {
			my $prot = shift (@proteins);
			$group{$prot} = $count;
			$count{$count}++;
			}
		if ($max < $count) {$max=$count};
		}

	if (@proteins) {print STDERR $#proteins," elements left in \@proteins after sorting??\n"}

	return \%group;

}








	
sub groupbyfamily {
	my %group;
	my $exc = $dbh->prepare("select count, family from phage");
	$exc->execute or die $dbh->errstr;
	my %family; my %genbyfam;
	while (my @ret = $exc->fetchrow_array) {$genbyfam{$ret[0]}=$ret[1]; $family{$ret[1]}=1}
	my @keys = sort {uc($a) cmp uc($b)} keys %family;
	my %familygroup;
	open FAM, ">family.groups" || &niceexit("Can't open family groups\n");
	foreach my $x (0 .. $#keys) {$familygroup{$keys[$x]} = $x+1; print FAM $x+1,"\t$keys[$x]\n"}
	close FAM;
	foreach my $genome (keys %genbyfam) {$group{$genome} = $familygroup{$genbyfam{$genome}}}
	return \%group;
	}
	
	
