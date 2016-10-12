#!/usr/bin/perl -w

=pod

This is scaffold_builder.pl. 

Written by Rob Edwards, 2011

scaffold_builder.pl takes contigs and assembles them into scaffolds based on their blast hits. The overall approach that 
scaffold_builder uses is to compare contigs to a reference genome using BLASTN, and build a simple scaffold. The scaffold is walked from 
beginning to end and gaps between contigs filled with an appropriate number of N's  


1. Concatenate the reference genome using 5000 x N between contigs (or something similar)
2. Use blastn to compare the assembled genome to the reference genome. I usually use a stringent cutoff here since we only parse out the best hit
3. Combine the contigs into scaffolds with scaffold_builder.pl


Version: $Id: scaffold_builder.pl,v 1.8 2011/08/24 15:41:58 linsalrob Exp $

=cut

use strict;
use Getopt::Std;
my %opts;
getopts('Da:b:c:d:f:g:hil:n:o:p:q:r:s:vx:z', \%opts);
die <<EOF if (!($opts{o} && $opts{q} && $opts{d}) || $opts{h});
$0
-q <file>	query fasta file
-d <file>	database fasta file
-b <file>	blast output file (if available). If not provided, the contigs in the database file will be concatenated and used for the blast
-o <file>	output file base name. We will make several output files with this name
-p <n>		minimum percent alignment length. To be considered a useful match, the alignment length (reported by blast) must be at least this % of the contig length 
-l <n>		minimum length of the contig to include. This allows another place to filter, eg, for longer contigs

Overlap alignment settings
-a <directory>	alignment directory. Default is alignment_dir
-x <n>		extra sequence to include in the alignment (bp on either side of the overlap region). Default = 5bp
-n <n>		maximum overlap to try and align. If there is too much overlap, it suggests either an assembly error or large duplicate regions (e.g. rrns). Limit the overlap.
-f <n>		if the overlap is this length or less, just force the alignment regardless of the percent similarity

Combine overlaps and insert gaps into the alignment
-g <n>		Insert gaps in the alignment provided the clustal alignment is > this percent identical. The percent does not include "-"'s inserted at the 5' or 3' end of the sequences
-c <n>		number of N's to insert between contigs that do not have an alignable overlap (default=0 i.e. join contigs together end-to-end)
-i 		include contigs in the final fasta file that were not included in the assembly (i.e. they are too short and/or they are not similar to the reference)

Other options
-D 		use some sensible default settings (currently: -a scaffold_alignments -p 60 -l 200 -n 100 -g 80 -x 10 -c 1000 -i -f 20. Other options will overwrite these if they are set).


Output
-s <n>		split the contigs back out following this many N's
-z 		after splitting the contigs, order by size. If size is not given (with -s) the default of 5000 is used. 
-v 		verbose output on STDERR. You should almost certainly redirect this to a file
-r 		directory to put all the output into (default = scaffold_builder)

EOF






#################################################################################
#                                                                               #
# DEFINE SOME LOCAL VARIABLES AND EXECUTABLES.                                  #
#                                                                               #
# You should probably change these things, but hopefully the default            #
# settings will work OK                                                         #
#                                                                               #
#################################################################################

my $blastExecutable   = "/usr/bin/blastall";
my $formatdbExecutable = "/usr/bin/formatdb";
my $clustalExecutable = "/usr/bin/clustalw";


#################################################################################
#                                                                               #
# This ends the local variables! Everything below should work as-is             #
#                                                                               #
#################################################################################

#################################################################################
#                                                                               #
# Set the default settings if we desire.                                        #
#                                                                               #
#################################################################################

if ($opts{D}) {
	$opts{p} = $opts{p} || 60;
	$opts{l} = $opts{l} || 200;
	$opts{n} = $opts{n} || 100;
	$opts{g} = $opts{g} || 80;
	$opts{x} = $opts{x} || 10;
	$opts{c} = $opts{c} || 1000;
	$opts{i} = 1;
	$opts{f} = $opts{f} || 20;
	$opts{o} = $opts{o} || "built_scaffolds";

	print STDERR "Running with these options:\n", (map {join("\n", " " . $_." : ".$opts{$_})} keys %opts), "\n\n";
}


# define the genetic code. Much of this is work by Gary Olsen, University of Illinois, Urbana Champaign. Thanks, Gary.
my %DNA_letter_can_be = (
		A => ["A"],                 a => ["a"],
		B => ["C", "G", "T"],       b => ["c", "g", "t"],
		C => ["C"],                 c => ["c"],
		D => ["A", "G", "T"],       d => ["a", "g", "t"],
		G => ["G"],                 g => ["g"],
		H => ["A", "C", "T"],       h => ["a", "c", "t"],
		K => ["G", "T"],            k => ["g", "t"],
		M => ["A", "C"],            m => ["a", "c"],
		N => ["A", "C", "G", "T"],  n => ["a", "c", "g", "t"],
		R => ["A", "G"],            r => ["a", "g"],
		S => ["C", "G"],            s => ["c", "g"],
		T => ["T"],                 t => ["t"],
		V => ["A", "C", "G"],       v => ["a", "c", "g"],
		W => ["A", "T"],            w => ["a", "t"],
		Y => ["C", "T"],            y => ["c", "t"]
		);
my %iupac;
map {$iupac{join("", @{$DNA_letter_can_be{$_}})}=$_} keys %DNA_letter_can_be;


unless (defined $opts{o}) {$opts{o}="-"}
unless (defined $opts{c}) {$opts{c}=0}
unless (defined $opts{a}) {$opts{a}="alignment_dir"}
unless (defined $opts{r}) {$opts{r}="scaffold_builder"}
unless (-d $opts{r}) {mkdir $opts{r}, 0755}


my %ignore;
my $cwd = `pwd`; chomp($cwd); # current working directory, in case we chdir around
if (index($opts{d}, "/") != 0) {$opts{d} = $cwd."/".$opts{d}}
if (index($opts{q}, "/") != 0) {$opts{q} = $cwd."/".$opts{q}}
if (defined $opts{b} && index($opts{b}, "/") != 0) {$opts{b} = $cwd."/".$opts{b}}
chdir($opts{r});
$cwd = `pwd`; chomp($cwd); # current working directory, in case we chdir around

# do we need to run the blast
unless (defined $opts{b}) {
	$opts{q} =~ /([\w\.\-]+)$/;
	$opts{b} = "$cwd/$1";
	$opts{d} =~ /([\w\.\-]+)$/;
	$opts{b} .= ".blastn.".$1;

	mkdir "blast", 0755;
	chdir("blast");

	my $dbs = &read_fasta($opts{d});
	open(BLASTDB, ">database_$$.fasta") || die "can't open blast/database_$$.fasta";
	print BLASTDB ">contig1\n";
	print BLASTDB join("N" x 5000, sort {length($b) <=> length($a)} values %$dbs), "\n";
	close BLASTDB;
	`$formatdbExecutable -i database_$$.fasta -p F`;
	chdir($cwd);
	`$blastExecutable -i $opts{q} -d blast/database_$$.fasta -p blastn -o $opts{b} -m 8 -e 1e-50`;
}

my $query; my $dblength;
{ # limit this scope of this reading


	my $queryInput = &read_fasta($opts{q});
	map {my $i = $_; s/\s+.*$//; $query->{$_}=$queryInput->{$i}} keys %$queryInput;
	if ($opts{l}) {
		map {$ignore{$_} = "contig too short (only ".length($query->{$_})." bp)" if (length($query->{$_}) < $opts{l})} keys %$query;
	}

	my $db = &read_fasta($opts{d});
	map {$dblength += length($db->{$_}); my $i = $_; s/\s+.*$//; $db->{$_}=$db->{$i}} keys %$db;
}


# convert the percent to a fraction so we can use it directly later
if ($opts{p} && $opts{p} > 1) {$opts{p} /= 100}


#################################################################################
#                                                                               #
# Figure out how big an array we need to use                                    #
#                                                                               #
#################################################################################


# We will use an array to hold all of the sequences (actually an array of hashes). 
# We are going to make a larger array, to allow for an offset if we have sequences that start before our database sequence does
# We also need to know how big the array should be, so we will scan through the blast results and calculate the dimensions of the array
# before actually populating it

my $max=0; my $min=1e6;
{ # scope this block
	open(IN, $opts{b}) || die "Can't open $opts{b}";
	my %seen;
	while (<IN>) {
		chomp;
		my ($q, $d, $p, $l, $m, $g, $qs, $qe, $ds, $de, $e, $b)=split /\t/;
		next if ($seen{$q});
		$seen{$q}=1;

		if ($opts{p} && $l/length($query->{$q}) < $opts{p}) {$ignore{$q} = "Alignment was only ". (100 * $l/length($query->{$q})). "% of the sequence length"}

		next if ($ignore{$q});


# where is the start of the contig
		($ds > $de) ? (($ds, $de)=($de, $ds)) : 1;
		($de > $max) ? ($max=$de) : 1;
# adjust the start for the offset at the beginning of the query
		$ds -= $qs;
		($ds < $min) ? ($min = $ds) : 1;

	}
	close IN;
}

# now we need an array from $min to $max, but $min can be <=0
# we also want to use a 1-based index so that things make sense to the biologists.
my $size = $max - $min;
my $sequence;
$#$sequence = $size + 1;

###########################################################################
#                                                                         #
# Read the blast file and store all the sequences in the array            #
#                                                                         #
###########################################################################

my %blastAlignmentResults; # the length of each contig that aligns.
{ # scope this block
	open(IN, $opts{b}) || die "Can't open $opts{b}";
	my %seen;
	while (<IN>) {
		chomp;
		my ($q, $d, $p, $l, $m, $g, $qs, $qe, $ds, $de, $e, $b)=split /\t/;
		next if ($seen{$q});
		$seen{$q}=1;
# next if ($ignore{$q});
		if ($ignore{$q}) {
			if ($opts{v}) {print STDERR "Contig $q ignored because $ignore{$q}\n"}
			next;
		}

		# where is the start of the alignment
		my $mylocation = $ds - $qs;

		# ends are not used, but I put them here as a reminder
# my $myend   = $de + (length($query->{$q}) - $qe); ## we don't use the real end -- we use the length of the query. But it should be the same
# $myend = $ds + $qs; ## this if if we reverse complement the sequence 

# reverse complement the query sequence if we need to
		if (($ds > $de && $qs < $qe) || ($ds < $de && $qs > $qe)) {
			$query->{$q} = &rc($query->{$q});
			$mylocation = $de - (length($query->{$q}) - $qe);
		}

# adjust the start for the offset at the beginning of the query, including any -ve numbers

		$mylocation -= $min; # adjust the start for the beginning of the array
		$mylocation++; # now make it a 1-based index: This is just to make the biologists happy!!

		my $posn=0;
		while ($posn < length($query->{$q})) {
			if (defined $sequence->[$mylocation]->{$q}) {
				die "we are attempting to overwrite a sequence. This should not happen";
			}
			$sequence->[$mylocation]->{$q} = [substr($query->{$q}, $posn, 1), $posn+1];
			$mylocation++;
			$posn++;
		}

		$blastAlignmentResults{$q} = [$p, $l];
	}
}


###########################################################################
#                                                                         #
# For each sequence, identify the sequence and the positions              #
#                                                                         #
###########################################################################

my $indices=[undef]; # an array of all contigs we are currently displaying
my %index; # a hash of where each contig will appear
my $locations; # a hash of contig as key and [start, stop] as value

# What we are going to print:
# This is an array of arrays. 
#	The first array is the sequence.
#	The second is the poistion in the original sequence
# 	The third is the contig name,
my $output; 
my $range; # what is the range of any given contig
{ # scope this block
	my %multipleContigsWarning;

	for (my $i=1; $i<=$max; $i++) {
#  Are we at the end of the previous contig?
		$indices = &check_ends($i, $indices);

# is this position even defined?
		next if (!defined $sequence->[$i]);


		foreach my $contig (@$indices) {
			if (defined $contig && defined $sequence->[$i]->{$contig}) {
				unless (defined $range->{$contig}->[0]) {$range->{$contig}->[0]=$i}
				$range->{$contig}->[1]=$i;
				push @{$output->[$i]->[0]}, $sequence->[$i]->{$contig}->[0];
				push @{$output->[$i]->[1]}, $sequence->[$i]->{$contig}->[1];
				push @{$output->[$i]->[2]}, $contig;
			} 
		}
		if ($output->[$i] && scalar(@{$output->[$i]->[2]}) > 2) {
			my $string = join("", sort {uc($a) cmp uc($b)} @{$output->[$i]->[2]});
			if (!$multipleContigsWarning{$string}) {
				print STDERR "WARNING: We have ", scalar(@{$output->[$i]->[2]}), " contigs  (",
				      join(" ", @{$output->[$i]->[2]}), ") starting at $i, results will be unpredictable in this region!\n";
				map {printf STDERR ("\t$_ alignment length: %d (%.3f%% of the query length) with %d%% similarity\n", 
				$blastAlignmentResults{$_}->[1], ($blastAlignmentResults{$_}->[1]/length($query->{$_}))*100, $blastAlignmentResults{$_}->[0])} @{$output->[$i]->[2]};
			}
			$multipleContigsWarning{$string}=1;
		}

	}

}

# make sure we end any contigs that are at the end of the sequence
foreach my $contig (@$indices) {
	next unless (defined $contig);
	$locations->{$contig}->[1] = $max;
}

# sort the contigs in the order that they start
my @contigs = sort {$locations->{$a}->[0] <=> $locations->{$b}->[0]} keys %$locations;

###########################################################################
#                                                                         #
# Align the overlaps if we need to                                        #
#                                                                         #
###########################################################################

## create the alignments of overlapping regions
my $clustalalignment; # a hash of the sequences that we align.
{ # scope this piece of work
	unless (defined $opts{x}) {$opts{x} = 5} # set the default extra DNA to extract
	unless (defined $opts{n}) {$opts{n} = 1e20}
	unless (-d "$opts{a}") {mkdir "$opts{a}", 0755}
	open(CLUSTAL, ">$opts{a}/clustal_progress.txt") || die "can't write to $opts{a}/clustal_progress.txt";
	my $overlap;
	for (my $i=1; $i<=$#contigs; $i++) {
		my $bp = $locations->{$contigs[$i-1]}->[1] - $locations->{$contigs[$i]}->[0];
		if ($bp > 0 && $bp <= $opts{n}) {
		# there is an overlap with the contig before. We need to get that many
		# bases from the beginning of this contig, and from the end of the previous
		# contig and align them
			$bp += $opts{x};
			if ($opts{v}) {print STDERR "Getting overlap between $contigs[$i] (",$locations->{$contigs[$i]}->[0],") and ", $contigs[$i-1], " (",$locations->{$contigs[$i-1]}->[1],") for a total of $bp base pairs\n"}
			$overlap++;
			my $outputfile = "overlap".$overlap."_".$contigs[$i]."_".$contigs[$i-1];
			open(FA, ">$opts{a}/$outputfile") || die "Can't open $opts{a}/$outputfile";
			print FA ">$contigs[$i-1] -$bp\n", substr($query->{$contigs[$i-1]}, -$bp), "\n";
			print FA ">$contigs[$i] 0 to $bp\n", substr($query->{$contigs[$i]}, 0, $bp), "\n";
			#$clustal->{$contigs[$i-1]}->{$contigs[$i]} = [length($query->{$contigs[$i-1]})-$bp, 0, $bp, "$opts{a}/$outputfile"];
			$clustalalignment->{$contigs[$i-1]}->{length($query->{$contigs[$i-1]})-$bp}="$opts{a}/$outputfile.aln";
			$clustalalignment->{$contigs[$i]}->{0}="$opts{a}/$outputfile.aln";


			close FA;
			if ($opts{v}) {print STDERR "Aligning overlap $overlap\n"}
			print CLUSTAL `$clustalExecutable $opts{a}/$outputfile 2>&1`;
		}
	}
	close CLUSTAL;
}


###########################################################################
#                                                                         #
# Insert gaps into our sequences based on the alignments                  #
#                                                                         #
###########################################################################
my $revOutput=[];
my $consensus=[];
$#$revOutput = $#$output;
$#$consensus = $#$output;
my $consensusPosition=1;
my $currentConsensusContig; # the contig we are reading for the consensus
my $nextConsensusContig; # the next one to be appended at the end of the current sequence
if ($opts{g}) {
	if ($opts{g} < 1) {$opts{g} *= 100} # make sure we are working with percents not fractions

# march down the list, and bubble as needed
	my $newx=1;
	my %newposition; # the new position adjusted for gaps
	my %lastposition; # the last position that we looked at
	for (my $curr=1; $curr <= $max; $curr++) {
		my $bases=[]; my $posns=[]; my $conts=[]; my %index;
		if (defined $output->[$curr]) {
			$bases = $output->[$curr]->[0];
			$posns = $output->[$curr]->[1];
			$conts = $output->[$curr]->[2];
			%index = map {($conts->[$_]=>$_)} (0 .. $#$conts);
		}
		
		my $skip; # if we have already written these bases (e.g. from an alignment) we want to skip out of them
		foreach my $contig (@$conts) {
			$skip =1 if (defined $lastposition{$contig} && $lastposition{$contig} >= $posns->[$index{$contig}]);
		}
		next if ($skip);

		# do we still have sequence for the currentConsensusContig
		if ($currentConsensusContig && !defined $index{$currentConsensusContig}) {
			if ($nextConsensusContig) {
				for (my $insertN=0; $insertN<=$opts{c}; $insertN++) {
					$consensus->[$consensusPosition++]->[0] = "n";
				}
				while (@{$nextConsensusContig->[0]}) {
					$consensus->[$consensusPosition]->[0] = shift @{$nextConsensusContig->[0]};
					$consensus->[$consensusPosition]->[1] = shift @{$nextConsensusContig->[1]};
					$consensusPosition++;
				}
			}
			undef $currentConsensusContig;
			undef $nextConsensusContig;
		}


# if we have an alignment here, read the alignment and insert it
		my $alignmentFile=undef;

		for (my $c=0; $c<=$#$conts; $c++) {
			if ($clustalalignment->{$conts->[$c]}->{$posns->[$c]}) {
				$alignmentFile = $clustalalignment->{$conts->[$c]}->{$posns->[$c]};
			}
		}
# do we actually want to correct this alignment, or not?

		my $clustal; # the clustal alignment object for this file
		if ($alignmentFile) {
			$clustal = &parse_file($alignmentFile);
			if (&alignment_length($clustal) < $opts{f}) {
				print STDERR "Aligning $alignmentFile anyway because trimmed alignment length was ", &alignment_length($clustal), "\n";
			}
			elsif (&trimmed_percent_identical($clustal) < $opts{g}) {
				if ($opts{v}) {printf STDERR "Did not correct %s because the similarity is too low - only %.3f%% identical\n", $alignmentFile, &trimmed_percent_identical($clustal)}
				undef $alignmentFile;
			}
		}

		###########
		# 
		# this is the part where we read the clustal file and add that into the sequence
		#
		###########

		if ($alignmentFile) {
			printf STDERR "Merging contigs : %s because they are %.3f%% identical\n", join(" ", @{$clustal->{'ids'}}), &trimmed_percent_identical($clustal);

			my $seqbases;
			foreach my $id (@{$clustal->{'ids'}}) {
				my $sequence = $clustal->{'sequences'}->{$id};
				@{$seqbases->{$id}} = split //, $sequence;
			}

			my $length = &alignment_length($clustal);
			for (my $alignmentposn = 0; $alignmentposn < $length; $alignmentposn++) {
				foreach my $contig (keys %$seqbases) {
					unless (defined $index{$contig}) {
						#push @$conts, $contig;
						$index{$contig}=$#$conts+1;
						if (!defined $newposition{$contig}) {$newposition{$contig}=1}
						if (!defined $lastposition{$contig}) {$lastposition{$contig}=0}
					}
					my $newbase = shift(@{$seqbases->{$contig}});
					$revOutput->[$newx]->[0]->[$index{$contig}] = $newbase;
					$revOutput->[$newx]->[1]->[$index{$contig}] = $newposition{$contig}++;
					$revOutput->[$newx]->[2]->[$index{$contig}] = $contig;
					$revOutput->[$newx]->[3] = "Merged by clustal";
					if ($newbase ne "-") {$lastposition{$contig}++}
					
					# for the consensus location
					my $locationString = "$contig:".$revOutput->[$newx]->[1]->[$index{$contig}]." ";
					$consensus->[$consensusPosition]->[1] .= $locationString;
				}
				$consensus->[$consensusPosition]->[0] = &iupac($revOutput->[$newx]->[0]);
				$newx++;
				$consensusPosition++;
			}
		}
		elsif (scalar(@$conts)) {
			
			################
			#
			# This is for those alignments that are too dissimilar to merge ... or there is one contig!
			#
			###############
			
			
			foreach my $contig (@$conts) {
				next unless (defined $contig);
				if (!defined $lastposition{$contig}) {$lastposition{$contig} = 0}
				if (!defined $newposition{$contig}) {$newposition{$contig} = $output->[$curr]->[1]->[$index{$contig}]}
				next if ($posns->[$index{$contig}] <= $lastposition{$contig});
				$lastposition{$contig}=$posns->[$index{$contig}];
				$revOutput->[$newx]->[0]->[$index{$contig}] = $output->[$curr]->[0]->[$index{$contig}];
				$revOutput->[$newx]->[1]->[$index{$contig}] = $newposition{$contig}++;
				$revOutput->[$newx]->[2]->[$index{$contig}] = $contig;
				$revOutput->[$newx]->[3] = "Original base";
				if (!defined $currentConsensusContig) {$currentConsensusContig=$contig}
				if ($currentConsensusContig eq $contig) {
					$consensus->[$consensusPosition]->[0] = $revOutput->[$newx]->[0]->[$index{$contig}];
					$consensus->[$consensusPosition]->[1] = $contig.":".$revOutput->[$newx]->[1]->[$index{$contig}];
					$consensusPosition++;
				}
				else {
					push @{$nextConsensusContig->[0]}, $revOutput->[$newx]->[0]->[$index{$contig}];
					push @{$nextConsensusContig->[1]}, $contig.":".$revOutput->[$newx]->[1]->[$index{$contig}];
				}
			}
			$newx++;
		}
		else {
			# we want to make sure that newx keeps up as we progress through the array
			$newx=$curr+1;
			undef $currentConsensusContig;
			if (defined $nextConsensusContig) {
				print STDERR "Oh crap, we have ", Dumper($nextConsensusContig), " as a nextcontig\n";
			}
			undef $nextConsensusContig;
			$consensus->[$consensusPosition++]->[0] = "N";
		}
	}
	$locations = &update_locations();
}









###########################################################################
#                                                                         #
# Print out the location file that has the starts and stops               #
#                                                                         #
###########################################################################

open(OUT, ">$opts{o}.locations.txt") || die "can't open $opts{o}.locations.txt: $!";
foreach my $contig (@contigs) {
	if ($locations->{$contig}) {print OUT join("\t", $contig, @{$locations->{$contig}}), "\n"}
	else {print OUT "$contig\n"}
}
close OUT;


###########################################################################
#                                                                         #
# Print out the map file that has the sequence, original, and new posns   #
#                                                                         #
###########################################################################
print STDERR "Printing out $opts{o}.map\n";
open(OUT, ">$opts{o}.map") || die "can't open $opts{o}.map for writing";
print OUT join("\t", "Position",  "Initial bases", "Initial Positions", "Initial Contigs");
if (defined $revOutput) {
	print OUT join("\t", "\tRevised bases", "Revised positions", "Revised contigs");
}
if (defined $consensus) {
	print OUT "\tConsensus base";
}
print OUT "\n";

for (my $i=1; $i<=$max; $i++) {
	if (!defined $output->[$i] || ref($output->[$i]) ne "ARRAY") {
		print OUT "$i\n";
		next;
	}

	my $sequences = join(" ", @{$output->[$i]->[0]}); 
	my $positions = join(" ", @{$output->[$i]->[1]});
	my $contigs   = join(" ", @{$output->[$i]->[2]});
	print OUT join("\t", $i, $sequences, $positions, $contigs);

	if (defined $revOutput->[$i]) {
		no warnings;
		my $rsequences = join(" ", @{$revOutput->[$i]->[0]});
		my $rpositions = join(" ", @{$revOutput->[$i]->[1]});
		my $rcontigs   = join(" ", @{$revOutput->[$i]->[2]});
		print OUT join("\t", "", $rsequences, $rpositions, $rcontigs, $revOutput->[$i]->[3]);
	}
	if (defined $consensus) {print OUT "\t", $consensus->[$i]->[0]}
	print OUT "\n";
}
# print out the rest of the consensus sequence
if (defined $consensus) {
	for (my $i=$max+1; $i<=$#$consensus; $i++) {
		print OUT "$i\t\t\t";
		if (defined $revOutput->[$i]) {print OUT "\t\t\t\t"}
		print OUT "\t", $consensus->[$i]->[0];
	}
}
close OUT;


###########################################################################
#                                                                         #
# Print out the consensus sequence                                        #
#                                                                         #
###########################################################################

{ # just scope this piece of work

	print STDERR "Printing the consensus\n";
	open(OUT, ">$opts{o}.consensus.locations") || die "Can't open $opts{o}.consensus";
	my $consensusSequence;
	my %consensusFasta;
	for (my $i=1; $i<=$#$consensus; $i++) {
		if (defined $consensus->[$i] && ref($consensus->[$i]) eq "ARRAY") {
			print OUT join("\t", $i, @{$consensus->[$i]}), "\n";
			$consensusSequence .= $consensus->[$i]->[0];
		} else {
			print OUT "$i\tN\n";
			$consensusSequence .= "N";
		}
	}

	my $consensusContig=0;
	if ($opts{z} && !defined $opts{s}) {$opts{s}=5000}
	if ($opts{s}) {
		map {$consensusContig++; $consensusFasta{$opts{o}."_Contig$consensusContig"}=$_} split /N{$opts{s},}/i, $consensusSequence;
	} else {
		$consensusContig++; 
		$consensusFasta{$opts{o}."_Contig$consensusContig"} = $consensusSequence;
	}


	if ($opts{i}) {
## Now we need to include all contigs that are not used in the scaffold
		my %index = map {($_=>1)} @contigs;
		foreach my $contig (sort {length($query->{$b}) <=> length($query->{$a})} keys %$query) {
			if (!$index{$contig}) {
				$consensusContig++; 
				$consensusFasta{$opts{o}."_Contig$consensusContig [original contig $contig]"}=$query->{$contig};
				if ($opts{v}) {print STDERR "$contig was not included in the scaffold or intermediate outputs but was included in the fasta file\n"}
			}
		}
	}
	
	open(FASTA, ">$opts{o}.consensus.fasta") || die "Can't open $opts{o}.consensus.fasta";
	my @contigsToPrint;
	if ($opts{z}) {
		@contigsToPrint = sort {length($consensusFasta{$b}) <=> length($consensusFasta{$a})} keys %consensusFasta;
	} else {
		@contigsToPrint = sort {&contignumberonly($a) <=> &contignumberonly($b)} keys %consensusFasta;
	}

	map {print FASTA ">$_\n$consensusFasta{$_}\n"} (sort {length($consensusFasta{$b}) <=> length($consensusFasta{$a})} keys %consensusFasta);
	close FASTA;
	close OUT;

	&print_statistics(\%consensusFasta);
}

					
###########################################################################
#                                                                         #
# Check whether we are at the end of a contig                             #
#                                                                         #
###########################################################################


sub check_ends {
	my ($i ,$inds)=@_;
	my %dontuse;
	foreach my $index (0 .. $#$inds) {
		if (defined $inds->[$index] && !(defined $sequence->[$i]->{$inds->[$index]})) {
			# we have a contig at this location, but we don't need the sequence for it any more
			# remember the stop position
			if ($inds->[$index] && !(defined $locations->{$inds->[$index]}->[1])) {
				$locations->{$inds->[$index]}->[1] = $i-1;
			}
			undef $inds->[$index]; # we don't need this position in the array
			$dontuse{$index}=$i; # and we don't want to use it again for at least 5 spaces
		}
	}

	return if (!defined $sequence->[$i]); # we only need to check if we are at an end here.

	# now we need to figure out where the contigs at this location should fit

	# delete any positions we are free to use
	map {delete $dontuse{$_} if (($i - $dontuse{$_}) > 5)} keys %dontuse;

	foreach my $contig (keys %{$sequence->[$i]}) {
		if (!defined $index{$contig}) {
			# we need to record this as the start of the contig
			$locations->{$contig}->[0] = $i;
			# we need to assign this a new index. Note that we use @$indices here, so we can assign it to the next free position
			for my $index (0 .. @$inds) {
				if (!defined $inds->[$index] && !$dontuse{$index}) {
					$index{$contig}=$index; 
					$inds->[$index]=$contig; 
					last;
				}
			}
		}
	}
	return $inds;
}



###########################################################################
#                                                                         #
# Update the locations if we edited the sequences                         #
#                                                                         #
###########################################################################

sub update_locations {
	my $newlocations;
	my %inlocs;
	for (my $curr=1; $curr <= $max; $curr++) {
		my $bases=[]; my $posns=[]; my $conts=[]; my %index;
		if (defined $revOutput->[$curr]) {
			$conts = $revOutput->[$curr]->[2];
		}
		my %theselocs;
		foreach my $contig (@$conts) {
			next if (!defined $contig);
			$theselocs{$contig}=1;
			unless (defined $newlocations->{$contig}->[0]) {$newlocations->{$contig}->[0] = $curr; $inlocs{$contig}=1}
		}
		foreach my $contig (keys %inlocs) {
			unless (defined $theselocs{$contig}) {$newlocations->{$contig}->[1] = $curr-1; delete $inlocs{$contig}}
		}
	}
	return $newlocations;
}
			
		
###########################################################################
#                                                                         #
# Convert an array of sequences to their relevant IUPAC codes             #
#                                                                         #
###########################################################################

sub iupac {
	my $array=shift;
	my %array;
	foreach my $base (@$array) {
		next if (!$base || $base eq " " || $base eq "-");
		$base=uc($base);
		if ($base eq "A" || $base eq "T" || $base eq "G" || $base eq "C") {$array{$base}=1}
		elsif ($DNA_letter_can_be{$base}) {
			map {$array{$_}=1} @{$DNA_letter_can_be{$base}};
		}
		else {
			die "Unknown base in alignment file: $base\n";
		}
	}

	my $code= join("", sort {$a cmp $b} grep {m/[GATC]+/} keys %array);
	return $iupac{$code} if (defined $iupac{$code});
	die "No IUPAC code for $code in ", Dumper(\%array);
}


###########################################################################
#                                                                         #
# Print some statstics about the data                                     #
#                                                                         #
###########################################################################


sub print_statistics {
	my $cfast = shift; # consensus fasta seqeunce
	open (STAT, ">$opts{o}.statistics.txt") || die "can't open $opts{o}.statistics.txt";
	print STAT join("\t", "Metric", "Original Sequence", "Scaffold"), "\n";
	print STAT join("\t", "Number of contigs", scalar(keys %$query), scalar(keys %$cfast)), "\n";
	my ($olen, $slen)=(0,0);
	map {$olen += length($query->{$_})} keys %$query;
	map {$slen += length($cfast->{$_})} keys %$cfast;
	print STAT join("\t", "Total length", $olen, $slen), "\n";
	print STAT join("\t", "N50", &N50($query), &N50($cfast)), "\n";

	($olen, $slen)=(0,0);
	map {my $s = $query->{$_}; $s =~ s/[N\-]+//g; $olen += length($s)} keys %$query;
	map {my $s = $cfast->{$_}; $s =~ s/[N\-]+//g; $slen += length($s)} keys %$cfast;
	print STAT join("\t", "Length without N's", $olen, $slen), "\n";
	my @qkeys = sort {length($query->{$b}) <=> length($query->{$a})} keys %$query;
	my @skeys = sort {length($cfast->{$b}) <=> length($cfast->{$a})} keys %$cfast;
	print STAT join("\t", "shortest contig", $qkeys[$#qkeys]." (".length($query->{$qkeys[$#qkeys]})." bp)", 
		$skeys[$#skeys]." (".length($cfast->{$skeys[$#skeys]})." bp)"), "\n";
	print STAT join("\t", "Longest contig", $qkeys[0]." (".length($query->{$qkeys[0]})." bp)", $skeys[0]." (".length($cfast->{$skeys[0]})." bp)"), "\n";
	close STAT;
}


###########################################################################
#                                                                         #
# Just remove the number from the end of the contig name                  #
#                                                                         #
###########################################################################

sub contignumberonly {
	my $d=shift;
	if ($d =~ m/(\d+)$/) {return $1} 
	elsif ($d =~ m/(\d+)\s+\[/) {return $1}
	else {return $d}
}



###########################################################################
#                                                                         #
# Read a fasta file and return a hash of the sequences                    #
#                                                                         #
###########################################################################

sub read_fasta {
	my ($file)=@_;
	if ($file =~ /\.gz$/) {open(IN, "gunzip -c $file|") || die "Can't open a pipe to $file"}
	elsif ($file =~ /\.zip$/) {open(IN, "unzip -p $file|") || die "Can't open a pipe to $file"}
	else {open (IN, $file) || die "Can't open $file"}
	my %f; my $t; my $s;
	while (<IN>) {
		if (/\r/) {s/\r/\n/g; s/\n\n/\n/g;}
		chomp;
		if (/^>/) {
			s#^>##;
			if ($t) {
				$s =~ s/\s+//g;
				$f{$t}=$s;
				undef $s;
			}
			$t=$_;
		}
		else {$s .= $_}
	}
	$s =~ s/\s+//g;
	$f{$t}=$s;
	close IN;
	return \%f;
}
	
###########################################################################
#                                                                         #
# Reverse complement a DNA sequence                                       #
#                                                                         #
###########################################################################


sub rc {
	my ($seq)=@_;
	$seq =~ tr/GATCgatc/CTAGctag/;
	$seq = reverse $seq;
	return $seq;
}

###########################################################################
#                                                                         #
# Calculate the N50 from a fasta sequence hash of [ids => sequences]      #
#                                                                         #
###########################################################################

sub N50 {
	my ($fa)=@_;
	my %length;
	my $total;
	foreach my $k (keys %$fa) {
		my $len = length($fa->{$k});
		$length{$len}++;
		$total += $len;
	}

	my @contigsizes = sort {$b <=> $a} keys %length;
	my $currsize=0;
	while ($currsize < int($total/2)) {
		my $l = shift @contigsizes;
		$currsize += ($l * $length{$l});
	}
	return $contigsizes[0];
}

###########################################################################
#                                                                         #
# Methods for parsing clustalw output files. These methods were written   #
# by RAE over some time, and part of a Clustal.pm module that is not      #
# provided with the scaffold_builder source code, but is available        #
# on request. Feel free to adapt these methods to your own needs.         # 
#                                                                         #
###########################################################################



###########################################################################
#                                                                         #
# Parse a clustalw output file and store the information as a hash.       #
#                                                                         #
###########################################################################


sub parse_file {
	my ($file)=@_;
	my $clustal;
	open(IN, $file) || die "Can't open $file";
	my $header = <IN>;
	my $currindex=0; my $spaces; 
	$clustal->{'identities'} = undef;
	$clustal->{'index'} = undef;
	$clustal->{'sequences'} = undef;

	while (<IN>) {
		chomp;
		next if (/^\s+$/);
		if (!$spaces && /^(\S+\s+)\S+/) {$spaces = length($1)} # how many spaces are there before the alignment
		if (/^[\s\:\*\.]+/) {
			# this is the identity line
			s/^.{$spaces}//;
			$clustal->{'identities'} .= $_;
			next;
		}
		if (/^(\S+)\s+(\S+)$/) {
			my ($id, $seq)=($1, $2);
			unless (defined $clustal->{'index'}->{$id}) {$clustal->{'index'}->{$id}=$currindex; $currindex++}
			$clustal->{'sequences'}->{$id}.=$seq;
		}
	}
	close IN;

	my @ids = sort {$clustal->{'index'}->{$a} <=> $clustal->{'index'}->{$b}} keys %{$clustal->{'index'}};
	$clustal->{'ids'}=\@ids;

	return $clustal; # we are going to pass around this pointer to the hash so we can use it when we need.
}



###########################################################################
#                                                                         #
# Create an array of the alignments                                       #
#                                                                         #
###########################################################################


sub alignments {
	my ($clustal) = @_;
	my @result;
	foreach my $id (@{$clustal->{'ids'}}) {
		my @seq = split //, $clustal->{'sequences'}->{$id};
		push @result, \@seq;
	}
	return \@result;
}



###########################################################################
#                                                                         #
# Generate the trimmed_alignments from the sequence                       #
#                                                                         #
###########################################################################


sub trimmed_alignments {
	my ($clustal) = @_;
	my @result;
	my $fivetrim = 0; my $threetrim = 0;
	# first iterate through and get the lengths to trim
	foreach my $id (@{$clustal->{'ids'}}) {
		my $seq = $clustal->{'sequences'}->{$id};
		if ($seq =~ /^(\-+)/) {
			(length($1) > $fivetrim) ? ($fivetrim = length($1)) : 1;
		}
		if ($seq =~ /(\-+)$/) {
			(length($1) > $threetrim) ? ($threetrim = length($1)) : 1;
		}
	}
	foreach my $id (@{$clustal->{'ids'}}) {
		my $seq = $clustal->{'sequences'}->{$id};
		if ($threetrim) {$seq = substr($seq, 0, length($seq)-$threetrim)}
		if ($fivetrim)  {$seq = substr($seq, $fivetrim)}
		push @result, [split //, $seq];
		$clustal->{'trimmed_sequence'}->{$id}=$seq;
	}
		
	my $idents = $clustal->{'identities'};
	if ($threetrim) {$idents = substr($idents, 0, length($idents)-$threetrim)}
	if ($fivetrim && length($idents) > $fivetrim)  {$idents = substr($idents, $fivetrim)}
	$clustal->{'trimmed_identities'}=$idents;

	return \@result;
}


###########################################################################
#                                                                         #
# Calculate the percent of the alignments that are similar                #
#                                                                         #
###########################################################################



sub calc_percent {
	my ($alignments)=@_;
	my $n=0; my $same=0;
	for my $j (0 .. $#{$alignments->[0]}) {
		$n++;
		my %bases;
		for my $i (0 .. $#$alignments) {$bases{uc($alignments->[$i]->[$j])}=1}
		(scalar(keys %bases) == 1) && ($same++);
	}
	return sprintf("%.3f", (($same/$n) *100));
}


###########################################################################
#                                                                         #
# The percent of the sequence that is identical, after trimming           #
#                                                                         #
###########################################################################


sub trimmed_percent_identical {
	my ($clustal) = @_;
	return &calc_percent(&trimmed_alignments($clustal));
}


###########################################################################
#                                                                         #
# The total length of the alignment, including gaps                       #
#                                                                         #
###########################################################################


sub alignment_length {
	my ($clustal) = @_;
	my $seq = $clustal->{'sequences'}->{$clustal->{'ids'}->[0]};
	return length($seq);
}



=pod

This software is released under the PERL Artistic License. 

If you use this software please cite:

A simple alignment-based scaffolder allows combined de novo and reference-based assembly.
Bas E. Dutilh, Keri Elkins, T. David Matthews, Anca M. Segall, Elizabeth Dinsdale and Robert A. Edwards

Sotware written by, and copyright, Rob Edwards, 2011.

The license is appended below.



=head1 The "Artistic License"

				Preamble

The intent of this document is to state the conditions under which a
Package may be copied, such that the Copyright Holder maintains some
semblance of artistic control over the development of the package,
while giving the users of the package the right to use and distribute
the Package in a more-or-less customary fashion, plus the right to make
reasonable modifications.

Definitions:

	"Package" refers to the collection of files distributed by the
	Copyright Holder, and derivatives of that collection of files
	created through textual modification.

	"Standard Version" refers to such a Package if it has not been
	modified, or has been modified in accordance with the wishes
	of the Copyright Holder as specified below.

	"Copyright Holder" is whoever is named in the copyright or
	copyrights for the package.

	"You" is you, if you're thinking about copying or distributing
	this Package.

	"Reasonable copying fee" is whatever you can justify on the
	basis of media cost, duplication charges, time of people involved,
	and so on.  (You will not be required to justify it to the
	Copyright Holder, but only to the computing community at large
	as a market that must bear the fee.)

	"Freely Available" means that no fee is charged for the item
	itself, though there may be fees involved in handling the item.
	It also means that recipients of the item may redistribute it
	under the same conditions they received it.

1. You may make and give away verbatim copies of the source form of the
Standard Version of this Package without restriction, provided that you
duplicate all of the original copyright notices and associated disclaimers.

2. You may apply bug fixes, portability fixes and other modifications
derived from the Public Domain or from the Copyright Holder.  A Package
modified in such a way shall still be considered the Standard Version.

3. You may otherwise modify your copy of this Package in any way, provided
that you insert a prominent notice in each changed file stating how and
when you changed that file, and provided that you do at least ONE of the
following:

    a) place your modifications in the Public Domain or otherwise make them
    Freely Available, such as by posting said modifications to Usenet or
    an equivalent medium, or placing the modifications on a major archive
    site such as uunet.uu.net, or by allowing the Copyright Holder to include
    your modifications in the Standard Version of the Package.

    b) use the modified Package only within your corporation or organization.

    c) rename any non-standard executables so the names do not conflict
    with standard executables, which must also be provided, and provide
    a separate manual page for each non-standard executable that clearly
    documents how it differs from the Standard Version.

    d) make other distribution arrangements with the Copyright Holder.

4. You may distribute the programs of this Package in object code or
executable form, provided that you do at least ONE of the following:

    a) distribute a Standard Version of the executables and library files,
    together with instructions (in the manual page or equivalent) on where
    to get the Standard Version.

    b) accompany the distribution with the machine-readable source of
    the Package with your modifications.

    c) give non-standard executables non-standard names, and clearly
    document the differences in manual pages (or equivalent), together
    with instructions on where to get the Standard Version.

    d) make other distribution arrangements with the Copyright Holder.

5. You may charge a reasonable copying fee for any distribution of this
Package.  You may charge any fee you choose for support of this
Package.  You may not charge a fee for this Package itself.  However,
you may distribute this Package in aggregate with other (possibly
commercial) programs as part of a larger (possibly commercial) software
distribution provided that you do not advertise this Package as a
product of your own.  You may embed this Package's interpreter within
an executable of yours (by linking); this shall be construed as a mere
form of aggregation, provided that the complete Standard Version of the
interpreter is so embedded.

6. The scripts and library files supplied as input to or produced as
output from the programs of this Package do not automatically fall
under the copyright of this Package, but belong to whoever generated
them, and may be sold commercially, and may be aggregated with this
Package.  If such scripts or library files are aggregated with this
Package via the so-called "undump" or "unexec" methods of producing a
binary executable image, then distribution of such an image shall
neither be construed as a distribution of this Package nor shall it
fall under the restrictions of Paragraphs 3 and 4, provided that you do
not represent such an executable image as a Standard Version of this
Package.

7. C subroutines (or comparably compiled subroutines in other
languages) supplied by you and linked into this Package in order to
emulate subroutines and variables of the language defined by this
Package shall not be considered part of this Package, but are the
equivalent of input as in Paragraph 6, provided these subroutines do
not change the language in any way that would cause it to fail the
regression tests for the language.

8. Aggregation of this Package with a commercial distribution is always
permitted provided that the use of this Package is embedded; that is,
when no overt attempt is made to make this Package's interfaces visible
to the end user of the commercial distribution.  Such use shall not be
construed as a distribution of this Package.

9. The name of the Copyright Holder may not be used to endorse or promote
products derived from this software without specific prior written permission.

10. THIS PACKAGE IS PROVIDED "AS IS" AND WITHOUT ANY EXPRESS OR
IMPLIED WARRANTIES, INCLUDING, WITHOUT LIMITATION, THE IMPLIED
WARRANTIES OF MERCHANTIBILITY AND FITNESS FOR A PARTICULAR PURPOSE.

				The End

=cut
