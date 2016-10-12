#!/usr/bin/perl -w

use strict;

use Rob;
my $rob=new Rob;
use Data::Dumper;
use Getopt::Std;
my %opts;
getopts('acf:k:m:p:q:r:n:v', \%opts);

die <<EOF unless ($opts{q});
$0
-q fastq file of sequences

-m minimum sequence length after primer removal

-p primer sequence (default: CTGCTGTACGGCCAAGGCGGATGTACGGTACAGCAG). You only need to specify one direction, the rc will be generated for you.

-k kmer size to split the primer up to look for substr matches (default = 6)

-n number of k-mer hits to still split on. Default is 16 6-mers which means there are three mismatches. We take the smallest index as the beginning, and then should adjust that incase the first base(s) are mis read, and then add the expected length of the primer (in case the last bases are misread).
	Default values depending on the k-mer size
		k=6 : n=14
		k=7 : n=13
		k=8 : n=12

-r reverse the left end of the sequence

-a write a random set of sequences out for alignment. This is really for development purposes!
-v verbose output. Explains the fate of each sequence

-c check the kmers produced for each primer and exit
EOF

unless (defined $opts{k}) {$opts{k} = 6}
unless (defined $opts{m}) {$opts{m} = 0}
unless (defined $opts{n}) {
	($opts{k} == 6) ? ($opts{n} = 14) : ($opts{k} == 7) ? ($opts{n} = 13) : ($opts{k} == 8) ? ($opts{n} = 12) : ($opts{n} = 16);
}


my @primers = ('CTGCTGTACGGCCAAGGCGGATGTACGGTACAGCAG', 'CTGCTGTACCGTACATCCGCCTTGGCCGTACAGCAG'); # IA_A and IA_B primers respectively
$opts{q} =~ /([\.\-\w]+)\.fastq/;
my $libname = $1;

if ($opts{p}) {
	@primers = (uc($opts{p}), $rob->rc(uc($opts{p})));
}
my $length = length($primers[0]);

my $allsequences;

# figure out the k-mers for each primer
my $kmers;
foreach my $p (@primers) {
	my $l=0;
	my @mers;
	while ($l < length($p) - $opts{k} + 1) {
		push @mers, substr($p, $l, $opts{k});
		$l++;
	}
	
	#my $c=0; print ">p\n$p\n", map {" " x $c++ . "$_\n"} @mers;
	#die "There are ", scalar(@mers), " primers", Dumper($p, \@mers), "\n";
	$kmers->{$p}=\@mers;
}

# check the kmers and exit
if ($opts{c}) {
	print Dumper($kmers);
	exit;
	foreach my $p (keys %$kmers) {
		print "$p\n";
		my $c=0;
		map {print " " x $c++, $_, "\n"} @{$kmers->{$p}};
		print "\n";
	}
	exit;
}

my $dest = "finding_lmp_primers.$$";
mkdir $dest, 0755;

open(XML, ">$dest/primer_finding.xml") || die "Can't write the parameter to an xml file";
my $xmlfile;
if ($opts{f}) {$xmlfile = "<filename format=\"fasta\">$opts{f}</filename>"}
elsif ($opts{q}) {$xmlfile = "<filename format=\"fastq\">$opts{q}</filename>"}
print XML <<XML;
<primerFinding>
	<primers>
		<primer1 direction="forward">$primers[0]</primer1>
		<primer2 direction="reverse">$primers[1]</primer2>
	</primers>
	<kmerSearches>
		<minKmer>$opts{k}</minKmer>
		<numKmers>$opts{n}</numKmers>
	</kmerSearches>
	<input>
		$xmlfile
	</input>
	<output>
		<minLengthToKeep>$opts{m}</minLengthToKeep>
	</output>
</primerFinding>
XML

close XML;

if ($opts{a}) {mkdir "$dest/alignments", 0755; mkdir "$dest/start_postitions", 0755}

open(CL, ">$dest/cleaned.fastq") || die "can't open cleaned.fastq";
open(UM, ">$dest/unmatched.fastq") || die "can't open unmatched.fastq";
open(PRIMER, ">$dest/primer_sequences.fasta") || die "can;t open primer sequences";
if ($opts{v}) {open(VERBOSE, ">$dest/verbose_output.txt") || die "Can't open verbose_output.txt"}
my $primercount=0;
my $sequences;
while ( $sequences = $rob->stream_fastq($opts{q}) ) {
#foreach my $id (keys %$allsequences) {#}
	foreach my $id (keys %$sequences) {
		my $seq;
		#if ($opts{q}) {$seq = uc($allsequences->{$id}->[0])}
		if ($opts{q}) {$seq = uc($sequences->{$id}->[0])}
		else {$seq = uc($allsequences->{$id})}
		my $idx = index($seq, $primers[0]);
		if ($idx > -1) {
			&printout($id, $seq, $idx, $length, "Exact");
			if ($opts{v}) {print VERBOSE "$id found exact match to $primers[0].\n"}
			next;
		}
		
		$idx = index($seq, $primers[1]);
		if ($idx > -1) {
			&printout($id, $seq, $idx, $length, "Exact (rc)");
			if ($opts{v}) {print VERBOSE "$id found exact match to $primers[1].\n"}
			next;
		}

# We need to see if there is a mismatch in the primer
# The approach we are going to take is to see if we have enough 6-mers of
# our primer to suggest that it is there, and then find the correct start
		my $primer_found=0;
		foreach my $key (@primers) {
			my $start = &primer_with_mismatch($id, $seq, $key);
			next unless (defined $start); # this can be 0!

# my $aseq = substr($seq, 0, $start) . "*" . substr($seq, $start,);
# print STDERR "\nFor $id and $key we think the start should be at ", $start , "\n";
# print STDERR "$seq\n$aseq\n";
			&printout($id, $seq, $start, $length, "Mismatch");
			if ($opts{v}) {print VERBOSE "$id found with mismatch\n"}
			$primer_found=1;
			last;

		}
		next if ($primer_found);
		if ($opts{v}) {print VERBOSE "$id Not found\n"}
		print UM "\@$id [unmatched]\n$seq\n+\n", $sequences->{$id}->[1], "\n";
	}
}


# To find a primer with a mismatch, we see if we have enough kmers that are
# exact matches, and then identify the start
sub primer_with_mismatch {
	my ($id, $seq, $primer) = @_;
	my @starts=(); 
	$#starts = length($seq);
	map {$starts[$_]=0} (0 .. $#starts);
	my $max= [0,0];

# first, see how many kmers match
# and record the start sites of all the kmers
	my $n=0;
	my %match;
	foreach my $i (0 .. $#{$kmers->{$primer}}) {
		my $subprimer = $kmers->{$primer}->[$i];
		my $index = index($seq, $subprimer);
		
		while ($index > -1) {
			$match{$subprimer}++;
			if ($index - $i >= 0) {$starts[$index - $i]++}
			else {$starts[0]++}
			$index = index($seq, $subprimer, $index+1);
		}
	}

	if (scalar(keys %match) < $opts{n}) {
# not enough matches to continue
		if ($opts{v}) {print VERBOSE "$id ", scalar(keys %match), " for $primer. Not enough - skipped\n"}
		return undef;

	}

# now sum up all the kmers using a window of 4bp either side
# this is to accomodate frameshifts.
	my $i=0;
	while ($i < $#starts-4) {
		my $x=0;
		map {$x+=$starts[$_] if (defined $starts[$_])} ($i-2 .. $i+2);
		($x > $max->[0]) ? ($max = [$x, $i]) : 1;
		$i++;
	}

# max->[1] should be very close to the start, but if there is a mismatch it might
# be one or two bases before the start, so we will use index once more to get the exact
# start.

	foreach my $i (0 .. $#{$kmers->{$primer}}) {
		my $subprimer = $kmers->{$primer}->[$i];
		my $index = index($seq, $subprimer);
		if ($index > -1) {
			$max->[1] = $index-$i;
			if ($max->[1] < 0) {$max->[1] = 0}
			last;
		}
	}

	if ((rand(1) > 0.95) && $opts{a}) {
		my $c=0;
		while (-e "$dest/alignments/$n.$c.fa") {$c++}
		open (TOALN, ">$dest/alignments/$n.$c.fa") || die "Can't open $dest/alignments/$n.$c.fa";
		print TOALN ">$id\n$seq\n>primer\n$primer\n";
		close TOALN;
		system("clustalw $dest/alignments/$n.$c.fa > /dev/null 2>&1 &");
		open(BAD, ">$dest/start_postitions/$n.$c.fa") || die "Can't open $dest/start_postitions/$n.$c.fa";
		print BAD "$id\n$seq\n$primer\n\n";
		print BAD "@starts\n";
		print BAD "Maximum hits: ", $max->[0], "\n";
		print BAD "Predicted start: ", $max->[1], "\n";
		close BAD;
	}

	if ($opts{v}) {print VERBOSE "$id Predicted mismatched start at ", $max->[1], "\n"}
	return $max->[1];
}


sub printout {
	my ($id, $seq, $start, $length, $evidence)=@_;
# print out the left end, the primer, and the right end of the sequence


	my ($left, $right, $primer)=("", "", "");
	my ($lqual, $rqual); # the quality scores
# sometimes the primer is at the end of the sequence, so we ignore those cases
	if (($start + $length + 1) < length($seq)) {
		$left = substr($seq, 0, $start);
		$right = substr($seq, ($start+$length+1), );
		$primer = substr($seq, $start, $length+1);

		$lqual = substr($sequences->{$id}->[1], 0, $start);
		$rqual = substr($sequences->{$id}->[1], ($start+$length+1), );
	} else {
		$left = substr($seq, 0, $start);
		$lqual = substr($sequences->{$id}->[1], 0, $start);
		$primer = substr($seq, $start, );
	}

	if ($opts{r}) {	
		$left = $rob->rc($left);
		$lqual = reverse $lqual;
	}
	print CL "\@${id} template=$id dir=F library=$libname\n$left\n+\n$lqual\n" if (length($left) > $opts{m});
	print CL "\@${id} template=$id dir=R library=$libname\n$right\n+\n$rqual\n" if (length($right) > $opts{m});
	print PRIMER ">primer_", $primercount++, " [start: $start] [evidence: $evidence]\n$primer\n";

}

