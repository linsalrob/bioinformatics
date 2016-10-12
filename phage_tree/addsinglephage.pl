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



# addsinglephage.pl

# this will add a single phage genome into the alignment. I am starting from having a directory of blast results,
# and taking it from there. We will run a clustal/protdist, as before, and then add the protdist results to
# the bottom line of the protdist matrix. Then we can simply run neighbor with only the lower half of the matrix
# need to check that tree against the good ones!!!

# this is also being converted to store everything in /tmp so it can be run as a cgi script

use strict;
use DBI;

my $dbh=DBI->connect("DBI:mysql:phage", "SQLUSER", "SQLPASSWORD") or die "Can't connect to database\n";

my $usage = "blastsqltoclustal.pl <original sequences dir> <input dir of blasts> <input file of original matrix>\n";
my $origseqdir = shift || &niceexit($usage); # for cgi this will be the /tmp separatemultifasta dir
my $dir = shift || &niceexit($usage); $dir =~ s/\/$//; # for cgi this will be generated from the blast results
my $matrixfile = shift || &niceexit($usage); # for cgi this will have to be stored somewhere and updated from time to time

my $outdir = "/tmp/addingphage".$$;
mkdir $outdir, 0755;

my $currentdir = `pwd`; chomp $currentdir;

my $numberofgenomes = &getnumgenomes;
&run_clustal; # this reads the blast files, makes the clustal input files, runs clustal, and then runs protdists
&add_newdata; # this reads the protdist files and outputs the matrix as infile
&run_neighbor; # this just runs neighbor on the protdists file
&cleanup;  # this takes the output file and prepares the colored html



&niceexit(0);


sub run_clustal { # this is the start of the old blastsqltoclustal program. This bracket will localize everything!
print STDERR "Parsing BLAST hits and getting sequences\n";

opendir (DIR, $dir) || &niceexit("Can't open $dir\n");
while (my $file = readdir(DIR)) {
	next if ($file =~ /^\./);
	open (IN, "$dir/$file") || &niceexit("Can't open $dir/$file\n");
	open (OUT, ">$outdir/$file") || &niceexit("Can't open $outdir/$file\n");
	my $get;
	while (my $line = <IN>) {
		if ($line =~ /Query=\s+(\S+)/) {
			my $querysequence=$1;
			open (SEQ, "$origseqdir/$querysequence") || die "Can't find $origseqdir/$querysequence for $file\n";
			while (my $seq = <SEQ>) {
				if ($seq =~ /^>/) {$seq =~ s/_/-/g} # this just removes _ which is used as a delimiter later.
				print OUT $seq;
				}
			close SEQ;
			}
		if ($line =~ /Sequences producing/) {$get = 1; next}
		next unless ($get);
		if ($line =~ /^>/) {last}
		next unless ($line =~ /\d_\d/);
		my @a = split (/\s+/, $line);
		my ($gene, $source) =split (/_/, $a[0]);
		unless ($gene && $source) {die "problem parsing $line\n"}
		my $exc = $dbh->prepare("SELECT translation from protein where count = $gene" ) or croak $dbh->errstr;
		$exc->execute or die $dbh->errstr;
		while (my @retrieved = $exc->fetchrow_array) {print OUT ">$a[0]\n$retrieved[0]\n"}
	}
	close IN; close OUT;
}
closedir(DIR);

print STDERR "Parsing complete. Results are in $outdir\nRunning Clustal\n";

unless (-e "$outdir.clustal") {mkdir "$outdir.clustal", 0755}
unless (-e "$outdir.protdist") {mkdir "$outdir.protdist", 0755}
unless (-e "/tmp/yes") {open TMP, ">/tmp/yes"; print TMP "y\n"; close TMP}


opendir(DIR, $outdir) || &niceexit("Can't open $outdir for reading\n");
while (my $file = readdir(DIR)) {
	next if ($file =~ /^\./);
	next if ($file =~ /\.dnd$/);
	system "/usr/local/genome/bin/clustalw -INFILE=$outdir/$file -OUTFILE=$outdir.clustal/infile -OUTPUT=PHYLIP";
	chdir "$outdir.clustal";
	system "/usr/local/genome/bin/protdist < /tmp/yes";
	system "mv infile $file";
	system "mv outfile $outdir.protdist/$file";
	chdir $currentdir;
	}

}





sub add_newdata {

# this is the end of the old blast to clustal program. Now we are going to try and merge in the combineprotdists
# program

# in this case, we need to read in each matrix file and figure out which sequence is our new one (won't be\d+_\d+)
# note for cgi have to make sure that this does not occur in tag line.

# then we need to figure out what the scores of ours against all the others were and store them in an array
# finally we need to add that to the bottom of the matrix.

my @matrix;
my @count;
my @average;
# read each file of the protdists one at a time, and add the data to an array
opendir(DIR, "$outdir.protdist") || &niceexit("Can't open $outdir.protdist\n");
while (my $file=readdir(DIR)) {
	next if ($file =~ /^\./);
	open (IN, "$outdir.protdist/$file") || &niceexit("Can't open $outdir.protdist/$file");
	my @protdists;
	my @genomes = ('0');
	# we need to know all the genomes before we can store the data. Therefore
	# read and store each line in @dists
	# then get all the genome numbers and store them in an array
	while (my $line = <IN>) {
		next if ($line =~ /^\s+/);
		chomp($line);
		# added loop to get the genome number from the database
		my @line = split (/\s+/, $line);
		unless ($line[0] =~ /_/) {
			# this is the line that has the data, we need to make sure the input file's
			# don't have _ in their names. (line 47)
			# we also need to make sure that 0's are saved for later
			foreach my $y (1 .. $#line) {unless ($line[$y]) {$line[$y] = -2}}
			@protdists=@line;
			push (@genomes, '0');
			next;	
		}
		else {
			my ($gene, $genome) = split (/_/, $line[0]);
			push (@genomes, $genome)
		}
		
	}
	close IN;
	# now we will average the matrix based on the count.
	foreach my $y (1 .. $numberofgenomes) {
		#if protdist is -2, the real value is 0, so we don't need to increment average, only count.
		#if protdist is -1, the distance was too large to calculate, therefore it is 100
		unless ($protdists[$y]) {$average[$genomes[$y]] += 100; $count[$genomes[$y]]++; next}
		if ($protdists[$y] == '-1') {$average[$genomes[$y]] += 100; $count[$genomes[$y]]++}
		elsif ($protdists[$y] == '-2') {$count[$genomes[$y]]++}
		elsif ($average[$genomes[$y]]) {
			$average[$genomes[$y]] += $protdists[$y];
			$count[$genomes[$y]] ++;
		}
		else {
			$average[$genomes[$y]] = $protdists[$y];
			$count[$genomes[$y]] ++;
		}
	}
	
} # end of the readdir

# now quickly loop through and average everything.
foreach my $y (1 .. $#average) {
	next unless ($average[$y]);
	$average[$y] = $average[$y]/$count[$y];
	}

# put the matrix into a temp file, and add the new data as the last line
open (IN, $matrixfile) || &niceexit("Can't open $matrixfile");
open (OUT, ">/tmp/infile") || &niceexit("Can't open /tmp/infile for writing");
my $line;
while (<IN>) {
	unless ($line) { # add one to the number of genomes on the first line
		my $in = $_;
		chomp($in);
		$in++;
		print OUT "$in\n";
		$line++;
		next;
		}
	chomp;
	print OUT;
	if ($average[$line]) {print OUT $average[$line], " ", $count[$line], "\n"} else {print OUT "100 0\n"}
	$line++;
	}
close IN;
print OUT "NEWGENOME ";
foreach my $y (0 .. $line-1) {
	if ($average[$y]) {print OUT $average[$y], " ", $count[$y], "  "} else {print OUT "100 0  "}
	}
print OUT "0 $line-1\n";
close OUT;
}

sub run_neighbor {
	chdir "/tmp/";
	open (OUT, ">neighbor.input");
	print OUT "s\nj\n133\ny\n";
	close OUT;
	system "/usr/local/genome/bin/neighbor < neighbor.input";
	}


sub cleanup {
	# this is a combination of treegenometohtml.pl and color.pl

	my %genomename;
	my %color;
	
	my $exc = $dbh->prepare("select count, organism, family from phage" ) or croak $dbh->errstr;
	$exc->execute or die $dbh->errstr;
	while (my @retrieved = $exc->fetchrow_array) {
		my $temp = "genome".$retrieved[0];
		$genomename{$temp} = $retrieved[1]." (".$retrieved[2].")";
		}
	

	$color{'Tectiviridae'} = "#000FF";
	$color{'Leviviridae'} = "#B22222";
	$color{'Plasmaviridae'} = "#9932CC";
	$color{'Inoviridae'} = "#32CD32";
	$color{'Fuselloviridae'} = "#FF7F50";
	$color{'Corticoviridae'} = "#006400";
	$color{'Myoviridae'} = "#778899";
	$color{'Podoviridae'} = "#FF0000";
	$color{'Microviridae'} = "#B8860B";
	$color{'Siphoviridae'} = "#20B2AA";

	open (IN, "/tmp/outfile") || &niceexit("Can't open /tmp/outfile for writing\n");
	open (OUT, ">/tmp/outfile.html")  || &niceexit("Can't open /tmp/outfile.html for writing\n");
	print OUT "<html><head><title>Phage trees</title></head><body bgcolor=\"#FFFFFF\">\n";
	my $check;
	while (my $line = <IN>) {
		chomp ($line);
		
		# htmlize
		if ($line =~ /populations/i) {$line .= "\n<pre>\n"}
		if ($line =~ /remember/i) {$line =  "</pre>\n$line"}
		if ($line =~ /between/i) {$check=1; print OUT "<table>\n"; $line =~ s/^/     /}
#		if (/\d-\d{3,}/) {print STDERR "Corrected $_"; s/(\d)(-\d{3,})/$1     $2/g}
		if ($check) {
			next if ($line =~ /\=\=/);
			$line =~ s/\s+/<tr><td>/;
			$line =~ s/\s*$/<\/td><\/tr>\n/;
			$line =~ s/\s\s+/<\/td><td>/g;
			}
		
		# change genome name
		if ($line =~ /genome\d+/) {$line =~ s/(genome\d+)/$genomename{$1}/g}
		
		#colorize
		foreach my $key (keys %color) {
			if ($line =~ /$key/) {
				$line =~ s/^/<font color=\"$color{$key}\">/;
				$line =~ s/$/<\/font>/;
				}
			}
		print OUT "$line\n";
	}
	
	print OUT "\n\n</table><p><p><hr><p><p>\n\nCOLOR CODES:\n<p>\n";
	foreach my $key (sort {$a cmp $b} keys %color) {print OUT "<font color=\"$color{$key}\">$key</font><br>\n"}
	print OUT "</body></html>";
}	

sub getnumgenomes {
	open (IN, $matrixfile) || die "Can't open $matrixfile for reading\n";
	my $number = <IN>;
	until ($number =~ /\d+/) {$number = <IN>}
	chomp ($number);
	$number =~ s/\s+//g;
	return $number;
	}
	

sub niceexit {
	my $reason = shift;
	$dbh->disconnect;
	if ($reason) {print STDERR $reason; exit(-1)}
	else {exit(0)}
	}
