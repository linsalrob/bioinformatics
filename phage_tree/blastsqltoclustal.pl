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



# blastsqltoclustal.pl

# get all the protein sequences from the trans-genome blast, write them into a single file, and then
# pipe them into clustalw to generate an alignment

use strict;
use DBI;

my $dbh=DBI->connect("DBI:mysql:phage", "SQLUSER", "SQLPASSWORD") or die "Can't connect to database\n";

my $usage = "blastsqltoclustal.pl <input dir of blasts> <output dir> <options>\nOPTIONS\n-r # fraction of complete to report\n";
my $indir = shift || &niceexit($usage);
my $outdir = shift || &niceexit($usage);
my $args = join (" ", @ARGV);
my $reportpercent;
if ($args =~ /-r (\d+)/) {$reportpercent=$1} else {$reportpercent=0}

if (-e $outdir) {&niceexit("$outdir already exists\n")} else {mkdir $outdir, 0755}

print STDERR "Parsing BLAST hits and getting sequences\n";

opendir (DIR, $indir) || &niceexit("Can't open $indir\n");
my @files=readdir(DIR); 
print STDERR "Doing ", $#files+1, " files\n";

my $filecount;
foreach my $file (@files) {
	next if ($file =~ /^\./);
	open (IN, "$indir/$file") || &niceexit("Can't open $indir/$file\n");
	my $get; my %output;
	while (my $line = <IN>) {
		if ($line =~ /Sequences producing/) {$get = 1; next}
		next unless ($get);
		if ($line =~ /^>/) {last}
		next unless ($line =~ /\d_\d/);
		my @a = split (/\s+/, $line);
		my ($gene, $source) =split (/_/, $a[0]);
		unless ($gene && $source) {die "problem parsing $line\n"}
		my $exc = $dbh->prepare("SELECT translation from protein where count = $gene" ) or croak $dbh->errstr;
		$exc->execute or die $dbh->errstr;
		while (my @retrieved = $exc->fetchrow_array) {$output{$a[0]} = $retrieved[0]}
	}
	my @keys = keys %output;
	next unless ($#keys);
	open (OUT, ">$outdir/$file") || &niceexit("Can't open $outdir/$file\n");
	foreach my $key (@keys) {print OUT ">$key\n$output{$key}\n"}
	close IN; close OUT;
}
closedir(DIR);

print STDERR "Parsing complete. Results are in $outdir\nRunning Clustal\n";

unless (-e "$outdir.clustal") {mkdir "$outdir.clustal", 0755}
unless (-e "$outdir.protdist") {mkdir "$outdir.protdist", 0755}
unless (-e "yes") {open YES, ">yes"; print YES "y\n"; close YES}

opendir(DIR, $outdir) || &niceexit("Can't open $outdir for reading\n");
while (my $file = readdir(DIR)) {
	next if ($file =~ /^\./);
	next if ($file =~ /\.dnd$/);
	$filecount++;
	if ($reportpercent) {unless ($filecount%$reportpercent) {print STDERR "$filecount done in ", time-$^T, " seconds\n"}}
	system "/usr/local/genome/bin/clustalw -INFILE=$outdir/$file -OUTFILE=$outdir.clustal/infile -OUTPUT=PHYLIP -MAXDIV=0";
	chdir "$outdir.clustal";
	system "/usr/local/genome/bin/protdist < ../yes";
	system "mv infile $file";
	system "mv outfile ../$outdir.protdist/$file";
	chdir "..";
	}

&niceexit();







sub niceexit {
	my $reason = shift;
	$dbh->disconnect;
	print STDERR "Done in ", time-$^T, " seconds\n";
	if ($reason) {print STDERR $reason; exit(-1)}
	else {exit(0)}
	}
