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



# clustalprotdist.pl

# once we have the fasta files run clustal and protdist only.

use strict;

my $usage = "clustalprotdist.pl <input dir of fastafiles> <options>\nOPTIONS\n-r # fraction of complete to report\n";
my $indir = shift || die $usage;
my $args = join (" ", @ARGV);
my $reportpercent;
if ($args =~ /-r (\d+)/) {$reportpercent=$1} else {$reportpercent=0}

unless (-e "$indir.clustal") {mkdir "$indir.clustal", 0755}
unless (-e "$indir.protdist") {mkdir "$indir.protdist", 0755}
unless (-e "yes") {open YES, ">yes"; print YES "y\n"; close YES}

my $filecount;
opendir(DIR, $indir) || die "Can't open $indir for reading\n";
while (my $file = readdir(DIR)) {
	next if ($file =~ /^\./);
	next if ($file =~ /\.dnd$/);
	$filecount++;
	if ($reportpercent) {unless ($filecount%$reportpercent) {print STDERR "$filecount done in ", time-$^T, " seconds\n"}}
	system "/usr/local/genome/bin/clustalw -INFILE=$indir/$file -OUTFILE=$indir.clustal/infile -OUTPUT=PHYLIP";
	chdir "$indir.clustal";
	unless (-e "infile") {die "INFILE does not exist before protdist while trying to parse $file\n"}
	system "/usr/local/genome/bin/protdist < ../yes";
	system "mv infile $file";
	system "mv outfile ../$indir.protdist/$file";
	chdir "..";
	}

