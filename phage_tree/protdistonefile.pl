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

my $protdistexe = '/home/redwards/bin/protdist';

my $usage = "$0 <file name>";
my $file = shift || die $usage;

my $outputfile = $file;
$outputfile =~ s/clustalw/protdist/;
if ($outputfile eq $file) {$outputfile .= ".protdist"}

# check if there is an infile, if there is just remove it :)

die "infile exists!" if (-e "infile"); # too nice

# but we need an outfile so that we are asked for the name of our file
unless (-e "outfile") {open OUT, ">outfile"; close OUT}

# fork a child to handle the processing
my $pid=fork();
if ($pid) {
# this is the parent
	waitpid($pid, 0);
}
elsif ($pid == 0) {
# this is the child
	open(PD, "|$protdistexe") || die "can't open a pipe to $protdistexe";
	print PD "$file\n";
	print PD "F\n";
	print PD "$outputfile\n";
	print PD "yes\n";
}
else {
	die "Not enough resources to fork";
}

print "\n";
exit;
