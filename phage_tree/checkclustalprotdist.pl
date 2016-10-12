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


use strict;
# count files in directories, and check for number of blast hits


my ($blast, $fasta, $clustal, $prot) = @ARGV || die "checkclustalprotdist.pl <blast dir> <fasta dir> <clustal dir> <protdist dir>\n";

opendir (DIR, $blast) || die "Can't open $blast\n";
my @blastfiles = readdir(DIR);
closedir DIR;

opendir (DIR, $fasta) || die "Can't open $fasta\n";
my @fastafiles = readdir(DIR);
closedir DIR;

opendir (DIR, $clustal) || die "Can't open $clustal\n";
my @clustalfiles = readdir(DIR);
closedir DIR;

opendir (DIR, $prot) || die "Can't open $prot\n";
my @protfiles = readdir(DIR);
closedir DIR;

print "TOTALS\nBLAST: ", $#blastfiles+1, "\nFasta: ", $#fastafiles+1, "\nClustal: ", $#clustalfiles+1, "\nProtdist: ", $#protfiles+1, "\n";

my %seen;
