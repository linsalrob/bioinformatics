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

my $protdistexe = '/home/redwards/bin/prodist';

my $usage = "protdistonly.pl <input dir of clustal files>\n";
my $indir = shift || die $usage;

unless (-e "$indir.protdist") {mkdir "$indir.protdist", 0755}
unless (-e "yes") {open YES, ">yes"; print YES "y\n"; close YES}

opendir(DIR, $indir) || die "Can't open $indir for reading\n";
chdir "$indir";
while (my $file = readdir(DIR)) {
	next if ($file =~ /^\./);
	next if ($file =~ /\.dnd$/);
	system "cp -f $file infile";
	system "/home/redwards/bin/prodist < ../yes";
	system "mv outfile ../$indir.protdist/$file";
}
system("rm -f infile");
