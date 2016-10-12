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

my $dir = shift || die "Need a dir to work with\n";
$dir =~ s/\/$//;

if (-e "$dir.treefiles") {die "$dir.treefiles exists\n"} else {mkdir "$dir.treefiles", 0755}
if (-e "$dir.outfiles") {die "$dir.outfiles exists\n"} else {mkdir "$dir.outfiles", 0755}
if (-e "$dir.infiles") {die "$dir.infiles exists\n"} else {mkdir "$dir.infiles", 0755}

unless (-e "yes") {open OUT, "yes"; print OUT "y\n"; close OUT}

opendir(DIR, $dir) || die "can't open $dir\n";
chdir "$dir.infiles";
while (my $file=readdir(DIR)) {
	next if ($file =~ /^\./);
	open (IN, "../$dir/$file") || die "Can't open ../$dir/$file\n";
	open (OUT, ">infile") || die "Can't open infile\n";; 
	while (<IN>) {if (/^\d+_\d/) {s/^\d+_(\d+)/genome$1/} print OUT} close IN; close OUT;
	system "neighbor < ../yes";
	system "mv outfile ../$dir.outfiles/$file.outfile";
	system "mv treefile ../$dir.treefiles/$file.treefile";
	system "mv infile ../$dir.infiles/$file.infile";
	}
		

