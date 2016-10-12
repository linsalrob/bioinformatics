#!/usr/bin/perl -w

# trim repeats. Save only those that start or stop between given coordinates

# Copyright Rob Edwards
#
# You may distribute this script under the GPL. See below.

=head1 NAME

trim_repeats.pl - remove repeats not in a specific region

=head1 SYNOPSIS

You can run this script from the command line with the command

perl -w trim_repeats.pl <filename> <from> <to>

=head1 DESCRIPTION

This script will read in a genbank file and remove any repeats that are outside of
a given range. I use it to cut down the number of repeats in a file.

=head1 AUTHOR -  Rob Edwards

Email redwards@utmem.edu

The author gratefully acknowledges all the help from the members of the BioPerl team.

=cut




use strict;

my ($file, $start, $stop)=@ARGV;
unless ($file && $start && $stop) {die "$0 <file name> <from> <to>"}

if ($start>$stop) {($start, $stop)=($stop, $start)}

open (IN, $file) || die "Can't open $file";

my $out=$file;
$out =~ s/.gbk$/.betwn.gbk/;

my $skip;
open (OUT, ">$out") || die "Can't open $out for writing";
while (my $l=<IN>) {
 if ($skip) {
  next if ($l =~ /^\s+\//);
  undef $skip;
 }
 unless ($l =~ /^\s+repeat/) {print OUT $l; next}
 if ($l =~ /repeat_region/) {print OUT $l; next}
 # repeat has the format join(411..423,13244..13256)
 $l=~m/(\d+)\.\.(\d+)\,(\d+)\.\.(\d+)/;
 my ($one, $two, $three, $four)=($1, $2, $3, $4);
 unless ($one && $two && $three && $four) {die "Can't parse numbers from $l\n"}
 my $keep;
 foreach my $pos ($one, $two, $three, $four) {
  if ($pos > $start && $pos < $stop) {$keep=1}
 }
 if ($keep) {print OUT $l; next}
 else {$skip=1}
}

#    Copyright 2002-2004 Rob Edwards
#    For updates, more information, or to discuss the scripts
#    please contact Rob Edwards at redwards@utmem.edu or via http://www.salmonella.org/
#
#    This file is part of Rob's Scripts developed by Rob Edwards.
#
#    These scripts are free software; you can redistribute and/or modify
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
