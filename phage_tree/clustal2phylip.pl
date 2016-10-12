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



$file = shift;
open (IN, $file) || die "Can't open $file";
open (OUT, ">$file.phylip") || die "Can't open $file.phylip";
$line=<IN>; chomp($line); print OUT $line;
while ($line=<IN>) {
	chomp($line);
	unless ($line =~ /^\s/) {
		print OUT "\n";
		@line = split /\s+/, $line;
		if (length($line[0])>10) {
			$line[0] =~ s/contig//i;
			$line[0] =~ s/_\d+\.\d+$//;
			}
		if (length($line[0])<10) {
			$spaces = " " x (10 - length($line[0]));
			$line[0] .= $spaces;
		}
		print OUT join (" ", @line);
		}
	else {$line =~ s/^\s+/ /; print OUT $line}
}

		
