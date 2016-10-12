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

my $file = shift || die "Need a file to work on\n";

open (IN, $file) || die "Can't open $file\n";
open (OUT, ">$file.txt") || die "Can't open $file.txt for writing\n"; 
 
my $offset;
while (<IN>) {
  if (/<BIG><STRONG>Viruses<\/STRONG>/) {$offset++; next}
  next unless ($offset);
  if (/<A HREF/i) {
 	/>(\w.*?)</;
	my $phage = $1;
	if (/( \(.*?\))/) {$phage .= $1}
	my $tab = "\t" x ($offset -1);
	print OUT $tab, $phage, "\n";
	}
 if (/<LI TYPE=circle/i) {$offset++}
 if (/<\/UL/i) {$offset--}
}


close IN; close OUT;
open (IN, "$file.txt") || die "Can't open $file.txt for reading\n";

my %phylo;
my @phylo = (' ');
while (<IN>) {
	chomp;
	my @line = split /\t/;
#	print " PHYLO: @phylo ELEMENTS: $#phylo\nLINE: @line ELEMENTS: $#line\n";
	unless ($#line > $#phylo) {splice (@phylo, $#line)}
	push (@phylo, $line[$#line]);
	$phylo{join (" --> ", @phylo)} = 1;
	}

open (OUT, ">$file.phylo") || die "Can't open $file.phylo\n";
print OUT join ("\n", sort {$a cmp $b} keys %phylo), "\n";
close IN; close OUT;	
		


