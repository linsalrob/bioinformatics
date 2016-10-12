#!/usr/bin/perl -w

#Copyright (C) 2005 Beltran Rodriguez-Brito,
#Pat MacNarnie, and Rober Edwards
#
#This program is free software; you can redistribute it and/or
#modify it under the terms of the GNU General Public License
#as published by the Free Software Foundation; either version 2
#of the License, or (at your option) any later version.
#
#This program is distributed in the hope that it will be useful,
#but WITHOUT ANY WARRANTY; without even the implied warranty of
#MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#GNU General Public License for more details.
#
#You should have received a copy of the GNU General Public License
#along with this program; if not, write to the Free Software
#Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301, USA.
#
#Released as part of the xipe package as supplemental online material
#to "A statistical approach for metagenome comparisons"
# by Beltran Rodriguez-Brito, Forest Rohwer, and Robert A. Edwards.


use strict;

my $file=shift || die "$0 <file>";
open(IN, $file) || die "can't open $file";
my $res;
while (<IN>) {
 chomp;
 my @a=split /\s+/;
 for (my $i=0; $i<=$#a;$i++) {
  $res->{$i}->{$a[$i]}++;
 }
}

foreach my $i (sort {$a <=> $b} keys %$res) {
 my $out;
 foreach my $a (sort {$a <=> $b} keys %{$res->{$i}}) {
  for (my $j=1; $j<=$res->{$i}->{$a}; $j++) {
   $out .= $a." ";
  }
 }
 $out =~ s/\s+$//;
 print $out, "\n";
}
