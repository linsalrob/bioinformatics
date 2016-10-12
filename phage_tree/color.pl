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




$color{'Tectiviridae'} = "#000FF";
$color{'Leviviridae'} = "#B22222";
$color{'Plasmaviridae'} = "#9932CC";
$color{'Inoviridae'} = "#32CD32";
$color{'Fuselloviridae'} = "#FF7F50";
$color{'Corticoviridae'} = "#006400";
$color{'Myoviridae'} = "#778899";
$color{'Podoviridae'} = "#FF0000";
$color{'Microviridae'} = "#B8860B";
$color{'Siphoviridae'} = "#20B2AA";

while (<>) {
	unless ($firstline) {unless (/html/) {print "<html><head><title>Phage tree</title></head><body>\n<pre>\n"} $firstline=1}
	foreach $key (keys %color) {
		if (/$key/) {s/^/<font color=\"$color{$key}\">/; s/$/<\/font>/}
		}
	print;
	}


print "\n\n<p><p><hr><p><p>\n\nCOLOR CODES:\n<p>\n";
foreach $key (sort {$a cmp $b} keys %color) {print "<font color=\"$color{$key}\">$key</font><br>\n"}

print "</pre></body></html>\n";


