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



# get similar phage

$dir = shift || die "Need a directory of blast results to work with\n";
$dir =~ s/\/$//;
opendir (DIR, $dir) || die "can't open $dir\n";
while ($file = readdir(DIR)) {
	next if ($file =~ /^\./);
	open (IN, "$dir/$file") || die "Can't open $dir/$file for reading\n";
	while ($line = <IN>) {
		if ($line =~ /Query=\s+(\S+)/) {($def, $gene) = split ("=", $1)}
		if ($line =~ /Sequences producing/) {
			$line = <IN>;
			until ($line =~ /^>/) {
				if ($line =~ /SWISS/) {$line = <IN>; next}
				unless ($line =~ /from/) {$line = <IN>; next}
				$line =~ /from(.*)/; 
				$match = $1;
				$match =~ s/\d+\s+\S+\s+\d+$//;
				push (@{$matches{$def}}, $match);
				$line = <IN>;
				}
			}
		if ($line =~ /^>/) {last}
		}
	close IN;
}

@keys = keys %matches;

foreach $key (@keys) {
	foreach $match (@{$matches{$key}}) {$count{$match}++}
	@matches = sort {$count{$b} <=> $count{$a}} keys %count;
	foreach $match (@matches) {print "$key\t$match\t$count{$match}\n"}
	}

				
