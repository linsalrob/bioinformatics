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




open (IN, $ARGV[0]) || die "Need a file\n";
open (OUT, ">$ARGV[0].pro") || die "Can't open output file\n";
while (<IN>) {
	if (/definition/i) {
		s/DEFINITION//i; s/Bacteriophage//gi; s/coliphage//gi;
		s/complete\s+genome//i; s/phage//gi;
		s/\,//g; s/\.//g; s/^\s+//; s/\s+$//; chomp; s/\s+/_/g;
		$definition = $_;
		}
	if (/CDS.*\d/) {/(\d+).*?(\d+)/; $start=$1."-".$2}
	if (/\/gene=\"(.*)\"/) {$gene = $1; $gene =~ s/^\s+//; $gene =~ s/\s+$//; $gene =~ s/\s/_/g;}
	if (/translation=/) {
		s/\/translation=\"//i;
		if (/\"/) {s/\"//; s/\s//g; chomp; $trans = $_}
		else {
			$trans = $_; $_ = <IN>;
			until (/\"/) {$trans .= $_; $_ = <IN>}
			$trans .= $_;
			$trans =~ s/\s//g; $trans =~ s/\n//g; $trans=~ s/\"//g;
			}
		if ($gene && $definition) {$tag = &check_dup($definition, $gene)}
		elsif ($start && $definition) {$tag = &check_dup($definition, $start)}
		else {die "Gene: $gene Definition: $definition Start: $start for\n$trans\n"}
		print OUT ">$tag\n$trans\n";
		undef $trans; undef $gene; undef $start; $count++; undef $tag;
		}
	if (/\/\//) {undef $trans}
	}

print "$count proteins found\n";

sub check_dup {
	$tag = join ("=", @_);
	if ($check{$tag}) {
		$check{$tag}++;
		$tag = $tag.".".$check{$tag}
		}
	else {
		$check{$tag}++
		}
	return $tag;
	}
