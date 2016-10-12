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



# blastsqltoclustal.pl

# get all the protein sequences from the trans-genome blast, write them into a single file, and then
# pipe them into clustalw to generate an alignment

use strict;
use DBI;

my $dbh=DBI->connect("DBI:mysql:phage", "apache") or die "Can't connect to database\n";

my $usage = "blastsqltoclustal.pl <input dir of blasts> <options>\nOPTIONS\n\t-o output by organism\n\t-w output by whole set (default)\n";
my $dir = shift || &niceexit($usage);
my $args = join (" ", @ARGV);

print STDERR "Parsing BLAST hits and getting sequences\n";

opendir (DIR, $dir) || &niceexit("Can't open $dir\n");
my @files=readdir(DIR); 
print STDERR "Doing ", $#files+1, " files\n";

my $filecount; my %count;
foreach my $file (@files) {
	next if ($file =~ /^\./);
	open (IN, "$dir/$file") || &niceexit("Can't open $dir/$file\n");
	my $get; my $query;
	while (my $line = <IN>) {
		if ($line =~ /Query=  (\d+)_\d+/) {$query=$1}
		if ($line =~ /Sequences producing/) {$get = 1; next}
		next unless ($get);
		if ($line =~ /^>/) {last}
		next unless ($line =~ /\d_\d/);
		$count{$query}++;
	}
}
closedir(DIR);

my @countsbyorder = sort {$count{$b} <=> $count{$a}} keys %count;


my %genome;
{
  my $exc = $dbh->prepare("select count,organism from phage" ) or croak $dbh->errstr;
  $exc->execute or die $dbh->errstr;
  while (my @ret = $exc->fetchrow_array) {$genome{$ret[0]}=$ret[1]}
}


my %org; my %func; my $maxorg=1;
foreach my $query (@countsbyorder) {
	my $exc = $dbh->prepare("select count,organism,gene,function,product,proteinid from protein where count = '$query'" ) or croak $dbh->errstr;
	$exc->execute or die $dbh->errstr;
	my @ret = $exc->fetchrow_array;
	$org{$query} = $ret[1];
	if ($ret[1] > $maxorg) {$maxorg=$ret[1]}
	foreach my $x (2 .. $#ret) {
		next unless ($ret[$x]); next if ($ret[$x] eq "NULL");
		$func{$query} = $ret[$x];
		last;
		}
	unless ($func{$query}) {&niceexit("Can't get a function for $query using ", join ("|", @ret))}
	}

if ($args =~ /-o/) { # sort by organism
  foreach my $x (1 .. $maxorg) {
   print "\n\n", $genome{$x}, "\n";
   foreach my $query (@countsbyorder) {
     next unless ($org{$query} == $x);
     print "$func{$query}\t", $count{$query}-1, "\n";
     }
  }
}
else { # just do the whole darn lot
  foreach my $query (@countsbyorder) {
     print "$func{$query}\t", $count{$query}-1, "\t", $genome{$org{$query}}, "\n";
     }
   }


&niceexit();







sub niceexit {
	my $reason = shift;
	$dbh->disconnect;
	print STDERR "Done in ", time-$^T, " seconds\n";
	if ($reason) {print STDERR $reason; exit(-1)}
	else {exit(0)}
	}
