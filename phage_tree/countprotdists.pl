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



# for every protein distance, count the number of genomes and the number of unique genomes

use DBI;
use strict;

my $dbh=DBI->connect("DBI:mysql:phage", "apache") or die "Can't connect to database\n";
my $dir= shift || &niceexit("Need a dir to work with\n");


my %total; my %count;
opendir(DIR, $dir) || &niceexit("Can't open $dir\n");
while (my $file = readdir(DIR)) {
	next if ($file =~ /^\./);
        open (IN, "$dir/$file") || &niceexit("Can't open $dir/$file");
	my ($target, $targetgenome) = split (/_/, $file);
	my %genomecount;
	while (<IN>) {
		next unless (/\d+_\d+/);
		/(\d+)_(\d+)/; my ($protein, $genome) = split /_/;
		$genomecount{$genome}=1;
		$total{$target}++;
		}
	my @keys = keys %genomecount;
	$count{$target} = $#keys +1;
	
}

my %func;
foreach my $query (keys %count) {
        my $exc = $dbh->prepare("select count,organism,gene,function,product,proteinid from protein where count = '$query'" ) or croak $dbh->errstr;
        $exc->execute or die $dbh->errstr;
        my @ret = $exc->fetchrow_array;
        foreach my $x (2 .. $#ret) {
                next unless ($ret[$x]); next if ($ret[$x] eq "NULL");
                $func{$query} = $ret[$x];
                last;
                }
        unless ($func{$query}) {&niceexit("Can't get a function for $query using ", join ("|", @ret))}
        }
	
print "Unique genomes\tTotal genomes\tProtein\tFunc\n";
foreach my $orf (keys %count) {
print  "$count{$orf}\t$total{$orf}\t$orf\t$func{$orf}\n";
}

sub niceexit {
        my $reason = shift;
        $dbh->disconnect;
        if ($reason) {print STDERR $reason; exit(-1)}
        else {exit(0)}
        }

