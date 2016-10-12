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
use DBI;


my $usage = "countblastspercent.pl <dir of blasts> <cutoff> <matrix file for output>\nCUTOFF: use 0.1 or 1e-10\n";
my $indir = shift || die $usage;
my $cutoff = shift || die $usage;
my$matrixfile = shift || die $usage;

if ($cutoff =~ /(\d*)e(-\d+)/i) {my $one; if ($1) {$one=$1} else {$one=1} $cutoff = $one*(10**$2)}
#print STDERR "Looking for less than $cutoff\n";

my $orfs = &getallorfs();

opendir (DIR, $indir) || die "Can't open $indir\n";

my @genomes;
while (my $file = readdir(DIR)) {
	next if ($file =~ /^\./);
	open (IN, "$indir/$file") || die "Can't open $indir/$file\n";
	my $count =0; my $file; my $query;
	while (<IN>) {
		if (/Query=\s+\d+_(\d+)/) {$query = $1; next}
		last if (/^>/);
		next unless (/\d+_\d+/);
		/\d+_(\d+)/;
		my $match=$1;
		my @part = split;
		if ($part[$#part-1]  =~ /(.*)e(-\d+)/i) {$part[$#part-1] = $1*(10**$2)}
		next unless ($part[$#part-1] < $cutoff);
		$genomes[$query][$match] ++;
		$count++
		}
	}


open (MAT, ">$matrixfile") || die "Can't open $matrixfile for writing\n";
print MAT "     ", $#genomes, "\n";
foreach my $y (1 .. $#genomes) {
	my $temp= "genome$y";
	my $spacestoadd = " " x (10 - length($temp));
	print MAT $temp, $spacestoadd;
	foreach my $x (1 .. $#genomes) {
		if ($x == $y) {print MAT "0 "; next}
		unless (defined $genomes[$y][$x]) {print MAT "100 "; next}
		unless ($genomes[$y][$x]) {print MAT "XXX "; next}
		my $norfs = $$orfs{$x};
		my $nmat = $genomes[$y][$x];
		my $percent = (1-($nmat/$norfs))*100;
if (($y == 1 && $x==7)||($x == 1 && $y ==7)) {print STDERR "y: $y, x: $x. genomes[$y][$x] = $genomes[$y][$x]; genomes[$y][$x]=$genomes[$x][$y], norfs: $norfs; nmat=$nmat; percent = $percent\n"}
		print MAT "$percent ";
		}
	print MAT "\n";
	}
		
sub getallorfs {
	my %orfs;
	my $dbh=DBI->connect('DBI:mysql:phage', "SQLUSER", "SQLPASSWORD") or die "Can't connect to database\n";
	my $exc = $dbh->prepare("SELECT organism from protein" ) or croak $dbh->errstr;
	$exc->execute or die $dbh->errstr;
	while (my @retrieved = $exc->fetchrow_array) {$orfs{$retrieved[0]}++}
	$dbh->disconnect;
	my @keys = keys %orfs;
	my $total;
	foreach my $key (@keys) {$total += $orfs{$key}}
	print STDERR $#keys+1, " Genomes checked and contain $total ORFS\n"; 
	return \%orfs;
	}
	
