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



# generate some data
use strict;
use DBI;

my $usage = "countoccurences.pl <protdist dir> <blast dir>\n";
my $dbh=DBI->connect("DBI:mysql:phage", "apache") or die "Can't connect to database\n";


my %filecount;
my $protdir= shift || die $usage;
my $blastdir= shift || die $usage;
opendir(DIR, $protdir) || &niceexit("Can't open $protdir");
my @protfiles = readdir(DIR);
closedir(DIR);
opendir(DIR, $blastdir) || &niceexit("Can't open $blastdir");
my @blastfiles = readdir(DIR);
closedir(DIR);

foreach my $file (@protfiles) {
	next if ($file =~ /^\./);
	open (IN, "$protdir/$file") || &niceexit("Can't open $protdir/$file");
	while (<IN>) {
		next unless (/\d+_\d+/);
		$filecount{$file}++;
		}
	}
	
my @files = sort {$filecount{$a} <=> $filecount{$b}} keys %filecount;

print "Most: ", $files[$#files], " : ", $filecount{$files[$#files]}, "\n";
print "Least: ", $files[0], " : ", $filecount{$files[0]}, "\n";

my $average; my $total; my %count; my %number;
foreach my $f (@files) {
	$number{$filecount{$f}}++;
	$average +=  $filecount{$f};
	$total++;
	$count{$filecount{$f}}++;
	}
$average = $average/$total;

my @mode = sort {$count{$b} <=> $count{$a}} keys %count;

print "Average: $average (from $total) files\n";
print "Mode: $mode[0] with $count{$mode[0]} files\n";


print "\n\n\nData for all the genomes\n\n";
foreach my $value (sort {$a <=> $b} keys %number) {print "$value\t$number{$value}\n"}

my %family;
my $exc = $dbh->prepare("select count, family from phage") or die $dbh->errstr;
$exc->execute or die $dbh->errstr;
while (my @retrieved = $exc->fetchrow_array) {$family{$retrieved[0]}=$retrieved[1]}

my %numberbyfam;
foreach my $f (@files) {
	my ($gene, $genome) = split /\_/, $f;
	${$numberbyfam{$family{$genome}}}{$filecount{$f}}++;
	}

my %prot; my %blast;
@blast{@blastfiles}=1;
foreach my $prot (@protfiles) {delete $blast{$prot}}
my @singles = keys %blast;
foreach my $single (@singles) {
	my ($gene, $genome) = split /\_/, $single;
	${$numberbyfam{$family{$genome}}}{'1'}++;
}

print "\n\n\nData by ICTV family 1\n\n";

my @family = (sort {$a cmp $b} keys %numberbyfam);
print "Count\t", join ("\t", @family), "\n";
foreach my $x (1 ..  $filecount{$files[$#files]}) {
	print "$x\t";
	foreach my $family (@family) {
	  if (exists ${$numberbyfam{$family}}{$x}) {print ${$numberbyfam{$family}}{$x},"\t"}
	  else {print "\t"}
	}
	print "\n";
}
 

print "\n\n\nData by ICTV family\n\n";
foreach my $family (@family) {
  print "$family\n";
  foreach my $value (sort {$a <=> $b} keys %{$numberbyfam{$family}}) {
    print "$value\t${$numberbyfam{$family}}{$value}\n";
    }
   print "\n\n\n";
 }


$dbh->disconnect;
