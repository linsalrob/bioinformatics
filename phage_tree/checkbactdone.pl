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

my ($dir1, $dir2, $dir3, $dir4)= @ARGV;

unless ($dir1 && $dir2 && $dir3 && $dir4) {die "checkbactdone.pl blastdir fastadir clustaldir protdistdir\n"}
 
opendir(DIR, $dir1) || die "Can't open $dir1\n";
my @blast = readdir(DIR);
closedir(DIR);
 
opendir(DIR, $dir2) || die "Can't open $dir2\n";
my @fasta = readdir(DIR);
closedir(DIR);
 
opendir(DIR, $dir3) || die "Can't open $dir3\n";
my @clustal = readdir(DIR);
closedir(DIR);

opendir(DIR, $dir4) || die "Can't open $dir4\n";
my @protdist = readdir(DIR);
closedir(DIR);

 
 
print STDERR $#blast+1, " files in $dir1; ", $#fasta+1, " files in $dir2; ", $#clustal+1, " files in $dir3; ", $#protdist+1, " files in $dir4\n";
 
my (%files1, %files2, %files3, %files4);
 
@files1{@blast}=(1);
@files2{@fasta}=(1);
@files3{@clustal}=(1);
@files4{@protdist}=(1);

 
open BLASTONLY, ">$dir1.no.fasta";
print "\n\n\n$dir1\n=========\n";
 
print STDERR "Checking blast hits for matching fasta\n";
# check for files that are BLAST results but no fasta file
foreach my $file (sort keys %files1) {
  next if (exists $files2{$file}); # we have a fasta file
  # countblast hits
  my $hits = countblasthits($file);
  next if ($hits == 1); # it only found itself, do not save
  print BLASTONLY "$file ";
  }
print BLASTONLY "\n";
close BLASTONLY;

print STDERR "Checking fasta files for matching clustal\n";
open OUT, ">$dir2.no.clustal"; 
 # check for files that have a fasta file but no clustal
 foreach my $file (sort keys %files2) {
  next if (exists $files3{$file}); # we have a fasta file
  #check how many sequences in the file
  open (FASTA, "$dir2/$file") || die "Can't open $dir2/$file while checking fasta\n";
  my $count;
  while (<FASTA>) {if (/^>/) {$count++}}; close FASTA;
  next if ($count==1);
  print OUT "$file ";
  }
print OUT "\n";
close OUT;


#check for clustal files that do not have a protdist file
print STDERR "Checking clustal files for matching protdist\n";
open (OUT, ">$dir3.no.protdist"); 
 # check for files that have a fasta file but no clustal
 foreach my $file (sort keys %files3) {
  next if (exists $files4{$file}); # we have a fasta file
  print OUT "$file ";
  }
print OUT "\n";
close OUT;
 







sub countblasthits {
  my $file=shift;
  open (BLAST, "$dir1/$file") || die "Can't open $dir1/$file while checking blast hits\n";
  my $get; my $count;
  while (<BLAST>) {
  	if (/^Sequences producing/) {$get=1}
	if (/^>/) {undef $get; last}
	if (/^\d+/) {$count++}
	}
 return $count;
 }
