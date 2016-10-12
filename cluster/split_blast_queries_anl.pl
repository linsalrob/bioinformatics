#!/usr/bin/perl -w

# split a fasta file into a specified number of sub files and blast it against a database

use strict;
#use POSIX;
use Schedule::SGE;

my ($blastexec)=('/vol/cee/bin/blastall');


my $usage=<<EOF;

$0 <options>

-f file to split
	-n number to break into OR 
	-max maximum number of sequences per file  Can Only be One Or The Other***
-d destination directory (default = ".")

-p blast program
-db blast database
-ex blast executable (default is $blastexec)

-N job name (default is blastfile) 

-rev reverse the order of files that are submitted to the BLAST queue. (i.e. so you can run twice and start from the end backwards!)

Other things will be used as blast options. Unless -p && -db, will split and stop

EOF

my ($file,$dest, $no, $blastp, $blastdb, $jobname, $rev, $required);
my $blastopt=" ";

while (@ARGV) {
 my $test=shift(@ARGV);
 if ($test eq "-f") {$file=shift @ARGV}
 elsif ($test eq "-d") {$dest=shift @ARGV}
 elsif($test eq "-n") {$no=shift @ARGV}
 elsif($test eq "-max") {$required=shift @ARGV}
 elsif($test eq "-p") {$blastp=shift @ARGV}
 elsif($test eq "-db") {$blastdb=shift @ARGV}
 elsif($test eq "-ex") {$blastexec=shift @ARGV}
 elsif($test eq "-N") {$jobname=shift @ARGV}
 elsif($test eq "-rev") {$rev=1}
 else {$blastopt .= " " . $test . " " . shift @ARGV}
}
die $usage unless ($file && ($no || $required));

$jobname="blast$file" unless ($jobname);
if ($dest) {unless (-e $dest) {mkdir $dest, 0755}}
else {$dest="."}
$dest =~ s/\/$//;


# read the file and see how many > we have
if ($file =~ /gz$/) {open(IN, "gunzip -c $file |") || die "Can't open a pipe to $file"}
else {open(IN, $file)|| die "Can't open $file"}

my $counttags;
while (<IN>) {$counttags++ if (/^>/)}
close IN;
#my $required=ceil($counttags/$no); # ceil rounds up so we should get less files than if we use int or real rounding.

if ($no) {$required=int($counttags/$no)+1}

print STDERR "There are $counttags sequences in $file and we are going to write $required per file\n";


my $filecount=1;
my @sourcefiles;
if ($file =~ /gz$/) {open(IN, "gunzip -c $file |") || die "Can't open a pipe to $file"}
else {open(IN, $file)|| die "Can't open $file"}
$file =~ s/^.*\///;
open (OUT, ">$dest/$file.$filecount") || die "Can't open $dest/$file.$filecount";
push @sourcefiles, "$file.$filecount";

my $sofar;
while (my $line=<IN>) {
 if ($line =~ /^>/) {$sofar++}
 if (($line =~ /^>/) && !($sofar % $required) && (($counttags - $sofar) > 20)) {
  # the last conditional is to make sure that we don't have a few sequences in a file at the end
  close OUT;
  $filecount++;
  open (OUT, ">$dest/$file.$filecount") || die "Can't open $dest/$file.$filecount";
  push @sourcefiles, "$file.$filecount";
 }
 print OUT $line;
}

print STDERR "Wrote to $filecount files in $dest\n";

unless ($blastp && $blastdb) 
{
	print STDERR "Can't blast because no program or database called\n";
	exit(0);
}

@sourcefiles=reverse @sourcefiles if ($rev);
my $pwd=`pwd`; chomp($pwd);

open(SH, ">bl$$.sh") || die "Can't write bl$$.sh";

print SH '#!/bin/bash';
print SH "\n/home/redwards/FIGdisk/env/cee/bin/blastall -p $blastp -d $blastdb -i $pwd/$dest/$file.\$SGE_TASK_ID -o $pwd/$dest/\$SGE_TASK_ID.$blastp $blastopt\n";

close SH;

if ($blastp eq "blastp" || $blastp eq "blastx" || $blastp eq "tblastx")
{
	if (-e "BLOSUM62")
	{
		print STDERR "great, found BLOSUM\n";
	}
	elsif (my $bls=`locate BLOSUM62`)
	{
		chomp($bls);
		print STDERR "Copying $bls .\n";
		`cp $bls .`;
	}
	else {die "Sorry, can't figure out where BLOSUM62 is and we need it"}
}

print STDERR "Submitting array job\nqsub -cwd -t 1-$filecount:1 -l low -e error -o error ./bl$$.sh\n";
unless (-e "error") {mkdir "error", 0755}
#`qsub -cwd -t 1-978:1 -l low -e error/ -o error/ ./blshit.sh`;
# Swap these two lines if you don't want to run on the low priority queue
print STDERR `qsub -cwd -t 1-$filecount:1 -e error -o error ./bl$$.sh`; 
#print STDERR `qsub -cwd -t 1-$filecount:1 -l low -e error -o error ./bl$$.sh`;

