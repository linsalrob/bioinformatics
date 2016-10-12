#!/usr/bin/perl -w

# split a fasta file into a specified number of sub files and blast it against a database

use strict;
use POSIX;

my $usage=<<EOF;

$0 <options>

-f file to split
-n number to break into
-d destination directory (default = ".")

-p blast program
-db blast database
-ex blast executable (default is /home/redwards/bin/blast/blastall)

-N job name (default is blastfile) 

-rev reverse the order of files that are submitted to the BLAST queue. (i.e. so you can run twice and start from the end backwards!)

Other things will be used as blast options. Unless -p && -db, will split and stop

EOF

my ($file,$dest, $no, $blastp, $blastdb, $jobname, $rev);
my $blastopt=" ";
my ($blastexec)=('/home/redwards/bin/blast/blastall');

while (@ARGV) {
 my $test=shift(@ARGV);
 if ($test eq "-f") {$file=shift @ARGV}
 elsif ($test eq "-d") {$dest=shift @ARGV}
 elsif($test eq "-n") {$no=shift @ARGV}
 elsif($test eq "-p") {$blastp=shift @ARGV}
 elsif($test eq "-db") {$blastdb=shift @ARGV}
 elsif($test eq "-ex") {$blastexec=shift @ARGV}
 elsif($test eq "-N") {$jobname=shift @ARGV}
 elsif($test eq "-rev") {$rev=1}
 else {$blastopt .= " " . $test . " " . shift @ARGV}
}

die $usage unless ($file && $no);

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
my $required=ceil($counttags/$no); # ceil rounds up so we should get less files than if we use int or real rounding.
#print STDERR "There are $counttags sequences in $file and we are going to write $required per file\n";


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

@sourcefiles=reverse @sourcefiles if ($rev);

#print STDERR "Wrote to $filecount files in $dest\n";

unless ($blastp && $blastdb) {die "Can't blast because either program or database were not specified"}

my $name="bl$file";
$name=substr($name, 0, 10);

my $pwd=`pwd`; chomp($pwd);

my $submitted;
foreach my $sf (@sourcefiles) {
 my $command = $blastexec . " -p $blastp -d $blastdb -i $pwd/$dest/$sf -o $pwd/$dest/$sf.$blastp";
 $command .= " $blastopt";
 #print STDERR "$command\n";
 
 my $qsub= "qsub ";
 #$stdoutdir && ($qsub .= " -o $stdoutdir");
 #$stderrdir && ($qsub .= " -e $stderrdir");
 $qsub .= " -o $dest -e $dest";

print STDERR "Submitting $command from $dest\n";

 open(QSUB, "| $qsub") || die "Can't open a pipe to qsub";
 print QSUB <<EOF;
#!/bin/sh
#PBS -V
#PBS -l nodes=1:ppn=2:mpi
#PBS -N $name
#PBS -q workq
mpiexec $command
EOF

 close QSUB;
 
 $submitted++;
}

print STDERR "$submitted jobs submitted\n";


