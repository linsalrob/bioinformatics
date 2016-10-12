#!/usr/bin/perl -w

# create a fasta file from some directory

use strict;

my $usage=<<EOF;
$0
-d list of source directories
-e list of directories to exclude
-o output file

EOF

my (@inc, @exc, $output);
my $exc;
while (@ARGV) {
 my $t=shift @ARGV;
 if ($t eq "-d") {$exc=1;next}
 elsif ($t eq "-e") {$exc=2;next}
 elsif ($t eq "-o") {$output=shift @ARGV; next}
 if ($exc==1) {push @inc, $t}
 elsif ($exc==2) {push @exc, $t}
}

die $usage unless ($inc[0] && $output);

my $all=0; my $kept=0;
foreach my $dir (@inc) {
 opendir(DIR, $dir) || die "Can't open $dir";
 while (my $f=readdir(DIR)) {
  next if ($f =~ /^\./);
  $all++;
  my $keep=1; 
  foreach my $e (@exc) {
   if (-e "$e/$f") {undef $keep; last}
  }
  next unless ($keep);
  `cat $dir/$f >> $output`;
  $kept++;
 }
}

print STDERR "Looked through $all files and kept $kept of them\n";

 

