#!/usr/bin/perl -w

# a quick gbk2fasta script that will put locus/definition on the ID line and DNA sequence

use strict;

my ($print, $locus, $def);

while (<>)
{
 if (/^LOCUS\s+(\S+)/) {$locus=$1; next}
 if (/^DEFINITION\s+(.*)/) {$def=$1; chomp($def)}
 if (/ORIGIN/)
 {
  if ($locus && $def) 
  {
   print ">$locus $def\n";
   $print=1;
   next;
  }
  else
  {
   print STDERR "No locus or definition we have |$locus| and |$def|\n";
   $print=2;
  }
  next;
 }
 if (m#^//#) {($locus, $def, $print)=(undef, undef, undef)}
 if ($print && $print ==2)
 {
  print STDERR "Next line is:\n$_";
  undef $print;
  next;
 }
 if ($print) 
 {
  chomp;
  s/\d//g;
  s/\s//g;
  print "$_\n";
 }
} 
 
