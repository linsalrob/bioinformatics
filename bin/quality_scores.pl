#!/usr/bin/perl -w

# count the quality scores for the files from Pedulla

use strict;
use lib '/home/rob/perl';
use Rob;
my $rob=new Rob;

my $usage=<<EOF;
$0 <options>

-f file
-c cutoff (default=20)

EOF

my ($file, $cutoff)=('', 20);
while (@ARGV) {
 my $test=shift;
 if ($test eq "-f") {$file=shift @ARGV}
 elsif ($test eq "-c") {$cutoff=shift @ARGV}
}

die $usage unless $file;

my $fasta=$rob->read_fasta($file);

print join "\t", "Sequence", "Number of bases", "Mean", "Median", "Standard Deviation", "Longest run of sequences greater than $cutoff", "\n"; 
foreach my $key (keys %$fasta) {
 my $max=0;
 my $current=0;
 my @values=split /\s+/, $fasta->{$key};
 #print STDERR join "|", $key, @values, "\n";
 my @temp;
 # this removes undefs from the array
 foreach (@values) {if (/\d/) {push @temp, $_}}
 @values=@temp;
 foreach my $val (@values) {
  if ($val > 99) {
   print STDERR "Problem with $val\n";
  }
  if ($val < $cutoff) {
   if ($current > $max) {$max=$current}
   $current=0;
  } 
  else {
   $current++;
  }
 }

 my $label=$key;
 $label =~ s/\s*PHD.*//;
 print join "\t", $label, scalar(@values), $rob->mean(\@values), $rob->median(\@values), $rob->stdev(\@values), $max, "\n";
}


