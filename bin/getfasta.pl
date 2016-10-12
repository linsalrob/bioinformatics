#!/usr/bin/perl -w

# get some fasta sequences from a file

use strict;
use lib '/clusterfs/home/rob/perl';
use Rob;

my $us=<<EOF;
$0

-f fasta file
-n needed sequences, one per line

EOF

my ($faf, $nf);

while (@ARGV) {
 my $t=shift;
 if ($t eq "-f") {$faf=shift}
 elsif ($t eq "-n") {$nf=shift}
}

die $us unless ($faf && $nf);

my $fasta = Rob->read_fasta($faf);
open(IN, $nf) || die "Can't opne $nf";
while (<IN>) {
 chomp;
 if ($fasta->{$_}) {
  $fasta->{$_} =~ s/\s+//g;
  print ">$_\n", $fasta->{$_}, "\n";
 } else {
  print STDERR "No sequence for $_\n";
 }
}



