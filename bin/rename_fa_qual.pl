#!/usr/bin/perl -w

# rename the sequences in a fasta file and a quality file begining with 1.

use strict;
use lib '/clusterfs/home/rob/perl';
use Rob;

my $usage=<<EOF;
$0 <options>
-f fasta file
-q quality file (can be ommitted if it is -f.qual)
-s start at (default=1)
-d destination fasta file (default = -f.renamed). Quality file is named from this 

EOF

my ($faf, $qaf, $daf, $start)=('','','',1);
while (@ARGV) {
 my $t=shift;
 if ($t eq "-f") {$faf=shift}
 elsif ($t eq "-q") {$qaf = shift}
 elsif ($t eq "-s") {$start = shift}
 elsif ($t eq "-d") {$daf=shift}
}

die $usage unless ($faf);
$qaf = $faf.".qual" unless ($qaf);
$daf = $faf.".renamed" unless ($daf);

die "Fasta file $faf not found" unless (-e $faf);
die "Quality file $qaf not found" unless (-e $qaf);


my $seqs=Rob->read_fasta($faf);
my $qual=Rob->read_fasta($qaf, 1);

open(FA, ">$daf") || die "Can't open $daf for writing";
open(QU, ">$daf.qual") || die "Can't open $daf.qual for writing";


foreach my $id (keys %$seqs) { 
 if ($qual->{$id}) {
  $seqs->{$id}=~s/\s+//g;
  print FA ">$start\n", $seqs->{$id}, "\n";
  print QU ">$start\n", $qual->{$id}, "\n";
  $start++;
 }
 else {
  print STDERR "No qualitites for $id. Skipped\n";
 }
}

foreach my $id (keys %$qual) {
 print STDERR "No fasta for quality $id\n" unless ($seqs->{$id});
}
 
print STDERR "$start sequences were renamed\n";

