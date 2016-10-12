#!/usr/bin/perl -w

# randomize a fasta file (and optionally a quality file)

use strict;
use lib '/clusterfs//home/rob/perl';
use Rob;

my $usage=<<EOF;
$0 
-f fasta file
-q quality file (optional)
-o output file (default = fasta file.rand)
EOF

my ($faf, $qaf, $oaf);
while (@ARGV) {
 my $t=shift;
 if ($t eq "-f") {$faf=shift}
 if ($t eq "-q") {$qaf=shift}
 if ($t eq "-o") {$oaf=shift}
}

die $usage unless ($faf);
$oaf = $faf.".rand" unless ($oaf);

my $fa=Rob->read_fasta($faf);
my $qual;
if ($qaf) {$qual=Rob->read_fasta($qaf, 1)} else {print STDERR "No quality file defined, just randomizing $faf\n"}


open(FA, ">$oaf") || die "Can't open $oaf for writing\n";
if ($qaf) {open(QU, ">$oaf.qual") || die "Can't open $oaf.qual for writing\n"}
my $allkeys=Rob->rand([keys %$fa]);
print STDERR "Writing ", scalar(@$allkeys), " sequences\n";
foreach my $key (@$allkeys)
{
 $fa->{$key}=~ s/\s+//g;
 print FA ">$key\n",$fa->{$key},"\n";
 if ($qaf) {
  if ($qual->{$key}) {
   print QU ">$key\n", $qual->{$key}, "\n";
  }
 }
}

print STDERR "DOne in ", $^T-time, " seconds\n";
 
 
