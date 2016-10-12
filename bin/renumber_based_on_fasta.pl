#!/usr/bin/perl -w

# renumber some sequences based on what is in the fasta file. We have renumbered the 454 sequences, and then lost them all in a harddrive failure
# but I have the fa (just not the qual files)

use lib '/clusterfs/home/rob/perl';
use Rob;

my $rob=new Rob;

my $usage=<<EOF;
$0 
-f new fasta file
-o original fasta file
-q original quality file
-x fasta ouptut file
-y quality output file

EOF

my ($fa, $of, $qf, $qo, $fo);
while (@ARGV)
{
 my $t=shift;
 if ($t eq "-f") {$fa=shift}
 elsif ($t eq "-o") {$of=shift}
 elsif ($t eq "-q") {$qf=shift}
 elsif ($t eq "-y") {$qo=shift}
 elsif ($t eq "-x") {$fo=shift}
}
die $usage unless ($fa && $of && $qf && $qo && $fo);
my $fasta=$rob->read_fasta($fa);
my $seq;
map {$fasta->{$_}=~s/\s+//g; $fasta->{$_}=uc($fasta->{$_}); $seq->{$fasta->{$_}}=$_} keys %$fasta;

my $orifa=$rob->read_fasta($of);
my $oriqf=$rob->read_fasta($qf, 1);


open(FA, ">$fo") || die "Can't open $fo for writing";
open(QU, ">$qo") || die "Can't open $qo for writing";

my %seen; my $skip; my $kept;
map 
{
 if ($seen{$orifa->{$_}}) {$skip++}
 else
 {
  $seen{$orifa->{$_}}=1;
  $orifa->{$_}=~s/\s+//g;
  $orifa->{$_}=uc($orifa->{$_});
  $kept++;
  if ($seq->{$orifa->{$_}})
  {
   print FA ">", $seq->{$orifa->{$_}}, "\n", $orifa->{$_}, "\n";
   print QU ">", $seq->{$orifa->{$_}}, "\n", $oriqf->{$_}, "\n";
  }
  else 
  {
   print STDERR "$_ not found!\n";
  }
 }
} keys %$orifa;

  
  
print STDERR "$kept sequences were kept (seen was: ", scalar(keys %seen), " and $skip sequences were ignored\n";
 
