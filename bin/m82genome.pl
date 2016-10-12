#!/usr/bin/perl -w

# convert BLASTs from M8 format to a summary of what and counts

use strict;
use DBI;
use lib '/home/rob/perl';
use Bio::SeqIO;

my $dbh=DBI->connect("DBI:mysql:ncbi", "rob", "forestry") or die "Can't connect to database\n";

my $usage=<<EOF;

$0 <options>

-d directory of blast results in m8 (tab) format
-c cutoff (E value) below which to ignore results
-z count zeroes only
-s directory of sequences. If this is included the lengths of each sequence will be used

EOF


my ($dir, $seqdir,$cutoff,$countzero)=('', '', 1000000,0);
while (@ARGV) {
 my $test=shift(@ARGV);
 if ($test eq "-d") {$dir=shift @ARGV}
 elsif ($test eq "-s") {$seqdir=shift @ARGV}
 elsif ($test eq "-c") {$cutoff = shift @ARGV}
 elsif ($test eq "-z") {$countzero=1}
}

die $usage unless $dir;

&count_zero($dir) if $countzero;

my $lengths;
if ($seqdir) {$lengths=seq_dir($seqdir)}

opendir(DIR, $dir) || die "Can't open $dir";
while (my $file=readdir(DIR)) {
 next if ($file =~ /^\./);
 my $len='xx';
 if ($lengths->{$file}) {$len=$lengths->{$file}}
 if (-z "$dir/$file") # file has zero size
 {
  print "$file\t$len\t0\n";
  next;
 }
 my %orghits; # hits by organism

 open (IN, "$dir/$file") || die "Can't open $dir/$file";
 my %hits;
 my %evalue;
 my $sawhit;
 my %totalhitlen;
 while (<IN>) {
  my @line=split /\t/;
  next unless ($line[10] < $cutoff);
  $sawhit++;
  my $id;
  if ($line[1] =~ /gi\|(\d+)\|/) {$id=$1}
  else {print STDERR "Don't know how to handle $line[1]\n"; next}
  $hits{$id}++;
  $evalue{$id}=$line[10];
  $totalhitlen{$id}=$line[3];
 }
 close IN;
 
 unless ($sawhit) {
  # no significant hit
  print "$file\t$len\t0\n";
  print STDERR "Not a zero size but no sig hits less than $cutoff for $dir/$file\n";
  next;
 }

 foreach my $gi (keys %hits) {
  my @row=$dbh->selectrow_array("select tax_id from gi_taxid_nucl where gi = $gi");
  my $tax=$row[0];
  unless ($tax) {
   print STDERR "No tax gi for $gi\n";
   print "$file\t$len\t$hits{$gi}\t$gi\t$evalue{$gi}\t\t$totalhitlen{$gi}\n";
   next;
  }
  my $exc=$dbh->prepare("select * from names where tax_id = $tax");
  $exc->execute || die $dbh->errstr;
  my $name=''; my $badname='';
  while (my @res=$exc->fetchrow_array) {
   if ($res[4] eq "scientific name") {$name=$res[2]}
   else {$badname=$res[2]}
  }
  if (!$name && $badname) {$name=$badname}
  unless ($name) {
   print STDERR "Couldn't parse tax id $tax for gi $gi\n";
   print "$file\t$len\t$hits{$gi}\t$gi\t$evalue{$gi}\t\t$totalhitlen{$gi}\n";
  }
  push @{$orghits{$name}}, $gi;
  unless (exists $evalue{$name}) {$evalue{$name}=$evalue{$gi}}
  if ($evalue{$name} < $evalue{$gi}) {$evalue{$name}=$evalue{$gi}}
  if ($gi == 45774915) {
   print STDERR "For 45774915 ($name) Had $totalhitlen{$name} Now added $totalhitlen{$gi}\n";
  }
  $totalhitlen{$name}+=$totalhitlen{$gi};
  if ($gi == 45774915) {
   print STDERR "For 45774915 ($name) Had something, now have $totalhitlen{$name}\n";
  }
 }
 foreach my $name (keys %orghits) {
  print  "$file\t$len\t", scalar(@{$orghits{$name}}), "\t", (join ", ", (@{$orghits{$name}})), "\t", $evalue{$name}, "\t$name\t$totalhitlen{$name}\n";
 }
}

  


sub count_zero {
 my $dir=shift;
 my $count=0;
 opendir(DIR, $dir) || die "Can't open %dir";
 while (my $file=readdir(DIR)) {
  next if ($file =~ /^\./);
  $count++ if (-z "$dir/$file");
 }
 print "$count files had zero size\n";
 exit(0);
}

sub seq_dir {
 my $seqdir=shift;
 opendir(DIR, $seqdir) || die "Can't open $seqdir";
 my $length;
 while (my $file = readdir(DIR)) {
  next if ($file =~ /^\./);
  my $sin=Bio::SeqIO->new(-file=>"$seqdir/$file", -format=>"fasta");
  while (my $seq=$sin->next_seq) {
   $length->{$file}=$seq->length();
   }
 }
 print STDERR "Read ", scalar(keys %$length), " lengths\n";
 return $length;
}
