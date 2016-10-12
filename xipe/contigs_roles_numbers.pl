#!/usr/bin/perl -w 

use strict;

my $usage=<<EOF;
$0
-r roles file with numbers (see /clusterfs/home/rob/seed_database/numbered_ss.txt which is the default)
-x correspondence file between contig/sequence name and number (see C29G3_reads_headers.txt)
-c cutoff
-i ignore file. List of sequences to ignore, one per line
* blast files

EOF

my $cutoff=100; my $rolesf="/clusterfs/home/rob/seed_database/numbered_ss.txt"; my @files; my $xfile;
my $ignoref;
while (@ARGV) {
 my $t=shift @ARGV;
 if ($t eq "-r") {$rolesf=shift @ARGV}
 elsif ($t eq "-c") {$cutoff=shift @ARGV}
 elsif ($t eq "-x") {$xfile=shift @ARGV}
 elsif ($t eq "-i") {$ignoref=shift @ARGV}
 elsif (-e $t) {push @files, $t}
 else {print STDERR "Can't figure out what $t is on command line\n"}
}

my $seqcount=1; my %seqcount;

die $usage unless ($rolesf && $files[0]);
print STDERR "Counting roles from $rolesf with cutoff of $cutoff\n";

my %ignore;
if ($ignoref) {
 print STDERR "Ignoring from $ignoref\n";
 open(IN, $ignoref) || die "Can't open $ignoref";
 while (<IN>) {
  chomp;
  $ignore{$_}=1;
 }
}

my %contigname;
if ($xfile) {
 open(IN, $xfile) || die "Can't open $xfile";
 while(<IN>) {
  chomp;
  my @a=split /\t/;
  $a[1]=~s/\s.*//;
  $contigname{$a[1]}=$a[0];
 }
}

open(IN, $rolesf) || die "Can't open $rolesf";
my $role;
while (<IN>) {
 chomp;
 my @line=split /\t/;
 push @{$role->{$line[2]}}, \@line;
}
close IN;

my $contig;
foreach my $file (@files) {
 open (IN, $file) || die "Can't open $file";
 print STDERR "READING: $file\n";
 while (<IN>) {
  my @l=split /\t/;
  next unless ($l[10] <= $cutoff);
  next if ($ignore{$l[0]});
  if ($role->{$l[1]}) {
   # we want to keep this hit because it is significant
   foreach my $depth (@{$role->{$l[1]}}) {
    my ($ss, $fun, $peg)=@$depth;
    unless ($seqcount{$l[0]}) {$seqcount{$l[0]}=$seqcount++}
    my $id=$seqcount{$l[0]};
    unless ($id) {print STDERR "no id from $_\n"; next}
    if ($contigname{$l[0]}) {$id=$contigname{$l[0]}}
    $contig->{$id}->{$ss}->{$fun}++;
   }
  }
 }
}

foreach my $cont (sort {$b cmp $a} keys %$contig) {
 foreach my $ss (sort keys %{$contig->{$cont}}) {
  foreach my $role (sort keys %{$contig->{$cont}->{$ss}}) {
   print "$cont\t$ss\t$role\n";
  }
 }
}


 
   
	  
