#!/usr/bin/perl -w 

use strict;

my $usage=<<EOF;
$0
-i ignore file. List of sequences to ignore, one per line
-r roles file (ss roles pegs)
-c cutoff
* blast files

EOF

my $cutoff=100; my $rolesf; my @files; my $ignoref;
while (@ARGV) {
 my $t=shift @ARGV;
 if ($t eq "-r") {$rolesf=shift @ARGV}
 elsif ($t eq "-c") {$cutoff=shift @ARGV}
 elsif ($t eq "-i") {$ignoref=shift @ARGV}
 elsif (-e $t) {push @files, $t}
 else {print STDERR "Can't figure out what $t is on command line\n"}
}

die $usage unless ($rolesf && $files[0]);

my %ignore;
if ($ignoref) {
 print STDERR "Ignoring from $ignoref\n";
 open(IN, $ignoref) || die "Can't open $ignoref";
 while (<IN>) {
  chomp;
  $ignore{$_}=1;
 }
}

print STDERR "Counting roles from $rolesf with cutoff of $cutoff\n";
open(IN, $rolesf) || die "Can't open $rolesf";
my $role;
while (<IN>) {
 chomp;
 my @line=split /\t/;
 if ($#line==4) {my $t=-shift @line; $t=shift @line} # get rid of classification if we have it
 push @{$role->{$line[2]}}, \@line;
}
close IN;

my $contig;
foreach my $file (@files) {
 open (IN, $file) || die "Can't open $file";
 while (<IN>) {
  my @l=split /\t/;
  next unless ($l[10] <= $cutoff);
  next if ($ignore{$l[0]});
  if ($role->{$l[1]}) {
   # we want to keep this hit because it is significant
   foreach my $depth (@{$role->{$l[1]}}) {
    my ($ss, $fun, $peg)=@$depth;
    $contig->{$l[0]}->{$ss}->{$fun}++;
   }
  }
 }
}

foreach my $cont (sort {$b cmp $a} keys %$contig) {
 foreach my $ss (sort keys %{$contig->{$cont}}) {
  foreach my $role (sort keys %{$contig->{$cont}->{$ss}}) {
   print "$cont\t$ss\t$role\t", $contig->{$cont}->{$ss}->{$role}, "\n";
  }
 }
}


 
   
	  
