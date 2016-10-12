#!/usr/bin/perl -w

use strict;

# convert the output from fig fig pegs_in_subsystems_by_homology to input for xipe. The output has three columns which are the peg id from the genome, the similar id, and the subsystem. We want to end up with a peg id, a subsystem, and a number, but all must just be numbers, so we need a key file too.

my $usage=<<EOF;

$0 
-k key file (probably ~/seed_database/ss_key.txt). This will be created or appended to with new ss
-i input file

This is expecting input in the form output by fig pegs_in_subsystems_by_homology

EOF

my ($key, $inf);
while (@ARGV) {
 my $t=shift @ARGV;
 if ($t eq "-k") {$key=shift @ARGV}
 elsif ($t eq "-i") {$inf=shift @ARGV}
}

die $usage unless ($key && $inf);


my %subsystem; my $max=1;
if (-e $key) {
 open(IN, $key) || die "Can't open $key";
 while (<IN>) {
  chomp;
  my @a=split /\t/;
  $a[0] = &clean_ss($a[0]);
  $subsystem{$a[0]}=$a[1];
  if ($a[1] > $max) {$max=$a[1]}
 }
}

my $found; my $peg;
my %new; my $rolecount; my $pegcount;
open(IN, $inf) || die "CAn't open $inf";
while (<IN>) {
 chomp;
 next unless ($_); # ignore blank lines
 next if (/^\s+$/);
 next if (/Looking for pegs/); # a comment line
 my @a=split /\t/;
 unless ($#a==2) {print STDERR "Not enough columns in $_\n"; next}
 $a[2]=&clean_ss($a[2]);
 next if ($found->{$a[2]}->{$a[0]}); # only count each peg once in each ss
 $found->{$a[2]}->{$a[0]}=++$rolecount;
 $peg->{$a[0]}=++$pegcount;
 unless ($subsystem{$a[2]}) {
  $subsystem{$a[2]}=++$max;
  $new{$a[2]}=$subsystem{$a[2]};
  print STDERR "Saving $a[2]\n";
 }
 print join("\t", $peg->{$a[0]}, $subsystem{$a[2]}, $found->{$a[2]}->{$a[0]}), "\n";
}

open(OUT, ">>$key") || die "CAn't append to $key";
print OUT map {"$_\t$new{$_}\n"} keys %new;
close OUT;

sub clean_ss {
 my $ss=shift;
 $ss =~ s/\_/ /g;
 $ss =~ s/  / /g;
 $ss =~ s/^\s+//;
 $ss =~ s/\s+$//;
 return $ss;
}
