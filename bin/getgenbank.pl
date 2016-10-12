#!/usr/bin/perl -w

# get genomes.
# first get accession numbers

use strict;
use lib '~/bioperl/1.3/bioperl-1.302';
use Bio::DB::GenBank;
use Bio::SeqIO;
my $gb = new Bio::DB::GenBank;
my $do = shift || die "$0 <-s -b> <file>\n\t-s for single, -b for batch\n";
my $file=shift || die "$0 <-s -b> <file>\n\t-s for single, -b for batch\n";

my @accs;
open (IN, $file) || die "$file";
while (<IN>) {
 chomp;
 next unless ($_);
 push (@accs, $_);
 }
close IN;

print STDERR $#accs+1, " sequences to get\n";


if ($do eq "-s") {singleget(\@accs)}
else {batchget(\@accs)}
exit(0);


sub singleget {
 my $accs = shift;
# my $seqout = new Bio::SeqIO(-fh => \*STDOUT, -format => 'genbank');
 foreach my $acc (@$accs) {
  print STDERR "Trying $acc\n";
  my $seqout = new Bio::SeqIO(-file => ">$acc", -format => 'genbank');
  my $seq = $gb->get_Seq_by_acc($acc);
  if ($seq) {$seqout->write_seq($seq)}
  else {print STDERR "ERROR with $acc\n"}
 }
}

sub batchget {
 my $accs=shift;
 my $seqout = new Bio::SeqIO(-fh => \*STDOUT, -format => 'genbank');
 my $seqio = $gb->get_Stream_by_batch($accs);

 my $seq;
 while( defined ($seq = $seqio->next_seq )) {
  if ($seq) {$seqout->write_seq($seq)}
  else {print STDERR "ERROR with a sequence?\n"}
 }
}
