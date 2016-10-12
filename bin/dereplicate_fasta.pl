#!/usr/bin/perl -w

# remove duplicate names and sequences

use lib '/home/rob/perl';
use strict;
use Bio::SeqIO;
use Getopt::Long;
my ($infile, $outfile);
GetOptions(
 'i|infile:s'   	=> \$infile,
 'o|outfile:s'		=> \$outfile,
);

die "$0 -i infile -o outfile" unless ($infile && $outfile);
my $sio=Bio::SeqIO->new(-file=>"$infile", -format=>'fasta');
my $sout=Bio::SeqIO->new(-file=>">$outfile", -format=>'fasta');

my %seq; my %id;
while (my $sin=$sio->next_seq) {
 next if ($seq{uc($sin->seq)});
 $seq{uc($sin->seq)}=1;
 if ($id{$sin->id}) {
  $id{$sin->id}++;
  $sin->id($sin->id . "." . $id{$sin->id});
 }
 else {$id{$sin->id}=1}
 
 $sout->write_seq($sin);
}

