#!/usr/bin/env perl

use strict;
use Bio::SeqIO;

my $file=shift;
my $dir=shift;
unless ($file && $dir) {die "$0 <genbank file> <dir to put separate files>"}

$dir =~ s/\/$//;
if (-e $dir) {die "$dir already exists"}
else {mkdir $dir, 0755}

my $si=Bio::SeqIO->new(-file=>$file, -format=>"genbank");
while (my $seq=$si->next_seq) {
 my $id=$seq->display_name;
 $id =~ s/\s+//g;
 my $so=Bio::SeqIO->new(-file=>">$dir/$id", -format=>"genbank");
 $so->write_seq($seq);
 $so->close;
}


