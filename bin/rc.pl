#!/usr/bin/perl -w
#

use strict;
use Rob;
my $rob=new Rob;

my $f=shift || die "fasta file?";
my $fa = $rob->read_fasta($f);
foreach my $id (keys %$fa) {
	print ">$id\n", $rob->rc($fa->{$id}), "\n";
}
