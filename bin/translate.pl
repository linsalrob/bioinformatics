#!/usr/bin/perl -w
#

use strict;
use Rob;
my $rob=new Rob;
use raeseqlib qw/translate_seq/;

my $f=shift || die "fasta file?";
my $fa = $rob->read_fasta($f);
map {print ">$_\n", translate_seq($fa->{$_}), "\n"} keys %$fa;
