#!/usr/bin/perl -w
#

use strict;
use Rob;

my $f=shift || die "text file to clean?";
open(IN, $f) || die "Can't open $f";
while (<IN>) {print Rob->ascii_clean($_)}
close IN;
