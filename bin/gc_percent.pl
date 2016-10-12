#!/usr/bin/perl -w

use strict;
use lib $ENV{HOME}."/bioinformatics/Modules";
use Rob;
my $rob=new Rob;

my $file=shift || die "fasta file?";
my $fa = $rob->read_fasta($file);

my ($gc, $at)=(0,0);
foreach my $i (keys %$fa)
{
	my $seq=$fa->{$i};

	while ($seq)
	{
		my $base = chop($seq);
		(lc($base) eq "c" || lc($base) eq "g") ? ($gc++) : ($at++);
	}
}
print "GC: $gc AT: $at Percent = ", int(($gc/($at+$gc))*10000)/100, "\n";
