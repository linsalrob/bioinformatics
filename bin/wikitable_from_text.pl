#!/usr/bin/perl -w
#

use strict;
use Getopt::Std;

my %opts;
getopts('hf:', \%opts);
unless ($opts{f}) {die "$0 \n-f text file \n-h if you have a header line\n"}
open(IN, $opts{f}) || die "Can't open $opts{f}";

print "{| class=\"wikitable\"  border=\"1\" cellpadding=\"2\"\n";

my $first=1;
while (<IN>)
{
	chomp;
	my @a=split /\t/;
	if ($opts{h}) {
		undef $opts{h};
		undef $first;
		map {print "! $_\n"} @a;
		next;
	}
	if ($first) {
		undef $first;
		map {print "| $_\n"} @a;
		next;
	}
	print "|-\n";
	map {print "| $_\n"} @a;
}
print "|}\n";


