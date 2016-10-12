#!/usr/bin/perl -w

# count jobs in the queue at specified intervals 

use strict;
$|=1;
my $time=shift || die "$0 <time delay in seconds>";

while (1) 
{
	my $in=`qstat | wc -l`;
	chomp($in);
	$in-=2;
	print scalar(localtime(time)), "\t$in\n";
	sleep($time);
}

