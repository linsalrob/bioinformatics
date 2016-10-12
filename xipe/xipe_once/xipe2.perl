#!/usr/bin/perl -s
#!This program made by Pat McNairnie

use strict;

if(!@ARGV){
	print "Usage: ./xipe2.cgi FILENAME\n";
	exit();
}
my $filename  = $ARGV[0];
my $j;
my $k;

open(FD, "<$filename") or die("Couldn't open $filename\n") ;
my @arr;
while(<FD>){
	my @lines = split(/ /, $_);
	push(@arr,\@lines);
}
close(FD);

for ($j=0; $j<=$#{$arr[0]}; $j++) {
		my @sorter;
       for ($k=0; $k<=$#arr; $k++) {
			push(@sorter,$arr[$k][$j]);
       }
		@sorter = sort {$a <=> $b} (@sorter);
		print "@sorter\n";
}
