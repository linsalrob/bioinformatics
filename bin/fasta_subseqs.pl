#!/usr/bin/perl -w

# grab some sequences from a file, at random

use strict;
use lib '/clusterfs/home/rob/perl';
use Rob;
my $rob=new Rob;

my $fasta=shift||die "$0 <fasta file><number of seqs>\n";
my $nos=shift||die "$0 <fasta file><number of seqs>\n"; 

my $fa=$rob->read_fasta($fasta);
my $arr=$rob->rand([keys %$fa]);
map {
 $fa->{$_}=~s/\s+//g; $fa->{$_}=~ s/(.{60})/$1\n/g; chomp($fa->{$_});
 print ">$_\n", $fa->{$_}, "\n";
} splice(@$arr, 0, $nos);

