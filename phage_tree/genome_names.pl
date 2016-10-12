#!/usr/bin/env perl 
#
# get all the proteins that belong to phages as a single fasta file

use strict;
use FIG;
my $fig=new FIG;
use Phage;
my $phage=new Phage;



foreach my $g ($phage->phages()) {
	my $n = my $m = $fig->genus_species($g);
	$n .= " ";
	$n =~ s/, complete genome\./i/;
	$n =~ s/ complete genome\./i/;
	$n =~ s/\s+genome\s+/ /; 
	$n =~ s/^unclassified\s+//i;
	$n =~ s/unclassified\.$//i;
	$n =~ s/^\S+viridae\s+//;
	$n =~ s/^\S+virus\s+//;
	$n =~ s/^\S+virales\s+//;
	$n =~ s/^\S+virinae\s+//;
	$n =~ s/^\S+\-like\s+//;
	$n =~ s/^Viruses\s+//i; $n =~ s/\s+Viruses\.\s+//i;
	$n =~ s/\s+prophage\s+/ /i;
	$n =~ s/\s+bacteriophage\s+/ /i;
	$n =~ s/\s+phage\s+/ /i;
	$n =~ s/\s+temperate\s+/ /i;
	$n =~ s/\s+complex\s+/ /i;
	$n =~ s/\s+filamentous virus\s+/ fv /;
	$n =~ s/\s+spindle-shaped virus\s+/ ssv /;
	$n =~ s/^\s*phages\s+//;
	$n =~ s/\s+phage\s+/ /;
	$n =~ s/\s+sensu lato\s+/ /;
	$n =~ s/\s+virus\s+/ /;
	$n =~ s/\s+$//;
	$n =~ s/^(\S{3})\S*\s+/$1\. /;
	print join("\t", $g, $m, $n), "\n"; 
}


