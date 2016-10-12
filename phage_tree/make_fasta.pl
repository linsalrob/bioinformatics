#!/usr/bin/env perl



use strict;
use Phage;
my $phage=new Phage;
use FIG;
my $fig=new FIG;

open(FASTA, ">phage_proteins.fasta") || die "can't open phage_proteins.fasta";

foreach my $genome ($phage->phages()) {
	foreach my $peg ($fig->pegs_of($genome)) {
		my $trans = $fig->get_translation($peg);
		print FASTA ">$peg\n$trans\n";
	}
}

		
