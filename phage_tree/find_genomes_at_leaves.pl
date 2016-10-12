#!/usr/bin/perl -w

use strict;
use lib '/home/redwards/bioinformatics/Modules';
use ParseTree;
use Getopt::Std;
my %opts;
getopts('f:p:', \%opts);
unless ($opts{f} && $opts{p}) {
	die <<EOF;
$0
-f tree file
-p phage_names.txt file

Print out the genome ids of all the genomes on the leaves.

EOF
}


my $pt = new ParseTree();

my $treein;

open(IN, $opts{f}) || die "Can't open $opts{f}";
while (<IN>) {chomp; $treein.=$_}
close IN;
my $tree = $pt->parse($treein);

my %nodes;
nodes($tree);

$nodes{"Ent_G4_clade_Enterobacteria_phage_ID2_Moscow/ID/2001"} = "Ent_G4_ID2_";
$nodes{"Str_mu1/6"}="Str_mu16";
$nodes{"Bac_Gamma_isolate_d'Herelle"}="Bac_Gamma_isolate_d_Herelle";
$nodes{"Bur_cepacia_complex_BcepC6B"}="Bur_cepacia_BcepC6B";

my %seen;

open(IN, $opts{p}) || die "can;t open $opts{p}";
while (<IN>) {
	chomp;
	my ($id, $name, $abbr)=split /\t/;
	$abbr =~ s/\-//g;
	$abbr =~ s/\.//g;
	$abbr =~ s/\ /\_/g;
	if ($nodes{$abbr}) {
		$seen{$abbr}=1;
		$seen{$nodes{$abbr}}=1;
		print join("\t", $id, $name, $nodes{$abbr}), "\n";
	} else {
		print STDERR "NO NODE FOR $abbr\n";
	}
}

print STDERR "MISSED:\n";
map {print STDERR $_,  "\n" unless ($seen{$_})} keys %nodes;




sub nodes {
	my ($node)=@_;
	if (ref($node) eq "ARRAY") {
		nodes($node->[0]);
		nodes($node->[1]);
	}
	else {
		$nodes{$node}=$node; 
	}
}



