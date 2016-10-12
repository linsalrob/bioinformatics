#!/usr/bin/perl -w
#
=pod

A program to convert a taxonomy id to the name at GenBank using eutils

=cut

use strict;
use LWP::Simple;

my $id=shift || die "Taxonomy id?";
my $oriid = $id;
if ($id =~ s/\.(\d+)$//) {
	print STDERR "The version number $1 was removed from the tax id\n";
}

my  $url = 'http://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=taxonomy&id=' . $id;
my $results = LWP::Simple::get($url);

$results =~ /<ScientificName>(.*)<\/ScientificName>/;
#print "$id\t$1\n";
my $sn=$1;
unless ($sn) {print STDERR "No scientific name for $id\n"}

$results =~ /<Lineage>(.*)<\/Lineage>/s;
my $ln = $1;
print join("\t", $oriid, $id, $sn, $ln), "\n";


