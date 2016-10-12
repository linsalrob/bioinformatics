#!/usr/bin/perl -w
#

=pod

Extract locus, definition, host, etc from phages, and then also get DNA md5 sum

=cut


use strict;
use Bio::SeqIO;

my $f=shift || die "Genbank file to parse";


my $sin=Bio::SeqIO->new(-file=>$f);
while (my $seq=$sin->next_seq()) {
	my $host;
	foreach my $feat ($seq->get_SeqFeatures()) {
		for my $ann ($feat->annotation()) {
			my %keys = map {($_=>1)} $ann->get_all_annotation_keys();
			if ($keys{'host'}) {
				my @ann = $ann->get_Annotations('host');
				$host = join(" :: ", @ann);
			}
		}
	}
	
	if (defined $host) {print join("\t", $seq->display_id(), $seq->desc(), $host), "\n"}
}



