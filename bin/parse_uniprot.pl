#__perl__

use strict;
use Digest::MD5;
use Data::Dumper;

# a simple uniprot parser because the parser in the version of bioperl  I have is broken.
# this is extremely simple, we are just printing out the id, accession number, go terms, md5 sum (so we can map to seed ids) and protein sequence

my $f = shift || die "Swiss Prot file?";
my $record;
my $inseq = 0;
# open the file, regardless of format
if ($f =~ /gz$/) {
	open(IN, "gunzip -c $f|") || die "Can't open a pipe to gunzip $f";
} elsif ($f =~ /zip$/) {
	open(IN, "unzip -p $f|") || die "can't open a pipe to unzip $f";
}	else {
	open(IN, $f) || die "can't open $f";
}

while (<IN>) {
	chomp;
	if (/^\/\/$/) {
		if ($record) {&output($record)};
		undef $record;
		$inseq=0;
		next;
	}
	if ($inseq) {
		s/^\s+//;
		$record->{'SQ'} .= $_;
		next;
	}
	m/^(\w\w)\s+(.*?)$/;
	my ($id, $data) = ($1, $2);
	if ($id eq "SQ") {$inseq = 1; next}
	push @{$record->{$id}}, $data;
}


sub output {
	my $rec = shift;
	#die Dumper($rec);

	my $trans = $rec->{'SQ'};
	$trans =~ s/\s//g;
	my $md5 = Digest::MD5::md5_hex(uc($trans));
	my $go = join("; ", grep {m/^GO/} @{$rec->{'DR'}});
	$rec->{'ID'}->[0] =~ m/^\s*(\S+)\s/;
	my $id = $1;

	foreach my $ac (map {split /\;\s*/, $_} @{$rec->{'AC'}}) {
		print join("\t", $id, $ac, $go, $md5, $trans), "\n";
	}
}

