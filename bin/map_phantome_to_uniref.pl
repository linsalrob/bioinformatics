#!/usr/bin/perl -w

use strict;

=pod

=head1

Map the phantome sequences to the UniRef100 sequences. For more information about this code, please see http://edwards.sdsu.edu/research/index.php/robs/348-mapping-uniref100-to-phantome

The code requires two arguments, the first is the fasta file from http://www.phantome.org/Downloads, the second is the fasta file from uniref. The code is written to handle zipped and gzipped files, so you may keep them compressed.

=cut

my ($phantomeF, $unirefF)=@ARGV;
die "$0 <phantome file> <uniref file>" unless (-e $phantomeF && -e $unirefF);

# first read the phantome file
my $fh = openfasta($phantomeF);
my %phantome;
my $id; my $seq=undef;
while (<$fh>) {
	chomp;
	if (s/^>//) {
		if ($seq) {
			$phantome{$seq}=$id;
			undef $seq;
		}
		s/\s+\[/\t/g; s/\]//g; # convert the [...] comments to tab separated text
		$id=$_;
	}
	else {
		$seq .= uc($_);
	}
}
$phantome{$seq}=$id; # don't forget the last sequence
close $fh; # close the filehandle so we can open the next one

($id, $seq)=(undef, undef); # clear these
$fh = openfasta($unirefF);
while (<$fh>) {
	chomp;
	if (s/^>//) {
		if ($seq && $phantome{$seq}) {&printout($id, $seq)}
		$id=$_;
		undef $seq;
	}
	else {
		$seq .= uc($_);
	}
}
if ($phantome{$seq}) {&printout($id, $seq)}
close $fh;
exit(0);







sub openfasta {
	my $file = shift;
	die "Fasta file not found" unless (-e $file);
	my $fh; # filehandle
		if ($file =~ /\.gz$/) {open($fh, "gunzip -c $file|") || die "Can't open a pipe to $file"}
	elsif ($file =~ /\.zip$/) {open($fh, "unzip -p $file|") || die "Can't open a pipe to $file"}
		else {open ($fh, $file) || die "Can't open $file"}
	return $fh;
}


sub printout {
	my ($id, $seq)=@_;
	$id =~ m/^(\S+)\s+(.*?)\s+n=\d+\s+Tax=(.*?)\s+RepID=(\S+)/;
	my ($uniid, $org, $repid)=($1, $2, $3);
	print join("\t", $uniid, $org, $repid, $phantome{$seq}, $seq), "\n";
}

