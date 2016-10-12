#!/usr/bin/perl -w
#

## Split sequences based on the tags in the sequences, and optionally remove those tags. 
# 
# You can specify a file of known tags with the -t option. Otherwise, we'll use the default 96 well Ion set.
# You can specify NOT removing the tag with the -k flag
# 
# If you provide a fastq file we will also remove the quality scores too.


use strict;
use Rob;
my $rob=new Rob;

unless (@ARGV) {
	die "$0 [-k keep tags] -[t file with barcode list] [-m minimum sequence length to keep] <list of fasta files>";
}

## check our options
my $barcodefile = "/home/redwards/bioinformatics/ion/barcodeList.txt";
my @fastafiles;
my @fastqfiles;
my $keeptag;
my $minimumLength = 0;
while (@ARGV) {
	my $arg=shift @ARGV;
	if ($arg eq "-t") {$barcodefile = shift @ARGV}
	elsif ($arg eq "-m") {$minimumLength = shift @ARGV}
	elsif ($arg eq "-k") {$keeptag = 1}
	elsif ($arg =~ /fna$/ || $arg =~ /fa$/) {push @fastafiles, $arg}
	elsif ($arg =~ /fq$/ || $arg =~ /fastq$/) {push @fastqfiles, $arg}
	else {print STDERR "Srry, we don't understand this argument: $arg\n";
	}
}

unless (@fastafiles || @fastqfiles) {die "You must provide some sequence files"}


my %bc; my %len;


open(IN, $barcodefile) || die "can't open $barcodefile";
while (<IN>) {
	next unless (/^barcode/);
	chomp;
	my ($number, $ie, $seq, $key, undef, undef, $length, undef)=split /\,/;
	unless ($length == length($seq)) {
		print STDERR "Expected length of $number ($length) is not the same as the length of the sequence. Are you sure that column is the length? (I just guestimated it was...)\n";
	}
	$bc{uc($seq)}=$ie;
	$len{uc($seq)}=$length;
}
close IN;


if (@fastafiles) {
	open(MISMATCH, ">not_mapped.fna") || die "Can't write to not_mapped";
	if ($minimumLength > 0) {
		open(MINIMUM, ">short_sequences.fna") || die "Can't write to short sequences";
	}

	my $fhs;
	foreach my $faf (@fastafiles) {
		print STDERR "$faf\n";
		my $fa = $rob->read_fasta($faf);
		my $did=0; my $didnot=0;
		foreach my $id (keys %$fa) {
			my $seq = $fa->{$id};
			if ($minimumLength > 0 && length($seq)+10 < $minimumLength) {
				print MINIMUM ">$id\n$seq\n";
				next;
			}
			my $matched;
			foreach my $bar (sort {$len{$b} <=> $len{$a}} keys %bc) {
				if (index(uc($seq), $bar) == 0) {
					if (!$fhs->{$bar}) {
						my $fh;
						open($fh, ">$bc{$bar}.fna") || die "Can't open $bc{$bar}.fna";
						$fhs->{$bar}=$fh;
					}
					my $fh = $fhs->{$bar};
					unless ($keeptag) {$seq =~ s/^$bar//}

					print $fh ">$id\n$seq\n";
					$matched=1;
					$did++;
					last;
				}
			}
			unless (defined $matched) {print MISMATCH ">$id\n$seq\n"; $didnot++}
		}
		print "$faf\t$did\t$didnot\n";
	}
}

if (@fastqfiles) {
	open(MISMATCH, ">not_mapped.fastq");
	if ($minimumLength > 0) {
		open(MINIMUM, ">short_sequences.fna") || die "Can't write to short sequences";
	}

	my $fhs;
	foreach my $fqf (@fastqfiles) {
		my $fq = $rob->read_fastq($fqf);
		my $did=0; my $didnot=0;
		foreach my $id (keys %$fq) {
			my $seq = $fq->{$id}->[0];
			my $qual = $fq->{$id}->[1];
			if ($minimumLength > 0 && length($seq)+10 < $minimumLength) {
				print MINIMUM "\@$id\n$seq\n+$id\n$qual\n";
				next;
			}
			my $matched;
			foreach my $bar (sort {$len{$b} <=> $len{$a}} keys %bc) {
				if (index(uc($seq), $bar) == 0) {
					if (!$fhs->{$bar}) {
						my $fh;
						open($fh, ">$bc{$bar}.fastq") || die "Can't open $bc{$bar}.fastq";
						$fhs->{$bar}=$fh;
					}
					my $fh = $fhs->{$bar};
					unless ($keeptag) {
						$seq = substr($fq->{$id}->[0], $len{$bar});
						$qual = substr($fq->{$id}->[0], $len{$bar});
					}
					if (length($seq) == 0) {
						$minimumLength && print MINIMUM "\@$id\n$seq\n+$id\n$qual\tAfter trimming $bar\n";
						next;
					}

					print $fh "\@$id\n$seq\n+$id\n$qual\n";
					$matched=1;
					$did++;
					last;
				}
			}
			unless (defined $matched) {print MISMATCH "\@$id\n$seq\n+$id\n$qual\n"; $didnot++}
		}
		print "$fqf\tmatched: $did (", (($did/($did+$didnot))*100), "%)\t not matched: $didnot\n";
	}
}






