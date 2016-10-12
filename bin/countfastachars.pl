#!/usr/bin/perl -w

# count fasta characters

use strict;

my $rewrite; # should we rewrite the file?
my $summary; #do we just want a summary
my $protein; # is it a protein
my $minsize=0; # minimum contig size to consider
unless (@ARGV) {die "$0: <-r rewrite file> <-m minimum contig size to include> <-s summary> <-p protein files> <fasta files or sequences to count>\n"}
while (@ARGV) {
	my $file = shift @ARGV;
	if ($file eq "-r") {$rewrite=1; next}
	if ($file eq "-s") {$summary=1; next}
	if ($file eq "-p") {$protein=1; next}
	if ($file eq "-m") {$minsize=shift @ARGV; next}
	if (-e $file) {
		my $runningtotal; my $number;
		my %long; my %short; $long{'length'}=1; $short{'length'}=1000000000; # use these to store shortest and longest
		if ($file =~ /\.gz$/) {open(IN, "gunzip -c $file|") || die "Can't open pipe to gunzip"}
		else {open (IN, $file)|| die "Can't open $file\n"}
		unless ($summary) {
			open (OUT, ">$file.count") || die "Can't open $file.count for writing\n";
			if ($protein) {print OUT "Sequence\tLength\tRunning Total\n"}
			else {print OUT "Sequence\tA\tG\tC\tT\tN\tLine total\tRunning total\n"}
		}
		if ($rewrite) {open REWRITE, ">$file.fasta" || die "Can't open $file.fasta for writing\n"}
		my $seq; my $tag; my @lengths;
		while (<IN>) {
			chomp;
			if (/^>/) {
				if ($seq) {
					my $length=length($seq);
					if ($length >= $minsize) {
						push @lengths, $length;
						if ($length > $long{'length'}) {%long=('length'=>$length, 'tag'=>$tag)}
						if ($length < $short{'length'}) {%short=('length'=>$length, 'tag'=>$tag)}
						if ($protein) {
							$runningtotal+=length($seq); $number++;
							unless ($summary) {print OUT "$tag\t", length($seq), "\t", $runningtotal, "\n"}
							undef $seq;
						}
						else {
							my $a = $seq =~ s/a/A/ig || '0';
							my $g = $seq =~ s/g/G/ig || '0';
							my $c = $seq =~ s/c/C/ig || '0';
							my $t = $seq =~ s/t/T/ig || '0';
							my $n = $seq =~ s/n/N/ig || '0';
							my $total = $a+$g+$c+$t+$n;
							$runningtotal+=$total; $number++;
							undef $seq;
							unless ($total == $length) {print STDERR "WARNING: Sequence contains letters that are not A, G, C, T, or N\n"}
							if ($rewrite) {$seq =~ s/(.{0,60})/$1\n/g; chomp($seq); print REWRITE "$tag\n$seq"}
							unless ($summary) {print OUT "$tag\t$a\t$g\t$c\t$t\t$n\t$total\t$runningtotal\n"}
						}
					}
				}
				$tag=$_;
			}
			else {s/\s//g; $seq.=$_}
		}
		my $length=length($seq);
		if ($length >= $minsize) {
			push @lengths, $length;
			if ($length > $long{'length'}) {%long=('length'=>$length, 'tag'=>$tag)}
			if ($length < $short{'length'}) {%short=('length'=>$length, 'tag'=>$tag)}
			if ($protein) {
				$runningtotal+=length($seq); $number++;
				unless ($summary) {print OUT "$tag\t", length($seq), "\t$runningtotal\n"}
			}
			else {
				my $a = $seq =~ s/a/A/ig || '0';
				my $g = $seq =~ s/g/G/ig || '0';
				my $c = $seq =~ s/c/C/ig || '0';
				my $t = $seq =~ s/t/T/ig || '0';
				my $n = $seq =~ s/n/N/ig || '0';
				my $total = $a+$g+$c+$t+$n;
				undef $seq;
				$runningtotal+=$total; $number++;
				unless ($total == $length) {print STDERR "WARNING: Sequence contains letters that are not A, G, C, T, or N\n"}
				if ($rewrite) {$seq =~ s/(.{0,60})/$1\n/g; chomp($seq); print REWRITE "$tag\n$seq"}
				unless ($summary) {print OUT "$tag\t$a\t$g\t$c\t$t\t$n\t$total\t$runningtotal\n"}
			}
		}
		if ($summary) {
			print "$file\n\tnumber of seqs:\t$number\n\ttotal:\t\t$runningtotal\n\taverage\t\t", $runningtotal/$number, "\n\tN50:\t\t", &N50(\@lengths, $runningtotal), "\n\tshortest:\t$short{'tag'} ($short{'length'})\n\tLongest:\t$long{'tag'} ($long{'length'})\n"}
	}
	else {
		my $seq = $file;
		my $length=length($seq);
		my $a = $seq =~ s/a/A/ig || '0';
		my $g = $seq =~ s/g/G/ig || '0';
		my $c = $seq =~ s/c/C/ig || '0';
		my $t = $seq =~ s/t/T/ig || '0';
		my $n = $seq =~ s/n/N/ig || '0';
		my $total = $a+$g+$c+$t+$n;
		unless ($total == $length) {print STDERR "WARNING: Sequence contains letters that are not A, G, C, T, or N\n"}
		print "A\tG\tC\tT\tN\tLine total\n";
		print "$a\t$g\t$c\t$t\t$n\t$total\n";
	}
}

sub N50 {
	my ($length, $total) = @_;
	my @contigsizes = sort {$b <=> $a} @$length;

	my $currsize=0;
	while ($currsize < int($total/2)) {
		my $l = shift @contigsizes;
		$currsize += $l;
	}
	return $contigsizes[0];
}

