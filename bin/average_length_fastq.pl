#!/usr/bin/perl -w

# count the longest, shortest, and average length for a set of fastq sequences.

my $f=shift || die "fastq file to work on?";
open(IN, $f) || die "Can't open $f";
my $max=[0, ""];
my $min=[1e6, ""];
my $total=0;
my $n=0;
my $line=0; 

# According to the fastq standard the ids, and sequences should only be on one line
# therefore we will just use a line count to see where we are:
# 0 -> ID
# 1 -> Sequence
# 2 -> Quality ID
# 3 -> Quality scores
my $posn=0;

while (<IN>) {
	chomp;
	$line++;
	if ($posn == 0) {$id = $_}
	elsif ($posn == 1) {&analyze_length(length($_), $id)}
	elsif ($posn == 2) {}
	elsif ($posn == 3) {$posn=-1}
	$posn++;
}

print "Shortest sequence: ", $min->[1], " (", $min->[0], " bp)\n";
print "Longest  sequence: ", $max->[1], " (", $max->[0], " bp)\n";
print "Average length: ", int(($total/$n)*100)/100, " from $n sequences and $total bp\n";

sub analyze_length {
	my ($len, $id) = @_;
	if ($len > $max->[0]) {$max=[$len, $id]}
	if ($len < $min->[0]) {$min=[$len, $id]}
	$total+=$len;
	$n++;
}
