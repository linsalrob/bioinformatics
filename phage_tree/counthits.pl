#!/usr/bin/perl -w

#    Copyright 2001, 20002 Rob Edwards
#    For updates, more information, or to discuss the scripts
#    please contact Rob Edwards at redwards@utmem.edu or via http://www.salmonella.org/
#
#    This file is part of The Phage Proteome Scripts developed by Rob Edwards.
#
#    Tnese scripts are free software; you can redistribute and/or modify
#    them under the terms of the GNU General Public License as published by
#    the Free Software Foundation; either version 2 of the License, or
#    (at your option) any later version.
#
#    They are distributed in the hope that they will be useful,
#    but WITHOUT ANY WARRANTY; without even the implied warranty of
#    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#    GNU General Public License for more details.
#
#    You should have received a copy of the GNU General Public License
#    in the file (COPYING) along with these scripts; if not, write to the Free Software
#    Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA
#
#



# count all the protdist and blast hits from the tree generation
# we will store all the blast scores, and all the protdist scores so we can correlate them

use strict;

my $usage = "counthits.pl <protdist dir> <blast dir> <output dir>\nOPTIONS\n\t-p protdist cutoff\n";
$usage .= "\t-e E value cutoff (eg -2)\n\n\t-v verbose output";


my $protdistdir=shift || die $usage; $protdistdir =~ s/\/$//;
my $blastdir = shift || die $usage; $blastdir =~ s/\/$//;
my $outdir = shift || die $usage; $blastdir =~ s/\/$//;

unless (-e $outdir) {mkdir $outdir, 0755}


my $args = join " ", @ARGV;
my $pcut = 10000;
my $ecut=10000;
my $verb;
if ($args =~ /-v/) {$verb = 1}
if ($args =~ /-p\s+(\d+)/) {$pcut = $1}
if ($args =~ /-e\s+(\S+)/) {$ecut = 10**$1}

my $prot = read_protdist($protdistdir); # $prot is a pointer to a hash with each protein's average protdist score
my $blast = read_blasts($blastdir); # blasts is a pointer to a hash  with all the e values

# print out all the e-value/protdist comparisons
print STDERR "Printing e-value/protdist comparisons\n";
open (OUT, ">$outdir/evalue.protdist.txt") || die "can't open evalue.protdist.txt";
print OUT "Protein Distance\tE value\n";
foreach my $x (keys %$prot) {
  foreach my $y (keys %{$$prot{$x}}) {
    next unless (exists ${$$prot{$x}}{$y});
    next unless (exists ${$$blast{$x}}{$y});
    next if ($x eq $y);
    print OUT "${$$prot{$x}}{$y}\t${$$blast{$x}}{$y}\n";
    }
  }
  

# figure out the protein distance distribution

count_protdist($prot);

count_blasts($blast);












sub read_protdist {
  my $dir = shift;
  print STDERR "Reading protdist\n";
  opendir(DIR, $dir) || die "Can't read $dir\n";
  my %match; my %count; my %similar; my $max=1; my $maxhit;
  # values used to count the protdist output
  my ($point01, $point1, $one, $two, $three, $four, $five, $six, $seven, $eight, $nine, $ten, $twenty);
  my ($thirty, $fourty, $fifty, $hundred, $unmatched);
  while (my $file = readdir(DIR)) {
    next if ($file =~ /^\./);
    open (IN, "$dir/$file") || print STDERR "Can't open $dir/$file\n";
    my %hits; my @seqs;
    while (<IN>) {
      # read all the protdists and proteins. Store them in a hash
      next unless (/\_/); # skip the first line
      my @line = split;
      my $seq = shift (@line);
      if ($seq eq $file) { # this is the driving protein. How many is it similar to
        my $similar = $#line+1;
	$similar{$similar}++;
	if ($similar > $max) {$max=$similar; $maxhit=$seq}
	}
      push (@seqs, $seq);
      @{$hits{$seq}} = @line;
      }
    # read all the hashes, and sum all the protein distances
    # also keep a count of it
    foreach my $seq (@seqs) {
      foreach my $x (0 .. $#seqs) {
        my $got;
        if (${$hits{$seq}}[$x] == -1) {${$hits{$seq}}[$x] = 100; $unmatched++; $got=1}
        $match{$seq}{$seqs[$x]} += ${$hits{$seq}}[$x];
	$count{$seq}{$seqs[$x]} ++;
	next if ($seq eq $seqs[$x]);
	next if ($got); # ignore unmatched hits
	# record the cumulative scores
	if (${$hits{$seq}}[$x] <= 0.01) {$point01++}
	elsif (${$hits{$seq}}[$x] <= 0.1) {$point1++}
	elsif (${$hits{$seq}}[$x] <= 1) {$one++}
	elsif (${$hits{$seq}}[$x] <= 2) {$two++}
	elsif (${$hits{$seq}}[$x] <= 3) {$three++}
	elsif (${$hits{$seq}}[$x] <= 4) {$four++}
	elsif (${$hits{$seq}}[$x] <= 5) {$five++}
	elsif (${$hits{$seq}}[$x] <= 6) {$six++}
	elsif (${$hits{$seq}}[$x] <= 7) {$seven++}
	elsif (${$hits{$seq}}[$x] <= 8) {$eight++}
	elsif (${$hits{$seq}}[$x] <= 9) {$nine++}
	elsif (${$hits{$seq}}[$x] <= 10) {$ten++}
	elsif (${$hits{$seq}}[$x] <= 20) {$twenty++}
	elsif (${$hits{$seq}}[$x] <= 30) {$thirty++}
	elsif (${$hits{$seq}}[$x] <= 40) {$fourty++}
	elsif (${$hits{$seq}}[$x] <= 50) {$fifty++}
	else {$hundred++}
	
	}
    }
  }
  # average the counts
  foreach my $x (keys %match) {
    foreach my $y (keys %{$match{$x}}) {$match{$x}{$y} = $match{$x}{$y}/$count{$x}{$y}} 
  }
  # print out all the protein distances
  open (OUT, ">$outdir/protdist.distrib.all.txt") || die "Can't open $outdir/protdist.distrib.txt";
  print OUT "Score\tOccurences\n";
  print OUT "0.01\t$point01\n0.1\t$point1\n1\t$one\n2\t$two\n3\t$three\n4\t$four\n5\t$five\n6\t$six\n7\t$seven\n";
  print OUT "8\t$eight\n9\t$nine\n10\t$ten\n20\t$twenty\n30\t$thirty\n40\t$fourty\n50\t$fifty\n100\t$hundred\n";
  print OUT "infinite\t$unmatched\n";
  close OUT;
  
  #print out the number of matches
  open (OUT, ">$outdir/prot.similars.txt")  || die "Can't open $outdir/prot.similars.txt";
  print OUT "Number of similarities\tNumber of proteins\n";
  foreach my $k (sort {$a <=> $b} keys %similar) {print OUT "$k\t$similar{$k}\n"}
  print OUT "\n\nMaximum hit: $maxhit ($max)\n";
  close OUT;
  return \%match;
}
      
      
sub read_blasts {
  my $dir = shift;
  # values to store the data
  my ($zero, $one, $two, $three, $four, $five, $six)=(0,0,0,0,0,0,0);
  my ($seven, $eight, $nine, $ten, $twenty)=(0,0,0,0,0);
  my ($thirty, $fourty, $fifty, $hundred)=(0,0,0,0);
  
  print STDERR "Reading BLASTs\n";
  opendir(DIR, $dir) || die "Can't read $dir\n";
  my %match; my %count;
  while (my $file = readdir(DIR)) {
    next if ($file =~ /^\./);
    open (IN, "$dir/$file") || print STDERR "Can't open $dir/$file\n";
    my $query;
    while (<IN>) {
      if (/Query=\s+(\d+_\d+)/) {$query = $1; next}
      last if (/^>/);
      next unless (/\d+_\d+/);
      my @line = split;
      my $seq = $line[0];
      my $e = $line[$#line-1];
      if ($e =~ /e/) {
        my @e = split /e/, $e;
	$e = $e[0]*(10**$e[$#e]);
	}
if ($verb) {print STDERR "E value: $e\n"}

  # count the raw blasts
      if ($e >= 1) {$zero++}
      if ($e >= 10**-1) {$one++}
      elsif ($e >= 10**-2) {$two++}
      elsif ($e >= 10**-3) {$three++}
      elsif ($e >= 10**-4) {$four++}
      elsif ($e >= 10**-5) {$five++}
      elsif ($e >= 10**-6) {$six++}
      elsif ($e >= 10**-7) {$seven++}
      elsif ($e >= 10**-8) {$eight++}
      elsif ($e >= 10**-9) {$nine++}
      elsif ($e >= 10**-10) {$ten++}
      elsif ($e >= 10**-20) {$twenty++}
      elsif ($e >= 10**-30) {$thirty++}
      elsif ($e >= 10**-40) {$fourty++}
      elsif ($e >= 10**-50) {$fifty++}
      else {$hundred++}



      $match{$query}{$seq}+=$e;
      $count{$query}{$seq}++;
      }
    }
  # average the counts
  foreach my $x (keys %match) {
    foreach my $y (keys %{$match{$x}}) {$match{$x}{$y} = $match{$x}{$y}/$count{$x}{$y}} 
  }
  # print out the data
  open (OUT, ">$outdir/blasts.raw.txt") || die "Can't open $outdir/protdist.distrib.txt";
  print OUT "Cutoff\tOccurences\n";
  print OUT "1\t$zero\n0.1\t$one\n1e-2\t$two\n1e-3\t$three\n1-e4\t$four\n1e-5\t$five\n1e-6\t$six\n1e-7\t$seven\n";
  print OUT "1e-8\t$eight\n1e-9\t$nine\n1e-10\t$ten\n1e-20\t$twenty\n1e-30\t$thirty\n1e-40\t$fourty\n1e-50\t$fifty\n0\t$hundred\n";
  close OUT;
  return \%match;
}
      





sub count_protdist {
  my $protdist=shift;
  # we want these numbers :
  # less than 
  # 0.01 0.1 1 2 3 4 5 6 7 8 9 10 20 30 40 50 100
  my ($point01, $point1, $one, $two, $three, $four, $five, $six, $seven, $eight, $nine, $ten, $twenty);
  my ($thirty, $fourty, $fifty, $hundred);
  foreach my $x (keys %$protdist) {
    foreach my $y (keys %{$$protdist{$x}}) {
      next if ($x eq $y);
      next unless (exists ${$$protdist{$x}}{$y});
      if (${$$protdist{$x}}{$y} <= 0.01) {$point01++}
      elsif (${$$protdist{$x}}{$y} <= 0.1) {$point1++}
      elsif (${$$protdist{$x}}{$y} <= 1) {$one++}
      elsif (${$$protdist{$x}}{$y} <= 2) {$two++}
      elsif (${$$protdist{$x}}{$y} <= 3) {$three++}
      elsif (${$$protdist{$x}}{$y} <= 4) {$four++}
      elsif (${$$protdist{$x}}{$y} <= 5) {$five++}
      elsif (${$$protdist{$x}}{$y} <= 6) {$six++}
      elsif (${$$protdist{$x}}{$y} <= 7) {$seven++}
      elsif (${$$protdist{$x}}{$y} <= 8) {$eight++}
      elsif (${$$protdist{$x}}{$y} <= 9) {$nine++}
      elsif (${$$protdist{$x}}{$y} <= 10) {$ten++}
      elsif (${$$protdist{$x}}{$y} <= 20) {$twenty++}
      elsif (${$$protdist{$x}}{$y} <= 30) {$thirty++}
      elsif (${$$protdist{$x}}{$y} <= 40) {$fourty++}
      elsif (${$$protdist{$x}}{$y} <= 50) {$fifty++}
      else {$hundred++}
      }
    }
  open (OUT, ">$outdir/protdist.distrib.av.txt") || die "Can't open $outdir/protdist.distrib.txt";
  print OUT "Score\tOccurences\n";
  print OUT "0.01\t$point01\n0.1\t$point1\n1\t$one\n2\t$two\n3\t$three\n4\t$four\n5\t$five\n6\t$six\n7\t$seven\n";
  print OUT "8\t$eight\n9\t$nine\n10\t$ten\n20\t$twenty\n30\t$thirty\n40\t$fourty\n50\t$fifty\n100\t$hundred\n";
  close OUT;
}



sub count_blasts {
  my $blasts=shift;
  # we want these numbers :
  # less than 
  # 0.01 0.1 1 2 3 4 5 6 7 8 9 10 20 30 40 50 100
  my ($zero, $one, $two, $three, $four, $five, $six)=(0,0,0,0,0,0,0);
  my ($seven, $eight, $nine, $ten, $twenty)=(0,0,0,0,0);
  my ($thirty, $fourty, $fifty, $hundred)=(0,0,0,0);
  foreach my $x (keys %$blasts) {
   foreach my $y (keys %{$$blasts{$x}}) {
     next if ($x eq $y);
     next unless (exists ${$$blasts{$x}}{$y});
      if (${$$blasts{$x}}{$y} >= 1) {$zero++}
      if (${$$blasts{$x}}{$y} >= 10**-1) {$one++}
      elsif (${$$blasts{$x}}{$y} >= 10**-2) {$two++}
      elsif (${$$blasts{$x}}{$y} >= 10**-3) {$three++}
      elsif (${$$blasts{$x}}{$y} >= 10**-4) {$four++}
      elsif (${$$blasts{$x}}{$y} >= 10**-5) {$five++}
      elsif (${$$blasts{$x}}{$y} >= 10**-6) {$six++}
      elsif (${$$blasts{$x}}{$y} >= 10**-7) {$seven++}
      elsif (${$$blasts{$x}}{$y} >= 10**-8) {$eight++}
      elsif (${$$blasts{$x}}{$y} >= 10**-9) {$nine++}
      elsif (${$$blasts{$x}}{$y} >= 10**-10) {$ten++}
      elsif (${$$blasts{$x}}{$y} >= 10**-20) {$twenty++}
      elsif (${$$blasts{$x}}{$y} >= 10**-30) {$thirty++}
      elsif (${$$blasts{$x}}{$y} >= 10**-40) {$fourty++}
      elsif (${$$blasts{$x}}{$y} >= 10**-50) {$fifty++}
      else {$hundred++}
      }
    }
  open (OUT, ">$outdir/blasts.distrib.av.txt") || die "Can't open $outdir/protdist.distrib.txt";
  print OUT "Cutoff\tOccurences\n";
  print OUT "1\t$zero\n0.1\t$one\n1e-2\t$two\n1e-3\t$three\n1-e4\t$four\n1e-5\t$five\n1e-6\t$six\n1e-7\t$seven\n";
  print OUT "1e-8\t$eight\n1e-9\t$nine\n1e-10\t$ten\n1e-20\t$twenty\n1e-30\t$thirty\n1e-40\t$fourty\n1e-50\t$fifty\n0\t$hundred\n";
  close OUT;
}





