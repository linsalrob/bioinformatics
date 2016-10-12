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



# count the number of times each protein occurs in protdist output as long as its score is less than the supplied cutoff.


use DBI;
use strict;

my $dbh=DBI->connect("DBI:mysql:phage", "apache") or die "Can't connect to database\n";
my $dir= shift || niceexit("countproteinmatches.byscore.pl <directory of protdists> <cutoff score>\n");
my $cutoff = shift || niceexit("countproteinmatches.byscore.pl <directory of protdists> <cutoff score>\n");



my %totalscore; my %count; my %matchcount;
opendir(DIR, $dir) || niceexit("Can't open $dir\n");
while (my $file = readdir(DIR)) {
	next if ($file =~ /^\./);
        open (IN, "$dir/$file") || niceexit("Can't open $dir/$file");
	my $match;
	my @proteins = '0';
	while (<IN>) {
		chomp;
		next unless (/(\d+)_(\d+)/);
		push (@proteins, $1);
		next unless (/^$file/);
		$match = $_;
		}
	my @match = split /\s+/, $match;
	
	my $proteintag = $match[0];
	$proteintag =~ s/_.*//;
	foreach $match (@match) {
		next if ($match eq $match[0]);
		if (($match < $cutoff) && ($match >0) ) {$matchcount{$proteintag}++}
		}
	
	my $skipped; my @skip;
	foreach my $x (1 .. $#match) {
		next if ($match[$x] < 0);
		if ($match[$x]<$cutoff) {
		unless ($proteins[$x]) {print STDERR "WARNING: ERROR with $x in $file\n"}
			$count{$proteins[$x]}++;
			$totalscore{$proteins[$x]}+=$match[$x];
			}
		else {$skipped++; push (@skip, $proteins[$x])}
		}
	if ($skipped) {print STDERR "$file, skipped $skipped: ", join (", ", @skip), "\n"}
}


print "Protein ID\tOccurences\tTotal score\tAverage Score\n";
my %order;

my @proteins = sort {$count{$b} <=> $count{$a}} keys %count;
foreach my $prot (@proteins) {
#	print "$prot\t$count{$prot}\t$totalscore{$prot}\t", $totalscore{$prot}/$count{$prot}, "\n";
$order{$count{$prot}}++;
}

my $allproteins=getallprots();
{
my @missed;
my $proteincount;
foreach my $prot (keys %$allproteins) {
	unless (exists $count{$prot}) {push (@missed, $prot)}
	$proteincount++;
	}

$order{0} = $#missed+1;
my @order = sort {$a <=> $b} keys %order;
print "Count\tNo. of occurences ($proteincount total proteins)\n";
foreach my $order (@order) {print "$order\t$order{$order}\n"}
}

{
print "\n\n\nMATCH COUNTS\n\n";
my @missed;
foreach my $prot (keys %$allproteins) {
	unless (exists $matchcount{$prot}) {push (@missed, $prot)}
	}

my %matchorder;
$matchorder{0}=$#missed+1;
foreach my $prot (keys %matchcount) {$matchorder{$matchcount{$prot}}++}
foreach my $match (sort {$a <=> $b} keys %matchorder) {
	print "$match\t$matchorder{$match}\n";
	}
}



sub getallprots {
   my %genomeproteins;
   my $exc = $dbh->prepare("select count,organism from protein");
   $exc->execute or die $dbh->errstr;
   while (my @ret = $exc->fetchrow_array) {$genomeproteins{$ret[0]}=$ret[1]}
   return \%genomeproteins;
   }

sub niceexit {
   my $reason = shift;
   $dbh->disconnect;
   if ($reason) {print STDERR $reason; exit(-1)}
   else {exit(0)}
   }
