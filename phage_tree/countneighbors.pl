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




use DBI;
use strict;
use OGD;
$| =1;


my $dbh=DBI->connect("DBI:mysql:phage", "SQLUSER", "SQLPASSWORD") or die "Can't connect to database\n";
my $file= shift || &niceexit("countneighbors <file of prot transfers>\n");

my %transfer;
open (IN, $file) || &niceexit("can't open $file\n");
while (<IN>) {
	my ($genome, $gene, $freq) = split /\t/;
	$transfer{$gene} = $genome;
	}

my ($name, $orf, $complement, $genomelength) = &getdata();


my %totaltransfer; # number of genes transfered per phage
my %one; # number of genes where one neighbor transfered
my %two; # number of genes where both neighbors transfered
my %none; #number of genes where both neighbors transfered
my (%same3of3, %same2of3, %same0of3);
my (%same2of2, %same0of2);

# we will also check the operon theory -- that operons are more likely to be transfered than non-operons

# for three genes transfered we can have the following options:
# ->->-> or <-<-<- all the same orintation - $same_threeofthree
# <-->-> or -><-<- or ->-><- or <-<--> two out of three the same $same_twoofthree
# <--><- or -><--> none the same $same_noneofthree

# for two transfered we can have
# ->-> or <-<- both the same: $same_twooftwo
# <--> or -><- none the same : $same_noneoftwo

foreach my $genome (keys %$orf) {
	my @allorfs = sort {$a <=> $b} @{$$orf{$genome}};
	foreach my $x (0 .. $#allorfs) {
		next unless ($transfer{$allorfs[$x]});
		$totaltransfer{$genome}++;
		
		my $prev; my $next; # previous and next ORFs to look at
		unless ($x) {$prev = $#allorfs} else {$prev = $x-1}
		if ($x == $#allorfs) {$next = 0} else {$next = $x+1}
		
		# check for transfers
		if ($transfer{$prev} && $transfer{$next}) {
			$two{$genome}++;
			# these two if loops do the operon test
			if ($$complement{$allorfs[$x]}) {
				if ($$complement{$allorfs[$prev]} && $$complement{$allorfs[$next]})
				  {$same3of3{$genome}++}
				elsif ($$complement{$allorfs[$prev]} || $$complement{$allorfs[$next]})
				  {$same2of3{$genome}++}
				else {$same0of3{$genome}++}
			}
			else {
				if ($$complement{$allorfs[$prev]} && $$complement{$allorfs[$next]})
				  {$same0of3{$genome}++}
				elsif ($$complement{$allorfs[$prev]} || $$complement{$allorfs[$next]})
				  {$same2of3{$genome}++}
				else {$same3of3{$genome}++}
			}
		}
		elsif ($transfer{$prev}) {
			$one{$genome}++;
			# this does the operon test
			if ($$complement{$allorfs[$x]}) {
				if ($$complement{$allorfs[$prev]}) {$same2of2{$genome}++}
				unless ($$complement{$allorfs[$prev]}) {$same0of2{$genome}++}
			}
			else {
				if ($$complement{$allorfs[$prev]}) {$same0of2{$genome}++}
				unless ($$complement{$allorfs[$prev]}) {$same2of2{$genome}++}
			}
		}
		elsif ($transfer{$next}) {
			$one{$genome}++;
			# this does the operon test
			if ($$complement{$allorfs[$x]}) {
				if ($$complement{$allorfs[$next]}) {$same2of2{$genome}++}
				unless ($$complement{$allorfs[$next]}) {$same0of2{$genome}++}
			}
			else {
				if ($$complement{$allorfs[$next]}) {$same0of2{$genome}++}
				unless ($$complement{$allorfs[$next]}) {$same2of2{$genome}++}
			}
			
		}
		else {$none{$genome}++}	
		}
	}

print "Name\tGenome Length\t# ORFS\tTotal transfer\tNone\tOne\tTwo\t0 of 2 in same dir\t2 of 2 in same dir\t0 of 3 in same dir\t2 of 3 in same dir\t3 of 3 in same dir\n";

foreach my $genome (sort {$a <=> $b} keys %$orf) {
	print "$$name{$genome}\t$$genomelength{$genome}\t";
	print $#{$$orf{$genome}}+1, "\t";
	if ($totaltransfer{$genome}) {print "$totaltransfer{$genome}\t"} else {print "0\t"}
	if ($none{$genome}) {print "$none{$genome}\t"} else {print "0\t"}
	if ($one{$genome}) {print "$one{$genome}\t"} else {print "0\t"}
	if ($two{$genome}) {print "$two{$genome}\t"} else {print "0\t"}
	if ($same0of2{$genome}) {print "$same0of2{$genome}\t"} else {print "0\t"}
	if ($same2of2{$genome}) {print "$same2of2{$genome}\t"} else {print "0\t"}
	if ($same0of3{$genome}) {print "$same0of3{$genome}\t"} else {print "0\t"}
	if ($same2of3{$genome}) {print "$same2of3{$genome}\t"} else {print "0\t"}
	if ($same3of3{$genome}) {print "$same3of3{$genome}\n"} else {print "0\n"}
	}




		
	










sub niceexit {
	my $reason = shift;
	$dbh->disconnect;
	if ($reason) {print STDERR $reason; exit(-1)}
	else {exit(0)}
	}

sub getdata {
print STDERR "Getting all data\n";
	my %name;
	my %orf;
	my %complement;
	my %genomelength;
	my $exc = $dbh->prepare("SELECT count, organism, sequence from phage" ) or croak $dbh->errstr;
	$exc->execute or die $dbh->errstr;
	while (my @retrieved = $exc->fetchrow_array) {$name{$retrieved[0]} = $retrieved[1]; $genomelength{$retrieved[0]}=length($retrieved[2])}
	
	
	foreach my $key (keys %name) {
		$exc = $dbh->prepare("SELECT count, organism, complement from protein where organism = $key" ) or croak $dbh->errstr;
		$exc->execute or die $dbh->errstr;
		while (my @retrieved = $exc->fetchrow_array) {
			push (@{$orf{$retrieved[1]}}, $retrieved[0]);
			$complement{$retrieved[0]}= $retrieved[2];
			}
		}
print STDERR "\tDone\n";

print "FOR WHOLE GENOME\nName\t#ORFs\tSame 0 of 3\tSame 2 of 3\tSame 3 of 3\n";
	
	foreach my $genome (sort keys %orf) {
	my ($same2of3, $same0of3, $same3of3) = (0, 0, 0);
		my @orfs = @{$orf{$genome}};
		foreach my $x (0 .. $#orfs) {
			my $prev; my $next; # previous and next ORFs to look at
			unless ($x) {$prev = $#orfs} else {$prev = $x-1}
			if ($x == $#orfs) {$next = 0} else {$next = $x+1}

			# these two if loops do the operon test for all the orfs in the genome
			if ($complement{$orfs[$x]}) {
				if ($complement{$orfs[$prev]} && $complement{$orfs[$next]})
				  {$same3of3++}
				elsif ($complement{$orfs[$prev]} || $complement{$orfs[$next]})
				  {$same2of3++}
				else {$same0of3++}
			}
			else {
				if ($complement{$orfs[$prev]} && $complement{$orfs[$next]})
				  {$same0of3++}
				elsif ($complement{$orfs[$prev]} || $complement{$orfs[$next]})
				  {$same2of3++}
				else {$same3of3++}
			}
		}

	print "$name{$genome}\t", $#{$orf{$genome}}+1, "\t$same0of3\t$same2of3\t$same3of3\n";
	}
	
	
	

	return \%name, \%orf, \%complement, \%genomelength;
}




