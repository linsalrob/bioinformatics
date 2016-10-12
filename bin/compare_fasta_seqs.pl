#!/usr/bin/perl -w

use strict;
use lib '/clusterfs/home/rob/perl';
use Rob;

my $f1 = shift || die "file1?";
my $f2 = shift || die "file2?";

print STDERR "Reading $f1\n";
my $fa1=Rob->read_fasta($f1);
my $seq1; my $count;
map {$fa1->{$_} =~ s/\s+//g; $count->{uc($fa1->{$_})}++; $seq1->{uc($fa1->{$_})}=$_} keys %$fa1;
print "Duplicates in $f1\n";
map {print  $seq1->{uc($fa1->{$_})}, "\n" if ($count->{$_} > 1)} keys %$count;
undef $count;

print STDERR "Reading $f2\n";
my $fa2=Rob->read_fasta($f2);
my $seq2;
map {$fa2->{$_} =~ s/\s+//g; $seq2->{uc($fa2->{$_})}=$_; $count->{uc($fa2->{$_})}++} keys %$fa2;
print STDERR "THis should be: ", $count->{'ZXD'}, "\n";
print "Duplicates in $f2\n";
map {print  $seq2->{uc($fa2->{$_})}, "\n" if ($count->{$_} > 1)} keys %$count;
undef $count;

print "Sequences in $f1 that are NOT in $f2\n";
map {print $seq1->{$_}, "\n" unless ($seq2->{$_})} keys %$seq1;

print "Sequences in $f2 that are NOT in $f1\n";
map {print $seq2->{$_}, "\n" unless ($seq1->{$_})} keys %$seq2;

print "Keys in $f1 that are not in $f2\n";
map {print $_, "\n" unless ($fa2->{$_})} keys %$fa1;

print "Keys in $f2 that are not in $f1\n";
map {print $_, "\n" unless ($fa1->{$_})} keys %$fa2;



