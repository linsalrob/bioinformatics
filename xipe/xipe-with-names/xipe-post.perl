#! /usr/bin/perl


# xipe-pre.pl
# Florent Angly
# Converts Xipe-redable data to Human-readable data


use warnings;
use strict;


# Usage
my $usage =<<STRING
$0 converts Xipe-readable data to Human-readable data
Usage is: $0 -i result_file -k keyfile

STRING
;


# Get arguments
my ($filetoconvert, $keyfile);
if (scalar @ARGV == 0) {die "$usage"};
while (@ARGV) {
   my $arg = shift @ARGV;
   if ("$arg" eq "-i") {
      $filetoconvert= shift @ARGV;
   } elsif ("$arg" eq "-k") {
      $keyfile = shift @ARGV;
   } else {
      die "$usage";
   }
}


# Make sure the presence of the mandatory arguments are enforced
if (not($filetoconvert) || not($keyfile)) {
   die "$usage";
}


# Load subsystem/role hashes
my %translhash;
open TRANSL, $keyfile || die "Cannot open file '$keyfile': $!\n\n";
while (<TRANSL>) {
   if (m/^(.*)\t(.*)$/) {
      $translhash{$2} = $1;
   } else {
      die "Cannot parse file '$keyfile': Bad file format\n\n";
   }
}
close TRANSL;


# Make the conversion
open RESULT, "$filetoconvert";
my @result = <RESULT>; 
close RESULT;
my $regexp = '^ (\d+) (\S+)$';
foreach my $line(@result) {
   if ($line =~ m/$regexp/) {
      $line =~ s/$regexp/$translhash{$1}\t$2/;
   }
   print "$line";
}


exit;
