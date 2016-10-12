#!/usr/bin/perl -w

# count all the files in directories in a directory

unless (@ARGV) {die "$0 <directories?>\n"}

foreach $dir (@ARGV) {
 unless (-d $dir) {print STDERR "$dir is not a dir. Ignoring\n"; next}
 opendir(DIR, $dir) or die "Can't open $dir\n";
 @files = readdir(DIR);
 print "$dir\t", $#files-1,"\n"; # note -1 because it should be +1 but then there is . and ..
 $total += $#files-1;
}

print "\nTOTAL:\t$total\n";
