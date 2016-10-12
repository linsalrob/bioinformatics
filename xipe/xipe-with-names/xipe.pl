#!/usr/bin/perl

#Copyright (C) 2005 Beltran Rodriguez-Brito,
#Florent Angly, Pat MacNarnie, and Rober Edwards
#
#This program is free software; you can redistribute it and/or
#modify it under the terms of the GNU General Public License
#as published by the Free Software Foundation; either version 2
#of the License, or (at your option) any later version.
#
#This program is distributed in the hope that it will be useful,
#but WITHOUT ANY WARRANTY; without even the implied warranty of
#MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#GNU General Public License for more details.
#
#You should have received a copy of the GNU General Public License
#along with this program; if not, write to the Free Software
#Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301, USA.
#
#Released as part of the xipe package as supplemental online material
#to "A statistical approach for metagenome comparisons"
# by Beltran Rodriguez-Brito, Forest Rohwer, and Robert A. Edwards.

use strict;
use warnings;


my $execloc='.';


# Usage message
my $usage=<<EOF;
$0 usage:

Mandatory arguments
  -f Blast input file for sample 1
  -g Blast input file for sample 2
  -k keyfile
  -o output file name

  The Blast file should be in the tab-delimited format.
  The keyfile gives the correspondance between fig_ID and subsystem/role name.

Options
  -h name of sample 1 (no space!)
  -i name of sample 2 (no space!)
  -c confidence level (default: 0.98)
  -s sample size (usually between 1 and 10,000; default: 5,000)
  -r number of repeats (default: 20,000)
  -t temporary output filename (will be a unique id if not provided)
  -l level of analysis, i.e. 'subsystem' or 'role' (default: subsystem)
  -p filter out the Blast lines not containing the PEGs listed in this file (one per line)

NOTE: YOU SHOULD RUN THIS PROGRAM IN THE FOLDER CONTAINING ALL YOUR INPUT FILES!

EOF


# Argument processing
# Get args
my ($inf1, $inf2, $sample, $rep, $tempf, $name1, $name2, $conf, $outputf,$keyfile, $level, $pegfile);
while (@ARGV) {
 my $t=shift;
 if ($t eq "-f") {$inf1=shift}
 elsif ($t eq "-g") {$inf2=shift}
 elsif ($t eq "-k") {$keyfile=shift}
 elsif ($t eq "-o") {$outputf=shift}
 elsif ($t eq "-h") {$name1=shift}
 elsif ($t eq "-i") {$name2=shift}
 elsif ($t eq "-c") {$conf=shift}
 elsif ($t eq "-s") {$sample=shift}
 elsif ($t eq "-r") {$rep=shift}
 elsif ($t eq "-t") {$tempf=shift}
 elsif ($t eq "-l") {$level=shift}
 elsif ($t eq "-p") {$pegfile=shift}
 else {die "$usage"}
}
# Enforce mandatory arguments and set up default values for the others
die "$usage" unless ($inf1 && $inf2 && $keyfile && $outputf);
$name1=$inf1 unless (defined $name1);
$name2=$inf2 unless (defined $name2);
$conf=0.98 unless (defined $conf);
$sample=5000 unless ($sample);
$rep=20000 unless ($rep);
$level='subsystem' unless ($level);
$tempf=$$ unless ($tempf);
if (-e $tempf) {
 my $c=1;
 while (-e $tempf.$c) {$c++}
 $tempf=$tempf.$c;
}


# Xipe start message
print STDOUT "Started running Xipe at ", scalar(localtime(time)), "\n";
my $starttime=time;
my $time=$starttime;
open(OUT, ">$tempf") || die "Can't open $tempf";
print OUT "Running Xipe with the following arguments:\n  File 1: $inf1\n  File 2: $inf2\n  Keyfile: $keyfile\n  Sample: $sample\n  Repeats: $rep\n  Confidence: $conf\n  Temp files: ${tempf}_medians.txt and ${tempf}_ranges.txt\n";
close OUT;


# Pre-process the data: from human-readable to xipe-readable
my $command = "$execloc/xipe-pre.perl -i $inf1 -j $inf2 -k $keyfile -f $pegfile";
print STDOUT "- Running '$command'\n";
system("$command");
# input files and key file for the next steps
my $inf11;
my $inf22;
my $key;
if ($level eq 'subsystem') {
 $key = "$keyfile.$inf1.$inf2.sub";
 $inf11 = "$inf1.sub.pre";
 $inf22 = "$inf2.sub.pre";
} elsif ($level eq 'role') {
 $key = "$keyfile.$inf1.$inf2.role";
 $inf11 = "$inf1.role.pre";
 $inf22 = "$inf2.role.pre"; 
} else {
 die $usage;
}
print STDOUT "  xipe-pre took ", time-$time, " seconds\n"; $time=time;


# Xipe part 1
$command = "$execloc/xipe -r $rep -s $sample -f $inf11 -g $inf22 -o ${tempf}_medians.txt -p ${tempf}_ranges.txt";
print STDOUT "- Running '$command'\n";
system("$command");
system("rm -f $inf1.sub.pre");
system("rm -f $inf1.role.pre");
system("rm -f $inf2.sub.pre");
system("rm -f $inf2.role.pre");
print STDOUT "  xipe took ", time-$time, " seconds\n"; $time=time;


# Xipe part 2
$command = "$execloc/xipe2.perl ${tempf}_medians.txt > ${tempf}_medians_reassigned.txt";
print STDOUT "- Running '$command'\n";
system("$command");
system("rm -f ${tempf}_medians.txt");
print STDOUT "  xipe2-medians took ", time-$time, " seconds\n"; $time=time;


# Xipe part 2 (bis)
$command = "$execloc/xipe2.perl ${tempf}_ranges.txt > ${tempf}_ranges_reassigned.txt";
print STDOUT "- Running '$command'\n";
system("$command");
system("rm -f ${tempf}_ranges.txt");
print STDOUT "  xipe2-ranges took ", time-$time, " seconds\n"; $time=time;


# Xipe part 3
$command = "$execloc/xipe3.perl ${tempf}_medians_reassigned.txt ${tempf}_ranges_reassigned.txt $rep $conf $name1 $name2 > $outputf.post";
print STDOUT "- Running '$command'\n";
system("$command");
system("rm -f ${tempf}_medians_reassigned.txt");
system("rm -f ${tempf}_ranges_reassigned.txt");
print STDOUT "  xipe3 took ", time-$time, " seconds\n"; $time=time;


# Post-process the data:  from xipe-readable to human-readable
$command = "$execloc/xipe-post.perl -i $outputf.post -k $key > $outputf";
print STDOUT "- Running '$command'\n";
system("$command");
system("rm -f $outputf.post");
system("rm -f $keyfile.$inf1.$inf2.role");
system("rm -f $keyfile.$inf1.$inf2.sub");
system("rm -f $tempf");
print STDOUT "  xipe-post took ", time-$time, " seconds\n"; $time=time;


# Final message
print STDOUT "Xipe successfully finished in ", time-$starttime, " seconds\n\n";


exit;
