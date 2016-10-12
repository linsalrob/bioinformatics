#!/usr/bin/perl -w

# another wrapper around Beltran's program xipe that allows everything to run on the SGE
# (c) Rob Edwards 6/27/05

use strict;

my $usage=<<EOF;
$0

Files. Only these three file names are required.
-f input file 1
-g input file 2
-o output file name

These should be in the format 
sequence_no	subsys_no	role_no

-c confidence level (default = 0.98)
-s sample size (usually between 1 and 10,000; default 5,000)
-r number of repeats (default=20,000)

-t temporary output filename (will be a unique id if not provided)


EOF

my $execloc='/home/rob/beltran';

my ($inf1, $inf2, $sample, $rep, $tempf, $conf, $outputf);
while (@ARGV) {
 my $t=shift;
 if ($t eq "-f") {$inf1=shift}
 elsif ($t eq "-g") {$inf2=shift}
 elsif ($t eq "-c") {$conf=shift}
 elsif ($t eq "-s") {$sample=shift}
 elsif ($t eq "-r") {$rep=shift}
 elsif ($t eq "-t") {$tempf=shift}
 elsif ($t eq "-o") {$outputf=shift}
}

die $usage unless ($inf1 && $inf2 && $outputf);
$conf=0.98 unless (defined $conf);
$rep=20000 unless ($rep);
$sample=5000 unless ($sample);
$tempf=$$ unless ($tempf);
if (-e $tempf) {
 my $c=1;
 while (-e $tempf.$c) {$c++}
 $tempf=$tempf.$c;
}

open(OUT, ">$tempf") || die "Can't open $tempf";
print OUT "Running xipe with the following commands:\nSample: $sample\nRepeats: $rep\nConfidence: $conf\nFile 1: $inf1\nFile 2: $inf2\nTemp files: ${tempf}_medians.txt and ${tempf}_ranges.txt\n";
close OUT;

# run the first part of xipe
print STDERR "Started running xipe at ", scalar(localtime(time)), "\n";
my $time=time;

print STDERR "Running $execloc/xipe -r $rep -s $sample -f $inf1 -g $inf2 -o ${tempf}_medians.txt -p ${tempf}_ranges.txt\n";
system("$execloc/xipe -r $rep -s $sample -f $inf1 -g $inf2 -o ${tempf}_medians.txt -p ${tempf}_ranges.txt");
print STDERR "xipe took ", time-$time, " seconds\n"; $time=time;

print STDERR "Running $execloc/xipe2.perl ${tempf}_medians.txt > ${tempf}_medians_reassigned.txt\n";
system("$execloc/xipe2.perl ${tempf}_medians.txt > ${tempf}_medians_reassigned.txt");
system("rm -f ${tempf}_medians.txt");
print STDERR "xipe2 - medians took ", time-$time, " seconds\n"; $time=time;

print STDERR "Running $execloc/xipe2.perl ${tempf}_ranges.txt > ${tempf}_ranges_reassigned.txt\n";
system("$execloc/xipe2.perl ${tempf}_ranges.txt > ${tempf}_ranges_reassigned.txt");
system("rm -f ${tempf}_ranges.txt");
print STDERR "xipe2 - ranges took ", time-$time, " seconds\n"; $time=time;

print STDERR "Running $execloc/xipe3.perl ${tempf}_medians_reassigned.txt ${tempf}_ranges_reassigned.txt $rep $conf $inf1 $inf2 > $outputf\n";
system("$execloc/xipe3.perl ${tempf}_medians_reassigned.txt ${tempf}_ranges_reassigned.txt $rep $conf $inf1 $inf2 > $outputf");
system("rm -f ${tempf}_medians_reassigned.txt ${tempf}_ranges_reassigned.txt");
print STDERR "xipe3 - the final xipe took ", time-$time, " seconds\n"; $time=time;



