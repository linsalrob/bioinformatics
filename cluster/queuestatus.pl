#!/usr/bin/perl -w

use strict;
use SGE;

my $sge=SGE->new(-verbose=>1);
my $stat=$sge->status();
print map {"$_ ==> @{$stat->{$_}}\n"} keys %$stat;
print "\n\n";
my @jobs=$sge->all_jobs;
print "There are ", scalar @jobs,  " jobs running (", join " ", @jobs, ")\n";
foreach my $job (@jobs) {
 my $bjs=$sge->brief_job_stats($job);
 print "$job is running on ", $bjs->[0], "\n";
}
