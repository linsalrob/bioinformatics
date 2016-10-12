#!/usr/bin/perl -w

=pod

=head1 split_files_by_column.pl

split a file into seperate files each with a different numbmer of columns. This allows you to take a big file and open it in excel where you are only allowed 255 columns.

=cut

# split into different numbers of columns

use strict;
my $usage=<<EOF;
$0

-c number of columns to split the files into
-f file to split
-k number of columns to keep and put at the beginning of each file. To keep just the first column (default) use -k 1. To keep the first two columns use -k 2

EOF

my ($file, $cols, $keep)=(undef, undef, 1);
while (@ARGV)
{
    my $t=shift;
    if ($t eq "-c") {$cols=shift}
    if ($t eq "-f") {$file=shift}
    if ($t eq "-k") {$keep=shift}
}


die $usage unless ($cols && $file);

my $filehandles;
open(IN, $file) || die "Can't open $file";
while (<IN>)
{
    chomp;
    my @a=split /\t/;
    my $posn=$keep;
    my $thisfile=0;
    my $end=$#a-$cols;
    while ($posn < $end)
    {
        $thisfile++;
        my $fh;
        if ($filehandles->{$thisfile}) {$fh=$filehandles->{$thisfile}}
        else
        {
            my $out=$file;
            ($out =~ /.txt$/) ? ($out =~ s/.txt$/.$thisfile.txt/) : ($out .= ".$thisfile");
            open($fh, ">$out") || die "Can't write to $out";
            $filehandles->{$thisfile}=$fh;
        }
        print $fh join("\t", @a[0..($keep-1)], @a[$posn .. $posn+$cols-1]), "\n";
        $posn+=$cols;
    }
    unless ($posn == $#a)
    {
        $thisfile++;my $fh;
        if ($filehandles->{$thisfile}) {$fh=$filehandles->{$thisfile}}
        else
        {
            my $out = $file;
            ($out =~ /.txt$/) ? ($out =~ s/.txt$/.$thisfile.txt/) : ($out .= ".$thisfile");
            open($fh, ">$out") || die "Can't write to $out";
            $filehandles->{$thisfile}=$fh;
        }
        print $fh join("\t", @a[0..($keep-1)], @a[$posn .. $#a]), "\n";
    }
}
