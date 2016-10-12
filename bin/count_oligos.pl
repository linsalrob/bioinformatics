#!/usr/bin/perl -w

=pod count_oligos.pl

=head1 count_oligos.pl

Count the occurence of different oligos in a sequence. This is the quickest way to count oligos upto about 11 or 12 nt. After that it is quicker to get the appropriate subsequences and hash them, because you can't make a hash of all possible 14-mers without running out of memory!

=cut


use strict;
use Rob;

my ($f, $min, $max)=@ARGV;

die "$0 <fasta file> <min default=3> <max default=10>" unless ($f);
defined $min or ($min=3);
defined $max or ($max=11);


my $oligos=&all_oligos($min, $max);

my $fa=Rob->read_fasta($f);

my %count;
foreach my $k (keys %$fa)
{
    my $seq=uc($fa->{$k});
    foreach my $o (keys %$oligos)
    {
        my $posn=index($seq, $o);
        while ($posn>-1)
        {
            $count{$o}++;
            $posn++;
            $posn=index($seq, $o, $posn);
        }
    }
}

print map {"$_\t$count{$_}\n"} sort {$count{$b} <=> $count{$a}} keys %count;





















sub all_oligos {
    my ($min, $max)=@_;
    my $time=time;
    my @oligos;
    my @nt=qw[G A T C];
    foreach my $o (@nt)
    {
        push @oligos, $o;
        foreach my $t (@nt)
        {
            push @oligos, $o.$t;
            foreach my $th (@nt)
            {
                push @oligos, $o.$t.$th;
                foreach my $f (@nt)
                {
                    push @oligos, $o.$t.$th.$f;
                    foreach my $fi (@nt)
                    {
                        push @oligos, $o.$t.$th.$f.$fi;
                        foreach my $s (@nt)
                        {
                            push @oligos, $o.$t.$th.$f.$fi.$s;
                            foreach my $se (@nt)
                            {
                                push @oligos, $o.$t.$th.$f.$fi.$s.$se;
                                foreach my $e (@nt)
                                {
                                    push @oligos, $o.$t.$th.$f.$fi.$s.$se.$e;
                                    foreach my $n (@nt)
                                    {
                                        push @oligos, $o.$t.$th.$f.$fi.$s.$se.$e.$n;
                                        foreach my $t (@nt)
                                        {
                                            push @oligos, $o.$t.$th.$f.$fi.$s.$se.$e.$n.$t;
                                            foreach my $el (@nt)
                                            {
                                                push @oligos, $o.$t.$th.$f.$fi.$s.$se.$e.$n.$t.$el;
                                                if ($max == 12)
                                                {
                                                    foreach my $tv (@nt)
                                                    {
                                                        push @oligos, $o.$t.$th.$f.$fi.$s.$se.$e.$n.$t.$el.$tv;
                                                    }
                                                }
                                            }
                                        }
                                    }
                                }
                            }
                        }
                    }
                }
            }
        }
    }
    my %want;
    my $totallen=0;
    map {if (length($_) >= $min && length($_) <= $max) {$want{$_}=1; $totallen+=length($_)}} @oligos;
    print STDERR "Returning ", scalar(keys %want), " oligos between $min and $max\nTotal length is $totallen\n";

    return \%want;
}

