#!/usr/bin/env /home/fig/FIGdisk/env/linux-debian-x86_64/bin/perl

BEGIN {
    @INC = qw(
              /home/fig/FIGdisk/dist/releases/current/linux-debian-x86_64/lib
              /home/fig/FIGdisk/dist/releases/current/linux-debian-x86_64/lib/FigKernelPackages
              /home/fig/FIGdisk/dist/releases/current/linux-debian-x86_64/lib/WebApplication
              /home/fig/FIGdisk/dist/releases/current/linux-debian-x86_64/lib/FortyEight
              /home/fig/FIGdisk/dist/releases/current/linux-debian-x86_64/lib/PPO
              /home/fig/FIGdisk/dist/releases/current/linux-debian-x86_64/lib/RAST
              /home/fig/FIGdisk/dist/releases/current/linux-debian-x86_64/lib/MGRAST
              /home/fig/FIGdisk/dist/releases/current/linux-debian-x86_64/lib/SeedViewer
              /home/fig/FIGdisk/dist/current/linux-debian-x86_64/lib
              /home/fig/FIGdisk/dist/current/linux-debian-x86_64/lib/FigKernelPackages
              /home/fig/FIGdisk/env/linux-debian-x86_64/lib/perl5/5.10.0/x86_64-linux
              /home/fig/FIGdisk/env/linux-debian-x86_64/lib/perl5/5.10.0
              /home/fig/FIGdisk/env/linux-debian-x86_64/lib/perl5/site_perl/5.10.0/x86_64-linux
              /home/fig/FIGdisk/env/linux-debian-x86_64/lib/perl5/site_perl/5.10.0
              .
              /home/fig/FIGdisk/config
 
);
}
use Data::Dumper;
use Carp;
use FIG_Config;
$ENV{'BLASTMAT'} = "/home/fig/FIGdisk/BLASTMAT";
$ENV{'FIG_HOME'} = "/home/fig/FIGdisk";
# end of tool_hdr
########################################################################

# get all the proteins that belong to phages as a single fasta file

use strict;
use FIG;
my $fig=new FIG;



foreach my $g (map {$_->[0]} $fig->get_attributes(undef, 'virus_type', 'Phage')) {
	print STDERR $fig->genus_species($g), " ($g)\n";
	foreach my $peg ($fig->pegs_of($g)) {
		print ">$peg\n", $fig->get_translation($peg), "\n";
	}
}


