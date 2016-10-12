#__perl__
#

use strict;
use lib '/usr/local/bioperl-dev/';
use Bio::TreeIO;

my $inf = shift || die "tree in newick format";
my $ouf = shift || die "output file to write phyloxml to";
my $input = new Bio::TreeIO(-file => $inf, -format => 'newick');
my $output = new Bio::TreeIO(-file => ">$ouf", -format => 'phyloxml');
while (my $tree = $input->next_tree()) {$output->write_tree($tree)}
