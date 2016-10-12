#__perl__
#
# Print a tuple of protein ID, LOCUS and phage name

my $gbk=shift || die "$0 <genbank file to parse?>";
open(IN, $gbk) || die $!;
my ($locus, $def)=("", "");
while (<IN>) {
	chomp;
	if (/^LOCUS\s+(\S+)/) {$locus = $1}
	if (/^DEFINITION\s+(\S+.*)$/) {$def = $1}
	if (/\/protein_id="(.*?)"/) {
		print join("\t", $1, $locus, $def), "\n";
	}
	if (m#^//$#) {($locus, $def)=("", "")}
}
close IN;
