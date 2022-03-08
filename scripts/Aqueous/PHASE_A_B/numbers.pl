#!/usr/bin/perl -w
################################################################################

if(@ARGV < 1) {
        print "Usage: numbers.pl file.pdb\n";
        exit;
        }

$pdb = shift @ARGV;

open(PDB,"$pdb") or die "failed to open $pdb, ";


@array = ();

my $natom = 0;
my $ntype = 0;

while (<PDB>) {

	   next unless /^ATOM/;

	   @array = split (" ",$_);

	
	   if($array[2] eq 'CA'){
	   		$natom = $natom + 1;
	   		$ntype = $ntype + 1;
	   }

	   if($array[2] eq 'CB'){
	   		$ntype = $ntype + 1;
	   }
}
print "$natom\n";
print "$ntype\n";

