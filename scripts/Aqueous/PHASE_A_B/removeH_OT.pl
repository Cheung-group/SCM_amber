#!/usr/bin/perl -w
#
# remove H in the pdb file

if(@ARGV < 1) {
        print "Usage: removeH.pl file.pdb\n";
        exit;
        }
  
$pdb = shift @ARGV;

open(PDB,"$pdb") or die "failed to open $pdb, ";

while(<PDB>){
  next unless /^ATOM/;
  $atom = substr($_,13,1);
  next if $atom =~ "H";
  $atom = substr($_,12,1);
  next if $atom =~ "H";
  $atom = substr($_,13,2);
  next if $atom =~ "OT";
  $atom = substr($_,12,2);
  next if $atom =~ "OT";
  $atom = substr($_,13,2);
  next if $atom =~ "OX";
  $atom = substr($_,12,2);
  next if $atom =~ "OX";
  print $_;
}

