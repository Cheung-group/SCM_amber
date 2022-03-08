#!/usr/bin/perl 
#	hbond-dssp.pl out.dssp (dssp output file)
#         2005
#	  Margaret S. Cheung
#
# extract HB pairs from DSSP output 

#$breakaa=32;

unless(@ARGV == 1){
    die "usage: hbond-dssp.pl out1.dssp\n";
}

$dssp = shift @ARGV; 
open(INFILE,$dssp) or die "failed to open $dssp,";

$p=0;
while(<INFILE>) {

 	if (/#  RESIDUE/) {$p=1; }

	$head =substr($_,5,5); 
	# 32 is the fake residue created by DSSP
#	if($head > $breakaa) {$head = $head - 1;}

#	print $head;	 
    if ($p==1 && $head =~ /\d+/ )
      { 
	$sep1 = substr($_,41,4);
	$ene1 = substr($_,46,4);
	$pair1 = $head+$sep1;
	
        if($sep1>3) {
	 	#if($head>$breakaa){$head -= 1; }
		#if($pair1>$breakaa){$pair1 -= 1;}
	        print "$head  $sep1 $pair1 $ene1 \n"; 
	}

	$sep2 = substr($_,52,4);
	$ene2 = substr($_,57,4);
	$pair2 = $head+$sep2;
        if($sep2>3) {
	 	#if($head>$breakaa){$head -= 1;}
		#if($pair2>$breakaa){$pair2 -= 1;}
	        print "$head  $sep2 $pair2 $ene2 \n"; 
	}

	$sep3 = substr($_,63,4);
	$ene3 = substr($_,68,4);
	$pair3 = $head+$sep3;
        if($sep3>3) {
	 	#if($head>$breakaa){$head -= 1;}
		#if($pair3>$breakaa){$pair3 -= 1;}
	        print "$head  $sep3 $pair3 $ene3 \n"; 
	}

	$sep4 = substr($_,74,4);
	$ene4 = substr($_,79,4);
	$pair4 = $head+$sep4;
        if($sep4>3) {
	 	#if($head>$breakaa){$head -= 1;}
		#if($pair4>$breakaa){ $pair4 -= 1;}
		print "$head  $sep4 $pair4 $ene4 \n";
	}
     }
} #end of while
close INFILE;

