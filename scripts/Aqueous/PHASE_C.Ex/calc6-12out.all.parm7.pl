#!/usr/bin/perl -w
# Script for calculating the 6-12 parameters in amber6 format.
# Assigning BT L-J interactions to ALL the contacts.
# Modified by Antonios Samiotakis in 11/21/09 in line 209
# so that   $BondEn = (1+$BondEn) * $GoEn; is changed to:
#           $BondEn = $GoEn + $BondEn; 
# Modified by Dirar (7/13/2011) to print the LJ Coefficient
# to two files to be compatibale with AMBER PARM7
#	Modified by Qian and Dirar (7.20.2012) to combine all calc_6-12 scripts  
#
### 	calc6-12out.BT.antonis.pl atomindex BTmap.dat pairindex crdfile

use strict;

my $ATOM = 245;
my $GoEn = 0.6;
my $RESNUM = 128;
my $MAXLENGTH = $ATOM;
my $FACTOR = 0.9; 
my $MAXCONTACT = 275;
my $NTYPE = $ATOM;
my $PAIRS = ($NTYPE*($NTYPE+1)/2);
my $MAXPAIRS = ($PAIRS+1);
my $SCALE = 0;
my $RG = 11.485;

my $sizefile = '';
my $hamfile = '';
my $mapfile = '';
my $crdfile = '';
my @coorx = ();
my @coory = ();
my @coorz = ();
my @contactmaph = ();
my @contactmapb = ();
my %hammap = ();
my @ressize = ();
my %restype = ();

my @godisth = ();
my @godistb = ();
my @repcn1 = (); #size dependent rep
my @cn1 = (); #size dependent rep
my @attcn2 = (); #size dependent epsilon
my @cn2 = (); #size dependent rep

my @convertindex = ();


if (@ARGV < 7){die "usage: 6-12out.pl ressize.dat ham.dat contactmap.dat crd htype utype crowdtype"; }

$sizefile = shift @ARGV;
$hamfile = shift @ARGV;
$mapfile = shift @ARGV;
$crdfile = shift @ARGV;
my $htype = shift @ARGV;
my $utype = shift @ARGV;
my $crowdtype = shift @ARGV;

open (INSIZE,$sizefile) or die "cant open $sizefile \n";
open (INHAM,$hamfile) or die "cant open $hamfile \n";
open (INMAP,$mapfile) or die "cant open $mapfile \n";
open (INCRD,$crdfile) or die "cant open $crdfile \n";

open (OUTr,">6-12out.dat.rep") or die "cant open 6-12out.dat.rep";
open (OUTa,">6-12out.dat.att") or die "cant open 6-12out.dat.att";
open (OUTNHB,">nnhb.dat") or die "cant open nnhb.dat";
open (OUTNSC,">nnsc.dat") or die "cant open nnhb.dat";

$NTYPE = $NTYPE+$crowdtype;
$PAIRS = ($NTYPE*($NTYPE+1)/2);
$MAXPAIRS = ($PAIRS+1);

#initialize
  for(my $j=1;$j<$MAXPAIRS;$j++)
{
	$repcn1[$j]=0;
	$attcn2[$j]=0;
}

& read_size ();
& read_crd ();
& read_ham ();
& read_map ();
& calc_dist();

if ($utype eq 'nourea')
{
if ($htype eq 'HOMO') # homo
	{
	$FACTOR = 1.0;
	$SCALE = 0.001;
	& calc_nn_homo_Go ();
	& calc_go_homo ();
	}
else {if ($htype eq 'GO') # Go
	{
	$FACTOR = 0.9;
	$SCALE = 1.0;
        & calc_nn_homo_Go ();
        & calc_go_Go ();
	}
else {if ($htype eq 'BT') # BT
	{
	$FACTOR = 0.9;
	$SCALE = 1.0;
        & calc_nn_BT ();
	}
        }}
}
else
{
if ($htype eq 'HOMO') # homo
        {
        $FACTOR = 1.0;
        $SCALE = 0.001;
        & calc_nn_homo_Go ();
        & calc_go_homo ();
        }
else {if ($htype eq 'GO') # Go
        {
        $FACTOR = 0.9;
        $SCALE = 1.0;
        & calc_nn_homo_Go ();
        & calc_go_Go_urea ();
        }
else {if ($htype eq 'BT') # BT
        {
        $FACTOR = 0.9;
        $SCALE = 1.0;
        & calc_nn_BT_urea ();
        }
        }}
}

& calc_xy (); #crowding agent and protein
& print_out ();


close OUTr;
close OUTa;
close INCRD;
close INMAP;
close INHAM;
close INSIZE;

sub print_out
{

printf(OUTr "FORMAT\n");
my $count = 0;
#fortran index

print "\n\n";

  for(my $j=1;$j<$MAXPAIRS;$j++)
    {
     printf(OUTr "  %.8le",$repcn1[$j]);
#	 print $j, " " ,$repcn1[$j],"\n";
      $count++;
     if( ( ($count%5) eq 0 ) || ($MAXPAIRS-$j)<2) {printf(OUTr "\n");}
    }




printf(OUTa "FORMAT\n");
$count = 0;
  for(my $j=1;$j<$MAXPAIRS;$j++)
    {

     if($attcn2[$j] < 0) 
	{ printf(OUTa " %.8le",$attcn2[$j]);
	$count++}
     else{
          printf(OUTa "  %.8le",$attcn2[$j]);
      $count++;}

     if( ( ($count%5) eq 0 ) || ($MAXPAIRS-$j)<2) {printf(OUTa "\n");}
    }

}


##### subroutines

sub calc_dist{
my $bond = '';
my $a1 = 0;
my $a2 = 0;
my $i=0;
my $j=0;

       for (  $i = 0 ; $i < ($RESNUM) ; $i ++)
       { for ( $j = 0 ; $j < ($RESNUM) ; $j ++){
        if (defined $contactmaph[$i][$j] eq  1 ){

	# hb;

	$a1= $convertindex[$i][0];
	$a2= $convertindex[$j][0];
	$bond=($coorx[$a1]-$coorx[$a2])**2+ ($coory[$a1]-$coory[$a2])**2+($coorz[$a1]-$coorz[$a2])**2;
	$bond = sqrt($bond);
	$godisth[$i][$j]=$bond;
	if ($i<$j ) {print "h ",$i," ",$a1," ",$j," ",$a2," ",$bond," ",$godisth[$i][$j],"\n";}
	   }
	if( defined $contactmapb[$i][$j] eq  1   )
	{
	#sidechain
        $a1= $convertindex[$i][1];
        $a2= $convertindex[$j][1];
        $bond=($coorx[$a1]-$coorx[$a2])**2+ ($coory[$a1]-$coory[$a2])**2+($coorz[$a1]-$coorz[$a2])**2;
        $bond = sqrt($bond);
        $godistb[$i][$j]=$bond;
#	if ($i<$j ) {print "b ",$i," ",$a1," ",$j," ",$a2," ",$bond," ",$godistb[$i][$j]," ", $coorx[$a1]," ",$coorx[$a2],"\n";}


	}

	  }
	}

}

sub calc_go_homo {

my $posen = 0;
my $negen = 0;
my $attr = '';
my $count=0; 
my $BondEn;
my $m = '';
my $n = '';
my $bond = '';
my $bond6 = '';
my $bond12 = '';
my $i = '';
my $j = '';
# evaluate the coefficients for 6-12 potential on specific pairwise go contacts	# get mapping fortran index
  
print "\n"; 

       for (  $i = 0 ; $i < ($RESNUM) ; $i ++)
       { for ( $j = 0 ; $j < ($RESNUM); $j ++){
        if ( defined ($contactmapb[$i][$j]) eq 1  ){

#sidechain
	 $m = $convertindex[$i][1] + 1;
	 $n = $convertindex[$j][1] + 1;
	$attr= $n*($n-1)/2;
 	$attr += $m;
     
	 $bond= $godistb[$i][$j];
	 $bond6 = $bond**6;
	 $bond12 = $bond6**2;
	 $BondEn= $hammap{$restype{$i}}{$restype{$j}}; 

        $BondEn=$GoEn;
	$attcn2[$attr]=2*$BondEn*$bond6; 
	$repcn1[$attr]=$BondEn*$bond12; 

	$negen += $BondEn; 

 if ($i<$j){print "Neg Fortran index", $i," ",$j," ",$bond," ",$restype{$i}," ",$restype{$j}, $BondEn ," ",$negen," ","cn2 ", $attcn2[$attr], "\n"; }


      } #end of if

        if ( defined ($contactmaph[$i][$j]) eq 1 ){
# go hb
	 $m = $convertindex[$i][0] + 1;
	 $n = $convertindex[$j][0] + 1;
	$attr= $n*($n-1)/2;
 	$attr += $m;

	 $bond= $godisth[$i][$j];
	 $bond6 = $bond ** 6;
	 $bond12 = $bond6 ** 2;
	 $BondEn= $GoEn; 
	$attcn2[$attr]=2*$BondEn*$bond6; 
	$repcn1[$attr]=$BondEn*$bond12; 
 if ($i<$j){print "HB Fortran index", $i," ",$j," ",$bond," ",$restype{$i}," ",$restype{$j}, $BondEn ,"\n"; }

	} # end of if
	} # end of for
       	 } #end of for
}
sub calc_go_Go {

my $posen = 0;
my $negen = 0;
my $attr = '';
my $count=0; 
my $BondEn;
my $m = '';
my $n = '';
my $bond = '';
my $bond6 = '';
my $bond12 = '';
my $i = '';
my $j = '';
# evaluate the coefficients for 6-12 potential on specific pairwise go contacts	# get mapping fortran index
  
print "\n"; 

       for (  $i = 0 ; $i < ($RESNUM) ; $i ++)
       { for ( $j = 0 ; $j < ($RESNUM); $j ++){
        if ( defined ($contactmapb[$i][$j]) eq 1  ){

#sidechain
	 $m = $convertindex[$i][1] + 1;
	 $n = $convertindex[$j][1] + 1;
	$attr= $n*($n-1)/2;
 	$attr += $m;
     
	 $bond= $godistb[$i][$j];
	 $bond6 = $bond**6;
	 $bond12 = $bond6**2;
	 $BondEn= $hammap{$restype{$i}}{$restype{$j}}; 

	if($BondEn le 0) { 
	# amber only reads in + value
	$BondEn = abs ($BondEn);
	$BondEn = $GoEn + $BondEn;
	$attcn2[$attr]=2*$BondEn*$bond6; 
	$repcn1[$attr]=$BondEn*$bond12; 

	$negen += $BondEn; 

 if ($i<$j){print "Neg Fortran index", $i," ",$j," ",$bond," ",$restype{$i}," ",$restype{$j}, $BondEn ," ",$negen," ","cn2 ", $attcn2[$attr], "\n"; }

		}
	else{
	
	$BondEn = $GoEn - $BondEn;
        $attcn2[$attr]= 2*$BondEn*$bond6;
        $repcn1[$attr]=$BondEn*$bond12;

	$posen += $BondEn; 

 if ($i<$j){print "Pos Fortran index", $i," ",$j," ",$bond," ",$restype{$i}," ",$restype{$j}, $BondEn ," ",$posen," ","cn2 ", $attcn2[$attr],"\n"; }

	   }

      } #end of if

        if ( defined ($contactmaph[$i][$j]) eq 1 ){
# go hb
	 $m = $convertindex[$i][0] + 1;
	 $n = $convertindex[$j][0] + 1;
	$attr= $n*($n-1)/2;
 	$attr += $m;

	 $bond= $godisth[$i][$j];
	 $bond6 = $bond ** 6;
	 $bond12 = $bond6 ** 2;
	 $BondEn= $GoEn; 
	$attcn2[$attr]=2*$BondEn*$bond6; 
	$repcn1[$attr]=$BondEn*$bond12; 
 if ($i<$j){print "HB Fortran index", $i," ",$j," ",$bond," ",$restype{$i}," ",$restype{$j}, $BondEn ,"\n"; }

	} # end of if
	} # end of for
       	 } #end of for
}

sub calc_nn_homo_Go {

my $attr;
my $bond;
my $bond6;
my $bond12;
my $r1;
my $r2;
my  $negen=0;
my  $posen=0;
my  $BondEn  = 0;
#fortran index

#cb-cb
#include seq dep interaction
	for(my $i=0;$i<($RESNUM);$i++){
	for(my $j=$i+2;$j<($RESNUM);$j++)
	{
#	if (!defined $contactmapb[$i][$j]   )
	{
	#fortran index;
	$r1= $convertindex[$i][1]+1;
	$r2= $convertindex[$j][1]+1;

	$attr= $r2*($r2-1)/2;
	$attr += $r1;
     
	$bond= $ressize[$i]+$ressize[$j];
	$bond=$bond*$FACTOR*$SCALE;
	print "cb-cb ",$attr," ", $i," ",$r1," ", $j, " ",$r2," ",  $bond, "\n"; 

	 $bond6 = $bond**6;
	 $bond12 = $bond6**2;

	print OUTNSC "sc $i $j \n ";

        $BondEn=$GoEn;
	$repcn1[$attr]=$BondEn*$bond12; 
	$attcn2[$attr]=0;

	$negen += $BondEn; 

 if ($i<$j){print "NN Neg Fortran index", $i," ",$j," ",$bond," ",$restype{$i}," ",$restype{$j}, $BondEn ," ",$negen," ","cn2 ", $attcn2[$attr], "\n"; }
	}#end of not defined
	}}

#ca-ca
#add non-specific pairwise interactions.

	for(my $i=0;$i<($RESNUM);$i++){
	for(my $j=$i+4;$j<($RESNUM);$j++)
	{
#	if (!defined $contactmaph[$i][$j]   )
	{
	#fortran index;
	$r1= $convertindex[$i][0]+1;
	$r2= $convertindex[$j][0]+1;

	$attr= $r2*($r2-1)/2;
	$attr += $r1;
	$bond= 1.0*$SCALE*$FACTOR; #calpha-calpha 1.22 (4.65AA/3.8)
         $bond6 = $bond**6;
         $bond12 = $bond6**2;

	print "aa " , $attr," ", $i," ",$r1," ", $j, " ",$r2," ",  $bond, "\n"; 

	my $BondEn=$GoEn;

        $repcn1[$attr]=$BondEn*$bond12;
        $attcn2[$attr]=0;
 if ($i<$j){print "nn HB Fortran index", $i," ",$j," ",$bond," ",$restype{$i}," ",$restype{$j}, $BondEn ,"\n"; }

	print OUTNHB "h $i $j\n";
	}#end of undefined
	}}

#ca-cb
	for(my $i=0;$i<($RESNUM);$i++){
	for(my $j=$i+2;$j<($RESNUM);$j++)
	{
	#fortran index;
	$r1= $convertindex[$i][0]+1;
	$r2= $convertindex[$j][1]+1;

	$attr= $r2*($r2-1)/2;
	$attr += $r1;

	$bond= 0.5+$ressize[$j];
	$bond=$bond*$FACTOR*$SCALE;
#	print $attr," ", $i," ",$r1," ", $j, " ",$r2," ",  $bond, "\n"; 

	my	 $BondEn=$GoEn;
		$bond = $bond ** 12;
		$repcn1[$attr]=$BondEn*$bond; 
		$attcn2[$attr]=0;    
	}}

#cb-ca
	for(my $i=0;$i<($RESNUM);$i++){
	for(my $j=$i+2;$j<($RESNUM);$j++)
	{
	#fortran index;
	$r1= $convertindex[$i][1]+1;
	$r2= $convertindex[$j][0]+1;

	$attr= $r2*($r2-1)/2;
	$attr += $r1;

	$bond= 0.5+$ressize[$i];
	$bond=$bond*$FACTOR*$SCALE;
#	print $attr," ", $i," ",$r1," ", $j, " ",$r2," ",  $bond, "\n"; 

		my $BondEn=$GoEn;
		$bond = $bond ** 12;
		$repcn1[$attr]=$BondEn*$bond; $cn1[$attr]=1;
		$attcn2[$attr]=0;    
	}}
}
sub calc_nn_BT {

my $attr;
my $bond;
my $bond6;
my $bond12;
my $r1;
my $r2;
my  $negen=0;
my  $posen=0;
my  $BondEn  = 0;
#fortran index

#cb-cb
#include seq dep interaction
	for(my $i=0;$i<($RESNUM);$i++){
	for(my $j=$i+2;$j<($RESNUM);$j++)
	{
#	if (!defined $contactmapb[$i][$j]   )
	{
	#fortran index;
	$r1= $convertindex[$i][1]+1;
	$r2= $convertindex[$j][1]+1;

	$attr= $r2*($r2-1)/2;
	$attr += $r1;
     
	$bond= $ressize[$i]+$ressize[$j];
	$bond=$bond*$FACTOR*$SCALE;
	print "cb-cb ",$attr," ", $i," ",$r1," ", $j, " ",$r2," ",  $bond, "\n"; 

	 $bond6 = $bond**6;
	 $bond12 = $bond6**2;
	 $BondEn= $hammap{$restype{$i}}{$restype{$j}}; 

	print OUTNSC "sc $i $j \n ";

	if($BondEn le 0) { 
	# amber only reads in + value
	$BondEn = abs ($BondEn);
	$BondEn = $GoEn + $BondEn;
	$attcn2[$attr]=2*$BondEn*$bond6; 
	$repcn1[$attr]=$BondEn*$bond12; 

	$negen += $BondEn; 

 if ($i<$j){print "NN Neg Fortran index", $i," ",$j," ",$bond," ",$restype{$i}," ",$restype{$j}, $BondEn ," ",$negen," ","cn2 ", $attcn2[$attr], "\n"; }
		}
	else{
	
	$BondEn = $GoEn - $BondEn;
        $attcn2[$attr]= 2*$BondEn*$bond6;
        $repcn1[$attr]=$BondEn*$bond12;

	$posen += $BondEn; 

 if ($i<$j){print "NN Pos Fortran index", $i," ",$j," ",$bond," ",$restype{$i}," ",$restype{$j}, $BondEn ," ",$posen," ","cn2 ", $attcn2[$attr],"\n"; }

	   }

	}#end of not defined
	}}

#ca-ca
#add non-specific pairwise interactions.

	for(my $i=0;$i<($RESNUM);$i++){
	for(my $j=$i+4;$j<($RESNUM);$j++)
	{
#	if (!defined $contactmaph[$i][$j]   )
	{
	#fortran index;
	$r1= $convertindex[$i][0]+1;
	$r2= $convertindex[$j][0]+1;

	$attr= $r2*($r2-1)/2;
	$attr += $r1;
	$bond= 1.22;
         $bond6 = $bond**6;
         $bond12 = $bond6**2;

	print "aa " , $attr," ", $i," ",$r1," ", $j, " ",$r2," ",  $bond, "\n"; 

	my $BondEn=$GoEn;

        $attcn2[$attr]=2*$BondEn*$bond6;
        $repcn1[$attr]=$BondEn*$bond12;
 if ($i<$j){print "nn HB Fortran index", $i," ",$j," ",$bond," ",$restype{$i}," ",$restype{$j}, $BondEn ,"\n"; }

	print OUTNHB "h $i $j\n";
	}#end of undefined
	}}

#ca-cb
	for(my $i=0;$i<($RESNUM);$i++){
	for(my $j=$i+2;$j<($RESNUM);$j++)
	{
	#fortran index;
	$r1= $convertindex[$i][0]+1;
	$r2= $convertindex[$j][1]+1;

	$attr= $r2*($r2-1)/2;
	$attr += $r1;

	$bond= 0.5+$ressize[$j];
	$bond=$bond*$FACTOR*$SCALE;
#	print $attr," ", $i," ",$r1," ", $j, " ",$r2," ",  $bond, "\n"; 

	my	 $BondEn=$GoEn;
		$bond = $bond ** 12;
		$repcn1[$attr]=$BondEn*$bond; 
		$attcn2[$attr]=0;    
	}}

#cb-ca
	for(my $i=0;$i<($RESNUM);$i++){
	for(my $j=$i+2;$j<($RESNUM);$j++)
	{
	#fortran index;
	$r1= $convertindex[$i][1]+1;
	$r2= $convertindex[$j][0]+1;

	$attr= $r2*($r2-1)/2;
	$attr += $r1;

	$bond= 0.5+$ressize[$i];
	$bond=$bond*$FACTOR*$SCALE;
#	print $attr," ", $i," ",$r1," ", $j, " ",$r2," ",  $bond, "\n"; 

		my $BondEn=$GoEn;
		$bond = $bond ** 12;
		$repcn1[$attr]=$BondEn*$bond; $cn1[$attr]=1;
		$attcn2[$attr]=0;    
	}}
}
sub calc_go_Go_urea {

my $posen = 0;
my $negen = 0;
my $attr = '';
my $count=0; 
my $BondEn;
my $m = '';
my $n = '';
my $bond = '';
my $bond6 = '';
my $bond12 = '';
my $i = '';
my $j = '';
# evaluate the coefficients for 6-12 potential on specific pairwise go contacts	# get mapping fortran index
  
print "\n"; 

       for (  $i = 0 ; $i < ($RESNUM) ; $i ++)
       { for ( $j = 0 ; $j < ($RESNUM); $j ++){
        if ( defined ($contactmapb[$i][$j]) eq 1  ){

#sidechain
	 $m = $convertindex[$i][1] + 1;
	 $n = $convertindex[$j][1] + 1;
	$attr= $n*($n-1)/2;
 	$attr += $m;
     
	 $bond= $godistb[$i][$j];
	 $bond6 = $bond**6;
	 $bond12 = $bond6**2;
	 $BondEn= $hammap{$restype{$i}}{$restype{$j}}; 

	# amber only reads in + value
	$BondEn = abs ($BondEn);
	$attcn2[$attr]=2*$BondEn*$bond6; 
	$repcn1[$attr]=$BondEn*$bond12; 

	$negen += $BondEn; 

 if ($i<$j){print "Neg Fortran index", $i," ",$j," ",$bond," ",$restype{$i}," ",$restype{$j}, $BondEn ," ",$negen," ","cn2 ", $attcn2[$attr], "\n"; }

      } #end of if

        if ( defined ($contactmaph[$i][$j]) eq 1 ){
# go hb
	 $m = $convertindex[$i][0] + 1;
	 $n = $convertindex[$j][0] + 1;
	$attr= $n*($n-1)/2;
 	$attr += $m;

	 $bond= $godisth[$i][$j];
	 $bond6 = $bond ** 6;
	 $bond12 = $bond6 ** 2;
	 $BondEn= $GoEn; 
	$attcn2[$attr]=2*$BondEn*$bond6; 
	$repcn1[$attr]=$BondEn*$bond12; 
 if ($i<$j){print "HB Fortran index", $i," ",$j," ",$bond," ",$restype{$i}," ",$restype{$j}, $BondEn ,"\n"; }

	} # end of if
	} # end of for
       	 } #end of for
}
sub calc_nn_BT_urea {

my $attr;
my $bond;
my $bond6;
my $bond12;
my $r1;
my $r2;
my  $negen=0;
my  $posen=0;
my  $BondEn  = 0;
#fortran index

#cb-cb
#include seq dep interaction
	for(my $i=0;$i<($RESNUM);$i++){
	for(my $j=$i+2;$j<($RESNUM);$j++)
	{
#	if (!defined $contactmapb[$i][$j]   )
	{
	#fortran index;
	$r1= $convertindex[$i][1]+1;
	$r2= $convertindex[$j][1]+1;

	$attr= $r2*($r2-1)/2;
	$attr += $r1;
     
	$bond= $ressize[$i]+$ressize[$j];
	$bond=$bond*$FACTOR*$SCALE;
	print "cb-cb ",$attr," ", $i," ",$r1," ", $j, " ",$r2," ",  $bond, "\n"; 

	 $bond6 = $bond**6;
	 $bond12 = $bond6**2;
	 $BondEn= $hammap{$restype{$i}}{$restype{$j}}; 

	print OUTNSC "sc $i $j \n ";

	# amber only reads in + value
	$BondEn = abs ($BondEn);
	$attcn2[$attr]=2*$BondEn*$bond6; 
	$repcn1[$attr]=$BondEn*$bond12; 

	$negen += $BondEn; 

 if ($i<$j){print "NN Neg Fortran index", $i," ",$j," ",$bond," ",$restype{$i}," ",$restype{$j}, $BondEn ," ",$negen," ","cn2 ", $attcn2[$attr], "\n"; }
	}#end of not defined
	}}

#ca-ca
#add non-specific pairwise interactions.

	for(my $i=0;$i<($RESNUM);$i++){
	for(my $j=$i+4;$j<($RESNUM);$j++)
	{
#	if (!defined $contactmaph[$i][$j]   )
	{
	#fortran index;
	$r1= $convertindex[$i][0]+1;
	$r2= $convertindex[$j][0]+1;

	$attr= $r2*($r2-1)/2;
	$attr += $r1;
	$bond= 1.22;
         $bond6 = $bond**6;
         $bond12 = $bond6**2;

	print "aa " , $attr," ", $i," ",$r1," ", $j, " ",$r2," ",  $bond, "\n"; 

	my $BondEn=$GoEn;

        $attcn2[$attr]=2*$BondEn*$bond6;
        $repcn1[$attr]=$BondEn*$bond12;
 if ($i<$j){print "nn HB Fortran index", $i," ",$j," ",$bond," ",$restype{$i}," ",$restype{$j}, $BondEn ,"\n"; }

	print OUTNHB "h $i $j\n";
	}#end of undefined
	}}

#ca-cb
	for(my $i=0;$i<($RESNUM);$i++){
	for(my $j=$i+2;$j<($RESNUM);$j++)
	{
	#fortran index;
	$r1= $convertindex[$i][0]+1;
	$r2= $convertindex[$j][1]+1;

	$attr= $r2*($r2-1)/2;
	$attr += $r1;

	$bond= 0.5+$ressize[$j];
	$bond=$bond*$FACTOR*$SCALE;
#	print $attr," ", $i," ",$r1," ", $j, " ",$r2," ",  $bond, "\n"; 

	my	 $BondEn=$GoEn;
		$bond = $bond ** 12;
		$repcn1[$attr]=$BondEn*$bond; 
		$attcn2[$attr]=0;    
	}}

#cb-ca
	for(my $i=0;$i<($RESNUM);$i++){
	for(my $j=$i+2;$j<($RESNUM);$j++)
	{
	#fortran index;
	$r1= $convertindex[$i][1]+1;
	$r2= $convertindex[$j][0]+1;

	$attr= $r2*($r2-1)/2;
	$attr += $r1;

	$bond= 0.5+$ressize[$i];
	$bond=$bond*$FACTOR*$SCALE;
#	print $attr," ", $i," ",$r1," ", $j, " ",$r2," ",  $bond, "\n"; 

		my $BondEn=$GoEn;
		$bond = $bond ** 12;
		$repcn1[$attr]=$BondEn*$bond; $cn1[$attr]=1;
		$attcn2[$attr]=0;    
	}}
}

sub read_map{
my $count = 0;
my ($type1,$type2,$r1,$r2);
my $i= '';
my $j= '';
my $index = 0;



	while (<INMAP>)
	{
	chomp ();
	$count ++ ;
	($type1,$type2,$r1,$r2) = split (' ',$_);

# r1,r2 is resnum

	if($type2 eq 'b'){ $contactmapb[$r1][$r2] = 1;

#sidechain
#	print "read b", $r1," ",$r2," ", $contactmapb[$r1][$r2],"\n";

	} 
	if($type2 eq 'h'){ 
#hb
	$contactmaph[$r1][$r2] = 1;
#	print "read h", $r1," ",$r2," ", $contactmaph[$r1][$r2],"\n";

}


	}
	
#	for ( my $i = 0 ; $i < $count ; $i ++)
#	{ for (my $j = 0 ; $j < $count ; $j ++){
#	if ( defined $contactmap[$i][$j] > 0 ){ print $i," ",$j," ",$contactmap[$i][$j] ; }
#	}}


	if($count ne $MAXCONTACT) {die "$count ne MAXCONTACT\n";}
}

sub read_ham
{
my @list = ();
my ($res1,$res2,$ener);

	while (<INHAM>)
	{
	chomp ($_);
	($res1, $res2, $ener) = split (' ',$_); 		
	$hammap{$res1}{$res2} = $ener;	
	$hammap{$res2}{$res1} = $ener;	
	}

        foreach my $r1 (sort keys %hammap){
	 foreach my $r2 (keys %{$hammap{$r1}}){
#	print $r1," ",$r2, " ",$hammap{$r1}{$r2},"\n";

	}}

}

sub read_size
{

my @list = ();
my $cindex = 0;
my $resnum = 0;
my $num = 0;

# cindex
	while (<INSIZE>)
	{

	@list = split (' ', $_);
	push @ressize, $list[2]; 	
	$restype{$list[0]} = lc $list[1];


	$num = 2;
	if( $list[1] eq 'GLY') {$num = 1;}
		for (my $i = 0 ; $i < $num; $i ++)
		{
		$convertindex[$resnum][$i]= $cindex;
	        print $restype{$list[0]}," ",$list[0]," ",$list[1], " ",$cindex,"\n";
		$cindex ++ ;	
		}	
	if( $list[1] eq 'GLY') 
	{$convertindex[$resnum][1]= $convertindex[$resnum][0];}

	$resnum ++ ;	

	}


#	print @ressize;

}
sub read_crd
{

my $count = 0;
my $in = '';
my @list = ();

$in = <INCRD>; #first line
chomp($_ =  <INCRD>); #second line

 s/ //g;

        if($_ ne $ATOM){ print "$_ is not eq to $ATOM";}

        while (<INCRD>)
        {
        chomp ($_);
	@list = split ;
	push @coorx, $list[0],$list[3];	
	push @coory, $list[1],$list[4];	
	push @coorz, $list[2],$list[5];	
	}
#	print @coorx;
}

sub calc_xy
{
my @crowd = ();
my $attr = '';
my $i = '';
my $a1 = '';
my $a2 = '';
my $bond = '';
my $BondEn='';
my $n = '';
my $m= '';
my $screen= 1.0;

for(my $i=0;$i<$crowdtype;$i++)
{
	$crowd[$i] = $ATOM + $i + 1;
}

#crowder and crowder
for(my $i=0;$i<$crowdtype;$i++)
    {
        for(my $j=$i;$j<$crowdtype;$j++)
        {
        #$n>$m

        $m=$crowd[$i];
        $n=$crowd[$j];
        $attr= $n*($n-1)/2;
        $attr += $m;

        $bond = ($RG+$RG)*$screen;
        $a = $bond;
                $BondEn=$GoEn;
                $bond = $bond ** 12;
                $repcn1[$attr]=$BondEn*$bond;
                $attcn2[$attr]=0;
        print "xx $n, yy $m, attr = $attr, bond=$a, repcn1= $repcn1[$attr]\n";
        }
     }

#crowder and protein
for ( my $i = 0 ; $i < ($RESNUM) ; $i ++)
        {
          for( my $j=0;$j<$crowdtype;$j++)
            {

#fortran index
$a1= $convertindex[$i][0]+1; #ca
$a2= $convertindex[$i][1]+1; #cb

# ca to xx 
                $n=$crowd[$j];
                $attr= $n*($n-1)/2;
                $attr += $a1;
                $bond=(0.5+$RG)*$screen;
                $BondEn=$GoEn;
                $bond = $bond ** 12;
                $repcn1[$attr]=$BondEn*$bond;
                $attcn2[$attr]=0;

 print "ca $i fortranInd $a1, xx $n, attr = $attr, repcn1= $repcn1[$attr]\n";
        unless ($restype{$i} eq 'gly')
        {
# cb to xx 
                $bond = $ressize[$i]+$RG;
                $bond=$bond*$screen;
                $attr= $n*($n-1)/2;
                $attr += $a2;
                $BondEn=$GoEn;
                $bond = $bond ** 12;
                $repcn1[$attr]=$BondEn*$bond;
                $attcn2[$attr]=0;
 print "cb $i fortranInd $a2, xx $j attr = $attr , repcn1=  $repcn1[$attr]\n";

        }

        } #end of i crowder
   } 

}
