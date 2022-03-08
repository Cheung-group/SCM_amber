#!/usr/bin/perl -w
################################################################################
#		A script to build a Calpha-Cbeta model in amber10              #
#		#######      Dirar Homouz Jan 2009        #######              #
#		#######   Antonios Samiotakis 4.13.2010   #######              #
#          #######   Modified to include Disulfid bonds   #######              #
#               #######  Dirar Homouz 7.19.2012     #######                    #
################################################################################

if(@ARGV < 6) {
        print "Usage: tlp.pl file.pdb dsfile crowdtype crowdnum RG box\n";
        exit;
        }

$pdb = shift @ARGV;
$dsfile = shift @ARGV;
$ctype = shift @ARGV;
$NC = shift @ARGV;
$RG = shift @ARGV;
$boxside = shift @ARGV;

open(PDB,"$pdb") or die "failed to open $pdb, ";
open(DSF,"$dsfile") or die "failed to open $dsfile, ";

@cbeta = ();
@array = ();
@ltrca = ();
@ltrcb = ();

for(my $j=0;$j<500;$j++){
  my $i= int $j/36;
  my $ja = ($j%36);
  my $i1= 65 + $i;
  if ($ja<10){$i2 = 48 + $ja;}
  else {$i2 = 65 +$ja -10;}
	
  my $c= sprintf("%c%c",$i1,$i2);
#  printf("%c%c\n",$i1,$i2);
  if($c eq 'DU'){$c = 'Z9';}
  if($c eq 'EP'){$c = 'ZP';}
  push (@ltrca, $c);
}

for(my $j=0;$j<500;$j++){
  my $i= int $j/36;
  my $ja = ($j%36);
  my $i1= 78 + $i;
  if ($ja<10){$i2 = 48 + $ja;}
  else {$i2 = 65 +$ja -10;}

  my $c= sprintf("%c%c",$i1,$i2);
#  printf("%c%c\n",$i1,$i2);
  if($c eq 'DU'){$c = 'Z9';}
  if($c eq 'EP'){$c = 'ZP';}
  push (@ltrcb, $c);
}



printf "cacb_unit = createUnit CACB\n";




my $count = 0;
my $res = 0;
my $atom = 0;

while (<PDB>) {


           $array[$count][0] = substr($_,0,4);
           $array[$count][1] = substr($_,4,7);
           $array[$count][2] = substr($_,11,5);
           $array[$count][3] = substr($_,16,1);
           $array[$count][4] = substr($_,17,3);
           $array[$count][5] = substr($_,20,3);
           $array[$count][6] = substr($_,23,3);
           $array[$count][7] = substr($_,26,12);
           $array[$count][8] = substr($_,38,8);
           $array[$count][9] = substr($_,46,8);
           $array[$count][10] = substr($_,54,6);
           $array[$count][11] = substr($_,60,6);
           $array[$count][12] = substr($_,66,12);
                
#          printf "$array[$count][4]\n";              

        $count++;
        }



#exit;

my $gly_flag = 0;
for ($l = 0; $l < $count; $l++) {
	

	
	   if(($array[$l][2] eq '  CA ') && ($array[$l][4] ne 'GLY')){
			printf "mon$atom = createAtom $ltrca[$res] $ltrca[$res] 0.0\n";
			$x = $array[$l][7]/3.8; $y = $array[$l][8]/3.8; $z=$array[$l][9]/3.8;
			printf "set mon$atom position {$x\t$y\t$z}\n";
			printf "cacb_res$res = createResidue R$res\n";
			printf "add cacb_res$res mon$atom\n";
			if ($atom >0){
					if ($gly_flag eq 1){
							$k = $atom-1;
							printf "bond mon$atom mon$k\n";
							}
					else {
						$k = $atom - 2;
						printf "bond mon$atom mon$k\n";
			                     }
					} 
			$atom = $atom +1;
			$gly_flag = 0;
				}


	    if(($array[$l][2] eq '  CA ') && ($array[$l][4] eq 'GLY')){
                        printf "mon$atom = createAtom $ltrca[$res] $ltrca[$res] 0.0\n";
			$x = $array[$l][7]/3.8; $y = $array[$l][8]/3.8; $z=$array[$l][9]/3.8;
			printf "set mon$atom position {$x\t$y\t$z}\n";
                        printf "cacb_res$res = createResidue R$res\n";
                        printf "add cacb_res$res mon$atom\n";
			printf "add cacb_unit cacb_res$res\n";
			if ($atom >0){
                                        if ($gly_flag eq 1){
                                                        $k = $atom-1;
                                                        printf "bond mon$atom mon$k\n";
                                                        }
                                        else {
                                                $k = $atom - 2;
                                                printf "bond mon$atom mon$k\n";
                                             }
                                        } 
                        $atom = $atom +1;
			$res = $res +1;
			$gly_flag = 1;
                                }


	   if($array[$l][2] eq '  CB '){
			$cbeta[$res] = $atom;
                        printf "mon$atom = createAtom $ltrcb[$res] $ltrcb[$res] 0.0\n";
			$x = $array[$l][7]/3.8; $y = $array[$l][8]/3.8; $z=$array[$l][9]/3.8;
			printf "set mon$atom position {$x\t$y\t$z}\n";
                        printf "add cacb_res$res mon$atom\n";
                        printf "add cacb_unit cacb_res$res\n";
			$k = $atom - 1;
                        printf "bond mon$atom mon$k\n";
			$res = $res + 1;
			$atom = $atom + 1;
                                }
}


#Disulfide bonds
@dsbond = ();
my $nds = <DSF>;
if($nds > 0){
  while(<DSF>){
	@dsbond = split (" ",$_);
	$aa = $cbeta[$dsbond[0]];	
	$bb = $cbeta[$dsbond[1]];	
	printf "bond mon$aa mon$bb\n";
  }
}


#Spherical Crowders
if($ctype ==1){
  for(my $n=0;$n<$NC;$n++){
        printf "cr_mon$n = createAtom ZZ ZZ 0.0\n";
        my $x = $boxside*(rand()-.5);
        my $y = $boxside*(rand()-.5);
        my $z = $boxside*(rand()-.5);
        printf "set cr_mon$n position { $x\t$y\t$z}\n";
        printf "cr_res$n = createResidue CRWD\n";
        printf "add cacb_unit cr_res$n\n";
        #printf "set cacb_unit.$n restype solvent\n";
        printf "add cr_res$n cr_mon$n\n";
  }
printf "set cacb_unit box $boxside\n";
}

#Dumbbell Crowders
if($ctype ==2){
  for(my $n=0;$n<$NC;$n++){
	$atom1 = 2*$n;
	$atom2 = 2*$n + 1;
        printf "cr_mon$atom1 = createAtom ZZ ZZ 0.0\n";
        my $x1 = $boxside*(rand()-.5);
        my $y1 = $boxside*(rand()-.5);
        my $z1 = $boxside*(rand()-.5);
        printf "set cr_mon$atom1 position { $x1\t$y1\t$z1}\n";
        printf "cr_mon$atom2 = createAtom ZY ZY 0.0\n";
        my $x2 = $x1 + $RG*(2.0/sqrt(3.0)+0.01*(rand()-.5));
        my $y2 = $y1 + $RG*(2.0/sqrt(3.0)+0.01*(rand()-.5));
        my $z2 = $z1 + $RG*(2.0/sqrt(3.0)+0.01*(rand()-.5));
        printf "set cr_mon$atom2 position { $x2\t$y2\t$z2}\n";
        printf "cr_res$n = createResidue CRWD\n";
        printf "add cacb_unit cr_res$n\n";
        #printf "set cacb_unit.$n restype solvent\n";
        printf "add cr_res$n cr_mon$atom1\n";
        printf "add cr_res$n cr_mon$atom2\n";
        printf "bond cr_mon$atom1 cr_mon$atom2\n";
  }
  printf "set cacb_unit box $boxside\n";
}



printf "myforce = loadamberparams frc.go\n";
printf "saveAmberParm cacb_unit cacb.prmtop cacb.crd\n";
printf "quit\n";
