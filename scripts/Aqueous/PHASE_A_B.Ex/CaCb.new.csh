#! /bin/csh -f

if ($#argv != 12) then
        echo "Usage:   		./CaCb.csh  pdb       name chain model map disulfile ureatype crowdtype RG MG crowdnum box"
        echo "Example: 		./CaCb.csh  1CFD.pdb  cfd   A    GO    BTmap.dat ds.inp nourea 0 11.485 3250.0 0 300"
	goto done
    endif
set pdb = $argv[1] 
set name = $argv[2]
set chain = $argv[3]
set model = $argv[4]
set map = $argv[5]
set disulfile = $argv[6]
set ureatype = $argv[7]
set crowtype = $argv[8]
set rg = $argv[9]
set mg = $argv[10]
set nc = $argv[11]
set box = $argv[12]


###################   (0.a) DEFINITIONS  ###########################

set PDB = $name.pdb
set NATIVE = $name.ab.crd
set ATOM = INDEX.$name.ATOM
set NUM = $name.num
set HB = INDEX.$name.HB
set SC = INDEX.$name.SC
set PAIR = INDEX.$name.PAIR
set TOP = $name.prmtop
set leapscript = leap_script

###################   (0.b) CLEANING UP PDB FILE  ###########################

./removeH_OT.pl $pdb > $PDB	#removing H atoms and OT's
./numbers.pl $PDB > $NUM      #finding natom ntype
sed -i "s/HSE/HIS/g" $PDB
set natom = `head -1 $NUM | awk '{print $1}'`	
set ntype = `tail -1 $NUM | awk '{print $1}'`	

###################   (1.a) HB NATIVE CONTACTS  ###########################

./dsspcmbi $PDB dssptmp				#using DSSP to generate HB list
./hbond-dssp.pl dssptmp > tmplist		#extracting HB pairs
./clean_HBlist.pl tmplist > HBondlist		#Cleaning the list with two rules
						# 1) remove duplicated entry. 2) ensure only a maximum of *2* HB pairs per residue. 
						# If there are more than two, delete the pair that has least bonding energy.

awk '{print NR,"h",$1-1,$2-1}' HBondlist > $HB


###################  (1.b) CLEANING  ###########################
/bin/rm dssptmp tmplist HBondlist


###################   (2.a)) SC NATIVE CONTACTS  ###########################

./makeSUN 	#compile resc.c

set Max = $natom
@ Max = $Max + 1

@ count = 1 

echo END > tmp
cat $PDB tmp > tmpPDB
while ($count < $Max)
./resc tmpPDB $count $chain  >> a
@ count = $count + 1
end
perl -wnle 'BEGIN {$p=0}; $p=1 if /Table II\b/; $p=0 if /Table III/; print $_ if $p' a >> b
perl -wnle 'BEGIN {$p=0}; $p=1 if /Table III\b/; $p=0 if /Table IV/; print $_ if $p' a >> c
perl -wnle 'BEGIN {$p=0}; $p=1 if /Table IV\b/; $p=0 if /Table I/; print $_ if $p' a >> d

./makemap.ab.pl a

awk '{print NR, "b", $3-1,$4-1}'  MapSC > tmpSC
./removeDScontact.pl  $disulfile tmpSC > $SC
/bin/rm tmpSC

###################   (2.b) JOINING ALL CONTACT PAIRS   ###########################


cat $SC $HB | awk '{print NR,$2,$3,$4}' > $PAIR


###################   (2.c) UPDATING HEADER FILES   ###########################

set pair = `wc -l $PAIR | awk '{print $1}'`	
set hb = `wc -l $HB | awk '{print $1}'`	

sed "s/define MAXLENGTH  NTYPE/define MAXLENGTH  $ntype/g;s/define NATOM  NATOM/define NATOM  $natom/g;s/define NTYPE  NTYPE/define NTYPE  $ntype/g;s/define MAXPAIR 500/define MAXPAIR $pair/g;s/define MAXHBPAIR 200/define MAXHBPAIR $hb/g"  IOab.header > IOab.h
sed "s/define MAXLENGTH NTYPE/define MAXLENGTH $ntype/g;s/define NATOM NATOM/define NATOM $natom/g;s/define NTYPE NTYPE/define NTYPE $ntype/g;s/define MAXPAIR 500/define MAXPAIR $pair/g;s/define MAXHBPAIR 200/define MAXHBPAIR $hb/g"  IOab.1.header > IOab.1.h
sed "s/define MAXLENGTH NTYPE/define MAXLENGTH $ntype/g;s/define NATOM NATOM/define NATOM $natom/g;s/define NTYPE NTYPE/define NTYPE $ntype/g;s/define MAXPAIR NPAIR/define MAXPAIR $pair/g;s/define MAXHBPAIR HBPAIR/define MAXHBPAIR $hb/g;s/define UREA UREA/define UREA $ureatype/g;s/define RG RG/define RG $rg/g;s/define MG MG/define MG $mg/g"  IOablinux.header > IOablinux.h

###################  (2.d) CLEANING  ###########################
/bin/rm a b c d MapSC tmp tmpPDB

###################   (3.a) COMPILING beta_coor.c  ##########################

cc beta_coor.c -lm -w -o beta

###################   (3.b) CREATING COARSE GRAINED CRD & INDEX files  ###########################

echo $PDB > input.beta
./beta < input.beta > output.beta

echo $name > HEAD
echo $ntype >> HEAD
awk '/./ {print $1/3.8, $2/3.8, $3/3.8}'  < $PDB.beta > b.crd
cat HEAD b.crd > a.crd

#formating the crd to AMBER crd
cc format.c -lm -w -o format
./format a.crd > $NATIVE

awk '{print $1,$2, $3/3.8}' $PDB.atom > $ATOM

###################  (3.c) CLEANING  ###########################
rm -f a.crd b.crd HEAD input.beta output.beta $PDB.atom $PDB.beta


###################   (4.a) FORCE FIELD  ###########################


cc frcfield_ab.amber10.c -lm -w -o frc	#compiling frcfield
sed 's/\(.*\)/\U\1/' $map > tmpMAP
./frc $NATIVE $ATOM $PAIR tmpMAP $disulfile		#creating the force field file

###################   (4.b) PRMTOP  ###########################

setenv AMBERHOME ~/amber10
./tlp.cacb.pl $PDB $disulfile $crowtype $nc $rg $box > $leapscript 	#creating tleap script
$AMBERHOME/exe/tleap -f $leapscript


###################   (4.c) NONBONDED INTERACTION  ###########################

sed "s/ATOM = TYPE/ATOM = $ntype/g;s/RESNUM = ATOM/RESNUM = $natom/g;s/MAXCONTACT = PAIR/MAXCONTACT = $pair/g;s/RG = RG/RG = $rg/g" calc6-12out.all.parm7.tmp > calc6-12out.all.parm7.pl

chmod +x calc6-12out.all.parm7.pl
./calc6-12out.all.parm7.pl $ATOM $map $PAIR $NATIVE $model $ureatype $crowtype 
sed -i "s/FORMAT/%FORMAT(5E16.8)/g" 6-12out.dat.att
sed -i "s/FORMAT/%FORMAT(5E16.8)/g" 6-12out.dat.rep

sed '/%FLAG LENNARD_JONES_ACOEF/,/%FLAG LENNARD_JONES_BCOEF/c\%FLAG LENNARD_JONES_ACOEF\n%FLAG LENNARD_JONES_BCOEF' cacb.prmtop > 0.top
sed '/%FLAG LENNARD_JONES_BCOEF/,/%FLAG BONDS_INC_HYDROGEN/c\%FLAG LENNARD_JONES_BCOEF\n%FLAG BONDS_INC_HYDROGEN' 0.top > 1.top
sed '/%FLAG LENNARD_JONES_ACOEF/ r  6-12out.dat.rep'  1.top > 0.top
sed '/%FLAG LENNARD_JONES_BCOEF/ r  6-12out.dat.att'  0.top > $TOP

###################   (4.d) PRODUCING psuedo.inp trip.inp  ###########################
if($model == "GO")  then
   sed "s/NONNATIVE XXX/NONNATIVE NO/g" countHB.tmp > countHB.c
endif
if ($model == "BT") then 
   sed "s/NONNATIVE XXX/NONNATIVE YES/g" countHB.tmp > countHB.c 
endif
cc countHB.c -lm -w -o countHB.x	#compile countHB.c
cc triple.c -lm -w -o triple	#compile triple.c
sed 2d $NATIVE > vmd.crd
echo $NATIVE $ATOM $PAIR 1 vmd.crd > countHB.in
echo $NATIVE $ATOM $PAIR 20 > triple.in
./countHB.x < countHB.in
./triple < triple.in

###################   (4.3) CLEANING UP  ###########################

rm -f cacb.prmtop cacb.crd 0.top 1.top
rm -f IO.h frc tmpMAP

###################   EXIT  ###########################
done:
     exit 0
