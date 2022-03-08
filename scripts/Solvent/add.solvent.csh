#! /bin/csh -f

if ($#argv != 12) then
        echo "Usage:   		./add.solvent.csh  pdb       name chain model map disulfile ureatype crowdtype RG MG crowdnum box"
        echo "Example: 		./add.solvent.csh  1CFD.pdb  cfd   A    GO    BTmap.dat ds.inp nourea 1 11.485 3250.0 470 300"
	goto done
    endif
set pdb = $argv[1] 
set name = $argv[2]
set chain = $argv[3]
set model = $argv[4]
set map = $argv[5]
set disulfile = $argv[6]
set ureatype = $argv[7]
set crowdtype = $argv[8]
set rg = $argv[9]
set mg = $argv[10]
set nc = $argv[11]
set box = $argv[12]

###################   (1) DEFINITIONS  ###########################

set PDB = $name.pdb
set ATOM = INDEX.$name.ATOM
set NUM = $name.num
set HB = INDEX.$name.HB
set PAIR = INDEX.$name.PAIR
set TOP = $name.prmtop
set leapscript = leap_script
set MINIZ = $name.min.crd
set OUTCRD = $name.solvent.crd

set natom = `head -1 $NUM | awk '{print $1}'`	
set ntype = `tail -1 $NUM | awk '{print $1}'`	
set pair = `wc -l $PAIR | awk '{print $1}'`	
set hb = `wc -l $HB | awk '{print $1}'`	

###################   (2) FORCE FIELD  ###########################


sed "s/define MAXLENGTH NTYPE/define MAXLENGTH $ntype/g;s/define NATOM NATOM/define NATOM $natom/g;s/define NTYPE NTYPE/define NTYPE $ntype/g;s/define MAXPAIR NPAIR/define MAXPAIR $pair/g;s/define MAXHBPAIR HBPAIR/define MAXHBPAIR $hb/g;s/define UREA UREA/define UREA $ureatype/g;s/define RG RG/define RG $rg/g;s/define MG MG/define MG $mg/g"  IOablinux.header > IOablinux.h
cc frcfield_ab.amber10.c -lm -w -o frc	#compiling frcfield
sed 's/\(.*\)/\U\1/' $map > tmpMAP
./frc $MINIZ $ATOM $PAIR tmpMAP $disulfile		#creating the force field file

###################   (3) PRMTOP  ###########################

setenv AMBERHOME ~/amber10
./tlp.cacb.pl $PDB $disulfile $crowdtype $nc $rg $box > $leapscript 	#creating tleap script
$AMBERHOME/exe/tleap -f $leapscript
./replaceCRD.pl $MINIZ cacb.crd > $OUTCRD

###################   (4) NONBONDED INTERACTION  ###########################

sed "s/ATOM = TYPE/ATOM = $ntype/g;s/RESNUM = ATOM/RESNUM = $natom/g;s/MAXCONTACT = PAIR/MAXCONTACT = $pair/g;s/RG = RG/RG = $rg/g" calc6-12out.all.parm7.tmp > calc6-12out.all.parm7.pl

chmod +x calc6-12out.all.parm7.pl
./calc6-12out.all.parm7.pl $ATOM $map $PAIR $MINIZ $model $ureatype $crowdtype 
sed -i "s/FORMAT/%FORMAT(5E16.8)/g" 6-12out.dat.att
sed -i "s/FORMAT/%FORMAT(5E16.8)/g" 6-12out.dat.rep

sed '/%FLAG LENNARD_JONES_ACOEF/,/%FLAG LENNARD_JONES_BCOEF/c\%FLAG LENNARD_JONES_ACOEF\n%FLAG LENNARD_JONES_BCOEF' cacb.prmtop > 0.top
sed '/%FLAG LENNARD_JONES_BCOEF/,/%FLAG BONDS_INC_HYDROGEN/c\%FLAG LENNARD_JONES_BCOEF\n%FLAG BONDS_INC_HYDROGEN' 0.top > 1.top
sed '/%FLAG LENNARD_JONES_ACOEF/ r  6-12out.dat.rep'  1.top > 0.top
sed '/%FLAG LENNARD_JONES_BCOEF/ r  6-12out.dat.att'  0.top > $TOP

###################   (5) CLEANING UP  ###########################

rm -f cacb.prmtop cacb.crd 0.top 1.top
rm -f IO.h frc tmpMAP

###################   EXIT  ###########################
done:
     exit 0
