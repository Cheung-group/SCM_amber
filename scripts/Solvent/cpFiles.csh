#! /bin/csh -f

if ($#argv != 2) then
        echo "Usage:   		./cpFiles.csh  name dir"
        echo "Example: 		./cpFiles.csh cfd ../One_Step.EX"
	goto done
    endif

set name = $argv[1]
set dir = $argv[2]

###################   (1) DEFINITIONS  ###########################

set PDB = $name.pdb
set MINIZ = $name.min.crd
set ATOM = INDEX.$name.ATOM
set PAIR = INDEX.$name.PAIR
set HB = INDEX.$name.HB
set NUM = $name.num

###################   (2) COPYInG  ###########################

#copying files from the bulk directory

cp $dir/$PDB .
cp $dir/$MINIZ .
cp $dir/$ATOM .
cp $dir/$PAIR .
cp $dir/$HB .
cp $dir/$NUM .

###################   EXIT  ###########################
done:
     exit 0
