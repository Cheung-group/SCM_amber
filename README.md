This readme is about the information and installation of this new version of amber10 (begins from a clean version!!) in XANADU
Qian Wang 7/17/2012

***If you want to trace the differences compared with the original AMBER10 code, type
   grep Antonios *.f
   grep Antonios *.h
   grep Qian *.f
   grep Qian *.h

1. Introduction
2. Installation


1. Introduction
   1) sander_source.pseudo.triple.2012.July --------------- thermodynamics, bulk, without debye-huckel potential
   2) sander_source.pseudo.triple.debye.2012.July --------- thermodynamics, bulk, with debye-huckel potential
	* the main change is in "short_ene.f"
	* to use this code, you need another input file named debye.inp, there are two numbers in the file:
	  the first is the relative dielectric constant, eg, 80.0
	  the second is the ionic strength, eg, 0.1 (M)
   3) sander_source.pseudo.triple.crowd.2012.July --------- thermodynamics, crowder (only suitable for ficoll), without debye-huckel potential
	* the reason for only suitable for ficoll is in "short_ene.f", NPP  = ntypes -1, this "-1" indicates that there is only one type for crowders
   4) sander_source.pseudo.triple.debye.crowd.2012.July --- thermodynamics, crowder (only suitable for ficoll), with debye-huckel potential
   5) sander_source.pseudo.triple.dumbbell.2012.July ----- thermodynamics, dextran model, without debye-huckel potential
	* in "short_ene.f", NPP  = ntypes - 2
	* we also need to record the bond energy of dextran, check the changes in "ene.f" and "force.f"
	* the bond energy of dextran will be printed out in md**en, line L3, the third column
   6) sander_source.pseudo.triple.HS.2012.July ------ thermodynamics, Hard Sphere model for cytoplasm (for azurin), without debye-huckel potential 
	* in "short_ene.f", NPP  = ntypes - 50
   7) sander_source.pseudo.triple.PD.2012.July ------ thermodynamics, PolyDisperse model for cytoplasm (for azurin), without debye-huckel potential
	* in "short_ene.f", NPP  = ntypes - 82
	* we also need to record the bond, angle, dihedral energy of crowders, check the changes in "ene.f" and "force.f"
	* the bond, angle, dihedral energy of crowders will be printed out in md**en, line L3, the third, fourth, fifth column
   8) sander_source.pseudo.triple.gamma.Brown.end.2012.July --- protein folding kinetics, bulk, without debye-huckel potential
	* you need prmtop, debye.inp, pseudo.inp, triple.inp, olap.inp, gamma.inp
   9) sander_source.pseudo.triple.gamma.Brown.end.debye.associate.2012.July  --- protein association kinetics, bulk, with debye-huckel potential
	* the main change is in "runmd.f"
	* you need prmtop, debye.inp, pseudo.inp, triple.inp, olap.inp, gamma.inp, associate.inp
  10) sander_source.pseudo.triple.crowd.gamma.Brown.end.2012.July --- protein folding kinetics, crowder (only suitable for ficoll), without debye-huckel potential
  11) sander_source.pseudo.triple.crowd.gamma.Brown.end.debye.2012.July --- protein folding kinetics, crowder (only suitable for ficoll), with debye-huckel potential 
  12) sander_source.pseudo.triple.crowd.gamma.Brown.end.debye.associate.2012.July --- protein association kinetics, crowder (only suitable for ficoll), with debye-huckel potential

*** sometimes you may see the following error when running MPI version: 
   Assertion 'ier == 0' failed in sander.f
   when this happens, it is often related to memory problem. Check locmem.f, there are two commands:
   if( numtasks <= 8 ) maxpr = maxpr/numtasks
   if( numtasks >  8 ) maxpr = 4*maxpr/(3*numtasks)
   play with the coefficient, eg. change 4 to 40, it can solve the problem

2. Installation

1) After log into XANADU
   module load intel/compiler/64/11.1/073
   module load openmpi/intel/64/1.4.3

2) cd amber10
   patch -p0 -N -r patch-rejects < bugfix.all

3) Install the AmberTool
   cd amber10/src
   ./configure_at gcc
   make -f Makefile_at

4) Install the serial version first
   cd amber10/src
   ./configure_amber ifort
   make

5) After this is done, in order to compile our in-house code one needs to enter, eg. 
   amber10/src/sander_source.pseudo.triple.2012.July
   and "make" there.

6) Installation for serial code is done
   No matter which code you compile, the name of the executable file is "sander" in the directory amber10/exe, change it to another name, eg.
   cp amber10/exe/sander amber10/sander.bulk.2012.July

7) In order to install the mpi version:
   begin from a clean version of amber10 (do not mix with the seiral code, name it differently, eg. amber10.parallel)
   repeat from step 1) to step 6)

8) in the directory amber10.parallel/src
   make clean
   ./configure_amber -openmpi ifort
   make parallel

9) After this is done, in order to compile our in-house code one needs to enter, eg.
   amber10.parallel/src/sander_source.pseudo.triple.2012.July
   and "make parallel" there.

10) Installation for MPI code is done
   No matter which code you compile, the name of the executable file is "sander.MPI" in the directory amber10.parallel/exe,
   change it to another name, eg., 
   cp amber10.parallel/exe/sander.MPI amber10.parallel/sander.bulk.2012.July.MPI
