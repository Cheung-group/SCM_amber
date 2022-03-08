This readme is about the information and installation of this modified version of amber10. 
Please start from the original version of amber/10 before compile this version.

***If you want to trace the differences compared with the original AMBER10 code, type
   grep Antonios *.f
   grep Antonios *.h
   grep Qian *.f
   grep Qian *.h
   grep Pengzhi *.f
   grep Pengzhi *.h
   grep Jianfa *.f
   grep Jianfa *.h


***According to License agreement, we only upload the sander source code. Place those sander_source directories in the $AMBERHOME/src and compile.



***Description***

Below are descriptions about versions of sander codes for various types of the Sidechain-Calpha Model (SCM).
   1) sander_source.pseudo.triple --------------- thermodynamics, bulk, without debye-huckel potential
   2) sander_source.pseudo.triple.debye --------- thermodynamics, bulk, with debye-huckel potential
		* The main change is in "short_ene.f"
		* To use this code, you need another input file named debye.inp, there are two numbers in the file in one line:
	  		The first is the relative dielectric constant, eg, 80.0
	  		The second is the ionic strength in M, eg, 0.1
   3) sander_source.pseudo.triple.crowd --------- thermodynamics, crowder (only suitable for ficoll), without debye-huckel potential
		* It is only suitable for ficoll because in "short_ene.f", NPP  = ntypes -1, this "-1" indicates that there is only one type of crowders
   4) sander_source.pseudo.triple.debye.crowd --- thermodynamics, crowder (only suitable for ficoll), with debye-huckel potential
   5) sander_source.pseudo.triple.dumbbell ----- thermodynamics, dextran model, without debye-huckel potential
		* In "short_ene.f", NPP  = ntypes - 2
		* We also need to record the bond energy of dextran, check the changes in "ene.f" and "force.f"
		* The bond energy of dextran will be printed out in md**en, line L3, the third column
   6) sander_source.pseudo.triple.HS ------ thermodynamics, Hard Sphere model for cytoplasm (for azurin), without debye-huckel potential 
		* In "short_ene.f", NPP  = ntypes - 50
   7) sander_source.pseudo.triple.PD ------ thermodynamics, PolyDisperse model for cytoplasm (for azurin), without debye-huckel potential
		* In "short_ene.f", NPP  = ntypes - 82
		* We also need to record the bond, angle, dihedral energy of crowders, check the changes in "ene.f" and "force.f"
		* The bond, angle, dihedral energy of crowders will be printed out in md**en, line L3, the third, fourth, fifth column
   8) sander_source.pseudo.triple.gamma.Brown.end --- protein folding kinetics, bulk, without debye-huckel potential
		* Input files required: prmtop, debye.inp, pseudo.inp, triple.inp, olap.inp, gamma.inp
   9) sander_source.pseudo.triple.gamma.Brown.end.debye.associate --- protein association kinetics, bulk, with debye-huckel potential
		* The main change is in "runmd.f"
		* Input files required: prmtop, debye.inp, pseudo.inp, triple.inp, olap.inp, gamma.inp, associate.inp
  10) sander_source.pseudo.triple.crowd.gamma.Brown.end --- protein folding kinetics, crowder (only suitable for ficoll), without debye-huckel potential
  11) sander_source.pseudo.triple.crowd.gamma.Brown.end.debye --- protein folding kinetics, crowder (only suitable for ficoll), with debye-huckel potential 
  12) sander_source.pseudo.triple.crowd.gamma.Brown.end.debye.associate --- protein association kinetics, crowder (only suitable for ficoll), with debye-huckel potential

*** sometimes you may see the following error when running MPI version: 
   Assertion 'ier == 0' failed in sander.f
   When this happens, it is often related to a memory problem. Check locmem.f, there are two commands:
   if( numtasks <= 8 ) maxpr = maxpr/numtasks
   if( numtasks >  8 ) maxpr = 4*maxpr/(3*numtasks)
   play with the coefficient, eg. change 4 to 40, it can solve the problem



***Installation***
1) cd amber10
   patch -p0 -N -r patch-rejects < bugfix.all

2) Install the AmberTool
   cd amber10/src
   ./configure_at gcc
   make -f Makefile_at

3) Install the serial version first
   cd amber10/src
   ./configure_amber ifort
   make

4) After this is done, in order to compile our in-house code one needs to enter, eg. 
   amber10/src/sander_source.pseudo.triple
   and "make" there.

5) Installation for serial code is done
   No matter which code you compile, the name of the executable file is "sander" in the directory amber10/exe, change it to another name, eg.
   cp amber10/exe/sander amber10/sander.bulk

6) In order to install the mpi version:
	In the directory amber10.parallel/src
   	make clean
   	./configure_amber -openmpi ifort
   	make parallel

9) After this is done, in order to compile our in-house code one needs to enter, eg.
   amber10.parallel/src/sander_source.pseudo.triple
   and "make parallel" there.

10) Installation for MPI code is done
   No matter which code you compile, the name of the executable file is "sander.MPI" in the directory amber10.parallel/exe,
   change it to another name, eg., 
   cp amber10.parallel/exe/sander.MPI amber10.parallel/sander.bulk.MPI
