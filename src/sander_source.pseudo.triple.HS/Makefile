#
#************************************************************************
#                              AMBER                                   **
#                                                                      **
#               Copyright (c) 1986, 1991, 1995, 1997, 1999             **
#                Regents of the University of California               **
#                       All Rights Reserved.                           **
#                                                                      **
#  This software provided pursuant to a license agreement containing   **
#  restrictions on its disclosure, duplication, and use. This software **
#  contains confidential and proprietary information, and may not be   **
#  extracted or distributed, in whole or in part, for any purpose      **
#  whatsoever, without the express written permission of the authors.  **
#  This notice, and the associated author list, must be attached to    **
#  all copies, or extracts, of this software. Any additional           **
#  restrictions set forth in the license agreement also apply to this  **
#  software.                                                           **
#************************************************************************
#
# Makefile for Amber 10   SANDER 
#
include ../config_amber.h

# Sander Fortran source files are free format
LOCALFLAGS = $(FREEFORMAT_FLAG)

SHELL=/bin/sh
VPATH=

NCSU_OBJECTS = mt19937.o ncsu-abmd-hooks.o ncsu-bbmd-ctxt.o ncsu-bbmd-hooks.o \
ncsu-cftree.o ncsu-colvar-math.o ncsu-colvar-type.o ncsu-colvar-utils.o \
ncsu-colvar.o ncsu-constants.o ncsu-cv-ANGLE.o ncsu-cv-DISTANCE.o \
ncsu-cv-HANDEDNESS.o ncsu-cv-MULTI_RMSD.o ncsu-cv-N_OF_BONDS.o \
ncsu-cv-N_OF_STRUCTURES.o ncsu-cv-R_OF_GYRATION.o ncsu-cv-TORSION.o \
ncsu-lexer.o ncsu-parser.o ncsu-pmd-hooks.o ncsu-rmsd.o \
ncsu-sander-hooks.o ncsu-sander-proxy.o ncsu-smd-hooks.o ncsu-umbrella.o \
ncsu-utils.o ncsu-value.o


EVB_BASE = evb_vars.o evb_mathf.o evb_init.o evb_bcast.o \
    evb_dealloc.o check_input.o evb_io.o evb_keywrd.o \
    exchange_gauss.o exchange_warshel.o react_flux_init.o \
    react_path.o react_flux.o kappa_keck.o \
    schlegel_dg.o schlegel_full.o schlegel_diis.o schlegel_gmres.o \
    schlegel_gbasis.o schlegel_pgbasis.o \
    bnd_strtch.o ang_bend.o dihed_tors.o torinternal.o \
    evb_gfchk.o bmat_inv.o bmat_grad.o \
    unitv.o crossv.o cart2internal.o wdc_bmat.o select_coord.o \
    schlegel_uff.o schlegel_poh_uff.o schlegel_poh_irc.o schlegel_hcn_cart.o \
    schlegel_hcn_int_irc.o schlegel_hcn_cart_irc.o evb_outoptions.o

EVB_DEBUG = evb_2stdebug.o egap_umb_2stdebug.o dg_debug.o \
    evb_io_debug.o bmat_debug.o bmat_inv_debug.o evb_io_debug2.o \
    morse_anal2num.o xwarshel_anal2num.o xgauss_anal2num.o \
    vdw_anal2num.o rcdbonds_anal2num.o rcbond_anal2num.o \
    dg_grad_anal2num.o

EVBOBJ =  $(EVB_BASE) $(EVB_DEBUG) evb_input.o evb_alloc.o evb_matrix.o \
    morsify.o mod_vdw.o evb_ntrfc.o evb_force.o out_evb.o \
    dbonds_anal2num.o bond_anal2num.o evb_umb.o evb_umb_primitive.o

EVBPIMD =  $(EVB_BASE) $(EVB_DEBUG) evb_input.LES.o evb_alloc.LES.o evb_matrix.LES.o \
    morsify.LES.o mod_vdw.LES.o evb_ntrfc.LES.o evb_force.LES.o \
    evb_pimd_init.LES.o out_evb.LES.o dbonds_anal2num.o bond_anal2num.o \
    evb_umb.LES.o evb_umb_primitive.LES.o

AMOEBAOBJ= \
    amoeba_valence.o amoeba_multipoles.o amoeba_recip.o amoeba_interface.o\
    amoeba_direct.o amoeba_mdin.o amoeba_adjust.o amoeba_self.o\
    amoeba_vdw.o amoeba_induced.o amoeba_runmd.o 

APBSOBJ= \
	apbs_vars.o apbs.o

QMOBJ= \
    qm_mm.o qm_link_atoms.o \
    qm_ewald.o qm_gb.o qm_zero_charges.o qm_print_info.o \
    qm_assign_atom_types.o

QM2OBJ= \
	qm2_allocate_e_repul.o qm2_calc_charges.o \
	qm2_calc_rij_and_eqns.o qm2_dihed.o qm2_energy.o qm2_fock.o \
	qm2_get_qm_forces.o qm2_get_qmmm_forces.o qm2_h1elec.o qm2_hcore_qmqm.o \
	qm2_hcore_qmmm.o qm2_identify_peptide_links.o \
	qm2_load_params_and_allocate.o qm2_repp.o qm2_rotate_qmqm.o \
	qm2_scf.o qm2_setup_orb_exp.o qm2_smallest_number.o

QMMMOBJ_DFTB= \
	qm2_dftb_broyden.o \
	qm2_dftb_dispersion_egr.o qm2_dftb_dispersion_params.o \
	qm2_dftb_dispersionread.o qm2_dftb_scf.o qm2_dftb_energy.o \
	qm2_dftb_ewevge.o qm2_dftb_externalshift.o \
	qm2_dftb_fermi.o qm2_dftb_get_qm_forces.o qm2_dftb_gamma.o \
	qm2_dftb_gammamat.o qm2_dftb_get_qmmm_forces.o qm2_dftb_gettab.o \
	qm2_dftb_load_params.o \
	qm2_dftb_mulliken.o \
	qm2_dftb_repulsiv.o qm2_dftb_self.o qm2_dftb_shift.o \
	qm2_dftb_short_range.o qm2_dftb_skpar.o qm2_dftb_slkode.o \
	qm2_dftb_slktrafo.o qm2_dftb_read_cm3.o qm2_dftb_cm3.o \
	qm2_dftb_ewald.o qm2_dftb_gb.o 

MMOBJ=  constants.o \
    nose_hoover.o nose_hoover_vars.o nose_hoover_init.o \
    stack.o qmmm_module.o qm2_dftb_module.o trace.o lmod.o decomp.o \
    icosasurf.o egb.o remd.o findmask.o pb_force.o pb_exmol.o pb_mpfrc.o pb_direct.o\
    pb_list.o np_force.o sa_driver.o relax_mat.o nmr.o \
    multisander.o sander.o trajene.o cshf.o nmrcal.o pearsn.o \
    printe.o runmin.o rdparm.o \
    mdread.o locmem.o runmd.o getcor.o degcnt.o decnvh.o \
    fastwt.o parallel.o shake.o ene.o mdwrit.o minrit.o \
    set.o dynlib.o mdfil.o nmlsrc.o \
    ew_force.o ew_setup.o ew_box.o ew_bspline.o\
    ew_fft.o nonbond_list.o short_ene.o ew_recip.o pcshift.o align.o csa.o \
    rfree.o rgroup.o random.o amopen.o \
    debug.o ew_recip_reg.o ew_handle_dips.o ew_dipole_recip.o \
    mexit.o new_time.o extra_pts.o thermo_int.o \
    matinv.o assert.o mmtsb.o mmtsb_client.o erfcfun.o veclib.o \
    is_init.o constantph.o prn_dipoles.o ips.o sglds.o \
    $(AMOEBAOBJ) bintraj.o \
    spatial_recip.o spatial_fft.o parms.o softcore.o molecule.o xref.o dssp.o

FULLPIMDOBJ = pimd_vars.o pimd_force.o pimd_init.o cmd_vars.o cmd_matrix.o 

PARTPIMDOBJ = pimd_vars.o pimd_force.LES.o pimd_init.LES.o cmd_vars.o cmd_matrix.LES.o

PUPILOBJ=  constants.o \
	nose_hoover.o nose_hoover_vars.o nose_hoover_init.o \
	stack.o qmmm_module.o qm2_dftb_module.o trace.o lmod.o decomp.o \
	pupildata.o icosasurf.o egb.o findmask.o pb_force.o pb_exmol.o \
        pb_mpfrc.o pb_direct.o\
	pb_list.o np_force.o sa_driver.o relax_mat.o nmr.o \
	multisander.o sander.PUPIL.o trajene.o cshf.o nmrcal.o pearsn.o \
	printe.PUPIL.o runmin.o force.PUPIL.o rdparm.o  \
	mdread.o locmem.o runmd.o getcor.o degcnt.o decnvh.o \
	fastwt.o parallel.o shake.o ene.o mdwrit.o minrit.o \
	set.o dynlib.PUPIL.o mdfil.PUPIL.o nmlsrc.o \
	ew_force.o ew_setup.o ew_box.o ew_bspline.o\
	ew_fft.o nonbond_list.o short_ene.o ew_recip.o pcshift.o align.o csa.o \
	rfree.o rgroup.o random.o amopen.o \
	debug.o ew_recip_reg.o ew_handle_dips.o ew_dipole_recip.o \
	mexit.PUPIL.o new_time.o extra_pts.o thermo_int.o \
	matinv.o assert.o mmtsb.o mmtsb_client.o erfcfun.o veclib.o \
	is_init.o constantph.o prn_dipoles.o ips.o sglds.o \
	$(AMOEBAOBJ) bintraj.o \
	spatial_recip.o spatial_fft.o parms.o softcore.o molecule.o xref.o \
	putvalues.o fixport.o  dssp.o

LSCIVROBJ = lscivr_vars.o lsc_init.o lsc_xp.o

LESOBJ=             \
    constants.o nose_hoover.o nose_hoover_vars.o nose_hoover_init_LES.o \
    qmmm_module.o qm2_dftb_module.o trace.o lmod.o decomp.o icosasurf.o \
    egb.LES.o remd.LES.o findmask.o pb_force.o pb_exmol.o pb_mpfrc.o \
    pb_direct.o \
    pb_list.o np_force.o sa_driver.o relax_mat.o nmr.LES.o \
    multisander.LES.o sander.LES.o trajene.o cshf.o nmrcal.o pearsn.o \
    printe.o runmin.o nonbond_list.LES.o force.LES.o rdparm.LES.o \
    mdread.LES.o locmem.LES.o runmd.LES.o getcor.o degcnt.LES.o decnvh.o \
    fastwt.o parallel.LES.o shake.o ene.LES.o mdwrit.o minrit.o \
    set.o dynlib.LES.o mdfil.o nmlsrc.o \
    ew_force.LES.o ew_setup.LES.o ew_box.o ew_bspline.o\
    ew_fft.o short_ene.LES.o ew_recip.LES.o pcshift.o \
    align.o csa.o stack.o rfree.o rgroup.o \
    random.o amopen.o \
    debug.o ew_recip_reg.o ew_handle_dips.o ew_dipole_recip.o \
    mexit.o new_time.o extra_pts.LES.o thermo_int.o matinv.o assert.o \
    mmtsb.o mmtsb_client.o erfcfun.o veclib.o \
    is_init.o constantph.o prn_dipoles.o ips.o sglds.o \
    $(QMOBJ) $(QM2OBJ) $(QMMMOBJ_DFTB) $(AMOEBAOBJ) bintraj.o \
    spatial_recip.o spatial_fft.o parms.o softcore.o molecule.o xref.o dssp.o

install:  serial

SERIALPROGS = sander$(SFX) sander.LES$(SFX) ambmask$(SFX)

serial:	checkserconf $(SERIALPROGS)
	/bin/mv $(SERIALPROGS) ../../exe

PARALLELPROGS = sander.MPI$(SFX) sander.LES.MPI$(SFX)

parallel:	checkparconf $(PARALLELPROGS)
	/bin/mv $(PARALLELPROGS) ../../exe

all_serial_programs: $(SERIALPROGS) 
	make sander.PUPIL$(SFX)
	make sander.APBS$(SFX)

#---------------------------------------------------------------------------
sander$(SFX): $(MMOBJ) $(QMOBJ) $(QM2OBJ) $(QMMMOBJ_DFTB) $(FULLPIMDOBJ) \
              $(LSCIVROBJ) qm_div.o force.o syslib \
              nxtsec lapack blas lmod $(INCDIVCON) checkserconf $(NCSU_OBJECTS)
	$(LOAD) -o sander$(SFX)  $(MMOBJ) $(QMOBJ) $(QM2OBJ) $(QMMMOBJ_DFTB) $(FULLPIMDOBJ) \
                $(LSCIVROBJ) qm_div.o force.o \
		../lmod/lmod.a $(LIBDIVCON) $(NCSU_OBJECTS) \
		../lapack/lapack.a ../blas/blas.a \
		../lib/nxtsec.o ../lib/sys.a $(NETCDFLIB) $(LOADLIB)

#---------------------------------------------------------------------------
sander.MPI$(SFX): $(MMOBJ) $(QMOBJ) $(QM2OBJ) $(QMMMOBJ_DFTB) $(FULLPIMDOBJ) \
              $(LSCIVROBJ) $(EVBOBJ) qm_div.o force.o \
              syslib nxtsec lapack blas lmod checkparconf $(NCSU_OBJECTS)
	$(LOAD) -o sander.MPI$(SFX)  $(MMOBJ) $(QMOBJ) $(QM2OBJ)  $(QMMMOBJ_DFTB) \
              $(FULLPIMDOBJ) $(LSCIVROBJ) $(EVBOBJ) qm_div.o force.o \
		../lmod/lmod.a ../lapack/lapack.a ../blas/blas.a \
		../lib/nxtsec.o ../lib/sys.a $(NCSU_OBJECTS) $(NETCDFLIB) $(LOADLIB)

#---------------------------------------------------------------------------
sander.PUPIL$(SFX): $(PUPILOBJ) $(QMOBJ) $(QM2OBJ) $(QMMMOBJ_DFTB) $(FULLPIMDOBJ) \
	$(LSCIVROBJ)  qm_div.o syslib \
	nxtsec lapack blas lmod $(INCDIVCON) checkserconf $(NCSU_OBJECTS)
	$(LOAD) -o sander.PUPIL$(SFX) $(PUPILOBJ) $(QMOBJ) $(QM2OBJ) $(QMMMOBJ_DFTB) $(FULLPIMDOBJ) \
	qm_div.o $(LSCIVROBJ) \
	../lmod/lmod.a $(LIBDIVCON) \
	../lapack/lapack.a ../blas/blas.a \
	../lib/nxtsec.o ../lib/sys.a $(NCSU_OBJECTS) $(NETCDFLIB) $(LOADLIB) \
	-lrt -lm -lc -L${PUPIL_PATH}/lib -lPUPIL -lPUPILBlind 

#---------------------------------------------------------------------------
sander.LES$(SFX): $(LESOBJ) $(PARTPIMDOBJ) syslib nxtsec lapack blas lmod \
		$(LSCIVROBJ) $(NCSU_OBJECTS) checkserconf
	$(LOAD) -o sander.LES$(SFX) $(LESOBJ) $(PARTPIMDOBJ) $(LSCIVROBJ) \
		../lmod/lmod.a ../lapack/lapack.a ../blas/blas.a \
		../lib/nxtsec.o ../lib/sys.a $(NCSU_OBJECTS) $(NETCDFLIB) $(LOADLIB)

#---------------------------------------------------------------------------
sander.LES.MPI$(SFX): $(LESOBJ) $(EVBPIMD) $(PARTPIMDOBJ) syslib nxtsec lapack blas lmod \
		$(LSCIVROBJ) $(NCSU_OBJECTS) checkparconf
	$(LOAD) -o sander.LES.MPI$(SFX) $(LESOBJ) $(EVBPIMD) $(PARTPIMDOBJ) $(LSCIVROBJ) \
		../lmod/lmod.a ../lapack/lapack.a ../blas/blas.a \
		../lib/nxtsec.o ../lib/sys.a $(NCSU_OBJECTS) $(NETCDFLIB) $(LOADLIB)

#---------------------------------------------------------------------------
sander.APBS$(SFX): $(APBSOBJ) $(MMOBJ) $(QMOBJ) $(QM2OBJ) \
              $(LSCIVROBJ) $(QMMMOBJ_DFTB) $(FULLPIMDOBJ) $(NCSU_OBJECTS) \
              qm_div.o force.o syslib \
              nxtsec lapack blas lmod $(INCDIVCON) checkserconf
	$(LOAD) -o sander.APBS$(SFX)  $(APBSOBJ) $(MMOBJ) $(QMOBJ) $(QM2OBJ) \
              $(QMMMOBJ_DFTB) $(FULLPIMDOBJ) $(LSCIVROBJ) \
              qm_div.o force.o \
              -L$(APBS_LIB) -liapbs -lapbsmainroutines -lapbs -lmaloc \
              ../lmod/lmod.a $(LIBDIVCON) \
              ../lapack/lapack.a ../blas/blas.a \
              ../lib/nxtsec.o ../lib/sys.a $(NCSU_OBJECTS) $(NETCDFLIB) $(LOADLIB)

#---------------------------------------------------------------------------
ambmask$(SFX): 	ambmask.o findmask.o amopen.o mexitl nxtsec
	$(LOAD) -o ambmask$(SFX) ambmask.o findmask.o \
		amopen.o ../lib/mexit.o ../lib/nxtsec.o $(LOADLIB)

#-----------LIBS

syslib:
	cd ../lib; make sys.a

lapack:  
	cd ../lapack; make

blas:  
	cd ../blas; make

nxtsec:  
	cd ../lib; make nxtsec.o

mexitl:
	cd ../lib; make mexit.o

lmod:  
	cd ../lmod; make
	cd ../lmod; cp lmodprmtop ../../exe; chmod 0755 ../../exe/lmodprmtop

divcon:	
	cd ../dcqtp; make libdivcon.a

netcdf.mod:
	cd ../netcdf/src; make install
	cp ../netcdf/include/*.mod .
	cp ../netcdf/bin/ncdump$(SFX) ../../exe

clean:
	-/bin/rm $(SERIALPROGS) $(PARALLELPROGS)
	-/bin/rm *.o
	-/bin/rm _*.f
	-/bin/rm *.mod
	-/bin/rm *.d
	-/bin/rm *nbflag
	cd ../dcqtp; make clean

depend::
	./makedepend > depend

checkparconf::
	./checkparconf

checkserconf::
	./checkserconf


#----------HEADER DEPENDENCIES (constructed from ./makedepend > depend)
include depend
