#include "copyright.h"
#include "dprec.h"
#include "assert.h"

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!+ [Enter a one-line description of subroutine rdparm1 here]
subroutine rdparm1(nf)

   use parms, only: numbnd,numang,nptra,nphb,MAX_ATOM_TYPE,MAX_BOND_TYPE, &
                    charmm
   use molecule

   implicit none
   
#  include "md.h"
#  include "memory.h"
#  include "files.h"
#  include "box.h"
#  include "nmr.h"
#  include "extra_pts.h"
#ifdef LES
#  include "les.h"
#endif
   integer nf
   integer i,nspsol,iok
   integer nhparm,idum,nttyp
   integer mbper,mgper,mdper,mbona,mtheta,mphia ! read but ignored
   character(len=80) fmt
   character(len=80) fmtin,ifmt,afmt,rfmt,type, line
   character(len=80), allocatable, dimension(:) :: ffdesc
   integer :: nlines, iscratch, ier
   ifmt = '(12I6)'
   afmt = '(20A4)'
   rfmt = '(5E16.8)'
   
   !     ----- READ THE MOLECULAR TOPOLOGY -----
   
   nspsol = 0
   
   !     ----- FORMATTED INPUT -----
   
   fmtin = '(A80)'
   type = 'TITLE'
   call nxtsec(nf,  6,  0,fmtin,  type,  fmt,  iok)
   read(nf,fmtin) title
   
   fmtin = ifmt
   type = 'POINTERS'
   call nxtsec(nf,  6,  0,fmtin,  type,  fmt,  iok)
   read(nf,fmt) natom,ntypes,nbonh,mbona,ntheth,mtheta,nphih,mphia, &
         nhparm,nparm,nnb,nres,nbona,ntheta,nphia, &
         numbnd,numang,nptra,natyp,nphb,ifpert,nbper,ngper,ndper, &
         mbper,mgper,mdper,ifbox,nmxrs,ifcap,numextra,ncopy
#ifdef LES
   if (nparm /= 1) then
      write(6,*) ' *** THIS VERSION ONLY ACCEPTS TOPOLOGY FILES'
      write(6,*) '     THAT WERE CREATED BY ADDLES, WITH NPARM=1'
      write(6,*) '     USE A VERSION COMPILED WITHOUT -DLES '
      call mexit(6,1)
   end if
#else
   if (nparm == 1) then
      write(6,*) ' *** THIS VERSION WILL NOT ACCEPT TOPOLOGY FILES'
      write(6,*) '     THAT WERE CREATED BY ADDLES, WITH NPARM=1'
      write(6,*) '     USE A VERSION COMPILED WITH -DLES '
      call mexit(6,1)
   end if
#endif
   write(6,8118) natom,ntypes,nbonh,mbona,ntheth,mtheta,nphih,mphia, &
         nhparm,nparm,nnb,nres,nbona,ntheta,nphia,numbnd, &
         numang,nptra,natyp,nphb,ifbox,nmxrs,ifcap,numextra,ncopy
   8118 format(t2, &
         'NATOM  = ',i7,' NTYPES = ',i7,' NBONH = ',i7,' MBONA  = ',i7, &
         /' NTHETH = ',i7,' MTHETA = ',i7,' NPHIH = ',i7,' MPHIA  = ',i7, &
         /' NHPARM = ',i7,' NPARM  = ',i7,' NNB   = ',i7,' NRES   = ',i7, &
         /' NBONA  = ',i7,' NTHETA = ',i7,' NPHIA = ',i7,' NUMBND = ',i7, &
         /' NUMANG = ',i7,' NPTRA  = ',i7,' NATYP = ',i7,' NPHB   = ',i7, &
         /' IFBOX  = ',i7,' NMXRS  = ',i7,' IFCAP = ',i7,' NEXTRA = ',i7 &
         ,/' NCOPY  = ',i7/)
   
   !     --- make sure we do not exceed memory limits in commons ---
   
   nttyp = ntypes*(ntypes+1)/2
   if (numbnd > MAX_BOND_TYPE .or. &
       numang > 5000 .or. &
       nptra  > 3000 .or. &
       nphb   > 480000 .or. &
       natyp  > MAX_ATOM_TYPE .or. &
       nttyp  > 400000                 ) then
         write(6,'(/,5x,a)') 'rdparm: a parameter array overflowed'
         write(6,'(/,5x,a)') '       (e.g. the table of dihedral params)'
         write(6,8119) numbnd,MAX_BOND_TYPE, &
                       numang,nptra, &
                       natyp,MAX_ATOM_TYPE, &
                       nphb,nttyp
 8119    format( &
          ' NUMBND = ',i7,'  max is',i5, &
         /' NUMANG = ',i7, '  max is 5000', &
         /' NPTRA  = ',i7, '  max is 3000', &
         /' NATYP  = ',i7, '  max is',i5,  &
         /' NPHB   = ',i7, '  max is 480000', &
         /' NTTYP  = ',i7, '  max is 400000'  &
          )
   
      call mexit(6, 1)
   end if
   
   if(nbona /= mbona .or. ntheta /= mtheta .or. nphia /= mphia) then
      write(6,*) 'Sander no longer allows constraints in prmtop'
      write(6,*) '...must have nbona=mbona, ntheta=mtheta, nphi=mphi'
      call mexit(6,1)
   end if

   ! Write implicit solvent radius and screening info to mdout
   if (( igb /= 0 .and. (ifcap == 0 .or. ifcap == 5)).or.hybridgb>0) then
      fmtin = afmt
      type = 'RADIUS_SET'
      call nxtsec(nf,  6,  1,fmtin,  type,  fmt,  iok)
      if (iok == 0) then         ! Allow failure to support pre-AMBER9 prmtop
         read(nf,fmt) type      ! Reuse type var to avoid declaring a new one
         write(6,'(A,A)') ' Implicit solvent radii are ',type
      end if
      if ( igb == 7 ) then
         write(6,'(A)') ' Replacing prmtop screening parameters with GBn (igb=7) values'
      end if
   end if

   !    --- Read the force field information from the prmtop if available
   fmtin = '(i,a)'
   type = 'FORCE_FIELD_TYPE'
   call nxtsec(nf,  6,  1,fmtin,  type,  fmt,  iok)
   if(iok == 0) then
     ! We found a force field description. Should be 1 or more lines
     ! of text which will be echoed back to the user. In some case
     ! we also take some decisions based on this.
     write (6,'(a)') '| Force field information read from topology file: '
     ! format should be (nlines, text) where we read the first lines nlines
     ! to know how many more lines to read - we ignore nlines on subsequent
     ! lines but they still need to be there.
     ! e.g.
     ! 2 PARM99
     ! 2 GAFF
     read(nf,fmt) nlines,line
     allocate (ffdesc(nlines), stat = ier)
     REQUIRE(ier==0)
     ffdesc(1) = line
     write(6,'(a,a)') '| ',ffdesc(1)
     do i = 2, nlines
       read(nf,fmt) iscratch,ffdesc(i)
       write(6,'(a,a)') '| ',ffdesc(i)
     end do

     !Test to see if any special force fields are in use.
     
     !1) Check for CHARMM
     charmm = .false.
     do i = 1, nlines
       if (index(ffdesc(i),'CHARMM') /= 0) then
         write(6,'("|")')
         write(6,'(a,a)') '| CHARMM force field in use. '
         charmm = .true.

         !AMBER 10 does not support CHARMM prmtop files quit
         !with an error message.
         call sander_bomb('rdparm','CHARMM prmtop file detected.', &
                          'AMBER 10 does not support CHARMM prmtop files.')

         !test to see if scee and scnb are set to 1.0
         if (scee /= 1.0d0) then
           write(6,'(a)') ' WARNING: CHARMM force field is in use but scee /= 1.0.'
           write(6,'(a,f4.2)') '          scee = ',scee
           write(6,'(a)') '          results will be INCORRECT!'
           REQUIRE(scee==1.0d0)
         end if
         if (scnb /= 1.0d0) then
           write(6,'(a)') ' WARNING: CHARMM force field is in use but scnb /= 1.0.'
           write(6,'(a,f4.2)') '          scnb = ',scnb
           write(6,'(a)') '          results will be INCORRECT!'
           REQUIRE(scnb==1.0d0)
         end if
       end if
     end do

     if(charmm)then
        type = 'CHARMM_NUM_IMPROPERS'
        call nxtsec(nf,  6,  1,fmtin,  type,  fmt,  iok)
        read(nf,fmt) nimphi
        write(0,*)"Got number of impropers: ",nimphi
        allocate(im(nimphi),jm(nimphi),km(nimphi),lm(nimphi),imp(nimphi), &
             stat=ier)
        REQUIRE(ier==0)
     endif


     !End Test
     deallocate (ffdesc, stat = ier)
     REQUIRE(ier==0)

   end if 

   return
end subroutine rdparm1 

!=======================================================================

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!+ [Enter a one-line description of subroutine rdparm2 here]
subroutine rdparm2(x,ix,ih,ipairs,nf)
   use parms
   use nblist, only: a,b,c
   use amoeba_mdin,only: iamoeba
   use pimd_vars, only: ipimd, dmdlm, itimass
#ifdef DSSP
   use dssp, only: ipepc, npepc
#endif
   implicit none
   integer nf
   _REAL_ x(*)
   integer ix(*),ipairs(*)
   character(len=4) ih(*)
   
#  include "md.h"
#  include "memory.h"
#  include "files.h"
#  include "box.h"
#  include "nmr.h"
#ifdef LES
#  include "les.h"
   integer iexcl,numex,k1,j1
#endif
   integer nttyp,ntype,i,iok
   integer i1,i2,i3,i4,j,jj,k,l,n,nn
   integer iptres,nspsol,natsm,idum,ip14
   integer l_ib,l_jb,l_bt
   integer l_it,l_jt,l_kt,l_tt
   integer l_id,l_jd,l_kd,l_ld,l_dt
   integer bp,ibp,jbp,btp
   _REAL_ dumd,oldbeta,duma,dumb,dumc
   character(len=4) dumchar
   integer dumint
   _REAL_ dumfloat
   _REAL_ massdiff   ! = mass[perturbed] - mass[original] for TI w.r.t. mass.
   integer irotat( natom )  ! dummy: disappears when leaving subroutine

   integer allocate_err
   integer ierr ! Allocation status.

   character(len=80) fmt
   character(len=80) fmtin,ifmt,afmt,rfmt,type

   ifmt = '(12I6)'
   afmt = '(20A4)'
   rfmt = '(5E16.8)'
   nttyp = ntypes*(ntypes+1)/2
   ntype = ntypes*ntypes
   
   !     ----- READ THE SYMBOLS AND THE CHARGES AND THE MASSES -----
   
   fmtin = afmt
   type = 'ATOM_NAME'
   call nxtsec(nf,  6,  0,fmtin,  type,  fmt,  iok)
   read(nf,fmt) (ih(m04+i-1),i = 1,natom)
   
   fmtin = rfmt
   type = 'CHARGE'
   call nxtsec(nf,  6,  0,fmtin,  type,  fmt,  iok)
   read(nf,fmt) (x(l15+i-1),i = 1,natom)
   
   fmtin = rfmt
   type = 'MASS'
   call nxtsec(nf,  6,  0,fmtin,  type,  fmt,  iok)
   read(nf,fmt) (x(lwinv+i-1),i = 1,natom)
   
   fmtin = ifmt
   type = 'ATOM_TYPE_INDEX'
   call nxtsec(nf,  6,  0,fmtin,  type,  fmt,  iok)
   read(nf,fmt) (ix(i04+i-1),i = 1,natom)
   
   fmtin = ifmt
   type = 'NUMBER_EXCLUDED_ATOMS'
   call nxtsec(nf,  6,  0,fmtin,  type,  fmt,  iok)
   read(nf,fmt) (ix(i+i08-1),i = 1,natom)
   
   fmtin = ifmt
   type = 'NONBONDED_PARM_INDEX'
   call nxtsec(nf,  6,  0,fmtin,  type,  fmt,  iok)
   read(nf,fmt) (ix(i+i06-1),i = 1,ntype)
   
   fmtin = afmt
   type = 'RESIDUE_LABEL'
   call nxtsec(nf,  6,  0,fmtin,  type,  fmt,  iok)
   read(nf,fmt) (ih(i+m02-1),i=1,nres)
   
   fmtin = ifmt
   type = 'RESIDUE_POINTER'
   call nxtsec(nf,  6,  0,fmtin,  type,  fmt,  iok)
   read(nf,fmt) (ix(i+i02-1),i=1,nres)
   ix(i02+nres) = natom+1
   
   !     ----- READ THE PARAMETERS -----
   
   fmtin = rfmt
   type = 'BOND_FORCE_CONSTANT'
   call nxtsec(nf,  6,  0,fmtin,  type,  fmt,  iok)
   read(nf,fmt) (rk(i),    i = 1,numbnd)
   
   fmtin = rfmt
   type = 'BOND_EQUIL_VALUE'
   call nxtsec(nf,  6,  0,fmtin,  type,  fmt,  iok)
   read(nf,fmt) (req(i),   i = 1,numbnd)
   
   fmtin = rfmt
   type = 'ANGLE_FORCE_CONSTANT'
   call nxtsec(nf,  6,  0,fmtin,  type,  fmt,  iok)
   read(nf,fmt) (tk(i),    i = 1,numang)
   
   fmtin = rfmt
   type = 'ANGLE_EQUIL_VALUE'
   call nxtsec(nf,  6,  0,fmtin,  type,  fmt,  iok)
   read(nf,fmt) (teq(i),   i = 1,numang)
   
   fmtin = rfmt
   type = 'DIHEDRAL_FORCE_CONSTANT'
   call nxtsec(nf,  6,  0,fmtin,  type,  fmt,  iok)
   read(nf,fmt) (pk(i),    i = 1,nptra)
   
   fmtin = rfmt
   type = 'DIHEDRAL_PERIODICITY'
   call nxtsec(nf,  6,  0,fmtin,  type,  fmt,  iok)
   read(nf,fmt) (pn(i),    i = 1,nptra)
   
   fmtin = rfmt
   type = 'DIHEDRAL_PHASE'
   call nxtsec(nf,  6,  0,fmtin,  type,  fmt,  iok)
   read(nf,fmt) (phase(i), i = 1,nptra)
   
   fmtin = rfmt
   type = 'SOLTY'
   call nxtsec(nf,  6,  0,fmtin,  type,  fmt,  iok)
   read(nf,fmt) (solty(i), i = 1,natyp)
   
   fmtin = rfmt
   type = 'LENNARD_JONES_ACOEF'
   call nxtsec(nf,  6,  0,fmtin,  type,  fmt,  iok)
   read(nf,fmt) (cn1(i),   i = 1,nttyp)
   
   fmtin = rfmt
   type = 'LENNARD_JONES_BCOEF'
   call nxtsec(nf,  6,  0,fmtin,  type,  fmt,  iok)
   read(nf,fmt) (cn2(i),   i = 1,nttyp)
   
   !     ----- READ THE BONDING INFORMATIONS -----
   
   fmtin = ifmt
   type = 'BONDS_INC_HYDROGEN'
   call nxtsec(nf,  6,  0,fmtin,  type,  fmt,  iok)
   read(nf,fmt) (ix(i+iibh-1),ix(i+ijbh-1),ix(i+iicbh-1), &
         i = 1,nbonh)
   
   fmtin = ifmt
   type = 'BONDS_WITHOUT_HYDROGEN'
   call nxtsec(nf,  6,  0,fmtin,  type,  fmt,  iok)
   read(nf,fmt)(ix(i+iiba-1),ix(i+ijba-1),ix(i+iicba-1),i = 1,nbona)
   
   fmtin = ifmt
   type = 'ANGLES_INC_HYDROGEN'
   call nxtsec(nf,  6,  0,fmtin,  type,  fmt,  iok)
   read(nf,fmt) (ix(i+i24-1),ix(i+i26-1),ix(i+i28-1),ix(i+i30-1), &
         i = 1,ntheth)
   
   fmtin = ifmt
   type = 'ANGLES_WITHOUT_HYDROGEN'
   call nxtsec(nf,  6,  0,fmtin,  type,  fmt,  iok)
   read(nf,fmt) (ix(i+i32-1),ix(i+i34-1),ix(i+i36-1),ix(i+i38-1), &
         i = 1,ntheta)
   
   fmtin = ifmt
   type = 'DIHEDRALS_INC_HYDROGEN'
   call nxtsec(nf,  6,  0,fmtin,  type,  fmt,  iok)
   read(nf,fmt) (ix(i+i40-1),ix(i+i42-1),ix(i+i44-1),ix(i+i46-1), &
         ix(i+i48-1),i = 1,nphih)
   
   fmtin = ifmt
   type = 'DIHEDRALS_WITHOUT_HYDROGEN'
   call nxtsec(nf,  6,  0,fmtin,  type,  fmt,  iok)
   read(nf,fmt) (ix(i+i50-1),ix(i+i52-1),ix(i+i54-1),ix(i+i56-1), &
         ix(i+i58-1),i = 1,nphia)
   
   fmtin = ifmt
   type = 'EXCLUDED_ATOMS_LIST'
   call nxtsec(nf,  6,  0,fmtin,  type,  fmt,  iok)
   read(nf,fmt) (ix(i+i10-1),i=1,nnb)
   
   !     ----- READ THE H-BOND PARAMETERS -----
   
   fmtin = rfmt
   type = 'HBOND_ACOEF'
   call nxtsec(nf,  6,  0,fmtin,  type,  fmt,  iok)
   read(nf,fmt) (asol(i),i=1,nphb)
   
#ifndef HAS_10_12
   do i=1,nphb
      if( asol(i) /= 0.d0 ) then
         write(6,*) 'Found a non-zero 10-12 coefficient, but source', &
               ' was not compiled with -DHAS_10_12.'
         write(6,*) 'If you are using a pre-1994 force field, you', &
               ' will need to re-compile with this flag.'
         call mexit(6,1)
      end if
   end do
#endif
   
   fmtin = rfmt
   type = 'HBOND_BCOEF'
   call nxtsec(nf,  6,  0,fmtin,  type,  fmt,  iok)
   read(nf,fmt) (bsol(i),i=1,nphb)
   
   fmtin = rfmt
   type = 'HBCUT'
   call nxtsec(nf,  6,  0,fmtin,  type,  fmt,  iok)
   read(nf,fmt) (hbcut(i),i=1,nphb)
   
   !     ----- READ ISYMBL,ITREE,JOIN AND IROTAT ARRAYS -----
   
   fmtin = afmt
   type = 'AMBER_ATOM_TYPE'
   call nxtsec(nf,  6,  0,fmtin,  type,  fmt,  iok)
   read(nf,fmt) (ih(i+m06-1),i=1,natom)
   
   fmtin = afmt
   type = 'TREE_CHAIN_CLASSIFICATION'
   call nxtsec(nf,  6,  0,fmtin,  type,  fmt,  iok)
   read(nf,fmt) (ih(i+m08-1),i=1,natom)
   
   fmtin = ifmt
   type = 'JOIN_ARRAY'
   call nxtsec(nf,  6,  0,fmtin,  type,  fmt,  iok)
   read(nf,fmt) (ix(i+i64-1),i=1,natom)
   
   fmtin = ifmt
   type = 'IROTAT'
   call nxtsec(nf,  6,  0,fmtin,  type,  fmt,  iok)
   read(nf,fmt) (irotat(i), i=1,natom)
   
   !     ----- READ THE BOUNDARY CONDITION STUFF -----
   
   nspm = 1
   ix(i70) = natom
   if (ifbox > 0) then
      
      fmtin = ifmt
      type = 'SOLVENT_POINTERS'
      call nxtsec(nf,  6,  0,fmtin,  type,  fmt,  iok)
      read(nf,fmt) iptres,nspm,nspsol
      
      fmtin = ifmt
      type = 'ATOMS_PER_MOLECULE'
      call nxtsec(nf,  6,  0,fmtin,  type,  fmt,  iok)
      read(nf,fmt) (ix(i+i70-1),i=1,nspm)
      
      fmtin = rfmt
      type = 'BOX_DIMENSIONS'
      call nxtsec(nf,  6,  0,fmtin,  type,  fmt,  iok)
      read(nf,fmt) oldbeta,duma,dumb,dumc
      
      !       ---(above values are still read, for backward compatibility, but
      !           ignored. Box info must come from the coord. file or from
      !           the &ewald namelist of the input file)
      
      if( igb /= 0  .or.  ntb == 0 )then
         box(1)=0.0d0
         box(2)=0.0d0
         box(3)=0.0d0
      else
         box(1)=a
         box(2)=b
         box(3)=c
      end if
      
   end if  ! (ifbox > 0)
   
   !     ----- LOAD THE CAP INFORMATION IF NEEDED -----
   
   if(ifcap == 1) then
      fmtin = '(I6)'
      type = 'CAP_INFO'
      call nxtsec(nf,  6,  0,fmtin,  type,  fmt,  iok)
      read(nf,fmt) natcap
      
      fmtin = '(4E16.8)'
      type = 'CAP_INFO2'
      call nxtsec(nf,  6,  0,fmtin,  type,  fmt,  iok)
      read(nf,fmt) cutcap,xcap,ycap,zcap
   end if
   
   if( igb /= 0 .and. ifcap == 0 .and. iok == -1 ) then
      write(0,*) 'GB/PB calculations now require a new-style prmtop file'
      write(6,*) 'GB/PB calculations now require a new-style prmtop file'
      call mexit(6,1)
   end if
   
   if (( igb /= 0 .and. (ifcap == 0 .or. ifcap == 5)).or.hybridgb>0) then
      fmtin = rfmt
      type = 'RADII'
      call nxtsec(nf,  6,  0,fmtin,  type,  fmt,  iok)
      read(nf,fmt) (x(l97+i-1),i=1,natom)
      type = 'SCREEN'
      call nxtsec(nf,  6,  0,fmtin,  type,  fmt,  iok)
      read(nf,fmt) (x(l96+i-1),i=1,natom)
   end if
   
   if (ipol > 0) then
      fmtin = rfmt
      type = 'POLARIZABILITY'
      call nxtsec(nf,  6,  0,fmtin,  type,  fmt,  iok)
      read(nf,fmt) (x(lpol+i-1),i=1,natom)
   end if

   !     ----- READ THE PERTURBED MASSES IF NEEDED  -----
   if (itimass > 0) then
      ! JVAN: dmdlm must be allocated here, not in pimd_init.
      allocate( dmdlm(1:natom),stat=ierr )
      REQUIRE( ierr == 0 )
      fmtin = rfmt
      type = 'TI_MASS'
      call nxtsec(nf,  6,  0,fmtin,  type,  fmt,  iok)
      read(nf,fmt) (dmdlm(i) ,i = 1,natom)
      do i=1,natom
         massdiff = (dmdlm(i) - x(lwinv+i-1))
         x(lwinv+i-1) = x(lwinv+i-1) + clambda * massdiff
         dmdlm(i) = massdiff/x(lwinv+i-1)
      end do
   end if
   
#ifdef DSSP
   !   ----- construct an array containing the atom numbers of the carbon atoms of all peptide
   !         groups
   allocate( ipepc(1:natom),stat=ierr )
   REQUIRE( ierr == 0 )
   k = 0
   do i=1,natom
      if( ih(m04+i-1) == 'C   ' ) then
         k = k+1
         ipepc(k) = i
      end if
   end do
   npepc = k
#endif
#ifdef LES

   if (nparm == 1.and.iamoeba.eq.0) then
      fmtin = ifmt
      type = 'LES_NTYP'
      call nxtsec(nf,  6,  0,fmtin,  type,  fmt,  iok)
      read (nf,fmt) nlesty
      lestmp=nlesty*nlesty
      
      ! check the array sizes to make sure we do not overflow
      
      ! LES types
      
      if (nlesty > maxlestyp) then
         write (6,*) 'Exceeded MAXLESTYP',nlesty
         stop
      end if
      
      ! LES atoms
      
      if (natom > maxles) then
         write (6,*) 'Exceeded MAXLES',natom
         stop
      end if

      fmtin = ifmt
      type = 'LES_TYPE'
      call nxtsec(nf,  6,  0,fmtin,  type,  fmt,  iok)
      read (nf,fmt) (lestyp(i),i=1,natom)
      fmtin = rfmt
      type = 'LES_FAC'
      call nxtsec(nf,  6,  0,fmtin,  type,  fmt,  iok)
      read (nf,fmt) (lesfac(i),i=1,lestmp)
      fmtin = ifmt
      type = 'LES_CNUM'
      call nxtsec(nf,  6,  0,fmtin,  type,  fmt,  iok)
      read (nf,fmt) (cnum(i), i=1,natom)
      fmtin = ifmt
      type = 'LES_ID'
      call nxtsec(nf,  6,  0,fmtin,  type,  fmt,  iok)
      read (nf,fmt) (subsp(i), i=1,natom)
      write (6,*) 'LES parameters were found'
      
      ! now create the list of atoms that have non-unitary scaling factors.
      ! this will be used in the Ewald calculation to correct for the
      ! lack of use of the intra-copy scaling factor in the charge grid.
      ! all of these pairs will need correction. The list will not change and
      ! is therefore only calculated once (here).
      
      nlesadj=0
      
      iexcl=0
      
      ! pairs are listed in two arrays, for i and j, rather than using
      ! a set of pointers like the nonbond and exclusion lists. This is 
      ! since many atoms will not have any correction partners (since 
      ! they are not in LES).

      if( ipimd.eq.0 ) then
      ! pimd are not going to use nb_adjust_les, so we do not need to generate
      ! les adjust list
      !
      do 6500 k=1,natom
         
         lestmp=nlesty*(lestyp(k)-1)
         !         write (6,*) 'atom1 : ',k,lestmp
         
         ! need to sum all f the number of exclusions even if non-LES atoms.
         ! see below.
         
         numex=ix(k+i08-1)
         
         do 6510 j=k+1,natom
            lfac=lesfac(lestmp+lestyp(j))
            
            ! check for non-zero scaling factor (meaning a correction will 
            ! be required)
            
            if (abs(lfac-1.0d0) > 0.01) then
               
               ! check to make sure these aren't excluded atoms (since then 
               ! no correction is wanted)
               
               !  FORMAT(12I6)  (NATEX(i), i=1,NEXT)
               !  the excluded atom list.  To get the excluded list for atom
               !  "i" you need to traverse the NUMEX list, adding up all
               !  the previous NUMEX values, since NUMEX(i) holds the number
               !  of excluded atoms for atom "i", not the index into the
               !  NATEX list.  Let IEXCL = SUM(NUMEX(j), j=1,i-1), then
               !  excluded atoms are NATEX(IEXCL)+1 to NATEX(IEXCL+NUMEX(i)).
               
               
               do k1=iexcl+1,iexcl+numex
                  
                  ! get exclusion
                  
                  j1=ix(k1+i10-1)
                  
                  ! check against atom j
                  
                  if (j1 == j) then
                     ! excluded, get next atom j
                     ! write (6,*) 'Exclusion list match'
                     goto 6510
                  end if
                  
                  ! check next entry
                  
               end do
               
               ! if we arrived here, the atom was not in the exclusion list
               ! so this pair will need correction
               ! (should add boundary checking for variables here)
               
               if (nlesadj == maxlesadj) then
                  write (6,*) 'EXCEEDED MAXLESADJ!'
                  stop
               end if
               
               nlesadj=nlesadj+1
               ileslst(nlesadj)=k
               jleslst(nlesadj)=j
            end if
            
            ! next j
            
         6510 continue
         
         ! increment the exclusion list pointer for atom i
         iexcl=iexcl+numex
         
         ! next i
         
      6500 continue
      
      end if !(ipimd == 0 )

      write (6,6520) nlesadj
      6520 format (1x,i7,' LES atom pairs require adjustment')
      
      ! end creation of LES adjustment list and reading LES info

   end if  ! (nparm==1)

#endif /* LES */
   
   !     ----- CALCULATE INVERSE, TOTAL MASSES -----
   !       -- save the masses for removal of mass weighted velocity,
   !          leaving the inverse masses in the legacy, Lwinv area
   
   tmass = 0.0d0
   !     -- index over molecules
   j = l75-1
   jj = i70-1
   !     -- index over mass->invmass
   k = lwinv-1
   !     -- index over saved mass
   l = lmass-1
   do n = 1,nspm
      j = j + 1
      jj = jj + 1
      x(j) = 0.0d0
      natsm = ix(jj)
      do nn = 1,natsm
         k = k+1
         l = l+1
         
         ! -- sum molecule
         x(j) = x(j) + x(k)
         
         ! -- save mass in "new" Lmass area
         x(l) = x(k)
         
         ! -- make inverse in "old" Lwinv area
         if( x(k) /= 0.d0 ) x(k) = 1.0d0 / x(k)
      end do
      tmass = tmass + x(j)
   end do
   tmassinv = 1.0d0 / tmass

 
   !     ----- SCALE THE CHARGES IF DIELC.GT.1.0E0 -----
   
   if (dielc /= 1.0e0 .and. igb == 0) then
      dumd = sqrt(dielc)
      do i = 1,natom
         x(i+l15-1) = x(i+l15-1)/dumd
      end do
   end if
   
   !     ----- INVERT THE HBCUT ARRAY -----
   
   do i = 1,nphb
      if(hbcut(i) <= 0.001e0) hbcut(i) = 1.0d-10
      hbcut(i) = 1.0e0/hbcut(i)
   end do
   
   !     ----- duplicate dihedral pointers for vector ephi -----
   
   call dihdup(nphih,ix(i40),ix(i42),ix(i44),ix(i46),ix(i48),pn)
   call dihdup(nphia,ix(i50),ix(i52),ix(i54),ix(i56),ix(i58),pn)
   
   !     --- pre-calculate some parameters for vector ephi ---
   
   call dihpar(nptra,pk,pn,phase,gamc,gams,ipn,fmn)
   
   if (charmm) then
     !    ---read in   1-4 parameters
   
!   allocate(cn114(nttyp),cn214(nttyp),stat=allocate_err)
!   if(allocate_err /= 0) &
!        call sander_bomb('rdparm2<rdparm.f>', &
!        'cannot allocate urey-bradley arrays',' BOMBING')
     fmtin = rfmt
     type = 'LENNARD_JONES_14_ACOEF'
     call nxtsec(nf,  6,  0,fmtin,  type,  fmt,  iok)
     read(nf,fmt) (cn114(i),   i = 1,nttyp)
   
     fmtin = rfmt
     type = 'LENNARD_JONES_14_BCOEF'
     call nxtsec(nf,  6,  0,fmtin,  type,  fmt,  iok)
     read(nf,fmt) (cn214(i),   i = 1,nttyp)

     !   
     !    --- read in   urey-bradley parameters:
 
  !   allocate(rub(numang),rkub(numang),stat=allocate_err)
  !   if(allocate_err /= 0) &
  !        call sander_bomb('rdparm2<rdparm.f>', &
  !        'cannot allocate urey-bradley arrays',' BOMBING')
     fmtin =  rfmt
     type = 'UREY_BRADLEY_EQUIL_VALUE'
     call nxtsec(nf,  6,  0,fmtin,  type,  fmt,  iok)
     read(nf,fmt) (rub(i),   i = 1,numang)
     fmtin = rfmt
     type = 'UREY_BRADLEY_FORCE_CONSTANT'
     call nxtsec(nf,  6,  0,fmtin,  type,  fmt,  iok)
     read(nf,fmt) (rkub(i),    i = 1,numang)

     fmtin =  '(i8)'
     type = 'CHARMM_NUM_IMPR_TYPES'
     call nxtsec(nf,  6,  0,fmtin,  type,  fmt,  iok)
     read(nf,fmt) nimprtyp
     write(0,*)"Got number of improper types: ",nimprtyp
     fmtin = rfmt
     type = 'CHARMM_IMPROPER_FORCE_CONSTANT'
     call nxtsec(nf,  6,  0,fmtin,  type,  fmt,  iok)
     read(nf,fmt) (pk_impr(i),    i = 1,nimprtyp)
     write(0,fmt)(pk_impr(i),    i = 1,nimprtyp)
     type = 'CHARMM_IMPROPER_PHASE'
     call nxtsec(nf,  6,  0,fmtin,  type,  fmt,  iok)
     read(nf,fmt) (phase_impr(i),    i = 1,nimprtyp)
     write(0,fmt)(phase_impr(i),    i = 1,nimprtyp)


   end if
   
   return
end subroutine rdparm2 


!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!+ [Enter a one-line description of subroutine istuff here]
subroutine istuff(i,j,iarray,k)
   
   ! routine to correctly load a strange shaped 2 dim matrix
   
   dimension iarray(15,*)
   iarray(i,j) = k
   return
end subroutine istuff 
