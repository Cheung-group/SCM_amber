#include "dprec.h"

module pupildata
   logical pupactive
   integer puperror,iPup,jPup,bs1,bs2,iresPup,l_puptmp
   integer pupLevelData                 ! Data to pass throught PUPIL Interface
   integer pupStep
   integer pupQZchange                  ! is != 0 if quantum zone has changed
   integer pupqatoms                    ! total number of PUPIL quantum atoms
   integer pupnumnb14                   ! total number of initial nonboded 14 pairs
   integer pupnbonh                     ! total number of initial boded H  pairs
   integer pupnbona                     ! total number of initial boded    pairs
   integer pupntheth                    ! total number of initial angled H  tern.
   integer pupntheta                    ! total number of initial angled    tern.
   integer pupnphih                     ! total number of initial dihed. H  quat.
   integer pupnphia                     ! total number of initial dihed.    quat.
   integer, allocatable :: pupnb14  (:) ! initial nonbonded 14 pair list  
   integer, allocatable :: pupbonh  (:) ! initial boded H  pairs
   integer, allocatable :: pupbona  (:) ! initial boded    pairs
   integer, allocatable :: puptheth (:) ! initial angled H  tern.
   integer, allocatable :: puptheta (:) ! initial angled    tern.
   integer, allocatable :: pupphih  (:) ! initial dihed. H  quat.
   integer, allocatable :: pupphia  (:) ! initial dihed.    quat.
   integer, allocatable :: pupmask  (:) ! mask: 0 -> classic 1-> quantum
   integer, allocatable :: pupqlist (:) ! atom number for quantum list
   integer, allocatable :: pupatm   (:)
   integer, allocatable :: pupres   (:) ! residue pointers
   integer, allocatable :: ixStack  (:) ! Captioning the ix pointer
   _REAL_,  allocatable :: pupchg   (:) ! Initial charges over each atoms in the system
   _REAL_,  allocatable :: qcdata   (:)
   _REAL_,  allocatable :: qfpup    (:)
   _REAL_,  allocatable :: qcell    (:)
   _REAL_,  allocatable :: realStack(:) ! Captioning the x pointer
   character(len=4), dimension(:), allocatable :: ihStack  ! Captioning the ih pointer
   character(len=10), dimension(:),allocatable :: keyres   ! to keep the MM key  RESIDUE_NAME
   character(len=20),dimension(:), allocatable :: keyMM    ! to keep the MM key particles RESIDUE_NAME.ATOM_NAME
   character*20  strAux
   _REAL_   qmEnergy
end module pupildata

! ***********************************************************************
! ***********************************************************************

subroutine get_atomic_number(name,atomic_number)
!
!
!     This subroutine assigns atomic number based upon the first and second
!     letter of the atom symbol. 

!     name ==>   character containing the atomic name
!     iqm_atomic_numbers ==>   integer array of atomic numbers assigned to atoms
!
   implicit none

   character(len=4) :: name
   integer :: atomic_number
  
   integer :: i, j, k, coincidence
   character*2, dimension(111):: elemnt=(/                &
    'H ','HE',                                            &
    'LI','BE','B ','C ','N ','O ','F ','NE',              &
    'NA','MG','AL','SI','P ','S ','CL','AR',              &
    'K ','CA',                                            &
    'SC','TI','V ','CR','MN','FE','CO','NI','CU','ZN',    &
              'GA','GE','AS','SE','BR','KR',              &
    'RB','SR',                                            & 
    'Y ','ZR','NB','MO','TC','RU','RH','PD','AG','CD',    &
              'IN','SN','SB','TE','I ','XE',              &
    'CS','BA',                                            &
    'LA','CE','PR','ND','PM','SM','EU',                   &
    'GD','TB','DY','HO','ER','TM','YB',                   &
    'LU','HF','TA','W ','RE','OS','IR','PT','AU','HG',    &
              'TL','PB','BI','PO','AT','RN',              &
    'FR','RA',                                            &
    'AC','TH','PA','U ','NP','PU','AM',                   &
    'CM','BK','CF','ES','FM','MD','NO',                   &          
    'LR','RF','DB','SG','BH','HS','MT','DS','RG'/)
!
!      compare values in name to those in elemnt and assign atomic
!      numbers accordingly
!
   coincidence = 0
   do j=1,111
     if(name(1:1) .eq. elemnt(j)(1:1)) then
        if(elemnt(j)(2:2) .eq. " ") then
           if(name(2:2) .eq. " ") then
             atomic_number = j
             coincidence = coincidence + 1
             !write(6,*) name,' COINCIDENCE .... with ',elemnt(j)
           else
             do k=j+1,111
               if((name(2:2) .eq. elemnt(k)(1:1)) .and.   &
                  (name(2:2) .eq. elemnt(k)(2:2)))then
                 atomic_number = j
                 coincidence = coincidence + 1
                 !write(6,*) name,' COINCIDENCE .... with ',elemnt(k)
               endif
             enddo
           endif
        else
           if(name(2:2) .eq. elemnt(j)(2:2)) then
              atomic_number = j
              coincidence = coincidence + 1
              !write(6,*) name,' COINCIDENCE .... with ',elemnt(j)
           endif
        endif
     end if
   enddo

   if(coincidence .gt. 1) then
     write(6,*) 'PUPIL: Unable to correctly identify element ', name
     call mexit(6,1)
!      else
!        write(6,*) 'FOUND: ',elemnt(atomic_number),' related to ',name
   endif

   return

end subroutine get_atomic_number

! ***********************************************************************
! ***********************************************************************

subroutine deleting_qm_atoms()
   use pupildata, x=>realStack , ix=>ixStack, ih=>ihStack
   use parms, only: req, fmn
   
#include "memory.h"
#include "extra_pts.h"
  
  integer replicates
  integer i
  replicates  = 1

!  Initializing the list of structures for a new QM zone
  call init_extra_pts( &
         ix(iibh),ix(ijbh),ix(iicbh), &
         ix(iiba),ix(ijba),ix(iicba), &
         ix(i24),ix(i26),ix(i28),ix(i30), &
         ix(i32),ix(i34),ix(i36),ix(i38), &
         ix(i40),ix(i42),ix(i44),ix(i46),ix(i48), &
         ix(i50),ix(i52),ix(i54),ix(i56),ix(i58), &
         ih(m06),ix,x,ix(i08),ix(i10),fmn, &
         nspm,ix(i70),x(l75),tmass,tmassinv,x(lmass),x(lwinv),req)

  if(nbonh.gt.0) call setbon(nbonh,ix(iibh),ix(ijbh),ix(iicbh), &
    ix(ibellygp), replicates) ! remove bonds between QM atoms from list

  if(nbona.gt.0) call setbon(nbona,ix(iiba),ix(ijba),ix(iicba), &
    ix(ibellygp),replicates) ! remove bonds between QM atoms from list

  if(ntheth.gt.0) call setang(ntheth,ix(i24),ix(i26),ix(i28),ix(i30), &
    ix(ibellygp)) ! remove angles between QM atoms from list

  if(ntheta.gt.0) call setang(ntheta,ix(i32),ix(i34),ix(i36),ix(i38),&
    ix(ibellygp)) ! remove angles between QM atoms from list

  if(nphih.gt.0) call setdih(nphih,ix(i40),ix(i42),ix(i44),ix(i46), &
    ix(i48), ix(ibellygp)) ! remove dihedrals between QM atoms from list

  if(nphia.gt.0) call setdih(nphia,ix(i50),ix(i52),ix(i54),ix(i56), &
    ix(i58), ix(ibellygp)) ! remove dihedrals between QM atoms from list

  !if(pupQZChange .ne. 0) then
  !  write(6,*) 'NUMBER OF H BONDS            = ',nbonh
  !  write(6,*) 'NUMBER OF   BONDS            = ',nbona
  !  write(6,*) 'NUMBER OF H ANGLES           = ',ntheth
  !  write(6,*) 'NUMBER OF   ANGLES           = ',ntheta
  !  write(6,*) 'NUMBER OF H DIHEDRALS        = ',nphih
  !  write(6,*) 'NUMBER OF   DIHEDRALS        = ',nphia
  !  write(6,*) 'NUMBER OF NONBONDED 14 PAIRS = ',numnb14
  !  do i=1,numnb14
  !     write(6,*) ix(inb_14+(i-1)*2),ix(inb_14+(i-1)*2+1)
  !  enddo
  !endif
  
!  write(6,*) 'nbonh',nbonh,'iibh',iibh,'nbona',nbona,'ntheth',ntheth
!  write(6,*) 'ntheta',ntheta,'nphih',nphih,'nphia',nphia
!  write(6,*) 'i26',i26,'i56',i56
!  write(6,*) 'IGB',igb

  return
end subroutine deleting_qm_atoms

! ***********************************************************************
! ***********************************************************************

subroutine write_energies(str,e)
  _REAL_ e(25)
  character*4  str
  integer i,k
  do i=1,5
    write(6,"(a4,2x,'ENERG',5(2x,i2,2x,d15.8))")   &
            str,( (i-1)*5+k,e((i-1)*5+k), k=1,5 )
  enddo
  write(6,*) " ----------------------------------"
end subroutine write_energies


!***********************************************************************

