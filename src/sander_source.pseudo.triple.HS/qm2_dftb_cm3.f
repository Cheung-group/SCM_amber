! <compile=optimized> !  -*- mode: f90; coding: iso-8859-15; -*-
#include "dprec.h"

subroutine qm2_dftb_cm3(scf_mchg)

!! Calculates the CM3 charges, according to 
!! Kalinowsli et al, JPC-A, 108, 2545 (2004).
!!
!!
!! After returning, the cm3%qcm3(j) array will contain
!! CM3 charges.
!!
!! -------------------------------------------------
!! Implementation in Amber/SCC-DFTB by
!!
!! GUSTAVO SEABRA
!! Quantum Theory Project, University of Florida
!! e-mail: seabra@qtp.ufl.edu
!! url:    http://www.qtp.ufl.edu/~seabra
!! -------------------------------------------------

   use qmmm_module, only: qmmm_struct, qmmm_nml, qmmm_mpi
   use qm2_dftb_module, only: cm3, mol, ks_struct, izp_str, mcharge, NDIM, dacc

   implicit none

!! Passed in
   _REAL_, intent(in) :: scf_mchg(qmmm_struct%nquant_nlink) !Mulliken charge per atom.

!! Pointers
   integer, pointer :: ind(:)
   integer, pointer :: izp(:)
   _REAL_ , pointer :: qzero(:)     ! (NNDIM)     Atomic electron population of the neutral atom
   _REAL_ , pointer :: qmat(:)      ! (NNDIM)     Atomic electron population
   _REAL_ , pointer :: overl(:,:)   ! (MDIM,MDIM) Overlap Matrix
   _REAL_ , pointer :: density(:,:) ! (NDIM,NDIM) Density Matrix

!! Locals
   integer :: i,j, izpi, izpj, lambda, omega, adim
   integer :: lambda_begin, lambda_end, omega_begin, omega_end, l
   _REAL_  :: PS1, PS2, temp_b, temp_c, temp_d,temp_b2

!! For Matrix Printing   
   integer :: tot_n
   integer :: n_blocks_print, n_extra_columns,begin
   _REAL_  :: tr

!--

!! Set pointers
   ind     => ks_struct%ind
   overl   => ks_struct%overl
   density => ks_struct%density
   izp     => izp_str%izp
   qmat    => mol%qmat
   qzero   => mcharge%qzero
!!

   cm3%b = 0.0d0
   cm3%t = 0.0d0

!! Find the number of OCCUPIED MO's (adim)
!! MOVED HERE BY ROSS WALKER. It used to be stored in the module but this was
!! causing strange performance behaviour the mulliken calculation. So for the
!! time being we just recalculate it here and pass it to qm2_dftb_dens_matrix.
   do  i = 1,NDIM
      if(ks_struct%occ(i) < dacc) exit
   end do

   ! The number of occupied orbitals will get stored into adim
   adim = i-1

   ! Calculates the density matrix (DFTB doesn't by default)
   call qm2_dftb_dens_matrix(adim, ndim)

! write(6,*) "DENSITY MATRIX"
!write(6,*) density
!write(6,*)
!write(6,*) "OVERLAP MATRIX"
!write(6,*) overl
!! -------------------------------------------------------------
!! ---  Build Mayer bond order matrix, B: (the CM3 B-Matrix) ---
!! -------------------------------------------------------------
   do i = 1, qmmm_struct%nquant_nlink      ! (k)
      lambda_begin = ind(i)+1
      lambda_end   = ind(i+1)
!write(6,*) "i = ", i,"     lambda begin:",lambda_begin,"       end:",lambda_end
      do j = 1, qmmm_struct%nquant_nlink
         omega_begin = ind(j)+1
         omega_end   = ind(j+1)
         temp_b = 0.0d0
!write(6,*) "j = ", j,"      omega begin:", omega_begin,"       end:", omega_end
         do lambda = lambda_begin, lambda_end
            do omega = omega_begin, omega_end
               PS1 = 0.0d0
               PS2 = 0.0d0
!write(6,*) "   l begin: ",1, "   l end: ", NDIM
               do l = 1, NDIM
                  PS1 = PS1 + Density(omega ,l) * Overl(lambda,l)
                  PS2 = PS2 + Density(lambda,l) * Overl(omega ,l)
               end do
               temp_b = temp_b +  PS1 * PS2
!write(6,*) "TEMP_B = ", temp_b
            end do ! omega
         end do ! lambda
         cm3%b(i,j) = temp_b
      end do ! j
   end do ! i

!!! Print Mayer bond order matrix
!
!      n_blocks_print  = int(qmmm_struct%nquant_nlink / 5)
!      n_extra_columns = qmmm_struct%nquant_nlink - (5 * n_blocks_print)
!
!      ! Trace of Mayer BO matrix. Significance??
!      tr = 0.0d0
!      do i=1,qmmm_struct%nquant_nlink
!         tr = tr + cm3%b(i,i)
!         write(6,*)cm3%b(i,i),tr
!      end do
!      write(6,*) "TRACE of Mayer BO matrix= ", tr
!
!      begin = 1
!      do i = 1, n_blocks_print
!         write(6,*)
!         write(6,'(7X,5(3X,I3,"(",A2,")"))') (l,mol%atyp(izp(l)),l=begin,begin+4)
!         do j = 1, qmmm_struct%nquant_nlink
!            write(6,'(I3,"(",A2,")",5(F10.5))')j,mol%atyp(izp(j)),(cm3%b(l,j),l=begin,begin+4)
!         end do
!         begin = begin + 5
!      end do
!
!      write(6,*)
!      write(6,'(7X,5(3X,I3,"(",A2,")"))') (l,mol%atyp(izp(l)),l=begin,begin+n_extra_columns-1)
!      do j = 1, qmmm_struct%nquant_nlink
!         write(6,'(I3,"(",A2,")",5(F10.5))')j,mol%atyp(izp(j)),(cm3%b(l,j),l=begin,begin+n_extra_columns-1)
!      end do
!



!! -------------------------------------
!! ----  Calculates the CM3 T-Matrix ---
!! ------------------------------------- 
   do i = 1, qmmm_struct%nquant_nlink
      izpi = izp(i)
      do j = 1, qmmm_struct%nquant_nlink
         izpj = izp(j)
         temp_b  = cm3%b(i,j)
         temp_c  = cm3%c(izpi,izpj)
         temp_d  = cm3%d(izpi,izpj)
         temp_b2 = temp_b * temp_b
         cm3%t(i,j) =  temp_d * temp_b + temp_c * temp_b2 
      end do
   end do


!      write(6,*) " CM3: T-MATRIX"
!
!      begin = 1
!      do i = 1, n_blocks_print
!         write(6,*)
!         write(6,'(7X,5(3X,I3,"(",A2,")"))') (l,mol%atyp(izp(l)),l=begin,begin+4)
!         do j = 1, qmmm_struct%nquant_nlink
!            write(6,'(I3,"(",A2,")",5(F10.5))')j,mol%atyp(izp(j)),(cm3%t(l,j),l=begin,begin+4)
!         end do
!         begin = begin + 5
!      end do
!
!      write(6,*)
!      write(6,'(7X,5(3X,I3,"(",A2,")"))') (l,mol%atyp(izp(l)),l=begin,begin+n_extra_columns-1)
!      do j = 1, qmmm_struct%nquant_nlink
!         write(6,'(I3,"(",A2,")",5(F10.5))')j,mol%atyp(izp(j)),(cm3%t(l,j),l=begin,begin+n_extra_columns-1)
!      end do

!! ----------------------------------
!! ---  Calculate the CM3 charges ---
!! ----------------------------------
   do i = 1, qmmm_struct%nquant_nlink
      cm3%qcm3(i) = scf_mchg(i)
      do j = 1, qmmm_struct%nquant_nlink
         if (i /= j) then
            cm3%qcm3(i) = cm3%qcm3(i) + cm3%t(i,j)
         end if
      end do
   end do

!write(6,*) " Here: Done calculating CM3 charges..."

   ! DEBUG: WRITE OUTPUT
   if (qmmm_nml%verbosity > 3 .and. qmmm_mpi%commqmmm_master) then
      write(6,*)
      write(6,'(" AT  (TYP)  Symb      Qzero        Qmat     QMullik        QCM3    new Qmat")')
      !         " 12  (  1)   H      1.00000     0.91905     0.08095     0.11405     0.88595
      do j = 1, qmmm_struct%nquant_nlink
         write(6,'(I3,2X,"(",I3,")",3X,A2,5(2X,F10.5))') &
               j, izp(j),mol%atyp(izp(j)),qzero( izp(j) ), qmat(j), scf_mchg(j),cm3%qcm3(j),(qzero( izp(j) ) - cm3%qcm3(j))
      end do
   end if


   ! Update scf_mchg to contain the CM3 charges
   ! do j = 1, qmmm_struct%nquant_nlink
   !    scf_mchg(j) = cm3%qcm3(j)
   ! end do

!   return

!end subroutine qm2_dftb_cm3



end subroutine qm2_dftb_cm3

!subroutine print_matrix(matrix, m, n)
!
!   _REAL_ , intent(in) :: matrix
!   integer, intent(in) :: m
!   integer, intent(in) :: n
!
!
!      n_blocks_print  = int(qmmm_struct%nquant_nlink / 5)
!      n_extra_columns = qmmm_struct%nquant_nlink - (5 * n_blocks_print)
!
!      ! Trace of Mayer BO matrix. Significance??
!      tr = 0.0d0
!      do i=1,qmmm_struct%nquant_nlink
!         tr = tr + cm3%b(i,i)
!         write(6,*)cm3%b(i,i),tr
!      end do
!      write(6,*) "TRACE of Mayer BO matrix= ", tr
!
!      begin = 1
!      do i = 1, n_blocks_print
!         write(6,*)
!         write(6,'(7X,5(3X,I3,"(",A2,")"))') (l,mol%atyp(izp(l)),l=begin,begin+4)
!         do j = 1, qmmm_struct%nquant_nlink
!            write(6,'(I3,"(",A2,")",5(F10.5))')j,mol%atyp(izp(j)),(cm3%b(l,j),l=begin,begin+4)
!         end do
!         begin = begin + 5
!      end do
!
!      write(6,*)
!      write(6,'(7X,5(3X,I3,"(",A2,")"))') (l,mol%atyp(izp(l)),l=begin,begin+n_extra_columns-1)
!      do j = 1, qmmm_struct%nquant_nlink
!         write(6,'(I3,"(",A2,")",5(F10.5))')j,mol%atyp(izp(j)),(cm3%b(l,j),l=begin,begin+n_extra_columns-1)
!      end do
!
!
!
!end subroutine print_matrix
