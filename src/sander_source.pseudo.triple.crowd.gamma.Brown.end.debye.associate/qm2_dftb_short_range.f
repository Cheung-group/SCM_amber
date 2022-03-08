! <compile=optimized>
!-*- mode: f90; coding: iso-8859-15; -*-

#include "dprec.h"


!==============================================================================
! evaluate short range expression i.e. sumR (gamma - 1/R)
!
! INPUT:
! _REAL_    rh(3)       vector between rmu and rnu
! _REAL_    umu         hubbard parameter of orbital mu
! _REAL_    unu         hubbard parameter of orbital nu
! _REAL_    basis(3,3)  basis of cell(unusual!!!:one line is a basis vector)
! _REAL_    tol         convergence tolerance (contribution of last shell)
!
! OUTPUT:
! _REAL_    value       value of the short range sum
!==============================================================================

subroutine shortrange(rh,umu,unu,tol,value)

   use nblist, only: ucell

   implicit none
   external gam12

   _REAL_ , intent(in ) :: rh(3)
   _REAL_ , intent(in ) :: umu
   _REAL_ , intent(in ) :: unu
   _REAL_ , intent(in ) :: tol
   _REAL_ , intent(out) :: value

!! Locals
   integer :: i,j,k,nreal,nmax,nmin
   _REAL_  :: rvec(3),R(3),lastshell,tmp,norm
   _REAL_  :: gval
   _REAL_  :: gamma_tmp

!! --

   gamma_tmp = 0.0d0
   nmax = 50   ! Sum at most 50 shells. If not converged, stop.
   nmin = 3    ! Sum at least 3 shells
   nreal = 0   ! Shell number
   lastshell = tol+1e-8 ! Contribution from the last shell calculated
   ! /* sum over R until tolerance is reached */
   do while ( (nreal .le. nmax) .and. &
              ((abs(lastshell) .gt. tol) .or. (nreal .le. nmin))  )
      lastshell = 0.0d0
      do i = -nreal,nreal
         do j = -nreal,nreal
            do k = -nreal,nreal
               ! /*only R belonging to outer shells are new ones */
               if((nreal == abs(i)) .or. (nreal == abs(j))  .or. &
                     (nreal == abs(k)) ) then

                  R(1)=i*ucell(1,1)+j*ucell(2,1)+k*ucell(3,1)
                  R(2)=i*ucell(1,2)+j*ucell(2,2)+k*ucell(3,2)
                  R(3)=i*ucell(1,3)+j*ucell(2,3)+k*ucell(3,3)

                  rvec(1) = rh(1) - R(1)
                  rvec(2) = rh(2) - R(2)
                  rvec(3) = rh(3) - R(3)

                  norm   = sqrt(rvec(1)**2+rvec(2)**2+rvec(3)**2)

                  ! get value for Gamma
                  call GAM12(norm,umu,unu,gval)

                  ! subtract long range part 1/R and multiply by Z(nu)

                  if (norm < 1.0e-6) then
                     tmp =  gval
                  else
                     tmp =  ( gval  - 1/norm )
                  endif

                  gamma_tmp = gamma_tmp + tmp
                  lastshell = lastshell + tmp

                  !write(6,'(4(I5),3(E20.10))') nreal, i,j,k,tmp,gamma_tmp,lastshell

               end if
            end do
         end do
      end do
      nreal = nreal + 1
   end do

   if(abs(lastshell) > tol) then
      write(6,'("lastshell = ",e20.10,10X,"tol = ",e20.10)') lastshell, tol
      call sander_bomb("shortrange <qm2_dftb_short_range.f>",&
         "Tolerance in subroutine shortrange not reached.",&
         "Exiting.")
   end if
   value = gamma_tmp

end subroutine SHORTRANGE




!=============================================================================
! evaluate derivative of short range expression: sumR (d gamma/dR - (-1/R^2))
!
! INPUT:
! _REAL_    rh(3)       vector between rmu and rnu
! _REAL_    umu         hubbard parameter of orbital mu
! _REAL_    unu         hubbard parameter of orbital nu
! _REAL_    basis(3,3)  basis of cell(unusual!!!:one line is a basis vector)
! _REAL_    tol         convergence tolerance (contribution of last shell)
!
! OUTPUT:
! _REAL_    deriv(3)    derivative of the short range sum
!==============================================================================

subroutine shortrange1(rh,umu,unu,tol,deriv)

   use nblist, only: ucell

   implicit none
   external gam121

   _REAL_ :: rh(3),umu,unu, tol, deriv(3)
   integer :: i,j,k,nreal,nmax,nmin
   _REAL_ :: rvec(3),R(3),lastshell,tmp,norm
   _REAL_ :: gdrv

   deriv(1) = 0.0
   deriv(2) = 0.0
   deriv(3) = 0.0
   nmax = 100
   nmin = 3
   nreal = 0

   lastshell = tol+1.0e-8
   ! /* sum over R until tolerance is reached */
   do while ((nreal .le. nmax) .and. ( (abs(lastshell) .gt. tol) &
         .or. (nreal .le. nmin)) )
      lastshell = 0.0
      do i = -nreal,nreal
         do j = -nreal,nreal
            do k = -nreal,nreal
               ! /*only R belonging to outer shells are new ones */
               if((nreal == abs(i)) .or. (nreal == abs(j))  .or. &
                     (nreal == abs(k)) ) then

                  R(1)=i*ucell(1,1)+j*ucell(2,1)+k*ucell(3,1)
                  R(2)=i*ucell(1,2)+j*ucell(2,2)+k*ucell(3,2)
                  R(3)=i*ucell(1,3)+j*ucell(2,3)+k*ucell(3,3)


                  rvec(1) = rh(1) - R(1)
                  rvec(2) = rh(2) - R(2)
                  rvec(3) = rh(3) - R(3)

                  norm   = sqrt(rvec(1)**2+rvec(2)**2+rvec(3)**2)

                  ! get derivative of gamma

                  call GAM121(norm,umu,unu,gdrv)

                  ! subtract long range -1/R^2/R (see def. of GAM121)
                  if (norm < 1.0e-6) then
                     tmp = gdrv
                  else
                     tmp = gdrv - (-1.0/(norm**3))
                  endif

                  deriv(1) = deriv(1) + tmp*rvec(1)
                  deriv(2) = deriv(2) + tmp*rvec(2)
                  deriv(3) = deriv(3) + tmp*rvec(3)

                  lastshell = lastshell + tmp
               end if
            end do
         end do
      end do
      nreal = nreal + 1
   end do

   if(abs(lastshell) > tol) then
      write(6,'("lastshell = ",e20.10,10X,"tol = ",e20.10)') lastshell, tol
      call sander_bomb("shortrange1 <qm2_dftb_short_range.f>",&
         "Tolerance in subroutine shortrange1 not reached.",&
         "Exiting.")
   end if

end subroutine SHORTRANGE1
