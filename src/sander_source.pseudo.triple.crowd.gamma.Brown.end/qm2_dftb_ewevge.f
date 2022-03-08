! <compile=optimized>

!***********************************************************************
! Parameters:
!
! NA      (I) :  Dimension of A
! NB      (I) :  Dimension of B
! N       (I) :  Dimension of Problem
! A       (I) :  Matrix A
!         (O) :  Eigenvector matrix
! B       (I) :  Matrix B
! EW      (O) :  Eigenvalues
! H       (-) :  Dummy vector
! AUX     (-) :  Auxiliary vector
! IEV     (I) :  0: No eigenvectors
! IER     (O) :  Dummy variable
!
! **********************************************************************
#include "dprec.h"
SUBROUTINE EWEVGE (NA,NB,N,A,B,EW,IEV,IORD,IER)

   use qm2_dftb_module, only: scr_space
   
   implicit none 

   !INTEGER :: NA, NB, N, IEV,IORD,IER
   !_REAL_  :: A(NA,N),B(NB,N),EW(N),H(N),AUX(3*NA)

!! Passed in
   integer, intent(in)    :: NA, NB, N, IEV,IORD
   integer, intent(out)   :: IER
   _REAL_ , intent(inout) :: A(NA,N)
   _REAL_ , intent(inout) :: B(NB,N)
   _REAL_ , intent(out)   :: EW(N)     ! output: eigenvalues

!! In DFTB,
!! IEV = 1

   CALL DSYGV(IEV,'V','L',N,A,NA,B,NB,EW,scr_space,3*NA,IER)
   
   RETURN
end SUBROUTINE EWEVGE

    ! ** To use DSYGVD (divide-and-conquer) **
    !DIMENSION  A(NA,N),B(NB,N),EW(N),H(N),AUX(3*NA)
    !integer :: IWORK(5*N+3)
    !integer :: LWORK, LIWORK

    !LWORK = 1 + 6*N + 2*N*N
    !CALL DSYGVD(IEV,'V','L',N,A,NA,B,NB,EW,AUX,LWORK,IWORK,5*N+3,IER)
