!-----------  begin ew_mpole.h --------------------------------------------
!     indmeth......<&ewald variable> order of extrapolation in 1st estimate
!                   of the iterative process. DEFAULT = 3
!                   Doc says [0,1,2] for 1st,2nd,3rd order, mfc does not
!                       know yet what 3 is....





!   ***********************************************************
!   BC_MULTPOLE needs to be set to the size of the common block:
!   BC_INDDIPR
!   BC_INDDIPI

#define BC_MULTPOLE 13
integer &
      ifirst,  imiddle, ithird,  lfixdip,  linddip, &
      ldipole, lquad,   lfield,  ltorque,  leold1, &
      leold2,  leold3,  ldipvel
common/multpole/ &
      ifirst,  imiddle, ithird,  lfixdip,  linddip, &
      ldipole, lquad,   lfield,  ltorque,  leold1, &
      leold2,  leold3,  ldipvel

#define BC_INDDIPR 4
_REAL_ diptol,dipmass, diptau,  diptemp
common/inddipr/ &
      diptol,dipmass, diptau,  diptemp

#define BC_INDDIPI 6
integer maxiter, indmeth, irstdip, iquench, nquench, &
      nttdip
common/inddipi/ &
      maxiter, indmeth, irstdip, iquench, nquench, &
      nttdip

!-----------  END   ew_mpole.h --------------------------------------------
