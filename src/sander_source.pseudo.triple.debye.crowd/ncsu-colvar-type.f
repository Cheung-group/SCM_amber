! <compile=optimized>

#include "ncsu-utils.h"
#include "ncsu-config.h"

module ncsu_colvar_type

implicit none

private

integer, public, parameter :: COLVAR_ANGLE           = 1
integer, public, parameter :: COLVAR_TORSION         = 2
integer, public, parameter :: COLVAR_DISTANCE        = 3
integer, public, parameter :: COLVAR_MULTI_RMSD      = 4
integer, public, parameter :: COLVAR_R_OF_GYRATION   = 5
integer, public, parameter :: COLVAR_HANDEDNESS      = 6
integer, public, parameter :: COLVAR_N_OF_BONDS      = 7
integer, public, parameter :: COLVAR_N_OF_STRUCTURES = 8

type, public :: colvar_t

   integer :: type = -1

   integer,   pointer :: i(:) => null()
   NCSU_REAL, pointer :: r(:) => null()

   integer :: tag ! (see ncsu-cv-priv.*)

end type colvar_t

end module ncsu_colvar_type
