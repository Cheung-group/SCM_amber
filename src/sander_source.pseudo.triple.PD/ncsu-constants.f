#include "ncsu-config.h"

module ncsu_constants

implicit none

private

#ifdef NCSU_REAL_IS_DOUBLE
NCSU_REAL, public, parameter :: PI = 3.14159265358979323846D0
NCSU_REAL, public, parameter :: ZERO  = 0.d0
NCSU_REAL, public, parameter :: ONE   = 1.d0
NCSU_REAL, public, parameter :: TWO   = 2.d0
NCSU_REAL, public, parameter :: THREE = 3.d0
NCSU_REAL, public, parameter :: FOUR  = 4.d0
#else
NCSU_REAL, public, parameter :: PI = 3.1415926535
NCSU_REAL, public, parameter :: ZERO  = 0.0
NCSU_REAL, public, parameter :: ONE   = 1.0
NCSU_REAL, public, parameter :: TWO   = 2.0
NCSU_REAL, public, parameter :: THREE = 3.0
NCSU_REAL, public, parameter :: FOUR  = 4.0
#endif /* NCSU_REAL_IS_DOUBLE */

integer, public, parameter :: STRING_LENGTH = 256

integer, public, parameter :: ERR_UNIT = SANDER_STDERR_UNIT
integer, public, parameter :: OUT_UNIT = SANDER_STDOUT_UNIT

integer, public, parameter :: ABMD_MONITOR_UNIT = SANDER_LAST_UNIT + 1
integer, public, parameter :: SMD_OUTPUT_UNIT = SANDER_LAST_UNIT + 2
integer, public, parameter :: PMD_OUTPUT_UNIT = SANDER_LAST_UNIT + 3

#ifdef MPI
integer, public, parameter :: REM_MDIN_UNIT = SANDER_LAST_UNIT + 4
integer, public, parameter :: PMD_REMLOG_UNIT = SANDER_LAST_UNIT + 5
integer, public, parameter :: ABMD_REMLOG_UNIT = SANDER_LAST_UNIT + 6
integer, public, parameter :: BBMD_MONITOR_UNIT = SANDER_LAST_UNIT + 7
integer, public, parameter :: BBMD_LOG_UNIT = SANDER_LAST_UNIT + 8
#endif /* MPI */

! and so forth

end module ncsu_constants
