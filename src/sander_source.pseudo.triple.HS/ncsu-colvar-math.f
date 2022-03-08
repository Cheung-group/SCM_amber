! <compile=optimized>

#include "ncsu-config.h"

module ncsu_colvar_math

implicit none

private

public :: distance
public :: distance_d

public :: angle
public :: angle_d

public :: torsion
public :: torsion_d

private :: dot3
private :: norm3
private :: cross3

private :: cosine
private :: cosine_d

private :: torsion3
private :: torsion3_d

!=============================================================================

contains

!=============================================================================

pure NCSU_REAL function distance(r1, r2)

   implicit none

   NCSU_REAL, intent(in) :: r1(3), r2(3)

   distance = norm3(r1 - r2)

end function distance

!=============================================================================

subroutine distance_d(r1, r2, d1, d2)

   implicit none

   NCSU_REAL, intent(in)  :: r1(3), r2(3)
   NCSU_REAL, intent(out) :: d1(3), d2(3)

   NCSU_REAL :: dr(3)

   dr = r1 - r2

   d1 = dr/norm3(dr)
   d2 = - d1

end subroutine distance_d

!=============================================================================

pure NCSU_REAL function angle(r1, r2, r3)

   implicit none

   NCSU_REAL, intent(in) :: r1(3), r2(3), r3(3)

   angle = acos(cosine(r1 - r2, r3 - r2))

end function angle

!=============================================================================

subroutine angle_d(r1, r2, r3, d1, d2, d3)

   use ncsu_constants, only : ONE

   implicit none

   NCSU_REAL, intent(in)  :: r1(3), r2(3), r3(3)
   NCSU_REAL, intent(out) :: d1(3), d2(3), d3(3)

   NCSU_REAL :: c, d

   c = cosine_d(r1 - r2, r3 - r2, d1, d3)
   d = -ONE/sqrt(ONE - c*c)

   d1 = d*d1
   d3 = d*d3

   d2 = - (d1 + d3)

end subroutine angle_d

!=============================================================================

!
! for A == B == C == D, a torsion along B == C is arccos([ABxBC]*[CDx(-BC)])
!        and its sign is given by sign(BC*[[ABxBC]x[CDx(-BC)]])
!

pure NCSU_REAL function torsion(r1, r2, r3, r4)

   implicit none

   NCSU_REAL, intent(in) :: r1(3), r2(3), r3(3), r4(3)

   torsion = torsion3(r2 - r1, r3 - r2, r4 - r3)

end function torsion

!=============================================================================

subroutine torsion_d(r1, r2, r3, r4, d1, d2, d3, d4)

   implicit none

   NCSU_REAL, intent(in)  :: r1(3), r2(3), r3(3), r4(3)
   NCSU_REAL, intent(out) :: d1(3), d2(3), d3(3), d4(3)

   NCSU_REAL :: t(3)

   call torsion3_d(r2 - r1, r3 - r2, r4 - r3, d1, t, d4)

   d2 = d1 - t
   d1 = - d1
   d3 = t - d4

end subroutine torsion_d

!=============================================================================

pure NCSU_REAL function dot3(v1, v2)

   implicit none

   NCSU_REAL, intent(in) :: v1(3), v2(3)

   dot3 = v1(1)*v2(1) + v1(2)*v2(2) + v1(3)*v2(3)

end function dot3

!=============================================================================

pure NCSU_REAL function norm3(v)

   implicit none

   NCSU_REAL, intent(in) :: v(3)

   norm3 = sqrt(v(1)**2 + v(2)**2 + v(3)**2)

end function norm3

!=============================================================================

pure function cross3(v1, v2) result(cross)

   implicit none

   NCSU_REAL             :: cross(3)
   NCSU_REAL, intent(in) :: v1(3), v2(3)

   cross(1) = v1(2)*v2(3) - v1(3)*v2(2)
   cross(2) = v1(3)*v2(1) - v1(1)*v2(3)
   cross(3) = v1(1)*v2(2) - v1(2)*v2(1)

end function cross3

!=============================================================================

pure NCSU_REAL function cosine(v1, v2)

   use ncsu_constants, only : ONE

   implicit none

   NCSU_REAL, intent(in) :: v1(3), v2(3)

   cosine = dot3(v1, v2)/(norm3(v1)*norm3(v2))

   if (cosine > ONE) cosine = ONE
   if (cosine < -ONE) cosine = -ONE

end function cosine

!=============================================================================

NCSU_REAL function cosine_d(v1, v2, d1, d2)

   use ncsu_constants, only : ONE

   implicit none

   NCSU_REAL, intent(in)  :: v1(3), v2(3)
   NCSU_REAL, intent(out) :: d1(3), d2(3)

   NCSU_REAL :: n1, n2, p12

   n1  = norm3(v1)
   n2  = norm3(v2)
   p12 = dot3(v1, v2)

   cosine_d = p12/(n1*n2)

   d1 = (v2 - v1*p12/(n1**2))/(n1*n2)
   d2 = (v1 - v2*p12/(n2**2))/(n1*n2)

   if (cosine_d > ONE) cosine_d = ONE
   if (cosine_d < -ONE) cosine_d = -ONE

end function cosine_d

!=============================================================================

pure NCSU_REAL function torsion3(v1, v2, v3)

   implicit none

   NCSU_REAL, intent(in) :: v1(3), v2(3), v3(3)

   NCSU_REAL :: n1(3), n2(3)

   n1 = cross3(v1, v2)
   n2 = cross3(v2, v3)

   torsion3 = sign(acos(cosine(n1, n2)), dot3(v2, cross3(n1, n2)))

end function torsion3

!=============================================================================

subroutine torsion3_d(v1, v2, v3, d1, d2, d3)

   use ncsu_constants, only : ONE

   implicit none

   NCSU_REAL, intent(in)  :: v1(3), v2(3), v3(3)
   NCSU_REAL, intent(out) :: d1(3), d2(3), d3(3)

   NCSU_REAL, parameter :: TINY = (ONE/10)**8

   NCSU_REAL :: n1(3), n2(3), dc1(3), dc2(3), c, c2, s, d

   n1 = cross3(v1, v2)
   n2 = cross3(v2, v3)

   c = cosine_d(n1, n2, dc1, dc2)
   s = dot3(v2, cross3(n1, n2))

   c2 = c*c - TINY
   d = sign(ONE/sqrt(ONE - c2), s)

   d1 = d*cross3(dc1, v2)
   d2 = d*(cross3(v1, dc1) + cross3(dc2, v3))
   d3 = d*cross3(v2, dc2)

end subroutine torsion3_d

!=============================================================================

end module ncsu_colvar_math
