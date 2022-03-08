! <compile=optimized>

#include "ncsu-utils.h"
#include "ncsu-config.h"

module ncsu_colvar

implicit none

private

!=============================================================================

public :: colvar_value
public :: colvar_force

public :: colvar_difference

public :: colvar_is_periodic

public :: colvar_has_min
public :: colvar_min

public :: colvar_has_max
public :: colvar_max

public :: colvar_print
public :: colvar_cleanup

public :: colvar_mdread
public :: colvar_bootstrap

!=============================================================================

contains

!=============================================================================

! the value is needed only on master
function colvar_value(cv, x) result(value)

   NCSU_USE_AFAILED

   use ncsu_colvar_type

   use ncsu_cv_ANGLE,           only : v_ANGLE           => colvar_value
   use ncsu_cv_TORSION,         only : v_TORSION         => colvar_value
   use ncsu_cv_DISTANCE,        only : v_DISTANCE        => colvar_value
   use ncsu_cv_MULTI_RMSD,      only : v_MULTI_RMSD      => colvar_value
   use ncsu_cv_R_OF_GYRATION,   only : v_R_OF_GYRATION   => colvar_value
   use ncsu_cv_HANDEDNESS,      only : v_HANDEDNESS      => colvar_value
   use ncsu_cv_N_OF_BONDS,      only : v_N_OF_BONDS      => colvar_value
   use ncsu_cv_N_OF_STRUCTURES, only : v_N_OF_STRUCTURES => colvar_value

   implicit none

   NCSU_REAL :: value

   type(colvar_t) :: cv ! mutable
   NCSU_REAL, intent(in) :: x(*)

   select case(cv%type)
      case(COLVAR_ANGLE)
         value = v_ANGLE(cv, x)
      case(COLVAR_TORSION)
         value = v_TORSION(cv, x)
      case(COLVAR_DISTANCE)
         value = v_DISTANCE(cv, x)
      case(COLVAR_MULTI_RMSD)
         value = v_MULTI_RMSD(cv, x)
      case(COLVAR_R_OF_GYRATION)
         value = v_R_OF_GYRATION(cv, x)
      case(COLVAR_HANDEDNESS)
         value = v_HANDEDNESS(cv, x)
      case(COLVAR_N_OF_BONDS)
         value = v_N_OF_BONDS(cv, x)
      case(COLVAR_N_OF_STRUCTURES)
         value = v_N_OF_STRUCTURES(cv, x)
      case default
         ncsu_assert_not_reached()
         value = NCSU_TO_REAL(0)
   end select

end function colvar_value

!=============================================================================

subroutine colvar_force(cv, x, fcv, f)

   NCSU_USE_AFAILED

   use ncsu_colvar_type

   use ncsu_cv_ANGLE,           only : f_ANGLE           => colvar_force
   use ncsu_cv_TORSION,         only : f_TORSION         => colvar_force
   use ncsu_cv_DISTANCE,        only : f_DISTANCE        => colvar_force
   use ncsu_cv_MULTI_RMSD,      only : f_MULTI_RMSD      => colvar_force
   use ncsu_cv_R_OF_GYRATION,   only : f_R_OF_GYRATION   => colvar_force
   use ncsu_cv_HANDEDNESS,      only : f_HANDEDNESS      => colvar_force
   use ncsu_cv_N_OF_BONDS,      only : f_N_OF_BONDS      => colvar_force
   use ncsu_cv_N_OF_STRUCTURES, only : f_N_OF_STRUCTURES => colvar_force

   implicit none

   type(colvar_t) :: cv ! mutable

   NCSU_REAL, intent(in) :: x(*), fcv
   NCSU_REAL, intent(inout) :: f(*)

   select case(cv%type)
      case(COLVAR_ANGLE)
         call f_ANGLE(cv, x, fcv, f)
      case(COLVAR_TORSION)
         call f_TORSION(cv, x, fcv, f)
      case(COLVAR_DISTANCE)
         call f_DISTANCE(cv, x, fcv, f)
      case(COLVAR_MULTI_RMSD)
         call f_MULTI_RMSD(cv, x, fcv, f)
      case(COLVAR_R_OF_GYRATION)
         call f_R_OF_GYRATION(cv, x, fcv, f)
      case(COLVAR_HANDEDNESS)
         call f_HANDEDNESS(cv, x, fcv, f)
      case(COLVAR_N_OF_BONDS)
         call f_N_OF_BONDS(cv, x, fcv, f)
      case(COLVAR_N_OF_STRUCTURES)
         call f_N_OF_STRUCTURES(cv, x, fcv, f)
      case default
         ncsu_assert_not_reached()
   end select

end subroutine colvar_force

!=============================================================================

function colvar_difference(cv, v1, v2) result(diff)

   NCSU_USE_AFAILED

   use ncsu_colvar_type
   use ncsu_constants, only : ZERO, PI

   implicit none

   NCSU_REAL :: diff
   type(colvar_t), intent(in) :: cv
   NCSU_REAL, intent(in) :: v1, v2

   NCSU_REAL :: t1, t2

   t1 = fix_value(v1)
   t2 = fix_value(v2)

   diff = t1 - t2

   select case(cv%type)
      case(COLVAR_TORSION)
         ncsu_assert(- PI.le.t1.and.t1.le.PI)
         ncsu_assert(- PI.le.t2.and.t2.le.PI)
         if (diff.gt.PI) then
            diff = diff - PI - PI
         else if (diff.lt.-PI) then
            diff = diff + PI + PI
         end if
      case default
         continue
   end select

contains

function fix_value(v) result(t)

   implicit none

   NCSU_REAL :: t
   NCSU_REAL, intent(in) :: v

   select case(cv%type)
      case(COLVAR_ANGLE)
         if (ZERO.le.v.and.v.le.PI) then
            t = v
         else
            t = acos(cos(v))
         end if
      case(COLVAR_TORSION)
         if (- PI.le.v.and.v.le.PI) then
	         t = v
         else
            t = atan2(sin(v), cos(v))
         endif
      case default
         t = v
   end select

end function fix_value

end function colvar_difference

!=============================================================================

logical function colvar_is_periodic(cv)

   use ncsu_colvar_type

   implicit none

   type(colvar_t), intent(in) :: cv

   select case(cv%type)
      case(COLVAR_TORSION)
         colvar_is_periodic = .true.
      case default
         colvar_is_periodic = .false.
   end select

end function colvar_is_periodic

!=============================================================================

logical function colvar_has_min(cv)

   use ncsu_colvar_type

   implicit none

   type(colvar_t), intent(in) :: cv

   select case(cv%type)
      case(COLVAR_ANGLE:COLVAR_R_OF_GYRATION)
         colvar_has_min = .true.
      case(COLVAR_N_OF_BONDS:COLVAR_N_OF_STRUCTURES)
         colvar_has_min = .true.
      case default
         colvar_has_min = .false.
   end select

end function colvar_has_min

!=============================================================================

NCSU_REAL function colvar_min(cv)

   NCSU_USE_AFAILED

   use ncsu_constants
   use ncsu_colvar_type

   implicit none

   type(colvar_t), intent(in) :: cv

   ncsu_assert(colvar_has_min(cv))

   select case(cv%type)
      case(COLVAR_ANGLE)
         colvar_min = ZERO
      case(COLVAR_TORSION)
         colvar_min = -PI
      case(COLVAR_DISTANCE:COLVAR_R_OF_GYRATION)
         colvar_min = ZERO
      case(COLVAR_N_OF_BONDS:COLVAR_N_OF_STRUCTURES)
         colvar_min = ZERO
      case default
         ncsu_assert_not_reached()
         colvar_min = ZERO
   end select

end function colvar_min

!=============================================================================

logical function colvar_has_max(cv)

   use ncsu_colvar_type

   implicit none

   type(colvar_t), intent(in) :: cv

   select case(cv%type)
      case(COLVAR_ANGLE:COLVAR_TORSION)
         colvar_has_max = .true.
      case default
         colvar_has_max = .false.
   end select

end function colvar_has_max

!=============================================================================

NCSU_REAL function colvar_max(cv)

   NCSU_USE_AFAILED

   use ncsu_constants
   use ncsu_colvar_type

   implicit none

   type(colvar_t), intent(in) :: cv

   ncsu_assert(colvar_has_max(cv))

   select case(cv%type)
      case(COLVAR_ANGLE:COLVAR_TORSION)
         colvar_max = PI
      case default
         ncsu_assert_not_reached()
         colvar_max = ZERO
   end select

end function colvar_max

!=============================================================================

subroutine colvar_print(cv, lun)

   use ncsu_utils
   use ncsu_constants
   use ncsu_colvar_type
   use ncsu_sander_proxy

   use ncsu_cv_ANGLE,           only : p_ANGLE           => print_details
   use ncsu_cv_TORSION,         only : p_TORSION         => print_details
   use ncsu_cv_DISTANCE,        only : p_DISTANCE        => print_details
   use ncsu_cv_MULTI_RMSD,      only : p_MULTI_RMSD      => print_details
   use ncsu_cv_R_OF_GYRATION,   only : p_R_OF_GYRATION   => print_details
   use ncsu_cv_HANDEDNESS,      only : p_HANDEDNESS      => print_details
   use ncsu_cv_N_OF_BONDS,      only : p_N_OF_BONDS      => print_details
   use ncsu_cv_N_OF_STRUCTURES, only : p_N_OF_STRUCTURES => print_details

   implicit none

   type(colvar_t), intent(in) :: cv
   integer,        intent(in) :: lun

   write (unit = lun, fmt = '(a,a)', advance = 'NO') NCSU_INFO, '  type = '''

   select case(cv%type)
      case(COLVAR_ANGLE)
         write (unit = lun, fmt = '(a)') 'ANGLE'''
         call p_ANGLE(cv, lun)
      case(COLVAR_TORSION)
         write (unit = lun, fmt = '(a)') 'TORSION'''
         call p_TORSION(cv, lun)
      case(COLVAR_DISTANCE)
         write (unit = lun, fmt = '(a)') 'DISTANCE'''
         call p_DISTANCE(cv, lun)
      case(COLVAR_MULTI_RMSD)
         write (unit = lun, fmt = '(a)') 'MULTI_RMSD'''
         call p_MULTI_RMSD(cv, lun)
      case(COLVAR_R_OF_GYRATION)
         write (unit = lun, fmt = '(a)') 'R_OF_GYRATION'''
         call p_R_OF_GYRATION(cv, lun)
      case(COLVAR_HANDEDNESS)
         write (unit = lun, fmt = '(a)') 'HANDEDNESS'''
         call p_HANDEDNESS(cv, lun)
      case(COLVAR_N_OF_BONDS)
         write (unit = lun, fmt = '(a)') 'N_OF_BONDS'''
         call p_N_OF_BONDS(cv, lun)
      case(COLVAR_N_OF_STRUCTURES)
         write (unit = lun, fmt = '(a)') 'N_OF_STRUCTURES'''
         call p_N_OF_STRUCTURES(cv, lun)
      case default
         ncsu_assert_not_reached()
         continue
   end select

end subroutine colvar_print

!=============================================================================

subroutine colvar_cleanup(cv)

   use ncsu_colvar_type

   use ncsu_cv_MULTI_RMSD,      only : c_MULTI_RMSD      => colvar_cleanup
   use ncsu_cv_R_OF_GYRATION,   only : c_R_OF_GYRATION   => colvar_cleanup
   use ncsu_cv_N_OF_STRUCTURES, only : c_N_OF_STRUCTURES => colvar_cleanup

   implicit none

   type(colvar_t), intent(inout) :: cv

   select case(cv%type)
      case(COLVAR_MULTI_RMSD)
         call c_MULTI_RMSD(cv)
      case(COLVAR_R_OF_GYRATION)
         call c_R_OF_GYRATION(cv)
      case(COLVAR_N_OF_STRUCTURES)
         call c_N_OF_STRUCTURES(cv)
      case default
         continue
   end select

   if (associated(cv%i)) &
      deallocate(cv%i)

   if (associated(cv%r)) &
      deallocate(cv%r)

   cv%type = -1

end subroutine colvar_cleanup

!=============================================================================

subroutine colvar_mdread(cv, node, cvno)

   use ncsu_utils
   use ncsu_value
   use ncsu_cftree
   use ncsu_constants
   use ncsu_colvar_type
   use ncsu_sander_proxy

   implicit none

   type(colvar_t), intent(inout) :: cv
   type(node_t),   intent(in)    :: node
   integer,        intent(in)    :: cvno

   integer :: n, error

   type(value_node_t), pointer :: alist, aiter

   character(len = STRING_LENGTH) :: type

   ncsu_assert(is_master())

   ncsu_assert(.not. associated(cv%i))
   ncsu_assert(.not. associated(cv%r))

   ncsu_assert(node_title(node) == 'variable')

   !
   ! type
   !

   if (.not.node_lookup_string(node, 'type', type)) then
      write (unit = ERR_UNIT, fmt = '(/a,a,'//pfmt(cvno)//'/)') &
            NCSU_ERROR, 'type is not specified for CV #', cvno
      call terminate()
   end if

   if (type == 'ANGLE') then
      cv%type = COLVAR_ANGLE
   else if (type == 'TORSION') then
      cv%type = COLVAR_TORSION
   else if (type == 'DISTANCE') then
      cv%type = COLVAR_DISTANCE
   else if (type == 'MULTI_RMSD') then
      cv%type = COLVAR_MULTI_RMSD
   else if (type == 'R_OF_GYRATION') then
      cv%type = COLVAR_R_OF_GYRATION
   else if (type == 'HANDEDNESS') then
      cv%type = COLVAR_HANDEDNESS
   else if (type == 'N_OF_BONDS') then
      cv%type = COLVAR_N_OF_BONDS
   else if (type == 'N_OF_STRUCTURES') then
      cv%type = COLVAR_N_OF_STRUCTURES
   else
      write (unit = ERR_UNIT, fmt = '(/a,a,a,a,'//pfmt(cvno)//',a/)') &
            NCSU_ERROR, 'CV type ''', trim(type), &
            ''' is not supported so far (CV #', cvno, ')'
      call terminate()
   end if

   !
   ! cv%i
   !

   if (node_lookup_list(node, 'i', alist)) then

      n = 0
      aiter => alist

      do while (associated(aiter))
         n = n + 1
         if (.not. value_is_integer(aiter%value)) then
            write (unit = ERR_UNIT, fmt = '(/a,a,'//pfmt(cvno)//',a/)') &
               NCSU_ERROR, 'CV #', cvno, ' : unexpected &
               &(not an integer) element of ''i'' list'
            call terminate()
         end if
         aiter => aiter%next
      end do

      if (n > 0) then
         allocate(cv%i(n), stat = error)
         if (error /= 0) &
            NCSU_OUT_OF_MEMORY

         n = 0
         aiter => alist
         do while (associated(aiter))
            n = n + 1
            cv%i(n) = value_get_integer(aiter%value)
            aiter => aiter%next
         end do
      end if

   end if ! node_lookup_list(vnode, 'i', alist))

   !
   ! cv%r
   !

   if (node_lookup_list(node, 'r', alist)) then

      n = 0
      aiter => alist

      do while (associated(aiter))
         n = n + 1
         if (.not. value_is_real(aiter%value)) then
            write (unit = ERR_UNIT, fmt = '(/a,a,'//pfmt(cvno)//',a/)') &
               NCSU_ERROR, 'CV #', cvno, ' : unexpected &
               &(not a real number) element of ''r'' list'
            call terminate()
         end if
         aiter => aiter%next
      end do

      if (n > 0) then
         allocate(cv%r(n), stat = error)
         if (error /= 0) &
            NCSU_OUT_OF_MEMORY

         n = 0
         aiter => alist
         do while (associated(aiter))
            n = n + 1
            cv%r(n) = value_get_real(aiter%value)
            aiter => aiter%next
         end do
      end if

   end if ! node_lookup_list(vnode, 'r', alist))

end subroutine colvar_mdread

!=============================================================================

subroutine colvar_bootstrap(cv, cvno, amass)

   use ncsu_utils
   use ncsu_colvar_type

   use ncsu_cv_ANGLE,           only : b_ANGLE           => colvar_bootstrap
   use ncsu_cv_TORSION,         only : b_TORSION         => colvar_bootstrap
   use ncsu_cv_DISTANCE,        only : b_DISTANCE        => colvar_bootstrap
   use ncsu_cv_MULTI_RMSD,      only : b_MULTI_RMSD      => colvar_bootstrap
   use ncsu_cv_R_OF_GYRATION,   only : b_R_OF_GYRATION   => colvar_bootstrap
   use ncsu_cv_HANDEDNESS,      only : b_HANDEDNESS      => colvar_bootstrap
   use ncsu_cv_N_OF_BONDS,      only : b_N_OF_BONDS      => colvar_bootstrap
   use ncsu_cv_N_OF_STRUCTURES, only : b_N_OF_STRUCTURES => colvar_bootstrap

   implicit none

   type(colvar_t), intent(inout) :: cv
   integer,        intent(in)    :: cvno
   NCSU_REAL,      intent(in)    :: amass(*)

#ifdef MPI
   integer :: bcastdata(3), ierr
#include "ncsu-mpi.h"

   !
   ! bcast type/i/r first
   !

#ifndef NCSU_DISABLE_ASSERT
   if (sanderrank == 0) then
      ncsu_assert(cv%type > 0)
   else
      ncsu_assert(.not. associated(cv%i))
      ncsu_assert(.not. associated(cv%r))
   end if
#endif /* NCSU_DISABLE_ASSERT */

   if (sanderrank == 0) then
      bcastdata(1) = cv%type

      bcastdata(2) = 0
      if (associated(cv%i)) &
         bcastdata(2) = size(cv%i)

      bcastdata(3) = 0
      if (associated(cv%r)) &
         bcastdata(3) = size(cv%r)
   end if ! sanderrank == 0

   call mpi_bcast(bcastdata, size(bcastdata), MPI_INTEGER, 0, commsander, ierr)
   ncsu_assert(ierr == 0)

   if (sanderrank /= 0) &
      cv%type = bcastdata(1)

   !
   ! cv%i
   !

   if (bcastdata(2) > 0) then
      if (.not. associated(cv%i)) then
         ncsu_assert(sanderrank > 0)
         allocate(cv%i(bcastdata(2)), stat = ierr)
         if (ierr /= 0) &
            NCSU_OUT_OF_MEMORY
      end if ! .not. associated(cv%i)

      call mpi_bcast(cv%i, bcastdata(2), MPI_INTEGER, 0, commsander, ierr)
      ncsu_assert(ierr == 0)
   else
      nullify(cv%i)
   end if ! bcastdata(2) > 0

   !
   ! cv%r
   !

   if (bcastdata(3) > 0) then
      if (.not. associated(cv%r)) then
         ncsu_assert(sanderrank > 0)
         allocate(cv%r(bcastdata(3)), stat = ierr)
         if (ierr /= 0) &
            NCSU_OUT_OF_MEMORY
      end if ! .not. associated(cv%r)

      call mpi_bcast(cv%r, bcastdata(3), MPI_DOUBLE_PRECISION, &
                     0, commsander, ierr)
      ncsu_assert(ierr == 0)
   else
      nullify(cv%r)
   end if ! bcastdata(3) > 0
#endif /* MPI */

   !
   ! dispatch according to the type
   !

   select case(cv%type)
      case(COLVAR_ANGLE)
         call b_ANGLE(cv, cvno, amass)
      case(COLVAR_TORSION)
         call b_TORSION(cv, cvno, amass)
      case(COLVAR_DISTANCE)
         call b_DISTANCE(cv, cvno, amass)
      case(COLVAR_MULTI_RMSD)
         call b_MULTI_RMSD(cv, cvno, amass)
      case(COLVAR_R_OF_GYRATION)
         call b_R_OF_GYRATION(cv, cvno, amass)
      case(COLVAR_HANDEDNESS)
         call b_HANDEDNESS(cv, cvno, amass)
      case(COLVAR_N_OF_BONDS)
         call b_N_OF_BONDS(cv, cvno, amass)
      case(COLVAR_N_OF_STRUCTURES)
         call b_N_OF_STRUCTURES(cv, cvno, amass)
      case default
         ncsu_assert_not_reached()
         continue
   end select

end subroutine colvar_bootstrap

!=============================================================================

end module ncsu_colvar
