#include "ncsu-utils.h"
#include "ncsu-config.h"

module ncsu_colvar_utils

implicit none

private

!=============================================================================

public :: check_i
public :: print_i

!=============================================================================

contains

!=============================================================================

subroutine check_i(cvi, cvno, cvtype, expected_isize)

   use ncsu_utils
   use ncsu_constants
   use ncsu_sander_proxy

   implicit none

   integer, pointer :: cvi(:)
   integer, intent(in) :: cvno
   character(*), intent(in) :: cvtype
   integer, optional, intent(in) :: expected_isize

   integer :: a, b

#include "ncsu-mpi.h"

   ncsu_assert(cvno > 0)

   if (.not.associated(cvi)) then
      NCSU_MASTER_ONLY_BEGIN
         write (unit = ERR_UNIT, fmt = '(/a,a,'//pfmt(cvno)//',a,a,a/)') &
            NCSU_ERROR, 'CV #', cvno, ' (', cvtype, &
            ') : no integers found'
      NCSU_MASTER_ONLY_END
      call terminate()
   end if ! .not. associated(cvi)

   if (present(expected_isize)) then
      if (size(cvi) /= expected_isize) then
         NCSU_MASTER_ONLY_BEGIN
            write (unit = ERR_UNIT, &
            fmt = '(/a,a,'//pfmt(cvno)//',a,a,a,'//pfmt &
            (size(cvi))//',a,'//pfmt(expected_isize)//',a/)') &
               NCSU_ERROR, 'CV #', cvno, &
               ' (', cvtype, ') : unexpected number of integers (', &
               size(cvi), ' instead of ', expected_isize, ')'
         NCSU_MASTER_ONLY_END
         call terminate()
      end if ! size(cvi) /= isize
   end if ! present(expected_isize)

   do a = 1, size(cvi)
      if (cvi(a) < 1 .or. cvi(a) > sander_natoms()) then
         NCSU_MASTER_ONLY_BEGIN
            write (unit = ERR_UNIT, &
               fmt = '(/a,a,'//pfmt(cvno)//',a,a,a,'//pfmt(a)//',a,'//pfmt &
               (cvi(a))//',a,'//pfmt(sander_natoms())//',a/)') &
               NCSU_ERROR, 'CV #', cvno, &
               ' (', cvtype, ') : integer #', a, ' (', cvi(a), &
               ') is out of range [1, ', sander_natoms(), ']'
         NCSU_MASTER_ONLY_END
         call terminate()
      end if
   end do

   ! check for duplicates

   do a = 1, size(cvi)
      do b = a + 1, size(cvi)
         if (cvi(a) == cvi(b)) then
            NCSU_MASTER_ONLY_BEGIN
               write (unit = ERR_UNIT, &
                  fmt = '(/a,a,'//pfmt(cvno)//',a,a,a,'//pfmt &
                  (a)//',a,'//pfmt(b)//',a,'//pfmt(cvi(a))//',a/)') &
                  NCSU_ERROR, 'CV #', cvno, ' (', cvtype, &
                  ') : integers #', a, ' and #', b, ' are equal (', cvi(a), ')'
            NCSU_MASTER_ONLY_END
            call terminate()
         end if ! cvi(a) == cvi(b)
      end do
   end do

end subroutine check_i

!=============================================================================

subroutine print_i(cvi, lun)

   use ncsu_utils
   use ncsu_sander_proxy

   implicit none

   integer, pointer :: cvi(:)
   integer, intent(in) :: lun

   integer :: a
   character(4) :: aname

   ncsu_assert(is_master())
   ncsu_assert(associated(cvi))
   ncsu_assert(size(cvi) > 0)

   write (unit = lun, fmt = '(a,a)', advance = 'NO') NCSU_INFO, '  atoms = ('

   do a = 1, size(cvi)

      ncsu_assert(cvi(a) > 0 .and. cvi(a) <= sander_natoms())
      aname = sander_atom_name(cvi(a))

      write (unit = lun, fmt = '('//pfmt(cvi(a))//',a,a,a)', advance = 'NO') &
         cvi(a), ' [', trim(aname), ']'

      if (a == size(cvi)) then
         write (unit = lun, fmt = '(a)') ')'
      else if (mod(a, 5) == 0) then
         write (unit = lun, fmt = '(a,/a,a)', advance = 'NO') &
            ',', NCSU_INFO, '          '
      else
         write (unit = lun, fmt = '(a)', advance = 'NO') ', '
      end if

   end do

end subroutine print_i

!=============================================================================

end module ncsu_colvar_utils
