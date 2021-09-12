      interface
      subroutine orbitals(lwan,lmwan,centwan,eulerwan,nwan,charr,ncharr)
      use, intrinsic :: iso_fortran_env !inc
      implicit none
      real_dp,      intent(out), allocatable :: eulerwan(:,:)
      integer,      intent(out), allocatable :: lwan(:),lmwan(:),centwan(:)
      integer,      intent(out)              :: nwan
      integer,      intent(in)               :: ncharr
      character(*), intent(in)               :: charr(ncharr)
      end subroutine orbitals
      end interface
