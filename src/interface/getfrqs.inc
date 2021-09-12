      interface
      subroutine getfrqs(frqs,nfrqs,fspec,wmin,wmax,lwrite)
      use, intrinsic :: iso_fortran_env !inc
      implicit none
      integer, intent(out)              :: nfrqs
      real_dp, intent(out), allocatable :: frqs(:)
      real_dp, intent(in)               :: wmin,wmax,fspec(2)
      logical, intent(in)               :: lwrite
      end subroutine getfrqs
      end interface
