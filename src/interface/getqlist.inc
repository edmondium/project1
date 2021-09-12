      interface
      subroutine getqlist(qvec,nq,mkpt_def,qfile)
      use, intrinsic :: iso_fortran_env !inc
      implicit none
      integer, intent(in)               :: mkpt_def
      integer, intent(out)              :: nq
      real_dp, intent(out), allocatable :: qvec(:,:)
      character(*), intent(in)          :: qfile
      end subroutine getqlist
      end interface
