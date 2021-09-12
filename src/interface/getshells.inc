      interface
      subroutine getshells(ptsh,nptsh,radsh,nshell,rad,lat,shift,nshift,lwrite)

# ifdef MPI      
# endif

      use, intrinsic :: iso_fortran_env !inc
      implicit none
      logical,              intent(in)  :: lwrite
      integer,              intent(in)  :: nshift
      integer,              intent(out) :: nptsh,nshell
      integer, allocatable, intent(out) :: ptsh(:,:)
      real_dp, allocatable, intent(out) :: radsh(:)
      real_dp,              intent(in)  :: rad,lat(3,3),shift(3,nshift)
      end subroutine getshells
      end interface
