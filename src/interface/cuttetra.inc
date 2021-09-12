      interface
      subroutine cuttetra(wtetra,vtetra,ntetra,etetra_in,efermi,tetra)
      use, intrinsic :: iso_fortran_env !inc
      implicit none
      real_dp,           intent(in)  :: efermi,etetra_in(4)
      real_dp,           intent(out) :: wtetra(4,3),vtetra(3)
      real_dp, optional, intent(out) :: tetra(4,4,3)
      integer,           intent(out) :: ntetra
      end subroutine cuttetra
      end interface
