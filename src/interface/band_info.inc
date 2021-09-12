      interface
      subroutine band_info(charr,cmtin,spin,band,result,title)
      use, intrinsic :: iso_fortran_env !inc
      implicit none
      character(*), intent(in)                         :: charr(:)
      character(*), intent(in),               optional :: title      
      integer,      intent(in),               optional :: spin,band(:)
      real_dp,      intent(out), allocatable, optional :: result(:,:,:)
      complex_dp,   intent(in),               optional :: cmtin(:,:,:,:)
      end subroutine band_info
      end interface
