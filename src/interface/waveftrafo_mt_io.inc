      interface
      subroutine waveftrafo_mt_io(cmtout,cmtin,dimout,dimin,ic,isym,ikpt00)


      use, intrinsic :: iso_fortran_env !inc
      implicit none
      integer,    intent(in)         :: ikpt00,isym,ic,dimout,dimin
      complex_dp, intent(in), target :: cmtin(dimin,*)
      complex_dp, intent(inout)      :: cmtout(dimout,*)
      end subroutine waveftrafo_mt_io
      end interface
