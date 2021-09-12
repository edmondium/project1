      interface
      subroutine wannier_spatial(matr,matk,lpt,nlpt,ldeg,ndiff,ispin,mode)
      use, intrinsic :: iso_fortran_env !inc
      implicit none
      integer,                 intent(inout) :: nlpt,ndiff
      integer,    allocatable, intent(inout) :: lpt(:,:),ldeg(:,:)
      integer,                 intent(in)    :: ispin,mode
      complex_dp,              intent(in)    :: matk(*)
      complex_dp, allocatable, intent(inout) :: matr(:,:,:)
      end subroutine wannier_spatial
      end interface
