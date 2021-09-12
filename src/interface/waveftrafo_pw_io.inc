      interface
      subroutine waveftrafo_pw_io(cpwout,cpwin,isym,ikpt00)


      use, intrinsic :: iso_fortran_env !inc
      implicit none
      integer,       intent(in)         :: ikpt00,isym
      MCOMPLEX_dp,   intent(in), target :: cpwin(*)
      MCOMPLEX_dp,   intent(inout)      :: cpwout(*)
      end subroutine waveftrafo_pw_io
      end interface
