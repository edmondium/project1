      interface
      subroutine hdf_kptwrt(Hfile,eig,cmtin,cpwin,maxbnd,nbnd,ikpt,ikptp,ispin,gpt1)
      use global   !inc
      use hdf5     !inc
      implicit none
      integer,           intent(in) :: ikpt,ikptp,ispin,nbnd,maxbnd
      integer(HID_T),    intent(in) :: Hfile
      integer, optional, intent(in) :: gpt1(:,:,:)
      real_dp,           intent(in) :: eig(nbnd)
      complex_dp,        intent(in) :: cmtin(maxlmindx,ncent,maxbnd,nspin3)
      MCOMPLEX_dp,       intent(in) :: cpwin(maxgpt,maxbnd,nspin3)
      end subroutine hdf_kptwrt
      end interface
