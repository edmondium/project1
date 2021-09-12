      interface
      subroutine selfenergy(job1,ikpt,eval,nbas,coul,Win(ctrafo),Win(screen), head, wing,
     &                                                           Win(screenc),headc,wingc,plasma)

      use global !inc
      use arrays !inc

      implicit none
      type(jobtype), intent(in)           :: job1
      integer,       intent(in)           :: nbas,ikpt
      logical,       intent(in)           :: eval(2)
      real_dp,       intent(in), optional :: coul(nbas)
      complex_dp,    intent(in), optional :: headc(3,3,nfreqc),wingc(3,nbas,2,nfreqc)
      MCOMPLEX_dp,   intent(in), optional :: head(6,nfreq),wing(3,nbas,nfreq)
      MCOMPLEX_dp,               optional :: ctrafo(nbasm(ikpt),nbas)      ! can be modified but leaves subroutine unchanged
      MCOMPLEX_dp,               optional :: screen(nbas*(nbas+1)/2,nfreq) ! modified by cblock
      complex_dp,                optional :: screenc(nbas,nbas,nfreqc)     ! modified by cblock
      real_dp,                   optional :: plasma                        ! changed if plasma=-1d0 (PLASMA METAL)
# ifdef MPI
      integer,       intent(in), optional :: win_screen,win_screenc,win_ctrafo
# endif
      end subroutine selfenergy
      end interface
