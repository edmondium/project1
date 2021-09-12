      interface
      subroutine define_freq(job1,contour_par,ospin,ltrs,lwrite)
      use global !inc
      use arrays !inc
      implicit none
      type(jobtype), intent(in)   :: job1
      real_dp,       intent(out)  :: contour_par(4)
      logical,       intent(in)   :: lwrite,ltrs
      integer,       intent(in)   :: ospin
      end subroutine define_freq
      end interface
