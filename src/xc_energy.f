c--------------------------------------------------------------------------------
c Copyright (c) 2021 Peter Grünberg Institut, Forschungszentrum Jülich, Germany
c This file is part of SPEX and available as free software under the conditions
c of the MIT license as expressed in the LICENSE file in more detail.
c--------------------------------------------------------------------------------

c     Write out RPA/HF energy

# include "cppmacro.h"
# include "jobtype.h"

      subroutine xc_energy(job1)
      use global
      use arrays
      use, intrinsic :: iso_fortran_env
      implicit none
      type(jobtype), intent(in) :: job1
      complex_dp                :: cdum
      real_dp                   :: rdum,rdum1,rdum2
      integer                   :: i
      real_dp                   :: freqintegral0,exchange_energy_core
      complex_dp                :: freqintegral0_pade

      if(all(job1%type/=[J_RPA,J_HFE])) return
      Rbegin
      if     (job1%type==J_RPA) then ; write(6,'(/A)') '=== Exchange-correlation energy ==='
      else if(job1%type==J_HFE) then ; write(6,'(/A)') '=== Exchange energy ==='
      endif
      if(job1%type==J_RPA) then
        write(6,'(/A)') '  Frequency  Integrand'
        write(6,'(3F18.10)') (freq(i),rpa_int(i),(freq(min(i+1,nfreq))-freq(max(i-1,1)))/2,i=1,nfreq)
        write(6,*)
        rdum = freqintegral0(rpa_int,freq,nfreq,3,4)
        cdum = freqintegral0_pade(rpa_int*(1d0,0d0),freq,nfreq) ! Pade integration is not very accurate for rpa_int; should be removed perhaps
        write(6,'(/A)') '  RPA correlation energy:'
        write(6,'( F12.7,12X,A)')  rdum,'  (spline)'
        write(6,'(2F12.7,A)')      cdum,'  (Pade)'
        deallocate ( rpa_int )
      endif
      rdum1 = 0
      rdum2 = exchange_energy_core()
      do i = 1,size(job1%band)
        rdum1 = rdum1 + selfx(i) * count(kptp==job1%kpt(i)) * wintgr(job1%kpt(i),job1%band(i),job1%spin(i)) / nspin2
      enddo
      write(6,'(/A,/F12.7,A,F12.7,A)') '  HF exchange energy:',rdum1,'   (',rdum2,' )'
      if(job1%type==J_RPA) then ; write(6,'(/A,/F12.7)') '  Total XC energy:',rdum1+rdum2+rdum
      else                      ; write(6,'(/A,/F12.7)') '  Total HF energy:',rdum1+rdum2
      endif
      write(6,'(/A/)') '==='
      Rend

      end

