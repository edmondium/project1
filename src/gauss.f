c--------------------------------------------------------------------------------
c Copyright (c) 2021 Peter Grünberg Institut, Forschungszentrum Jülich, Germany
c This file is part of SPEX and available as free software under the conditions
c of the MIT license as expressed in the LICENSE file in more detail.
c--------------------------------------------------------------------------------

c     Initializes BZ integration using Gaussian functions for occupied states,
c     i.e. calculates (for any f) the weights w(k,n,s) in the sum
c
c      2                                  2  occ               3
c     SUM SUM wintgr(k,n,s) * f(k,n,s) = SUM SUM INT f(k,n,s) d k .
c     s=1 k,n                            s=1  n
c
c     SUM wintgr(k,n,s)  gives the number of states up to the Fermi energy (efermi).
c     This can be used to determine the Fermi energy.
c
c     Gaussian function : exp(-(e-ene)**2/(2*gauss(1)**2)) / ( sqrt(2*pi) * gauss(1) )

# include "cppmacro.h"

      subroutine gauss_init(wintgr_out,dim,minb,maxb,efermi,spin,lshift)

      use global, only: ene,nkpt,gauss,nspin1

      use, intrinsic :: iso_fortran_env
      implicit none
      logical, intent(in)  :: lshift
      integer, intent(in)  :: minb,maxb,spin,dim
      real_dp, intent(in)  :: efermi
      real_dp, intent(out) :: wintgr_out(dim,minb:maxb,*)
      real_dp              :: sqrt2,rdum
      integer              :: ikpt,ikpt1,iband,ispin,spin1,spin2
      if(spin==0) then ; spin1 = 1    ; spin2 = nspin1
      else             ; spin1 = spin ; spin2 = spin
      endif
      sqrt2                          = sqrt(2d0)
      wintgr_out(:,:,:spin2-spin1+1) = 0
      do ispin = spin1,spin2
        do ikpt = 1,nkpt ; if(ikpt>dim) exit
          ikpt1 = ikpt ; if(lshift) ikpt1 = ikpt + nkpt
          do iband = minb,maxb
            rdum                                 = (efermi-ene(iband,ikpt1,ispin)) / (sqrt2 * gauss(1))
            wintgr_out(ikpt,iband,ispin-spin1+1) = (1+erf(rdum)) / 2
          enddo
        enddo
      enddo
      wintgr_out(:,:,:spin2-spin1+1) = wintgr_out(:,:,:spin2-spin1+1) / nkpt
      end

c ------------

c     Initializes Gaussian integration on a surface given by h(k)=h0 for occupied and unoccupied states,
c     i.e. calculates (for any f) the weights w(n',n,k) in the sum
c
c      2                                           2  occ unocc                          3
c     SUM  SUM  wintgr3(n',n,k,s) * f(k,n,n',s) = SUM SUM  SUM  INT        f(k,n,n',s) d k .
c     s=1 k,n,n'                                  s=1  n    n'     h(k)=h0
c
      subroutine gauss3_init(wintgr3,h0,ikpt,ispin1,ispin2,bandi,bandf)

      use global

      use, intrinsic :: iso_fortran_env
      implicit none
      integer, intent(in)  :: ikpt,ispin1,ispin2,bandi,bandf
      real_dp, intent(in)  :: h0
      real_dp, intent(out) :: wintgr3(bandi:bandf,bando,nkpt)
      real_dp              :: sqrt2,rdum1,rdum2,ene1,ene2,g2,g2fac
      integer              :: ikpt1,ikpt2,iband1,iband2
      integer              :: kptsum
      sqrt2   = sqrt(2d0)
      g2      = 2*gauss(2)**2
      g2fac   = 1d0 / (sqrt(2*pi)*gauss(2))
      wintgr3 = 0
      do ikpt1 = 1,nkpt
        if(ikpt<=nkpt) then
          ikpt2 = kptsum(ikpt1,ikpt)
        else
          ikpt2 = nkpt + ikpt1
        endif
        do iband1 = 1,bando
          ene1  = ene(iband1,ikpt1,ispin1)
          rdum1 = ( efermi - ene1  ) / ( sqrt2 * gauss(1) )
          rdum1 = ( 1 + erf(rdum1) ) / 2
          do iband2 = bandi,bandf
            ene2  = ene(iband2,ikpt2,ispin2)
            rdum2 = ( ene2  - efermi ) / ( sqrt2 * gauss(1) )
            rdum2 = ( 1 + erf(rdum2) ) / 2
            wintgr3(iband2,iband1,ikpt1) = rdum1 * rdum2 * g2fac * exp(- (ene2-ene1-h0)**2 / g2 )
          enddo
        enddo
      enddo
      wintgr3 = wintgr3 / nkpt
      end

c ------------

c     Initializes Gaussian integration on a surface given
c
c      s'           s
c     e  (q+k)  -  e (q)  =  w(i)     for occupied (n) and unoccupied states (n') and frequencies w(i),
c      n'           n
c
c     i.e. calculates (for any f) the weights w(n',n,k) in the sum
c
c      2                                           2  occ unocc                          3
c     SUM  SUM  wintgr4(n',n,k,s) * f(k,n,n',s) = SUM SUM  SUM  INT        f(k,n,n',s) d k .
c     s=1 k,n,n'                                  s=1  n    n'     h(k)=h0
c
c     freq(:nfreq) - frequencies w(i)
c     ikpt         - k index
c     kpt1(:nkpt1) - q indices (currently unused)
c     kpt1p(:nkpt) - parents of kpoints in kpt1 (can be 0 if parent is not in kpt1)
c     kflag        - 1 (normal), 2 (q and q+k switched)
c     ispin1/2     - s/s'
c     bandi/f|1/2  - n=(bandi1..bandf1), n'=(bandi2..bandf2)
c
c     See tetrahedron5_init for details about the frequency mesh.
c
      subroutine gauss4_init(wintgr4,freq,nfreq,ikpt,kpt1,kpt1p,nkpt1,kflag,ispin1,ispin2,bandi1,bandf1,bandi,bandf)

      use global, only: nkpt,ene,efermi,pi,gauss

      use, intrinsic :: iso_fortran_env
      implicit none
      integer, intent(in)  :: nfreq,ikpt,nkpt1,ispin1,ispin2,bandi1,bandf1,bandi,bandf,kpt1(nkpt1),kpt1p(nkpt),kflag
      real_dp, intent(in)  :: freq(nfreq)
      real_dp, intent(out) :: wintgr4(bandi:bandf,bandi1:bandf1,nkpt1,nfreq)
      real_dp              :: sqrt2,rdum,rdum1,rdum2,ene1,ene2,g2fac,gfac,f1,f2,df,dene,factor,a
      integer              :: ikpt1,ikpt2,ikptp,iband1,iband2,ifreq
      integer              :: kptsum

      ! Check whether freq(:) is in ascending order and calculate factor a (see above).
      rdum = 0
      do ifreq = 2,nfreq
        if(freq(ifreq)<=freq(ifreq-1)) Error('input frequency mesh not properly ordered.')
        if(ifreq/=nfreq) then
          rdum  = rdum  + (freq(ifreq)-freq(ifreq-1)) / (freq(ifreq+1)-freq(ifreq-1))
          rdum1 = rdum / (ifreq-1)
          if(ifreq>2.and.abs(rdum1-a)>1d-10) Bug('Hilbert mesh not exponential.')
          a     = rdum1
        endif
      enddo

      ! Define weights
      sqrt2   = sqrt(2d0)
      gfac    = 1d0 / (sqrt(2d0)*gauss(2))
      g2fac   = 1d0 / (sqrt(2*pi)*gauss(2))
      factor  = sqrt(0.5d0*pi) * gauss(2)
      wintgr4 = 0
      do ikpt1 = 1,nkpt
        ikptp = kpt1p(ikpt1) ; if(ikptp==0) cycle
        ikpt2 = kptsum(ikpt1,ikpt)
        do iband1 = bandi1,bandf1
          if(iband1==0) then
            ene1  = efermi-1d-8 ! special case: dispersion-less "core" state at efermi
            rdum1 = 1           ! "core" state fully occupied
          else
            if(kflag==1) then ; ene1 = ene(iband1,ikpt1,ispin1)
            else              ; ene1 = ene(iband1,ikpt2,ispin1)
            endif
            rdum1 = ( efermi - ene1  ) / ( sqrt2 * gauss(1) )
            rdum1 = ( 1 + erf(rdum1) ) / 2
          endif
          do iband2 = bandi,bandf
            if(kflag==1) then ; ene2 = ene(iband2,ikpt2,ispin2)
            else              ; ene2 = ene(iband2,ikpt1,ispin2)
            endif
            rdum2 = ( ene2 - efermi ) / ( sqrt2 * gauss(1) )
            rdum2 = ( 1 + erf(rdum2) ) / 2
            dene  = ene2 - ene1
            do ifreq = 1,nfreq
              if(freq(ifreq)==0) cycle
              if(ifreq/=nfreq) then ; df =    a  * ( freq(ifreq+1)- freq(ifreq)   )
              else                  ; df = (1-a) * ( freq(ifreq)  - freq(ifreq-1) )
              endif
              f1 = freq(ifreq) - df
              f2 = freq(ifreq) + df
              wintgr4(iband2,iband1,ikptp,ifreq) = wintgr4(iband2,iband1,ikptp,ifreq) +
     &                                             rdum1 * rdum2 * g2fac * factor * (
     &                                             erf ( gfac * (f2-dene) ) - erf ( gfac * (f1-dene) ) ) / (2*df)
            enddo
          enddo
        enddo
      enddo
      wintgr4 = wintgr4 / nkpt
      end

