c--------------------------------------------------------------------------------
c Copyright (c) 2021 Peter Grünberg Institut, Forschungszentrum Jülich, Germany
c This file is part of SPEX and available as free software under the conditions
c of the MIT license as expressed in the LICENSE file in more detail.
c--------------------------------------------------------------------------------

c     Solves the homogeneous scalar-relativistic Dirac equations
c
c     p' =  p/r + 2Mq
c     q' = -q/r +  wp
c
c     with
c       M  = 1 + ( e - V ) / (2c^2)
c       w  = l(l+1)/(2Mr^2) + V - e.
c
c     Uses 4th-order Runge-Kutta.
c
c     The sign of itype0 used as a special flag for an extended MT sphere.

# include "cppmacro.h"

      subroutine dirac_hom(p,q,dp,l,e,itype0,ispin)

      use global, only: grid,maxgrid,vmt,clight

      use, intrinsic :: iso_fortran_env
      implicit none
      integer, intent(in)  :: l,itype0,ispin
      real_dp, intent(in)  :: e
      real_dp, intent(out) :: p(maxgrid),q(maxgrid),dp
      real_dp              :: pp(maxgrid),qp(maxgrid)
      real_dp, pointer_cnt :: v0(:)
      real_dp              :: inc,incB,inci,inciB
      real_dp              :: ci,cci
      real_dp              :: rA,rB,rC,rAi,rBi,rCi
      real_dp              :: dr,m,w,v,ph,qh
      real_dp              :: k1,k2,k3,k4,l1,l2,l3,l4
      real_dp, parameter   :: sixi = 0.166666666666667d0
      integer              :: itype
      integer              :: i,ll,n
      real_dp              :: intgrf

      itype = abs(itype0)
      inc   = exp(grid(itype)%increment)
      
      if(itype0>0) then ! point v0=>vmt
        v0 => vmt(:,1,itype,ispin)
      else              ! extended grid: define v0 beyond original MT
        allocate(v0(maxgrid))        
        n        = size(vmt,1)
        v0(:n)   = vmt(:,1,itype,ispin)
        v0(n+1:) = 0
        rA       = grid(itype)%first / inc
        do n = 1,maxgrid
          if(v0(n)==0) exit
          rA = rA * inc
        enddo
        w  = v0(n-1) * rA
        v  = v0(n-2) * (rA/inc)
        v  = log(w/v) / ( rA-rA/inc )
        w  = w / exp(v*rA)
        rA = grid(itype)%first
        do i = 1,maxgrid
          if(v0(i)==0) v0(i) = w * exp ( v*rA ) / rA
          rA = rA * inc
        enddo
      endif

      rA    = grid(itype)%first
      n     = grid(itype)%number
      incB  = (1+inc)/2
      inci  = 1/inc
      inciB = 1/incB
      ll    = l*(l+1)
      ci    = 1/clight
      cci   = ci*ci

      ! Start values are derived from p=r**a and q=Q*r**a (w=a used temporarily in code):
      ! For r->0, we can set v0(r)=Z/r, so 2M=-Z/(r*c^2) and w=Z/r-l(l+1)*c^2/(Z*r).
      w     = sqrt(1+ll-(v0(1)*rA*ci)**2)
      m     = rA * ( 2 + ( e - v0(1) ) * cci )
      p(1)  = rA**(l+1)
      q(1)  = p(1) * ( w - 1 ) / m
      rAi   = 1/rA
      pp(1) = w * p(1) * rAi
      qp(1) = w * q(1) * rAi

      do i = 1,n-1
        rC      = inc  * rA
        rB      = incB * rA
        dr      = rC - rA

        rCi     = inci  * rAi
        rBi     = inciB * rAi

        k1      = pp(i)
        l1      = qp(i)

        v       = ( v0(i) + v0(i+1)*inc ) * inciB / 2
        m       = 2 + ( e - v ) * cci
        w       = ll/(m*rB**2) + v - e
        ph      = p(i) + k1*dr/2
        qh      = q(i) + l1*dr/2
        k2      =   ph*rBi + m*qh
        l2      = - qh*rBi + w*ph

        ph      = p(i) + k2*dr/2
        qh      = q(i) + l2*dr/2
        k3      =   ph*rBi + m*qh
        l3      = - qh*rBi + w*ph

        v       = v0(i+1)
        m       = 2 + ( e - v ) * cci
        w       = ll/(m*rC**2) + v - e
        ph      = p(i) + k3*dr
        qh      = q(i) + l3*dr
        k4      =   ph*rCi + m*qh
        l4      = - qh*rCi + w*ph

        p(i+1)  = p(i) + dr * sixi * (k1+k2+k2+k3+k3+k4)
        q(i+1)  = q(i) + dr * sixi * (l1+l2+l2+l3+l3+l4)
        pp(i+1) =   p(i+1)*rCi + m*q(i+1)
        qp(i+1) = - q(i+1)*rCi + w*p(i+1)

        rA      = rC
        rAi     = rCi
      enddo

      q  = q * ci
      if(itype0>0) then ; w = 1/sqrt(intgrf(p**2+q**2,itype))
      else              ; w = 1
      endif
      p  = p * w
      q  = q * w
      dp = ( pp(n)*w - p(n)/rA ) / rA

      if(itype0>0) then
        nullify(v0)
      else
        deallocate(v0)
      endif

      end

c ---------------------------

c     Same for complex energies e

      subroutine cdirac_hom(p,q,dp,l,e,itype,ispin)

      use global, only: grid,maxgrid,vmt,clight

      use, intrinsic :: iso_fortran_env
      implicit none
      integer, intent(in)     :: l,itype,ispin
      complex_dp, intent(in)  :: e
      complex_dp, intent(out) :: p(maxgrid),q(maxgrid),dp
      complex_dp              :: pp(maxgrid),qp(maxgrid)
      complex_dp              :: m,w,ph,qh
      complex_dp              :: k1,k2,k3,k4,l1,l2,l3,l4
      real_dp,    pointer_cnt :: v0(:)
      real_dp                 :: inc,incB,inci,inciB
      real_dp                 :: ci,cci
      real_dp                 :: rA,rB,rC,rAi,rBi,rCi,dr,v
      real_dp, parameter      :: sixi = 0.166666666666667d0
      integer                 :: i,ll,n
      real_dp                 :: intgrf

      rA    = grid(itype)%first
      n     = grid(itype)%number
      inc   = exp(grid(itype)%increment)
      incB  = (1+inc)/2
      inci  = 1/inc
      inciB = 1/incB
      ll    = l*(l+1)
      ci    = 1/clight
      cci   = ci*ci
      v0    =>vmt(:,1,itype,ispin)

      w     = sqrt(1+ll-(v0(1)*rA*ci)**2)
      m     = rA * ( 2 + ( e - v0(1) ) * cci )
      p(1)  = rA**(l+1)
      q(1)  = p(1) * ( w - 1 ) / m
      rAi   = 1/rA
      pp(1) = w * p(1) * rAi
      qp(1) = w * q(1) * rAi

      do i = 1,n-1
        rC      = inc  * rA
        rB      = incB * rA
        dr      = rC - rA

        rCi     = inci  * rAi
        rBi     = inciB * rAi

        k1      = pp(i)
        l1      = qp(i)

        v       = ( v0(i) + v0(i+1)*inc ) * inciB / 2
        m       = 2 + ( e - v ) * cci
        w       = ll/(m*rB**2) + v - e
        ph      = p(i) + k1*dr/2
        qh      = q(i) + l1*dr/2
        k2      =   ph*rBi + m*qh
        l2      = - qh*rBi + w*ph

        ph      = p(i) + k2*dr/2
        qh      = q(i) + l2*dr/2
        k3      =   ph*rBi + m*qh
        l3      = - qh*rBi + w*ph

        v       = v0(i+1)
        m       = 2 + ( e - v ) * cci
        w       = ll/(m*rC**2) + v - e
        ph      = p(i) + k3*dr
        qh      = q(i) + l3*dr
        k4      =   ph*rCi + m*qh
        l4      = - qh*rCi + w*ph

        p(i+1)  = p(i) + dr * sixi * (k1+k2+k2+k3+k3+k4)
        q(i+1)  = q(i) + dr * sixi * (l1+l2+l2+l3+l3+l4)
        pp(i+1) =   p(i+1)*rCi + m*q(i+1)
        qp(i+1) = - q(i+1)*rCi + w*p(i+1)

        rA      = rC
        rAi     = rCi
      enddo

      q  = q * ci
      w  = 1/sqrt(intgrf(real(p)**2+imag(p)**2+real(q)**2+imag(q)**2,itype))
      p  = p * w
      q  = q * w
      dp = ( pp(n)*w - p(n)/rA ) / rA

      nullify(v0)

      end

c ---------------------------

c     Solves the inhomogeneous scalar-relativistic Dirac equations
c
c     p' =  p/r + 2Mq + c 2M'qi
c     q' = -q/r +  wp +    w'pi
c
c     with
c       M  = 1 + ( e - V ) / (2c^2)
c       w  = l(l+1)/(2Mr^2) + V - e
c       M' = dM/de = 1/(2*c^2)
c       w' = dw/de = -1 - l(l+1)/(2cMr)^2.
c
c     Uses 4th-order Runga-Kutta.

      subroutine dirac_inhom(p,q,dp,l,e,pi,qi,itype,ispin)

      use global, only: grid,maxgrid,vmt,clight

      use, intrinsic :: iso_fortran_env
      implicit none
      integer, intent(in)  :: l,itype,ispin
      real_dp, intent(out) :: p(maxgrid),q(maxgrid),dp
      real_dp, intent(in)  :: pi(maxgrid),qi(maxgrid),e
      real_dp              :: pp(maxgrid),qp(maxgrid)
      real_dp, pointer_cnt :: v0(:)
      real_dp              :: inc,incB,inci,inciB
      real_dp              :: ci,cci
      real_dp              :: rA,rB,rC,rAi,rBi,rCi
      real_dp              :: dr,m,w,wh,v,ph,qh,inh1,inh2
      real_dp              :: k1,k2,k3,k4,l1,l2,l3,l4
      real_dp, parameter   :: sixi = 0.166666666666667d0
      integer              :: i,ll,n

      rA    = grid(itype)%first
      n     = grid(itype)%number
      inc   = exp(grid(itype)%increment)
      incB  = (1+inc)/2
      inci  = 1/inc
      inciB = 1/incB
      ll    = l*(l+1)
      ci    = 1/clight
      cci   = ci*ci
      v0    =>vmt(:,1,itype,ispin)

      m     = 2 + ( e - v0(1) ) * cci
      p(1)  = 0
      q(1)  = 0
      pp(1) =  qi(1) * ci
      qp(1) = -pi(1) * ( 1 + ll*cci/(m*rA)**2 )

      rAi   = 1/rA

      do i = 1,n-1
        rC      = inc  * rA
        rB      = incB * rA
        dr      = rC - rA

        rCi     = inci  * rAi
        rBi     = inciB * rAi

        k1      = pp(i)
        l1      = qp(i)

        v       = ( v0(i) + v0(i+1)*inc ) * inciB / 2
        m       = 2 + ( e - v ) * cci
        wh      = ll/(m*rB**2)
        w       = wh + v - e
        ph      = p(i) + k1*dr/2  ; inh1 =   (qi(i)+qi(i+1))/2 * ci
        qh      = q(i) + l1*dr/2  ; inh2 = - (pi(i)+pi(i+1))/2 * ( 1 + wh*cci/m )
        k2      =   ph*rBi + m*qh + inh1
        l2      = - qh*rBi + w*ph + inh2

        ph      = p(i) + k2*dr/2
        qh      = q(i) + l2*dr/2
        k3      =   ph*rBi + m*qh + inh1
        l3      = - qh*rBi + w*ph + inh2

        v       = v0(i+1)
        m       = 2 + ( e - v ) * cci
        wh      = ll/(m*rC**2)
        w       = wh + v - e
        ph      = p(i) + k3*dr    ; inh1 =   qi(i+1) * ci
        qh      = q(i) + l3*dr    ; inh2 = - pi(i+1) * ( 1 + wh*cci/m )
        k4      =   ph*rCi + m*qh + inh1
        l4      = - qh*rCi + w*ph + inh2

        p(i+1)  = p(i) + dr * sixi * (k1+k2+k2+k3+k3+k4)
        q(i+1)  = q(i) + dr * sixi * (l1+l2+l2+l3+l3+l4)
        pp(i+1) =   p(i+1)*rCi + m*q(i+1) + inh1
        qp(i+1) = - q(i+1)*rCi + w*p(i+1) + inh2

        rA      = rC
        rAi     = rCi
      enddo

      q  = q * ci
      dp = ( pp(n) - p(n)/rA ) / rA

      nullify(v0)

      end

c ---------------------------

c     Same for complex energies e

      subroutine cdirac_inhom(p,q,dp,l,e,pi,qi,itype,ispin)

      use global, only: grid,maxgrid,vmt,clight

      use, intrinsic :: iso_fortran_env
      implicit none
      integer,    intent(in)  :: l,itype,ispin
      complex_dp, intent(in)  :: e
      complex_dp, intent(out) :: p(maxgrid),q(maxgrid),dp
      complex_dp, intent(in)  :: pi(maxgrid),qi(maxgrid)
      complex_dp              :: pp(maxgrid),qp(maxgrid)
      complex_dp              :: m,w,wh,ph,qh,inh1,inh2
      complex_dp              :: k1,k2,k3,k4,l1,l2,l3,l4
      real_dp,    pointer_cnt :: v0(:)
      real_dp                 :: inc,incB,inci,inciB
      real_dp                 :: ci,cci
      real_dp                 :: rA,rB,rC,rAi,rBi,rCi,dr,v
      real_dp,    parameter   :: sixi = 0.166666666666667d0
      integer                 :: i,ll,n

      rA    = grid(itype)%first
      n     = grid(itype)%number
      inc   = exp(grid(itype)%increment)
      incB  = (1+inc)/2
      inci  = 1/inc
      inciB = 1/incB
      ll    = l*(l+1)
      ci    = 1/clight
      cci   = ci*ci
      v0    =>vmt(:,1,itype,ispin)

      m     = 2 + ( e - v0(1) ) * cci
      p(1)  = 0
      q(1)  = 0
      pp(1) =  qi(1) * ci
      qp(1) = -pi(1) * ( 1 + ll*cci/(m*rA)**2 )

      rAi   = 1/rA

      do i = 1,n-1
        rC      = inc  * rA
        rB      = incB * rA
        dr      = rC - rA

        rCi     = inci  * rAi
        rBi     = inciB * rAi

        k1      = pp(i)
        l1      = qp(i)

        v       = ( v0(i) + v0(i+1)*inc ) * inciB / 2
        m       = 2 + ( e - v ) * cci
        wh      = ll/(m*rB**2)
        w       = wh + v - e
        ph      = p(i) + k1*dr/2  ; inh1 =   (qi(i)+qi(i+1))/2 * ci
        qh      = q(i) + l1*dr/2  ; inh2 = - (pi(i)+pi(i+1))/2 * ( 1 + wh*cci/m )
        k2      =   ph*rBi + m*qh + inh1
        l2      = - qh*rBi + w*ph + inh2

        ph      = p(i) + k2*dr/2
        qh      = q(i) + l2*dr/2
        k3      =   ph*rBi + m*qh + inh1
        l3      = - qh*rBi + w*ph + inh2

        v       = v0(i+1)
        m       = 2 + ( e - v ) * cci
        wh      = ll/(m*rC**2)
        w       = wh + v - e
        ph      = p(i) + k3*dr    ; inh1 =   qi(i+1) * ci
        qh      = q(i) + l3*dr    ; inh2 = - pi(i+1) * ( 1 + wh*cci/m )
        k4      =   ph*rCi + m*qh + inh1
        l4      = - qh*rCi + w*ph + inh2

        p(i+1)  = p(i) + dr * sixi * (k1+k2+k2+k3+k3+k4)
        q(i+1)  = q(i) + dr * sixi * (l1+l2+l2+l3+l3+l4)
        pp(i+1) =   p(i+1)*rCi + m*q(i+1) + inh1
        qp(i+1) = - q(i+1)*rCi + w*p(i+1) + inh2

        rA      = rC
        rAi     = rCi
      enddo

      q  = q * ci
      dp = ( pp(n) - p(n)/rA ) / rA

      nullify(v0)

      end

c ---------------------------

c     CORE STATES

c     Outward integration up to the classical turning point (rt).
c     Returns the homogeneous solution in pout/qout (not normalized) and the boundary values at rt in pt/qt.
c
c     rt and nt are returned if nt=0.
c
      subroutine core_dirac_hom_out(pout,qout,pt,qt,rt,nt,l,e,itype,ispin)

      use global, only: grid,maxgrid,vmt,clight

      use, intrinsic :: iso_fortran_env
      implicit none
      integer,    intent(in)  :: l,itype,ispin
      integer,    intent(in)  :: nt
      real_dp,    intent(in)  :: rt
      complex_dp, intent(in)  :: e
      complex_dp, intent(out) :: pout(maxgrid),qout(maxgrid),pt,qt
      complex_dp              :: m,w,ph,qh
      complex_dp              :: k1,k2,k3,k4,l1,l2,l3,l4
      complex_dp              :: pp,qp
      complex_dp, allocatable :: p(:),q(:)
      real_dp,    allocatable :: v0(:)
      real_dp                 :: inc,incB,inci,inciB
      real_dp                 :: ci,cci
      real_dp                 :: r,rA,rB,rC,rAi,rBi,rCi,dr,v,va,vb
      real_dp, parameter      :: sixi = 0.166666666666667d0
      integer                 :: i,ll,n

      if(nt==0) Error('turning point (rt,nt) undefined.')

      ! prepare potential continuation ( va * exp(vb*r)/r )
      inc = exp(grid(itype)%increment)
      rA  = grid(itype)%radius
      rB  = rA / inc
      n   = grid(itype)%number
      va  = vmt(n,  1,itype,ispin) * rA
      vb  = vmt(n-1,1,itype,ispin) * rB
      vb  = log(va/vb) / ( rA-rB )
      va  = va / exp(vb*rA)

      ! define potential up to rt
      ll = l*(l+1)
      n  = nt
      allocate ( v0(n),p(n),q(n) )
      r = grid(itype)%first
      do i = 1,n
        if(i<=grid(itype)%number) then ; v0(i) = vmt(i,1,itype,ispin)
        else                           ; v0(i) = va * exp ( vb*r ) / r
        endif
        r = r * inc
      enddo

      rA    = grid(itype)%first
      incB  = (1+inc)/2
      inci  = 1/inc
      inciB = 1/incB
      ci    = 1/clight
      cci   = ci*ci

      ! start values (p=r**w, q=b*r**w)
      w     = sqrt(1+ll-(v0(1)*rA*ci)**2)
      m     = rA * ( 2 + ( e - v0(1) ) * cci )
      p(1)  = rA**(l+1)
      q(1)  = p(1) * ( w - 1 ) / m
      rAi   = 1/rA
      pp    = w * p(1) * rAi
      qp    = w * q(1) * rAi

      ! Runge-Kutta
      do i = 1,n-1

        if(abs(real(p(i)))+abs(imag(p(i)))+abs(real(q(i)))+abs(imag(q(i)))>1d2) then
          p(:i) = p(:i) / 1d4
          q(:i) = q(:i) / 1d4
          pp    = pp    / 1d4
          qp    = qp    / 1d4
        endif

        rC      = inc  * rA
        rB      = incB * rA
        dr      = rC - rA

        rCi     = inci  * rAi
        rBi     = inciB * rAi

        k1      = pp
        l1      = qp

        v       = ( v0(i) + v0(i+1)*inc ) * inciB / 2
        m       = 2 + ( e - v ) * cci
        w       = ll/(m*rB**2) + v - e
        ph      = p(i) + k1*dr/2
        qh      = q(i) + l1*dr/2
        k2      =   ph*rBi + m*qh
        l2      = - qh*rBi + w*ph

        ph      = p(i) + k2*dr/2
        qh      = q(i) + l2*dr/2
        k3      =   ph*rBi + m*qh
        l3      = - qh*rBi + w*ph

        v       = v0(i+1)
        m       = 2 + ( e - v ) * cci
        w       = ll/(m*rC**2) + v - e
        ph      = p(i) + k3*dr
        qh      = q(i) + l3*dr
        k4      =   ph*rCi + m*qh
        l4      = - qh*rCi + w*ph

        p(i+1)  = p(i) + dr * sixi * (k1+k2+k2+k3+k3+k4)
        q(i+1)  = q(i) + dr * sixi * (l1+l2+l2+l3+l3+l4)
        pp      =   p(i+1)*rCi + m*q(i+1)
        qp      = - q(i+1)*rCi + w*p(i+1)

        rA      = rC
        rAi     = rCi
      enddo

      q        = q * ci
      pt       = p(n)
      qt       = q(n)
      n        = min(maxgrid,n)
      pout(:n) = p(:n)
      qout(:n) = q(:n)
      deallocate ( v0,p,q )

      end

c ---------------------------

c     Inward integration down to rc.

      subroutine core_dirac_hom_in(pout,qout,pt,qt,rt,nt,l,e,itype,ispin)

      use global, only: grid,maxgrid,vmt,clight

      use, intrinsic :: iso_fortran_env
      implicit none
      integer,    intent(in)  :: l,itype,ispin
      integer,    intent(in)  :: nt
      real_dp,    intent(in)  :: rt
      complex_dp, intent(in)  :: e
      complex_dp, intent(out) :: pout(maxgrid),qout(maxgrid),pt,qt
      complex_dp              :: m,w,ph,qh
      complex_dp              :: k1,k2,k3,k4,l1,l2,l3,l4
      complex_dp              :: pp,qp
      complex_dp, allocatable :: p(:),q(:)
      real_dp,    allocatable :: v0(:)
      real_dp                 :: inc,incB,inci,inciB
      real_dp                 :: ci,cci
      real_dp                 :: r,rA,rB,rC,rAi,rBi,rCi,dr,v,va,vb
      real_dp, parameter      :: sixi = 0.166666666666667d0
      integer                 :: i,j,ll,n

      if(nt==0) Error('turning point (rt,nt) undefined.')

      ! prepare potential continuation ( va * exp(vb*r)/r )
      inc = exp(grid(itype)%increment)
      rA  = grid(itype)%radius
      rB  = rA / inc
      n   = grid(itype)%number
      va  = vmt(n,  1,itype,ispin) * rA
      vb  = vmt(n-1,1,itype,ispin) * rB
      vb  = log(va/vb) / ( rA-rB )
      va  = va / exp(vb*rA)

      ! extend mesh to 15*rt
      n = nt
      r = rt
      do while(r<15*rt)!grid(itype)%radius+20)
        n = n + 1
        r = r * inc
      enddo

      ! define potential up to 15*rt
      allocate ( v0(n),p(n),q(n) )
      r = grid(itype)%first
      do i = 1,n
        if(i<=grid(itype)%number) then ; v0(i) = vmt(i,1,itype,ispin)
        else                           ; v0(i) = va * exp ( vb*r ) / r
        endif
        r = r * inc
      enddo
      r = r / inc

      ll    = l*(l+1)
      incB  = (1+inc)/2
      inci  = 1/inc
      inciB = 1/incB
      ci    = 1/clight
      cci   = ci*ci
      rA    = r
      rAi   = 1/rA

      ! start values ( q/p=-sqrt(-e/m) )
      m     = 2 + e*cci
      w     = ll/(m*r**2) - e
      p(n)  = 1d-8
      q(n)  = 1d-8 * (-sqrt(-e/m))
      pp    =  p(n)/r + m * q(n)
      qp    = -q(n)/r + w * p(n)

c      rewind(1000)

      ! Runge-Kutta
      do i = n,nt+1,-1 !; write(1000,*) i,real(p(i))

        if(abs(real(p(i)))+abs(imag(p(i)))>1d2) then
          j      = max(i,maxgrid)
          p(i:j) = p(i:j) / 1d4
          q(i:j) = q(i:j) / 1d4
          pp     = pp     / 1d4
          qp     = qp     / 1d4
        endif

        rC      = inci  * rA
        rB      = inciB * rA
        dr      = rC - rA

        rCi     = inc  * rAi
        rBi     = incB * rAi

        k1      = pp
        l1      = qp

        v       = ( v0(i) + v0(i-1)*inci ) * incB / 2 !; write(*,*) v0(i),v0(i-1),v
        m       = 2 + ( e - v ) * cci
        w       = ll/(m*rB**2) + v - e
        ph      = p(i) + k1*dr/2
        qh      = q(i) + l1*dr/2
        k2      =   ph*rBi + m*qh
        l2      = - qh*rBi + w*ph

        ph      = p(i) + k2*dr/2
        qh      = q(i) + l2*dr/2
        k3      =   ph*rBi + m*qh
        l3      = - qh*rBi + w*ph

        v       = v0(i-1)
        m       = 2 + ( e - v ) * cci
        w       = ll/(m*rC**2) + v - e
        ph      = p(i) + k3*dr
        qh      = q(i) + l3*dr
        k4      =   ph*rCi + m*qh
        l4      = - qh*rCi + w*ph

        p(i-1)  = p(i) + dr * sixi * (k1+k2+k2+k3+k3+k4) !; write(*,*) rA,p(i),p(i-1),k1,k2,k3,k4 ; read(*,*)
        q(i-1)  = q(i) + dr * sixi * (l1+l2+l2+l3+l3+l4)
        pp      =   p(i-1)*rCi + m*q(i-1)
        qp      = - q(i-1)*rCi + w*p(i-1) !; write(*,*) i-1,p(i-1),q(i-1);read(*,*)

        rA      = rC
        rAi     = rCi
      enddo

      q  = q * ci
      pt = p(nt)
      qt = q(nt)
      if(nt<=maxgrid) then
        i          = min(n,maxgrid)
        pout(nt:i) = p(nt:i)
        qout(nt:i) = q(nt:i)
        if(n<maxgrid) then
          pout(n+1:) = 0
          qout(n+1:) = 0
        endif
      endif

      deallocate ( v0,p,q )

      end

c ---------------------------

c     Outward integration up to rt.
c     Returns the inhomogeneous solution in pout/qout and pt/qt.

      subroutine core_dirac_inhom_out(pout,qout,pt,qt,rt,nt,alpha,l,e,pc,qc,bm,itype,ispin)

      use global, only: grid,maxgrid,vmt,clight

      use, intrinsic :: iso_fortran_env
      implicit none
      integer,    intent(in)  :: l,itype,ispin
      integer,    intent(in)  :: nt,alpha
      real_dp,    intent(in)  :: rt
      real_dp,    intent(in)  :: pc(maxgrid),qc(maxgrid),bm(maxgrid)
      complex_dp, intent(in)  :: e
      complex_dp, intent(out) :: pout(maxgrid),qout(maxgrid),pt,qt
      complex_dp              :: m,w,wh,ph,qh,inh1,inh2
      complex_dp              :: k1,k2,k3,k4,l1,l2,l3,l4
      complex_dp              :: k01,k02,k03,k04,l01,l02,l03,l04
      complex_dp              :: pp,qp,pp0,qp0
      complex_dp, allocatable :: p(:),q(:),p0(:),q0(:)
      real_dp,    allocatable :: v0(:),pi(:),qi(:)
      real_dp                 :: inc,incB,inci,inciB
      real_dp                 :: ci,cci
      real_dp                 :: r,rA,rB,rC,rAi,rBi,rCi,dr,v,va,vb,ca1,cb1,ca2,cb2,bm0,bm1
      real_dp, parameter      :: sixi = 0.166666666666667d0
      integer                 :: i,ll,n

      if(nt==0) Error('turning point (rt,nt) undefined.')

      ! prepare potential continuation ( va * exp(vb*r)/r )
      inc = exp(grid(itype)%increment)
      rA  = grid(itype)%radius
      rB  = rA / inc
      n   = grid(itype)%number
      va  = vmt(n,  1,itype,ispin) * rA
      vb  = vmt(n-1,1,itype,ispin) * rB
      vb  = log(va/vb) / ( rA-rB )
      va  = va / exp(vb*rA)

      ! continuation of bm(:)
      bm0 =   bm(n) / rA
      bm1 = ( bm(n) - bm(n-1) ) / (rA-rB)

      ! prepare core state (pc/qc) continuation ( r**alpha * ca1 * exp(cb1*r), ca2,cb2 for qc )
      ca1 = pc(n)   / rA**(alpha+1)
      cb1 = pc(n-1) / rB**(alpha+1)
      cb1 = log(ca1/cb1) / ( rA-rB )
      ca1 = ca1 / exp(cb1*rA)
      ca2 = qc(n)   / rA**(alpha+1)
      cb2 = qc(n-1) / rB**(alpha+1)
      cb2 = log(ca2/cb2) / ( rA-rB )
      ca2 = ca2 / exp(cb2*rA)

      ! define potential and pi/qi up to rt
      n = nt ; allocate ( v0(n),p(n),q(n),p0(n),q0(n),pi(n),qi(n) )
      r = grid(itype)%first
      do i = 1,n
        if(i<=grid(itype)%number) then
          v0(i) = vmt(i,1,itype,ispin)
          pi(i) = pc(i) * bm(i)/r
          qi(i) = qc(i) * bm(i)/r
        else
          v0(i) = va * exp ( vb*r ) / r
          pi(i) = r**(alpha+1) * ca1 * exp(cb1*r) * ( bm0 + bm1*(r-rA)/r )
          qi(i) = r**(alpha+1) * ca2 * exp(cb2*r) * ( bm0 + bm1*(r-rA)/r )
        endif
        r = r * inc
      enddo

      ll    = l*(l+1)
      rA    = grid(itype)%first
      incB  = (1+inc)/2
      inci  = 1/inc
      inciB = 1/incB
      ci    = 1/clight
      cci   = ci*ci

      ! start values (p=r**w, q=b*r**w)
      w     = sqrt(1+ll-(v0(1)*rA*ci)**2)
      m     = 2 + ( e - v0(1) ) * cci
      p0(1) = rA**(l+1)
      q0(1) = p0(1) * ( w - 1 ) / (m*rA)
      pp0   = w * p(1) / rA
      qp0   = w * q(1) / rA
      p(1)  = 0
      q(1)  = 0
      pp    =  qi(1) * ci
      qp    = -pi(1) * ( 1 + ll*cci/(m*rA)**2 )
      rAi   = 1/rA

      ! Runge-Kutta
      do i = 1,n-1

        if(abs(real(p0(i)))+abs(imag(p0(i)))+abs(real(q0(i)))+abs(imag(q0(i)))>1d2) then
          p0(:i) = p0(:i) / 1d4
          q0(:i) = q0(:i) / 1d4
          pp0    = pp0    / 1d4
          qp0    = qp0    / 1d4
        endif
        if(abs(real(p(i)))+abs(imag(p(i)))+abs(real(q(i)))+abs(imag(q(i)))>1d2) then
          w      = ( conjg(p0(i))*p(i)  + conjg(q0(i))*q(i)  )  /
     &             ( conjg(p0(i))*p0(i) + conjg(q0(i))*q0(i) )
          p(:i)  = p(:i) - w * p0(:i)
          q(:i)  = q(:i) - w * q0(:i)
          pp     = pp    - w * pp0
          qp     = qp    - w * qp0
        endif

        rC      = inc  * rA
        rB      = incB * rA
        dr      = rC - rA

        rCi     = inci  * rAi
        rBi     = inciB * rAi

        k1      = pp
        l1      = qp
        k01     = pp0
        l01     = qp0

        v       = ( v0(i) + v0(i+1)*inc ) * inciB / 2
        m       = 2 + ( e - v ) * cci
        wh      = ll/(m*rB**2)
        w       = wh + v - e
        ph      = p(i) + k1*dr/2  ; inh1 =   (qi(i)+qi(i+1))/2 * ci
        qh      = q(i) + l1*dr/2  ; inh2 = - (pi(i)+pi(i+1))/2 * ( 1 + wh*cci/m )
        k2      =   ph*rBi + m*qh + inh1
        l2      = - qh*rBi + w*ph + inh2
        ph      = p0(i) + k01*dr/2
        qh      = q0(i) + l01*dr/2
        k02     =   ph*rBi + m*qh
        l02     = - qh*rBi + w*ph

        ph      = p(i) + k2*dr/2
        qh      = q(i) + l2*dr/2
        k3      =   ph*rBi + m*qh + inh1
        l3      = - qh*rBi + w*ph + inh2
        ph      = p0(i) + k02*dr/2
        qh      = q0(i) + l02*dr/2
        k03     =   ph*rBi + m*qh
        l03     = - qh*rBi + w*ph

        v       = v0(i+1)
        m       = 2 + ( e - v ) * cci
        wh      = ll/(m*rC**2)
        w       = wh + v - e
        ph      = p(i) + k3*dr    ; inh1 =   qi(i+1) * ci
        qh      = q(i) + l3*dr    ; inh2 = - pi(i+1) * ( 1 + wh*cci/m )
        k4      =   ph*rCi + m*qh + inh1
        l4      = - qh*rCi + w*ph + inh2
        ph      = p0(i) + k03*dr
        qh      = q0(i) + l03*dr
        k04     =   ph*rCi + m*qh
        l04     = - qh*rCi + w*ph

        p(i+1)  = p(i) + dr * sixi * (k1+k2+k2+k3+k3+k4)
        q(i+1)  = q(i) + dr * sixi * (l1+l2+l2+l3+l3+l4)
        pp      =   p(i+1)*rCi + m*q(i+1) + inh1
        qp      = - q(i+1)*rCi + w*p(i+1) + inh2
        p0(i+1) = p0(i) + dr * sixi * (k01+k02+k02+k03+k03+k04)
        q0(i+1) = q0(i) + dr * sixi * (l01+l02+l02+l03+l03+l04)
        pp0     =   p0(i+1)*rCi + m*q0(i+1)
        qp0     = - q0(i+1)*rCi + w*p0(i+1)

        rA      = rC
        rAi     = rCi
      enddo

      q        = q * ci
      pt       = p(n)
      qt       = q(n)
      n        = min(grid(itype)%number,n)
      pout(:n) = p(:n)
      qout(:n) = q(:n)
      deallocate ( v0,p,q,p0,q0,pi,qi )

      end

c ---------------------------

c     Inward integration down to rc.

      subroutine core_dirac_inhom_in(pout,qout,pt,qt,rt,nt,alpha,l,e,pc,qc,bm,itype,ispin)

      use global, only: grid,maxgrid,vmt,clight

      use, intrinsic :: iso_fortran_env
      implicit none
      integer,    intent(in)  :: l,itype,ispin
      integer,    intent(in)  :: nt,alpha
      real_dp,    intent(in)  :: rt
      real_dp,    intent(in)  :: pc(maxgrid),qc(maxgrid),bm(maxgrid)
      complex_dp, intent(in)  :: e
      complex_dp, intent(out) :: pout(maxgrid),qout(maxgrid),pt,qt
      complex_dp              :: m,w,wh,ph,qh,inh1,inh2
      complex_dp              :: k1,k2,k3,k4,l1,l2,l3,l4
      complex_dp              :: k01,k02,k03,k04,l01,l02,l03,l04
      complex_dp              :: pp,qp,pp0,qp0
      complex_dp, allocatable :: p(:),q(:),pi(:),qi(:),p0(:),q0(:)
      real_dp,    allocatable :: v0(:)
      real_dp                 :: inc,incB,inci,inciB
      real_dp                 :: ci,cci
      real_dp                 :: r,rA,rB,rC,rAi,rBi,rCi,dr,v,va,vb,ca1,cb1,ca2,cb2,bm0,bm1 !,rAA
      real_dp, parameter      :: sixi = 0.166666666666667d0
      integer                 :: i,j,ll,n

      if(nt==0) Error('turning point (rt,nt) undefined.')

      ! prepare potential continuation ( va * exp(vb*r)/r )
      inc = exp(grid(itype)%increment)
      rA  = grid(itype)%radius
      rB  = rA / inc
      n   = grid(itype)%number
      va  = vmt(n,  1,itype,ispin) * rA
      vb  = vmt(n-1,1,itype,ispin) * rB

c      vb = (va-vb) / (rA-rB)
c      va = va - vb*rA
      vb  = log(va/vb) / ( rA-rB )
      va  = va / exp(vb*rA)

      ! continuation of bm(:)
      bm0 =   bm(n) / rA
      bm1 = ( bm(n) - bm(n-1) ) / (rA-rB)

      ! prepare core state (pc/qc) continuation ( r**alpha * ca1 * exp(cb1*r), ca2,cb2 for qc )
      ca1 = pc(n)   / rA**(alpha+1)
      cb1 = pc(n-1) / rB**(alpha+1)
      cb1 = log(ca1/cb1) / ( rA-rB )
      ca1 = ca1 / exp(cb1*rA)
      ca2 = qc(n)   / rA**(alpha+1)
      cb2 = qc(n-1) / rB**(alpha+1)
      cb2 = log(ca2/cb2) / ( rA-rB )
      ca2 = ca2 / exp(cb2*rA)

      ! extend mesh to 15*rt
      n = nt
      r = rt
      do while(r<15*rt)!grid(itype)%radius+10)
        n = n + 1
        r = r * inc
      enddo

      ! define potential and pi/qi up to 15*rt
      allocate ( v0(n),p(n),q(n),pi(n),qi(n),p0(n),q0(n) )
      bm0 = bm(grid(itype)%number) / grid(itype)%radius
      bm1 = ( bm(grid(itype)%number) - bm(grid(itype)%number-1) ) / (rA-rB)
      r   = grid(itype)%first
c      write(1004,*) '#',bm0,bm1
      do i = 1,n
        if(i<=grid(itype)%number) then
          v0(i) = vmt(i,1,itype,ispin)
          pi(i) = pc(i) * bm(i)/r
          qi(i) = qc(i) * bm(i)/r
c          write(1004,'(6ES40.20)') r,pi(i),qi(i),bm(i)
        else
c          v0(i) = va*r + vb
          v0(i) = va * exp ( vb*r ) / r
          pi(i) = r**(alpha+1) * ca1 * exp(cb1*r) * ( bm0 + bm1*(r-rA)/r )
          qi(i) = r**(alpha+1) * ca2 * exp(cb2*r) * ( bm0 + bm1*(r-rA)/r )
c          write(1004,'(5ES40.20)') r,pi(i),qi(i)
        endif
        r = r * inc
      enddo
      r = r / inc

      ll    = l*(l+1)
      incB  = (1+inc)/2
      inci  = 1/inc
      inciB = 1/incB
      ci    = 1/clight
      cci   = ci*ci
      rA    = r
      rAi   = 1/rA

      ! start values ( q/p=-sqrt(-e/m) )
      m     = 2 + e*cci
      w     = ll/(m*r**2) - e
      p0(n) = 1d-8
      q0(n) = 1d-8 * (-sqrt(-e/m))
      pp0   =  p0(n)/r + m * q0(n)
      qp0   = -q0(n)/r + w * p0(n)
      p(n)  = p0(n)
      q(n)  = q0(n)
      pp    = pp0 + qi(n) * ci
      qp    = qp0 - pi(n) * ( 1 + ll*cci/(m*rA)**2 )

c      rAA = r
c      rewind(1001)

      ! Runge-Kutta
      do i = n,nt+1,-1 !; write(1001,*) rA,real(p0(i)),real(p(i)) !; write(1001,*) i,abs(p(i)),abs(pi(i)),abs(v0(i))

        if(abs(real(p0(i)))+abs(imag(p0(i)))+abs(real(q0(i)))+abs(imag(q0(i)))>1d2) then
          p0(i:) = p0(i:) / 1d4
          q0(i:) = q0(i:) / 1d4
          pp0    = pp0    / 1d4
          qp0    = qp0    / 1d4
        endif
        if(abs(real(p(i)))+abs(imag(p(i)))+abs(real(q(i)))+abs(imag(q(i)))>1d2) then
          j       = max(i,maxgrid)
          w       = ( conjg(p0(i))*p(i)  + conjg(q0(i))*q(i)  )  /
     &              ( conjg(p0(i))*p0(i) + conjg(q0(i))*q0(i) ) !; write(1001,*) '#scale',i,w
          p(i:j)  = p(i:j) - w * p0(i:j) !; write(1001,*) '#->',p(i)
          q(i:j)  = q(i:j) - w * q0(i:j)
          pp      = pp    - w * pp0 !; write(1001,*) '#->',pp
          qp      = qp    - w * qp0 !; write(1001,*) '#->',qp
        endif

        rC      = inci  * rA
        rB      = inciB * rA
        dr      = rC - rA

        rCi     = inc  * rAi
        rBi     = incB * rAi

        k1      = pp
        l1      = qp
        k01     = pp0
        l01     = qp0

        v       = ( v0(i) + v0(i-1)*inci ) * incB / 2 !; write(*,*) v0(i),v0(i-1),v
c        v = exp ( ( log(-v0(i)) + log(-v0(i-1)*inci) ) / 2 )
        m       = 2 + ( e - v ) * cci
        wh      = ll/(m*rB**2)
        w       = wh + v - e
        ph      = p(i) + k1*dr/2  ; inh1 =   (qi(i)+qi(i-1))/2 * ci
        qh      = q(i) + l1*dr/2  ; inh2 = - (pi(i)+pi(i-1))/2 * ( 1 + wh*cci/m )
        k2      =   ph*rBi + m*qh + inh1
        l2      = - qh*rBi + w*ph + inh2
        ph      = p0(i) + k01*dr/2
        qh      = q0(i) + l01*dr/2
        k02     =   ph*rBi + m*qh
        l02     = - qh*rBi + w*ph

        ph      = p(i) + k2*dr/2
        qh      = q(i) + l2*dr/2
        k3      =   ph*rBi + m*qh + inh1
        l3      = - qh*rBi + w*ph + inh2
        ph      = p0(i) + k02*dr/2
        qh      = q0(i) + l02*dr/2
        k03     =   ph*rBi + m*qh
        l03     = - qh*rBi + w*ph

        v       = v0(i-1)
        m       = 2 + ( e - v ) * cci
        wh      = ll/(m*rC**2)
        w       = wh + v - e
        ph      = p(i) + k3*dr    ; inh1 =   qi(i-1) * ci
        qh      = q(i) + l3*dr    ; inh2 = - pi(i-1) * ( 1 + wh*cci/m )
        k4      =   ph*rCi + m*qh + inh1
        l4      = - qh*rCi + w*ph + inh2
        ph      = p0(i) + k03*dr
        qh      = q0(i) + l03*dr
        k04     =   ph*rCi + m*qh
        l04     = - qh*rCi + w*ph

        p(i-1)  = p(i) + dr * sixi * (k1+k2+k2+k3+k3+k4) !; write(*,*) i,abs(p(i)),k1,k2,k3,k4
        q(i-1)  = q(i) + dr * sixi * (l1+l2+l2+l3+l3+l4)
        pp      =   p(i-1)*rCi + m*q(i-1) + inh1
        qp      = - q(i-1)*rCi + w*p(i-1) + inh2
        p0(i-1) = p0(i) + dr * sixi * (k01+k02+k02+k03+k03+k04) !; ; write(*,*) i,abs(p0(i)),k1_,k2_,k3_,k4_ ; read(*,*)
        q0(i-1) = q0(i) + dr * sixi * (l01+l02+l02+l03+l03+l04)
        pp0     =   p0(i-1)*rCi + m*q0(i-1)
        qp0     = - q0(i-1)*rCi + w*p0(i-1)

        rA      = rC
        rAi     = rCi
      enddo

      q  = q * ci
      pt = p(nt)
      qt = q(nt)
      if(nt<=maxgrid) then
        i          = min(n,maxgrid)
        pout(nt:i) = p(nt:i)
        qout(nt:i) = q(nt:i)
        if(n<maxgrid) then
          pout(n+1:) = 0
          qout(n+1:) = 0
        endif
      endif

c      r = rAA
c      do i = n,nt+1,-1
c        write(1003,'(9ES40.18)') r,p(i),q(i),p0(i),q0(i)
c        r = r * inci
c      enddo

      deallocate ( v0,p,q,pi,qi,p0,q0 )

      end

c ---------------------------

c     Returns turning point rt=r(nt) and mesh end point rmax=r(nmax)

      subroutine core_turning_point(rt,nt,l,e,itype,ispin)

      use global, only: grid,vmt

      use, intrinsic :: iso_fortran_env
      implicit none
      integer,    intent(in)  :: l,itype,ispin
      integer,    intent(out) :: nt
      real_dp,    intent(out) :: rt
      real_dp,    intent(in)  :: e
      real_dp                 :: inc
      real_dp                 :: r,rA,rB,v,va,vb
      logical                 :: ldum
      integer                 :: ll,n

      ! prepare potential continuation ( va * exp(vb*r)/r )
      inc = exp(grid(itype)%increment)
      rA  = grid(itype)%radius
      rB  = rA / inc
      n   = grid(itype)%number
      va  = vmt(n,  1,itype,ispin) * rA
      vb  = vmt(n-1,1,itype,ispin) * rB
      vb  = log(va/vb) / ( rA-rB )
      va  = va / exp(vb*rA)

      ! determine classical turning point rt
      ll   = l*(l+1)
      ldum = ll==0
      n    = 1
      r    = grid(itype)%first
      do
        if(n<=grid(itype)%number) then ; v = vmt(n,1,itype,ispin)
        else                           ; v = -exp ( va + vb*r ) / r
        endif
        v = v + ll/(2*r**2)! ; if(n>1000) then ; write(*,*) r,v,real(e),ldum; read(*,*) ; endif
        if(v<real(e)) ldum = .true.
        if(v>real(e).and.ldum) exit
        n = n + 1
        r = r * inc
      enddo
      rt = r
      nt = n

c      ! determine mesh end point rmax
c      do
c        if(r>15*rt) exit
c        n = n + 1
c        r = r + inc
c      enddo
c      rmax = r
c      nmax = n

      end

c ---------------------------

c     Determines core states with scalar-relativistic approximation.

      subroutine corestate(pout,qout,e,l,itype,ispin)
      use global, only: maxgrid,gridtype,grid
      use, intrinsic :: iso_fortran_env
      implicit none
      real_dp,    intent(out)   :: pout(maxgrid),qout(maxgrid)
      real_dp,    intent(inout) :: e
      integer,    intent(in)    :: l,itype,ispin
      complex_dp, allocatable   :: p(:),q(:)
      complex_dp                :: pt,qt,pt1o,qt1o,pt1i,qt1i,pt2o,qt2o,pt2i,qt2i,pt3o,qt3o,pt3i,qt3i
      complex_dp                :: e1,e2,e3
      real_dp                   :: rt
      real_dp                   :: r,k1,k2,k3,dk
      integer                   :: nt
      integer                   :: maxgrid0
      type(gridtype)            :: grid1
      real_dp                   :: intgr

      ! define extended grid for intgr (normalization)
      call core_turning_point(rt,nt,l,e,itype,ispin)
      maxgrid0        = maxgrid
      grid1%number    = 1
      grid1%increment = grid(itype)%increment
      grid1%first     = grid(itype)%first
      r               = grid(itype)%first
      do while ( r<=15*rt )
        r             = r * exp(grid(itype)%increment)
        grid1%number  = grid1%number + 1
      enddo
      maxgrid = grid1%number
      allocate ( p(maxgrid),q(maxgrid) )

      ! determine core eigenvalue
      e1 = e
      call core_dirac_hom_out(p,q,pt,qt,rt,nt,l,e1,itype,ispin)
      e1 = e - 0.01d0
      e2 = e + 0.01d0
      call core_dirac_hom_out(p,q,pt1o,qt1o,rt,nt,l,e1,itype,ispin)
      call core_dirac_hom_in (p,q,pt1i,qt1i,rt,nt,l,e1,itype,ispin)
      call core_dirac_hom_out(p,q,pt2o,qt2o,rt,nt,l,e2,itype,ispin)
      call core_dirac_hom_in (p,q,pt2i,qt2i,rt,nt,l,e2,itype,ispin)
      qt1o = qt1o * pt / pt1o ! define all q's so that the p's are continuous
      qt2o = qt2o * pt / pt2o !
      qt1i = qt1i * pt / pt1i !
      qt2i = qt2i * pt / pt2i !
      k1   = qt1i - qt1o
      k2   = qt2i - qt2o
 1    dk   = (k2-k1) / (e2-e1)
      e3   = e1 - k1 / dk
      call core_dirac_hom_out(p,q,pt3o,qt3o,rt,nt,l,e3,itype,ispin)
      call core_dirac_hom_in (p,q,pt3i,qt3i,rt,nt,l,e3,itype,ispin)
      qt3o = qt3o * pt / pt3o
      qt3i = qt3i * pt / pt3i
      k3   = qt3i - qt3o
      if(abs(k3)<1d-8) then
        e        = real(e3)
        p(:nt-1) = p(:nt-1) * pt / pt3o
        q(:nt-1) = q(:nt-1) * pt / pt3o
        p(nt:)   = p(nt:)   * pt / pt3i
        q(nt:)   = q(nt:)   * pt / pt3i
        if(sum(abs(imag(p)))+sum(abs(imag(q)))>1d-12) Error('Nonzero imaginary part.')
        rt       = sqrt ( intgr(real(p)**2+real(q)**2,grid1) )
        p        = p / rt
        q        = q / rt
        maxgrid  = maxgrid0
        pout     = p(:maxgrid)
        qout     = q(:maxgrid)
        deallocate ( p,q )
        return
      endif
      if(abs(k3)>max(abs(k1),abs(k2))) Error('Not converged.')
      if(abs(k1)>abs(k2)) then
        k1 = k3
        e1 = e3
      else
        k2 = k3
        e2 = e3
      endif
      goto 1

      end


























































































      subroutine core_dirac_inhom_in1(pout,qout,pt,qt,rt,nt,alpha,l,e,pc,qc,bm,itype,ispin)

      use global, only: grid,maxgrid,vmt,clight

      use, intrinsic :: iso_fortran_env
      implicit none
      integer,    intent(in)  :: l,itype,ispin
      integer,    intent(in)  :: nt,alpha
      real_dp,    intent(in)  :: rt
      real_dp,    intent(in)  :: pc(maxgrid),qc(maxgrid),bm(maxgrid)
      complex_dp, intent(in)  :: e
      complex_dp, intent(out) :: pout(maxgrid),qout(maxgrid),pt,qt
      complex_dp              :: m,w,wh,ph,qh,inh1,inh2
      complex_dp              :: k1,k2,k3,k4,l1,l2,l3,l4
      complex_dp              :: k01,k02,k03,k04,l01,l02,l03,l04
      complex_dp              :: pp,qp,pp0,qp0,nodes
      complex_dp, allocatable :: p(:),q(:),pi(:),qi(:),p0(:),q0(:)
      real_dp,    allocatable :: v0(:)
      real_dp                 :: h,inc,incB,inci,inciB
      real_dp                 :: ci,cci
      real_dp                 :: r,rA,rB,rC,rAi,rBi,rCi,dr,v,va,vb,ca1,cb1,ca2,cb2,bm0,bm1 !,rAA
      real_dp, parameter      :: sixi = 0.166666666666667d0
      integer                 :: i,j,ll,n

      if(nt==0) Error('turning point (rt,nt) undefined.')

      ! prepare potential continuation ( va * exp(vb*r)/r )
      h   = grid(itype)%increment
      inc = exp(h)
      rA  = grid(itype)%radius
      rB  = rA / inc
      n   = grid(itype)%number
      va  = vmt(n,  1,itype,ispin) * rA
      vb  = vmt(n-1,1,itype,ispin) * rB
      vb  = log(va/vb) / ( rA-rB )
      va  = va / exp(vb*rA)

      ! continuation of bm(:)
      bm0 =   bm(n) / rA
      bm1 = ( bm(n) - bm(n-1) ) / (rA-rB)

      ! prepare core state (pc/qc) continuation ( r**alpha * ca1 * exp(cb1*r), ca2,cb2 for qc )
      ca1 = pc(n)   / rA**(alpha+1)
      cb1 = pc(n-1) / rB**(alpha+1)
      cb1 = log(ca1/cb1) / ( rA-rB )
      ca1 = ca1 / exp(cb1*rA)
      ca2 = qc(n)   / rA**(alpha+1)
      cb2 = qc(n-1) / rB**(alpha+1)
      cb2 = log(ca2/cb2) / ( rA-rB )
      ca2 = ca2 / exp(cb2*rA)

      ! extend mesh to 15*rt
      n = nt
      r = rt
      do while(r<15*rt)!grid(itype)%radius+10)
        n = n + 1
        r = r * inc
      enddo

      ! define potential and pi/qi up to 15*rt
      allocate ( v0(n),p(n),q(n),pi(n),qi(n),p0(n),q0(n) )
      bm0 = bm(grid(itype)%number) / grid(itype)%radius
      bm1 = ( bm(grid(itype)%number) - bm(grid(itype)%number-1) ) / (rA-rB)
      r   = grid(itype)%first
c      write(1004,*) '#',bm0,bm1
      do i = 1,n
        if(i<=grid(itype)%number) then
          v0(i) = vmt(i,1,itype,ispin) * r
          pi(i) = pc(i) * bm(i)/r
          qi(i) = qc(i) * bm(i)/r
c          write(1004,'(6ES40.20)') r,pi(i),qi(i),bm(i)
        else
c          v0(i) = va*r + vb
          v0(i) = va * exp ( vb*r )
          pi(i) = r**(alpha+1) * ca1 * exp(cb1*r) * ( bm0 + bm1*(r-rA)/r )
          qi(i) = r**(alpha+1) * ca2 * exp(cb2*r) * ( bm0 + bm1*(r-rA)/r )
c          write(1004,'(5ES40.20)') r,pi(i),qi(i)
        endif
        r = r * inc
      enddo
      r = r / inc

      ll    = l*(l+1)
      incB  = sqrt(inc)
      inci  = 1/inc
      inciB = 1/incB
      ci    = 1/clight
      cci   = ci*ci
      rA    = r
      rAi   = 1/rA

      ! start values ( q/p=-sqrt(-e/m) )
      m     = ( 2 + cci*(e-v0(n)/r) ) * r
      w     = ll/m + v0(n) - e*r
      p0(n) = 1d-8
      q0(n) = 1d-8 * (-sqrt(-e/m))
      pp0   =  p0(n) + m * q0(n)
      qp0   = -q0(n) + w * p0(n)
      p(n)  = p0(n)
      q(n)  = q0(n)
      pp    = pp0 + r * qi(n) * ci
      qp    = qp0 - r * pi(n) !* ( 1 + ll*cci/m**2 )

c      write(*,*) 'm',m
c      write(*,*) 'PP',pp
c      write(*,*) 'QP',qp
c      write(*,*) 't1',r*qi(n)*ci
c      write(*,*) 't2',r*pi(n)

c      write(*,*) 'comp',m,q0(n),p0(n),r*qi(n)*ci

c      rAA = r
      rewind(1001)

      ! Runge-Kutta
      do i = n,nt+1,-1 ; write(1001,*) rA,real(p0(i)),real(p(i)) !; write(1001,*) i,abs(p(i)),abs(pi(i)),abs(v0(i))

c        if(abs(real(p0(i)))+abs(imag(p0(i)))+abs(real(q0(i)))+abs(imag(q0(i)))>1d2) then
c          p0(i:) = p0(i:) / 1d4
c          q0(i:) = q0(i:) / 1d4
c          pp0    = pp0    / 1d4
c          qp0    = qp0    / 1d4
c        endif
        m = conjg(p(i)) * p(i) + conjg(q(i)) * q(i)
        if(real(m)>1d6) then
c        if(abs(real(p(i)))+abs(imag(p(i)))+abs(real(q(i)))+abs(imag(q(i)))>1d2.and..true.) then
          m = conjg(p0(i)) * p0(i) + conjg(q0(i)) * q0(i)
          j       = max(i,maxgrid)
          w = conjg(p(i)) * p0(i) + conjg(q(i)) * q0(i)
c          w       = (p(i)*conjg(p0(i))+q(i)*conjg(q0(i))) / (abs(p0(i))**2+abs(q0(i))**2) ; write(1001,*) '#scale',i,w
          w = -real(w) + (0d0,1d0) * imag(w)
          w = w / m
          p(i:j)  = p(i:j) - w * p0(i:j) ; write(1001,*) '#->',p(i)
          q(i:j)  = q(i:j) - w * q0(i:j)
          pp      = pp    - w * pp0 ; write(1001,*) '#->',pp
          qp      = qp    - w * qp0; write(1001,*) '#->',qp
        endif

        rC      = inci  * rA
        rB      = inciB * rA
        dr      = rC - rA

        rCi     = inc  * rAi
        rBi     = incB * rAi

        k1      = pp
        l1      = qp
        k01     = pp0
        l01     = qp0

        v       = ( v0(i) + v0(i-1) ) / 2 !; write(*,*) v0(i),v0(i-1),v
        m       = 2*rB + ( e*rB - v ) * cci
        wh      = ll/m
        w       = wh + v - e*rB
        ph      = p(i) - k1*h/2 ; inh1 =   (qi(i)+qi(i-1))/2 * ci               * rB
        qh      = q(i) - l1*h/2 ; inh2 = - (pi(i)+pi(i-1))/2 * rB!( 1 + wh*cci/m ) * rB
        k2      =   ph + m*qh + inh1
        l2      = - qh + w*ph + inh2
        ph      = p0(i) - k01*h/2
        qh      = q0(i) - l01*h/2
        k02     =   ph + m*qh
        l02     = - qh + w*ph

        ph      = p(i) - k2*h/2
        qh      = q(i) - l2*h/2
        k3      =   ph + m*qh + inh1
        l3      = - qh + w*ph + inh2
        ph      = p0(i) - k02*h/2
        qh      = q0(i) - l02*h/2
        k03     =   ph + m*qh
        l03     = - qh + w*ph

        v       = v0(i-1)
        m       = 2*rC + ( e*rC - v ) * cci
        wh      = ll/m
        w       = wh + v - e*rC
        ph      = p(i) - k3*h ; inh1 =   qi(i-1) * ci               * rC
        qh      = q(i) - l3*h ; inh2 = - pi(i-1) * rC!( 1 + wh*cci/m ) * rC
        k4      =   ph + m*qh + inh1
        l4      = - qh + w*ph + inh2
        ph      = p0(i) - k03*h
        qh      = q0(i) - l03*h
        k04     =   ph + m*qh
        l04     = - qh + w*ph

c        write(*,*) k1,k2,k3,k4
c        write(*,*)

        p(i-1)  = p(i) - h * sixi * (k1+k2+k2+k3+k3+k4) !; write(*,*) i,abs(p(i)),k1,k2,k3,k4
        q(i-1)  = q(i) - h * sixi * (l1+l2+l2+l3+l3+l4)
        pp      =   p(i-1) + m*q(i-1) + inh1
        qp      = - q(i-1) + w*p(i-1) + inh2
        p0(i-1) = p0(i) - h * sixi * (k01+k02+k02+k03+k03+k04) !; ; write(*,*) i,abs(p0(i)),k1_,k2_,k3_,k4_ ; read(*,*)
        q0(i-1) = q0(i) - h * sixi * (l01+l02+l02+l03+l03+l04)
        pp0     =   p0(i-1) + m*q0(i-1)
        qp0     = - q0(i-1) + w*p0(i-1)

        rA      = rC
        rAi     = rCi
      enddo

      q  = q * ci
      pt = p(nt)
      qt = q(nt)
      if(nt<=maxgrid) then
        i          = min(n,maxgrid)
        pout(nt:i) = p(nt:i)
        qout(nt:i) = q(nt:i)
        if(n<maxgrid) then
          pout(n+1:) = 0
          qout(n+1:) = 0
        endif
      endif

c      if(imag(e)==0) return
      return
      q = q / ci
c      write(*,*) imag(e)
      call cradsra_inward_inh1(
     &  e,l,v0,grid(itype)%first,grid(itype)%increment,n,
     &  nt,n,clight,
     &  p(n),q(n),[(1d0,i=1,n)],real(pi),real(qi),
     &  p0,q0*ci,
     &  nodes,p,q)

c      write(*,*) pt,p(nt)
c      write(*,*) qt,q(nt)

      pt = p(nt)
      qt = q(nt)
      if(nt<=maxgrid) then
        i          = min(n,maxgrid)
        pout(nt:i) = p(nt:i)
        qout(nt:i) = q(nt:i)
        if(n<maxgrid) then
          pout(n+1:) = 0
          qout(n+1:) = 0
        endif
      endif

c      r = rAA
c      do i = n,nt+1,-1
c        write(1003,'(9ES40.18)') r,p(i),q(i),p0(i),q0(i)
c        r = r * inci
c      enddo

      deallocate ( v0,p,q,pi,qi,p0,q0 )

      end


# define CORR 1


      SUBROUTINE cradsra_inward_inh1(
     >                             e,l,vr,r0,h,jri,ictp,jmtd,c,
     >                             p_in,q_in,inh,p0,q0,
     >                             p_hom,q_hom,
     <                             nodes,p,q)

      use, intrinsic :: iso_fortran_env
      IMPLICIT NONE

      ! - scalars -
      INTEGER, INTENT (IN) :: l
      INTEGER, INTENT (IN) :: jmtd,jri,ictp
      real_dp, INTENT (IN) :: h,r0
      real_dp, INTENT (IN) :: c
      complex_dp, INTENT (IN) :: e
      complex_dp, INTENT (OUT):: nodes
      complex_dp, INTENT (IN) :: p_in,q_in

      ! - arrays -
      real_dp, INTENT (IN) :: vr(jmtd)
      real_dp, INTENT (IN) :: inh(jmtd)
      real_dp, INTENT (IN) :: p0(jmtd),q0(jmtd)
      complex_dp, INTENT (IN) :: p_hom(jmtd),q_hom(jmtd)
      complex_dp, INTENT (inOUT):: p(jmtd),q(jmtd)

      ! - local scalars -
      INTEGER                 :: i,j
      real_dp                 :: dr,drh
      real_dp                 :: fl1
      real_dp                 :: r,rh,rh2
      real_dp                 :: norm0,norm
      real_dp                 :: cin,cin2
      complex_dp              :: sk1,sk2,sk3,sk4
      complex_dp              :: sl1,sl2,sl3,sl4
      complex_dp              :: t1,t2
      complex_dp              :: rm,rve
      complex_dp              :: yn,zn
      complex_dp              :: cfac
      complex_dp              :: cdum


      ! - local arrays -
      real_dp                 :: rmsh(jmtd)
      complex_dp              :: pp(jmtd),qp(jmtd)

      cin  = 1.0/c
      cin2 = cin*cin

      IF (jri.GT.jmtd) STOP 'radsra_inward_inh'

      ! setup radial mesh
      rmsh(1) = r0
      dr      = exp(h)
      DO i = 2,jri
        rmsh(i) = rmsh(i-1)*dr
      END DO

      ! set up initial conditions
      fl1     = l* (l+1)
      rm      = 2.*rmsh(jri) + cin2* (rmsh(jri)*e-vr(jri))

      p(jri)  = p_in !; write(*,*) p_in,p(jri)
      q(jri)  = q_in !; write(*,*) q_in,q(jri)

      t1      = rmsh(jri)*inh(jri)*q0(jri)*cin  ! a factor 1/c is already included in q0
      t2      = rmsh(jri)*inh(jri)*p0(jri)*CORR

      pp(jri) = rm*q(jri)                           + p(jri) + t1
      qp(jri) = (fl1/rm+vr(jri)-rmsh(jri)*e)*p(jri) - q(jri) - t2

c      write(*,*) 'RM',rm
c      write(*,*) 'PP',pp(jri)
c      write(*,*) 'QP',qp(jri)
c      write(*,*) 't1',t1
c      write(*,*) 't2',t2

c      write(*,*) 'comp',rm,q(jri),p(jri),t1

      ! use 4th order runge-kutta to get first few mesh points by
      ! inward integration
      nodes = 0
      dr    = exp(-h)
      drh   = sqrt(dr)
      r     = rmsh(jri)
      DO i = jri,ictp,-1
        rh2 = drh*r
        rh  = dr*r
        ! calculate k1
        sk1 = h*pp(i)
        sl1 = h*qp(i)
        ! calculate k2
        rve = 0.5* (vr(i-1)+vr(i)) - rh2*e
        rm  = 2.*rh2 - cin2*rve
        yn  = p(i) - 0.5*sk1
        zn  = q(i) - 0.5*sl1
        t1  = 0.5*rh2*(inh(i-1)*q0(i-1)+inh(i)*q0(i))*cin
        t2  = 0.5*rh2*(inh(i-1)*p0(i-1)+inh(i)*p0(i))*CORR
!         t2  = 0.5*rh2*( inh(i  )*p0(i  )*(1+cin2*fl1/(2*r +cin2*(r *e-vr(i  ))))
!      & +inh(i-1)*p0(i-1)*(1+cin2*fl1/(2*rh+cin2*(rh*e-vr(i-1)))) )
        sk2 = h* (rm*zn           + yn + t1)
        sl2 = h* ((fl1/rm+rve)*yn - zn - t2)
        ! calculate k3
        yn  = p(i) - 0.5*sk2
        zn  = q(i) - 0.5*sl2
        sk3 = h* (rm*zn           + yn + t1)
        sl3 = h* ((fl1/rm+rve)*yn - zn - t2)
        ! calculate k4
        rve = vr(i-1) - rh*e
        rm  = 2.*rh - cin2*rve
        yn  = p(i) - sk3
        zn  = q(i) - sl3
        t1  = rh*inh(i-1)*q0(i-1)*cin
        t2  = rh*inh(i-1)*p0(i-1)*CORR
        sk4 = h* (rm*zn           + yn + t1)
        sl4 = h* ((fl1/rm+rve)*yn - zn - t2)
        ! calculate p/q at i-1

c        write(*,*) sk1/h,sk2/h,sk3/h,sk4/h
c        write(*,*) p(i-1),q(i-1)

        p(i-1)  = p(i) - (sk1+2.*sk2+2.*sk3+sk4)/6.
        q(i-1)  = q(i) - (sl1+2.*sl2+2.*sl3+sl4)/6.

c        write(*,*) p(i-1),q(i-1) ; read(*,*)

        norm = conjg(p(i-1))*p(i-1) + conjg(q(i-1))*q(i-1)
        IF( sqrt(norm) > 1d6 .and..true.) THEN
          norm0   = conjg(p_hom(i-1))*p_hom(i-1) + c*c*conjg(q_hom(i-1))*q_hom(i-1)
          IF( norm0 > 1d-12 ) THEN

            cdum   = conjg(p(i-1))*p_hom(i-1) + conjg(q(i-1))*c*q_hom(i-1)

            cfac   = CMPLX( -real(cdum),imag(cdum),kind=real64 )
            cfac   = cfac/norm0

            DO j = i-1,jri

              p(j) = p(j) + cfac *   p_hom(j)
              q(j) = q(j) + cfac *c* q_hom(j)

            END DO
          ELSE
            DO j = i-1,jri
              p(j) = 0.
              q(j) = 0.
            END DO
          END IF
        END IF

        pp(i-1) = rm*q(i-1)           + p(i-1) + t1
        qp(i-1) = (fl1/rm+rve)*p(i-1) - q(i-1) - t2

        nodes = nodes +
     & CMPLX(int(0.501*abs(sign(1d0,real(p(i-1)))-sign(1d0,real(p(i))))),
     & int(0.501*abs(sign(1d0,imag(p(i-1)))-sign(1d0,imag(p(i))))))

        r       = rh
      END DO

      DO i = 1,jri
        q(i) = cin*q(i)
      END DO

      END SUBROUTINE cradsra_inward_inh1

