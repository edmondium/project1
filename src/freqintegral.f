c--------------------------------------------------------------------------------
c Copyright (c) 2021 Peter Grünberg Institut, Forschungszentrum Jülich, Germany
c This file is part of SPEX and available as free software under the conditions
c of the MIT license as expressed in the LICENSE file in more detail.
c--------------------------------------------------------------------------------

# include "cppmacro.h"

      module freq_integral
      use, intrinsic :: iso_fortran_env

      interface freqintegral
      module procedure  freqintegral_r,freqintegral_c
      end interface

      interface spline_init
      module procedure  spline_init_r,spline_init_c
      end interface

      contains


c     Initializes the evaluation of the integral
c
c            infinity
c     I(c) =   INT    y(w') / (w'+c) dw'
c           -infinity
c
c     c = farg(1:nfarg); y(w') is symmetric around w'=0
c
      subroutine freqintegral_init(wfreqintegral,freq,n,farg,narg)
      use global, only: img
      implicit none
      integer,    intent(in)  :: n,narg
      real_dp,    intent(in)  :: freq(n)
      complex_dp, intent(in)  :: farg(narg)
      complex_dp, intent(out) :: wfreqintegral(0:3,n,narg)
      complex_dp              :: f(0:3,n)
      complex_dp              :: c,c1,c1i,c2i,lnd,tail(2),lnc(n),integ1,integ2,integ3
      real_dp                 :: lnx(n-1),xi(n-1),xi2(n-1),xi3(n-1),mat(2,2),r,logr,x1,rr,r3i
      real_dp                 :: a0,a1,a2,r21,r22,r31,r32,ri21,ri22,ri31,ri32,ri41,ri42,ri51,ri52
      integer                 :: i,j,k,sgn
      if(freq(1)/=0) Bug('freq(1) not zero.')
      do i = 1,n-1
        lnx(i) = log(1+freq(i+1))       - log(1+freq(i))
        xi(i)  =    (1+freq(i+1))**(-1) -    (1+freq(i))**(-1)
        xi2(i) =  ( (1+freq(i+1))**(-2) -    (1+freq(i))**(-2) ) / 2
        xi3(i) =  ( (1+freq(i+1))**(-3) -    (1+freq(i))**(-3) ) / 3
      enddo
      r        = freq(n)
      mat(1,1) =    r**(-2) ; mat(1,2) =    r**(-4)
      mat(2,1) = -2*r**(-3) ; mat(2,2) = -4*r**(-5)
      rr       =   mat(1,1) * mat(2,2) - mat(1,2) * mat(2,1) ! -> mat = invert(mat)
      x1       =   mat(1,1)
      mat(1,1) =   mat(2,2) / rr
      mat(2,1) = - mat(2,1) / rr
      mat(1,2) = - mat(1,2) / rr
      mat(2,2) =         x1 / rr
      logr     = log(r)
      rr       = (1+r)**(-2)
      r3i      = 1/r**3
      x1       = r / (r+1)
      do i = 1,narg
        f    = 0
        sgn  = 1
        c    = farg(i)
        c2i  = c**(-2)
 1      c1   = c - 1
        c1i  = 1/c1
        do j = 1,n
          a0 = abs(freq(j)+c)
          if          (a0==0) then ; lnc(j) = 0
          else if(imag(c)==0) then ; lnc(j) = log(a0)
          else                     ; lnc(j) = log(a0) + img * atan2(imag(c),real(c)+freq(j)) ! complex log
          endif
        enddo
        ! intervals (freq(j)->freq(j+1))
        if(abs(c1)<1d-3) then ! Taylor
          r22  = 0
          r32  = 0
          ri22 = 1
          ri32 = 1
          ri42 = 1
          ri52 = 1
          do j = 1,n-1
            lnd    = lnc(j+1) - lnc(j)
            r21    = r22
            r31    = r32
            r22    = freq(j+1)**2
            r32    = freq(j+1)*r22
            ri21   = ri22
            ri31   = ri32
            ri41   = ri42
            ri51   = ri52
            ri52   = 1/(1+freq(j+1))
            ri22   = ri52**2
            ri32   = ri22*ri52
            ri42   = ri32*ri52
            ri52   = ri42*ri52
            f(0,j) = f(0,j) + sgn * lnd
            a0     = xi(j) + lnx(j)
            a1     = xi2(j) - xi(j)
            a2     = xi3(j) - xi2(j)
            f(1,j) = f(1,j) + sgn * ( a0 - c1 * ( a1 - c1*a2 ) )
            a0     = a0     - ( r22*ri22 - r21*ri21 ) / 2
            a1     = a1*2/3 - ( r22*ri32 - r21*ri31 ) / 3
            a2     = a2/2   - ( r22*ri42 - r21*ri41 ) / 4
            f(2,j) = f(2,j) + sgn * ( a0 - c1 * ( a1 - c1*a2 ) )
            a0     = a0     - ( r32*ri32 - r31*ri31 ) / 3
            a1     = a1*3/4 - ( r32*ri42 - r31*ri41 ) / 4
            a2     = a2*3/5 - ( r32*ri52 - r31*ri51 ) / 5
            f(3,j) = f(3,j) + sgn * ( a0 - c1 * ( a1 - c1*a2 ) )
          enddo
        else                     ! exact
          do j = 1,n-1
            lnd    = lnc(j+1) - lnc(j)
            integ1 =   ( lnx(j) - lnd    ) * c1i
            integ2 = - ( xi(j)  + integ1 ) * c1i
            integ3 = - ( xi2(j) + integ2 ) * c1i
            integ2 =   integ2 -   integ1
            integ3 = - integ3 + 2*integ2 + integ1
            if(sgn==1) then
              f(0,j) =          lnd ; lnd = lnd - integ1
              f(1,j) =          lnd ; lnd = lnd + integ2
              f(2,j) =          lnd ; lnd = lnd + integ3
              f(3,j) =          lnd
            else
              f(0,j) = f(0,j) - lnd ; lnd = lnd - integ1
              f(1,j) = f(1,j) - lnd ; lnd = lnd + integ2
              f(2,j) = f(2,j) - lnd ; lnd = lnd + integ3
              f(3,j) = f(3,j) - lnd
            endif
          enddo
        endif
        ! tail contribution (freq(n)->infinity)
        if(abs(c)<1d-1) then ! Taylor
          tail = 0
          do k = 10,0,-1
            tail(1) = -c/r * tail(1) + 1d0 / (k+2)    ! tail w**(-2)
            tail(2) = -c/r * tail(2) + 1d0 / (k+4)    ! tail w**(-4)
          enddo
          tail(1) = tail(1) / r**2
          tail(2) = tail(2) / r**4
        else                    ! exact
          c1      = c/r + logr - lnc(n)
          tail(1) =   c1 * c2i                           ! tail w**(-2)
          tail(2) = ( c1 * c2i + (c/3-r/2) * r3i ) * c2i ! tail w**(-4)
        endif
        c1      = tail(1) * mat(1,1) + tail(2) * mat(2,1)
        tail(2) = tail(1) * mat(1,2) + tail(2) * mat(2,2)
        tail(1) = c1
        f(0,n)  = f(0,n) + sgn *   tail(1)
        f(1,n)  = f(1,n) + sgn * ( tail(1) * x1 + tail(2) * rr   )
        f(2,n)  = f(2,n) + sgn * ( tail(1) * x1 + tail(2) * rr*2 ) * x1
        f(3,n)  = f(3,n) + sgn * ( tail(1) * x1 + tail(2) * rr*3 ) * x1**2
        ! add negative axis
        if(sgn==1) then
          c   = -c
          sgn = -1
          goto 1
        endif
        wfreqintegral(:,:,i) = f
      enddo
      end subroutine freqintegral_init

c ------

      function ln(z)
      implicit none
      complex_dp             :: ln
      complex_dp, intent(in) :: z
      if      (abs(z)==0) then ; ln = 0
      else if(imag(z)==0) then ; ln = log(abs(z))
      else                     ; ln = log(z)
      endif
      end function ln

c ------

      subroutine spline_init_r(aspline,x,y,n)
      use global, only: metal
      implicit none
      integer, intent(in)  :: n
      real_dp, intent(out) :: aspline(0:3,n)
      real_dp, intent(in)  :: x(n)
      real_dp              :: x1(n+1)
      real_dp, intent(in)  :: y(n)
      real_dp              :: y1(n+1),aa(n),bb(n)
      integer              :: i
      x1(:n)  = x/(1+x)
      x1(n+1) = 1d0
      y1(:n)  = y
      y1(n+1) = 0d0
      if(metal) then ; i = 4
      else           ; i = 1
      endif
      call cubicspline(aa,bb,y1,x1,n+1,0d0,i,0d0,1)
      do i = 1,n
        call cubiccoeff(aspline(:,i),aa(i),bb(i),x1(i),y1(i))
      enddo
      end subroutine spline_init_r

      subroutine spline_init_c(aspline,x,y,n)
      use global, only: metal
      implicit none
      integer,    intent(in)  :: n
      complex_dp, intent(out) :: aspline(0:3,n)
      real_dp,    intent(in)  :: x(n)
      real_dp                 :: x1(n+1)
      complex_dp, intent(in)  :: y(n)
      complex_dp              :: y1(n+1),aa(n),bb(n)
      integer                 :: i
      x1(:n)  = x/(1+x)
      x1(n+1) = 1d0
      y1(:n)  = y
      y1(n+1) = 0d0
      if(metal) then ; i = 4
      else           ; i = 1
      endif
      call cubicspline_c(aa,bb,y1,x1,n+1,(0d0,0d0),i,(0d0,0d0),1)
      do i = 1,n
        call cubiccoeff_c(aspline(:,i),aa(i),bb(i),x1(i),y1(i))
      enddo
      end subroutine spline_init_c

c ------

c     Calculates the mth derivative of the integral
c
c                      infinity
c     I(c) = 1/(2*pi*i)   INT    y(w') / (w'+c) dw'
c                     -infinity
c
c     m > 0 : Finite differences are used.
c
c     mode=0 : integrate only from 0 to infinity
c     mode=1 : y(-w) = y(w)  is assumed
c     mode=2 : y(-w) = y*(w) is assumed (only relevant for freqintegral_c)
c
c     Due to the branch cut of the complex logarithm (along the negative real axis),
c     there is a sign change: real[I(w+i*eta)] = -real[I(w-i*eta)], where w and eta are
c     real and eta is infinitesimal. We define real[I(w)] = 0 (i.e., to be the average).
c     (Note the imaginary prefactor!)
c
c     Slow routine. No initialization necessary.
c
      recursive function freqintegral_r(y,x,n,mode,c,m) result(freqintegral_out)
      use global, only: pi,img,metal
      use wrapper
      implicit none
      complex_dp              :: freqintegral_out
      complex_dp,  intent(in) :: c
      integer,     intent(in) :: n,m,mode
      real_dp,     intent(in) :: x(n)
      real_dp                 :: x1(n+1),xa(18),mat(2,2),r,r1
      real_dp,     intent(in) :: y(n)
      real_dp                 :: y1(n+1),aa(n),bb(n),a(0:3,n)
      complex_dp              :: cc,c1,f,integ(0:3,0:1),integ1(0:3),tail(2)
      integer                 :: i,j,sgn
      if(m>=1) then
        freqintegral_out = ( freqintegral_r(y,x,n,mode,c+1d-4,m-1) - freqintegral_r(y,x,n,mode,c-1d-4,m-1) ) / 2d-4
        return
      else if(m<0) then
        Bug('6th argument out of range.')
      endif
      x1(:n)  = x/(1+x)
      x1(n+1) = 1
      y1(:n)  = y
      y1(n+1) = 0
      if(metal) then ; i = 4
      else           ; i = 1
      endif
      call cubicspline(aa,bb,y1,x1,n+1,0d0,i,0d0,1)
      do i = 1,n
        call cubiccoeff(a(:,i),aa(i),bb(i),x1(i),y1(i))
      enddo
      sgn = 1
      f   = 0
      cc  = c
 1    c1  = cc - 1
      ! intervals (x(i)->x(i+1))
      do i = 1,n-1
        if(abs(c1)<2d-1) then
          r = 1 ; r1 = 1
          do j = 1,18
            r     = r  / (1+x(i))
            r1    = r1 / (1+x(i+1))
            xa(j) = r1 - r
          enddo
          integ(:,1) = 0
          integ(0,1) = ln(cc+x(i+1)) - ln(cc+x(i))
          do j = 15,0,-1
            integ(1,1) = - integ(1,1)*c1 - xa(j+1)/(j+1)
            integ(2,1) = - integ(2,1)*c1 - xa(j+2)/(j+2)
            integ(3,1) = - integ(3,1)*c1 - xa(j+3)/(j+3)
          enddo
        else
          integ(1,0) =   log(1+x(i+1))       - log(1+x(i))
          integ(2,0) =    - (1+x(i+1))**(-1) +    (1+x(i))**(-1)
          integ(3,0) =  ( - (1+x(i+1))**(-2) +    (1+x(i))**(-2) ) / 2
          integ(0,1) =   ln(cc+x(i+1))       - ln(cc+x(i))
          integ(1,1) = ( integ(1,0)-integ(0,1) ) / c1
          integ(2,1) = ( integ(2,0)-integ(1,1) ) / c1
          integ(3,1) = ( integ(3,0)-integ(2,1) ) / c1
        endif
        integ1(0) = integ(0,1)
        integ1(1) = integ(0,1) -   integ(1,1)
        integ1(2) = integ(0,1) - 2*integ(1,1) +   integ(2,1)
        integ1(3) = integ(0,1) - 3*integ(1,1) + 3*integ(2,1) - integ(3,1)
        f         = f + sgn * dot_product(a(:,i),integ1)
      enddo
      ! tail contribution (freq(n)->infinity)
      r        = x(n)
      mat(1,1) =    r**(-2) ; mat(1,2) =    r**(-4)
      mat(2,1) = -2*r**(-3) ; mat(2,2) = -4*r**(-5)
      mat      = invert(mat)
      if(abs(cc)<2d-1) then ! Taylor
        tail = 0
        do i = 15,0,-1
          tail(1) = tail(1)*(-cc/r) + 1/((i+2)*r**2)         ! tail w**(-2)
          tail(2) = tail(2)*(-cc/r) + 1/((i+4)*r**4)         ! tail w**(-4)
        enddo
      else                  ! exact
        c1      = cc/r + log(r) - ln(cc+r)
        tail(1) =   c1 / cc**2                               ! tail w**(-2)
        tail(2) = ( c1 / cc**2 + (cc/3-r/2) / r**3 ) / cc**2 ! tail w**(-4)
      endif
      tail = matmul ( tail,mat )
      r    = (1+r)**(-2)
      f    = f + sgn * ( a(0,n) *   tail(1)
     &                 + a(1,n) * ( tail(1) * x1(n)    + tail(2) * r            )
     &                 + a(2,n) * ( tail(1) * x1(n)**2 + tail(2) * r*2*x1(n)    )
     &                 + a(3,n) * ( tail(1) * x1(n)**3 + tail(2) * r*3*x1(n)**2 ) )
      ! add negative axis
      if(sgn==1.and.mode>0) then
        cc  = -cc
        sgn = -1
        goto 1
      endif
      if(real(c)==0) then ; freqintegral_out = imag ( f / (2*pi) )   ! argument c on imaginary axis: f is purely imaginary (fulfilled up to machine precision without this if)
      else                ; freqintegral_out =        f / (2*pi*img)
      endif
      end function freqintegral_r

      recursive function freqintegral_c(y,x,n,mode,c,m) result(freqintegral_out)
      use global, only: pi,img,metal
      use wrapper
      implicit none
      complex_dp              :: freqintegral_out
      complex_dp,  intent(in) :: c
      integer,     intent(in) :: n,m,mode
      real_dp,     intent(in) :: x(n)
      real_dp                 :: x1(n+1),xa(18),mat(2,2),r,r1
      complex_dp,  intent(in) :: y(n)
      complex_dp              :: y1(n+1),aa(n),bb(n),a(0:3,n)
      complex_dp              :: cc,c1,f,integ(0:3,0:1),integ1(0:3),tail(2)
      integer                 :: i,j,sgn
      if(m>=1) then
        freqintegral_out = ( freqintegral_c(y,x,n,mode,c+1d-4,m-1) - freqintegral_c(y,x,n,mode,c-1d-4,m-1) ) / 2d-4
        return
      else if(m<0) then
        Bug('6th argument out of range.')
      endif
      x1(:n)  = x/(1+x)
      x1(n+1) = 1
      y1(:n)  = y
      y1(n+1) = 0
      if(metal) then ; i = 4
      else           ; i = 1
      endif
      call cubicspline_c(aa,bb,y1,x1,n+1,(0d0,0d0),i,(0d0,0d0),1)
      do i = 1,n
        call cubiccoeff_c(a(:,i),aa(i),bb(i),x1(i),y1(i))
      enddo
      sgn = 1
      f   = 0
      cc  = c
 1    c1  = cc - 1
      ! intervals (x(i)->x(i+1))
      do i = 1,n-1
        if(abs(c1)<2d-1) then
          r = 1 ; r1 = 1
          do j = 1,18
            r     = r  / (1+x(i))
            r1    = r1 / (1+x(i+1))
            xa(j) = r1 - r
          enddo
          integ(:,1) = 0
          integ(0,1) = ln(cc+x(i+1)) - ln(cc+x(i))
          do j = 15,0,-1
            integ(1,1) = - integ(1,1)*c1 - xa(j+1)/(j+1)
            integ(2,1) = - integ(2,1)*c1 - xa(j+2)/(j+2)
            integ(3,1) = - integ(3,1)*c1 - xa(j+3)/(j+3)
          enddo
        else
          integ(1,0) =   log(1+x(i+1))       - log(1+x(i))
          integ(2,0) =    - (1+x(i+1))**(-1) +    (1+x(i))**(-1)
          integ(3,0) =  ( - (1+x(i+1))**(-2) +    (1+x(i))**(-2) ) / 2
          integ(0,1) =   ln(cc+x(i+1))       - ln(cc+x(i))
          integ(1,1) = ( integ(1,0)-integ(0,1) ) / c1
          integ(2,1) = ( integ(2,0)-integ(1,1) ) / c1
          integ(3,1) = ( integ(3,0)-integ(2,1) ) / c1
        endif
        integ1(0) = integ(0,1)
        integ1(1) = integ(0,1) -   integ(1,1)
        integ1(2) = integ(0,1) - 2*integ(1,1) +   integ(2,1)
        integ1(3) = integ(0,1) - 3*integ(1,1) + 3*integ(2,1) - integ(3,1)
        f         = f + sgn * sum(a(:,i)*integ1)
      enddo
      ! tail contribution (freq(n)->infinity)
      r        = x(n)
      mat(1,1) =    r**(-2) ; mat(1,2) =    r**(-4)
      mat(2,1) = -2*r**(-3) ; mat(2,2) = -4*r**(-5)
      mat      = invert(mat)
      if(abs(cc)<2d-1) then ! Taylor
        tail = 0
        do i = 15,0,-1
          tail(1) = tail(1)*(-cc/r) + 1/((i+2)*r**2)         ! tail w**(-2)
          tail(2) = tail(2)*(-cc/r) + 1/((i+4)*r**4)         ! tail w**(-4)
        enddo
      else                  ! exact
        c1      = cc/r + log(r) - ln(cc+r)
        tail(1) =   c1 / cc**2                               ! tail w**(-2)
        tail(2) = ( c1 / cc**2 + (cc/3-r/2) / r**3 ) / cc**2 ! tail w**(-4)
      endif
      tail = matmul ( tail,mat )
      r    = (1+r)**(-2)
      f    = f + sgn * ( a(0,n) *   tail(1)
     &                 + a(1,n) * ( tail(1) * x1(n)    + tail(2) * r            )
     &                 + a(2,n) * ( tail(1) * x1(n)**2 + tail(2) * r*2*x1(n)    )
     &                 + a(3,n) * ( tail(1) * x1(n)**3 + tail(2) * r*3*x1(n)**2 ) )
      ! add negative axis
      if(sgn==1.and.mode>0) then
        cc  = -cc
        sgn = -1
        if(mode==2) a = conjg(a)
        goto 1
      endif
      freqintegral_out = f / (2*pi*img)
      end function freqintegral_c

c ------

c     Returns the mth derivative of
c
c     I(c) = 1/(2*pi*i) * INT p(w) / (w+c) dw
c          = 1/(2*pi*i) * INT SUM(i=poles) resid(i) / (img*w-pole(i)) / (w+c) dw.
c
c     The function p(w) is defined only for positive w.
c     mode=0 : integrate from 0 to infinity
c     mode=1 : integrate from -infinity to infinity; p(-w) = p(w)  is assumed
c     mode=2 : integrate from -infinity to infinity; p(-w) = p*(w) is assumed
c
c     See comments for freqintegral_r above.
c
      recursive function freqintegral_poles(pole,resid,npole,mode,c,m) result(freqintegral_out)
      use global, only: img,pi
      implicit none
      integer,    intent(in) :: npole,m,mode
      complex_dp, intent(in) :: pole(npole),resid(npole),c
      complex_dp             :: freqintegral_out,fi
      complex_dp             :: r,p,e,lnp
      integer                :: i,sgn
      if(m>=1) then
        freqintegral_out = ( freqintegral_poles(pole,resid,npole,mode,c+1d-4,m-1) -
     &                       freqintegral_poles(pole,resid,npole,mode,c-1d-4,m-1) ) / 2d-4
        return
      else if(m<0) then
        Bug('6th argument out of range.')
      endif
      if(any(pole==0)) Error('one of the poles is zero.')
      fi = 0
      do i = 1,npole
        p   = pole(i) ; if(abs(real(p))<1d-8) p = img*imag(p) ! single pole on imaginary axis: (numerically small but) finite real part would lead to delta function (zero real part corresponds to two
        r   = resid(i)
        e   = c
        sgn = 1
 1      lnp = ln(img*p) ; if(p==e) Error('Argument coincides with a pole.')
        fi  = fi + sgn * r * (lnp-ln(e)) / (e-img*p)
        if(sgn==1.and.mode>0) then
          if(mode==2) then
            p = -conjg(p)
            r = -conjg(r)
          endif
          e   = -e
          sgn = -1
          goto 1
        endif
      enddo
      freqintegral_out = fi / (2*pi)
      end function freqintegral_poles

c ------

c     Returns the mth derivative of
c
c     I(c) = 1/(2*pi*i) * INT de wght(e) INT dw p(iw) / (w-ic+ie)     [ i=img , the integration is along the imaginary axis ]
c
c     with the Pade approximant  p(iw) = SUM(n=poles) resid(n) / (iw-pole(n))
c     and the tetrahedron weight function wght(e) defined between mesh points frq(:nfrq) as cubic polynomials with coefficients pfrq (see tetrahedron_nonlin).
c
c     For a given pole and order m (m=0..3) : 1/(2*pi*i) * INT de e**m [ log(e-c) - log(pole) ] / (c+pole-e)
c
c     The function p(iw) is defined only for positive w.
c     mode=0 : integrate from 0 to infinity
c     mode=1 : integrate from -infinity to infinity; p(-iw) = p(iw)  is assumed
c     mode=2 : integrate from -infinity to infinity; p(-iw) = p*(iw) is assumed
c
c     Remark:
c     mode=1 : INT(-infty,0) p(-iw)  / (w-ic+ie) = INT(0,-infty)   r   / ((w-ic+ie)*(iw+b))  (same formula if b -> -b)
c     mode=2 : INT(-infty,0) p*(-iw) / (w-ic+ie) = INT(0,-infty) (-r*) / ((w-ic+ie)*(iw-b*)) (same formula if b -> b* and r -> -r*)
c
c     See comments for freqintegral_r above.
c
c     Taylor expansion around e2-e1 can be included as algorithm (has a computational overhead):
# define include_taylor
c
c     Routine can be tested:
c # define freqintegral_test
      recursive function freqintegral_poles_nonlin(pole,resid,npole,mode,cin,m,e0,nfrq,frq,pfrq) result(freqintegral_out)
      use global, only: img,pi
      implicit none
      integer,    intent(in) :: npole,m,mode,nfrq
      real_dp,    intent(in) :: frq(nfrq),pfrq(0:3,nfrq-1),e0
      complex_dp, intent(in) :: pole(npole),resid(npole),cin
      complex_dp             :: freqintegral_out,fi
      complex_dp             :: c,r,l0,l1,l2,j1,j2,ec1,ec2,ec1i,ec2i,ec1n,ec2n,lec1,lec2,lec1p,lec2p,q(0:3),w,wn,wi,lw,lc
      complex_dp             :: p,lp,pp,pn,ecp,ecp1,ecp2,ecp1n,ecp2n,b,pc,pcn,decp(0:3),wcp,wcpn,wcpi,dwcp(0:3)
      real_dp                :: e1,e2,sgnpi,de
      logical                :: brcut,notdef(nfrq,2,npole),lser,larr(4),lstab
      integer                :: i,j,n,ifrq,loop,algo
      complex_dp             :: dlecp(nfrq,2,npole),lecp(nfrq,2,npole),lec(nfrq,2),lpp(npole,2),lcp(npole)
      complex_dp             :: dilog
      real                   :: acp(npole,2)
      real,      target      :: ecmin0,ecmax0,ecmin1,ecmax1
      real,      pointer     :: ecmin,ecmax
      real                   :: emax,emin,x,ap,ac,aec1,aec2,aw,xarr(4)
      real,      parameter   :: accur = 1e-14 ! accuracy
      real,      parameter   :: xcut  = 0.5   ! cut criterion for series expansion
      real,      parameter   :: fstab = 500.  ! stability criterion abs(w)<fstab*abs(e2-e1) for Jn(e2-w)-Jn(e1-w)
      integer,   parameter   :: naccur = ceiling(log(accur)/log(xcut)) + 1
      real,      parameter   :: xaccur(naccur)= [ (accur**(1.0/n),n=1,naccur) ]
      logical                :: isnan
# ifdef include_taylor
      real                   :: acpe(nfrq,2,npole),min4,der
# endif      

# ifdef freqintegral_test
      integer sgn,k,ntst,ntst0
      real_dp e
      real_dp, save :: n1=0d0
      complex_dp fj,qq(0:3)
      character(5) line
      integer, save :: nalgo(2,0:4) =0
      ntst = 2500
# endif
      
      if(m>=1) then
        freqintegral_out = ( freqintegral_poles_nonlin(pole,resid,npole,mode,cin+1d-4,m-1,e0,nfrq,frq,pfrq) -
     &                       freqintegral_poles_nonlin(pole,resid,npole,mode,cin-1d-4,m-1,e0,nfrq,frq,pfrq) ) / 2d-4
        return
      else if(m<0) then
        Bug('6th argument out of range.')
      endif
      if(imag(cin)<0)  Bug('Argument has negative imaginary part.')
      if(any(pole==0)) Error('one of the poles is zero.')
      if(npole==0)     Bug('Number of poles is zero.')
      if(nfrq==0)      Bug('Frequency dimension is zero.')

      c  = cin - e0
      aw = 0
      if(imag(c)>0) then
        w  = img * imag(c) ! position of branch cut (loop=2, so c->-c)
        wi = 1/w
        lw = 2*log(w)*img + pi
        aw = abs(w)
      endif
      lc = log(-img*c)
      ac = abs(c)

      do ifrq = 1,nfrq
        lec(ifrq,1) = log(img*(frq(ifrq)-c))
        if(mode>0) lec(ifrq,2) = log(-img*(frq(ifrq)-c))        
      enddo
      do i = 1,npole
        p = pole(i)
        if(abs(real(p))<1d-8) then ! single pole on imaginary axis: (numerically small but) finite real part would lead to delta function (zero real part corresponds to two poles on either sides of the axis)
          p        = img*imag(p)
          lpp(i,1) = log(abs(imag(p)))
          lpp(i,2) = lpp(i,1)
          lcp(i)   = 2 * log(abs(1-imag(c)/imag(p)))
          acp(i,1) = real(abs(c+p))
          acp(i,2) = real(abs(-c+p))
        else
          lpp(i,1) = log(img*p)
          acp(i,1) = real(abs(c+p))
          if(mode==2) then ; p = -conjg(p) ; lpp(i,2) = log(img*p)
          else             ;                 lpp(i,2) = lpp(i,1)
          endif
          lcp(i)   = 2 * log(1-img*imag(c)/p)
          acp(i,2) = real(abs(-c+p))
        endif
      enddo

# ifdef include_taylor
      do i = 1,npole
        p = pole(i)
        do ifrq = 1,nfrq
          e1             = frq(ifrq)
          acpe(ifrq,1,i) = abs(cmplx(c+p-e1))
          if     (mode==1) then ;                 acpe(ifrq,2,i) = abs(cmplx(-c+p+e1))
          else if(mode==2) then ; p = -conjg(p) ; acpe(ifrq,2,i) = abs(cmplx(-c+p+e1))
          endif
        enddo
      enddo
# endif      
      
      notdef  = .true.
      loop    = 1
      fi      = 0
      do ifrq = 1,nfrq-1
        e1     = frq(ifrq)
        e2     = frq(ifrq+1)
        ec1    = e1 - c
        ec2    = e2 - c
        aec1   = abs(ec1)
        aec2   = abs(ec2)
        ecmax0 = max(aec1,aec2)
        ecmax1 = max(ecmax0,aw)
        if(e1*e2<0) then ; ecmin0 = abs(imag(c))   ; emin = 0
        else             ; ecmin0 = min(aec1,aec2) ; emin = min(e1,e2)
        endif
        ecmin1 = min(ecmin0,aw)
        emax   = max(abs(e1),abs(e2))
        lstab  = ac < fstab*abs(e2-e1)
        
        do i = 1,npole
          r     = resid(i)
          p     = pole(i) ; if(abs(real(p))<1d-8) p = img*imag(p) ! single pole on imaginary axis: (numerically small but) finite real part would lead to delta function (zero real part corresponds to two poles on either sides of the axis)
          ap    = abs(p)
 1        lp    = lpp(i,loop)
          brcut = imag(c)<0 .and. real(c)>min(e1,e2) .and. real(c)<max(e1,e2)
          if(brcut) then ; ecmax => ecmax1 ; ecmin => ecmin1
          else           ; ecmax => ecmax0 ; ecmin => ecmin0
          endif
          x    = min(acp(i,loop),ac)
# ifdef include_taylor          
          de   = e2-e1
          der  = abs(real(de))
          min4 = minval( [ acpe(ifrq,loop,i) , acpe(ifrq+1,loop,i) , aec1 , aec2 ] )
          larr = [ emax<x , ap>max(emax,ac).and.lstab , ap<min(emin,ac).and.lstab , der<min4.and..not.brcut ] ! stability criteria for series expansions
# else
          larr = [ emax<x , ap>max(emax,ac).and.lstab , ap<min(emin,ac).and.lstab , .false. ]                 ! stability criteria for series expansions
# endif
          algo = 0
          if(any(larr)) then
# ifdef include_taylor            
            xarr = [ emax/x , ecmax/ap , ap/ecmin , der/min4 ]
# else
            xarr = [ emax/x , ecmax/ap , ap/ecmin , 0. ]
# endif            
            n    = minloc(xarr,1,larr)
            x    = xarr(n)
            if(x<xcut) then
              algo = n
              lser = .true.
              n    = 1 ; do while(x>xaccur(n)) ; n = n + 1 ; enddo ; if(n>naccur) Bug('naccur too small.')
            endif
          endif

c         Series for small e1,e2
          if(algo==1) then
            pc      = 1 + p/c
            pcn     = pc
            ecp1    = e1 / (c+p)
            ecp2    = e2 / (c+p)
            ecp1n   = 1
            ecp2n   = 1
            b       = lc - lp
            q       = 0
            if(.not.brcut) then
              if(loop==2) b = b + sign(1d0,real(ec1)) * pi * img
              do j = 0,3
                ecp1n   = ecp1n * ecp1
                ecp2n   = ecp2n * ecp2
                decp(j) = ( ecp2n - ecp1n ) / (j+1)
              enddo
              do j = 1,n
                q        = q + b * decp
                b        = b - pcn / j
                pcn      = pcn * pc
                ecp1n    = ecp1n * ecp1
                ecp2n    = ecp2n * ecp2
                decp(:2) = decp(1:)
                decp(3)  = ( ecp2n - ecp1n ) / (j+4)
              enddo
            else
              sgnpi = sign(1d0,real(ec1)) * pi
              wcp   = real(c) / (c+p)
              wcpn  = 2
              do j = 0,3
                ecp1n   = ecp1n * ecp1
                ecp2n   = ecp2n * ecp2
                wcpn    = wcpn  * wcp
                decp(j) = ( ecp2n - ecp1n )        / (j+1)
                dwcp(j) = ( ecp2n + ecp1n - wcpn ) / (j+1)
              enddo
              do j = 1,n
                q        = q + b*decp - sgnpi*dwcp*img
                b        = b - pcn / j
                pcn      = pcn * pc
                ecp1n    = ecp1n * ecp1
                ecp2n    = ecp2n * ecp2
                wcpn     = wcpn * wcp
                decp(:2) = decp(1:)
                decp(3)  = ( ecp2n - ecp1n )        / (j+4)
                dwcp(:2) = dwcp(1:)
                dwcp(3)  = ( ecp2n + ecp1n - wcpn ) / (j+4)
              enddo
            endif
            pc   = c + p
            pcn  = pc
            q(1) = pcn * q(1) ; pcn = pcn * pc
            q(2) = pcn * q(2) ; pcn = pcn * pc
            q(3) = pcn * q(3)

c         Series for large p          
          else if(algo==2) then
            pp   = 1/p
            pn   = pp
            ec1n = ec1
            ec2n = ec2
            if(.not.brcut) then
              lec1 = lec(ifrq,loop)
              lec2 = lec(ifrq+1,loop)
              j1   = lec1 - lp
              j2   = lec2 - lp
              q(0) = pn * ( ec2n * (j2-1)     - ec1n * (j1-1)     )     ; ec1n = ec1n * ec1 ; ec2n = ec2n * ec2 ; pn = pn * pp
              q(1) = pn * ( ec2n * (j2-1d0/2) - ec1n * (j1-1d0/2) ) / 2 ; ec1n = ec1n * ec1 ; ec2n = ec2n * ec2 ; pn = pn * pp
              q(2) = pn * ( ec2n * (j2-1d0/3) - ec1n * (j1-1d0/3) ) / 3
              q(3) = 0
              do j = 4,n+3
                ec1n = ec1n * ec1
                ec2n = ec2n * ec2
                pn   = pn * pp
                q(3) = q(3) + pn * ( ec2n * (j2-1d0/j) - ec1n * (j1-1d0/j) ) / j
              enddo
            else
              lec1 = lec(ifrq,3-loop)
              lec2 = lec(ifrq+1,3-loop)
              sgnpi= sign(1d0,real(ec1)) * pi
              j1   = lec1 - lp + sgnpi*img
              j2   = lec2 - lp - sgnpi*img
              wn   = 2*w * sgnpi*img
              q(0) = pn * ( ec2n * (j2-1)     - ec1n * (j1-1)     + wn )     ; ec1n=ec1n*ec1 ; ec2n=ec2n*ec2 ; pn=pn*pp ; wn=wn*w
              q(1) = pn * ( ec2n * (j2-1d0/2) - ec1n * (j1-1d0/2) + wn ) / 2 ; ec1n=ec1n*ec1 ; ec2n=ec2n*ec2 ; pn=pn*pp ; wn=wn*w
              q(2) = pn * ( ec2n * (j2-1d0/3) - ec1n * (j1-1d0/3) + wn ) / 3
              q(3) = 0
              do j = 4,n+3
                ec1n = ec1n * ec1
                ec2n = ec2n * ec2
                wn   = wn * w
                pn   = pn * pp
                q(3) = q(3) + pn * ( ec2n * (j2-1d0/j) - ec1n * (j1-1d0/j) + wn ) / j
              enddo
            endif
            q(2) = q(2) + q(3)
            q(1) = q(1) + q(2)
            q(0) = q(0) + q(1)
            q(1) = p    * q(1)
            q(2) = p**2 * q(2)
            q(3) = p**3 * q(3) + c * ( q(2)*3 + c * ( q(1)*3 + c * q(0) ) )
            q(2) =        q(2) + c * ( q(1)*2 + c *   q(0)                )
            q(1) =        q(1) + c *   q(0)

c         Series for small p
          else if(algo==3) then
            ec1i = 1/ec1
            ec2i = 1/ec2
            if(.not.brcut) then
              lec1 = lec(ifrq,loop)
              lec2 = lec(ifrq+1,loop)
              j1   = lec1 - lp
              j2   = lec2 - lp
              q(0) = ( lec2*(lec2/2-lp)  - lec1*(lec1/2-lp)  ) * p      ; ec1n = ec1        ; ec2n = ec2
              q(1) =   ec2n * (j2-1)     - ec1n * (j1-1)                ; ec1n = ec1n * ec1 ; ec2n = ec2n * ec2
              q(2) = ( ec2n * (j2-1d0/2) - ec1n * (j1-1d0/2) ) / (2*p)  ; ec1n = ec1n * ec1 ; ec2n = ec2n * ec2 ; pn = p**2
              q(3) = ( ec2n * (j2-1d0/3) - ec1n * (j1-1d0/3) ) / (3*pn) ; ec1n = 1          ; ec2n = 1
              do j = 1,n
                ec1n = ec1n * ec1i
                ec2n = ec2n * ec2i
                q(0) = q(0) - pn * ( ec2n * (j2+1d0/j) - ec1n * (j1+1d0/j) ) / j
                pn   = pn   * p
              enddo              
            else
              lec1 = lec(ifrq,3-loop)
              lec2 = lec(ifrq+1,3-loop)
              sgnpi= sign(1d0,real(ec1)) * pi
              wn   = 2*w * sgnpi*img
              j1   = lec1/2 - lp + sgnpi*img
              j2   = lec2/2 - lp - sgnpi*img
              q(0) = ( lec2*j2         - lec1*j1   + sgnpi*lw ) * p      ; ec1n = ec1        ; ec2n = ec2
              j1   = j1 + lec1/2
              j2   = j2 + lec2/2
              q(1) =   ec2n*(j2-1)     - ec1n*(j1-1)     + wn            ; ec1n = ec1n * ec1 ; ec2n = ec2n * ec2 ; wn = wn * w
              q(2) = ( ec2n*(j2-1d0/2) - ec1n*(j1-1d0/2) + wn ) / (2*p)  ; ec1n = ec1n * ec1 ; ec2n = ec2n * ec2 ; wn = wn * w
              pn   = p*p
              q(3) = ( ec2n*(j2-1d0/3) - ec1n*(j1-1d0/3) + wn ) / (3*pn) ; ec1n = 1          ; ec2n = 1          ; wn = 2*sgnpi*img
              do j = 1,n
                ec1n = ec1n * ec1i
                ec2n = ec2n * ec2i
                wn   = wn   * wi
                q(0) = q(0) - pn * ( ec2n * (j2+1d0/j) - ec1n * (j1+1d0/j) + wn ) / j
                pn   = pn   * p
              enddo
            endif
            q(1) = q(1) + q(0)
            q(2) = q(2) + q(1)
            q(3) = q(3) + q(2)
            q(0) = q(0) / p
            q(2) = q(2) * p
            q(3) = p**2 * q(3) + c * ( q(2)*3 + c * ( q(1)*3 + c * q(0) ) )
            q(2) =        q(2) + c * ( q(1)*2 + c *   q(0)                )
            q(1) =        q(1) + c *   q(0)
            q    = -q

c         Taylor wrt e2-e1
          else if(algo==4) then
            lec1     = lec(ifrq,loop)
            wcp      = c+p-e1
            wcpi     = de/wcp
            dwcp(0)  = ( lec1 - lp ) * wcpi
            dwcp(1:) = 0
            q        = 0
            ec1i     = de/ec1
            ec1n     = -1
            do j = 1,n
              dwcp(1) = de * dwcp(1) + e1 * dwcp(0)
              dwcp(2) = de * dwcp(2) + e1 * dwcp(1)
              dwcp(3) = de * dwcp(3) + e1 * dwcp(2)
              q       = q + dwcp / j
              ec1n    = - ec1n * ec1i
              dwcp(3) = dwcp(2)
              dwcp(2) = dwcp(1)
              dwcp(1) = dwcp(0)
              dwcp(0) = wcpi * ( dwcp(0) + ec1n / j )
            enddo
            
c         General algorithm            
          else

            do j = ifrq,ifrq+1
              if(notdef(j,loop,i)) then
                notdef(j,loop,i)  = .false. ; if(j==ifrq) then ; ecp = 1-ec1/p ; else ; ecp = 1-ec2/p ; endif
                lecp (j,loop,i)   = log(ecp)
                dlecp(j,loop,i)   = dilog(ecp)
              endif
            enddo

            lec1p = lecp(ifrq,loop,i)
            lec2p = lecp(ifrq+1,loop,i)
            if(.not.brcut) then
              lec1 = lec(ifrq,loop)
              lec2 = lec(ifrq+1,loop)
              l0   =       ( ( lec2 - 1     ) * ec2      - ( lec1 - 1     ) * ec1      )
              j1   = img * ( ( lec2 - 1d0/2 ) * ec2**2/2 - ( lec1 - 1d0/2 ) * ec1**2/2 )
              j2   =     - ( ( lec2 - 1d0/3 ) * ec2**3/3 - ( lec1 - 1d0/3 ) * ec1**3/3 )
              l1   = c*   l0 -   img*j1
              l2   = c**2*l0 - 2*img*j1*c - j2
              q(0) = dlecp(ifrq,loop,i) - dlecp(ifrq+1,loop,i) + (lp-lec2) * lec2p - (lp-lec1) * lec1p
              q(1) = (c+p)*q(0) - l0 + lp * (e2   -e1   )
              q(2) = (c+p)*q(1) - l1 + lp * (e2**2-e1**2) / 2
              q(3) = (c+p)*q(2) - l2 + lp * (e2**3-e1**3) / 3
            else
              if(lser) then
                lec1 = lec(ifrq,3-loop)
                lec2 = lec(ifrq+1,3-loop)
                l0   =       ( ( lec2 - 1     ) * ec2      - ( lec1 - 1     ) * ec1      )
                j1   = img * ( ( lec2 - 1d0/2 ) * ec2**2/2 - ( lec1 - 1d0/2 ) * ec1**2/2 )
                j2   =     - ( ( lec2 - 1d0/3 ) * ec2**3/3 - ( lec1 - 1d0/3 ) * ec1**3/3 )
                l1   = c*   l0 -   img*j1
                l2   = c**2*l0 - 2*img*j1*c - j2
              else ; l0 = -l0 ; l2 = -l2
              endif
              sgnpi = sign(1d0,real(ec1)) * pi
              q(0)  = dlecp(ifrq,loop,i) - dlecp(ifrq+1,loop,i) + (lp-lec2) * lec2p - (lp-lec1) * lec1p +
     &                sgnpi * ( lec1p + lec2p - lcp(i) ) * img
              q(1)  = (c+p)*q(0) - l0 + lp * (e2   -e1   )     + sgnpi * ( e1    + e2    - 2*real(c)    )     * img
              q(2)  = (c+p)*q(1) - l1 + lp * (e2**2-e1**2) / 2 + sgnpi * ( e1**2 + e2**2 - 2*real(c)**2 ) / 2 * img
              q(3)  = (c+p)*q(2) - l2 + lp * (e2**3-e1**3) / 3 + sgnpi * ( e1**3 + e2**3 - 2*real(c)**3 ) / 3 * img
            endif
            lser = .false.
          endif
          if(loop==2) then ; q(1) = -q(1) ; q(3) = -q(3) ; endif

          fi = fi + r * sum( pfrq(:,ifrq) * q )

# ifdef freqintegral_test
          if(isnan(real(fi))) then
            write(0,'(A)') 'NaN detected, algo=',algo
            Bug('NaN')
          endif
          if(brcut) then ; nalgo(2,algo) = nalgo(2,algo) + 1
          else           ; nalgo(1,algo) = nalgo(1,algo) + 1
          endif
          if(algo>0) n1 = n1+n
          k = sum(nalgo)
          e = n1 / sum(nalgo(:,1:))
c          if     (algo==1) then ; write(*,'(A,8I9,''  %'',5F6.1,I3,F10.5)') 'small e',nalgo,(100.*sum(nalgo(:,j))/k,j=0,3),e,n,x
c          else if(algo==2) then ; write(*,'(A,8I9,''  %'',5F6.1,I3,F10.5)') 'large p',nalgo,(100.*sum(nalgo(:,j))/k,j=0,3),e,n,x
c          else if(algo==3) then ; write(*,'(A,8I9,''  %'',5F6.1,I3,F10.5)') 'small p',nalgo,(100.*sum(nalgo(:,j))/k,j=0,3),e,n,x
c          else                  ; write(*,'(A,8I9,''  %'',5F6.1)')          'general',nalgo,(100.*sum(nalgo(:,j))/k,j=0,3),e
c          endif
          if     (algo==1) then ; write(*,'(A,10I9,''  %'',6F6.1,I3,F10.5)') 'small e',nalgo,(100.*sum(nalgo(:,j))/k,j=0,4),e,n,x
          else if(algo==2) then ; write(*,'(A,10I9,''  %'',6F6.1,I3,F10.5)') 'large p',nalgo,(100.*sum(nalgo(:,j))/k,j=0,4),e,n,x
          else if(algo==3) then ; write(*,'(A,10I9,''  %'',6F6.1,I3,F10.5)') 'small p',nalgo,(100.*sum(nalgo(:,j))/k,j=0,4),e,n,x
          else if(algo==4) then ; write(*,'(A,10I9,''  %'',6F6.1,I3,F10.5)') 'Taylore',nalgo,(100.*sum(nalgo(:,j))/k,j=0,4),e,n,x
          else                  ; write(*,'(A,10I9,''  %'',6F6.1)')          'general',nalgo,(100.*sum(nalgo(:,j))/k,j=0,4),e
          endif

          sgn = 1 ; if(loop==2) sgn = -1
          ntst0 = ntst

          do k = 0,3

            do
            
            fj = 0
            do j = 0,ntst-1
              e  = ( e1*(ntst-(j+0.5d0)) + e2*(j+0.5d0) ) / ntst
              fj = fj + (e2-e1)/ntst * (log(img*(e-c))-lp) / (c-e+p) * (sgn*e)**k
c              fj = fj + (e2-e1)/ntst * (log(img*(e-c))-lp) / (c-e+p) * (sgn*e)**k
c              fj = fj + (e2-e1)/ntst * (log(1-e/c)) / (c-e+p) * (sgn*e)**k              
            enddo
            if(abs(q(k)-fj)/abs(fj)>10d0**(-(5-k)).and.abs(q(k)-fj)>1d-13) then
              lec1p = lecp(ifrq,loop,i)
              lec2p = lecp(ifrq+1,loop,i)
              brcut = imag(c)<0 .and. real(c)>min(e1,e2) .and. real(c)<max(e1,e2)
              if(.not.brcut) then
                lec1  = lec(ifrq,loop)
                lec2  = lec(ifrq+1,loop)
                l0    =       ( ( lec2 - 1     ) * ec2      - ( lec1 - 1     ) * ec1      )
                j1    = img * ( ( lec2 - 1d0/2 ) * ec2**2/2 - ( lec1 - 1d0/2 ) * ec1**2/2 )
                j2    =     - ( ( lec2 - 1d0/3 ) * ec2**3/3 - ( lec1 - 1d0/3 ) * ec1**3/3 )
                l1    = c*   l0 -   img*j1
                l2    = c**2*l0 - 2*img*j1*c - j2
                qq(0) = dilog(1-ec1/p) - dilog(1-ec2/p) + (lp-lec2) * log(1-ec2/p) - (lp-lec1) * log(1-ec1/p)
                qq(1) = (c+p)*qq(0) - l0 + lp * (e2   -e1   )
                qq(2) = (c+p)*qq(1) - l1 + lp * (e2**2-e1**2) / 2
                qq(3) = (c+p)*qq(2) - l2 + lp * (e2**3-e1**3) / 3
              else
                lec1  = lec(ifrq,3-loop)
                lec2  = lec(ifrq+1,3-loop)
                l0    = -     ( ( lec2 - 1     ) * ec2      - ( lec1 - 1     ) * ec1      )
                j1    = img * ( ( lec2 - 1d0/2 ) * ec2**2/2 - ( lec1 - 1d0/2 ) * ec1**2/2 )
                j2    =     - ( ( lec2 - 1d0/3 ) * ec2**3/3 - ( lec1 - 1d0/3 ) * ec1**3/3 )
                l1    =    c*   l0 -   img*j1
                l2    = -( c**2*l0 - 2*img*j1*c - j2 )
                sgnpi = sign(1d0,real(e1-c)) * pi
                qq(0) = dilog(1-ec1/p) - dilog(1-ec2/p) + (lp-lec2) * log(1-ec2/p) - (lp-lec1) * log(1-ec1/p) +
     &            sgnpi * ( log(1-ec1/p) + log(1-ec2/p) - lcp(i) ) * img
                qq(1) = (c+p)*qq(0) + l0 + lp * (e2   -e1   )     + sgnpi * ( e1    + e2    - 2*real(c)    )     * img
                qq(2) = (c+p)*qq(1) - l1 + lp * (e2**2-e1**2) / 2 + sgnpi * ( e1**2 + e2**2 - 2*real(c)**2 ) / 2 * img
                qq(3) = (c+p)*qq(2) + l2 + lp * (e2**3-e1**3) / 3 + sgnpi * ( e1**3 + e2**3 - 2*real(c)**3 ) / 3 * img
              endif
              if(loop==2) then ; qq(1) = -qq(1) ; qq(3) = -qq(3) ; endif
              if(ntst>10000*ntst0) then
                write(*,*) 'e1,e2',e1,e2
                write(*,*) 'c,p',c,p
                write(*,'(4F20.10,I2,L2)') r * sum( pfrq(:,ifrq) * q ),fi,loop,brcut!r,q,loop
                write(*,*) k,ntst
                write(*,'(2F25.15)') fj,q(k),qq(k)
                write(*,'(A,$)') 'c=continue:'
                read(*,'(A)') line ; if(line=='c') then ; ntst = ntst0 ; exit ; endif
              endif
              ntst = ntst + ntst
            else
              ntst = ntst0
              exit
            endif

            enddo

          enddo          
# endif

          if(mode>0) then
            c   = -c
            e1  = -e1 ; ec1 = -ec1
            e2  = -e2 ; ec2 = -ec2
            if(loop==1) then
              loop = 2
              if(mode==2) then
                p = -conjg(p)
                r = -conjg(r)
              endif
              goto 1
            endif
            loop = 1
          endif

        enddo
      enddo

# ifdef freqintegral_test
      write(*,*) 'notdef:',1d0*count(notdef)/size(notdef)      
# endif      

      freqintegral_out = fi / (2*pi*img)

      if(isnan(real(fi))) Bug('Detected NaN.')
      
      end function freqintegral_poles_nonlin

      end module freq_integral

c ------

c     Calculates the integral
c
c                    infinity
c     I(c) = 1/(2*pi)   INT    y(w') dw'
c                        0
c
c     The asymptote y(w)~w^(-m1)+O[w^-m2] is assumed.
c
      function freqintegral0(y,x,n,m1,m2)
      use global, only: pi,metal
      use wrapper
      use, intrinsic :: iso_fortran_env
      implicit none
      real_dp                 :: freqintegral0
      integer,     intent(in) :: n,m1,m2
      real_dp,     intent(in) :: x(n)
      real_dp                 :: x1(n+1),mat(2,2),r
      real_dp,     intent(in) :: y(n)
      real_dp                 :: y1(n+1),aa(n),bb(n),a(0:3,n)
      complex_dp              :: f,tail(2)
      integer                 :: i
      x1(:n)  = x/(1+x)
      x1(n+1) = 1
      y1(:n)  = y
      y1(n+1) = 0
      if(metal) then ; i = 4
      else           ; i = 1
      endif
      call cubicspline(aa,bb,y1,x1,n+1,0d0,i,0d0,1)
      do i = 1,n
        call cubiccoeff(a(:,i),aa(i),bb(i),x1(i),y1(i))
      enddo
      f = 0
      ! intervals (x(j)->x(j+1))
      do i = 1,n-1
        f = f + a(0,i) * ( x(i+1) - x(i) )
     &    + a(1,i) * ( x(i+1)                                         -   log(1+x(i+1)) -
     &               ( x(i)                                           -   log(1+x(i))   ) )
     &    + a(2,i) * ( x(i+1) -                        1 / (1+x(i+1)) - 2*log(1+x(i+1)) -
     &               ( x(i)   -                        1 / (1+x(i))   - 2*log(1+x(i))   ) )
     &    + a(3,i) * ( x(i+1) - ( 3 - 0.5d0/(1+x(i+1)) ) / (1+x(i+1)) - 3*log(1+x(i+1)) -
     &               ( x(i)   - ( 3 - 0.5d0/(1+x(i))   ) / (1+x(i))   - 3*log(1+x(i))   ) )
      enddo
      ! tail contribution (freq(n)->infinity)
      r        = x(n)
      mat(1,1) =     r**(-m1)   ; mat(1,2) =     r**(-m2)
      mat(2,1) = -m1*r**(-m1)/r ; mat(2,2) = -m2*r**(-m2)/r
      mat      = invert(mat)
      tail(1)  = y(n)
      tail(2)  = (a(1,n)+a(2,n)*2*x1(n)+a(3,n)*3*x1(n)**2) / (1+r)**2
      tail     = matmul ( mat,tail )
      f        = f + tail(1) / r**(m1-1) / (m1-1) + tail(2) / r**(m2-1) / (m2-1)
c      mat(1,1) =    r**(-4) ; mat(1,2) =    r**(-6)
c      mat(2,1) = -4*r**(-5) ; mat(2,2) = -6*r**(-7)
c      mat      = invert(mat)
c      tail(1)  = y(n)
c      tail(2)  = (a(1,n)+a(2,n)*2*x1(n)+a(3,n)*3*x1(n)**2) / (1+r)**2
c      tail     = matmul ( mat,tail )
c      f        = f + tail(1) / r**3 / 3 + tail(2) / r**5 / 5
      freqintegral0 = f / (2*pi)
      end

c ------

c     Returns the mth derivative of
c
c                      infinity
c     I(c) = 1/(2pi) *   INT   SUM(i=poles) resid(i) / (img*w-pole(i)) dw.
c                         0
c
      function freqintegral0_pade(y,x,n)
      use global, only: img,pi,smooth
      use, intrinsic :: iso_fortran_env
      implicit none
      integer,    intent(in) :: n
      real_dp,    intent(in) :: x(n)
      complex_dp, intent(in) :: y(n)
      complex_dp             :: freqintegral0_pade
      complex_dp             :: cpade(n+1),pole((n+1)/2*smooth(2)),resid((n+1)/2*smooth(2))
      integer                :: i,npole
      call pade_init(cpade,img*x,y,n,-1)
      call pade_poles(pole,resid,npole,img*x,y,cpade,n,smooth(2),-1,.true.)
      freqintegral0_pade = sum( [ ( resid(i) * ( img*pi/2 - log(-pole(i)) ), i=1,npole) ] ) / (2*pi*img)
      end

c -------------------------------------


c backup
# if 0
      function freqintegral1(y,x,n,c,m)
      use global, only: pi,img
      use, intrinsic :: iso_fortran_env
      implicit none
      complex_dp              :: freqintegral1
      complex_dp,  intent(in) :: c
      integer,     intent(in) :: n,m
      real_dp,     intent(in) :: x(n)
      real_dp                 :: x1(n+1),xa(n,7)
      MCOMPLEX_dp, intent(in) :: y(n)
      MCOMPLEX_dp             :: y1(n+1),aa(n),bb(n),a(0:3,n)
      complex_dp              :: cc,c1,f,integ(0:3,0:m+1),integ1(0:3)
      integer                 :: i,j,sgn
      x1(:n)  = x/(1+x)
      x1(n+1) = 1
      y1(:n)  = y
      y1(n+1) = 0
# ifdef INV
      call cubicspline  (aa,bb,y1,x1,n+1,0d0,1,0d0,1)
# else
      call cubicspline_c(aa,bb,y1,x1,n+1,(0d0,0d0),1,(0d0,0d0),1)
# endif
      do i = 1,n
        call cubiccoeff(a(:,i),aa(i),bb(i),x1(i),y1(i))
      enddo
      sgn = 1
      f   = 0
      cc  = c
 1    c1  = cc - 1
      if(abs(c1)<1d-3) then
        do i = 1,n
          do j = 1,7
            xa(i,j) = (1+x(i))**(-j)
          enddo
        enddo
      endif
      ! intervals (x(j)->x(j+1))
      do i = 1,n-1
        if(abs(c1)<1d-3) then
          if(m==0) then
            integ(0,1) = ln(cc+x(i+1)) - ln(cc+x(i))
            integ(1,1) = - (xa(i+1,1)-xa(i,1))   + (xa(i+1,2)-xa(i,2))*c1/2 - (xa(i+1,3)-xa(i,3))*c1**2/3
            integ(2,1) = - (xa(i+1,2)-xa(i,2))/2 + (xa(i+1,3)-xa(i,3))*c1/3 - (xa(i+1,4)-xa(i,4))*c1**2/4
            integ(3,1) = - (xa(i+1,3)-xa(i,3))/3 + (xa(i+1,4)-xa(i,4))*c1/4 - (xa(i+1,5)-xa(i,5))*c1**2/5
          else if(m==1) then
            integ(0,2) = - (cc+x(i+1))**(-1) + (cc+x(i))**(-1)
            integ(1,2) = - (xa(i+1,2)-xa(i,2))/2 + (xa(i+1,3)-xa(i,3))*c1*2/3 - (xa(i+1,4)-xa(i,4))*c1**2*3/4
            integ(2,2) = - (xa(i+1,3)-xa(i,3))/3 + (xa(i+1,4)-xa(i,4))*c1/2   - (xa(i+1,5)-xa(i,5))*c1**2*3/5
            integ(3,2) = - (xa(i+1,4)-xa(i,4))/4 + (xa(i+1,5)-xa(i,5))*c1*2/5 - (xa(i+1,6)-xa(i,6))*c1**2/2
          else
            integ(0,3) = ( (cc+x(i+1))**(-2) - (cc+x(i))**(-2) )/2
            integ(1,3) = - (xa(i+1,3)-xa(i,3))/3 + (xa(i+1,4)-xa(i,4))*c1*3/4 - (xa(i+1,5)-xa(i,5))*c1**2*6/5
            integ(2,3) = - (xa(i+1,4)-xa(i,4))/4 + (xa(i+1,5)-xa(i,5))*c1*3/5 - (xa(i+1,6)-xa(i,6))*c1**2
            integ(3,3) = - (xa(i+1,5)-xa(i,5))/5 + (xa(i+1,6)-xa(i,6))*c1/2   - (xa(i+1,7)-xa(i,7))*c1**2*6/7
          endif
        else
          integ(1,0) =   log(1+x(i+1))       - log(1+x(i))
          integ(2,0) =    - (1+x(i+1))**(-1) +    (1+x(i))**(-1)
          integ(3,0) =  ( - (1+x(i+1))**(-2) +    (1+x(i))**(-2) ) / 2
          integ(0,1) =   ln(cc+x(i+1))       - ln(cc+x(i))
          if(m>=1) integ(0,2) =   - (cc+x(i+1))**(-1) +  (cc+x(i))**(-1)
          if(m==2) integ(0,3) = ( - (cc+x(i+1))**(-2) +  (cc+x(i))**(-2) ) / 2
          integ(1,1) = ( integ(1,0)-integ(0,1) ) / c1
          integ(2,1) = ( integ(2,0)-integ(1,1) ) / c1
          integ(3,1) = ( integ(3,0)-integ(2,1) ) / c1
          if(m>=1) then
            integ(1,2) = ( integ(1,1)-integ(0,2) ) / c1
            integ(2,2) = ( integ(2,1)-integ(1,2) ) / c1
            integ(3,2) = ( integ(3,1)-integ(2,2) ) / c1
          endif
          if(m==2) then
            integ(1,3) = ( integ(1,2)-integ(0,3) ) / c1
            integ(2,3) = ( integ(2,2)-integ(1,3) ) / c1
            integ(3,3) = ( integ(3,2)-integ(2,3) ) / c1
          endif
        endif
        integ1(0) = integ(0,1+m)
        integ1(1) = integ(0,1+m) -   integ(1,1+m)
        integ1(2) = integ(0,1+m) - 2*integ(1,1+m) +   integ(2,1+m)
        integ1(3) = integ(0,1+m) - 3*integ(1,1+m) + 3*integ(2,1+m) - integ(3,1+m)
        f         = f + sgn**(m+1) * dot_product(a(:,i),integ1)
      enddo
      ! tail contribution (freq(n)->infinity)

      ! add negative axis
      if(sgn==1) then
        cc  = -cc
        sgn = -1
        goto 1
      endif
      if     (m==0) then ; freqintegral1 =   f / (2*pi*img)
      else if(m==1) then ; freqintegral1 =  -f / (2*pi*img)
      else               ; freqintegral1 = 2*f / (2*pi*img)
      endif
      end function freqintegral1

c     OLD ROUTINE (BACKUP):

c
c     Calculates the mth derivative of the integral
c
c                      infinity
c     I(c) = 1/(2*pi*i)   INT    y(w') / (w'+c) dw'
c                     -infinity
c
c     (If one of the mesh points coincides with c, instead c+dc and c-dc are used, and the results are averaged.)
c
      function freqintegral0(y,x,n,c,m)
      use global, only: pi,img
      use wrapper
      use, intrinsic :: iso_fortran_env
      implicit none
      complex_dp              :: freqintegral0
      complex_dp,  intent(in) :: c
      integer,     intent(in) :: n,m
      real_dp,     intent(in) :: x(n)
      MCOMPLEX_dp, intent(in) :: y(n)
      MCOMPLEX_dp             :: aa(n),bb(n),y1(n+1),a(4)
      real_dp                 :: dc,ac,ac1,ac1_ ,xx,dx,x1(n+1),amat2(2,2),xpower(n,0:9)
      real_dp                 :: d0,d1,d2,d3,d4,d5,d6,d7,d8,d9
      complex_dp              :: result,cdum,d,lnx(n,-1:1) , hlp
      complex_dp              :: e0,e0_,e1,e1_,e2,e2_
      complex_dp              :: cc,cc2,cc3,cc4,c1,c1_,c12,c12_,c13,c13_,c14,c14_,c15,c15_
      complex_dp              :: p0,p1,p2,p3
      integer                 :: i,j,idum
      logical                 :: tay
      freqintegral0 = 0
      result       = 0
      x1(:n)       = x/(1+x)
      x1(n+1)      = 1d0
      y1(:n)       = y
      y1(n+1)      = 0d0
      do i = 1,n
        xx          = 1d0/(1+x(i))
        xpower(i,0) = log(1+x(i))
        xpower(i,1) = xx ; xx = xx * xx
        xpower(i,2) = xx ; xx = xx * xx
        xpower(i,3) = xx ; xx = xx * xx
        xpower(i,4) = xx ; xx = xx * xx
        xpower(i,5) = xx ; xx = xx * xx
        xpower(i,6) = xx ; xx = xx * xx
        xpower(i,7) = xx ; xx = xx * xx
        xpower(i,8) = xx ; xx = xx * xx
        xpower(i,9) = xx
      enddo
      if( m<0.or.m>2) Error('Power out of range (bug).')
      if((m==0.or.m==2).and.c==0) return
      if     (m==0) then ; dc = 1d-12
      else if(m==1) then ; dc = 1d-6
      else               ; dc = 1d-3
      endif
      if(any(abs(c-x)<dc)) then ; cc = c + dc ! If c is identical to any mesh point, use c+/-dc and average.
      else                      ; cc = c
      endif
      hlp = 0 !; rewind(123) ; if(all(abs(y)<1d-6)) return
 1    do i = 1,n
        lnx(i, 1) = ln(x(i)+cc)
        lnx(i,-1) = ln(x(i)-cc)
      enddo
      cc2  = cc  * cc
      cc3  = cc2 * cc
      cc4  = cc2 * cc2
      c1   =  cc - 1
      c1_  = -cc - 1
      ac   = abs(real(cc))  + abs(imag(cc))
      ac1  = abs(real(c1))  + abs(imag(c1))
      ac1_ = abs(real(c1_)) + abs(imag(c1_))
      tay  = .false.
      if     (ac1 <1d-3) then ; c1  = 1d30 ; tay = .true. ! if we come too close to 1 (or -1), we use Taylor expansions
      else if(ac1_<1d-3) then ; c1_ = 1d30 ; tay = .true. !
      endif
      c12  = c1 **2 ; c13  = c12 *c1
      c12_ = c1_**2 ; c13_ = c12_*c1_
c      if(min(abs(c1),abs(c1_))<1d-6) Error('argument c too close to 1. Modify MESH parameters slightly.')
      if(m==0) then
        p0  = -(3*cc2-3*cc+1)/c13 + (3*cc2+3*cc+1)/c13_
        p1  = -(3*cc-2)      /c12 + (-3*cc-2)     /c12_
        p2  =               1/c1  -              1/c1_
        p3  = -(2*cc-1)      /c12 + (-2*cc-1)     /c12_
      else if(m==1) then
        c14 = c13 *c1 ; c14_ = c13_*c1_
        p0  = 1/c14 + 1/c14_
        p1  = 1/c13 - 1/c13_
        p2  = 1/c12 + 1/c12_
        p3  = (3*cc-1)/c13 + (-3*cc-1)/c13_
      else
        c14 = c13 *c1 ; c14_ = c13_*c1_
        c15 = c14 *c1 ; c15_ = c14_*c1_
        p0  = 1/c13 - 1/c13_
        p1  = 1/c14 + 1/c14_
      endif
# ifdef INV
      call cubicspline  (aa,bb,y1,x1,n+1,0d0,1,0d0,1)
# else
      call cubicspline_c(aa,bb,y1,x1,n+1,(0d0,0d0),1,(0d0,0d0),1)
# endif
      do i = n,1,-1
        call cubiccoeff0(a,aa(i),bb(i),x1(i),y1(i))
        if(i==n) then ! tail contribution (find expansion a(1)/x**2 + a(2)/x**4 and integrate)
          amat2(1,1) =    x(n)**(-2) ; amat2(1,2) =    x(n)**(-4)
          amat2(2,1) = -2*x(n)**(-3) ; amat2(2,2) = -4*x(n)**(-5)
          a(:2)      = matmul(invert(amat2),[y(n),(3*x1(n)**2*a(1)+2*x1(n)*a(2)+a(3))/(1+x(n))**2])
          if(ac<1d-1) then ! cc approx. 0
            cdum = cc/x(n)
            if     (m==0) then
              result = result - 2 * a(1) * ( cdum/3 + cdum**3/5 + cdum**5/7 + cdum**7/9  + cdum**9/11 ) / x(n)**2
     &                        - 2 * a(2) * ( cdum/5 + cdum**3/7 + cdum**5/9 + cdum**7/11 + cdum**9/13 ) / x(n)**4
            else if(m==1) then
              result = result + 2 * a(1) * ( 1d0/3 + 3*cdum**2/5 + 5*cdum**4/7 + 7*cdum**6/9  + 9*cdum**8/11 ) / x(n)**3
     &                        + 2 * a(2) * ( 1d0/5 + 3*cdum**2/7 + 5*cdum**4/9 + 7*cdum**6/11 + 9*cdum**8/13 ) / x(n)**5
            else
              result = result - 2 * a(1) * ( 3*cdum/5 + 10*cdum**3/7 + 21*cdum**5/9  + 36*cdum**7/11 + 55*cdum**9/13 ) / x(n)**4
     &                        - 2 * a(2) * ( 3*cdum/7 + 10*cdum**3/9 + 21*cdum**5/11 + 36*cdum**7/13 + 55*cdum**9/15 ) / x(n)**6
            endif
          else                ! otherwise
            if     (m==0) then
              result = result + (  a(1)+   a(2)/cc2) * (lnx(n,-1)-lnx(n,1)+2*cc/x(n)) / cc2
     &                        + 2*a(2)/(cc*x(n)**3)/3
            else if(m==1) then
              result = result + (2*a(1)+ 4*a(2)/cc2) * (lnx(n,-1)-lnx(n,1)) / cc3
     &                        + 2*(a(1)+   a(2)*(3/cc2+1/x(n)**2/3)) / (cc2*x(n))
     &                        +   (a(1)+   a(2)/cc2)                 /  cc2 * ((x(n)+cc)**(-1)+(x(n)-cc)**(-1))
            else
              result = result + (3*a(1)+10*a(2)/cc2) * (lnx(n,-1)-lnx(n,1)) / cc4
     &                        + 2*(a(1)+   a(2)*(6/cc2+1/x(n)**2/3)) / (cc3*x(n))
     &                        + (2*a(1)+ 4*a(2)/cc2)                 /  cc3 * ((x(n)+cc)**(-1)+(x(n)-cc)**(-1))
     &                        + (  a(1)+   a(2)/cc2)                 /  cc2 * ((x(n)+cc)**(-2)-(x(n)-cc)**(-2)) / 2
            endif
          endif
# if 0
          ! testing
          hlp = 0
          xx  = x(n)-0.005d0
          do idum = 0,1000
            xx = xx + 0.01d0
            hlp = hlp + (1/(xx+cc)**(m+1)-(-1)**m*1/(xx-cc)**(m+1))*(a(1)/xx**2+a(2)/xx**4)*0.01d0
            write(123,*) xx,real(a(1)/xx**2+a(2)/xx**4)
          enddo
c          write(*,*) result
c          write(*,*) hlp
c          write(*,*)
# endif
        else ! contribution of interval [x(i),x(i+1)]
          idum = max(1,i-1)
          d0   = xpower(i+1,0) - xpower(i,0)
          d1   = xpower(i+1,1) - xpower(i,1)
          d2   = xpower(i+1,2) - xpower(i,2)
          e0   = lnx(i+1, 1) - lnx(i, 1)
          e0_  = lnx(i+1,-1) - lnx(i,-1)
          if(m/=0) then
            if(m>=1) then ; e1 = (x(i+1)+cc)**(-1) - (x(i)+cc)**(-1) ; e1_ = (x(i+1)-cc)**(-1) - (x(i)-cc)**(-1) ; endif
            if(m==2) then ; e2 = (x(i+1)+cc)**(-2) - (x(i)+cc)**(-2) ; e2_ = (x(i+1)-cc)**(-2) - (x(i)-cc)**(-2) ; endif
          endif
          if(m==0) then
            result = result
     &        + a(1) * ( cc3*(e0/c13+e0_/c13_) + p0*d0 + p1*d1 + p2*d2/2 )
     &        + a(2) * ( cc2*(e0/c12-e0_/c12_) + p3*d0 - p2*d1 )
     &        + a(3) * ( cc *(e0/c1 +e0_/c1_ ) - p2*d0 )
     &        + a(4) *       (e0-e0_)
          else if(m==1) then
            result = result
     &        + a(1) * ( 3*cc2*(e0/c14+e0_/c14_-d0*p0) - cc3*(e1/c13-e1_/c13_) + p2*d2/2 - p3*d1 )
     &        + a(2) * ( 2*cc *(e0/c13-e0_/c13_-d0*p1) - cc2*(e1/c12+e1_/c12_) - p2*d1 )
     &        + a(3) * (        e0/c12+e0_/c12_-d0*p2  - cc *(e1/c1 -e1_/c1_ ) )
     &        + a(4) * (-e1-e1_)
          else
            result = result
     &        + a(1) * ( (3*cc2+3*cc)*(e0-d0)/c15-(3*cc2-3*cc)*(e0_-d0)/c15_ - 3*cc2*(e1/c14-e1_/c14_) - cc3*(e2/c13+e2_/c13_)/2
     &        - 3*cc*p1*d1 + p0*d2/2 )
     &        + a(2) * (     (2*cc+1)*(e0-d0)/c14+    (2*cc-1)*(e0_-d0)/c14_ - 2*cc *(e1/c13+e1_/c13_) - cc2*(e2/c12-e2_/c12_)/2
     &        - p0*d1 )
     &        + a(3) * (              (e0-d0)/c13-             (e0_-d0)/c13_ -       (e1/c12-e1_/c12_) - cc *(e2/c1 +e2_/c1_ )/2 )
     &        + a(4) * (-e2+e2_)/2
          endif
c          write(*,*) result
          if(tay) then ! Taylor expansion for the case cc-1 or -cc-1 approx. 0
            if     (c1 ==1d30) then ; d =  cc - 1 ; idum =  1
            else if(c1_==1d30) then ; d = -cc - 1 ; idum = -1
            else                      ; Error('bug.')
            endif
            d3 = xpower(i+1,3) - xpower(i,3)
            d4 = xpower(i+1,4) - xpower(i,4)
            d5 = xpower(i+1,5) - xpower(i,5)
            d6 = xpower(i+1,6) - xpower(i,6)
            d7 = xpower(i+1,7) - xpower(i,7)
            d8 = xpower(i+1,8) - xpower(i,8)
            d9 = xpower(i+1,9) - xpower(i,9)
            if(m==0) then
              result = result + idum * (
     &          + a(1) * (d0+3*d1-3*d2/2+d3/3 - d*(-d1+3*d2/2-d3+d4/4       - d*(-d2/2+d3-3*d4/4+d5/5
     &                                        - d*(-d3/3+3*d4/4-3*d5/5+d6/6 - d*(-d4/4+3*d5/5-d6/2+d7/7
     &                                        - d*(-d5/5+d6/2-3*d7/7+d8/8   - d*(-d6/6+3*d7/7-3*d8/8+d9/9)))))))
     &          + a(2) * (d0+2*d1-d2/2        - d*(-d1+d2-d3/3              - d*(-d2/2+2*d3/3-d4/4
     &                                        - d*(-d3/3+d4/2-d5/5          - d*(-d4/4+2*d5/5-d6/6
     &                                        - d*(-d5/5+d6/3-d7/7          - d*(-d6/6+2*d7/7-d8/8)))))))
     &          + a(3) * (d0+d1               - d*(-d1+d2/2                 - d*(-d2/2+d3/3
     &                                        - d*(-d3/3+d4/4               - d*(-d4/4+d5/5
     &                                        - d*(-d5/5+d6/6               - d*(-d6/6+d7/7))))))))
            else if(m==1) then
              result = result
     &          + a(1) * (-d1+3*d2/2-d3+d4/4 - 2*d   *(-d2/2+d3-3*d4/4+d5/5)   + 3*d**2*(-d3/3+3*d4/4-3*d5/5+d6/6)
     &                                       - 4*d**3*(-d4/4+3*d5/5-d6/2+d7/7) + 5*d**4*(-d5/5+d6/2-3*d7/7+d8/8)
     &                                       - 6*d**5*(-d6/6+3*d7/7-3*d8/8+d9/9))
     &          + a(2) * (-d1+d2-d3/3        - 2*d*   (-d2/2+2*d3/3-d4/4)      + 3*d**2*(-d3/3+d4/2-d5/5)
     &                                       - 4*d**3*(-d4/4+2*d5/5-d6/6)      + 5*d**4*(-d5/5+d6/3-d7/7)
     &                                       - 6*d**5*(-d6/6+2*d7/7-d8/8))
     &          + a(3) * (-d1+d2/2           - 2*d   *(-d2/2+d3/3)             + 3*d**2*(-d3/3+d4/4)
     &                                       - 4*d**3*(-d4/4+d5/5)             + 5*d**4*(-d5/5+d6/6)
     &                                       - 6*d**5*(-d6/6+d7/7))
            else
              result = result + idum * (
     &          + a(1) * (-d2/2+d3-3*d4/4+d5/5 -  3*d   *(-d3/3+3*d4/4-3*d5/5+d6/6) +  6*d**2*(-d4/4+3*d5/5-d6/2+d7/7)
     &                                         - 10*d**3*(-d5/5+d6/2-3*d7/7+d8/8)   + 15*d**4*(-d6/6+3*d7/7-3*d8/8+d9/9))
     &          + a(2) * (-d2/2+2*d3/3-d4/4    -  3*d   *(-d3/3+d4/2-d5/5)          +  6*d**2*(-d4/4+2*d5/5-d6/6)
     &                                         - 10*d**3*(-d5/5+d6/3-d7/7)          + 15*d**4*(-d6/6+2*d7/7-d8/8))
     &          + a(3) * (-d2/2+d3/3           -  3*d   *(-d3/3+d4/4)               +  6*d**2*(-d4/4+d5/5)
     &                                         - 10*d**3*(-d5/5+d6/6)               + 15*d**4*(-d6/6+d7/7)))
            endif
          endif
# if 0
          ! testing
          dx = (x(i+1)-x(i))/1000
          xx = x(i) - dx/2
          do j = 1,1000
            xx = xx + dx
            hlp = hlp + (-(-1)**m/(xx-cc)**(m+1)) *dx *
     &        (a(1)*(xx/(1+xx))**3+a(2)*(xx/(1+xx))**2+a(3)*(xx/(1+xx))+a(4))+
     &        a(4)*dx /(xx+cc)**(m+1)
            if(j==1000) then
              write(123,'(3F20.10)') xx,real(a(1)*(xx/(1+xx))**3+a(2)*(xx/(1+xx))**2+a(3)*(xx/(1+xx))+a(4)),dble(y(i+1))
            else
              write(123,'(3F20.10)') xx,real(a(1)*(xx/(1+xx))**3+a(2)*(xx/(1+xx))**2+a(3)*(xx/(1+xx))+a(4))
            endif
          enddo
c          write(*,*) hlp
          xx = x(i) - dx/2
          do j = 1,1000
            xx = xx + dx
            hlp = hlp + 1/(xx+cc)**(m+1) * dx *
     &        (a(1)*(xx/(1+xx))**3+a(2)*(xx/(1+xx))**2+a(3)*(xx/(1+xx)))!+a(4))
          enddo
c          write(*,*) hlp
c          write(*,*)
# endif
        endif
      enddo

      if     (real(cc)>real(c)) then ; cc = c - dc ; goto 1
      else if(real(cc)<real(c)) then ; result = result / 2 !; hlp = hlp / 2
      endif
      if     (m==0) then ; freqintegral0 =   result / (2*pi*img)
      else if(m==1) then ; freqintegral0 =  -result / (2*pi*img)
      else               ; freqintegral0 = 2*result / (2*pi*img)
      endif

      contains
      function ln(z)
      use, intrinsic :: iso_fortran_env
      implicit none
      complex_dp             :: ln
      complex_dp, intent(in) :: z
      if   (abs(z)<1d-12) then ; ln = 0
      else if(imag(z)==0) then ; ln = log(abs(z))
      else                     ; ln = log(z)
      endif
      end function ln
      end

      subroutine cubiccoeff0(c,a,b,x,y)
      use, intrinsic :: iso_fortran_env
      implicit none
      MCOMPLEX_dp, intent(out) :: c(4)
      real_dp,     intent(in)  :: x(2)
      MCOMPLEX_dp, intent(in)  :: a,b,y(2)
      real_dp                  :: dx
      dx   = ( x(2) - x(1) )**3
      c(1) = ( -           y(1) +           y(2) -                      a +                      b ) / dx
      c(2) = (   3*x(2)*   y(1) - 3*x(1)*   y(2) +      (2*x(1)+x(2)) * a -      (2*x(2)+x(1)) * b ) / dx
      c(3) = ( - 3*x(2)**2*y(1) + 3*x(1)**2*y(2) - x(1)*(2*x(2)+x(1)) * a + x(2)*(2*x(1)+x(2)) * b ) / dx
      c(4) = (     x(2)**3*y(1) -   x(1)**3*y(2) + x(1)**2*x(2)       * a - x(2)**2*x(1)       * b ) / dx
      end
# endif
