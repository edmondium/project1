c--------------------------------------------------------------------------------
c Copyright (c) 2021 Peter Grünberg Institut, Forschungszentrum Jülich, Germany
c This file is part of SPEX and available as free software under the conditions
c of the MIT license as expressed in the LICENSE file in more detail.
c--------------------------------------------------------------------------------

c     Integrates function f numerically (Lagrange and Simpson integration).
c     The one-dimensional grid is defined in "grid" (see global.f and below).
c
c     normally  : outward integration from origin (0) to r1*exp(n*h) (normally but not necessarily = grid%radius)
c                 integral from 0 to r1 approximated by leading term in power series expansion  c * r**x  with x > -1
c
c     if  h < 0 : inward integration from r1 to r1*exp(n*h) with r1*exp(n*h) < r1
c
c     Function intgrf provides a much faster integration but must be
c     initialized with intgrf_init (see below).

# include "cppmacro.h"

      function intgr(f,grid)

      use global, only: gridtype
      use util,   only: chr

      use, intrinsic :: iso_fortran_env
      implicit none
      real_dp              :: intgr
      real_dp, intent(in)  :: f(*)
      type(gridtype)       :: grid
      real_dp              :: r(7),fr(7),r1,h,dr,rdum,x
      integer              :: n,nstep,n0=6,i,j
      integer, parameter   :: simpson(7)   = [41,216,27,272,27,216,41]
      real_dp, parameter   :: lagrange(7,6)= reshape(
     &                        [ 19087d0,65112d0,-46461d0, 37504d0,-20211d0, 6312d0, -863d0,
     &                           -863d0,25128d0, 46989d0,-16256d0,  7299d0,-2088d0,  271d0,
     &                            271d0,-2760d0, 30819d0, 37504d0, -6771d0, 1608d0, -191d0,
     &                           -191d0, 1608d0, -6771d0, 37504d0, 30819d0,-2760d0,  271d0,
     &                            271d0,-2088d0,  7299d0,-16256d0, 46989d0,25128d0, -863d0,
     &                           -863d0, 6312d0,-20211d0, 37504d0,-46461d0,65112d0,19087d0 ], ! The last row is actually never used.
     &                        [7,6])

      n  = grid%number
      r1 = grid%first
      h  = grid%increment

      ! integral from 0 to r1 approximated by leading term in power series expansion
      ! (if the increment is negative (h<0), then we obviously perform inward integration
      ! and no approximation is needed)
c      if (h>0.and.f(1)*f(2)>1d-10) then
c        if(f(2)==f(1)) then
c          intgr = r1*f(1)
c        else
c          x     = (f(3)-f(2))/(f(2)-f(1))
c          a     = (f(2)-x*f(1)) / (1-x)
c          x     = log(x)/h
c          if(x<0) then
c            if(x>-1) write(6,'(A,ES9.1)') 'intgr: Warning! Negative exponent x in extrapolation a+c*r**x:',x
c            if(x<=-1) write(6,'(A,ES9.1)') 'intgr: Negative exponent x in extrapolation a+c*r**x:',x
c            if(x<=-1) Error('Negative exponent x in extrapolation a+c*r**x')
c          endif
c          intgr = r1*(f(1)+x*a) / (x+1)
c        endif
c      else
c        intgr = 0
c      endif
      if (h>0.and.f(1)*f(2)>1d-20) then
        x = log ( f(2)/f(1) ) / h
        if(x<-1) then ; Error('Exponent x in extrapolation c*r**x below -1: '//Chf(x,'ES9.1')) ; endif
        if(x<-0.5)      Warn('Large negative exponent x in extrapolation c*r**x: '//Chf(x,'ES9.1'))
        intgr = r1*f(1) / (x+1)
      else
        intgr = 0
      endif

      nstep = (n-1)/6
      n0    = n-6*nstep
      dr    = exp(h)

      ! Lagrange integration from r(1) to r(n0)
      r(1)=r1
      if(n0>1) then
        do i=2,7
          r(i) = r(i-1)*dr
        enddo
        fr = f(:7) * r
c        do i=n+1,7                                      ! If n<7 we still need fr(1:7) for the Lagrange integration. Therefore, points
c          fr(i) = fr(i-1) * (f(n)/f(n-1)) * r(i)/r(i-1) ! fr(n+1:7) are added using the leading term in the power series expansion.
c        enddo                                           !
        intgr = intgr + h/60480 *
     &                  sum ( [(dot_product(lagrange(:,i),fr),i=1,n0-1)] )
        r(1)  = r(n0)
      endif

      ! Simpson integration from r(n0) to r(n)
      rdum = 0
      do i = 1,nstep
        do j = 2,7
          r(j) = r(j-1)*dr
        enddo
        fr   = f(n0:n0+6) * r
        rdum = rdum + dot_product(simpson,fr)
        n0   = n0 + 6
        r(1) = r(7)
      enddo
      intgr = intgr + h/140 * rdum

      end

c     -----------------------------

c     Initializes fast numerical integration intgrf (see below).

      subroutine intgrf_init

      use global, only: ntype, grid, gridf

      use, intrinsic :: iso_fortran_env
      implicit none
      real_dp              :: r(7),r1,h,dr
      integer              :: n,nstep,n0=6,i,j,itype
      integer, parameter   :: simpson(7)   = [41,216,27,272,27,216,41]
      real_dp, parameter   :: lagrange(7,6)= reshape(
     &                        [ 19087d0,65112d0,-46461d0, 37504d0,-20211d0, 6312d0, -863d0,
     &                           -863d0,25128d0, 46989d0,-16256d0,  7299d0,-2088d0,  271d0,
     &                            271d0,-2760d0, 30819d0, 37504d0, -6771d0, 1608d0, -191d0,
     &                           -191d0, 1608d0, -6771d0, 37504d0, 30819d0,-2760d0,  271d0,
     &                            271d0,-2088d0,  7299d0,-16256d0, 46989d0,25128d0, -863d0,
     &                           -863d0, 6312d0,-20211d0, 37504d0,-46461d0,65112d0,19087d0 ], ! The last row is actually never used.
     &                        [7,6])

      n = maxval(grid(:)%number)
      allocate ( gridf(n,ntype) )

      gridf = 0

      do itype = 1,ntype

        n  = grid(itype)%number
        r1 = grid(itype)%first
        h  = grid(itype)%increment

        nstep = (n-1)/6
        n0    = n-6*nstep
        dr    = exp(h)

        ! Calculate Lagrange-integration coefficients from r(1) to r(n0)
        r(1)=r1
        if(n0>1) then
          do i=2,7
            r(i) = r(i-1)*dr
          enddo
          do i=1,7
            gridf(i,itype) = h/60480 * r(i) * sum(lagrange(i,1:n0-1))
          enddo
          r(1)  = r(n0)
        endif

        ! Calculate Simpson-integration coefficients from r(n0) to r(n)
        do i = 1,nstep
          do j = 2,7
            r(j) = r(j-1)*dr
          enddo
          do j = n0,n0+6
            gridf(j,itype) = gridf(j,itype) + h/140 * r(j-n0+1) * simpson(j-n0+1)
          enddo
          n0   = n0 + 6
          r(1) = r(7)
        enddo

      enddo

      end

c     -----------------------------

c     Integrates function f numerically (Lagrange and Simpson integration) on grid(itype)
c     and is much faster than intgr. (Only normal outward integration.)
c
c     It assumes a power-law behavior for small radii: f(r) = c*r**x with x > 0.
c     If x < 0, a warning is issued and a linear function f(r) = c*r is used instead.
c
c     Before first use of this function it has to be initialized with intgrf_init.
c
c# define dumpradial
      function intgrf(f,itype)

      use util,   only: chr
      use file
# ifdef dumpradial
#   warning dumpradial defined
      use global, only: grid,gridf,rgrid
# else
      use global, only: grid,gridf
# endif

      use, intrinsic :: iso_fortran_env
      implicit none
      real_dp             :: intgrf
      integer, intent(in) :: itype
      real_dp, intent(in) :: f(*)
      integer             :: n
      real_dp             :: r1,h,x
      integer, save       :: iwarn = 0
# ifdef dumpradial
      integer             :: iunit,i
# endif

      n     = grid(itype)%number
      r1    = grid(itype)%first
      h     = grid(itype)%increment

      ! integral from 0 to r1 approximated by leading term in power series expansion
c      ! f(r) = a+c*r**x
c      if (f(1)*f(2)>1d-12) then
c        if(f(2)==f(1)) then
c          intgrf = r1*f(1)
c        else
c          x      = (f(3)-f(2)) / (f(2)-f(1))
c          a      = (f(2)-x*f(1)) / (1-x)
c          x      = log(x)/h
c          if(x<0) then
c            if(x>-1) write(6,'(A,ES9.1)') 'intgrf: Warning! Negative exponent x in extrapolation a+c*r**x:',x
c            if(x<=-1) write(6,'(A,ES9.1)') 'intgrf: Negative exponent x in extrapolation a+c*r**x:',x
c            if(x<=-1) Error('Negative exponent x in extrapolation a+c*r**x')
c          endif
c          intgrf = r1*(f(1)+x*a) / (x+1)
c        endif
c      else
c        intgrf = 0
c      endif
      ! With f(r) = r**x
      if( f(1)>1d-10 .and. f(2)>1d-10 ) then
        x = log ( f(2)/f(1) ) / h
c        if(x<0) then
c          Warn('Negative exponent x in extrapolation c*r**x:')
c          intgrf = r1*f(1) / 2
c        else
c          intgrf = r1*f(1) / (x+1)
c        endif
        if(x<=-1) then
          iwarn = iwarn + 1
          if(iwarn<=10) then
            Warn('Exponent x in extrapolation c*r**x below -1: '//Chf(x,'ES9.1')//'  Atom type: '//Chr(itype)) !Integrand dumped to "spex.radial"
            if(iwarn==10) Warn('Reached 10 warnings in this routine. Further warnings not reported anymore...')
          endif
# ifdef dumpradial
          iunit = fopen('spex.radial',status='unknown')
          write(iunit,'(F15.10,F25.10)') (rgrid(i,itype),f(i),i=1,n)
          call fclose(iunit)
          Error('Exponent below -1. Integrand written to spex.radial.')
# endif
          intgrf = r1*f(1) / 2
        else
          if(x<-0.5) then
            iwarn = iwarn + 1
            if(iwarn<=10) then
              Warn('Large negative exponent x in extrapolation c*r**x: '//Chf(x,'ES9.1')//'  Atom type: '//Chr(itype))
              if(iwarn==10) Warn('Reached 10 warnings in this routine. Further warnings not reported anymore...')
            endif
          endif
          intgrf = r1*f(1) / (x+1)
        endif
      else
        intgrf = 0
      endif

      ! integrate from r(1) to r(n) by multiplying with gridf
      intgrf = intgrf + dot_product(gridf(:n,itype),f(:n))

      end

c     -----------------------------

c     Same for complex integrands

      function cintgrf(f,itype)

      use global, only: grid,gridf
      use util,   only: chr

      use, intrinsic :: iso_fortran_env
      implicit none
      complex_dp             :: cintgrf
      integer,    intent(in) :: itype
      complex_dp, intent(in) :: f(*)
      integer                :: n
      complex_dp             :: x
      real_dp                :: r1,h

      n  = grid(itype)%number
      r1 = grid(itype)%first
      h  = grid(itype)%increment

      ! With f(r) = r**x
      if (real(f(1)*f(2))>1d-20) then
        x = log ( f(2)/f(1) ) / h
c        if(real(x)<0) then
c          Warn('Negative exponent x in extrapolation c*r**x: '//Chf(x,'ES9.1'))
c          cintgrf = r1*f(1) / 2
c        else
c          cintgrf = r1*f(1) / (x+1)
c        endif
        if(real(x)<-1) then
          Warn('Exponent x in extrapolation c*r**x below -1: '//Chf(real(x),'ES9.Ã1'))
          cintgrf = r1*f(1) / 2
        else
          if(real(x)<-0.5) Warn('Large negative exponent x in extrapolation c*r**x: '//Chf(real(x),'ES9.1'))
          cintgrf = r1*f(1) / (x+1)
        endif
      else
        cintgrf = 0
      endif

      ! integrate from r(1) to r(n) by multiplying with gridf
      cintgrf = cintgrf + dot_product(gridf(:n,itype),f(:n))

      end

c     -----------------------------

c     Calculates the primitive of f, on grid(itypein):
c
c                   r
c     primf(r) = integral f(r') dr'   ( r = grid point )
c                   0
c
c     If itypein is negative, the primitive
c
c                   R
c     primf(r) = integral f(r') dr'   ( R = MT sphere radius )
c                   r
c
c     is calculated instead.
c
c     The routine primitivef is much faster.

      subroutine primitive(primf,f,itypein)

      use global, only: grid,gridtype

      use, intrinsic :: iso_fortran_env
      implicit none
      integer, intent(in)  :: itypein
      real_dp, intent(out) :: primf(grid(abs(itypein))%number)
      real_dp, intent(in)  ::     f(grid(abs(itypein))%number)
      real_dp              :: h
      integer              :: itype,n,i
      type(gridtype)       :: grid1
      real_dp              :: intgr

      itype = abs(itypein)

      n = grid(itype)%number
      h = grid(itype)%increment

      if(itypein>0) then
        grid1%first     = grid(itype)%first  ! perform outward integration
        grid1%increment = h                  ! (from 0 to r)
      else
        grid1%first     = grid(itype)%radius ! perform inward integration
        grid1%increment = -h                 ! (from MT sphere radius to r)
      endif

      do i=1,n
        grid1%number = i
        if(itypein>0) then
          primf(i)     =  intgr(f,grid1)
        else
          Error('check whether inward integration works here. Better use primitivef.')
          primf(n-i+1) = -intgr(f(n:1:-1),grid1)
        endif
      enddo

      end

c     -----------------------------

c     Fast calculation of primitive.
c     Integration is done differently from routine primitive. Only Lagrange integration is used
c     Usage is identical to routine primitive.
c     Uses power-law behavior for small radii: c*r**x with x > -1.
      subroutine primitivef(primf,fin,itypein)

      use global, only: grid
      use util,   only: chr
      use, intrinsic :: iso_fortran_env
      implicit none
      integer, intent(in)  :: itypein
      real_dp, intent(out) :: primf(grid(abs(itypein))%number)
      real_dp, intent(in)  ::   fin(grid(abs(itypein))%number)
      real_dp              ::     f(grid(abs(itypein))%number)
      real_dp              ::     r(grid(abs(itypein))%number)
      real_dp, parameter   :: lagrange(7,6)= reshape(
     &                        [ 19087d0,65112d0,-46461d0, 37504d0,-20211d0, 6312d0, -863d0,
     &                           -863d0,25128d0, 46989d0,-16256d0,  7299d0,-2088d0,  271d0,
     &                            271d0,-2760d0, 30819d0, 37504d0, -6771d0, 1608d0, -191d0,
     &                           -191d0, 1608d0, -6771d0, 37504d0, 30819d0,-2760d0,  271d0,
     &                            271d0,-2088d0,  7299d0,-16256d0, 46989d0,25128d0, -863d0,
     &                           -863d0, 6312d0,-20211d0, 37504d0,-46461d0,65112d0,19087d0 ], ! The last row is actually never used.
     &                        [7,6])
      real_dp              :: h,x,h1
      integer              :: itype,n,i,n0
      real_dp              :: intgr,r1,dr,fr(7)

      itype = abs(itypein)

      primf = 0

      n = grid(itype)%number
      h = grid(itype)%increment

      if(itypein>0) then
        r1 = grid(itype)%first  ! perform outward integration
        f  = fin                ! (from 0 to r)
      else
        r1 = grid(itype)%radius ! perform inward integration
        h  = -h                 ! (from MT sphere radius to r)
        f  = fin(n:1:-1)        !
      endif

      ! integral from 0 to r1 approximated by leading term in power series expansion (only if h>0)
c      if(h>0.and.f(1)*f(2)>1d-20) then
c        if(f(2)==f(1)) then
c          intgr = r1*f(1)
c        else
c          x     = (f(3)-f(2))/(f(2)-f(1))
c          a     = (f(2)-x*f(1)) / (1-x)
c          x     = log(x)/h
c          if(x<0) then
c            if(x>-1) write(6,'(A,ES9.1)') 'intgr: Warning! Negative exponent x in extrapolation a+c*r**x:',x
c            if(x<=-1) write(6,'(A,ES9.1)') 'intgr: Negative exponent x in extrapolation a+c*r**x:',x
c            if(x<=-1) Error('Negative exponent x in extrapolation a+c*r**x')
c          endif
c          intgr = r1*(f(1)+x*a) / (x+1)
c        endif
c      else
c        intgr = 0
c      endif
      if (h>0.and.f(1)*f(2)>1d-20) then
        x = log ( f(2)/f(1) ) / h
        if(x<-1)   write(6,'(A,ES9.1)') 'primitivef: Exponent x in extrapolation c*r**x below -1: ',x
        if(x<-1)   Error('Exponent x in extrapolation c*r**x below -1.')
        if(x<-0.5) Warn('Large negative exponent x in extrapolation c*r**x: '//Chf(x,'ES9.1'))
        intgr = r1*f(1) / (x+1)
      else
        intgr = 0
      endif

      primf(1) = intgr
      dr       = exp(h)
      r(1)     = r1
      n0       = 1
      h1       = h/60480

      ! Lagrange integration from r(n0) to r(n0+5)
 1    do i=2,7
        r(i) = r(i-1)*dr
      enddo
      fr = f(n0:n0+6) * r(:7)
      do i=1,6
        intgr = intgr + h1 * dot_product(lagrange(:,i),fr)
        if(primf(n0+i)==0) primf(n0+i) = intgr ! avoid double-definition
      enddo
      if(n0+12<=n) then
        r(1) = r(7)
        n0   = n0 + 6
        goto 1
      else if(n0+6<n) then
        r(1)  = r(n-5-n0)
        n0    = n-6
        intgr = primf(n-6)
        goto 1
      endif

      if(itypein<0) then       !
        primf = -primf(n:1:-1) ! Inward integration
      endif                    !

      end

c     -----------------------------

c     Returns derivative of f in df.
      subroutine derivative(df,f,itype)
      use global, only: grid
      use, intrinsic :: iso_fortran_env
      implicit none
      integer, intent(in)  :: itype
      real_dp, intent(in)  :: f(grid(itype)%number)
      real_dp, intent(out) :: df(grid(itype)%number)
      real_dp              :: h,r,d21,d32,d43,d31,d42,d41,df1,df2,s
      real_dp              :: y0,y1,y2
      integer              :: i,n
      n = grid(itype)%number
      h = grid(itype)%increment
      r = grid(itype)%first
      ! use Lagrange interpolation of 3rd order (and averaging) for points 3 to n
      d21 = r * (exp(h)-1) ; d32 = d21 * exp(h) ; d43 = d32 * exp(h)
      d31 = d21 + d32      ; d42 = d32 + d43
      d41 = d31 + d43
      df(1) =   d31*d41 / (d21*d32*d42) * f(2) + ( -1d0/d21 - 1d0/d31 - 1d0/d41) * f(1)
     &        - d21*d41 / (d31*d32*d43) * f(3) + d21*d31 / (d41*d42*d43) * f(4)
      df(2) = - d32*d42 / (d21*d31*d41) * f(1) + (  1d0/d21 - 1d0/d32 - 1d0/d42) * f(2)
     &        + d21*d42 / (d31*d32*d43) * f(3) - d21*d32 / (d41*d42*d43) * f(4)
      df1   =   d32*d43 / (d21*d31*d41) * f(1) - d31*d43 / (d21*d32*d42) * f(2) +
     &  ( 1d0/d31 + 1d0/d32 - 1d0/d43 ) * f(3) + d31*d32 / (d41*d42*d43) * f(4)
      do i = 3,n-2
        d21 = d32 ; d32 = d43 ; d43 = d43 * exp(h)
        d31 = d42 ; d42 = d42 * exp(h)
        d41 = d41 * exp(h)
        df2   = - d32*d42 / (d21*d31*d41) * f(i-1) + ( 1d0/d21 - 1d0/d32 - 1d0/d42) * f(i) +
     &            d21*d42 / (d31*d32*d43) * f(i+1) - d21*d32 / (d41*d42*d43) * f(i+2)
        df(i) = ( df1 + df2 ) / 2
        df1   =   d32*d43 / (d21*d31*d41) * f(i-1) - d31*d43 / (d21*d32*d42) * f(i) +
     &    ( 1d0/d31 + 1d0/d32 - 1d0/d43 ) * f(i+1) + d31*d32 / (d41*d42*d43) * f(i+2)
      enddo
      df(n-1) = df1
      df(n)   = - d42*d43 / (d21*d31*d41) * f(n-3) + d41*d43 / (d21*d32*d42) * f(n-2) -
     &            d41*d42 / (d31*d32*d43) * f(n-1) + ( 1d0/d41 + 1d0/d42 + 1d0/d43 ) * f(n)
      ! for first two points use Lagrange interpolation of second order for log(f(i))
      ! or, as a fall-back, Lagrange interpolation with the conditions f(1), f(2), f(3), f'(3).
      s = sign(1d0,f(1))
      if(sign(1d0,f(2))/=s.or.sign(1d0,f(3))/=s.or.any(f(:3)==0)) then
        d21   = r * (exp(h)-1)
        d32   = d21 * exp(h)
        d31   = d21 + d32
        s     = df(3) / (d31*d32) - f(1) / (d21*d31**2) + f(2) / (d21*d32**2) - f(3) / (d31**2*d32) - f(3) / (d31*d32**2)
        df(1) = - (d21+d31) / (d21*d31) * f(1) + d31 / (d21*d32) * f(2) - d21 / (d31*d32) * f(3) + d21*d31 * s

        df(2) = - (d21-d32) / (d21*d32) * f(2) - d32 / (d21*d31) * f(1) + d21 / (d31*d32) * f(3) - d21*d32 * s
      else
        y0    = log(abs(f(1)))
        y1    = log(abs(f(2)))
        y2    = log(abs(f(3)))
        df(1) = ( - 3*y0/2 + 2*y1 - y2/2 ) * f(1) / (h*r)
        df(2) = (y2-y0)/2                  * f(2) / (h*r*exp(h))
      endif
      end

c     -----------------------------

c     Returns the Legendre polynomials P_l(x) for l=0,...,ll.
      subroutine legendre(leg,x,ll)
      
      use, intrinsic :: iso_fortran_env

      implicit none
      integer, intent(in)  :: ll
      real_dp, intent(in)  :: x
      real_dp, intent(out) :: leg(0:ll)
      integer              :: l      
      if(ll<0) Bug('Input l quantum number negative.')
      leg(0) = 1 ; if(ll==0) return
      leg(1) = x
      do l = 1,ll-1
        leg(l+1) = ( (2*l+1)*x*leg(l) - l*leg(l-1) ) / (l+1)
      enddo
      end

c     -----------------------------

c     Returns the spherical harmonics Y_lm(theta,phi)
c     for l = 0,...,ll in Y(1,...,(ll+1)**2).
      subroutine harmonics(Y,theta,phi,ll)

      use global, only: img
      use, intrinsic :: iso_fortran_env
      implicit none
      integer,    intent(in)  :: ll
      real_dp,    intent(in)  :: theta,phi
      complex_dp, intent(out) :: Y((ll+1)**2)
      real_dp                 :: stheta,ctheta

      if(ll<0) return

      stheta = sin(theta)
      ctheta = cos(theta)

      call harmonicsr(Y,[stheta*cos(phi),stheta*sin(phi),ctheta],ll)

      end

c     -----------

c     Returns the spherical harmonics Y_lm(^rvec)
c     for l = 0,...,ll in Y(1,...,(ll+1)**2).
      subroutine harmonicsr(Y,rvec,ll)

      use global, only: img

      use, intrinsic :: iso_fortran_env
      implicit none
      integer,    intent(in)  :: ll
      real_dp,    intent(in)  :: rvec(3)
      complex_dp, intent(out) :: Y((ll+1)**2)
      complex_dp              :: c,cxy
      integer                 :: l,m,lm,lm1,lm2
      integer,    parameter   :: msqrtn0 = 31
      real_dp,    parameter   :: sqrtn0(msqrtn0) = [ (sqrt(dble(l)),l=1,msqrtn0) ]      
      real_dp,    parameter   :: thrsh = 1d-16
      real_dp                 :: sqrtn(2*ll+1)
      real_dp                 :: ctheta,r

      if(ll<0) return

      Y(1) = 0.282094791773878d0
      if(ll==0) return

      ctheta = 0
      r      = sqrt( rvec(1)**2 + rvec(2)**2 + rvec(3)**2 )
      if(r>thrsh) then
        cxy    = ( rvec(1) + img*rvec(2) ) / r
        ctheta = rvec(3) / r
      else
        Y(2:) = 0d0
        return
      endif

c     precalculate sqrt(n), n = 1,...,2*ll+1
      l         = min(2*ll+1,msqrtn0)
      sqrtn(:l) = sqrtn0(:l)
      do l = msqrtn0+1,2*ll+1
        sqrtn(l) = sqrt(dble(l))
      enddo

c     define Y,l,-l and Y,l,l
      r   = Y(1)
      c   = 1
      m   = -1
      do l=1,ll
        r           = r * sqrtn(2*l+1) / sqrtn(2*l)
        c           = c * cxy
        Y(l**2+1)   = r * conjg(c)  ! l,-l
        Y((l+1)**2) = m * r * c     ! l,l
        m           = -m
      enddo

c     define Y,l,-l+1 and Y,l,l-1
      Y(3) = 0.48860251190292d0 * ctheta
      do l = 2,ll
        r          = sqrtn(2*l+1) * ctheta
        Y(l**2+2)  = r * Y((l-1)**2+1) ! l,-l+1
        Y(l*(l+2)) = r * Y(l**2)       ! l,l-1
      enddo

c     define Y,l,m, |m|<l-1
      do l = 2,ll
        lm  = l**2 + 2
        lm1 = lm - 2*l
        lm2 = lm - 4*l + 2
        do m = -l+2,l-2
          lm    = lm + 1
          lm1   = lm1 + 1
          lm2   = lm2 + 1
          Y(lm) = sqrtn(2*l+1) / sqrtn(l+m) / sqrtn(l-m) * (
     &            sqrtn(2*l-1) * ctheta                      * Y(lm1)-
     &            sqrtn(l+m-1) * sqrtn(l-m-1) / sqrtn(2*l-3) * Y(lm2) )
        enddo
      enddo

      end

c     -----------

c     Calculates the Wigner 3j symbols using Racah's formula
      function wigner3j(l1,l2,l3,m1,m2,m3)
      use global, only: fac, sfac
      use, intrinsic :: iso_fortran_env
      implicit none
      real_dp             :: wigner3j
      integer, intent(in) :: l1,l2,l3,m1,m2,m3
      integer             :: tmin,tmax,t,f1,f2,f3,f4,f5
      wigner3j = 0
c     The following IF clauses should be in the calling routine and commented here.
c      if(-m3/=m1+m2)  return
c      if(abs(m1)>l1) return
c      if(abs(m2)>l2) return
c      if(abs(m3)>l3) return
c      if(l3<abs(l1-l2).or.l3>l1+l2) return
      f1   = l3-l2+m1
      f2   = l3-l1-m2
      f3   = l1+l2-l3
      f4   = l1-m1
      f5   = l2+m2
      tmin = max(0,-f1,-f2) ! The arguments to fac (see below)
      tmax = min(f3,f4,f5)  ! must not be negative.
      ! The following line is only for testing and should be removed at a later time.
      if(tmax-tmin/=min(l1+m1,l1-m1,l2+m2,l2-m2,l3+m3,l3-m3,l1+l2-l3,l1-l2+l3,-l1+l2+l3))
     &  Error('Number of terms incorrect.')
      if(tmin<=tmax) then
        do t = tmin,tmax
          wigner3j = wigner3j + (-1)**t /
     &                          ( fac(t) * fac(f1+t) * fac(f2+t) * fac(f3-t) * fac(f4-t) * fac(f5-t) )
        enddo
        wigner3j = wigner3j * (-1)**(l1-l2-m3) *
     &                        sfac(l1+l2-l3) * sfac(l1-l2+l3) * sfac(-l1+l2+l3) / sfac(l1+l2+l3+1) *
     &                        sfac(l1+m1)   * sfac(l1-m1) *
     &                        sfac(l2+m2)   * sfac(l2-m2) *
     &                        sfac(l3+m3)   * sfac(l3-m3)
      endif
      end

c     ----------------

c     Calculates Gaunt coefficients, i.e. the integrals of three spherical harmonics
c     integral ( conjg(Y(l1,m1)) * Y(l2,m2) * conjg(Y(l3,m3)) )
c     They are also the coefficients C(l1,l2,l3,m1,m2,m3) in
c     conjg(Y(l1,m1)) * Y(l2,m2) = sum(l3,m3) C(l1,l2,l3,m1,m2,m3) Y(l3,m3)
      function gaunt(l1,l2,l3,m1,m2,m3)
      use global, only: pi
      use, intrinsic :: iso_fortran_env
      implicit none
      real_dp             :: gaunt
      integer, intent(in) :: l1,l2,l3,m1,m2,m3
      real_dp             :: wigner3j
      gaunt = 0
      if(m3/=m2-m1)  return
      if(abs(m1)>l1) return
      if(abs(m2)>l2) return
      if(abs(m3)>l3) return
      if(mod(l1+l2+l3,2)/=0) return
      if(l3<abs(l1-l2).or.l3>l1+l2) return
      gaunt = (-1)**(m1+m3) *
     &        sqrt((2*l1+1)*(2*l2+1)*(2*l3+1)/pi/4)*
     &        wigner3j(l1,l2,l3,-m1,m2,-m3)*
     &        wigner3j(l1,l2,l3, 0,  0,  0)
      end

c     --------------

c     see Varshalovich, Moskalev, Khersonskii, "Quantum Theory of Angular Momentum"
      subroutine wignerrot(mat,alpha,beta,gamma,l) ! not tested yet
      use global, only: fac,sfac,img
      use, intrinsic :: iso_fortran_env
      implicit none
      integer,    intent(in)  :: l
      complex_dp, intent(out) :: mat(-l:l,-l:l)
      real_dp,    intent(in)  :: alpha,beta,gamma
      real_dp                 :: cosb2,sinb2,rdum
      integer                 :: m1,m2,f1,f2,f3,tmin,tmax,t
      cosb2 = cos(beta/2)
      sinb2 = sin(beta/2)
      do m1=-l,l
        do m2=-l,l
          rdum = 0
          f1   = l-m1
          f2   = l-m2
          f3   = m1+m2
          tmin = max(0,-f3) ! The arguments to fac (see below)
          tmax = min(f1,f2) ! must be non-negative.
          do t = tmin,tmax
            rdum = rdum + (-1)**t * cosb2**(m1+m2+2*t) * sinb2**(2*l-m1-m2-2*t) /
     &                    ( fac(t) * fac(f1-t) * fac(f2-t) * fac(f3+t) )
          enddo
          if(rdum/=0) then
            rdum = rdum * (-1)**(l-m2) * sfac(l+m1) * sfac(l-m1) * sfac(l+m2) * sfac(l-m2)
            mat(m1,m2) = exp(-img*m1*alpha) * rdum * exp(-img*m2*gamma)
          else
            mat(m1,m2) = 0
          endif
        enddo
      enddo
      end

c     -------------

c     Calculates the Bessel functions (not used)
      subroutine bessel(bes,x,l)
      use global, only: pi
      use, intrinsic :: iso_fortran_env
      implicit none
      integer, intent(in)  :: l
      real_dp, intent(out) :: bes(0:l)
      real_dp, intent(in)  :: x
      real_dp              :: fac1,fac2,y,y2
      integer              :: i
      if(l<=0) Error('l<1!')
      bes = 0
      ! Evaluate j_0 and j_1 (this is somewhat unstable)
      if(x>30d0) then
        y    = x-pi/4
        fac1 = sqrt(2d0/pi/x)
        bes(0) = fac1 * cos(y)
        bes(1) = fac1 * sin(y)
      else
        fac1      = 1
        fac2      = 1
        i         = 0
        y         = 0.25d0*x**2
        y2        = y**2
        bes(0)    = 1-y
        bes(1)    = 1-y/2
        do while (max(abs(fac1),abs(fac2))>1d-10)
          i    = i + 2
          fac1 = fac1/(i**2-i)**2        * y2
          fac2 = fac2/((i-1)*i**2*(i+1)) * y2
          bes(0) = bes(0) + fac1 * (1-y/(i+1)**2)
          bes(1) = bes(1) + fac2 * (1-y/((i+1)*(i+2)))
        enddo
        bes(1) = x/2 * bes(1)
      endif
      ! Evaluate j_l (l>1) by recursion
      do i=1,l-1
        bes(i+1) = 2*i*bes(i)/x - bes(i-1)
      enddo
      end

c     ------------

c     Calculates the spherical Bessel functions of orders 0 to l at x
c     by backward recurrence using j_l(x) = (2l+3)/x j_l+1(x) - j_l+2(x) .
c     (Starting points are calculated according to Zhang, Min, "Computation of Special Functions".)
      subroutine sphbessel(sphbes,x,l)
      use, intrinsic :: iso_fortran_env
      implicit none
      integer, intent(in)  :: l
      real_dp, intent(in)  :: x
      real_dp, intent(out) :: sphbes(0:l)
      real_dp              :: s0,s1,f,f0,f1,cs
      integer              :: ll,lsta,lmax
      if(x>100) call sphbessel_simple(sphbes,x,l)
      if(l<0) return
      if(x<0) then
        Bug('negative argument.')
      else if(x==0) then
        sphbes(0)  = 1d0
        do ll = 1,l
          sphbes(ll) = 0d0
        enddo
        return
      endif
      sphbes(0) = sin(x) / x
      sphbes(1) = ( sphbes(0) - cos(x) ) / x
      if(l<=1) return
      s0   = sphbes(0)
      s1   = sphbes(1)
      lsta = lsta1(x,200)      !
      lmax = l                 !
      if(lsta<l) then          !
        lmax            = lsta ! determine starting point lsta
        sphbes(lmax+1:) = 0d0  ! for backward recurrence
      else                     !
        lsta = lsta2(x,l,15)   !
      endif                    !
      f0 = 0d0                                                    !
      f1 = 1d-200                                                 !
      do ll = lsta,0,-1                                           ! backward recurrence
        f  = f1 / x * (2*ll+3) - f0 ; if(ll<=lmax) sphbes(ll) = f ! with arbitrary start values
        f0 = f1                                                   !
        f1 = f                                                    !
        if(ll>lmax) then ! avoids NAN, should be improved ...
          f0 = f0 / f
          f1 = f1 / f
          f  = 1d0
        endif
      enddo                                   !
      if(abs(s0)>abs(s1)) then ; cs = s0 / f  !
      else                     ; cs = s1 / f0 ! scale to correct values
      endif                                   !
      sphbes = cs * sphbes                    !
      contains
c     Test starting point
      function lsta0(x,mp)
      implicit none
      integer             :: lsta0
      integer, intent(in) :: mp
      real_dp, intent(in) :: x
      real_dp             :: f,lgx
      lgx   = log10(x)
      lsta0 = 0
      f     = lgx
      do while(f>-mp)
        lsta0 = lsta0 + 1
        f     = f + lgx - log10(2d0*lsta0+1)
      enddo
      end function lsta0
c     Returns starting point lsta1 for backward recurrence such that sphbes(lsta1) approx. 10^(-mp).
      function lsta1(x,mp)
      implicit none
      integer             :: lsta1
      integer, intent(in) :: mp
      real_dp, intent(in) :: x
      real_dp             :: f0,f1,f
      integer             :: n0,n1,nn,it
      n0 = int(1.1*x) + 1
      f0 = envj(n0,x) - mp
      n1 = n0 + 5
      f1 = envj(n1,x) - mp
      do it = 1,20
        nn = n1 - (n1-n0) / (1d0-f0/f1)
        f  = envj(nn,x) - mp
        if(abs(nn-n1)<1) exit
        n0 = n1
        f0 = f1
        n1 = nn
        f1 = f
      enddo
      lsta1 = nn
      end function lsta1
c     Returns the starting point lsta2 for backward recurrence such that all sphbes(l) have mp significant digits.
      function lsta2(x,l,mp)
      implicit none
      integer             :: lsta2
      integer, intent(in) :: l,mp
      real_dp, intent(in) :: x
      real_dp             :: f0,f1,f,hmp,ejn,obj
      integer             :: n0,n1,nn,it
      hmp = 0.5d0 * mp
      ejn = envj(l,x)
      if(ejn<=hmp) then
        obj = mp
        n0  = int(1.1*x) + 1
      else
        obj = hmp + ejn
        n0  = l
      endif
      f0 = envj(n0,x) - obj
      n1 = n0 + 5
      f1 = envj(n1,x) - obj
      do it = 1,20
        nn = n1 - (n1-n0) / (1d0-f0/f1)
        f  = envj(nn,x) - obj
        if(abs(nn-n1)<1) exit
        n0 = n1
        f0 = f1
        n1 = nn
        f1 = f
      enddo
      lsta2 = nn + 10
      end function lsta2
      function envj(l,x)
      implicit none
      real_dp             :: envj
      real_dp, intent(in) :: x
      integer, intent(in) :: l
      envj = 0.5d0 * log10(6.28d0*l) - l*log10(1.36d0*x/l)
      end function envj
      end

c     -----------

c     Calculates the spherical Bessel functions of orders 0 to l at x
c     (less robust than sphbessel but faster)
      subroutine sphbessel_simple(sphbes,x,l)
      use, intrinsic :: iso_fortran_env
      implicit none
      integer, intent(in)  :: l
      real_dp, intent(in)  :: x
      real_dp, intent(out) :: sphbes(0:l)
      integer              :: ll
      if(l<0) return
      if(x<0) then
        Bug('negative argument.')
      else if(x==0) then
        sphbes(0)  = 1d0
        do ll = 1,l
          sphbes(ll) = 0d0
        enddo
        return
      endif
      sphbes(0) = sin(x) / x
      sphbes(1) = ( sphbes(0) - cos(x) ) / x
      if(l<=1) return
      do ll = 2,l
        sphbes(ll) = (2*ll-1)*sphbes(ll-1)/x - sphbes(ll-2)
      enddo
      end

c     -----------

c     Calculates the integral of sphbes(y,l)*y**l from zero to x  ! not used
      subroutine sphbespow_integral(sphbesp,sphbes,x,l)
      use, intrinsic :: iso_fortran_env
      implicit none
      integer, intent(in)  :: l
      real_dp, intent(out) :: sphbesp(0:l)
      real_dp, intent(in)  :: sphbes(0:l),x
      real_dp              :: Si
      integer              :: ll
      sphbesp(0) = Si(x)
      do ll=1,l
        sphbesp(ll) = (2*ll-1)*sphbesp(ll-1) - x**ll * sphbes(ll-1)
      enddo
      end
c     Calculates the integral of sphbes(y,l)/y**(l+1) from infinity to x  ! not used
      subroutine sphbespowi_integral(sphbesp,sphbes,x,l)
      use, intrinsic :: iso_fortran_env
      implicit none
      integer, intent(in)  :: l
      real_dp, intent(out) :: sphbesp(0:l)
      real_dp, intent(in)  :: sphbes(0:l),x
      real_dp              :: Ci
      integer              :: ll
      sphbesp(0) = Ci(x) - sphbes(0)
      do ll=1,l
        sphbesp(ll) = ( sphbesp(ll-1) - sphbes(ll)/x**ll ) / (2*ll+1)
      enddo
      end

c     -----------

c     Calculates the spherical Neumann functions
      subroutine sphneumann(sphneu,x,l)
      use, intrinsic :: iso_fortran_env
      implicit none
      integer, intent(in)  :: l
      real_dp, intent(out) :: sphneu(0:l)
      real_dp, intent(in)  :: x
      integer              :: ll
      sphneu(0) = -cos(x)/x
      if(l==0) return
      sphneu(1) = sphneu(0)/x - sin(x)/x
      do ll=1,l-1
        sphneu(ll+1) = (2*ll+1) * sphneu(ll) / x - sphneu(ll-1)
      enddo
      end

c     -----------

c     Calculates the sine integral ( i.e. integral of sin(y)/y from zero to x)
      function Si(x)
      use global, only: pi,img
      use, intrinsic :: iso_fortran_env
      implicit none
      real_dp             :: Si
      real_dp, intent(in) :: x
      real_dp             :: fac,x2,x1
      complex_dp          :: cdum,c,d,delta
      integer             :: i
      x2 = x**2
      if(abs(x)>2) then ! continued fraction using Lentz's method
        x1   = abs(x)
        cdum = 1 + img*x1
        c    = cdum
        d    = 0
        i    = 0
        do
          i     = i + 1
          d     = ( 2*i+1 + img*x1 - i**2 * d )**(-1)
          c     =   2*i+1 + img*x1 - i**2 / c
          delta = c*d
          cdum  = cdum * delta
          if(abs(1-delta)<1d-15) exit
        enddo
        cdum = exp(-img*x1) / cdum
        Si   = imag(cdum) + pi/2
        if(x<0) Si = -Si
      else ! McLaurin series
        fac = x
        Si  = x
        i   = 1
        do while(abs(fac)>1d-16)
          i   = i + 2
          fac = -fac * x2 / ((i-1)*i)
          Si  = Si + fac/i
        enddo
      endif
      end

c     -------------

c     Calculates the cosine integral ( i.e. integral of cos(y)/y from +infinity to x)
      function Ci(x)
      use global, only: img
      use, intrinsic :: iso_fortran_env
      implicit none
      real_dp             :: Ci
      real_dp, intent(in) :: x
      real_dp             :: x2,fac
      real_dp, parameter  :: euler = 0.57721566490153286060651209008240243d0
      complex_dp          :: cdum,c,d,delta
      integer             :: i
      if(x<=0) Error('Argument is non-positive')
      x2 = x**2
      if(abs(x)>2) then ! continued fraction using Lentz's method
        cdum = 1 + img*x
        c    = cdum
        d    = 0
        i    = 0
        do
          i     = i + 1
          d     = ( 2*i+1 + img*x - i**2 * d )**(-1)
          c     =   2*i+1 + img*x - i**2 / c
          delta = c*d
          cdum  = cdum * delta
          if(abs(1-delta)<1d-15) exit
        enddo
        cdum = exp(-img*x) / cdum
        Ci   = - real(cdum)
      else ! McLaurin series
        fac = 1
        Ci  = 0
        i   = 0
        do while(abs(fac)>1d-16)
          i   = i + 2
          fac = -fac * x2 / ((i-1)*i)
          Ci  = Ci + fac/i
        enddo
        Ci = euler + log(x) + Ci
      endif
      end

c     --------------

c     Calculates the complex E1 function ( i.e. integral of exp(-u)/u from z to infinity which is related to the
c     complex exponential function Ei by E1(z) = -Ei(-z) )
c     Either an expansion or a continued fraction is evaluated.
      function cE1(z)
      use global, only: pi,img
      use, intrinsic :: iso_fortran_env
      implicit none
      complex_dp             :: cE1
      complex_dp, intent(in) :: z
      complex_dp             :: cdum
      real_dp                :: x,y,znorm
      real_dp, parameter     :: euler = 0.57721566490153286060651209008240243d0
      integer                :: i,nfrac
      x     = real(z)
      y     = imag(z)
      znorm = abs(z)
      if(znorm==0) then
        cE1 = 0
        return
        Error('lower bound of the E1-function integral must not be zero.')
      else if(znorm<=10.or.x<0.and.abs(y)<=2) then ! This is somewhat experimental.
        cE1  = 1d0
        cdum = 1d0
        do i = 1,10000
          cdum = - cdum * i * z / (i+1d0)**2
          cE1  = cE1 + cdum
          if(abs(cdum)<=abs(cE1)*1d-15) goto 1
        enddo
        write(6,'(A'NoA) 'cEi: expansion for complex E1 function with argument' ; write(6,*) z
        write(6,'(A)')   '     not converged.'
        Error('expansion for complex E1 function not converged.')
 1      cE1 = - euler - log(z) + z * cE1
        write(*,*) 'cE1-iter',i
      else
        nfrac = 120
        cdum  = 0d0
        do i = nfrac,1,-1
          cdum = i / ( 1 + i / (z+cdum) )
        enddo
        cdum = 1d0 / (z+cdum)
        cE1  = exp(-z) * cdum
      endif
      if(x<=0.and.y==0) cE1 = real(cE1) - pi*img
      end

c     --------------

c     Calculates
c               z2
c     exp(z) * INT exp(-u)/u du  with  z = (z1+z2)/2 .
c               z1
c     ( M.S. Kluskens, IEEE Trans. Antennas Propagat. 47, 803 (1999) )
      recursive function cE1_diff(z1,z2) result(output)
      use, intrinsic :: iso_fortran_env
      implicit none
      complex_dp             :: output
      complex_dp, intent(in) :: z1,z2
      complex_dp             :: dz,z,a,b,dz2,dz2z2
      real_dp                :: h,zn,dzn
      integer                :: i,nstep
      if(z1==z2) then ; output = 0 ; return ; endif
      z   = (z1+z2)/2 ; zn  = abs(z)
      dz  = z2-z      ; dzn = abs(dz)
      dz2 = dz**2
      h   = zn / dzn
      if(dzn>20) then
c        write(*,'(A'NoA) '#'
        output = cE1_diff(z1,z) * exp( z - (z1+dz/2) ) +
     &           cE1_diff(z,z2) * exp( z - (z2-dz/2) )
      else if(h>=1) then
        dz2z2  = dz2/z**2
        a      = dz/z
        b      = dz/2
        output = a
        do i = 3,1001,2
          a      = dz2z2 / i * ( (z+i-1)*b + (i-2)*a )
          output = output + a ; if(abs(a)<=abs(output)*1d-15) goto 1
          b      = b * dz2 / (i*(i+1))
        enddo
        write(*,*) z1,z2
c        Error('difference of complex E1 function not converged')
 1      output = 2 * output
      else if(h<1) then
        nstep  = 51 ! should be enough
        a      = 0
        output = 0
        do i = nstep,1,-2
          a      = ( 1 + z*(1+z*a)/(i+1) ) / i
          output = dz2 * ( output + a ) / max(1,(i*(i-1)))
        enddo
        output = exp(z) * log(z2/z1) - 2 * output/dz
      endif

      end

c     --------------

c     Returns the complex error function.
      function cerf(z)
      use global, only: pi
      use, intrinsic :: iso_fortran_env
      implicit none
      complex_dp             :: cerf
      complex_dp, intent(in) :: z
      complex_dp             :: z1,z2,c,d,delta
      integer                :: i
      z1 = z ; if(real(z)<0) z1 = -z1
      if(real(z1)<2d0) then ! McLaurin series
        z2   = z1**2
        i    = 0
        c    = z1
        cerf = z1
        do
          i    = i + 1
          c    = -c * z2 / i
          cerf = cerf + c/(2*i+1)
          if(abs(c/(2*i+1))<1d-20) exit
        enddo
        cerf = cerf * 2/sqrt(pi)
      else ! continued fraction using Lentz's method
        d    = 0d0
        c    = z1
        cerf = z1
        i    = 0
        do
          i     = i + 1
          c     =   2*z1 + i/c
          d     = ( 2*z1 + i*d )**(-1)
          delta = c*d
          cerf  = cerf * delta ; if(abs(1-delta)<1d-15) exit
          i     = i + 1
          c     =   z1 + i/c
          d     = ( z1 + i*d )**(-1)
          delta = c*d
          cerf  = cerf * delta ; if(abs(1-delta)<1d-15) exit
          if(i==10000) Error('Lentz method not converged after 10000 steps.')
        enddo
        cerf = 1 - exp(-z1**2) / cerf / sqrt(pi)
      endif
      if(real(z)<0) cerf = -cerf
      end

c     --------------

c     Calculates the Wigner-D matrices for the rotation ROT and stores them in DWGN.
c     (More generally, the rotation can also be an improper rotation.)
c
c     DWGN is the transformation matrix for spherical harmonics wrt a rotation of the coordinate system.
c     The rotation matrix is given by ROT. Then,
c            -1                     l
c     Y  (ROT  r) = SUM Y   (r) DWGN   (ROT)
c      lm            m'  lm'        m'm
c
c     or, equivalently,
c                          l *
c     Y  (ROT*r) = SUM DWGN   (ROT) Y   (r) ,
c      lm           m'     mm'       lm'
c
c     since DWGN is unitary.
c
c     (Note that the coordinates of a vector r transform as ROT^(-1)*r, if the coordinate system rotates by ROT.)
c
c     (See D.M. Brink, G.R. Satchler, "Angular Momentum", Clarendon Press, Oxford 1968)
      subroutine dwigner(dwgn,rot,lcut)
      use global, only: img,fac,sfac
      use, intrinsic :: iso_fortran_env
      implicit none
      integer,    intent(in)  :: lcut
      complex_dp, intent(out) :: dwgn(-lcut:lcut,-lcut:lcut,0:lcut)
      real_dp,    intent(in)  :: rot(3,3)
      real_dp                 :: rot1(3,3)
      real_dp                 :: det,cosa,sina,cosb,sinb,cosc,sinc,cosb2,sinb2,nom,denom,cosb2_,sinb2_
      complex_dp              :: cdum
      integer                 :: l,m1,m2,t,tmin,tmax,nexpo
      real_dp                 :: determinant
      dwgn = 0
      det  = determinant(rot)
      if(abs(1-abs(det))>1d-10) then
        write(6,'(A)') 'dwigner: determinant neither 1 nor -1'
        Error('determinant neither 1 nor -1')
      endif
      rot1 = rot / det ! rot1 is proper rotation
      ! Determine sines and cosines of Euler angles a,b,c
      cosb  = rot1(3,3)
      if(abs(cosb)>1) then ! avoid rounding errors
        if(1-abs(cosb)>1d-8) Error('rounding error.')
        cosb = sign(1d0,cosb)
      endif
      if(1-abs(cosb)<1d-12) then ! avoid rounding errors
        sinb = 0d0
        if(cosb>0) then ; cosb2 = 1d0 ; sinb2 = 0d0
        else            ; cosb2 = 0d0 ; sinb2 = 1d0
        endif
      else
        sinb  = sqrt(1-cosb**2) ! 0 <= b <= 180 deg.
        cdum  = sqrt(cosb+img*sinb)
        cosb2 = real(cdum) ! cosb2 = cos(b/2)
        sinb2 = imag(cdum) ! sinb2 = sin(b/2)
      endif
      if(sinb>1d-8) then
        cosa =  rot1(1,3) / sinb
        sina =  rot1(2,3) / sinb
        cosc = -rot1(3,1) / sinb
        sinc =  rot1(3,2) / sinb
      else ! in this case we can set angle c to zero
        cosa = rot1(1,1) / cosb
        sina = rot1(2,1) / cosb
        cosc = 1d0
        sinc = 0d0
        sinb = 0d0
      endif
      if(abs(cosa)>1) cosa = sign(1d0,cosa) ! avoid rounding errors
      if(abs(sina)>1) sina = sign(1d0,sina) ! avoid rounding errors
      if(abs(cosc)>1) cosc = sign(1d0,cosc) ! avoid rounding errors
      if(abs(sinc)>1) sinc = sign(1d0,sinc) ! avoid rounding errors
      if( any( abs ( [ cosa,sina,cosb,sinb,cosc,sinc ] ) > 1 )) then ! The 1d-10 is a dirty fix!
        write(6,'(A)') 'dwigner: Defective symmetry transformation matrix.'
        write(6,'(A)') '         Not a proper rotation matrix:'
        write(6,'(9X,3F10.5)') rot1
        Error('Defective symmetry transformation matrix.')
      endif

      ! Calculate Wigner D-matrix elements
      do l=0,lcut
        do m1=-l,l
          do m2=-l,l
            nom  = sfac(l+m1) * sfac(l-m1) * sfac(l+m2) * sfac(l-m2)
            tmin = max(0,m1-m2)
            tmax = min(l+m1,l-m2)
            do t=tmin,tmax
              nexpo         = 2*l+m1-m2-2*t ; if(nexpo/=0) then ; cosb2_ = cosb2**nexpo ; else ; cosb2_ = 1 ; endif
              nexpo         = 2*t+m2-m1     ; if(nexpo/=0) then ; sinb2_ = sinb2**nexpo ; else ; sinb2_ = 1 ; endif
              denom         = fac(l+m1-t) * fac(l-m2-t) * fac(t) * fac(t+m2-m1)
              dwgn(m1,m2,l) = dwgn(m1,m2,l) + (-1)**t * nom / denom * cosb2_ * sinb2_
            enddo
            dwgn(m1,m2,l) = (cosa-img*sina)**m1 * dwgn(m1,m2,l) * (cosc-img*sinc)**m2 * nint(det)**l
            if(abs(dwgn(m1,m2,l))<1d-10) dwgn(m1,m2,l) = 0
          enddo
        enddo
      enddo
      end

c     --------------

c     Returns the coefficients H of
c
c     SUM G   Y  (r)
c      lm  lm  lm
c     -------------- = SUM H   Y  (r) .
c              *        lm  lm  lm
c     SUM L   Y  (r)
c      lm  lm  lm
c
c     The expansion becomes exact in the limit cuth -> infinity.
c     Note that the expression on the left-hand side must not have a pole. (Otherwise the result is meaningless.)
c
      subroutine expandfrac(hcoeff,gcoeff,lcoeff,cuth,cutg,cutl)
      use wrapper
      use, intrinsic :: iso_fortran_env
      implicit none
      integer,    intent(in)  :: cuth,cutg,cutl
      complex_dp, intent(in)  :: gcoeff((cutg+1)**2),lcoeff((cutl+1)**2)
      complex_dp, intent(out) :: hcoeff((cuth+1)**2)
      complex_dp              :: matrix((cuth+1)**2,(cuth+1)**2)
      integer                 :: ll,lh,lg,ml,mh,mg,lml,lmh,lmg
      real_dp                 :: gaunt
      hcoeff = 0 ; if(all(abs(real(gcoeff))+abs(imag(gcoeff))<1d-12)) return
      matrix = 0
      lml    = 0
      do ll = 0,cutl
        do ml = -ll,ll
          lml = lml + 1
          lmh = 0
          do lh = 0,cuth
            do mh = -lh,lh
              lmh = lmh + 1
              mg  = mh - ml
              do lg = max(abs(mg),abs(ll-lh)),min(cuth,ll+lh)
                lmg = lg**2 + (lg+mg) + 1
                matrix(lmg,lmh) = matrix(lmg,lmh) + gaunt(ll,lh,lg,ml,mh,mg) * lcoeff(lml)
              enddo
            enddo
          enddo
        enddo
      enddo
      call solve(hcoeff,matrix,gcoeff)
      end

c     ------------

c     Returns matrix A that fulfills  r^T*A*r = 1/(r^T*B*r) up to quadratic order in r.
c
c     Cases of bad convergence: It seems that matY(1) converges exponentially as a function of l, i.e., f(l) = a + b * exp(-c*x).
c     Given three value pairs (l,f1),(l+2,f2),(l+4,f3), we can derive c = log[(f2-f1)/(f3-f2)]/2 and
c     a = ( f2**2 - f1*f3 ) / ( 2*f2 - f1 - f3 ). a is the limit l->\infty.
c     boost: if defined, this extrapolated value (if l>=4) is used to replace matY(1).
# define boost
      subroutine invert_angular(matrix)
      use global, only: pi
      use util,   only: chr
      use, intrinsic :: iso_fortran_env
      implicit none
      complex_dp, intent(inout) :: matrix(3,3)
      complex_dp, allocatable   :: matY(:)
      complex_dp                :: denom(9),old,last(3),extra=0
      integer                   :: l
      if(abs(matrix(1,1))+abs(matrix(2,2))+abs(matrix(3,3))>huge(0d0)) then
        matrix = 0
        return
      endif
      call matrix2harmonics(denom(1),denom(5),matrix) ; denom(2:4) = 0
      call conjgharmonics(denom,2)
      old = -1000
      l   = 0
      do
        allocate ( matY((l+1)**2) )
        call expandfrac(matY,[(1d0,0d0)*sqrt(4*pi)],denom,l,0,2)
        last(1) = last(2)
        last(2) = last(3)
        last(3) = matY(1)
        if(l>=4) extra = ( last(2)**2 - last(1)*last(3) ) / ( 2*last(2) - last(1) - last(3) )
c        write(60,*) l,real(matY(1)),real(extra)
        if(abs(old-matY(1))<1d-6) then ; exit
        else if(l==20)            then ; Info('Inversion not converged with l=20.')
        else if(l==40)            then ; Warn('Convergence not reached with l=40. Last error: '//Chf(abs(old-matY(1)),'F10.6'))
                                         goto 1
        endif
        old = matY(1)
        deallocate ( matY )
        l = l + 2
      enddo
      if(l>20) Info('Inversion converged with l='//chr(l))
 1    continue
# ifdef boost
      if(l>=4) matY(1) = extra
# endif
      call harmonics2matrix(matrix,matY(1),matY(5))
      deallocate ( matY )
      end

c     ------------

c     Returns an interpolated function value and derivative (in y and dy) at x
c     of a function given by n value pairs (xx,yy) using spline interpolation
c     (or linear interpolation if LINEAR is defined).
c     Outside the range [x(1),x(n)] the function is extrapolated linearly.
c     x, xx : real valued ; y, yy : complex valued.
c# define LIN_INTERPOL
      recursive subroutine interpolate(y,dy,x,yy,xx,n)
      use, intrinsic :: iso_fortran_env
      implicit none
      integer,    intent(in)  :: n
      real_dp,    intent(in)  :: xx(n),x
      complex_dp, intent(in)  :: yy(n)
      complex_dp, intent(out) :: y,dy
      complex_dp              :: a(n-1),b(n-1)
      real_dp                 :: u,v
      integer                 :: i
      if(n<=1) Error('Frequency mesh contains a single point.')
      if(xx(1)>xx(n)) then
        call interpolate(y,dy,-x,yy,-xx,n)
        dy = -dy
        return
      endif
# ifdef LIN_INTERPOL
#   warning LIN_INTERPOL interpolation defined for testing
      if(x>=xx(n)) then
        i = n - 1
      else
        i = 1
        do while(xx(i+1)<x)
          i = i + 1
        enddo
      endif
      dy = (yy(i+1)-yy(i)) / (xx(i+1)-xx(i))
      y  = yy(i) + dy*(x-xx(i))
# else
      call cubicspline_c(a,b,yy, xx,n,(0d0,0d0),2,(0d0,0d0),2)
      if(x<=xx(1)) then
        dy = (b(1)-3*yy(1)) / (xx(2)-xx(1))
        y  = yy(1) + dy*(x-xx(1))
      else if(x>=xx(n)) then
        dy = (3*yy(n)-a(n-1)) / (xx(n)-xx(n-1))
        y  = yy(n) + dy*(x-xx(n))
      else
        i = 1
        do while(x>xx(i+1))
          i = i + 1
        enddo
        u  = (x-xx(i)) / (xx(i+1)-xx(i))
        v  = 1 - u
        y  = yy(i)*v**3 + yy(i+1)*u**3 + a(i)*u**2*v + b(i)*u*v**2
        dy = (b(i)-3*yy(i)) * v**2 + (3*yy(i+1)-a(i)) * u**2 + 2*(a(i)-b(i)) * u*v
        dy = dy / (xx(i+1)-xx(i))
      endif
# endif
      end

c     ------------

c
c     Same as interpolate, but uses precalculated a(:) and b(:) (as returned by routine cubicspline)
c     (Does not return gradient.)      
      subroutine interpolatef(a,b,y,x,yy,xx,n)
      use, intrinsic :: iso_fortran_env
      implicit none
      integer, intent(in)  :: n
      real_dp, intent(in)  :: a(n-1),b(n-1),xx(n),yy(n),x
      real_dp, intent(out) :: y
      real_dp              :: u,v
      integer              :: j
      j = 1
      do while(xx(j+1)<x)
        j = j + 1 ; if(j+1==n) exit
      enddo
      u = (x-xx(j))/(xx(j+1)-xx(j))
      v = 1 - u
      y = yy(j)*v**3 + yy(j+1)*u**3 + a(j)*u**2*v + b(j)*u*v**2
      end

c     ------------      

c     Returns the coefficients a(i) and b(i), i=1..n-1, of the interpolation s(x) with cubic splines of the set
c     {x(i),y(i); i=1..n}. Within the interval [x(i),x(i+1)] the interpolation is given by
c
c     s(x) = y(i)*v**3 + y(i+1)*u**3 + a(i)*u**2*v + b(i)*u*v**2
c
c     with u=(x-x(i))/(x(i+1)-x(i)) and v=1-u.
c
c     There are two free parameters, which are fixed with the last four arguments:
c
c     (d/dx)^m1 s(x(1)) = y1 and
c     (d/dx)^mn s(x(n)) = yn for m1 = 1,2,3,4 and mn = 1,2,3,
c     and (d/dw)^3 s(x(w(1))) = y1 with x(w) = w/(1+w) for m1 = 4.
c
      subroutine cubicspline(a,b,y,x,n,y1,m1,yn,mn)
      use, intrinsic :: iso_fortran_env
      implicit none
      integer, intent(in)  :: n,m1,mn
      real_dp, intent(in)  :: x(n)
      real_dp, intent(in)  :: y(n),y1,yn
      real_dp, intent(out) :: a(n-1),b(n-1)
      real_dp              :: h(n-1)
      real_dp              :: d(n),dl(n-1),du(n-1),p(n)
      integer              :: i
      if(n<=1) Error('at least two points are needed.')
      do i = 1,n-1
        if(x(i+1)<=x(i)) Bug('x array not ordered according to size.')
        h(i) = x(i+1) - x(i)
      enddo
      ! n-2 tridiagonal system of linear equations
      do i = 2,n-1
        d(i)    = 4/h(i) + 4/h(i-1)
        du(i)   = 2/h(i)
        dl(i-1) =          2/h(i-1)
        p(i)    = 6 * ( (y(i)-y(i-1)) / h(i-1)**2 + (y(i+1)-y(i)) / h(i)**2 )
      enddo
      ! imposing the additional boundary conditions
      if(m1==1) then
        d(1)  = 1
        du(1) = 0
        p(1)  = y1
      else if(m1==2) then
        d(1)  = 4/h(1)
        du(1) = 2/h(1)
        p(1)  = 6 * (y(2)-y(1)) / h(1)**2 - y1
      else if(m1==3) then
        d(1)  = 6/h(1)**2
        du(1) = 6/h(1)**2
        p(1)  = 12 * (y(2)-y(1)) / h(1)**3 + y1
      else if(m1==4) then
        d(1)  = 6/h(1)**2 + 24/h(1) + 6
        du(1) = 6/h(1)**2 + 12/h(1)
        p(1)  = 12*(3*h(1)+1) * (y(2)-y(1)) / h(1)**3 + y1
      else
        Bug('wrong m1 argument.')
      endif
      if(mn==1) then
        d(n)    = 1
        dl(n-1) = 0
        p(n)    = yn
      else if(mn==2) then
        d(n)    = 4/h(n-1)
        dl(n-1) = 2/h(n-1)
        p(n)    = 6 * (y(n)-y(n-1)) / h(n-1)**2 + yn
      else if(mn==3) then
        d(n)    = 6/h(n-1)**2
        dl(n-1) = 6/h(n-1)**2
        p(n)    = 12 * (y(n)-y(n-1)) / h(n-1)**3 + yn
      else
        Bug('wrong mn argument.')
      endif
      ! solve matrix equation
      call dgtsv(n,1,dl,d,du,p,n,i)
      if(i/=0) Error('dgtsv failed.')
      do i = 1,n-1
        a(i) = 3*y(i+1) - h(i)*p(i+1)
        b(i) = 3*y(i)   + h(i)*p(i)
      enddo
      end

c     ------------

c     Complex version of cubicspline.
      subroutine cubicspline_c(a,b,y,x,n,y1,m1,yn,mn)
      use global, only: img
      use, intrinsic :: iso_fortran_env
      implicit none
      integer,    intent(in)  :: n,m1,mn
      real_dp,    intent(in)  :: x(n)
      complex_dp, intent(in)  :: y(n),y1,yn
      complex_dp, intent(out) :: a(n-1),b(n-1)
      real_dp                 :: ar(n-1),br(n-1)
      real_dp                 :: ai(n-1),bi(n-1)
      call cubicspline(ar,br,real(y),x,n,real(y1),m1,real(yn),mn)
      call cubicspline(ai,bi,imag(y),x,n,imag(y1),m1,imag(yn),mn)
      a = ar + img*ai
      b = br + img*bi
      end

c     ------------

c     Transforms from the representation
c     s(x) = y1*v**3 + y2*u**3 + a*u**2*v + b*u*v**2, u=(x-x1)/(x2-x1), v=1-u
c     to
c     s(x) = c(1)*x**3 + c(2)*x**2 + c(3)*x + c(4)
      subroutine cubiccoeff(c,a,b,x,y)
      use, intrinsic :: iso_fortran_env
      implicit none
      real_dp, intent(out) :: c(0:3)
      real_dp, intent(in)  :: x(2)
      real_dp, intent(in)  :: a,b,y(2)
      real_dp                  :: dx
      dx   = ( x(2) - x(1) )**3
      c(0) = (     x(2)**3*y(1) -   x(1)**3*y(2) + x(1)**2*x(2)       * a - x(2)**2*x(1)       * b ) / dx
      c(1) = ( - 3*x(2)**2*y(1) + 3*x(1)**2*y(2) - x(1)*(2*x(2)+x(1)) * a + x(2)*(2*x(1)+x(2)) * b ) / dx
      c(2) = (   3*x(2)*   y(1) - 3*x(1)*   y(2) +      (2*x(1)+x(2)) * a -      (2*x(2)+x(1)) * b ) / dx
      c(3) = ( -           y(1) +           y(2) -                      a +                      b ) / dx
      end

      subroutine cubiccoeff_c(c,a,b,x,y)
      use, intrinsic :: iso_fortran_env
      implicit none
      complex_dp, intent(out) :: c(0:3)
      real_dp,    intent(in)  :: x(2)
      complex_dp, intent(in)  :: a,b,y(2)
      real_dp                  :: dx
      dx   = ( x(2) - x(1) )**3
      c(0) = (     x(2)**3*y(1) -   x(1)**3*y(2) + x(1)**2*x(2)       * a - x(2)**2*x(1)       * b ) / dx
      c(1) = ( - 3*x(2)**2*y(1) + 3*x(1)**2*y(2) - x(1)*(2*x(2)+x(1)) * a + x(2)*(2*x(1)+x(2)) * b ) / dx
      c(2) = (   3*x(2)*   y(1) - 3*x(1)*   y(2) +      (2*x(1)+x(2)) * a -      (2*x(2)+x(1)) * b ) / dx
      c(3) = ( -           y(1) +           y(2) -                      a +                      b ) / dx
      end

c     ------------

c     Returns an interpolated function value and derivative (in y and dy) at x
c     of a function given by n value pairs (xx,yy).
c     x, xx : real valued ; y, yy : complex valued.
c     The xx array must be ordered according to size.
      subroutine interpolate_old(y,dy,x,yy,xx,n) ! replaced by cubic interpolation now
      use, intrinsic :: iso_fortran_env
      implicit none
      integer,    intent(in)  :: n
      real_dp,    intent(in)  :: xx(n),x
      complex_dp, intent(in)  :: yy(n)
      complex_dp, intent(out) :: y,dy
      complex_dp              :: c(0:2),d(0:2)
      integer                 :: i
      real_dp                 :: qq,dqq,h
      real_dp                 :: q,dq
      q(h)  = -2*h**3 + 3*h**2 ! interpolates smoothly between 0 and 1
      dq(h) = -6*h**2 + 6*h
      do i = 2,n
        if(xx(i)<=xx(i-1)) Bug('xx array not ordered according to size.')
      enddo
      if(n==2) then
        dy = (yy(2)-yy(1)) / (xx(2)-xx(1))
        y  = yy(1) + dy * (x-xx(1))
      else if(n==3.or.x<xx(2)) then
        call polynomial(c,xx,yy)
        if(x>=xx(1)) then
          dy =        c(1)   + c(2)*x*2
          y  = c(0) + c(1)*x + c(2)*x**2
        else
          dy =        c(1)   + c(2)*xx(1)*2
          y  = yy(1) + dy * (x-xx(1))
        endif
      else if(x>xx(n-1)) then
        call polynomial(c,xx(n-2),yy(n-2))
        if(x<=xx(n)) then
          dy =        c(1)   + c(2)*x*2
          y  = c(0) + c(1)*x + c(2)*x**2
        else
          dy =        c(1)   + c(2)*xx(n)*2
          y  = yy(n) + dy * (x-xx(n))
        endif
      else
        do i = 3,n-1
          if(x<=xx(i)) then
            call polynomial(c,xx(i-2),yy(i-2))
            call polynomial(d,xx(i-1),yy(i-1))
            h   = (x-xx(i-1)) / (xx(i)-xx(i-1))
            qq  = q(h)
            dqq = dq(h) / (xx(i)-xx(i-1))
            y   = (1-qq) * ( c(0) + c(1)*x + c(2)*x**2 ) +
     &               qq  * ( d(0) + d(1)*x + d(2)*x**2 )
            dy  = (1-qq) * (        c(1)   + c(2)*x*2  ) - dqq * ( c(0) + c(1)*x + c(2)*x**2 ) +
     &               qq  * (        d(1)   + d(2)*x*2  ) + dqq * ( d(0) + d(1)*x + d(2)*x**2 )
            exit
          endif
        enddo
      endif

      contains
      subroutine polynomial(c,xx,yy)
      use, intrinsic :: iso_fortran_env
      implicit none
      complex_dp, intent(out) :: c(0:2)
      complex_dp, intent(in)  :: yy(3)
      real_dp,    intent(in)  :: xx(3)
      c(0) =    xx(2)*xx(3)  / ( (xx(1)-xx(2)) * (xx(1)-xx(3)) ) * yy(1) +
     &          xx(1)*xx(3)  / ( (xx(2)-xx(1)) * (xx(2)-xx(3)) ) * yy(2) +
     &          xx(1)*xx(2)  / ( (xx(3)-xx(1)) * (xx(3)-xx(2)) ) * yy(3)
      c(1) = - (xx(2)+xx(3)) / ( (xx(1)-xx(2)) * (xx(1)-xx(3)) ) * yy(1) -
     &         (xx(1)+xx(3)) / ( (xx(2)-xx(1)) * (xx(2)-xx(3)) ) * yy(2) -
     &         (xx(1)+xx(2)) / ( (xx(3)-xx(1)) * (xx(3)-xx(2)) ) * yy(3)
      c(2) =             1d0 / ( (xx(1)-xx(2)) * (xx(1)-xx(3)) ) * yy(1) +
     &                   1d0 / ( (xx(2)-xx(1)) * (xx(2)-xx(3)) ) * yy(2) +
     &                   1d0 / ( (xx(3)-xx(1)) * (xx(3)-xx(2)) ) * yy(3)
      end subroutine polynomial

      end

c     ------------

c
c     Returns dilogarithm of z (Spence's function)
c
c     Li2(z) = INT(1,z) ln(z) / (1-t) dt
c
c     Accuracy can be set with parameter "accur".
c     Algorithm i or number of iterations n can be returned with (z,i) and (z,n) as argument.
      function dilog(z)
      use global, only: pi,img
      use, intrinsic :: iso_fortran_env
      implicit none
      complex_dp             :: dilog
      complex_dp, intent(in) :: z
      complex_dp             :: zz(6),zq,series
      real_dp                :: az(6)
      real_dp,    parameter  :: pi26 = pi**2/6 , pi3 = pi/3 , accur = 1d-14
      complex_dp, parameter  :: q = (1+img*sqrt(3d0))/2 , qc = conjg(q)
      complex_dp, parameter  :: dilogq = pi**2/36 - img*1.01494160640965362502d0
      integer                :: i,n
      zz = [ z , 1-z , 1/z , 1/(1-z) , 1-1/z , 1/(1-1/z) ]
      az = abs(1-zz)
      i  = minloc(az,1)
      zq = z-q
      if  (abs(zq)<az(i)) then ; series = dilog_series2(zq,n) ; i = 7 ; zq = z
      else ; zq = conjg(z)-q
        if(abs(zq)<az(i)) then ; series = dilog_series2(zq,n) ; i = 8 ; zq = conjg(z)
        else                   ; series = dilog_series(zz(i),n)
        endif
      endif
      select case(i)
        case(1)      ; dilog =  series
        case(2)      ; dilog = -series + pi26               ; if(az(i)>accur) dilog = dilog - log(z)*log(1-z)
        case(3)      ; dilog = -series        - log(z)**2/2
        case(4)      ; dilog =  series + pi26 + log(1-z) * ( log(1-z)/2 - log(z) )
        case(5)      ; dilog =  series - pi26 - log(z)**2/2 ; if(az(i)>accur) dilog = dilog - log(z)*log(1-1/z)
        case(6)      ; dilog = -series - pi26 - log(z-1)**2/2
        case default ; dilog =  dilogq + pi3 * (pi3-img*log(1-zq)) - series ; if(i==8) dilog = conjg(dilog)
      end select

      contains

      ! Returns Li2(z) for |1-z|<=1. (i is the number of terms included.)
      function dilog_series(z,i)
      implicit none
      complex_dp              :: dilog_series
      complex_dp, intent(in)  :: z
      integer,    intent(out) :: i
      complex_dp              :: z1,sum,fac,add
      sum = 0
      z1  = 1-z ; if(abs(z1)>1) Error('Input parameter z does not fulfill |1-z|<=1.')
      fac = z1
      i   = 0
      do
        i   = i + 1
        add = fac / i**2
        sum = sum + add ; if(abs(real(add))+abs(imag(add))<accur) exit
        fac = fac * z1
      enddo
      dilog_series = sum
      end function dilog_series

      ! Returns auxiliary series SUM(m) Qm x**m with Qm = [ SUM(n=1,m-1) q**m/m ] / (m*conjg(q)**m). (i is the number of terms included.)
      function dilog_series2(x,i)
      implicit none
      complex_dp              :: dilog_series2
      complex_dp, intent(in)  :: x
      integer,    intent(out) :: i
      complex_dp              :: sum,xq,xqi,qi,add,qq
      if(abs(x)>1) Error('Input parameter x does not fulfill |x|<=1.')
      sum = 0
      qq  = 0
      xq  = x*q
      xqi = xq**2
      qi  = q
      i   = 1
      do
        qq  = qq + qi/i
        i   = i + 1
        add = qq/i * xqi
        sum = sum + add ; if(abs(real(add))+abs(imag(add))<accur) exit
        xqi = xqi * xq
        qi  = qi * q
      enddo
      dilog_series2 = sum
      end function dilog_series2

      end

c     ------------
      
c     Returns least common multiple of the integers iarr(1:n).
      function kgv(iarr,n)
      implicit none
      integer              :: kgv
      integer, intent(in)  :: n,iarr(n)
      logical              :: lprim(2:maxval(iarr))
      integer, allocatable :: prim(:),expo(:)
      integer              :: nprim,marr
      integer              :: i,j,ia,k
      ! Determine prime numbers
      marr  = maxval(iarr)
      lprim = .true.
      do i = 2,marr
        j = 2
        do while (i*j<=marr)
          lprim(i*j) = .false.
          j          = j + 1
        enddo
      enddo
      nprim = count(lprim)
      allocate ( prim(nprim),expo(nprim) )
      j = 0
      do i = 2,marr
        if(lprim(i)) then
          j       = j + 1
          prim(j) = i
        endif
      enddo
      ! Determine least common multiple
      expo = 0
      do i = 1,n
        ia = iarr(i)
        if(ia==0) cycle
        do j = 1,nprim
          k = 0
          do while(ia/prim(j)*prim(j)==ia)
            k  = k + 1
            ia = ia / prim(j)
          enddo
          expo(j) = max(expo(j),k)
        enddo
      enddo
      kgv = 1
      do j = 1,nprim
        kgv = kgv * prim(j)**expo(j)
      enddo
      deallocate ( prim,expo )
      end

c     ----------------

c     Returns prime factorization of |n0| in fac(:nfac)
c     n0>0 : prime factorization of  n0 in the form  n0 = fac(1) * fac(2) * ... * fac(nfac); fac(*) are in ascending order.
c     n0<0 : prime factorization of -n0 in the form -n0 = a(1)**b(1) * a(2)**b(2) * ... * a(nfac)**b(nfac); a(*) are in ascending order.
c            with fac(:2*nfac) = [ a(:nfac) , b(:nfac) ]; can be given as 2D array fac(:nfac,2) in the calling routine.
c     Dimension of fac must be large enough, i.e., size(fac) >=     ceiling( log(1d0*abs(n0)) / log(2d0) ) for n0>0
c                                              and size(fac) >= 2 * ceiling( log(1d0*abs(n0)) / log(2d0) ) for n0<0
c     Returns nfac=0 if |n0|=1.
      subroutine primfac(fac,nfac,n0)
      use util, only: chr
      implicit none
      integer, intent(out)  :: nfac,fac(*)
      integer, intent(in)   :: n0
      logical               :: lprim(2:abs(n0))
      integer, allocatable  :: prim(:),pow(:)
      integer               :: nprim
      integer               :: n,n1,i,j
      if(n0==0) Bug('Argument must not be zero.')
      n = abs(n0)
      ! Determine prime numbers
      lprim = .true.
      do i = 2,n
        j = 2
        do while (i*j<=n)
          lprim(i*j) = .false.
          j          = j + 1
        enddo
      enddo
      nprim = count(lprim)
      allocate ( prim(nprim) )
      j = 0
      do i = 2,n
        if(lprim(i)) then
          j       = j + 1
          prim(j) = i
        endif
      enddo
      ! Determine factorization
      n1 = n
      j  = 0
      do i = 1,nprim
        do while(mod(n1,prim(i))==0)
          j      = j + 1
          fac(j) = prim(i)
          n1     = n1 / prim(i)
        enddo
      enddo
      deallocate(prim)
      nfac = j
      if(n1/=1) Bug('Prime factorization could not be determined.')
      if(n0<0.and.nfac>0) then
        allocate(pow(nfac))
        pow = 1
        j   = 1
        do i = 2,nfac
          if(fac(i-1)==fac(i)) then ; pow(j) = pow(j) + 1
          else                      ; fac(j) = fac(i-1) ; j = j + 1
          endif
        enddo
        fac(j)       = fac(nfac)
        fac(j+1:2*j) = pow(:j)        
        nfac         = j
        deallocate(pow)
      endif
      end

c     ----------------

c     Returns greatest common divisor of n1 and n2
      function ggt(n1,n2)
      implicit none
      integer              :: ggt
      integer, intent(in)  :: n1,n2
      integer, allocatable :: p1(:),p2(:)
      integer              :: np1,np2,i1,i2
      allocate ( p1(ceiling(log(1d0*n1)/log(2d0))) )
      allocate ( p2(ceiling(log(1d0*n2)/log(2d0))) )
      call primfac(p1,np1,n1)
      call primfac(p2,np2,n2)
      ggt = 1
      do i1 = 1,np1
        do i2 = 1,np2
          if(p1(i1)==p2(i2)) then
            ggt    = ggt * p1(i1)
            p2(i2) = 0
            exit
          endif
        enddo
      enddo
      end

c     ----------------

c     Returns all integer products that yield "result" (e.g., 2*3*7=42).
c     result           : Result of products.
c     nprd             : Number of products found.
c     dim              : Number of integer factors.
c     dim<0            : Only return number of products nprd; prd is not referenced.
c     prd(:dim,:nprod) : Integer factors of all products (if dim>0).
c     
      subroutine iproducts(prd,nprd,dim,result)
      implicit none
      integer, intent(in)  :: dim,result
      integer, intent(out) :: prd(abs(dim),*)
      integer, intent(out) :: nprd
      integer              :: fac( 2 * ceiling( log(1d0*result) / log(2d0) ) )
      integer              :: n,p(abs(dim)),d
      logical              :: def
      def = dim>0
      d   = abs(dim)
      if(dim==0)    Bug('Dimension "dim" must not be zero.')
      if(result<=0) Bug('result must be positive.')
      if(result==1) then
        nprd = 1
        if(def) prd(:,1) = 1
        return
      endif
      call primfac(fac,n,-result)
      nprd = 0
      call iproducts1(p,prd,nprd,d,d,fac,fac(n+1),n,result,def)

      contains

      recursive subroutine iproducts1(p,prd,nprd,d,dim,fac,pow,n,result,def)
      implicit none
      logical, intent(in)    :: def
      integer, intent(in)    :: d,dim,n,fac(n),pow(n),result
      integer, intent(inout) :: p(dim),prd(dim,*),nprd
      integer, allocatable   :: pw(:)
      integer                :: nn,nn1,i
      nn = product([fac**pow])
      if(d==1) then
        p(1) = nn
        nprd = nprd + 1
        if(def) prd(:,nprd) = p
        if(product(p)/=result) Bug('Product incorrect.')
        return
      endif
      allocate(pw(n))
      pw = 0
      do
        nn1  = product([fac**pw])
        p(d) = nn1
        call iproducts1(p,prd,nprd,d-1,dim,fac,pow-pw,n,result,def)
        if(nn1==nn) exit
        call nextpw(pw,pow)
      enddo
      deallocate(pw)
      end subroutine iproducts1

      recursive subroutine nextpw(pw,pow)
      implicit none
      integer, intent(inout) :: pw(*)
      integer, intent(in)    :: pow(*)
      pw(1) = pw(1) + 1
      if(pw(1)>pow(1)) then
        pw(1) = 0
        call nextpw(pw(2),pow(2))
      endif
      end subroutine nextpw

      end
