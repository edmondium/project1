c--------------------------------------------------------------------------------
c Copyright (c) 2021 Peter Grünberg Institut, Forschungszentrum Jülich, Germany
c This file is part of SPEX and available as free software under the conditions
c of the MIT license as expressed in the LICENSE file in more detail.
c--------------------------------------------------------------------------------

c     Calculates the coefficients cpade(:) for the Pade interpolation pade(x) (Thiele's algorithm).
c
c     If n odd, an additional constraint must be imposed. Here,
c     pade(x)*x**2 = finite for x->infinity  if constr = -1
c     no constraint                          if constr =  0,
c     Re pade'(x(n)) = 0                     if constr =  1,
c     Im pade'(x(n)) = 0                     if constr =  2,
c        pade'(x(n)) = 0                     if constr =  3.

# include "cppmacro.h"

      recursive subroutine pade_init(cpade,x,y,n,constr)
      use global, only: img
      use, intrinsic :: iso_fortran_env
      implicit none
      integer,    intent(in)  :: n,constr
      complex_dp, intent(out) :: cpade(*) ! dimension n (n+1) if n even (odd)
      complex_dp, intent(in)  :: x(n),y(n)
      complex_dp              :: dpade,dpade_func,pade(n)
      integer                 :: i,j
      if(all(constr/=[-1,0,1,2,3])) Bug('constr must be -1, 0, 1, 2, or 3.')
      if(all(abs(real(y))<1d-10).and.all(abs(imag(y))<1d-10)) then
        cpade(1) = huge(1d0)
        return
      endif
      if(mod(n,2)/=0.and.constr/=-1.and.constr/=3) then
        call pade_init(cpade,x(2),y(2),n-1,0)
        dpade = dpade_func(x(n),x(2),cpade,n-1)
      endif
      do j = 1,n
        cpade(j) = 1 / y(j)
        do i = 1,j-1
          cpade(j) = ( x(j) - x(i) ) / ( cpade(j) - cpade(i) )
        enddo
      enddo
      if(mod(n,2)/=0) then
        if(constr==-1) then
          cpade(n+1) = -sum(cpade(2:n-1:2))
        else
          if(x(n)/=0) Bug('Check code. Arrays x(:)/y(:) must probably be reversed.')
          if     (constr==1) then ; dpade = imag(dpade)*img
          else if(constr==2) then ; dpade = real(dpade)
          else if(constr==3) then ; dpade = 0d0
          endif
          pade(n) = cpade(n)
          do i = n-1,1,-1
            pade(i) = cpade(i) + (x(n)-x(i)) / pade(i+1)
          enddo
          dpade = -dpade * pade(1)**2
          do i = 1,n-1
            dpade = ( dpade - 1/pade(i+1) ) / (x(i)-x(n)) * pade(i+1)**2
          enddo
          cpade(n+1) = 1/dpade
        endif
      endif
      end

c     ------------

      function pade_inter(z,pole,resid,npol)
      use, intrinsic :: iso_fortran_env
      implicit none
      integer,    intent(in) :: npol
      complex_dp             :: pade_inter
      complex_dp, intent(in) :: z,pole(npol),resid(npol)
      integer                :: i
      pade_inter = 0
      do i = 1,npol
        pade_inter = pade_inter + resid(i)/(z-pole(i))
      enddo
      end function

c     ------------

c     Returns the interpolated value using the coefficients cpade(:), i.e., the continued fraction
c
c     f(z) = 1 / cpade(1) + (z-x(1)) / cpade(2) + (z-x(2)) / cpade(3) + ... + (z-x(n-1)) / cpade(n).
      function pade_func(z,x,cpade,n)
      use, intrinsic :: iso_fortran_env
      implicit none
      integer,    intent(in) :: n
      complex_dp             :: pade_func,pade
      complex_dp, intent(in) :: z,x(n),cpade(*)
      integer                :: i,nn
      if(cpade(1)==huge(1d0)) then
        pade_func = 0
        return
      endif
      if(mod(n,2)==0) then ; nn = n
      else                 ; nn = n + 1
      endif
      pade = 0
      do i = nn,2,-1
        pade = (z-x(i-1)) / (cpade(i)+pade)
      enddo
      pade_func = 1 / (cpade(1)+pade)
      end

c     ------------

c     Returns the derivative.
      function dpade_func(z,x,cpade,n)
      use, intrinsic :: iso_fortran_env
      implicit none
      integer,    intent(in) :: n
      complex_dp, intent(in) :: z,x(n),cpade(*)
      complex_dp             :: dpade_func,pade,dpade,den
      integer                :: i,nn
      if(cpade(1)==huge(1d0)) then
        dpade_func = 0
        return
      endif
      if(mod(n,2)==0) then ; nn = n
      else                 ; nn = n + 1
      endif
      pade  = 0
      dpade = 0
      do i = nn,2,-1
        den   = (cpade(i)+pade)
        pade  = (z-x(i-1)) / den
        dpade = ( 1 - pade*dpade ) / den
      enddo
      pade       = cpade(1) + pade
      dpade_func = -dpade/pade**2
      end

c     ------------

c     Returns the dominant poles of the Pade approximation defined as a continued fraction in cpade(:), i.e.,
c
c     f(z) ~= SUM(i=1,npole) resid(i) / (z-pole(i)).
c
c     f(iz) is real valued (lreal): then, if z is a pole with residue r, -conjg(z) is a pole, too, with residue -conjg(r).
      recursive subroutine pade_poles(pole,resid,npole,x,y,cpade,n,nsmooth,constr,lwrite)
      use util, only: chr
      use, intrinsic :: iso_fortran_env
      implicit none
      integer,    intent(in)  :: n,nsmooth,constr
      integer,    intent(out) :: npole
      logical,    intent(in)  :: lwrite
      complex_dp, intent(in)  :: x(n),y(n),cpade(*)
      complex_dp, intent(out) :: pole(*),resid(*)
      complex_dp, allocatable :: x_rand(:)
      complex_dp              :: cpad(n+1)
      real_dp,    parameter   :: minstep   = 1d-12 ! minimum step size (convergence criterion)
      real_dp,    parameter   :: error0    = 1d-6  ! maximum permissible error
      real_dp,    parameter   :: ransmooth = 1d-8  ! random noise range for smoothing
      real_dp,    parameter   :: ranrepeat = 1d-10 ! random noise range for repeating
      integer,    parameter   :: maxrepeat = 9     ! maximum number of repeated attempts with random noise
      complex_dp              :: z,dpade,pade,den,sum,dsum,step
      integer                 :: i,j,nn,istep,irand,npol
      logical                 :: purify,lreal,isnan
      real_dp                 :: r,error
      integer,    save        :: nrepeat = 0
      if(cpade(1)==huge(1d0)) then
        npole = 0
        return
      endif
      if(nsmooth>1) then
        call pade_poles(pole,resid,npole,x,y,cpade,n,1,constr,.false.)
        allocate(x_rand(n))
        do i = 2,nsmooth
          do j = 1,n
            call random_number(r)
            r         = (r-0.5d0) * ransmooth
            x_rand(j) = x(j) + r * imag(x(j)) * (0d0,1d0)
          enddo
          call pade_init(cpad,x_rand,y,n,constr)
          call pade_poles(pole(npole+1),resid(npole+1),npol,x_rand,y,cpad,n,1,constr,.false.)
          npole = npole + npol
        enddo
        resid(:npole) = resid(:npole) / nsmooth
        deallocate(x_rand)
        return
      endif
      if(mod(n,2)==0) then ; nn = n
      else                 ; nn = n + 1
      endif
      if(lwrite) then
        write(6,'(6X,A)')                  'Dominant poles'
        write(6,'(25X,A,25X,A,5X,A,3X,A)') 'position','weight','error','st'
      endif
      lreal  = all(imag(y)==0)
      purify = .true.
      irand  = 0
      npole  = 0
      do while(npole<nn/2)
        call random_number(r) ; z = r
        call random_number(r) ; z = z + (0d0,1d0) * r
        z     = ( z - (0.5d0,0.5d0) ) * 2 ; if(imag(x(2))*imag(z)>0) z = conjg(z)
        istep = 0
 1      pade  = 0
        dpade = 0
        do i = nn,2,-1
          den   = cpade(i) + pade
          pade  = (z-x(i-1)) / den
          dpade = ( 1 - pade*dpade ) / den
        enddo
        resid(npole+1) = 1/dpade ! Defined here before subtraction of poles leads to additional numerical inaccuracy.
        pade = pade + cpade(1)
        if(.not.purify) then
          sum  = 0
          dsum = 0
          do i = 1,npole
            den  = 1 / (z-pole(i))
            sum  = sum  + resid(i) * den
            dsum = dsum - resid(i) * den**2
          enddo
          den   = 1 / ( 1 - pade*sum )
          dpade = ( dpade + pade**2*dsum ) * den**2
          pade  = pade * den
        endif
        step  = pade / dpade / (1+istep/100) ! division by a factor to force convergence in problematic cases.
        z     = z - step
        istep = istep + 1
        if(isnan(real(z))) then
          if(lwrite) write(6,'(6X,A)') 'Detected NAN. Search stopped.'
          exit
        endif
        if(any(abs(z-pole(:npole))<1d-8)) then
          if(error>error0) then
            istep = 10001 ! Trigger restart with random start point
          else
            if(lwrite) write(6,'(6X,A)') 'Recurring pole. Search stopped.'
            exit
          endif
        endif
        if(istep>10000) then
          irand = irand + 1
          if(irand>20) then
            if(lwrite) write(6,'(6X,A)') 'Not converged after 20 restarts with random points.'
            exit
          endif
          call random_number(r) ; z = r
          call random_number(r) ; z = z + (0d0,1d0) * r
          z      = ( z - (0.5d0,0.5d0) ) * 2 ; if(imag(x(2))*imag(z)>0) z = conjg(z)
          purify = .false.
          istep  = 0
          goto 1
        endif
        if(purify) then ; if(abs(step)>minstep)           goto 1
        else            ; if(abs(step)>max(1d-8,minstep)) goto 1
        endif
        if(.not.purify) then
          purify = .true.
          goto 1
        endif
        purify      = .false.
        npole       = npole + 1 ; if(npole>nn/2) Bug('More poles than expected.')
        pole(npole) = z
        if(lreal) then
          if(abs(real(pole(npole)))>1d-8) then
            npole        = npole + 1
            pole(npole)  = -conjg(z)
            resid(npole) = -conjg(resid(npole-1))
          else
            if(abs(real(resid(npole)))>1d-8)
     &        Error('large residue real part for pole on imaginary axis: '//Chf(real(resid(npole)),'ES8.1'))
            resid(npole) = (0d0,1d0) * imag(resid(npole))
          endif
        endif
        error = 0
        do i = 1,n
          pade = 0
          do j = 1,npole
            pade = pade + resid(j) / (x(i)-pole(j))
          enddo
          error = error + abs(y(i)-pade)
        enddo
        if(lwrite) write(6,'(I3,2F15.10,'' '',2F15.10,''  '',ES8.1,I5)') npole,pole(npole),resid(npole),error,istep
        if(error<1d-11.and.npole/=nn/2) then
          if(lwrite) write(6,'(6X,A)') 'Error below 1e-12. Search stopped.'
          exit
        endif
      enddo
      if(error>error0) then
        nrepeat = nrepeat + 1
        if(nrepeat>maxrepeat) then
          write(0,'(A)') 'Could not determine poles of Pade approximant. Current data (with added noise):'
          write(0,'(39X,A,39X,A)') 'x','y'
          do i = 1,n
            write(0,'(2F20.13,2ES20.10)') x(i),y(i)
          enddo
          Error('Large error ('//Chf(error,'ES8.1')//') in pole detection. Giving up after '//Chr(nrepeat)//' attempts.' )
        else
          if(nrepeat>1) then
            Warn('Large error ('//Chf(error,'ES8.1')//') in pole detection (attempt '//Chr(nrepeat))
     &                                                                //'). Attempt with random noise ...'
          else
            Warn('Large error ('//Chf(error,'ES8.1')//') in pole detection. Attempt with random noise ...')
          endif
        endif
        allocate(x_rand(n))
        do j = 1,n
          call random_number(r)
          r         = (r-0.5d0) * ranrepeat
          x_rand(j) = x(j) + r * imag(x(j)) * (0d0,1d0)
        enddo
        call pade_init(cpad,x_rand,y,n,constr)
        call pade_poles(pole,resid,npol,x_rand,y,cpad,n,1,constr,.false.)
        deallocate(x_rand)
      endif
      nrepeat = 0
      end

c     ------------

# if 0

c     Returns the polynomial coefficients of the numerator and denominator, i.e.,
c
c     f(z) = ( num(1) + num(2)*z + ... ) / ( den(1) + den(2) * z + ... ).
c
c     The algorithm is unstable for larger n.
      subroutine pade_coeff(num,den,x,y,cpade,n)
      use, intrinsic :: iso_fortran_env
      implicit none
      integer,    intent(in)  :: n
      complex_dp, intent(in)  :: x(n),y(n),cpade(n)
      complex_dp, intent(out) :: num(n/2),den(n/2+1)
      complex_dp              :: g(n,n),a(n)
      integer                 :: i
      if(mod(n,2)/=0) Error('Not implemented for odd number of mesh points.')
      g(:,1) = y
      do i = 2,n
        g(:,i) = ( g(i-1,i-1) - g(:,i-1) ) / ( (x-x(i-1)) * g(:,i-1) )
      enddo
      a      = [ (g(i,i),i=1,n) ]
      g      = 0
      g(1,1) = a(1)
      g(1,2) = a(1)
      do i = 3,n
        g( :,i) = g(:,i-1) - x(i-1) * a(i) * g(:   ,i-2)
        g(2:,i) = g(2:,i)  +          a(i) * g(:n-1,i-2)
      enddo
      num    = g(:n/2,n)
      g      = 0
      g(1,1) = 1
      g(1,2) = 1 - x(1) * a(2)
      g(2,2) =            a(2)
      do i = 3,n
        g( :,i) = g(:,i-1) - x(i-1) * a(i) * g(:   ,i-2)
        g(2:,i) = g(2:,i)  +          a(i) * g(:n-1,i-2)
      enddo
      den = g(:n/2+1,n)
      end

# endif




# if 0
      subroutine pade_init_old(cpade,x,y,n)
      use, intrinsic :: iso_fortran_env
      implicit none
      integer,    intent(in)  :: n
      complex_dp, intent(out) :: cpade(n)
      complex_dp, intent(in)  :: x(n),y(n)
      complex_dp              :: rho(n,-1:n),cpade1(n)
      integer                 :: istride,i1,i2
      integer i,j
      character(10) line
      if(mod(n,2)/=0) Error('Number of frequencies must be even.')
      rho = 0
      do istride = 0,n-1
        do i1 = 1,n-istride
          i2 = i1 + istride
          if(i1==i2) then
            rho(i1,i1) = y(i1)**(-1)
          else
            rho(i1,i2) = (x(i1)-x(i2)) / (rho(i1,i2-1)-rho(i1+1,i2)) + rho(i1+1,i2-1)
          endif
        enddo
      enddo
      do istride = 1,n
        cpade(istride) = rho(1,istride) - rho(1,istride-2)
      enddo
      do j = 1,n
        cpade1(j) = 1 / y(j)
        do i = 1,j-1
          cpade1(j) = ( x(j) - x(i) ) / ( cpade1(j) - cpade1(i) )
        enddo
      enddo
      cpade = cpade1
c      read(*,'(A)') line
c      if(line=='a') cpade = cpade1
c      if(i==1) cpade = cpade1
c      do i = 1,n
c        write(*,*) i
c        write(*,'(2F30.15)') cpade(i)
c        write(*,'(2F30.15)') cpade1(i)
c      enddo
c      read(*,*)
      end

      subroutine pade_poles1(pole,resid,x,y,cpade,n)
      use, intrinsic :: iso_fortran_env
      implicit none
      integer,    intent(in)  :: n
      complex_dp, intent(in)  :: x(n),y(n),cpade(n)
      complex_dp, intent(out) :: pole(n/2),resid(n/2)
      complex_dp              :: z,dpade,pade,den,sum,dsum,pade_func,step,dpad
      integer                 :: i,ipole,j
      do ipole = 1,n/2
        z     = 0
 1      pade  = 0
        dpade = 0
        do i = n,2,-1
          den   = 1 + pade
          pade  = cpade(i) * (z-x(i-1)) / den
          dpade = ( cpade(i) - pade*dpade ) / den
        enddo
        dpade = dpade / cpade(1)
        pade = (1+pade) / cpade(1)
        write(*,*)
        write(*,*) pade,dpade
        sum  = 0
        dsum = 0
        do i = 1,ipole-1
          den  = 1 / (z-pole(i))
          sum  = sum  + resid(i) * den
          dsum = dsum - resid(i) * den**2
        enddo
        dpad  = dpade
        dpade = dpade / pade**2
        pade  = ( pade**(-1) - sum )**(-1)
        dpade = (dsum+dpade) * pade**2
        step  = pade / dpade
        z     = z - step
        write(*,*) z
        if(abs(step)>1d-15) goto 1
        pole(ipole) = z
        resid(ipole) = 1/dpad
        write(*,*) 'Found pole ',pole(ipole)
        write(*,*) 'with residue ',resid(ipole)
        dpade = 0
        do i = 1,n
          pade = 0
          do j = 1,ipole
            pade = pade + resid(j) / (x(i)-pole(j))
          enddo
c          write(*,'(5F20.15)') y(i),pade,abs(y(i)-pade)
          dpade = dpade + abs(y(i)-pade)
        enddo
        write(*,*) ipole,dpade
        do i = 1,ipole
          write(*,'(4F20.15)') pole(i),resid(i)
        enddo
        read(*,*)
      enddo
      end
c#########
      subroutine pade_poles2(pole,resid,npole,x,y,cpade,nn,lwrite,lpur)
      use, intrinsic :: iso_fortran_env
      implicit none
      integer,    intent(in)  :: nn
      integer,    intent(out) :: npole
      logical,    intent(in)  :: lwrite,lpur
      complex_dp, intent(in)  :: x(nn),y(nn),cpade(nn)
      complex_dp, intent(out) :: pole(nn/2),resid(nn/2)
      complex_dp              :: c(nn),c1(nn)
      complex_dp              :: z,dpade,pade,den,sum,dsum,step
      integer                 :: i,j,n,istep,irand
      real_dp                 :: error,r
      logical :: purify
      if(lwrite) then
        write(6,'(6X,A)')        'Dominant poles'
        write(6,'(35X,A,36X,A)') 'position','weight'
      endif
      purify = .true.
      n     = nn
      c     = cpade
      irand = 0
      npole = 0
      do while(npole<nn/2)
        z     = 0
        istep = 0
 1      pade  = 0
        dpade = 0
        do i = n,2,-1
          den   = c(i) + pade
          pade  = (z-x(i-1)) / den
          dpade = ( 1 - pade*dpade ) / den
        enddo
        pade = pade + c(1)
c        sum  = 0
c        dsum = 0
c        do i = 1,npole
c          den  = 1 / (z-pole(i))
c          sum  = sum  + resid(i) * den
c          dsum = dsum - resid(i) * den**2
c        enddo
c        den   = 1 / ( 1 - pade*sum )
c        dpade = ( dpade + pade**2*dsum ) * den**2
c        pade  = pade * den
        step  = pade / dpade / (1+istep/100) ! division by a factor to force convergence in problematic cases.
        z     = z - step
        istep = istep + 1
        if(isnan(real(z))) then
          if(lwrite) write(6,'(6X,A)') 'Detected NAN. Search stopped'
          exit
        endif
        if(any(abs(z-pole(:npole))<1d-13)) then
          write(6,'(6X,A)') 'Recurring pole. Search stopped.'
          exit
        endif
        if(istep>1000) then
          irand = irand + 1
          if(irand>20) then
            write(6,'(6X,A)') 'Not converged after 20 restarts with random points.'
            exit
          endif
          call random_number(r) ; z = r
          call random_number(r) ; z = z + (0d0,1d0) * r
          istep = 0
          goto 1
        endif
        if(abs(step)>1d-15) goto 1
        if(lpur.and..not.purify) then
          c = cpade
          write(*,'(A$)') 'p'
          purify = .true.
          goto 1
        endif
        purify = .false.
        npole        = npole + 1
        pole(npole)  = z
        resid(npole) = 1/dpade
        error = 0
        c1 = 0
        do i = 1,nn
          pade = 0
          do j = 1,npole
            pade = pade + resid(j) / (x(i)-pole(j))
          enddo
          c1(i) = y(i) - pade
          error = error + abs(y(i)-pade)
c          write(*,'(I3,4F20.10)') i,y(i),c1(i)
        enddo
        if(lwrite) write(6,'(I3,2F20.15,''  '',3F20.15,I3)') npole,pole(npole),resid(npole),error,istep
        if(error<1d-12) then
c          write(6,'(6X,A)') 'Error below 1e-12. Search stopped.'
c          exit
        endif
        if(n>2) then
c          n = n - 2
          call pade_init(c,x,c1,n)
        endif
c        do j = 1,n-2
c          pade = 0
c          do i = j,2,-1
c            pade = (x(j)-x(i-1)) / (c(i)+pade)
c          enddo
c          pade  = 1 / (c(1)+pade) - resid(npole) / (x(j)-pole(npole))
c          c1(j) = 1 / pade
c          do i = 1,j-1
c            c1(j) = ( x(j) - x(i) ) / ( c1(j) - c1(i) )
c          enddo
c        enddo
c        write(*,*) c(:n)
c        c = c1
c        n = n - 2
c        write(*,*)
c        write(*,*) c1(:n)
c        read(*,*)
c        do j = 1,nn
c          pade = 0
c          do i = n,2,-1
c            pade = (x(j)-x(i-1)) / (c(i)+pade)
c          enddo
c          pade = c(1) + pade
c          write(*,*) x(j),1/pade,y(j)-resid(npole)/(x(j)-pole(npole))
c          write(*,*) x(j),1/pade-(y(j)-resid(npole)/(x(j)-pole(npole)))
c          read(*,*)
c        enddo
      enddo
c      if(lwrite)
      write(6,'(6X,A,F21.15)') 'Final deviation:',error
      if(error>1d-6) Error('Large deviation in pole detection.')
      if(error>1d-10) Warn('Large deviation in pole detection.')
c      read(*,*)
      end

      subroutine pade_init1(cpade,x,y,n)
      use, intrinsic :: iso_fortran_env
      implicit none
      integer,    intent(in)  :: n
      complex_dp, intent(in)  :: x(n),y(n)
      complex_dp, intent(out) :: cpade(n)
      complex_dp              :: g(n,n)
      integer                 :: i
      if(mod(n,2)/=0) Error('Number of frequencies must be even.')
      g(:,1) = y
      do i = 2,n
        g(:,i) = ( g(i-1,i-1) - g(:,i-1) ) / ( (x-x(i-1)) * g(:,i-1) )
      enddo
      cpade = [ (g(i,i),i=1,n) ]
      end
# endif


