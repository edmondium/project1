c--------------------------------------------------------------------------------
c Copyright (c) 2021 Peter Grünberg Institut, Forschungszentrum Jülich, Germany
c This file is part of SPEX and available as free software under the conditions
c of the MIT license as expressed in the LICENSE file in more detail.
c--------------------------------------------------------------------------------

c     Hilbert mesh

# include "cppmacro.h"

c
c     Get Hilbert-mesh definition (fspec) from spex.inp (HILBERT or FSPEC) or auto-adjusted.
      subroutine getfspec(fspec,nfrqh,nfrqc,frqc,lwrite)
      use global, only: inp,metal,egap
      use key
      use, intrinsic :: iso_fortran_env
      implicit none
      real_dp,       intent(out) :: fspec(2)
      integer,       intent(in)  :: nfrqh,nfrqc
      complex_dp,    intent(in)  :: frqc(nfrqc)
      logical,       intent(in)  :: lwrite
      character(:),  allocatable :: arg(:)
      real_dp                    :: f2,d
      real_dp                    :: ch2r
      integer                    :: n
      logical                    :: isinteger
      fspec = 0
      allocate ( character(30) :: arg(2) )
      call getkey(inp,'FSPEC', arg, section='SUSCEP', default=[' ',' '])
      if(arg(1)/=' ') then
        fspec(1) = ch2r(arg(1),'FSPEC') ; if(fspec(1)<=0) Error('First argument to FSPEC must be positive.')
        fspec(2) = ch2r(arg(2),'FSPEC')
        fspec(1) = -fspec(1) ! minus sign used as flag for FSPEC (backwards compatibility)
        if(fspec(2)<1) Error('Second argument to FSPEC must be greater than one.')
        if(lwrite) Info('Keyword FSPEC obsolete; use HILBERT.')
        call getkey(inp,'HILBERT', arg, section='SUSCEP', default=[' ',' '])
        if(arg(1)/=' ') Error('Keywords HILBERT and FSPEC given at the same time.')
      else ! new keyword HILBERT
        deallocate(arg)
        call getkey(inp,'HILBERT', arg, section='SUSCEP',  default=[' '] )
        if(arg(1)=='NONE') then   ! HILBERT NONE
          if(size(arg)>1) Error('Extra argument(s) after HILBERT NONE.')
        else if(arg(1)/=' ') then ! HILBERT explicit input
          if(size(arg)/=2) Error('Expected two arguments after HILBERT.')
          fspec(1) = ch2r(arg(1),'HILBERT') ; if(fspec(1)<0) Error('First HILBERT argument negative.')
          fspec(2) = ch2r(arg(2),'HILBERT') ; if(fspec(2)<1) Error('Second HILBERT argument smaller than one.')
          if(isinteger(arg(1))) then ; call      conversion ; if(lwrite) write(6,'(A,F13.8,F10.5)') 'Converted to HILBERT',fspec(:2)
          else                       ; call back_conversion ; if(lwrite) write(6,'(A,I5,F13.5)')    'Corresponds to HILBERT',n,f2
          endif
        else                      ! HILBERT default
          if(nfrqh==0.and.nfrqc>1.and.all(imag(frqc)==0d0).and.count(abs(real(frqc))>=egap)>2) then ! define according to frqc real mesh
            !n  = minloc(abs(real(frqc)),1,abs(real(frqc))>=egap) ; if(n==0) Bug('Error in minloc.')
            n  = maxloc(abs(real(frqc)),1,abs(real(frqc))>=egap) ; if(n==0) Bug('Error in maxloc.')
            f2 = frqc(n) - egap
            if(n-1<1) then ; d = abs( frqc(n+1) - frqc(n) )
            else           ; d = abs( frqc(n-1) - frqc(n) )
            endif
            fspec(1) = min(d,0.01d0)     ! assume frqs-increment to be min(d,0.01) at egap (or 0) and max(5*fspec(1),d) at f2
            fspec(2) = max(d,5*fspec(1)) !
            fspec(2) = ( f2 - fspec(1) ) / ( f2 - fspec(2) ) ! solution x of: d + dx + dx**2 + ... + dx**n = f2   and  D = dx**n  with d=fspec(1), D=fspec(2)
            if(lwrite) write(6,'(A,F13.8,F10.5)') 'Auto-adjusted  HILBERT',fspec(:2)
            if(fspec(1)<=0.or.fspec(2)<1) Error('Auto-adjusted parameters out of reasonable range.')
          else
            if(metal) then ; fspec = [ 150,50 ]
            else           ; fspec = [  50,30 ]
            endif
            call conversion
            write(6,'(A,F13.8,F10.5)') 'Default mesh   HILBERT',fspec(:2)
          endif
          call back_conversion ; if(lwrite) write(6,'(A,I5,F13.5)') 'Corresponds to HILBERT',n,f2
        endif
      endif
      deallocate ( arg )

      contains
      subroutine conversion ! fspec -> fspec
      implicit none
      integer  :: n
      n        = fspec(1)
      fspec(2) = fspec(2)**(1/(n+0.5d0))
      if(abs(1-fspec(2))<1d-6) then ; fspec(1) = 5d0 / n
      else                          ; fspec(1) = 5d0 / (1-fspec(2)**n) * (1-fspec(2))
      endif
      end subroutine conversion
      subroutine back_conversion ! fspec -> n,f2
      implicit none
      if(fspec(2)==1) then ; n = 5/fspec(1)
      else                 ; n = log(1-5*(1-fspec(2))/fspec(1)) / log(fspec(2))
      endif
      f2 = fspec(2)**(n+0.5d0)
      end subroutine back_conversion
      end

c     ---------------

c
c     Define Hilbert mesh frqs(1:nfrqs) in the frequency range wmin..wmax
c
c     frqs(1) , frqs(2)=frqs(1)+winc , frqs(3)=frqs(2)+winc*h , frqs(4)=frqs(3)+winc*h*h , ... , frqs(nfrqs) > wmax-a*[frqs(nfrqs)-frqs(nfrqs-1)]
c
c     frqs(1)=wmin+a*[frqs(2)-frqs(1)] (so that the BZ integration region starts at wmin)
c
c     with winc = fspec(1), h = fspec(2), a = 1/(1+h)
c
c begin interface
      subroutine getfrqs(frqs,nfrqs,fspec,wmin,wmax,lwrite)
      use util, only: chr
      use file
      use key
      use, intrinsic :: iso_fortran_env !inc
      implicit none
      integer, intent(out)              :: nfrqs
      real_dp, intent(out), allocatable :: frqs(:)
      real_dp, intent(in)               :: wmin,wmax,fspec(2)
      logical, intent(in)               :: lwrite
c end interface
      logical                           :: define
      integer                           :: i,n
      real_dp                           :: winc,wdif,wmin1,winc1,a
      real_dp                           :: h,rdum

      define = .false.
      nfrqs  = 0

c     Define mesh frqs(1:N) for Hilbert trafo (integration over spectral function).
      n = 0
      if(fspec(1)/=0) then
        ! Hilbert mesh defined:
        ! Define first stride winc and factorial increment h [ w_i+1-w_i = h * (w_i-w_i-1) ]
        wmin1 = wmin
        if(fspec(1)<0) wmin1 = 0 ! Backwards compatibility: FSPEC
        wdif = wmax - wmin1
        winc = abs(fspec(1))
        h    =     fspec(2)
        if(fspec(1)<0) then ! Backwards compatibility: FSPEC
          if(lwrite) Info('FSPEC used. Hilbert mesh starts at zero.')
          n = nint ( -fspec(1) )
          if(abs(n+fspec(1))<1d-12) then
            if(h==1) then ; winc = wdif / (n-1)
            else          ; winc = wdif * (h-1) / (h**(n-1)-1)
            endif
            if(lwrite) Info('First FSPEC argument interpreted as integer. First increment is '//Chf(winc,'F14.10'))
          else
            if(lwrite) Info('First FSPEC argument interpreted as real.')
          endif
        endif
        ! Define Hilbert mesh
        a     = 1/(1+h)
 1      rdum  = wmin1 + winc * a ; if(fspec(1)<0) rdum = 0 ! Backwards compatibility: FSPEC
        winc1 = winc
        i     = 0
        do while(rdum-winc1*a<wmax)
          i     = i + 1 ; if(define) frqs(i) = rdum
          rdum  = rdum + winc1
          winc1 = winc1 * h
        enddo
        if(.not.define) then
          define = .true.
          nfrqs  = i
          allocate(frqs(nfrqs))
          goto 1
        endif
        ! Test Hilbert mesh
        do i = 1,nfrqs-1
          if(frqs(i)<0)          Error('Detected an error in the spectral-function mesh (negative value).')
          if(frqs(i)>=frqs(i+1)) Error('Detected an error in the spectral-function mesh (not ordered).')
          if(i>1) then
            if(abs(a-(frqs(i)-frqs(i-1))/(frqs(i+1)-frqs(i-1)))>1d-13)
     &                           Error('Detected an error in the spectral-function mesh (stretching factor inconsistent).')
          endif
        enddo
      endif

      end

c     ---------------

c     Calculate integration weights (wghtSr and wghtSi) for Hilbert trafo
c
c     suscep(j) = SUM (wghtSr(i,j)+img*wghtSi(i,j)) S(i)
c
c     where S(w) = sgn(w) Im suscep(w) is the spectral function. The weights are obtained from
c
c     INT(w1,w2) S(w') / (w-w'+i*disorder) dw'.
c
c     For the wings and head, it makes sense to replace S(w)->S'(w)=S(w)*w and S(w)->S'(w)=S(w)*w**2, respectively (multdiff1/=0).
c     Then, the weights are obtained from
c
c     INT(w1,w2) S'(w')/w'    / (w-w'+i*disorder) dw'   if md=1 (wings)
c     INT(w1,w2) S'(w')/w'**2 / (w-w'+i*disorder) dw'   if md=2 (head)
c
c     iterm = 0,1,2 : full/positive/negative frequency branch
      subroutine getwghtS(wghtSr,wghtSi,frq,nfrq,frqs,nfrqs,disorder,md,iterm)
      use global, only: img
      use, intrinsic :: iso_fortran_env
      implicit none
      integer,    intent(in)  :: nfrq,nfrqs,iterm,md
      complex_dp, intent(in)  :: frq(nfrq)
      real_dp,    intent(in)  :: frqs(nfrqs),disorder
      real_dp,    intent(out) :: wghtSr(nfrqs,nfrq),wghtSi(nfrqs,nfrq)
      real_dp                 :: winc,w1,w2,a,rdum,rdum1,rdum2
      complex_dp              :: w,cdum,int1,int2
      integer                 :: j,isign,i1,i2
      if(nfrq==0) return
      winc   = frqs(2) - frqs(1) ; a = 0.5d0 ; if(nfrqs>=3) a = winc / (frqs(3)-frqs(1))
      wghtSr = 0
      wghtSi = 0
c      Info('Hilbert mesh written to fort.111 for testing.')
c      write(111,'(F12.8)') frqs
      isign  = 1 ; if(iterm==2) isign = -1
 1    do j = 1,nfrq
        w  = frq(j) * isign
        if(frqs(1)==0) then ; i2 = 1 !; write(*,*) 'set1'
        else                ; i2 = 0 !; write(*,*) 'set0'
        endif
        ! Loop over intervals [i1,i2]
        do while(i2<nfrqs)
          i1 = i2
          ! avoid Hilbert frequencies that are too close to the pole at freq
          rdum = abs(w-frqs(i1+1))
          if(rdum>min(winc/3,1d-6)) then ; i2 = i1 + 1
          else                           ; i2 = i1 + 2 ; if(i2>nfrqs) exit
          endif
          if(i1==0) then ; w1 = max( 0d0 , frqs(i1+1)-a*winc )
          else           ; w1 = frqs(i1)
          endif
          w2   = frqs(i2)
          int1 = 0
          int2 = 0
          ! positive branch:  INT(w1,w2)   S(w') /(w-w'+iD) dw'                                     (D=disorder)
          ! negative branch: -INT(-w2,-w1) S(-w')/(w-w'-iD) dw' = INT(w1,w2) S(w')/(-w-w'+iD) dw'   (so, w -> -w)
          cdum = 0 ; rdum = 0 ; rdum1 = 0 ; rdum2 = 0
          if(w+w1-img*disorder/=0)
     &      cdum  = ln ( w-w1+img*disorder ) - ln ( w-w2+img*disorder )
          if(w1/=0.and.md>0) then
            rdum  = log( w2     ) - log( w1 )
            rdum1 =    ( w2**(-1) -      w1**(-1) )
            rdum2 =    ( w2**(-2) -      w1**(-2) ) / 2
          endif
          if(md==0) then      ! S(w') = S'(w) =  S0 + S1*w'
            int1 = cdum
            int2 = cdum * (w+img*disorder) - (w2-w1)
          else if(md==1) then ! S(w') = S'(w)/w = ( S0 + S1*w' ) / w'
            int2 = cdum
            if(w-img*disorder/=0) then ; int1 = (cdum+rdum) / (w+img*disorder)
            else                       ; int1 = rdum1
            endif
          else if(md==2) then ! S(w') = S'(w)/w**2 = ( S0 + S1*w' ) / w'**2
            if(w-img*disorder/=0) then ; int2 = (cdum+rdum)  / (w+img*disorder)
                                         int1 = (int2+rdum1) / (w+img*disorder)
            else                       ; int2 = rdum1
                                         int1 = rdum2
            endif
          endif
          if(i2-i1==1) then
            cdum         = ( w2 * int1 - int2 ) / (w2-w1)
            if(i1>0) then
            wghtSr(i1,j) = wghtSr(i1,j) + real(cdum)
            wghtSi(i1,j) = wghtSi(i1,j) + imag(cdum)
            endif
            cdum         = ( int2 - w1 * int1 ) / (w2-w1)
            wghtSr(i2,j) = wghtSr(i2,j) + real(cdum)
            wghtSi(i2,j) = wghtSi(i2,j) + imag(cdum)
          else if(i2-i1==2) then
            if(i1>0)
     &      wghtSr(i1,j) = wghtSr(i1,j) + isign
            i1           = i1 + 1
            wghtSr(i1,j) = wghtSr(i1,j) + real(int1)
            wghtSi(i1,j) = wghtSi(i1,j) + imag(int1)
            wghtSr(i2,j) = wghtSr(i1,j) - isign
          else
            Bug('i2-i1 neither 1 nor 2.')
          endif
c          if(j/=1) then ; write(*,*) frqs(i1),i1,j,wghtSr(i1,j) ; read(*,*) ; endif
        enddo
      enddo
      if(iterm==0.and.isign==1) then
        isign = -1
        goto 1
      endif

      contains

c     Principal value of complex logarithm. (The built-in function was changed in the Fortran 2003 standard, which distinguishes between +0 and -0!)
      function ln(z)
      use global, only: img,pi
      implicit none
      complex_dp             :: ln
      complex_dp, intent(in) :: z
      real_dp                :: a,b
      a = real(z)
      b = imag(z)
      if(a<0.and.b==0) then
        ln = log(-a) + img*pi
      else
        ln = log(z)
      endif
      end function ln

      end

