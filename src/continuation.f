c--------------------------------------------------------------------------------
c Copyright (c) 2021 Peter Grünberg Institut, Forschungszentrum Jülich, Germany
c This file is part of SPEX and available as free software under the conditions
c of the MIT license as expressed in the LICENSE file in more detail.
c--------------------------------------------------------------------------------

c     Fit to the model function
c
c     f(x) = SUM_i a(i) / (x-b(i))
c
c     scratch = 0: start with one pole and add one-by-one up to npole poles
c     scratch = 1: start directly with npole poles with the input a and b
c     scratch = 2: start directly with npole random poles and random weights.

# include "cppmacro.h"
      subroutine continuation(x,y,n,a,b,c,npole,scratch,allowup,fitconst,fitconstr,lwrite)
      use global, only: img MpiC(Mrank)
      use util,   only: chr
      use, intrinsic :: iso_fortran_env
      implicit none
      integer,    intent(in)    :: n,scratch
      integer,    intent(inout) :: npole
      logical,    intent(in)    :: allowup,fitconst,fitconstr,lwrite
      complex_dp, intent(in)    :: x(n),y(n)
      complex_dp, intent(inout) :: a(npole),b(npole),c
      integer                   :: np,nc,i,j,k,npole1,niter,npole0
      real_dp                   :: p(4*npole+2),error
      if(npole<=0) Error('at least one pole is required.')
      if(all(abs(y)<1d-8)) then
        npole = 1 ; a = 0 ; b = -img*1d30 ; return
      endif
      nc     = 0
      k      = 0
      npole1 = 0
      if(scratch==0) then
        npole0 = 1
        b(1)   = y(n)*x(n) / (y(n)-y(1))
        a(1)   = -b(1)*y(1)              ; if(imag(b(1))>0) b(1) = conjg(b(1))
        p      = 0
        p(1)   = real(b(1)) ; p(2) = imag(b(1))
        p(3)   = real(a(1)) ; p(4) = imag(a(1))
        c      = 0
      else if(scratch==1) then
        npole0 = npole
        j      = 0
        do i = 1,npole
          j = j + 1 ; p(j) = real(b(i))
          j = j + 1 ; p(j) = imag(b(i))
          j = j + 1 ; p(j) = real(a(i))
          j = j + 1 ; p(j) = imag(a(i))
        enddo
        if(.not.fitconst) c = 0
      else
        npole0 = npole
        call random_number(p) ; c = 0
      endif
      if(lwrite) write(6,'(A'NoA) '  Number of iterations:'
 1    do i = npole0,npole
        np = 4*i
 2      if(fitconst) then ; call fit(p,np+2,y,2*n,nc,niter,error)
        else              ; call fit(p,np  ,y,2*n,nc,niter,error)
        endif
        if(lwrite) then
          if(i/=npole0.or.k/=0) write(6,'(22X,'' '''NoA)
          write(6,'(I5,'' ('',ES7.1,'')'''NoA) niter,sqrt(error)
        endif
        if(niter==10000) then
          if(lwrite) then
            write(6,'(A)')     '  Fit not converged after 10000 iterations'
            write(6,'(38X,A)') '  Retry with random start parameters'
          endif
          call random_number(p) ; c = 0 ; k = k + 1
          if(k<=10) then ; goto 2
          else           ; exit
          endif
        endif
        k = 0
        if(any(p(2:np:4)>=0).and..not.allowup) then
c          write(*,'(4F10.5)') (p(j*4+1:j*4+4),j=0,i-1)
          if(lwrite) write(6,'(A)') '  Fit produced poles in upper half plane'
        else
          if(lwrite) write(6,'(A)') '  Fit successful'
          if(nc==2) then ; npole1 = -i ! sign here used as a flag
          else           ; npole1 =  i
          endif
          if(fitconst) c = p(np+1) + img * p(np+2)
          b(:i) = p(1:np:4) + img * p(2:np:4)
          a(:i) = p(3:np:4) + img * p(4:np:4)
        endif
        if(i/=npole) then
          if(fitconst) p(4*i+5:4*i+6) = p(4*i+1:4*i+2)
          j = maxloc(abs((p(1:np:4)+img*p(2:np:4))/(p(3:np:4)+img*p(4:np:4))),1)
          p(4*i+1:4*i+4) = p(4*j-3:4*j) * [2d0,2d0,0.5d0,0.5d0]
        endif
      enddo
      if(npole1==0) then
        if(scratch==0.or.scratch==2) then
          write(6,'(2F10.5,2F15.10)') (x(i),y(i),i=1,n)
          Error('Fit not successful.')
        else
          RWarn('Fit not successful. Return to input parameters.')
          return
        endif
      endif
      if(fitconstr.and.nc==0) then
        if(lwrite) write(6,'(A'NoA) '  Fit with constraints:'
        nc        = 2
        npole0    = abs(npole1)
        p(1:np:4) = real(b) ; p(2:np:4) = imag(b)
        p(3:np:4) = real(a) ; p(4:np:4) = imag(a)
        goto 1
      endif
      if(npole/=abs(npole1)) then
        RWarn('Number of poles reduced to: '//chr(abs(npole1)))
        npole = abs(npole1)
      endif
      if(fitconstr.and.npole1>0) then
        RWarn('Fit with constraints not successful.')
      endif
      if(lwrite) then
        write(6,'(A)') '            Position            Weight'
        do i = 1,npole
          write(6,'(A,I2,A,4F10.5)') '  Pole',i,':',b(i),a(i)
        enddo
        if(fitconst) write(6,'(A,2F10.5)') '  Constant shift:',c
      endif

      contains

      subroutine fit(p,np,y,nf,nc,niter,error)
      use wrapper
      use, intrinsic :: iso_fortran_env
      implicit none
      integer,    intent(in)    :: nf,np,nc
      integer,    intent(out)   :: niter
      real_dp,    intent(inout) :: p(np)
      real_dp,    intent(out)   :: error
      complex_dp, intent(in)    :: y(*)
      real_dp                   :: p0(np),mat(nf+nc,np),vec(nf+nc),dp(np),x1,x2,y1,y2,yy
      real_dp,    parameter     :: threshold = 1d-10
      integer                   :: i,j
      p0 = p
      i  = 0
      do
        vec   = f(p,np,y,nf,nc)
        mat   = df(p,np,nf,nc)
        error = sum(vec**2) ; if(error<threshold**2) exit
        dp    = matmul(invert(matmul(transpose(mat),mat)),matmul(transpose(mat),vec))
        if(maxval(abs(dp))>1) dp = dp/maxval(abs(dp)) ! restrict maximal step to 1
        x1 = -0.2d0 ; y1 = sum(f(p-dp*x1,np,y,nf,nc)**2)
        x2 =  1d0   ; y2 = sum(f(p-dp*x2,np,y,nf,nc)**2)
        do while(x2-x1>max(min(1d-6,error),1d-10))
          yy = sum(f(p-dp*(x1+x2)/2,np,y,nf,nc)**2)
          if(y1<y2) then
            x2 = (x1+x2)/2 ; y2 = yy
          else
            x1 = (x1+x2)/2 ; y1 = yy
          endif
        enddo
        p = p - dp * (x1+x2)/2 !; write(*,*) 'Current error:     ',error,x1
        i = i + 1
        if(maxval(abs(p-p0))<threshold) then ; j = j + 1 ; if(j==5) exit
        else                                 ; j = 0
        endif
        if(i==10000)                    then ; niter = i ; return
        endif
        p0 = p
      enddo
      niter = i
c      write(6,'(A,I5,A)') 'Minimum found after',i,' iterations.'
      end subroutine fit

      function f(p,np,y,nf,nc)
      use, intrinsic :: iso_fortran_env
      implicit none
      integer,    intent(in) :: np,nf,nc
      real_dp                :: f(nf+nc)
      real_dp,    intent(in) :: p(np)
      complex_dp, intent(in) :: y(nf/2)
      complex_dp, parameter  :: img = (0d0,1d0)
      complex_dp             :: a,b,sum
      integer                :: i,j
      do i = 1,nf/2
        if(fitconst) then ; sum = p(np-1) + img * p(np)
        else              ; sum = 0
        endif
        do j = 1,np/4*4,4
          b   = p(j)   + img * p(j+1)
          a   = p(j+2) + img * p(j+3)
          sum = sum + a * (x(i)-b)**(-1)
        enddo
        f(2*i-1) = real ( sum - y(i) )
        f(2*i)   = imag ( sum - y(i) )
      enddo
      if(nc==1.or.nc==2) then ! constraints
        if(nc==2)
     &  f(nf+2) = 0
        f(nf+1) = 0
        do j = 1,np/4*4,4
          b       = p(j)   + img * p(j+1)
          a       = p(j+2) + img * p(j+3)
          if(nc==2)
     &    f(nf+2) = f(nf+2) + imag(a/b**2) ! continuity of derivative at w=0
          f(nf+1) = f(nf+1) + imag(a)      ! integrability
        enddo
      else if(nc/=0) then
        Error('number of constraints neither 0 nor 2.')
      endif
      end function f

      function df(p,np,nf,nc)
      use, intrinsic :: iso_fortran_env
      implicit none
      integer,    intent(in) :: np,nf,nc
      real_dp                :: df(nf+nc,np)
      real_dp,    intent(in) :: p(np)
      complex_dp, parameter  :: img = (0d0,1d0)
      complex_dp             :: a,b,cdum
      integer                :: i,j
      do i = 1,nf/2
        do j = 1,np/4*4,4
          b             = p(j)   + img * p(j+1)
          a             = p(j+2) + img * p(j+3)
          cdum          = a * (x(i)-b)**(-2)
          df(2*i-1,j)   = real(cdum)     ; df(2*i,j)   = imag(cdum)
          df(2*i-1,j+1) = real(cdum*img) ; df(2*i,j+1) = imag(cdum*img)
          cdum          = (x(i)-b)**(-1)
          df(2*i-1,j+2) = real(cdum)     ; df(2*i,j+2) = imag(cdum)
          df(2*i-1,j+3) = real(cdum*img) ; df(2*i,j+3) = imag(cdum*img)
        enddo
      enddo
      if(fitconst) then
        do i = 1,nf/2
          df(2*i-1,np-1) = 1d0 ; df(2*i,np-1) = 0d0
          df(2*i-1,np)   = 0d0 ; df(2*i,np)   = 1d0
        enddo
      endif
      if(nc/=0) then
        do j = 1,np/4*4,4
          b              = p(j)   + img * p(j+1)
          a              = p(j+2) + img * p(j+3)
          df(nf+1,j:j+2) = 0d0
          df(nf+1,  j+3) = 1d0
          if(nc==2) then
            df(nf+2,j  ) = imag ( -2*a/b**3 )
            df(nf+2,j+1) = imag ( -2*a/b**3 * img )
            df(nf+2,j+2) = imag (    1/b**2 )
            df(nf+2,j+3) = imag (  img/b**2 )
          endif
        enddo
      endif
      end function df

      end

c     ----------

      subroutine write_seinfo(x,y,n,a,b,c,npole,state1,state2)
      use global, only: img,statetype
      use file
      use, intrinsic :: iso_fortran_env
      implicit none
      integer,         intent(in) :: n,npole
      complex_dp,      intent(in) :: x(n),y(n),a(npole),b(npole),c
      type(statetype), intent(in) :: state1,state2
      real_dp                     :: rdum
      integer                     :: iunit1,iunit2
      integer                     :: i,j
      iunit1 = fopen('spex.sew',status='new',numbered=.true.)
      iunit2 = fopen('spex.sec',status='new',numbered=.true.)
      if(state1%k/=state2%k) Bug('Self-energy matrix element for different kpoints.')
      if(state1%s/=state2%s) Bug('Self-energy matrix element for different spins.')
      if(state1%b==state2%b) then
        write(iunit1,'(A)') '# Expectation value of the self-energy on the imaginary frequency axis.'
        write(iunit2,'(A)') '# Expectation value of the self-energy continued to the real frequency axis.'
        write(iunit1,'(A/A,2I5)') '#','# kpt/spin',state1%k,state1%s
        write(iunit2,'(A/A,2I5)') '#','# kpt/spin',state1%k,state1%s
        write(iunit1,'(A,I5/)')       '# band    ',state1%b
        write(iunit2,'(A,I5/)')       '# band    ',state1%b        
      else
        write(iunit1,'(A)') '# Off-diagonal matrix element of the self-energy on the imaginary frequency axis.'
        write(iunit2,'(A)') '# Off-diagonal matrix element of the self-energy continued to the real frequency axis.'
        write(iunit1,'(A/A,2I5)')  '#','# kpt/spin',state1%k,state1%s
        write(iunit2,'(A/A,2I5)')  '#','# kpt/spin',state1%k,state1%s
        write(iunit1,'(A,2I5/)')       '# bands   ',state1%b,state2%b
        write(iunit2,'(A,2I5/)')       '# bands   ',state1%b,state2%b        
      endif
      do i = -n,n-1
        if(i==-1.or.i==0) cycle
        do j = 0,100
          rdum = ( (100-j) * imag(x(abs(i))) + j * imag(x(abs(i+1))) ) / 100
          if(i<0) then ; write(iunit1,'(3F20.10'NoA) -rdum, conjg ( sum(a(:npole)/(img*rdum-b(:npole))) ) + c
          else         ; write(iunit1,'(3F20.10'NoA)  rdum,         sum(a(:npole)/(img*rdum-b(:npole)))   + c
          endif
          if(j==0)   write(iunit1,'(2F20.10'NoA) real(y(abs(i)))   , sign(1,i) * imag(y(abs(i)))
          if(j==100) write(iunit1,'(2F20.10'NoA) real(y(abs(i+1))) , sign(1,i) * imag(y(abs(i+1)))
          write(iunit1,*)
        enddo
      enddo
      do i = -1000,1000
        if(i==0) cycle
        rdum = i * imag(x(n)) / 1000
        if(i<0) then ; write(iunit2,'(3F20.10)')  rdum, conjg ( sum(a(:npole)/(rdum-b(:npole))) ) + c
        else         ; write(iunit2,'(3F20.10)')  rdum,         sum(a(:npole)/(rdum-b(:npole)))   + c
        endif
      enddo
      call fclose(iunit1)
      call fclose(iunit2)
      end

