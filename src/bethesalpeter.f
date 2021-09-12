c--------------------------------------------------------------------------------
c Copyright (c) 2021 Peter Grünberg Institut, Forschungszentrum Jülich, Germany
c This file is part of SPEX and available as free software under the conditions
c of the MIT license as expressed in the LICENSE file in more detail.
c--------------------------------------------------------------------------------

c     Solves Bethe-Salpeter equation
c
c     R = K * ( 1 - W K )**(-1) = ( 1 - K W )**(-1) * K
c
c     Input:
c       K : Two-particle Green function -iGG
c       W : Screened interaction (vertical indexing according to reorder4)
c
c     Output suscepR:
c       X0        : Free susceptiblity          (contracted K) projected onto exp(ikr)/sqrt(vol)
c       X         : Renormalized susceptibility (contracted R) projected onto exp(ikr)/sqrt(vol)
c       tr(LOSS)  : trace of loss matrix [ LOSS = ImR = ( R - R^T ) / (2i) ]
c       eig(LOSS) : largest four eigenvalues of loss matrix
c
c
c     Goldstone mode at k=0 and w=0
c
c     In the absence of spin-orbit coupling, the magnetization m(r) can be rotated rigidly without a
c     cost of energy. For the rotation to be rigid, the change dm(r) of the magnetization must be
c     proportional to m(r). Let us consider the function R^(-1)=dB(r)/dm(r), giving the field dB(r)
c     that has to be applied to get the change dm(r) in m(r). Since there is no cost of energy,
c     R^(-1)*m(r)=0, so m(r) is an eigenvector of R^(-1) with eigenvalue 0. (In fact, if dm(r) is
c     not proportional to m(r), the rotation is not rigid, and a finite dB(r) is be required.)
c     Since K^(-1) is never zero (because K does not diverge around w=0), it follows that m(r) is an
c     eigenvector of K*W with eigenvalue 1. This criterion is tested by the code. The relative
c     contribution of each Wannier orbital density to m(r) as percentage is written to the output
c     and can be compared with the orbital contributions to m(r) from DFT.
c
c     Interestingly, this means, conversely, that the applied field dB(r) that would rotate the
c     magnetization rigidly, should be proportional to m(r) too. (Although the argument is a bit
c     questionable mathematically because one would invert a matrix with an eigenvalue 0, but one
c     should be able to take a non-divergent limit.) This, for example, guarantees that the applied
c     field dB(r) has the same sign as m(r) so that regions with opposite magnetization rotate in
c     the same direction. Apart from this, I frankly do not understand this point fully at the moment.
c     Anyway, it leads to the fact that m(r) is an eigenfunction of R in the limit k->0, w->0 with
c     diverging eigenvalue.
c
c     Around the Goldstone peak, the eigenvalue <m|ImR(w)|m> has the same form as <1|ImR(w)|1> because
c     <m|1> is nonzero and |m> is the dominant eigenfunction. They are related by
c     <1|ImR(w)|1> = |<m|1>|^2 <m|ImR(w)|m>. The prefactor is written.
c
c     Mathematical note:
c     K(w=0) and W(w=0) are hermitian, but neither K*W nor W*K is. The symmetrized product
c     sqrt(K)*W*sqrt(K) is hermitian, though, because K(w=0) is positive definite (W is not).
c     Therefore, this symmetrized product but also K*W and W*K have real (and the same) eigenvalues.
c     Since R^(-1) = K^(-1) * ( 1 - K*W), the product K*W should have an eigenvalue close to 1 with
c     the eigenvector m(r).
c
# include "cppmacro.h"

      subroutine bethesalpeter(suscepw,screenw,ikpt,frq,nfrq,spin,cjob)

      use global
      use wrapper
      use file
      use, intrinsic :: iso_fortran_env
      Mpi( use Mwrapper )

      implicit none
      integer,     intent(in)  :: ikpt,nfrq,spin
      character(7),intent(in)  :: cjob
      complex_dp,  intent(in)  :: frq(nfrq)
      complex_dp,  intent(in)  :: suscepw(nwan*nwan,nwan*nwan,nfrq)
      complex_dp,  intent(in)  :: screenw(nwan*nwan,nwan*nwan)
      complex_dp               :: kernel(nwan*nwan,nwan*nwan),w(nwan*nwan,nwan*nwan)
      complex_dp               :: suscepwr(nwan*nwan,nwan*nwan),orthog(nwan*nwan,nwan*nwan),invers(nwan*nwan,nwan*nwan)
      complex_dp               :: evec(nwan*nwan,nwan*nwan)
      complex_dp               :: projm(nbasm(ikpt),nwan*nwan)
      complex_dp               :: proje(nwan*nwan)
      complex_dp               :: kern(nfrq),sus(nfrq)
      real_dp                  :: eval(nwan*nwan),eig(4,nfrq),spur(nfrq)
      real                     :: cputime
      integer                  :: iunit,ikindx,k
      integer                  :: pointer(nwan*nwan)
      integer                  :: ifrq,ispin1,ispin2,nfrq1
      character(14), parameter :: dots = '..............'

      if     (spin==0) then ; ispin1 = 1 ; ispin2 = 1
      else if(spin==1) then ; ispin1 = 1 ; ispin2 = 1
      else if(spin==2) then ; ispin1 = 2 ; ispin2 = 2
      else if(spin==3) then ; ispin1 = 1 ; ispin2 = 2
      else if(spin==4) then ; ispin1 = 2 ; ispin2 = 1
      else                  ; Error('unknown spin index.')
      endif

      Rwrite(6,'(/A)') 'Solution of Bethe-Salpeter equation (BSE)'

      ! Magnetic response function will be calculated in two different ways:
      ! a) Projecting onto plane waves, which is more expensive (bare susceptibility)
      ! b) By diagonalizing the loss matrix in Wannier basis

      ! Calculate overlap of Wannier products |n1,n2*> (for the contraction of vertices)
      call overlap_wan (orthog,invers,projm,proje)

c # define testmode
# ifdef testmode
#   warning testmode!
      rewind(144)
      read(144,*) proje
# endif

      ! Test Goldstone condition and set WSCALE if requested
      w = screenw
      if(ikpt==1) then
        Rcall goldstone(wscale,suscepw,w,frq,nfrq) ; Mpi( call Mcast(wscale) ; call Mcast(w) )
      endif
      w = w * wscale

      Rwrite(6,'(A'NoA) '  Solve BSE and diagonalize loss matrix '

      call cpu_time(cputime)
      sus  = 0
      kern = 0
      eval = 0
      spur = 0
      eig  = 0
      do ifrq = Mrange(1,nfrq)
        kernel     = suscepw(:,:,ifrq)
        ! solve BSE
        suscepwr   = matmat ( kernel , invert ( identity(nwan*nwan) - matmat(w,kernel) ) )
        kern(ifrq) = dotprod ( proje , matmul ( kernel,  proje ) ) ! NAG produces internal compiler error
        sus(ifrq)  = dotprod ( proje , matmul ( suscepwr,proje ) ) ! for dot_product instead of dotprod!
        ! contract vertices and orthonormalize
        suscepwr   = matmat ( orthog , matmat ( suscepwr , orthog ) )
        ! Construct LOSS Matrix (imaginary part of suscepwr)
        suscepwr   = (suscepwr - transpose(conjg(suscepwr))) / (2*img)
        ! Take trace of LOSS Matrix
        spur(ifrq) = trace(suscepwr)
        ! Diagonalize LOSS Matrix
        call diagonalize (evec, eval, suscepwr)
        ! Reorder eigenvalues (largest first)
        call rorderp(pointer,abs(eval),nwan*nwan)
        evec        = evec(:,pointer(nwan*nwan:1:-1))
        eval        = eval(  pointer(nwan*nwan:1:-1))
        eig(:,ifrq) = eval(:4)
        ! Progress bar
        nfrq1 = nfrq Mpi(/Msize)
        Rif(14*ifrq/nfrq1>14*(ifrq-1)/nfrq1) write(6,'(A'NoA) dots(:14*ifrq/nfrq1-14*(ifrq-1)/nfrq1)
      enddo
      Rcall cpu_done(cputime)
      Mpi( call Msum(kern); call Msum(sus); call Msum(eig) ; call Msum(spur) )

      Rbegin
      if(ikpt<=nkpt) then ; ikindx = ikpt
      else                ; ikindx = 0 ; do k = 1,nkpt ; if(sum(abs(kpt(:,ikpt)-kpt(:,k)))<1d-8) ikindx = -k ; enddo
      endif
      iunit = fopen('suscepR'//cjob,status='unknown')
      write(iunit,'(A)')        '# Magnetic susceptibility X'
      write(iunit,'(A,9F10.5)') '# lattvec:',lat
      write(iunit,'(A,3F10.5)') '# k point:',kpt(:,ikpt)
      write(iunit,'(A,I6'NoA)   '# k index:',abs(ikindx) ; if(ikindx<=0) write(iunit,'(A'NoA) '  (additional)'
      write(iunit,'(/A,I2)')    '# spin:',spin
      write(iunit,'(A'NoA)      '# Freq.[meV]'; if(any(imag(frq)/=0)) write(iunit,'(''            '''NoA)
      write(iunit,'(16X,A,31X,A,24X,A,7X,A)') 'X0','X','tr(LOSS)','eig(LOSS) [(vol*eV)^-1; vol=unitcell volume]'
      do ifrq = 1,nfrq
        if(any(imag(frq)/=0)) then ; write(iunit,'(2F12.4'NoA)      frq(ifrq) *hartree*1000
        else                       ; write(iunit,'( F12.4'NoA) real(frq(ifrq))*hartree*1000
        endif
        write(iunit,'(''  '',9F16.6)') kern(ifrq)*vol/hartree,sus(ifrq)*vol/hartree,spur(ifrq)*vol/hartree,eig(:,ifrq)*vol/hartree
      enddo
      call fclose(iunit)
      write(6,'(2X,A)') 'Magnetic response function written to suscepR'//trim(cjob)//' .'
      if(ikpt==1.and.all(imag(frq)==0)) then
        ifrq     = maxloc( abs(imag(sus)),1 )
        kernel   = suscepw(:,:,ifrq)
        suscepwr = matmat ( kernel , invert ( identity(nwan*nwan) - matmat(w,kernel) ) )
        suscepwr = matmat ( orthog , matmat ( suscepwr , orthog ) )
        suscepwr = (suscepwr - transpose(conjg(suscepwr))) / (2*img)
        call diagonalize (evec, eval, suscepwr)
        write(6,'(2X,A,F10.5)') 'Projection <1|m> squared:',abs(dotprod(proje,matmul(matmul(invers,orthog),evec(:,nwan**2))))**2
      endif
      write(6,*)

      Rend

      contains

c --------------

c     orthog = overlap**(1/2)
c     invers = overlap**(-1)
c     projm  = < M(ikpt) w | w >
c     proje  = < w | w exp(iqr) > / sqrt(vol)
      subroutine overlap_wan ( orthog,invers,projm,proje )
      Mpi( use Mwrapper )
      implicit none
      complex_dp,  intent(out) :: orthog(nwan*nwan,nwan*nwan)
      complex_dp,  intent(out) :: invers(nwan*nwan,nwan*nwan)
      complex_dp,  intent(out) :: projm(nbasm(ikpt),nwan*nwan)
      complex_dp,  intent(out) :: proje(nwan*nwan)
      complex_dp,  allocatable :: cprod(:,:,:)
      MCOMPLEX_dp              :: olapm(nbasm(ikpt),nbasm(ikpt)),olap(ngptm(ikpt),ngptm(ikpt)),ctrafo(nbasm(ikpt))
      complex_dp               :: olap_wan(nwan*nwan,nwan*nwan),evec_wan(nwan*nwan,nwan*nwan)
      real_dp                  :: eval_wan(nwan*nwan),diag1(nwan*nwan),diag2(nwan*nwan)
      real_dp                  :: rdum
      real_dp,     parameter   :: wtol = 1d-4
      real                     :: cputime
      integer                  :: k,nkpack,nk,k1,k2,kk1,kk2
      integer                  :: i,n1,n2

      call olap_pw ( olap,gptm(:,pgptm(:ngptm(ikpt),ikpt)),ngptm(ikpt) )
      if(fullpw) call inverse(olap)
      olapm = blockmat ( MUNIT*identity(nbasp),olap )

      Load( call prepare_wbloch )

      ! projm(I,n1n2) = SUM(q) < M(I,k) | w(n1,q+k) w(n2,q)^*> -> < M(I,k) | w(n1) w(n2)^* >
      Rwrite(6,'(A'NoA) '  Project Wannier products onto mixed basis ..........'
      call cpu_time(cputime)
      projm = 0
      MrangeDef1(kk1,kk2,nkpt)
      rdum   = ( maxmem - mem ) Mpi(/Nsize)
      nkpack = int ( rdum / (nwan**2*nbasm(ikpt)*16d0) ) !; if(nkpack<kk2-kk1+1) call teststop('nkpack<kk2-kk1+1')
      do k1 = kk1,kk2,nkpack
        k2 = min(k1+nkpack-1,kk2)
        nk = k2 - k1 + 1
        allocate ( cprod(nwan*nwan,nk,nbasm(ikpt)) )
        call wavefproducts3_mt(cprod,              [(k,k=k1,k2)],nk,ikpt,ispin2,ispin1,1,nwan,1,nwan)
        call wavefproducts3_pw(cprod(:,:,nbasp+1:),[(k,k=k1,k2)],nk,ikpt,ispin2,ispin1,1,nwan,1,nwan)
        do i = 1,nwan*nwan
          do k = 1,nk
            projm(:,i) = projm(:,i) + cprod(i,k,:)
          enddo
        enddo
        deallocate ( cprod )
      enddo
      projm = projm / nkpt
      Mpi( call Msum(projm) )
      Rcall cpu_done(cputime)

      Load( Deallocate_(cmtu) )
      Load( Deallocate_(cpwu) )

      Rbegin

      ! Symmetrize in case of inversion symmetry
      Inv( call symmetrize(projm,nbasm(ikpt),nwan*nwan,1) )

      ! Calculate the overlap of Wannier functions <n1n2|n3n4>
      write(6,'(A'NoA) '  Calculate and invert overlap of Wannier products ...'
      do n2 = 1,nwan*nwan
        do n1 = 1,nwan*nwan
          olap_wan(n1,n2) = dotprod ( projm(:,n1), matmul(olapm,projm(:,n2)) )
        enddo
      enddo

      if(.not.fullpw) projm = matmul(olapm,projm) ! return <M|...> instead of <M~|...> for APPROXPW

      call getctrafo(ctrafo MpiC(MPI_UNDEFINED) ,1,[1],ikpt)
      proje = conjg ( matmul(MCONJG(ctrafo),projm) )

      call diagonalize ( evec_wan,eval_wan,olap_wan )
      if(any(eval_wan<-1d-10)) Warn('Overlap eigenvalue(s) negative.')
      diag1 = 0
      diag2 = 0
      where(eval_wan>wtol)
        diag1 = sqrt(eval_wan)
        diag2 = 1/eval_wan
      end where
      orthog = matmat(evec_wan,diagmat(diag1,conjg(transpose(evec_wan))))
      invers = matmat(evec_wan,diagmat(diag2,conjg(transpose(evec_wan))))
      call cpu_done(cputime)
      if(nbasm(ikpt)<count(diag1/=0)) Warn('Mixed basis smaller than Wannier product basis.')
      write(6,'(A,I8)')   '  Number of Wannier products:   ',count(diag2/=0)
      write(6,'(A,F8.5)') '  Ability to represent exp(iqr):',
     &  real ( dotprod ( proje , matmul ( invers , proje ) ) )

      Rend

# ifdef MPI
      call Mcast(proje)
      call Mcast(orthog)
      call Mcast(invers)
# endif

      end subroutine overlap_wan

c --------------

      end

c --------------

      ! Reorder from 13 to 13  (r)
      !              42    24  (r')
      subroutine reorder_matrix(matout,matin)
      use global, only: nwan
      use, intrinsic :: iso_fortran_env
      implicit none
      complex_dp, intent(out)  :: matout(nwan,nwan,nwan,nwan)
      complex_dp, intent(in)   :: matin(nwan,nwan,nwan,nwan)
      integer                  :: n1,n3
      do n1 = 1,nwan
        do n3 = 1,nwan
          matout(n1,:,n3,:) = transpose(matin(n1,:,n3,:))
        enddo
      enddo
      end

      subroutine reorder4(mat)
      use global, only: nwan
      use, intrinsic :: iso_fortran_env
      implicit none
      complex_dp, intent(inout) :: mat(nwan,nwan,nwan,nwan)
      integer                   :: n1,n3
      do n1 = 1,nwan
        do n3 = 1,nwan
          mat(n1,:,n3,:) = transpose(mat(n1,:,n3,:))
        enddo
      enddo
      end

c --------------

      ! Tests Goldstone condition and
      ! - sets wscale if requested (WSCALE 0)
      ! - modifies w by cutting Goldstone mode if requested (CUTGOLD)

      subroutine goldstone(wscale,suscepw,w,frq,nfrq)

      use global,  only: nwan,inp
      use wrapper, only: diagonalize,diagonalize_gen,matmat,invert,dotprod
      use key
      use, intrinsic :: iso_fortran_env

      implicit none
      integer,    intent(in)    :: nfrq
      real_dp,    intent(inout) :: wscale
      complex_dp, intent(in)    :: suscepw(nwan**2,nwan**2,nfrq),frq(nfrq)
      complex_dp, intent(inout) :: w(nwan**2,nwan**2)
      complex_dp                :: evec(nwan**2,nwan**2),carr(nwan**2)
      real_dp                   :: eval(nwan**2),error
      logical                   :: cutgold,ldum
      integer                   :: ifrq,ieig,i,j,k

      ifrq = minloc(abs(frq),1)
      if(abs(frq(ifrq))>1d-10) then
        if(wscale==0) Error('Missing zero frequency.')
        if(abs(real(frq(ifrq)))>1d-10) then
          write(6,'(A)') 'Cannot test Goldstone condition as Re(w)!=0.'
          return
        endif
        write(6,'(2X,A,F10.5,F13.10)') 'Test Goldstone condition (1-KW)*m=0 at Im(w)=',imag(frq(ifrq))
      else
        write(6,'(2X,A,F10.5)')        'Test Goldstone condition (1-KW)*m=0'
      endif

      call diagonalize_gen(evec,carr,matmat(suscepw(:,:,ifrq),w))
      call diagonalize(eval,suscepw(:,:,ifrq))
      if(all(eval>=0).or.all(eval<=0)) then
        if     (all(eval>=0)) then ; write(6,'(2X,A)') 'K(w=0) is positive definite.'
        else if(all(eval<=0)) then ; write(6,'(2X,A)') 'K(w=0) is negative definite.'
        else                       ; write(6,'(2X,A)') 'K(w=0) is indefinite.'
        endif
        if(any(abs(imag(carr))>1d-10)) then
          call getkey(inp,'DISORDER',ldum,section='SUSCEP',default=.false.)
          if(ldum) then ; write(6,'(A)') 'K*W has complex eigenvalues because of DISORDER (?) Goldstone test skipped.' ; return
          else          ; Error('K*W should have real eigenvalues but has not.')
          endif
        endif
        ieig = minloc(abs(1-carr),1)
        write(6,'(2X,A,F10.5)') 'Eigenvalue of K*W closest to 1:',real(carr(ieig))
      else
        write(6,'(2X,A)') 'K(w=0) is indefinite.'
        ieig = minloc(abs(1-carr),1)
        write(6,'(2X,A,2F10.5)') 'Eigenvalue of K*W closest to 1:',carr(ieig)
      endif
      if(wscale==0) then
        wscale = 1/real(carr(ieig))
        write(6,'(2X,A,F11.5)') 'Scaling factor for W (1/eigv):',wscale
      endif
      call getkey(inp,'CUTGOLD',cutgold,section='WANNIER',default=.false.)
      if(cutgold) then
        carr = carr(ieig) * matmul(invert(suscepw(:,:,ifrq)),evec(:,ieig))
        do j = 1,nwan**2 ; do i = 1,nwan**2
          w(i,j) = w(i,j) - carr(i) * conjg(evec(j,ieig))
        enddo ; enddo
      endif
      carr  = evec(:,ieig)
      carr  = carr / sum(carr)
      error = 0
      write(6,'(2X,A,F9.5,A)') 'Goldstone eigenfunction m(r) at w=0'
      do i = 1,nwan
        do j = 1,nwan
          k = i+(j-1)*nwan
          if(i==j) then
            write(6,'(2X,I3,F10.5,A'NoA) i,real(carr(k))*100,' %'
            if(abs(imag(carr(k)))>1d-12) write(6,'(A,F10.5'NoA) '   imag:',imag(carr(k))*100
            write(6,*)
          else
            error = error + abs(carr(k))
          endif
        enddo
      enddo
      write(6,'(2X,A,F10.5)') 'Sum of off-diagonal elements:',error

      ! CUTGOLD:
      ! Expectation value of K(w)*W'with |c>=evec(:,ieig) must be zero for all frequencies w  ( W' = W - <c|K(0)W|c> K(0)^(-1)|c><c| )
      ! <c|K(w)W'|c> = <c|K(w)W|c> - <c|K(w)K(0)^(-1)|c> <c|K(0)W|c> = 0  ( because <c|K(0)W|c> is diagonal )
      if(cutgold) then
        carr = matmul(w,evec(:,ieig))
        do ifrq = 1,nfrq
          if(abs(dotprod(evec(:,ieig),matmul(suscepw(:,:,ifrq),carr)))>1d-10) Error('CUTGOLD test failed.')
        enddo
      endif

      end

c --------------

# define LEADINGORDER 3
# define RENORMALIZE
# if LEADINGORDER != 3
#   warning LEADINGORDER is not 3
# endif
      subroutine bethesalpeter_t(matrix,screenw,frq,nfrq,ikpt)
      use global, only: nwan,wscale MpiC(Mrank) MpiC(Msize)
      use util
      use wrapper
      use, intrinsic :: iso_fortran_env
      Mpi( use Mwrapper )
      implicit none
      integer,    intent(in)    :: nfrq,ikpt
      complex_dp, intent(in)    :: frq(nfrq)
      complex_dp, intent(in)    :: screenw(nwan*nwan,nwan*nwan)
      complex_dp, intent(inout) :: matrix(nwan*nwan,nwan*nwan,nfrq)
      complex_dp                :: w(nwan*nwan,nwan*nwan),t(nwan*nwan,nwan*nwan),kw(nwan*nwan,nwan*nwan)
      integer                   :: ifrq,i
      real                      :: cputime

      w = screenw
      if(ikpt==1) then
        Rcall goldstone(wscale,matrix,w,frq,nfrq) ; Mpi( call Mcast(wscale) ; call Mcast(w) )
      endif
      w = w * wscale

      Rwrite(6,'(/A'NoA) 'Calculation of T matrix... '
      call cpu_time(cputime)

      do ifrq = Mrange1(nfrq)
        kw = matmat(matrix(:,:,ifrq),w)
        t  = w
        do i = 1,LEADINGORDER-1
          t = matmat(t,kw)
        enddo
# ifdef RENORMALIZE
        t = matmat(t,invert(identity(nwan*nwan) - kw ) )
# else
#   warning no renormalization in bethesalpeter
# endif
        matrix(:,:,ifrq) = t
      enddo

      MrangeDistr(matrix(:,:,McolD1(nfrq,i)),i)

      Rcall cpu_done(cputime)

c      write(*,*) 'updated T:',sum(matrix),sum(abs(matrix))
c      read(*,*)

      end

c --------------


