c--------------------------------------------------------------------------------
c Copyright (c) 2021 Peter Grünberg Institut, Forschungszentrum Jülich, Germany
c This file is part of SPEX and available as free software under the conditions
c of the MIT license as expressed in the LICENSE file in more detail.
c--------------------------------------------------------------------------------

c Common subroutines for selfenergy and selfenergy_wannier

# include "cppmacro.h"
# include "jobtype.h"
# include "restype.h"

c -----------------------------

c
c     Add contribution from W contraction to self-energy (called after k loop)
      subroutine selfenergy_contract(eval,job1)
      use global
      use arrays
      implicit none
      type(jobtype), intent(in) :: job1
      logical,       intent(in) :: eval(2)
      real                      :: cputime
      Mpi( integer              :: Merr )
      Rwrite(6,*)
c     (1) IBC: Contraction contribution
      if(eval(2).and.oibc/=0) then
        call ibc2_contraction(job1,freq*img,nfreq)
        ! integrate missing ...
      endif
c     (2) CORES: core contribution
      if(eval(2).and.any(cores)) then
        Rwrite(6,'(A'NoA) 'Core contribution to self-energy... '
        Nfence(selfc)
        Rcall cpu_time(cputime)
        Rcall selfenergy_core(job1)
        Rcall cpu_done(cputime)
        Nfence(selfc)
      endif
c     (3) COSX: Coulomb-hole contribution
      if(eval(1).and.job1%type==J_COSX) then
        Rwrite(6,'(A'NoA) 'Coulomb-hole contribution to self-energy... '
        Rcall cpu_time(cputime)
        call coulombhole(job1)
        Rcall cpu_done(cputime)
      endif
      end

c ------------------------------

c
c Calculates the core contribution to the self-energy.
c
c Implemented in an analogous (but much simpler) way to above.
c
c The projections (cprod) calculated below are
c              s        c,s
c < M       phi    | phi       > = delta    exp(-ikR ) / sqrt(N)
c    k,aLMP    n,q      n',q+k          aa'         a
c
c                                    L+l            s *                         R_MT  s         c,s
c                                *   SUM    SUM  cmt              C              INT u     (r) u   (r) M   (r) dr
c                                  l'=|L-l|  p'     nq,al'(m-M)p'  l'(m-M),lm,LM  0   al'p'     alp     aLP
c
c with the core state n'=(a'lmp).
c The array screenk_mt defined in the parent routine is used.
c
c CORESOC: SOC split states are averaged and treated non-relativistically
c
      subroutine selfenergy_core(job1)
      use wrapper
      use global
      use arrays
      use freq_integral
      implicit none
      type(jobtype), intent(in) :: job1
      complex_dp, allocatable   :: cprod(:,:),mat(:,:,:)
      complex_dp                :: cprod1(maxlmindxm)
      complex_dp                :: cvec(size(selfc,2)),cpade(nfreq+1),pole(nfreq*smooth(2)),resid(nfreq*smooth(2)),cdum
      complex_dp                :: wfreqintegral(0:3,nfreq,size(selfc,2))
      complex_dp                :: aspline(0:3,nfreq)
      real_dp                   :: freqr1(nfreqr)
      real_dp                   :: integral(maxindx,0:maxlcut,maxindxm,0:maxlcutm),gnt,ecor,enediff,rdum
      integer                   :: i,j
      integer                   :: itype,ieq,ic,lc,nc,mc,ll,nn,mm,l,n,m,llm,lm,s
      integer                   :: ifreq,ifreqr,iself,iselfc,npole1
      integer                   :: iblock,ibandq,ikptq,ispin,nb
      logical                   :: ldum
      real_dp                   :: intgrf,gaunt

      if(l_soc) call teststop('CORES & GW & SOC')
      Load( Error('Not implemented for LOAD yet.') )

      iselfc = 0
      do iblock = 1,nblock
        i      = block(1,iblock) ! first band of block
        ibandq = job1%band(i)
        ikptq  = job1%kpt(i) ; if(storeibz) ikptq = kptp(ikptq)
        ispin  = job1%spin(i)
        nb     = sizeblock(iblock)
        allocate ( cprod(maxlmindxm,nb),mat(nb,nb,nfreq) )

        ! Define frequency arguments w of self-energy SIGMA(w)  (-> freqr1)
        if(oselfc==1.or.oselfc==4) then
          freqr1 = freqr
        else
          if(ene(ibandq,ikptq,ispin)<=efermi) then ; freqr1 = ene(ibandq,ikptq,ispin) - efermi - freqr
          else                                     ; freqr1 = ene(ibandq,ikptq,ispin) - efermi + freqr
          endif
        endif

        ! Loop over core states
        do itype = 1,ntype ; if(.not.any(cores(:,:,itype))) cycle
          do lc = 0,lcutc(itype)
            do nc = 1,nindxc(lc,itype)
              ldum = cores(nc,lc,itype)
              if(lcore_soc.and.lc/=0) then
                if(mod(nc,2)==0) then ; cycle
                else                  ; ldum = ldum.or.cores(nc+1,lc,itype)
                endif
              endif

              if(ldum) then

                ! Calculate radial integrals (-> integral)
                do ll = 0,lcutm(itype)
                  do nn = 1,nindxm(ll,itype)
                    do l = abs(ll-lc),min(ll+lc,lcut(itype)),2
                      do n = 1,nindx(l,itype)
                        if(lcore_soc.and.lc/=0) then
                          integral(n,l,nn,ll) = intgrf ( basm(:,nn,ll,itype) / rgrid(:,itype) / (2*lc+1) *
     &                                       ( matmul(core1(:,nc:nc+1,lc,itype,ispin),[lc+1,lc]) * bas1(:,n,l,itype,ispin) +
     &                                         matmul(core2(:,nc:nc+1,lc,itype,ispin),[lc+1,lc]) * bas2(:,n,l,itype,ispin) ) ,itype)
                        else
                          integral(n,l,nn,ll) = intgrf ( basm(:,nn,ll,itype) / rgrid(:,itype) *
     &                                                 ( core1(:,nc,lc,itype,ispin) * bas1(:,n,l,itype,ispin) +
     &                                                   core2(:,nc,lc,itype,ispin) * bas2(:,n,l,itype,ispin) ) ,itype)
                        endif
                      enddo
                    enddo
                  enddo
                enddo

                ! Loop over equivalent atoms and mc (ecor identical!) and sum cprod*screen*cprod (-> mat)
                mat  = 0
                ic   = sum(neq(:itype-1))
                ecor = ecore(nc,lc,itype,ispin) - efermi
                if(lcore_soc.and.lc/=0) ecor = ( (lc+1)*ecor + lc*(ecore(nc+1,lc,itype,ispin)-efermi) ) / (2*lc+1)

                do ieq = 1,neq(itype)
                  ic = ic + 1
                  do mc = -lc,lc

                    ! Calculate projection < M nq | core > (-> cprod)
                    cprod = 0
                    do s = ispin,ispin+nspin3-1 ! SOC loop
                    llm = 0
                    do ll = 0,lcutm(itype)
                      do mm = -ll,ll
                        m = mc - mm
                        do nn = 1,nindxm(ll,itype)
                          llm = llm + 1
                          do l = max(abs(ll-lc),abs(m)),min(ll+lc,lcut(itype))
                            gnt = gaunt(l,lc,ll,m,mc,mm) ; if(gnt==0) cycle
                            lm  = sum ( [ ((2*i+1)*nindx(i,itype),i=0,l-1) ] ) + (l+m) * nindx(l,itype)
                            do n = 1,nindx(l,itype)
                              lm   = lm + 1
                              rdum = integral(n,l,nn,ll) * gnt ; if(abs(rdum)<1d-10) Bug('rdum=0')
                              do i = 1,nb
                                cprod(llm,i) = cprod(llm,i) + rdum * conjg(cmt(lm,ic,job1%band(block(i,iblock)),ikptq,s))
                              enddo
                            enddo
                          enddo
                        enddo
                      enddo
                    enddo
                    enddo ! SOC loop
                    if(llm/=sum([((2*i+1)*nindxm(i,itype),i=0,lcutm(itype))])) Bug('Count error.')

                    ! Accumulate cprod(i)*screen*cprod(j) (-> mat)
                    do ifreq = 1,nfreq
                      do i = 1,nb
                        cprod1(:llm) = matvec( screenk_mt(:llm,:llm,ic,ifreq) , cprod(:llm,i) )
                        do j = 1,i
                          mat(i,j,ifreq) = mat(i,j,ifreq) + dotprod( cprod(:llm,j) , cprod1(:llm) )
                          mat(j,i,ifreq) = conjg(mat(i,j,ifreq))
                        enddo
                      enddo
                    enddo

                  enddo ! mc
                enddo ! ieq

                ! Frequency convolutions (-> selfc)
                if(freqint<=1) then
                  if(oselfc==1) then ; cvec = freq + img * ecor
                  else               ; cvec = img * (ecor-freqr1)
                  endif
                  call freqintegral_init(wfreqintegral,freq,nfreq,cvec,size(cvec))
                endif
                do j = 1,nb
                  do i = 1,j
                    if(maxval(abs(mat(i,j,:)))>1d-10) then
                      if(freqint==2) then
                        call pade_init(cpade,img*freq,mat(i,j,:),nfreq,-1)
                        call pade_poles(pole,resid,npole1,img*freq,mat(i,j,:),cpade,nfreq,smooth(2),-1,.false.)
                      else
                        call spline_init(aspline,freq,mat(i,j,:),nfreq)
                      endif
                      do ifreq = 1,size(selfc,2)
                        if(oselfc==1) then ; ifreqr = 1     ; rdum = freq(ifreq)
                        else               ; ifreqr = ifreq ; rdum = 0
                        endif
                        enediff = ecor - efermi - freqr1(ifreqr)
                        if(freqint<=1) then ; cdum = sum(wfreqintegral(:,:,ifreq)*aspline) / (2*pi*img)
                        else                ; cdum = freqintegral_poles(pole,resid,npole1,1,rdum+img*enediff,0)
                        endif
                        iself              = iselfc + i + (j-1)*nb
                        selfc(iself,ifreq) = selfc(iself,ifreq) - cdum / nkpt
                        if(i/=j) then
# ifndef INV
                          if(freqint<=1) then ; cdum = sum(wfreqintegral(:,:,ifreq)*conjg(aspline)) / (2*pi*img)
                          else                ; cdum = freqintegral_poles(-conjg(pole),-conjg(resid),npole1,1,rdum+img*enediff,0)
                          endif
# endif
                          iself              = iselfc + j + (i-1)*nb
                          selfc(iself,ifreq) = selfc(iself,ifreq) - cdum / nkpt
                        endif
                      enddo
                    endif
                  enddo
                enddo

              endif

            enddo
          enddo
        enddo

        deallocate ( cprod,mat )
        iselfc = iselfc + nb**2

      enddo ! blocks

      end 

c ------------------------------

      subroutine coulombhole(job1)
      use global
      use arrays
      Mpi( use Mwrapper )
      Load( use readwrite )
      implicit none
      type(jobtype), intent(in) :: job1
      MCOMPLEX_dp,  allocatable :: screenk(:,:,:),stepfunc(:,:,:)
      complex_dp                :: selfx1(size(selfx))
      MCOMPLEX_dp,  allocatable :: cpw1(:,:,:),cpw2(:)
      complex_dp,   allocatable :: cmt1(:,:,:)
      complex_dp                :: cmt2(maxlmindx),coh(maxgrid,(2*maxlcutm+1)**2)
      complex_dp                :: integral(maxindx,maxindx,(2*maxlcutm+1)**2,0:maxlcut,0:maxlcut,nspin2)
      real_dp                   :: basbas(maxgrid,maxindx,maxindx),gnt
      integer                   :: itype,ieq,ic,ng
      integer                   :: l,l1,l2,m,m1,m2,n1,n2,lm,lm1,lm2
      integer                   :: iblock,ikptq,ispin,i1,i2,ib1,ib2,s,nb
      integer                   :: ig1,ig2,ig3
      integer                   :: jg1,jg2,jg3
      integer                   :: g(3),g1(3),g2(3)
      integer                   :: iselfx,iselfx0
      MCOMPLEX_dp               :: stepfunction
      real_dp                   :: gaunt
      complex_dp                :: cintgrf
# ifdef MPI
      logical :: storeibz1
      storeibz1 = storeibz
      storeibz  = .false.
# endif

c      write(*,*) sum(screenk_mt),sum(screenk_pw)

c
c     MT
c
      selfx1 = 0
      ic     = 0
      do itype = 1,ntype
        ng = grid(itype)%number
        do ieq = 1,neq(itype)
          ic  = ic + 1

          coh = 0

          lm1 = 0
          do l1 = 0,lcutm(itype)
            do m1 = -l1,l1
              lm2 = 0
              do l2 = 0,lcutm(itype)
                do m2 = -l2,l2
                  do l = abs(l1-l2),l1+l2,2
                    m   = m1 - m2 ; if(abs(m)>l) cycle
                    lm  = l*(l+1) + m + 1 ; Mcycle(lm)
                    gnt = gaunt(l,l1,l2,m,m1,m2)

                    if(gnt/=0) then
                      do n1 = 1,nindxm(l1,itype)
                        do n2 = 1,nindxm(l2,itype)
                          coh(:ng,lm) = coh(:ng,lm) + gnt * screenk_mt(lm1+n1,lm2+n2,ic,1) *
     &                                  basm(:ng,n1,l1,itype) * basm(:ng,n2,l2,itype) / rgrid(:ng,itype)**2
                        enddo
                      enddo
                    endif

                  enddo
                  lm2 = lm2 + nindxm(l2,itype)
                enddo
              enddo
              lm1 = lm1 + nindxm(l1,itype)
            enddo
          enddo

          integral = 0

          do ispin = 1,nspin2
            do l1 = 0,lcut(itype)
              do l2 = 0,lcut(itype)

                do n2 = 1,nindx(l2,itype)
                  do n1 = 1,nindx(l1,itype)
                    basbas(:ng,n1,n2) = bas1(:ng,n1,l1,itype,ispin)*bas1(:ng,n2,l2,itype,ispin) +
     &                                  bas2(:ng,n1,l1,itype,ispin)*bas2(:ng,n2,l2,itype,ispin)
                  enddo
                enddo

                do l = abs(l1-l2),min(l1+l2,2*lcutm(itype)),2
                  lm = l**2
                  do m = -l,l
                    lm = lm + 1 ; Mcycle(lm)

                    do n2 = 1,nindx(l2,itype)
                      do n1 = 1,nindx(l1,itype)
                        integral(n1,n2,lm,l1,l2,ispin) = cintgrf( coh(:ng,lm) * basbas(:ng,n1,n2) , itype )
                      enddo
                    enddo

                  enddo
                enddo

              enddo
            enddo
          enddo

          iselfx0 = 0
          do iblock = 1,nblock
            ib1    = block(1,iblock) ! first band of block
            ikptq  = job1%kpt(ib1)
            ispin  = job1%spin(ib1)
            nb     = sizeblock(iblock)

            if(l_soc) then ; allocate ( cmt1(maxlmindx,nb,2)           )
            else           ; allocate ( cmt1(maxlmindx,nb,ispin:ispin) )
            endif
# ifdef LOAD
            allocate ( cmt(maxlmindx,ncent,nb,1,nspin3) )
            call read_wavef2(job1%band(block(:nb,iblock)),nb,ikptq,ispin,cmt)
            kindx(ikptq) = 1
# endif

            do i2 = 1,nb
              ib2 = ifLoad( i2 , job1%band(block(i2,iblock)) )

# ifndef old_trafo
              call wavefunction_mt(cmt1(:,i2,:),maxlmindx,ic,ib2,ikptq,ispin)
# else
              if(storeibz.and.kptp(ikptq)/=ikptq) then
                if(l_soc) then
                  call waveftrafo_mt(cmt1(:,i2,1),maxlmindx,ikptq,ib2,1,ic,.false.)
                  call waveftrafo_mt(cmt1(:,i2,2),maxlmindx,ikptq,ib2,2,ic,.false.)
                  call waveftrafo_soc(cmt1(:,i2,:),maxlmindx,symkpt(ikptq))
                else
                  call waveftrafo_mt(cmt1(:,i2,ispin),maxlmindx,ikptq,ib2,ispin,ic,.false.)
                endif
              else
                if(l_soc) then ; cmt1(:,i2,:)     = cmt(:,ic,ib2,kindx(ikptq),:)
                else           ; cmt1(:,i2,ispin) = cmt(:,ic,ib2,kindx(ikptq),ispin)
                endif
              endif
# endif
              do s = 1,2 ; if(.not.l_soc.and.s/=ispin) cycle ! SOC spin

              cmt2 = 0
              lm1  = 0
              do l1 = 0,lcut(itype)
                do m1 = -l1,l1
                  lm2 = 0
                  do l2 = 0,lcut(itype)
                    do m2 = -l2,l2

                      do l = abs(l1-l2),min(l1+l2,2*lcutm(itype)),2
                        m   = m1 - m2 ; if(abs(m)>l) cycle
                        lm  = l*(l+1) + m + 1 ; Mcycle(lm)
                        gnt = gaunt(l,l1,l2,m,m1,m2)

                        do n2 = 1,nindx(l2,itype)
                          do n1 = 1,nindx(l1,itype)
                            cmt2(lm1+n1) = cmt2(lm1+n1) + gnt * cmt1(lm2+n2,i2,s) * integral(n1,n2,lm,l1,l2,s)
                          enddo
                        enddo

                      enddo

                      lm2 = lm2 + nindx(l2,itype)
                    enddo
                  enddo
                  lm1 = lm1 + nindx(l1,itype)
                enddo
              enddo

              iselfx = iselfx0

              do i1 = 1,i2
                iselfx         = iselfx + 1
                selfx1(iselfx) = selfx1(iselfx) + dot_product(cmt1(:,i1,s),cmt2) / (2*nkpt)
              enddo

              enddo ! SOC spin

              iselfx0 = iselfx

            enddo
            deallocate ( cmt1 )
            Load( deallocate ( cmt ) )

          enddo

        enddo
      enddo

c
c     PW
      ! Dimension screenk (G vectors needed for multiplication with phi from both sides)
      g1 =  huge(0)
      g2 = -huge(0)
      do ikptq = 1,nkpt
        do ig1 = 1,ngpt(ikptq)
          do ig2 = 1,ngpt(ikptq)
            g  = gpt(:,pgpt(ig2,ikptq)) - gpt(:,pgpt(ig1,ikptq))
            g1 = min(g1,g)
            g2 = max(g2,g)
          enddo
        enddo
      enddo
      allocate ( screenk(g1(1):g2(1),g1(2):g2(2),g1(3):g2(3)) )
      screenk = 0
      ! Dimension step function
      g1 =  huge(0)
      g2 = -huge(0)
      do ig1 = lbound(screenk,1),ubound(screenk,1)
      do ig2 = lbound(screenk,2),ubound(screenk,2)
      do ig3 = lbound(screenk,3),ubound(screenk,3)
        do jg1 = lbound(screenk_pw,1),ubound(screenk_pw,1)
        do jg2 = lbound(screenk_pw,2),ubound(screenk_pw,2)
        do jg3 = lbound(screenk_pw,3),ubound(screenk_pw,3)
          g  = [ ig1,ig2,ig3 ] - [ jg1,jg2,jg3 ]
          g1 = min(g1,g)
          g2 = max(g2,g)
        enddo
        enddo
        enddo
      enddo
      enddo
      enddo
      allocate ( stepfunc(g1(1):g2(1),g1(2):g2(2),g1(3):g2(3)) )
      ! Define step function
      do ig1 = g1(1),g2(1)
      do ig2 = g1(2),g2(2)
      do ig3 = g1(3),g2(3)
        g                     = [ig1,ig2,ig3]
        stepfunc(ig1,ig2,ig3) = stepfunction(g)
      enddo
      enddo
      enddo
      ! Multiply step function with screenk_pw (->screenk)
      do ig1 = lbound(screenk,1),ubound(screenk,1)
      do ig2 = lbound(screenk,2),ubound(screenk,2)
      do ig3 = lbound(screenk,3),ubound(screenk,3)
        do jg1 = lbound(screenk_pw,1),ubound(screenk_pw,1)
        do jg2 = lbound(screenk_pw,2),ubound(screenk_pw,2)
        do jg3 = lbound(screenk_pw,3),ubound(screenk_pw,3)
          g                    = [ ig1,ig2,ig3 ] - [ jg1,jg2,jg3 ]
          screenk(ig1,ig2,ig3) = screenk(ig1,ig2,ig3) + stepfunc(g(1),g(2),g(3)) * screenk_pw(jg1,jg2,jg3)
        enddo
        enddo
        enddo
      enddo
      enddo
      enddo
      ! Multiply with wave functions
      iselfx0 = 0
      do iblock = 1,nblock
        ib1    = block(1,iblock) ! first band of block
        ikptq  = job1%kpt(ib1)
        ispin  = job1%spin(ib1)
        nb     = sizeblock(iblock)

        if(l_soc) then ; allocate ( cpw2(ngpt(ikptq)),cpw1(ngpt(ikptq),nb,2)           )
        else           ; allocate ( cpw2(ngpt(ikptq)),cpw1(ngpt(ikptq),nb,ispin:ispin) )
        endif

# ifdef LOAD
        allocate ( cpw(maxgpt,nb,1,nspin2) )
        call read_wavef2(job1%band(block(:nb,iblock)),nb,ikptq,ispin,cpwout=cpw)
        kindx(ikptq) = 1
# endif

        do i2 = 1,nb
          ib2 = ifLoad( i2 , job1%band(block(i2,iblock)) )
# ifndef old_trafo
          call wavefunction_pw(cpw1(:,i2,:),ngpt(ikptq),ib2,ikptq,ispin)
# else
          if(storeibz.and.kptp(ikptq)/=ikptq) then
            if(l_soc) then
              call waveftrafo_pw(cpw1(:,i2,1),ikptq,ib2,1,.false.)
              call waveftrafo_pw(cpw1(:,i2,2),ikptq,ib2,2,.false.)
              call waveftrafo_soc(cpw1(:,i2,:),ngpt(ikptq),symkpt(ikptq))
            else
              call waveftrafo_pw(cpw1(:,i2,ispin),ikptq,ib2,ispin,.false.)
            endif
          else
            if(l_soc) then ; cpw1(:,i2,:)     = cpw(:ngpt(ikptq),ib2,kindx(ikptq),:)
            else           ; cpw1(:,i2,ispin) = cpw(:ngpt(ikptq),ib2,kindx(ikptq),ispin)
            endif
          endif
# endif

          do s = 1,2 ; if(.not.l_soc.and.s/=ispin) cycle ! SOC spin

          cpw2 = 0
          do ig1 = 1,ngpt(ikptq) ; Mcycle(ig1)
            do ig2 = 1,ngpt(ikptq)
              g         = gpt(:,pgpt(ig1,ikptq)) - gpt(:,pgpt(ig2,ikptq))
              cpw2(ig1) = cpw2(ig1) + cpw1(ig2,i2,s) * screenk(g(1),g(2),g(3))
            enddo
          enddo

          iselfx = iselfx0

          do i1 = 1,i2
            iselfx         = iselfx + 1
            selfx1(iselfx) = selfx1(iselfx) + dot_product(cpw1(:,i1,s),cpw2) / (2*nkpt*vol)
          enddo

          enddo ! SOC spin

          iselfx0 = iselfx

        enddo
        deallocate ( cpw1,cpw2 )
        Load( deallocate ( cpw ) )

      enddo
      Mpi( call Msum(selfx1,0) )
# ifdef INV
      Rif(any(abs(imag(selfx1))>1d-10)) Error('Nonzero imaginary element in coulombhole term.')
# endif
      ifR selfx = selfx + selfx1

      deallocate ( screenk,stepfunc )

      Mpi( storeibz = storeibz1 )

      end 

c ------------------------------

c Returns coefficients for interpolation (linear interpolation at the moment). Assumes arr(:) to be size-ordered.
      subroutine getcoeff(coeff,x,arr,n)
      use util, only: chr
      use, intrinsic :: iso_fortran_env
      implicit none
      integer, intent(in)  :: n
      real_dp, intent(out) :: coeff(n)
      real_dp, intent(in)  :: x,arr(n)
      real_dp, save        :: maxdev = 0.05d0
      integer              :: i,i1,i2
      coeff = 0
      if(x<=arr(1).or.x>=arr(n)) then
        if(x<=arr(1)) then ; i = 1
        else               ; i = n
        endif
        coeff(i) = 1d0
        if(abs(x-arr(i))>maxdev) then
          Warn('Frequency argument of W exceeds frequency mesh by '//trim(chr(abs(x-arr(i))))//)
     &         '. You might want to remove spex.cor/t or rerun without RESTART.'
          maxdev = abs(x-arr(i))
        endif
      else
        i1 = 1
        i2 = n
        do while(i1+1<i2)
          i = (i1+i2)/2
          if(x<arr(i)) then ; i2 = i
          else              ; i1 = i
          endif
        enddo
        if(i1+1/=i2) Bug('Interval not found.')
        coeff(i1) = (arr(i2)-x) / (arr(i2)-arr(i1))
        coeff(i2) = (x-arr(i1)) / (arr(i2)-arr(i1))
      endif
      end

c ------------------------------

# ifdef maxfreqc_test
c Determine the maximal argument to W. (This function has been replaced by the routine freqc_range.)
c
      function maxfreqc(job1)
      use global
      use arrays
      use, intrinsic :: iso_fortran_env
      implicit none
      real_dp                    :: maxfreqc
      type(jobtype), intent(in)  :: job1
      real_dp,       allocatable :: wintgrc(:,:,:,:)
      real_dp                    :: eneq,enek,enediff,weight,rdum1,rdum2
      real_dp                    :: freqr1(nfreqr)
      integer                    :: minb(nfreqr,nspin1),maxb(nfreqr,nspin1)
      integer                    :: iblock,ib,ibandq,iband,ikptq,ispin,s,ifreqr
      integer                    :: i,j,k,l,ikpt1,ikpt2
      integer                    :: kptsum

      if(oselfc==1) Error('wrong self-energy modus.')

      maxfreqc = -1d30

c
c     Loop over bands

      do iblock = 1,nblock

        ib     = block(1,iblock)
        ibandq = job1%band(ib)
        ikptq  = job1%kpt(ib)
        ispin  = job1%spin(ib)

        eneq   = ene(ibandq,ikptq,ispin)
        if(oselfc==4) then
          freqr1 = freqr
        else
          if(eneq<=efermi) then ; freqr1 = eneq - efermi - freqr
          else                  ; freqr1 = eneq - efermi + freqr
          endif
        endif

c
c       Define tetrahedron weights (only once for oselfc=1 & 4, every time for oselfc=2 & 3)
        if(oselfc/=4.or.oselfc==4.and..not.allocated(wintgrc)) then
          ! determine minimum and maximum bands
          do ifreqr = 1,nfreqr
            do s = 1,nspin1 ; if(oselfc/=4.and.s/=ispin) cycle
              ! (1) minimum
              i = 1
              do while(all(ene(i,:,s)-freqr1(ifreqr)<=efermi).and.i/=maxeband) ! all states below efermi(+freqr1) will get weight 1/nkpt
                i = i + 1
              enddo
              j = i
              do k = 1,size(ene,2) ! be sure to catch all degenerate states
                if(j<=nband(k,s)) i = min(i,deg(j,k,s))
              enddo
              minb(ifreqr,s) = i
              ! (2) maximum
              do while(any(ene(i,:,s)-freqr1(ifreqr)<efermi).and.i/=maxeband) ! all states above efermi(+freqr1) will get weight 0
                i = i + 1
              enddo
              i = i - 1
              if(i>0) then
                j = i
                do k = 1,size(ene,2) ! be sure to catch all degenerate states
                  if(j<=nband(k,s)) then
                    l = deg(j,k,s) ; if(l<j) l = deg(l,k,s)
                    i = max(i,l)
                  endif
                enddo
              endif
              maxb(ifreqr,s) = i
            enddo
          enddo
          ! calculate weights
          i = maxval(maxb-minb)+1
          if(oselfc/=4) then ; allocate ( wintgrc(size(ene,2),i,nfreqr,ispin:ispin) )
          else               ; allocate ( wintgrc(size(ene,2),i,nfreqr,nspin1) )
          endif
          if(i>0) then
            wintgrc = 0
            rdum1   = efermi
            do s = 1,nspin1 ; if(oselfc/=4.and.s/=ispin) cycle
              do ifreqr = 1,nfreqr
                efermi = rdum1 + freqr1(ifreqr)
                call tetrahedron_init(wintgrc(:nkpt,:,ifreqr,s),nkpt,minb(ifreqr,s),maxb(ifreqr,s),efermi,s,.false.)
                if(lkptadd)
     &          call tetrahedron_init(wintgrc(nkpt+1:,:,ifreqr,s),nkpt,minb(ifreqr,s),maxb(ifreqr,s),efermi,s,.true.)
              enddo
            enddo
            efermi = rdum1
          endif
        endif

        do ikpt1 = 1,nkpt
          ikpt2 = kptsum(ikpt1,ikptq)
          ! loop over frequencies and search for maximum W argument
          do ifreqr = 1,nfreqr
            do iband = 1,nband(ikpt2,ispin)
              if(iband>bando) then ; rdum1 = 0
              else                 ; rdum1 = wintgr(ikpt2,iband,ispin)
              endif
              if     (iband<minb(ifreqr,ispin)) then ; rdum2 = 1d0/nkpt
              else if(iband>maxb(ifreqr,ispin)) then ; rdum2 = 0
              else                                   ; rdum2 = wintgrc(ikpt2,iband-minb(ifreqr,ispin)+1,ifreqr,ispin)
              endif
              weight = rdum2 - rdum1
              if(abs(weight)>1d-10) then
                enek    = ene(iband,ikpt2,ispin) - efermi
                enediff = enek - freqr1(ifreqr)
                if(abs(enediff)>maxfreqc) then
c                  write(101,'(F10.5,4I3,3F15.10)') abs(enediff),ifreqr,iband,ikpt2,ispin,weight*nkpt,rdum1*nkpt,rdum2*nkpt
                  maxfreqc = abs(enediff)
                endif
              endif
            enddo
          enddo
        enddo
        if(oselfc==2.or.oselfc==3) deallocate ( wintgrc )
      enddo

      if(oselfc==4) deallocate ( wintgrc )

      end
# endif

c ------------------------------

c     Returns bounds [f1,f2] for freqc array used to interpolate W in contour integration.
c
c     ospin = 0 : Green function poles have same spin as self-energy.
c     ospin = 1 : Green function poles have opposite spin to self-energy; frequency range flipped for spin-down.
c     ospin = 2 : Green function poles have opposite spin to self-energy; frequency range flipped for spin-up.
c
c     Self-energy:  SIGMA(w) = INT G(w+w') W(-w') dw'  ! here, W is a general effective interaction
c     Residue part: SUM(km) W(w-e(km))
c     (Factors omitted)
c
      subroutine freqc_range(f1,f2,job1,ospin)
      use global
      use arrays
      use, intrinsic :: iso_fortran_env
      implicit none
      real_dp,       intent(out) :: f1,f2
      type(jobtype), intent(in)  :: job1
      integer,       intent(in)  :: ospin
      logical                    :: done(2,0:1)
      MCOMPLEX_dp                :: mom(3)
      real_dp                    :: eneq,enek,w,wght_ef,wght_w
      real_dp                    :: freqr1(nfreqr)
      real_dp                    :: logdiv1,logdiv2
      integer                    :: iblock,ib,ibandq,iband,ikptq,ispin,ispin2,ifreqr,k1,k2
      integer                    :: ikpt2
      real_dp                    :: tetrahedron_weight,logdiv

      if(oselfc==1)               Bug('wrong self-energy mode.')
      if(ospin<0.or.ospin>nspin1) Bug('ospin parameter out of range.')

      f1   =  huge(0d0)
      f2   = -huge(0d0)
      done = .false.

c
c     Loop over bands

      do iblock = 1,nblock

        ib     = block(1,iblock)
        ibandq = job1%band(ib)
        ikptq  = job1%kpt(ib)
        ispin  = job1%spin(ib)
        if(oselfc==4) then
          if(done(ispin,ikptq/(nkpt+1))) cycle
          done(ispin,ikptq/(nkpt+1)) = .true.
        endif
        if(ospin/=0) then ; ispin2 = 3 - ispin
        else              ; ispin2 =     ispin
        endif

        eneq = ene(ibandq,ikptq,ispin)
        if(oselfc==4) then
          freqr1 = freqr
        else
          if(job1%full) Error('FULL matrix calculations impossible with this self-energy type.')
          if(eneq<=efermi) then ; freqr1 = eneq - efermi - freqr
          else                  ; freqr1 = eneq - efermi + freqr
          endif
        endif

        if(ikptq==nkpt+1) then ; k1 = nkpt+1 ; k2 = nkpt+nkpti2
        else                   ; k1 = 1      ; k2 = nkpti
        endif

        if(ologdiv) then
          do ib = 1,size(block,1)
            if(block(ib,iblock)/=0) then
              ibandq = job1%band(block(ib,iblock))
              eneq   = ene(ibandq,ikptq,ispin) - efermi
              call momentum_matrix(mom,[ikptq],1,ispin,ispin,ibandq,ibandq,ibandq,ibandq,.false. MpiC(.false.) )
              logdiv1 = logdiv(eneq,dble(mom))
              do ifreqr = 1,nfreqr
                w       = freqr1(ifreqr)
                logdiv2 = logdiv(eneq-w,dble(mom))
                if( abs(logdiv1-logdiv2)>1d-8 ) then
                  f1 = min(f1,w-eneq)
                  f2 = max(f2,w-eneq)
                endif
              enddo
            endif
          enddo
        endif

        do ikpt2 = k1,k2
          do iband = 1,nband(ikpt2,ispin2)
            enek    = ene(iband,ikpt2,ispin2) - efermi
            wght_ef = 0
            if(iband<=bando) then
              wght_ef = wintgr(ikpt2,iband,ispin2) * nkpt
              if(min(wght_ef,1-wght_ef)<-1d-10) Bug('wght_ef or 1-wght_ef below zero.')
            endif
            do ifreqr = 1,nfreqr
              w = freqr1(ifreqr)
              if( w<=0.and.wght_ef>1d-8 .or. w>0.and.1-wght_ef>1d-8 ) then ! w<0 (w>0): only occ. (unocc.) states
                wght_w = tetrahedron_weight(iband,ikpt2,ispin2,w) * nkpt
                if(abs(wght_w-wght_ef)>1d-8) then
                  if(ospin==0.or.ospin==ispin) then
                    f1 = min( f1 , w - enek )
                    f2 = max( f2 , w - enek )
                  else
                    f1 = min( f1 , enek - w )
                    f2 = max( f2 , enek - w )
                  endif
                endif
              endif
            enddo
          enddo
        enddo

      enddo

      if(f1>f2) Error('Could not determine freqc range.')

      end

c ------------------------------

c     Returns integrated residues in res: INT(emin,emax) W(w') A(w') dw'.
c     A is defined on an equidistant size-ordered grid frqc(:nfrqc)
c     W is the tetrahedron weight function (from tetrahedron_nonlin) defined in (nfrq,frq,pfrq).
c
      function residues(emin,emax,ene0,nfrq,frq0,pfrq,nfrqc,frqc0,a0,linvert) result(res)
      use global, only: restart
      use, intrinsic :: iso_fortran_env
      implicit none      
      complex_dp              :: res
      logical,    intent(in)  :: linvert
      integer,    intent(in)  :: nfrqc,nfrq
      complex_dp, intent(in)  :: a0(nfrqc)
      real_dp,    intent(in)  :: emin,emax,ene0,frqc0(nfrqc),frq0(nfrq),pfrq(0:3,nfrq-1)
      complex_dp              :: a(nfrqc)
      real_dp                 :: frq(nfrq),frqc(nfrqc)
      complex_dp              :: t0,t1
      real_dp                 :: f1,f2
      integer                 :: i,j
      res = 0
      frq = frq0 + ene0
      if(linvert) then
        frqc = -frqc0(nfrqc:1:-1)
        a    =  a0   (nfrqc:1:-1)
      else
        frqc =  frqc0
        a    =  a0
      endif
      if(frq(1)>emax.or.frq(nfrq)<emin) return
      if(frqc(2)<frqc(1)) Bug('Array frqc not ordered according to size.')
      if(max(frq(1),emin)+1d-12<frqc(1)) then
        write(0,'(A)') 'residues: Out of bounds error of frequency mesh for W or T (lower bound).'
        goto 1000
      endif
      if(min(frq(nfrq),emax)-1d-12>frqc(nfrqc)) then
        write(0,'(A)') 'residues: Out of bounds error of frequency mesh for W or T (upper bound).'
        goto 1000
      endif
      res = 0
      i   = 1 ; do while(frq (i+1)<=emin)             ; i = i + 1 ; enddo
      j   = 1 ; do while(frqc(j+1)<=max(emin,frq(i))) ; j = j + 1 ; enddo
      t1  = (a(j+1)-a(j)) / (frqc(j+1)-frqc(j))
      t0  = a(j) - t1*(frqc(j)-ene0)
      do
        f1  = max( max ( frq(i),  frqc(j)   ) , emin ) - ene0
        f2  = min( min ( frq(i+1),frqc(j+1) ) , emax ) - ene0
        res = res + t0 * (    f2*(pfrq(0,i)   + f2*(pfrq(1,i)/2 + f2*(pfrq(2,i)/3 + f2*pfrq(3,i)/4 )))
     &                   -    f1*(pfrq(0,i)   + f1*(pfrq(1,i)/2 + f1*(pfrq(2,i)/3 + f1*pfrq(3,i)/4 ))) )
     &            + t1 * ( f2*f2*(pfrq(0,i)/2 + f2*(pfrq(1,i)/3 + f2*(pfrq(2,i)/4 + f2*pfrq(3,i)/5 )))
     &                   - f1*f1*(pfrq(0,i)/2 + f1*(pfrq(1,i)/3 + f1*(pfrq(2,i)/4 + f1*pfrq(3,i)/5 ))) )
        if(f1>f2) Bug('f1>f2.')
        if(frq(i+1)<frqc(j+1)) then ; i = i + 1 ; if(i==nfrq .or.frq(i) >=emax) exit
        else                        ; j = j + 1 ; if(j==nfrqc.or.frqc(j)>=emax) exit
          t1 = (a(j+1)-a(j)) / (frqc(j+1)-frqc(j))
          t0 = a(j) - t1*(frqc(j)-ene0)
        endif
      enddo
      return
 1000 if(iand(restart,R_cor)/=0.or.iand(restart,R_cot)/=0) then ! not clear at this point if this is a GW or GT calculation
        write(0,'(A)') '          This is a RESTART run. The frequency mesh was read from spex.cor (spex.cot) and is'
        write(0,'(A)') '          inappropriate for the present JOB and CONTOUR parameters.'
        Error('Out of bounds error.')
      else
        write(0,'(A)') '          This should not happen.'
        Bug('Out of bounds error.')
      endif
      end

c ------------------------------

c     Returns integrated residues in res: INT W(w') A(w') dw'.
c     A is defined by a Pade approximant A(w')=SUM resid / (w'-pole).
c     W is the tetrahedron weight function (from tetrahedron_nonlin) defined in (ene0,nfrq,frq,pfrq).
c     If lconjg, then A(-w')=conjg(A(w')), otherwise A(-w')=A(w')     
c
      function residues_pade(emin,emax,ene0,nfrq,frq0,pfrq,npol,pole,resid) result(res)
      use, intrinsic :: iso_fortran_env
      implicit none
      complex_dp              :: res
      integer,    intent(in)  :: npol,nfrq
      complex_dp, intent(in)  :: pole(npol),resid(npol)
      real_dp,    intent(in)  :: emin,emax,ene0,frq0(nfrq),pfrq(0:3,nfrq-1)
      real_dp                 :: frq(nfrq)
      complex_dp              :: q(0:3),p,r
      real_dp                 :: f1,f2
      integer                 :: i,j
      frq = frq0 + ene0
      if(frq(1)>emax.or.frq(nfrq)<emin) return
      if(emin*emax<0) Bug('integral boundaries have different sign.')
      res = 0
      i   = 1 ; do while(frq(i+1)<=emin) ; i = i + 1 ; enddo
      do
        f1 = max( frq(i)   , emin ) - ene0
        f2 = min( frq(i+1) , emax ) - ene0 ; if(f1>f2) Bug('f1>f2.')
        do j = 1,npol
          p    = pole(j) - ene0
          r    = resid(j)
          q(0) = log(f2-p) - log(f1-p)
          q(1) = q(0) * p + (f2-f1)
          q(2) = q(1) * p + (f2**2-f1**2)/2
          q(3) = q(2) * p + (f2**3-f1**3)/3
          res  = res + r * sum(pfrq(:,i)*q)
        enddo
        i = i + 1 ; if(i==nfrq.or.frq(i)>=emax) exit
      enddo
      end

c ------------------------------
