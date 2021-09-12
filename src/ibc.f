c--------------------------------------------------------------------------------
c Copyright (c) 2021 Peter Grünberg Institut, Forschungszentrum Jülich, Germany
c This file is part of SPEX and available as free software under the conditions
c of the MIT license as expressed in the LICENSE file in more detail.
c--------------------------------------------------------------------------------

# include "cppmacro.h"

c Incomplete basis-set correction for
c - susceptibility (routines ibc1_contraction, ibc1_doublecount)
c - self-energy    (routines ibc2_contraction, ibc2_doublecount1, ibc2_doublecount2)

c -----------------------

c Until we know what the optimal "gauge" for the response is, we allow for different choices with the preprocessor identifier GAUGE.
c In general, an ll perturbation of the function u(l,El) (El is the energy parameter) leads to a response in all l1 channels with l1=|l-ll|..l+ll .
c The response is not uniquely defined as we may always add a multiple of the homogeneous solution hom(l1,El-w) (w is the frequency).
c Furthermore, we want the response to go to zero (in value and radial derivative) at the MT sphere boundary. We use two gauge definitions:
c (A) Match with hom(l1,El-w) and u(l1,El1) [ use udot(l1,El1) if El-w=El1 ]
c (B) Match with hom(l1,El-w) and u(l1,El)  [ use udot(l1,El)  if    w=0   ] -- this requires u(l1,El) [and udot(l1,El)] to be calculated.
c
c Identifier GAUGE:
c 1 - Always use (A)
c 2 - Use (A) and (B) for the response of u/udot and ulo, respectively.
c 3 - Always use (B).

# define GAUGE 2

c -----------------------

c     FOR SUSCEPTIBILITY
c
c
c     (A) Contraction term
c
c     SUM(n=occ) <J n | tilde(nI,w) >
c
      subroutine ibc1_contraction(frq,nfrq,defsuscep)

      use global
      use, intrinsic :: iso_fortran_env
      implicit none
      integer,    intent(in) :: nfrq
      logical,    intent(in) :: defsuscep
      complex_dp, intent(in) :: frq(nfrq)
      complex_dp             :: efrq
      complex_dp             :: h1(maxgrid,maxindx),h2(maxgrid,maxindx),dh(maxindx)
      complex_dp             :: h1_(maxgrid,2),h2_(maxgrid,2),dh_(2)
      complex_dp             ::  hom1(maxgrid), hom2(maxgrid),dhom
      complex_dp             :: ihom1(maxgrid),ihom2(maxgrid),dihom
      complex_dp             ::  rhs1(maxgrid), rhs2(maxgrid)
      complex_dp             :: ihom1_(maxgrid),ihom2_(maxgrid)
      complex_dp             :: c0,c1,c2,c3
      complex_dp             :: cont(maxval(neq),maxlmindx,(2*maxlcut+1)*maxindx)
      complex_dp             :: cmt1(bando,nkpt),carr(maxval(neq))
      complex_dp             :: hom1t(2),hom2t(2),ihom1t(2),ihom2t(2)
      real_dp,   allocatable :: b1(:,:,:),b2(:,:,:),db(:,:)
      real_dp                :: b1_(maxgrid,2),b2_(maxgrid,2)
      real_dp                :: p1(maxgrid,maxindx),p2(maxgrid,maxindx),dp(maxindx),ep(maxindx)
      real_dp                :: pd1(maxgrid,maxindx),pd2(maxgrid,maxindx),dpd(maxindx)
      real_dp                :: tild1(maxgrid),tild2(maxgrid)
      real_dp                :: help1(maxgrid),help2(maxgrid)
      real_dp                :: cor1(maxgrid),cor2(maxgrid),ecor
      real_dp                :: wintgr1(bando,nkpt),int(max(maxindxm,maxindx))
      real_dp                :: gaunt1(0:maxlcutm,0:maxlcut,-maxlcutm:maxlcutm,-maxlcut:maxlcut)
      real_dp                :: debas,rdum,rdum1,squarec
      real_dp                :: rg,rt
      real_dp,   parameter   :: delta = 0.005d0
      integer                :: ng,nt
      integer                :: ifrq,ispin,itype,ieq,ic
      integer                :: l,l0,l1,ll,ll0,n,n0,n1,nn,nn0,m,m0,m1,mm,mm0,lm,lm0,llm,lm_0,l2
      integer                :: i,j,icount
      real                   :: time1,time2
      real_dp                :: intgrf,gaunt
      complex_dp             :: cintgrf

# ifdef MPI
      Error('IBC & MPI not implemented.')
# endif

      icount = 0

      if(oibc==0)   Error('IBC undefined.')
      if(storeibz)  Error('not implemented')
      if(l_soc)     Error('not implemented')
      if(defsuscep) then
        write(6,'(/A)') 'IBC susceptibility contraction'
      else
        write(6,'(A'NoA) 'IBC projections... '
      endif

c      write(140,'(2F20.10)') (rgrid(i,1),vmt(i,1,1,1),i=1,grid(1)%number)

c      write(*,*) 'out-projecting of core states switched off'

# if 0
      write(*,*) 'BASM ORDER CHANGED'

      do i = 1,grid(1)%number
        tild1(i) = 1
        if(rgrid(i,1)>1) then
          rdum     = (rgrid(i,1)-1)/(grid(1)%radius-1)
          tild1(i) = ( 2*rdum**3 - 3*rdum**2 + 1 )
        endif
      enddo

      do l = 0,maxlcutm
        do i = 1,maxgrid
c          basm(i,:nindxm(l,1),l,1) = basm(i,nindxm(l,1):1:-1,l,1) !* tild1(i)
          basm(i,:nindxm(l,1),l,1) = basm(i,:nindxm(l,1),l,1) * tild1(i)
c          write(777,'(F15.10'NoA) rgrid(i,1),basm(i,:nindxm(l,1),l,1)/rgrid(i,1)
c          write(777,*)
        enddo
      enddo
# endif

c
c     Allocate array for contraction contribution (ibc_suscep0) and array for double-counting (ibc_proj)
      lm = 0
      nn = 0
      n  = 0
      do itype = 1,ntype
        lm0 = sum ( [ ((2*l2+1)*nindxm(l2,itype),l2=0,lcutm(itype)) ] ) ; if(lm0>lm) lm = lm0
        lm0 = sum ( [ (         nindxm(l2,itype),l2=0,lcutm(itype)) ] ) ; if(lm0>nn) nn = lm0
        lm0 = sum ( [ (         nindx (l2,itype),l2=0,lcut (itype)) ] ) ; if(lm0>n)  n  = lm0
      enddo
      if(defsuscep) then
        allocate ( ibc_suscep0(lm,lm,ncent,nfrq) )
        ibc_suscep0 = 0
      endif
      allocate ( ibc_proj(n,n,nn,ntype,nfrq,nspin) )
      ibc_proj   = 0
      if(maxlmindxm/=lm) Error('!')

c
c     Core contribution

      if(oibc==2.and.defsuscep) then
        if(lcore_soc) Error('not implemented')
        call cpu_time(time1)
        do ispin = 1,nspin
          if(nspin==2) then
            if(ispin==1) write(6,'(A)') ' Spin up'
            if(ispin==2) write(6,'(A)') ' Spin down'
          endif
          do itype = 1,ntype
            do l = 0,lcutc(itype)
              do n = 1,nindxc(l,itype)
                ! SRA cores
                ecor = ecore(n,l,itype,ispin)
                call corestate(cor1,cor2,ecor,l,itype,ispin)
c                cor1 = core1(:,n,l,itype,ispin)
c                cor2 = core2(:,n,l,itype,ispin)
                squarec = intgrf ( cor1**2 + cor2**2 ,itype )
                write(6,'(1X,2I3,A,F13.6,A,F13.6,A,F9.7,A,F9.7,A)') itype,n+l,lchar(l),ecore(n,l,itype,ispin),'  ->',ecor,
     &            '    (',intgrf(core1(:,n,l,itype,ispin)**2+core2(:,n,l,itype,ispin)**2,itype),',',squarec,')'
                call core_turning_point(rt,nt,l,ecor,itype,ispin)
                write(*,*) nt,rt
                do l1 = 0,lcut(itype)
                  if(abs(l-l1)>lcutm(itype)) cycle
c                  rewind(120)
c                  write(120,'(3F20.15)') (rgrid(i,itype),core1(i,n,l,itype,ispin),
c     &              core2(i,n,l,itype,ispin),i=1,grid(itype)%number) ; read(*,*)
c                  cycle
                  do ifrq = 1,nfrq!,1,-1
                    ! homogeneous solution
                    efrq = ecor + frq(ifrq) ; write(*,*) 'efrq definition changed to ecor + frq(ifrq)'
                    call core_dirac_hom_out(hom1,hom2,hom1t(1),hom2t(1),rt,nt,l1,efrq,itype,ispin)
                    call core_dirac_hom_in (hom1,hom2,hom1t(2),hom2t(2),rt,nt,l1,efrq,itype,ispin)
c                    write(*,*) 'hom'
c                    write(*,*) hom1t
c                    write(*,*) hom2t
                    do ll = abs(l-l1),min(l+l1,lcutm(itype)),2
c                      write(*,*) ll
                      do nn = 1,nindxm(ll,itype)
                        ! inhomogeneous solution
                        if(l1==l) then
                          debas = intgrf ( ( cor1**2 + cor2**2 ) * basm(:,nn,ll,itype) / rgrid(:,itype) ,itype) / squarec
                          help1 = -basm(:,nn,ll,itype) + debas * rgrid(:,itype)
                        else
                          help1 = -basm(:,nn,ll,itype)
                        endif
                        call core_dirac_inhom_out(ihom1,ihom2,ihom1t(1),ihom2t(1),rt,nt,n+l-1,l1,efrq,cor1,cor2,help1,itype,ispin)
                        call core_dirac_inhom_in (ihom1,ihom2,ihom1t(2),ihom2t(2),rt,nt,n+l-1,l1,efrq,cor1,cor2,help1,itype,ispin)
c                        cycle
c                        write(*,*) 'ihom'
c                        write(*,*) ihom1t
c                        write(*,*) ihom2t
c                        read(*,*)
c                        write(1002,'(9E40.20)') (rgrid(i,itype),ihom1(i),ihom2(i),hom1(i),hom2(i),i=1,grid(itype)%number) !; Error(' ')
                        ! match outward and inward solution
                        if(l==l1.and.frq(ifrq)==0) then;c0 = 0
                          c1        = 1d0
                          c2        = ( ihom1t(1)+hom1t(1) - ihom1t(2) ) / hom1t(2)
c                          write(*,'(A$)') 'a'
                        else
                          ihom1t(1) = ihom1t(1) - ihom1t(2)
                          ihom2t(1) = ihom2t(1) - ihom2t(2)
                          c0        = - (  hom1t(1) * hom2t(2) -  hom2t(1) * hom1t(2) )
                          c1        =   ( ihom1t(1) * hom2t(2) - ihom2t(1) * hom1t(2) ) / c0
                          c2        =   ( ihom1t(1) * hom2t(1) - ihom2t(1) * hom1t(1) ) / c0
c                          write(*,'(A$)') 'b'
                        endif
                        n0         = min(nt-1,grid(itype)%number)
                        ihom1(:n0) = ihom1(:n0) + c1 * hom1(:n0)
                        ihom2(:n0) = ihom2(:n0) + c1 * hom2(:n0)
                        if(nt<=grid(itype)%number) then
                          ihom1(nt:) = ihom1(nt:) + c2 * hom1(nt:)
                          ihom2(nt:) = ihom2(nt:) + c2 * hom2(nt:)
                        endif
                        ! orthogonalize to core state (l=l1)
                        if(l==l1) then
                          c3    = cintgrf( ihom1*cor1 + ihom2*cor2 ,itype) / squarec
c                          write(1002,'(9F40.20)') (rgrid(i,itype),ihom1(i),ihom2(i),hom1(i),hom2(i),i=1,grid(itype)%number)
                          ihom1 = ihom1 - c3 * cor1
                          ihom2 = ihom2 - c3 * cor2!;stop
                          write(*,'(6F30.10)') imag(frq(ifrq)),abs(c3),abs(c1),abs(c2),abs(c0),intgrf(abs(ihom1)**2+abs(ihom2)**2,
     &                      itype)
                        endif
c                        rewind(120)
                        write((nindxm(ll,itype)-nn+1)*100+ll*10+l1,'(9F40.15)')
     &                    (rgrid(i,itype),ihom1(i),ihom2(i),hom1(i),hom2(i),i=1,grid(itype)%number)
                        write((nindxm(ll,itype)-nn+1)*100+ll*10+l1,*)
                        write((nindxm(ll,itype)-nn+1)*100+ll*10+l1,*)
                        ! add negative frequency (make real)
                        tild1 = 2*real(ihom1)
                        tild2 = 2*real(ihom2)
                        ! Calculate contraction SUM(n) <M(i)|cor(n)*tild(nj)> (->ibc_suscep0)
                        help1 = cor1 * tild1 / rgrid(:,itype)
                        help2 = cor2 * tild2 / rgrid(:,itype)
                        do ll0 = abs(l-l1),min(l+l1,lcutm(itype)),2
                          do nn0 = 1,nindxm(ll0,itype)
                            int(nn0) = intgrf(basm(:,nn0,ll0,itype)*(help1+help2),itype)
                          enddo
                          do mm = -min(ll0,ll),min(ll0,ll)
                            do m = -l,l
                              m1   = mm + m ; if(abs(m1)>l1) cycle
                              rdum = gaunt(ll0,l1,l,mm,m1,m) * gaunt(ll,l1,l,mm,m1,m)
                              i    = (ll0+mm)*nindxm(ll0,itype)      + sum( [ ((2*l2+1)*nindxm(l2,itype),l2=0,ll0-1) ] )
                              j    = (ll +mm)*nindxm(ll, itype) + nn + sum( [ ((2*l2+1)*nindxm(l2,itype),l2=0,ll -1) ] )
                              do nn0 = 1,nindxm(ll0,itype)
                                rdum1 = rdum * int(nn0)
                                i     = i + 1
                                ic    = sum(neq(:itype-1))
                                do ieq = 1,neq(itype)
                                  ic                        = ic + 1
                                  ibc_suscep0(i,j,ic,ifrq) = ibc_suscep0(i,j,ic,ifrq) + rdum1
                                enddo
                              enddo
                            enddo
                          enddo
                        enddo

                      enddo
                    enddo
                  enddo
                enddo
              enddo
            enddo
          enddo
        enddo
        call cpu_time(time2)
        write(6,'(A,F10.5)') '  Timing:',time2-time1
        write(*,*) sum(ibc_suscep0)
      endif
c      return
c      Error(' ')

c
c     Valence contribution

      if(defsuscep) write(6,'(A'NoA) '  Valence contribution... '
      call cpu_time(time1)

      do ispin = 1,nspin
        wintgr1 = transpose ( wintgr(:nkpt,:bando,ispin) )
        do itype = 1,ntype
          ng = grid(itype)%number
          rg = grid(itype)%radius

c         Recalculate radial functions (bas1/2) to adapt to our radial integrators (->b1,b2,db)
          allocate ( b1(maxgrid,maxval(nindx(:,itype)),0:lcut(itype)) )
          allocate ( b2(maxgrid,maxval(nindx(:,itype)),0:lcut(itype)) )
          allocate ( db(maxval(nindx(:,itype)),0:lcut(itype)) )
          do l = 0,lcut(itype)
            do n = 1,nindx(l,itype)
              if(n/=2) then
                call dirac_hom(b1(:,n,l),b2(:,n,l),db(n,l),l,ebas(n,l,itype,ispin),itype,ispin)
              else
                call dirac_inhom(b1(:,2,l),b2(:,2,l),db(2,l),l,ebas(1,l,itype,ispin),
     &                           b1(:,1,l),b2(:,1,l),itype,ispin)
                rdum      = intgrf(b1(:,1,l)*b1(:,2,l)+b2(:,1,l)*b2(:,2,l),itype)
                b1(:,2,l) = b1(:,2,l) - rdum * b1(:,1,l)
                b2(:,2,l) = b2(:,2,l) - rdum * b2(:,1,l)
                db(2,l)   = db(2,l)   - rdum * db(1,l)
              endif
            enddo
          enddo

          do l = 0,lcut(itype)

            ! calculate radial functions with "shifted" energy parameters (+-delta) for udot response
            call dirac_hom(b1_(:,1),b2_(:,1),c0,l,ebas(1,l,itype,ispin)+delta,itype,ispin) ! c0 is complex but real in dirac.f!
            call dirac_hom(b1_(:,2),b2_(:,2),c0,l,ebas(1,l,itype,ispin)-delta,itype,ispin)

            ! contract coefficients (->cont)
            lm_0 = sum ( [ ((2*l2+1)*nindx(l2,itype),l2=0,l-1) ] )
            ic   = sum(neq(:itype-1))
            do ieq = 1,neq(itype)
              ic = ic + 1
              lm = 0
              do m = -l,l
                do n = 1,nindx(l,itype)
                  lm   = lm + 1
                  cmt1 = wintgr1 * cmt(lm_0+lm,ic,:bando,:nkpt,ispin)
                  lm0  = 0
                  do l0 = 0,lcut(itype)
                    do m0 = -l0,l0
                      do n0 = 1,nindx(l0,itype)
                        lm0              = lm0 + 1
                        cont(ieq,lm0,lm) = sum(conjg(cmt(lm0,ic,:bando,:nkpt,ispin))*cmt1)
                      enddo
                    enddo
                  enddo
                enddo
              enddo
            enddo

            do l1 = 0,lcut(itype)
              if(abs(l-l1)>lcutm(itype)) cycle
              ! precalculate Gaunt coefficients (->gaunt1)
              gaunt1 = 0
              do l0 = 0,lcut(itype)
                do ll0 = abs(l1-l0),min(l1+l0,lcutm(itype)),2
                  do m0 = -l0,l0
                    do mm0 = -ll0,ll0
                      m1                    = mm0 + m0 ; if(abs(m1)>l1) cycle
                      gaunt1(ll0,l0,mm0,m0) = gaunt(ll0,l1,l0,mm0,m1,m0)
                    enddo
                  enddo
                enddo
              enddo

              ! construct the radial response functions (->tild1,tild2)

              do ifrq = 1,nfrq
                ! solve homogeneous differential equation (-> h1,h2,dh -> hom1,hom2)
                do n = 1,nindx(l,itype)
                  efrq = ebas(n,l,itype,ispin) + frq(ifrq)
                  call cdirac_hom(h1(:,n),h2(:,n),dh(n),l1,efrq,itype,ispin)
                enddo
                ! solve homogeneous differential equation for udot (->h1_,h2_,dh_)
                efrq = ebas(1,l,itype,ispin) + frq(ifrq)
                call cdirac_hom(h1_(:,1),h2_(:,1),dh_(1),l1,efrq+delta,itype,ispin)
                call cdirac_hom(h1_(:,2),h2_(:,2),dh_(2),l1,efrq-delta,itype,ispin)
# if   GAUGE == 1
#   define NCOPY nindx(l,itype)
# elif GAUGE == 2
#   define NCOPY 2
# else
#   define NCOPY 0
# endif
                ! define radial functions for matching (->p1,p2,dp,pd1...)
                if(ifrq==1) then ; if(frq(ifrq)/=0) Error('!')
                  do n = 1,NCOPY
                    p1(:,n) = b1(:,1,l1) ; pd1(:,n) = b1(:,2,l1)
                    p2(:,n) = b2(:,1,l1) ; pd2(:,n) = b2(:,2,l1)
                    dp(n)   = db(1,l1)   ; dpd(n)   = db(2,l1)
                    ep(n)   = ebas(1,l1,itype,ispin)
                  enddo
                  do n = NCOPY+1,nindx(l,itype)
                    p1(:,n) = h1(:,n)
                    p2(:,n) = h2(:,n)
                    dp(n)   = dh(n)
                    call dirac_inhom(pd1(:,n),pd2(:,n),dpd(n),l1,ebas(n,l,itype,ispin),p1(:,n),p2(:,n),itype,ispin)
                    rdum     = intgrf(p1(:,n)*pd1(:,n)+p2(:,n)*pd2(:,n),itype)
                    pd1(:,n) = pd1(:,n) - rdum * p1(:,n)
                    pd2(:,n) = pd2(:,n) - rdum * p2(:,n)
                    dpd(n)   = dpd(n)   - rdum * dp(n)
                    ep(n)    = ebas(n,l,itype,ispin)
                  enddo
                endif

                do ll = abs(l-l1),min(l+l1,lcutm(itype)),2 ! Gaunt condition l+ll+l1=even
                  do nn = 1,nindxm(ll,itype)

                    do n = 1,nindx(l,itype)
                      efrq = ebas(n,l,itype,ispin) + frq(ifrq)

                      if(n/=2) then
                        ! (A) RESPONSE FOR u AND ulo
                        ! right-hand side (->rhs1,rhs2)
                        rhs1  = -b1(:,n,l) * basm(:,nn,ll,itype) / rgrid(:,itype)
                        rhs2  = -b2(:,n,l) * basm(:,nn,ll,itype) / rgrid(:,itype)
                        if(l1==l) then
                          debas = intgrf( (b1(:,n,l)**2+b2(:,n,l)**2) * basm(:,nn,ll,itype) / rgrid(:,itype) ,itype)
                          rhs1  = rhs1 + b1(:,n,l) * debas
                          rhs2  = rhs2 + b2(:,n,l) * debas
                        endif
                        ! solve inhomogeneous differential equation (->ihom1,ihom2)
                        call cdirac_inhom(ihom1,ihom2,dihom,l1,efrq,rhs1,rhs2,itype,ispin)
                        ! match with uhom and u (or udot and u)
                        if(abs(efrq-ep(n))>1d-8) then
                          hom1 = h1(:,n)
                          hom2 = h2(:,n)
                          dhom = dh(n)
                          call continuous(c0,c1,ihom1(ng)/rg,dihom,hom1(ng)/rg,dhom,p1(ng,n)/rg*(1d0,0d0),dp(n)*(1d0,0d0))
                          ihom1 = ihom1 + c0 * hom1 + c1 * p1(:,n)
                          ihom2 = ihom2 + c0 * hom2 + c1 * p2(:,n)
                        else
                          call continuous(c0,c1,ihom1(ng)/rg,dihom,pd1(ng,n)/rg*(1d0,0d0),dpd(n)*(1d0,0d0),
     &                                                              p1(ng,n)/rg*(1d0,0d0), dp(n)*(1d0,0d0))
                          ihom1 = ihom1 + c0 * pd1(:,n) + c1 * p1(:,n)
                          ihom2 = ihom2 + c0 * pd2(:,n) + c1 * p2(:,n)
                        endif
                      else
                        ! (B) RESPONSE FOR udot (by finite difference)
                        rdum = delta
                        j    = 1
  1                     rhs1 = -b1_(:,j) * basm(:,nn,ll,itype) / rgrid(:,itype)
                        rhs2 = -b2_(:,j) * basm(:,nn,ll,itype) / rgrid(:,itype)
                        if(l1==l) then
                          debas = intgrf( (b1_(:,j)**2+b2_(:,j)**2) * basm(:,nn,ll,itype) / rgrid(:,itype) ,itype)
                          rhs1  = rhs1 + b1_(:,j) * debas
                          rhs2  = rhs2 + b2_(:,j) * debas
                        endif
                        call cdirac_inhom(ihom1,ihom2,dihom,l1,efrq+rdum,rhs1,rhs2,itype,ispin)
                        if(abs(efrq+rdum-ep(1))>1d-8) then
                          hom1 = h1_(:,j)
                          hom2 = h2_(:,j)
                          dhom = dh_(j)
                          call continuous(c0,c1,ihom1(ng)/rg,dihom,hom1(ng)/rg,dhom,p1(ng,1)/rg*(1d0,0d0),dp(1)*(1d0,0d0))
                          ihom1 = ihom1 + c0 * hom1 + c1 * p1(:,1)
                          ihom2 = ihom2 + c0 * hom2 + c1 * p2(:,1)
                        else
                          call continuous(c0,c1,ihom1(ng)/rg,dihom,pd1(ng,1)/rg*(1d0,0d0),dpd(1)*(1d0,0d0),
     &                                                              p1(ng,1)/rg*(1d0,0d0), dp(1)*(1d0,0d0))
                          ihom1 = ihom1 + c0 * pd1(:,1) + c1 * p1(:,1)
                          ihom2 = ihom2 + c0 * pd2(:,1) + c1 * p2(:,1)
                        endif
                        if(rdum>0) then
                          rdum   = -delta
                          j      = 2
                          ihom1_ = ihom1
                          ihom2_ = ihom2
                          goto 1
                        endif
                        ihom1 = (ihom1_-ihom1) / (2*delta)
                        ihom2 = (ihom2_-ihom2) / (2*delta)
                      endif
                      ! add negative frequency, make real (->tild1,tild2)
                      tild1 = 2*real(ihom1)
                      tild2 = 2*real(ihom2)

c                      tild1 = bas1(:,n,l,itype,ispin) * basm(:,nn,ll,itype) / rgrid(:,itype)
c                      tild2 = bas2(:,n,l,itype,ispin) * basm(:,nn,ll,itype) / rgrid(:,itype)
c                      help1 = tild1 ; tild1 = 0
c                      help2 = tild2 ; tild2 = 0
c                      if(l1==1) then
c                      do n0 = 1,1!nindx(l1,itype)
c                        c0    = intgrf(bas1(:,n0,l1,itype,ispin)*help1+bas2(:,n0,l1,itype,ispin)*help2,itype)
c                        c1    = intgrf(bas1(:,n0,l1,itype,ispin)**2   +bas2(:,n0,l1,itype,ispin)**2,   itype)
c                        tild1 = tild1 + bas1(:,n0,l1,itype,ispin) * c0/c1
c                        tild2 = tild2 + bas2(:,n0,l1,itype,ispin) * c0/c1
c                      enddo
c                      endif

                      ! Precalculate projections <u|u^tilde> for double-counting term (->ibc_proj)
                      lm  = sum(nindx (:l -1,itype))
                      lm0 = sum(nindx (:l1-1,itype))
                      llm = sum(nindxm(:ll-1,itype))
                      do n1 = 1,nindx(l1,itype)
                        ibc_proj(lm0+n1,lm+n,llm+nn,itype,ifrq,ispin) =
     &                    intgrf(bas1(:,n1,l1,itype,ispin)*tild1+bas2(:,n1,l1,itype,ispin)*tild2,itype)
                      enddo
                      if(.not.defsuscep) cycle

                      ! project core states out (->tild1,tild2)
                      if(l1<=lcutc(itype)) then
                        do n0 = 1,nindxc(l1,itype)
                          rdum  = intgrf ( core1(:,n0,l1,itype,ispin)*tild1 + core2(:,n0,l1,itype,ispin)*tild2 , itype)
                          tild1 = tild1 - rdum * core1(:,n0,l1,itype,ispin)
                          tild2 = tild2 - rdum * core2(:,n0,l1,itype,ispin)
                        enddo
                      endif

                      ! Calculate contraction <M(i)|sum(nk)phi(nk)phi^tilde(nkj)> (->ibc_suscep0)
                      do l0 = 0,lcut(itype)
                        do ll0 = abs(l0-l1),min(l0+l1,lcutm(itype)),2
                          do n0 = 1,nindx(l0,itype)
                            help1 = bas1(:,n0,l0,itype,ispin) * tild1 / rgrid(:,itype)
                            help2 = bas2(:,n0,l0,itype,ispin) * tild2 / rgrid(:,itype)
                            do nn0 = 1,nindxm(ll0,itype)
                              int(nn0) = intgrf(basm(:,nn0,ll0,itype)*(help1+help2),itype)
                              icount = icount + 1
                            enddo

                            do mm = -ll,ll
                              do mm0 = -ll0,ll0
                                carr = 0
                                do m1 = -l1,l1
                                  m    = m1 - mm  ; if(abs(m) >l ) cycle
                                  m0   = m1 - mm0 ; if(abs(m0)>l0) cycle
                                  lm0  = (l0+m0)*nindx(l0,itype) + n0 + sum( [ ((2*l2+1)*nindx(l2,itype),l2=0,l0-1) ] )
                                  lm   = (l +m )*nindx(l ,itype) + n
                                  rdum = gaunt1(ll,l,mm,m) * gaunt1(ll0,l0,mm0,m0)
                                  carr = carr + rdum*cont(:,lm0,lm)
                                enddo
                                if(all(carr==0)) cycle
                                i = (ll0+mm0)*nindxm(ll0,itype)      + sum( [ ((2*l2+1)*nindxm(l2,itype),l2=0,ll0-1) ] )
                                j = (ll +mm )*nindxm(ll, itype) + nn + sum( [ ((2*l2+1)*nindxm(l2,itype),l2=0,ll -1) ] )
                                do nn0 = 1,nindxm(ll0,itype)
                                  i  = i + 1
                                  ic = sum(neq(:itype-1))
                                  do ieq = 1,neq(itype)
                                    ic                       = ic + 1
                                    ibc_suscep0(i,j,ic,ifrq) = ibc_suscep0(i,j,ic,ifrq) + carr(ieq) * int(nn0)
                                  enddo
                                enddo
                              enddo
                            enddo

                          enddo
                        enddo
                      enddo

                    enddo
                  enddo
                enddo
              enddo
            enddo
          enddo
          deallocate ( b1,b2,db )
        enddo
      enddo

      call cpu_done(time1)
      if(defsuscep) write(*,*) sum(ibc_suscep0) !; Error(' ')

      if(nspin==1.and.defsuscep) ibc_suscep0 = ibc_suscep0 * 2 ! spin factor

# if 0
      write(*,*) 'change cmt cpw'
      cpw = 0
      cmt = 0
      j = 0
      ic = 0
      do itype = 1,ntype
        do ieq = 1,neq(itype)
          ic = ic + 1
          lm = sum ( [ ((2*l+1)*nindx(l,itype),l=0,lcut(itype)) ] )
          do i = 1,lm
            j = j + 1
            cmt(i,ic,j,:,1) = 1
          enddo
        enddo
      enddo
# endif

# if 0
c      return
      l = maxlmindxm
      do ifrq = 1,nfrq
        do j = 1,maxlmindxm
          do i = 1,maxlmindxm

c          do ifrq = 1,nfrq
            write(900,'(2I4,2F15.10)') i,j,ibc_suscep0(i,j,1,ifrq)
            write(901,'(2I4,2F15.10)') i,j,conjg(ibc_suscep0(j,i,1,ifrq))
c            write(10000+i*100+j,'(5F20.10)') imag(frq(ifrq)),ibc_suscep0(i,j,:,ifrq)
          enddo
        enddo
        write(900,*)
        write(901,*)
      enddo
c      write(999,'(2F20.10)') (rgrid(i,1),basm(i,2,0,1),i=1,grid(1)%number)
c      Error(' ')
# endif
      contains

c     Returns coefficients a and b so that uu + a*u + b*udot is continuous at the MT sphere boundary (uu and duu are boundary values)
      subroutine continuous(c1,c2,uu,duu,u1,du1,u2,du2)
      implicit none
      complex_dp, intent(out) :: c1,c2
      complex_dp, intent(in)  :: uu,duu,u1,du1,u2,du2
      complex_dp              :: wronsk
      wronsk = u1 * du2 - u2 * du1
      c1     = - ( uu * du2 - duu * u2 ) / wronsk
      c2     =   ( uu * du1 - duu * u1 ) / wronsk
      end subroutine continuous

      end
# undef NCOPY

c -----------------------

c     (B) Double-counting contribution
c
c     - SUM(m=all) SUM(n=occ.) <J n|m> <m|tilde(nI,w)>


# ifdef FIBC
      subroutine ibc1_doublecount(ibc_suscep1,cprod,dim,ikpt,kpt1,nkpt1,nkpts,nfrq,bo1,bo2,bu1,bu2)
# else
      subroutine ibc1_doublecount(ibc_suscep1,cprod,    ikpt,kpt1,nkpt1,nkpts,nfrq,bo1,bo2,bu1,bu2)
# endif

      use global
      use wrapper, only: dotprod
      use, intrinsic :: iso_fortran_env
      implicit none
      integer,     intent(in)    :: bo1,bo2,bu1,bu2,ikpt,nkpt1,kpt1(nkpt1),nkpts(nkpt1),nfrq
# ifdef FIBC
      integer,     intent(in)    :: dim
      complex_dp,  intent(inout) :: ibc_suscep1(dim,nbasp,nfrq)
      MCOMPLEX_dp, intent(in)    :: cprod(bu1:bu2,bo1:bo2,nkpt1,dim)
      complex_dp                 :: cprod0(bu1:bu2)
# else
      complex_dp,  intent(inout) :: ibc_suscep1(maxlmindxm,maxlmindxm,ncent,nfrq)
      MCOMPLEX_dp, intent(in)    :: cprod(bu1:bu2,bo1:bo2,nkpt1,nbasp)
      complex_dp                 :: cprod0(bu1:bu2,bo1:bo2,nkpt1,maxlmindxm)
# endif
      complex_dp,  allocatable   :: cprod1(:,:),cprod2(:,:)
      complex_dp                 :: cmt1((bo2-bo1+1)*nkpt1,maxindx)
      complex_dp                 :: cmt2(bu2-bu1+1,maxlmindx,nkpt1)
      complex_dp                 :: cexp
      real_dp                    :: wintgr1(bo1:bo2,nkpt)
      real_dp                    :: rdum
      integer                    :: ifrq,ispin,itype,ieq,ic,ikpt1,ikpt2,ib1,k
      integer                    :: l,l1,ll,n,n1,nn,m,m1,mm,lm,lm1,l2,nlm
      integer                    :: ln,ln1,lln,ln0
      integer                    :: i,j,i0,ni
      integer                    :: kptsum
      real_dp                    :: gaunt

      if(storeibz) Error('not implemented')

      do ispin = 1,nspin
        wintgr1 = transpose ( wintgr(:nkpt,bo1:bo2,ispin) )
        i0      = 0
        ic      = 0
        do itype = 1,ntype
          do ieq = 1,neq(itype)
            ic   = ic + 1
            ni   = sum ( [ ((2*l2+1)*nindxm(l2,itype),l2=0,lcutm(itype)) ] )
            ln   = sum ( [ (         nindx (l2,itype),l2=0,lcut(itype)) ] )
            nlm  = sum ( [ ((2*l2+1)*nindx (l2,itype),l2=0,lcut(itype)) ] )
            cexp = exp ( img * 2*pi * dot_product(kpt(:,ikpt),cent(:,ic)) )
            if(nspin==1) cexp = cexp * 2 ! spin factor
# ifndef FIBC
#   ifdef INV
            allocate ( cprod1(nbasp,1) )
            do k = 1,nkpt1
              do ib1 = bo1,bo2
                do i = bu1,bu2
                  cprod1(:,1)         = cprod(i,ib1,k,:nbasp)
                  call desymmetrize(cprod1,nbasp,1,1)
                  cprod0(i,ib1,k,:ni) = cprod1(i0+1:i0+ni,1)
                enddo
              enddo
            enddo
            deallocate ( cprod1 )
#   else
            cprod0(:,:,:,:ni) = cprod(:,:,:,i0+1:i0+ni)
#   endif
# endif
            allocate ( cprod1((bo2-bo1+1)*nkpt1,nlm) )
            allocate ( cprod2(ln,ln) )
            do k = 1,nkpt1
              ikpt2       = kptsum(kpt1(k),ikpt)
              cmt2(:,:,k) = transpose ( cmt(:,ic,bu1:bu2,ikpt2,ispin) )
            enddo

# ifdef FIBC
            do i = 1,dim
# else
            do i = 1,ni
# endif

              ! contract cmt(ib2,ikpt2)^* (->cprod1)
              cprod1 = 0
              j      = 0
              do k = 1,nkpt1
                ikpt1 = kpt1(k)
                do ib1 = bo1,bo2
# ifdef FIBC
                  cprod0 = cprod(:,ib1,k,i)
# endif
                  j = j + 1
                  do lm = 1,nlm
# ifdef FIBC
                    cprod1(j,lm) = nkpts(k) * wintgr1(ib1,ikpt1) * dotprod ( cmt2(:,lm,k) , cprod0 )
# else
                    cprod1(j,lm) = nkpts(k) * wintgr1(ib1,ikpt1) * dotprod ( cmt2(:,lm,k) , cprod0(:,ib1,k,i) )
# endif
                  enddo
                enddo
              enddo

# ifdef FIBC
              j = i0
# else
              j = 0
# endif
              do ll = 0,lcutm(itype)
                do mm = -ll,ll
                  ! contract cmt(ib1,ikpt1)*Gaunt (->cprod2)
                  cprod2 = 0
                  lm     = 0
                  do l = 0,lcut(itype)
                    ln0 = sum ( [ (nindx(l2,itype),l2=0,l-1) ] )
                    do m = -l,l
                      do n = 1,nindx(l,itype)
                        cmt1(:,n) = reshape ( conjg(cmt(lm+n,ic,bo1:bo2,kpt1,ispin)) , [ (bo2-bo1+1)*nkpt1 ] )
                      enddo
                      m1 = m + mm
                      do l1 = abs(ll-l),min(ll+l,lcut(itype)),2
                        if(abs(m1)>l1) cycle
                        rdum = gaunt(ll,l1,l,mm,m1,m)
                        lm1  = sum ( [ ((2*l2+1)*nindx(l2,itype),l2=0,l1-1) ] ) + (l1+m1)*nindx(l1,itype)
                        ln1  = sum ( [          (nindx(l2,itype),l2=0,l1-1) ] )
                        do n1 = 1,nindx(l1,itype)
                          ln1 = ln1 + 1
                          ln  = ln0
                          do n = 1,nindx(l,itype)
                            ln             = ln + 1
                            cprod2(ln1,ln) = cprod2(ln1,ln) + rdum * dotprod ( cmt1(:,n) , cprod1(:,lm1+n1) )
                          enddo
                        enddo
                      enddo
                      lm = lm + nindx(l,itype)
                    enddo
                  enddo
                  ! contract ibc_proj (->suscep1)
                  ln  = sum ( [ (nindx (l2,itype),l2=0,lcut(itype)) ] )
                  lln = sum ( [ (nindxm(l2,itype),l2=0,ll-1)        ] )
                  do nn = 1,nindxm(ll,itype)
                    j = j + 1
                    do ifrq = 1,nfrq
# ifdef FIBC
                      ibc_suscep1(i,j,   ifrq) = ibc_suscep1(i,j,   ifrq) - cexp *
# else
                      ibc_suscep1(i,j,ic,ifrq) = ibc_suscep1(i,j,ic,ifrq) - cexp *
# endif
     &                                            sum ( ibc_proj(:ln,:ln,lln+nn,itype,ifrq,ispin) * cprod2 )
                    enddo
                  enddo
                enddo
              enddo

            enddo
            deallocate ( cprod1,cprod2 )

            i0 = i0 + ni
          enddo
        enddo
      enddo

      if(.true.) then ! Pulay
        do ispin = 1,nspin
          do itype = 1,ntype
            do ll = 0,lcutm(itype)
              do nn = 1,nindxm(ll,itype)

              enddo
            enddo
          enddo
        enddo
      endif

# if 0
      l = maxlmindxm
      do ifrq = 1,nfrq
        do j = 1,maxlmindxm
          do i = 1,maxlmindxm

c          do ifrq = 1,nfrq
            write(910,'(2I4,2F15.10)') i,j,ibc_suscep1(i,j,1,ifrq)
            write(911,'(2I4,2F15.10)') i,j,conjg(ibc_suscep1(j,i,1,ifrq))
c            write(10000+i*100+j,'(5F20.10)') imag(frq(ifrq)),ibc_suscep0(i,j,:,ifrq)
          enddo
        enddo
        write(910,*)
        write(911,*)
      enddo

      ibc_suscep1 = ibc_suscep0 + ibc_suscep1

      l = maxlmindxm
      do ifrq = 1,nfrq
        do j = 1,maxlmindxm
          do i = 1,maxlmindxm

c          do ifrq = 1,nfrq
            write(920,'(2I4,2F15.10)') i,j,ibc_suscep1(i,j,1,ifrq)
            write(921,'(2I4,2F15.10)') i,j,conjg(ibc_suscep1(j,i,1,ifrq))
c            write(10000+i*100+j,'(5F20.10)') imag(frq(ifrq)),ibc_suscep0(i,j,:,ifrq)
            do ic = 1,ncent
              write(*,'(4F15.10'NoA) ibc_suscep1(i,j,ic,ifrq) , conjg(ibc_suscep1(j,i,ic,ifrq))
              ibc_suscep1(i,j,ic,ifrq) = avrge ( ibc_suscep1(i,j,ic,ifrq) , conjg(ibc_suscep1(j,i,ic,ifrq)) )
              ibc_suscep1(j,i,ic,ifrq) = conjg ( ibc_suscep1(i,j,ic,ifrq) )
              write(*,'(2F15.10'NoA) ibc_suscep1(i,j,ic,ifrq)
              if(abs(ibc_suscep1(i,j,ic,ifrq))>1d-8) then
                read(*,*)
              else
                write(*,*)
              endif
            enddo
          enddo
        enddo
        write(920,*)
        write(921,*)
      enddo

c      Error(' ')

      contains

      function avrge(z1,z2)
      use global, only: img
      use, intrinsic :: iso_fortran_env
      implicit none
      complex_dp             :: avrge
      complex_dp, intent(in) :: z1,z2
      real_dp                :: r1,r2,ra,rg
      r1 = real( z1 )
      r2 = real( z2 )
      ra = (r1+r2)/2
      if(r1*r2<-1d-15) then
        avrge = ra
      else
        rg = sqrt(abs(r1*r2))
        if(abs(rg-ra)>abs(-rg-ra)) rg = -rg
        avrge = rg
      endif
      r1 = imag( z1 )
      r2 = imag( z2 )
      ra = (r1+r2)/2
      if(r1*r2<-1d-15) then
        avrge = avrge + img * ra
      else
        rg = sqrt(abs(r1*r2))
        if(abs(rg-ra)>abs(-rg-ra)) rg = -rg
        avrge = avrge + img*rg
      endif
      end function avrge

# endif

      end

c -----------------------

c     FOR SELF-ENERGY
c
c
c     (A) Contraction term
c
c     SUM(IJ) INT dw <J n | tilde(n'I,w) > W(IJ,w)
c
      subroutine ibc2_contraction(job1,frq,nfrq)

      use global
      use, intrinsic :: iso_fortran_env
      implicit none
      type(jobtype)           :: job1
      integer,    intent(in)  :: nfrq
      complex_dp, intent(in)  :: frq(nfrq)
      complex_dp              :: efrq
      complex_dp              :: h1(maxgrid,maxindx),h2(maxgrid,maxindx),dh(maxindx)
      complex_dp              :: h1_(maxgrid,2),h2_(maxgrid,2),dh_(2)
      complex_dp              ::  hom1(maxgrid), hom2(maxgrid),dhom
      complex_dp              :: ihom1(maxgrid),ihom2(maxgrid),dihom
      complex_dp              ::  rhs1(maxgrid), rhs2(maxgrid)
      complex_dp              :: ihom1_(maxgrid),ihom2_(maxgrid)
      complex_dp              :: c0,c1
      complex_dp, allocatable :: screen1(:,:,:),contract(:,:,:,:),carr(:)
      real_dp,    allocatable :: b1(:,:,:),b2(:,:,:),db(:,:)
      real_dp                 :: b1_(maxgrid,2),b2_(maxgrid,2)
      real_dp                 :: p1(maxgrid,maxindx),p2(maxgrid,maxindx),dp(maxindx),ep(maxindx)
      real_dp                 :: pd1(maxgrid,maxindx),pd2(maxgrid,maxindx),dpd(maxindx)
      real_dp                 :: tild1(maxgrid),tild2(maxgrid)
      real_dp                 :: help1(maxgrid),help2(maxgrid)
      real_dp                 :: int(max(maxindxm,maxindx))
      real_dp                 :: gaunt1(0:maxlcutm,0:maxlcut,-maxlcutm:maxlcutm,-maxlcut:maxlcut)
      real_dp                 :: debas,rdum
      real_dp                 :: rg
      real_dp,    parameter   :: delta = 0.005d0
      integer                 :: ng
      integer                 :: ifrq,ispin,itype,ieq,ic
      integer                 :: l,l0,l1,ll,ll0,n,n0,nn,nn0,m,m0,m1,mm,mm0,lm,lm0,l2
      integer                 :: i,j
      integer                 :: nlmp,nlm
      integer,    allocatable :: band(:)
      integer                 :: ib,jb,iblock,nb,ikptq,iselfc
      real                    :: time1
      real_dp                 :: intgrf,gaunt,freqintegral0

      complex_dp ibc_selfc1(size(ibc_selfc),nfrq)
      ibc_selfc1 = ibc_selfc

      if(oibc==0) Error('IBC undefined.')
      if(storeibz)  Error('not implemented')
      if(l_soc)     Error('not implemented')
      write(6,'(/A)') 'IBC self-energy contraction'

      do ic = 1,ncent
        write(*,*) sum(screenk_mt(:,:,ic,:))
      enddo

# if 0
      write(*,*) 'BASM ORDER CHANGED'

      do i = 1,grid(1)%number
        tild1(i) = 1
        if(rgrid(i,1)>1) then
          rdum     = (rgrid(i,1)-1)/(grid(1)%radius-1)
          tild1(i) = ( 2*rdum**3 - 3*rdum**2 + 1 )
        endif
      enddo

      do l = 0,maxlcutm
        do i = 1,maxgrid
c          basm(i,:nindxm(l,1),l,1) = basm(i,nindxm(l,1):1:-1,l,1) !* tild1(i)
          basm(i,:nindxm(l,1),l,1) = basm(i,:nindxm(l,1),l,1) * tild1(i)
c          write(777,'(F15.10'NoA) rgrid(i,1),basm(i,:nindxm(l,1),l,1)/rgrid(i,1)
c          write(777,*)
        enddo
      enddo
# endif

      call cpu_time(time1)

      do ispin = 1,nspin
        do itype = 1,ntype
          ng = grid(itype)%number
          rg = grid(itype)%radius


c         Recalculate radial functions (bas1/2) to adapt to our radial integrators (->b1,b2,db)
          allocate ( b1(maxgrid,maxval(nindx(:,itype)),0:lcut(itype)) )
          allocate ( b2(maxgrid,maxval(nindx(:,itype)),0:lcut(itype)) )
          allocate ( db(maxval(nindx(:,itype)),0:lcut(itype)) )
          do l = 0,lcut(itype)
            do n = 1,nindx(l,itype)
              if(n/=2) then
                call dirac_hom(b1(:,n,l),b2(:,n,l),db(n,l),l,ebas(n,l,itype,ispin),itype,ispin)
              else
                call dirac_inhom(b1(:,2,l),b2(:,2,l),db(2,l),l,ebas(1,l,itype,ispin),
     &                           b1(:,1,l),b2(:,1,l),itype,ispin)
                rdum      = intgrf(b1(:,1,l)*b1(:,2,l)+b2(:,1,l)*b2(:,2,l),itype)
                b1(:,2,l) = b1(:,2,l) - rdum * b1(:,1,l)
                b2(:,2,l) = b2(:,2,l) - rdum * b2(:,1,l)
                db(2,l)   = db(2,l)   - rdum * db(1,l)
              endif
            enddo
          enddo

          nlm  = sum ( [ ((2*l2+1)*nindx (l2,itype),l2=0,lcut (itype)) ] )
          nlmp = sum ( [ ((2*l2+1)*nindxm(l2,itype),l2=0,lcutm(itype)) ] )
          allocate ( screen1(nlmp,neq(itype),-lcutm(itype):lcutm(itype)) )
          allocate ( contract(nlm,nlm,neq(itype),nfrq) )
          allocate ( carr(neq(itype)) )
          contract = 0

          do l = 0,lcut(itype)

            ! calculate radial functions with "shifted" energy parameters (+-delta) for udot response
            call dirac_hom(b1_(:,1),b2_(:,1),c0,l,ebas(1,l,itype,ispin)+delta,itype,ispin)
            call dirac_hom(b1_(:,2),b2_(:,2),c0,l,ebas(1,l,itype,ispin)-delta,itype,ispin)

            do l1 = 0,lcut(itype)
              if(abs(l-l1)>lcutm(itype)) cycle
              ! precalculate Gaunt coefficients (->gaunt1)
              gaunt1 = 0
              do l0 = 0,lcut(itype)
                do ll0 = abs(l1-l0),min(l1+l0,lcutm(itype)),2
                  do m0 = -l0,l0
                    do mm0 = -ll0,ll0
                      m1                    = mm0 + m0 ; if(abs(m1)>l1) cycle
                      gaunt1(ll0,l0,mm0,m0) = gaunt(ll0,l1,l0,mm0,m1,m0)
                    enddo
                  enddo
                enddo
              enddo
              ! construct the radial response functions (->tild1,tild2)
              do ifrq = 1,nfrq
c                if(ifrq>4) then
c                  if(ifrq<10) write(*,*) 'ifrq>4 cycled'
c                  cycle
c                endif
                ! solve homogeneous differential equation (-> h1,h2 -> hom1,hom2)
                do n = 1,nindx(l,itype)
                  efrq = ebas(n,l,itype,ispin) + frq(ifrq)
                  call cdirac_hom(h1(:,n),h2(:,n),dh(n),l1,efrq,itype,ispin)
                enddo
                ! solve homogeneous differential equation for udot (->h1_,h2_)
                efrq = ebas(1,l,itype,ispin) + frq(ifrq)
                call cdirac_hom(h1_(:,1),h2_(:,1),dh_(1),l1,efrq+delta,itype,ispin)
                call cdirac_hom(h1_(:,2),h2_(:,2),dh_(2),l1,efrq-delta,itype,ispin)

# if   GAUGE == 1
#   define NCOPY nindx(l,itype)
# elif GAUGE == 2
#   define NCOPY 2
# else
#   define NCOPY 0
# endif
                ! define radial functions for matching (->p1,p2,dp,pd1...)
                if(ifrq==1) then ; if(frq(ifrq)/=0) Error('!')
                  do n = 1,NCOPY
                    p1(:,n) = b1(:,1,l1) ; pd1(:,n) = b1(:,2,l1)
                    p2(:,n) = b2(:,1,l1) ; pd2(:,n) = b2(:,2,l1)
                    dp(n)   = db(1,l1)   ; dpd(n)   = db(2,l1)
                    ep(n)   = ebas(1,l1,itype,ispin)
                  enddo
                  do n = NCOPY+1,nindx(l,itype)
                    p1(:,n) = h1(:,n)
                    p2(:,n) = h2(:,n)
                    dp(n)   = dh(n)
                    call dirac_inhom(pd1(:,n),pd2(:,n),dpd(n),l1,ebas(n,l,itype,ispin),p1(:,n),p2(:,n),itype,ispin)
                    rdum     = intgrf(p1(:,n)*pd1(:,n)+p2(:,n)*pd2(:,n),itype)
                    pd1(:,n) = pd1(:,n) - rdum * p1(:,n)
                    pd2(:,n) = pd2(:,n) - rdum * p2(:,n)
                    dpd(n)   = dpd(n)   - rdum * dp(n)
                    ep(n)    = ebas(n,l,itype,ispin)
                  enddo
                endif

                do ll = abs(l-l1),min(l+l1,lcutm(itype)),2 ! Gaunt condition l+ll+l1=even
                  do nn = 1,nindxm(ll,itype)

                    lm = sum ( [ ((2*l2+1)*nindxm(l2,itype),l2=0,ll-1) ] ) + nn
                    ic = sum(neq(:itype-1))
                    do ieq = 1,neq(itype)
                      ic                    = ic + 1
                      screen1(:,ieq,-ll:ll) = transpose(screenk_mt(lm:lm+2*ll*nindxm(ll,itype):nindxm(ll,itype),:nlmp,ic,ifrq))
                    enddo

                    do n = 1,nindx(l,itype)
                      efrq = ebas(n,l,itype,ispin) + frq(ifrq)

                      if(n/=2) then
                        ! (A) RESPONSE FOR u AND ulo
                        ! right-hand side (->rhs1,rhs2)
                        rhs1 = -b1(:,n,l) * basm(:,nn,ll,itype) / rgrid(:,itype)
                        rhs2 = -b2(:,n,l) * basm(:,nn,ll,itype) / rgrid(:,itype)
                        if(l1==l) then
                          debas = intgrf( (b1(:,n,l)**2+b2(:,n,l)**2) * basm(:,nn,ll,itype) / rgrid(:,itype) ,itype)
                          rhs1  = rhs1 + b1(:,n,l) * debas
                          rhs2  = rhs2 + b2(:,n,l) * debas
                        endif
                        ! solve inhomogeneous differential equation (->ihom1,ihom2)
                        call cdirac_inhom(ihom1,ihom2,dihom,l1,efrq,rhs1,rhs2,itype,ispin)
                        ! match with uhom and u (or udot and u)
                        if(abs(efrq-ep(n))>1d-8) then
                          hom1 = h1(:,n)
                          hom2 = h2(:,n)
                          dhom = dh(n)
                          call continuous(c0,c1,ihom1(ng)/rg,dihom,hom1(ng)/rg,dhom,p1(ng,n)/rg*(1d0,0d0),dp(n)*(1d0,0d0))
                          ihom1 = ihom1 + c0 * hom1 + c1 * p1(:,n)
                          ihom2 = ihom2 + c0 * hom2 + c1 * p2(:,n)
                        else
                          call continuous(c0,c1,ihom1(ng)/rg,dihom,pd1(ng,n)/rg*(1d0,0d0),dpd(n)*(1d0,0d0),
     &                                                              p1(ng,n)/rg*(1d0,0d0), dp(n)*(1d0,0d0))
                          ihom1 = ihom1 + c0 * pd1(:,n) + c1 * p1(:,n)
                          ihom2 = ihom2 + c0 * pd2(:,n) + c1 * p2(:,n)
                        endif
                      else
                        ! (B) RESPONSE FOR udot (by finite difference)
                        rdum = delta
                        j    = 1
  1                     rhs1 = -b1_(:,j) * basm(:,nn,ll,itype) / rgrid(:,itype)
                        rhs2 = -b2_(:,j) * basm(:,nn,ll,itype) / rgrid(:,itype)
                        if(l1==l) then
                          debas = intgrf( (b1_(:,j)**2+b2_(:,j)**2) * basm(:,nn,ll,itype) / rgrid(:,itype) ,itype)
                          rhs1  = rhs1 + b1_(:,j) * debas
                          rhs2  = rhs2 + b2_(:,j) * debas
                        endif
                        call cdirac_inhom(ihom1,ihom2,dihom,l1,efrq+rdum,rhs1,rhs2,itype,ispin)
                        if(abs(efrq+rdum-ep(1))>1d-8) then
                          hom1 = h1_(:,j)
                          hom2 = h2_(:,j)
                          dhom = dh_(j)
                          call continuous(c0,c1,ihom1(ng)/rg,dihom,hom1(ng)/rg,dhom,p1(ng,1)/rg*(1d0,0d0),dp(1)*(1d0,0d0))
                          ihom1 = ihom1 + c0 * hom1 + c1 * p1(:,1)
                          ihom2 = ihom2 + c0 * hom2 + c1 * p2(:,1)
                        else
                          call continuous(c0,c1,ihom1(ng)/rg,dihom,pd1(ng,1)/rg*(1d0,0d0),dpd(1)*(1d0,0d0),
     &                                                              p1(ng,1)/rg*(1d0,0d0), dp(1)*(1d0,0d0))
                          ihom1 = ihom1 + c0 * pd1(:,1) + c1 * p1(:,1)
                          ihom2 = ihom2 + c0 * pd2(:,1) + c1 * p2(:,1)
                        endif
                        if(rdum>0) then
                          rdum   = -delta
                          j      = 2
                          ihom1_ = ihom1
                          ihom2_ = ihom2
                          goto 1
                        endif
                        ihom1 = (ihom1_-ihom1) / (2*delta)
                        ihom2 = (ihom2_-ihom2) / (2*delta)
                      endif
                      ! add negative frquency, make real (->tild1,tild2)
                      tild1 = 2*real(ihom1)
                      tild2 = 2*real(ihom2)

c                      tild1 = bas1(:,n,l,itype,ispin) * basm(:,nn,ll,itype) / rgrid(:,itype)
c                      tild2 = bas2(:,n,l,itype,ispin) * basm(:,nn,ll,itype) / rgrid(:,itype)

                      ! Calculate contraction   contract(lm0,lm) = SUM(IJ) <J lm0 | tilde(lmI) > W(IJ)
                      do l0 = 0,lcut(itype)
                        do ll0 = abs(l0-l1),min(l0+l1,lcutm(itype)),2
                          do n0 = 1,nindx(l0,itype)
                            help1 = bas1(:,n0,l0,itype,ispin) * tild1 / rgrid(:,itype)
                            help2 = bas2(:,n0,l0,itype,ispin) * tild2 / rgrid(:,itype)
                            do nn0 = 1,nindxm(ll0,itype)
                              int(nn0) = intgrf(basm(:,nn0,ll0,itype)*(help1+help2),itype)
                            enddo

                            do mm = -ll,ll
                              do mm0 = -ll0,ll0
                                j    = (ll0+mm0)*nindxm(ll0,itype) + sum( [ ((2*l2+1)*nindxm(l2,itype),l2=0,ll0-1) ] )
                                nn0  = nindxm(ll0,itype)
                                carr = matmul ( int(:nn0) , screen1(j+1:j+nn0,:,mm) )
                                do m1 = -l1,l1
                                  m    = m1 - mm  ; if(abs(m) >l ) cycle
                                  m0   = m1 - mm0 ; if(abs(m0)>l0) cycle
                                  lm0  = (l0+m0)*nindx(l0,itype) + n0 + sum( [ ((2*l2+1)*nindx(l2,itype),l2=0,l0-1) ] )
                                  lm   = (l +m )*nindx(l ,itype) + n  + sum( [ ((2*l2+1)*nindx(l2,itype),l2=0,l -1) ] )
                                  rdum = gaunt1(ll,l,mm,m) * gaunt1(ll0,l0,mm0,m0)
                                  contract(lm0,lm,:,ifrq) = contract(lm0,lm,:,ifrq) + carr * rdum
                                enddo
                              enddo
                            enddo

                          enddo
                        enddo
                      enddo

                    enddo
                  enddo
                enddo
              enddo
            enddo
          enddo

          deallocate ( carr )
          allocate ( carr(nlm) )

          iselfc = 0

          rdum = 0
          i = 0
          do l0 = 0,lcut(itype)
            do m0 = -l0,l0
              do n0 = 1,nindx(l0,itype)
                i = i + 1
                j = 0
                do l = 0,lcut(itype)
                  do m = -l,l
                    do n = 1,nindx(l,itype)
                      j = j + 1
c                      write(*,*) l0,m0,n0,l,m,n
c                      write(*,*) contract(i,j,:,1) * vol
                      if(l0==l.and.m0==m) then
c                        write(*,*) intgrf(bas1(:,n0,l,itype,ispin)*bas1(:,n,l,itype,ispin)+
c     &                    bas2(:,n0,l,itype,ispin)*bas2(:,n,l,itype,ispin),itype)
c                        read(*,*)
                        do ieq = 1,neq(itype)
                          rdum = rdum + abs ( contract(i,j,ieq,1) * vol -
     &                      intgrf(bas1(:,n0,l,itype,ispin)*bas1(:,n,l,itype,ispin)+
     &                      bas2(:,n0,l,itype,ispin)*bas2(:,n,l,itype,ispin),itype) )
                        enddo
                      else
                        do ieq = 1,neq(itype)
                          rdum = rdum + abs ( contract(i,j,ieq,1) * vol )
                        enddo
                      endif
                    enddo
                  enddo
                enddo
              enddo
            enddo
          enddo
          write(*,*) ' error',rdum,i,j,nlm
c          Error(' ')

c          ibc_selfc = 0
          do iblock = 1,nblock
            ib = block(1,iblock)
            nb = sizeblock(iblock)

            if(job1%spin(ib)==ispin) then
              allocate ( band(nb) )
              ikptq = job1%kpt(ib)
              band  = job1%band(block(:nb,iblock))

              ic = sum(neq(:itype-1))
              do ieq = 1,neq(itype)
                ic = ic + 1 ; write(*,*) ic
                do ifrq = 1,nfrq
                  i = iselfc
                  do jb = 1,nb
                    carr = matmul ( contract(:,:,ieq,ifrq) , cmt(:nlm,ic,band(jb),ikptq,ispin) )
                    do ib = 1,nb
                      i                 = i + 1
                      ibc_selfc(i,ifrq) = ibc_selfc(i,ifrq) + dot_product( cmt(:nlm,ic,band(ib),ikptq,ispin) , carr )
                    enddo
                  enddo
                enddo
              enddo
            endif

            ! integrate
            i = iselfc
            do jb = 1,nb
              do ib = 1,jb
                i = i + 1
                write(*,*) ib,jb,freqintegral0(real(ibc_selfc(i,:)),imag(frq),nfrq,4,6)
              enddo
            enddo

            iselfc = iselfc + nb**2
            deallocate ( band )
          enddo

          deallocate ( carr )
          deallocate ( contract )
          deallocate ( screen1 )

        enddo
      enddo

      write(*,'(3F20.10)') (imag(frq(ifrq)),ibc_selfc(1,ifrq)-ibc_selfc1(1,ifrq),ifrq=1,nfrq)
      write(801,'(3F20.10)') (imag(frq(ifrq)),ibc_selfc(1,ifrq)-ibc_selfc1(1,ifrq),ifrq=1,nfrq)
      write(0,'(F20.10)') real(ibc_selfc(1,:4))
C      write(0,'(F20.10)') ibc_selfc(1,2)

      call cpu_done(time1)

      contains

c     Returns coefficients a and b so that uu + a*u + b*udot is continuous at the MT sphere boundary (uu and duu are boundary values)
      subroutine continuous(c1,c2,uu,duu,u1,du1,u2,du2)
      implicit none
      complex_dp, intent(out) :: c1,c2
      complex_dp, intent(in)  :: uu,duu,u1,du1,u2,du2
      complex_dp              :: wronsk
      wronsk = u1 * du2 - u2 * du1
      c1     = - ( uu * du2 - duu * u2 ) / wronsk
      c2     =   ( uu * du1 - duu * u1 ) / wronsk
      end subroutine continuous

      end
# undef NCOPY
# undef GAUGE

c -----------------------

c     (B) Double-counting contribution
c
c     - SUM(m) SUM(IJ) INT dw <J n|m> <m|tilde(n'I,w)> W(IJ,w)
c
c     (a) Subroutine ibc2_doublecount1 contracts over m:                   cprod_ibc    = <J n|a> = SUM(m) <J n |m> <m|a>
c     (b) Subroutine ibc2_doublecount2 constructs the frequency integrand: ibc_selfc(w) = SUM(I) { SUM(a) [ SUM(J) W(IJ,w) <J n|a> ] <a|tild(n'I,w)> }
c
c     (a)
      subroutine ibc2_doublecount1(cprod_ibc,cprod,n2,cmt2)
      use global
      use wrapper
      complex_dp,  intent(inout) :: cprod_ibc(nbasp,maxlmindx)
      MCOMPLEX_dp, intent(in)    :: cprod(maxbasm,maxband)
      complex_dp,  intent(in)    :: cmt2(n2,maxlmindx,ncent)
      complex_dp                 :: cprod1(nbasp,n2)
      integer                    :: iband
      integer                    :: i0,ic,itype,ieq,nlmp,nlm,l
# ifdef INV
      do iband = 1,n2
        cprod1(:,iband) = cprod(:nbasp,iband)
        call desymmetrize(cprod1(1,iband),nbasp,1,1)
      enddo
# else
      cprod1 = cprod(:nbasp,:n2)
# endif
      i0 = 0
      ic = 0
      do itype = 1,ntype
        do ieq = 1,neq(itype)
          ic   = ic + 1
          nlmp = sum ( [ ((2*l+1)*nindxm(l,itype),l=0,lcutm(itype)) ] )
          nlm  = sum ( [ ((2*l+1)*nindx (l,itype),l=0,lcut (itype)) ] )
          cprod_ibc(i0+1:i0+nlmp,:nlm) = cprod_ibc(i0+1:i0+nlmp,:nlm) +
     &                                   matmat ( cprod1(i0+1:i0+nlmp,:n2), conjg(cmt2(:n2,:nlm,ic)) )
          i0   = i0 + nlmp
        enddo
      enddo
      end

c -----------------------
c
c     (b)
c
c     The factors nkpts1 and nsym1 are needed for adding the equivalent k points (see selfenergy)
      subroutine ibc2_doublecount2(self,band,nb,nfrq,screen,dim,ctrafo,cprod_ibc,ikpt1,nkpts1,nsym1,ispin)
      use global
      use wrapper, only: matmat
      use, intrinsic :: iso_fortran_env
      implicit none
      integer,     intent(in)    :: nb,nfrq,dim,band(nb),ikpt1,nkpts1,nsym1,ispin
      complex_dp,  intent(out)   :: self(nb,nb,nfrq)
      complex_dp,  intent(in)    :: cprod_ibc(nbasp,maxlmindx,nb)
      MCOMPLEX_dp, intent(in)    :: screen(dim*(dim+1)/2,nfrq),ctrafo(nbasm(ikpt1),dim)
      MCOMPLEX_dp                :: screen0(dim*(dim+1)/2)
      complex_dp,  allocatable   :: screen1(:,:),cprod1(:,:),cprod2(:,:,:),ctrafo1(:,:)
      complex_dp                 :: cmt1(nb,maxindx),proj(maxindx,maxindxm),cexp
      real_dp                    :: gnt
      integer                    :: ifrq,itype,ieq,ic,i0,i,j,nwblock,wblock
      integer                    :: nlmp,nlm,ln,ln1,lln,ib,jb
      integer                    :: l,m,n,l1,m1,n1,ll,mm,nn,l2,lm,lm1,lm1_0,llm
      real_dp                    :: gaunt!,intgrf
      if(symkpt(ikpt1)/=1) Error('!')

c      write(*,*) 'ibc_proj redefined'
      if(.not.allocated(ibc_proj)) Bug('ibc_proj not allocated.')
# if 0
      do itype = 1,ntype
        ln = 0
        do l = 0,lcut(itype)
          do n = 1,nindx(l,itype)
            ln  = ln + 1
            ln1 = 0
            do l1 = 0,lcut(itype)
              do n1 = 1,nindx(l1,itype)
                ln1 = ln1 + 1
                lln = 0
                do ll = 0,lcutm(itype)
                  do nn = 1,nindxm(ll,itype)
                    lln = lln + 1
                    ibc_proj(ln1,ln,lln,itype,:,ispin) = intgrf((
     &                bas1(:,n1,l1,itype,ispin)*bas1(:,n,l,itype,ispin)+
     &                bas2(:,n1,l1,itype,ispin)*bas2(:,n,l,itype,ispin))*basm(:,nn,ll,itype)/rgrid(:,itype),itype)
                  enddo
                enddo
              enddo
            enddo
          enddo
        enddo
      enddo
# endif

# ifdef INV
      allocate ( ctrafo1(nbasp,dim) )
      ctrafo1 = ctrafo(:nbasp,:)
      call desymmetrize(ctrafo1,nbasp,dim,1)
# endif
      do ifrq = 1,nfrq
        i0 = 0
        ic = 0

        ! change from block-diagonal to normal form (if cblock is defined) (->screen0)
        if(allocated(cblock)) then
          nwblock = maxval(cblock)
          screen0 = 0
          n       = 0
          m       = 0
          do ib = 1,nwblock
            wblock = count(cblock==ib)
            do j = 1,wblock ; do i = 1,j
              n           = n + 1
              nn          = (m+j)*(m+j-1)/2 + m + i
              screen0(nn) = screen(n,ifrq)
            enddo ; enddo
            m = m + wblock
          enddo
        else
          screen0 = screen(:,ifrq)
        endif

        do itype = 1,ntype
          nlmp = sum( [ ((2*l2+1)*nindxm(l2,itype),l2=0,lcutm(itype)) ] )
          nlm  = sum( [ ((2*l2+1)*nindx (l2,itype),l2=0,lcut (itype)) ] )
          allocate ( screen1(nlmp,nlmp) )
          allocate ( cprod1(nlmp,nlm) )
          allocate ( cprod2(nlm,maxindxm,nb) )
          do ieq = 1,neq(itype)
            ic      = ic + 1
            cexp    = exp ( img * 2*pi * dot_product(kpt(:,ikpt1),cent(:,ic)) )
# ifdef INV
            screen1 =        matmat ( ctrafo1(i0+1:i0+nlmp,:) , matmat ( screen0,
     &                transpose(conjg(ctrafo1(i0+1:i0+nlmp,:))) ) )
# else
            screen1 =        matmat ( ctrafo(i0+1:i0+nlmp,:) , matmat ( screen0,
     &                transpose(conjg(ctrafo(i0+1:i0+nlmp,:))) ) )
# endif

            do ib = 1,nb
              cprod1 = matmul(screen1,cprod_ibc(i0+1:i0+nlmp,:nlm,ib))

              llm = 0
              do ll = 0,lcutm(itype)
                lln = sum ( [ (nindxm(l2,itype),l2=0,ll-1) ] )
                do mm = -ll,ll

                  cprod2 = 0

                  lm = 0
                  do l = 0,lcut(itype)
                    ln = sum ( [ (nindx(l2,itype),l2=0,l-1) ] )
                    do m = -l,l

                      n          = nindx(l,itype)
                      cmt1(:,:n) = cmt(lm+1:lm+n,ic,band,ikpt1,ispin) * cexp

                      m1 = mm + m
                      do l1 = abs(l-ll),min(l+ll,lcut(itype)),2
                        if(abs(m1)>l1) cycle
                        ln1   = sum ( [ (         nindx(l2,itype),l2=0,l1-1) ] )
                        lm1_0 = sum ( [ ((2*l2+1)*nindx(l2,itype),l2=0,l1-1) ] ) + (m1+l1)*nindx(l1,itype)
                        gnt   = gaunt(ll,l1,l,mm,m1,m)

                        do jb = 1,nb
                          proj = 0
                          n1   = nindx(l1,itype)
                          nn   = nindxm(ll,itype)
                          do n = 1,nindx(l,itype)
                            proj(:n1,:nn) = proj(:n1,:nn) +
     &                                      ibc_proj(ln1+1:ln1+n1,ln+n,lln+1:lln+nn,itype,ifrq,ispin) * gnt * cmt1(jb,n)
                          enddo
                          do nn = 1,nindxm(ll,itype)
                            lm1 = lm1_0
                            do n1 = 1,nindx(l1,itype)
                              lm1               = lm1 + 1
                              cprod2(lm1,nn,jb) = cprod2(lm1,nn,jb) + proj(n1,nn)
                            enddo
                          enddo
                        enddo

                      enddo
                      lm = lm + nindx(l,itype)
                    enddo
                  enddo
                  if(lm/=nlm) Bug('Count error.')

                  do nn = 1,nindxm(ll,itype)
                    llm = llm + 1
                    do jb = 1,nb
                      self(ib,jb,ifrq) = self(ib,jb,ifrq) - sum ( cprod1(llm,:)*cprod2(:,nn,jb) ) * nkpts1 / nsym1
                    enddo
                  enddo

                enddo
              enddo

              if(llm/=nlmp) Bug('Count error.')

            enddo

            i0 = i0 + nlmp
          enddo
          deallocate ( screen1,cprod1,cprod2 )
        enddo
      enddo
# ifdef INV
      deallocate ( ctrafo1 )
# endif
      if(i0/=nbasp) Bug('Count error.')
      end



