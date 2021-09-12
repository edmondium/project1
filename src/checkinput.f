c--------------------------------------------------------------------------------
c Copyright (c) 2021 Peter Grünberg Institut, Forschungszentrum Jülich, Germany
c This file is part of SPEX and available as free software under the conditions
c of the MIT license as expressed in the LICENSE file in more detail.
c--------------------------------------------------------------------------------

# include "cppmacro.h"

      subroutine checkinput
      use global
      use util, only: chr
      use wrapper
      use file
      use key
      Load( use readwrite )
c      use m_radsra

      use, intrinsic :: iso_fortran_env
      implicit none
      integer                 :: ispin,itype,l,i,j,m,n,lm,icent,ikpt,ineq,lm1,igpt,iband,nk,iarr(2),ikindx
      real_dp,    allocatable :: olapcb(:),olapbb(:,:,:,:,:)
      real_dp,    allocatable :: olapcv_avg(:,:,:,:,:),olapcv_max(:,:,:,:,:)
      MCOMPLEX_dp,allocatable :: olappp(:,:)
      complex_dp, allocatable :: olapcv(:,:),olapww(:,:),carr(:),carr2(:,:)
      integer,    allocatable :: olapcv_loc(:,:,:,:,:,:)
      real_dp                 :: olapww_avg,olapww_max,rdum1
      real_dp                 :: sphbes(0:maxlcut),q(3),qnorm,rdum,rarr(maxband)
c      real_dp                 :: dcore1(maxgrid)
      real_dp                 :: e,emin,emax,de
      real_dp                 :: u1(maxgrid),u2(maxgrid),du,olap(maxindx,maxindx),vec(maxindx),evec(maxindx,maxindx)
      complex_dp              :: y((maxlcut+1)**2),cdum,cexp
      integer                 :: iunit
      logical                 :: ldum
      real_dp,    allocatable :: mtaccur(:)
      real                    :: time1,time2
      real_dp                 :: intgr,intgrf

      beginSingle

      write(6,'(//A)') '### subroutine: checkinput ###'

      write(6,'(/A)') 'Overlap <core|core>'
      do ispin = 1,nspin1
        if(nspin==2.and.ispin==1) write(6,'(A)') 'Spin up'
        if(nspin==2.and.ispin==2) write(6,'(A)') 'Spin down'
        do itype = 1,ntype
          if(all(nindxc(:lcutc(itype),itype)==0)) cycle
          if(ntype>1) write(6,'(A,I3)') 'Atom type',itype
          do l = 0,lcutc(itype)
            do i = 1,nindxc(l,itype)
              if(lcore_soc) then
                j = 2*l + 1 ; if(l/=0.and.mod(i,2)==0) j = 2*l-1
                n = i       ; if(l/=0)                 n = (i-1)/2 + l + 1
                write(6,'(I1,A,I1,A,1X'NoA) n,lchar(l),j,'/2'
              else
                write(6,'(I1,A,2X'NoA) i+l,lchar(l)
              endif
              do j = 1,i
                if(lcore_soc.and.l/=0.and.mod(i,2)/=mod(j,2)) cycle
                write(6,'(F10.6'NoA) intgr(
     &            core1(:,i,l,itype,ispin)*core1(:,j,l,itype,ispin)+
     &            core2(:,i,l,itype,ispin)*core2(:,j,l,itype,ispin),
     &            grid(itype))
              enddo
              write(6,*)
            enddo
          enddo
        enddo
      enddo

c      write(6,'(/A)') 'Wronskian (core/valence)'
c      do ispin=1,nspin
c        if(nspin==2.and.ispin==1) write(6,'(A)') 'Spin up'
c        if(nspin==2.and.ispin==2) write(6,'(A)') 'Spin down'
c        do itype=1,ntype
c          if(ntype>1) write(6,'(A,I3)') 'Atom type',itype
c          do l=0,lcutc(itype)
c            if(l>lcut(itype)) exit ! very improbable case
c            write(6,"(9X,'u(',A,')',4X,'udot(',A,')',:,3X,'ulo(',A,') ...')") (lchar(l),i=1,min(3,nindx(l,itype)))
c            do i=1,nindxc(l,itype)
c              call derivative(dcore1,core1(:,i,l,itype,ispin),itype)
c              dcore1 = dcore1 - core1(:,i,l,itype,ispin) / rgrid(:,itype)
c              write(6,'(I1,A,2X'NoA) i+l,lchar(l)
c              do j=1,nindx(l,itype)
c                write(6,'(F10.6'NoA) grid(itype)%radius / 2 *
c     &            ( ubas(j,l,itype,ispin) * dcore1(grid(itype)%number) -
c     &             dubas(j,l,itype,ispin) *  core1(grid(itype)%number,i,l,itype,ispin) )
c              enddo
c              write(6,*)
c            enddo
c          enddo
c        enddo
c      enddo
c      write(123,'(4F20.10)') (rgrid(i,1),core1(i,1,0,1,1),core1(i,2,0,1,1),core1(i,1,1,1,1),i=1,grid(1)%number)
c      Error(' ')

      write(6,'(/A)') 'Overlap <core|basis>'
      nk = size(cpw,3)
      allocate ( olapcb(maxindx),olapcv(maxband,nk),
     &           olapcv_avg(  -maxlcutc:maxlcutc,maxindxc,0:maxlcutc,ntype,nspin1),
     &           olapcv_max(  -maxlcutc:maxlcutc,maxindxc,0:maxlcutc,ntype,nspin1),
     &           olapcv_loc(2,-maxlcutc:maxlcutc,maxindxc,0:maxlcutc,ntype,nspin1) )
      do ispin = 1,nspin1
        if(nspin1==2.and.ispin==1) write(6,'(A)') 'Spin up'
        if(nspin1==2.and.ispin==2) write(6,'(A)') 'Spin down'
        do itype = 1,ntype
          if(all(nindxc(:lcutc(itype),itype)==0)) cycle
          if(ntype>1) write(6,'(A,I3)') 'Atom type',itype
          do l = 0,lcutc(itype)
            if(l>lcut(itype)) exit ! very improbable case
            write(6,"(9X,'u(',A,')',4X,'udot(',A,')',:,3X,'ulo(',A,') ...')") (lchar(l),i=1,min(3,nindx(l,itype)))
            do i = 1,nindxc(l,itype)
              if(lcore_soc) then
                j = 2*l + 1 ; if(l/=0.and.mod(i,2)==0) j = 2*l-1
                n = i       ; if(l/=0)                 n = (i-1)/2 + l + 1
                write(6,'(I1,A,I1,A,1X'NoA) n,lchar(l),j,'/2'
              else
                write(6,'(I1,A,2X'NoA) i+l,lchar(l)
              endif
              do j = 1,nindx(l,itype)
                olapcb(j) = intgr(
     &            core1(:,i,l,itype,ispin)*bas1(:,j,l,itype,ispin)+
     &            core2(:,i,l,itype,ispin)*bas2(:,j,l,itype,ispin),
     &            grid(itype) )
                write(6,'(F10.6'NoA) olapcb(j)
              enddo
# ifndef LOAD
              lm    = sum ( [ (nindx(j,itype)*(2*j+1),j=0,l-1) ] )
              icent = sum(neq(1:itype-1))+1 ! take first of group of equivalent atoms
              do m = -l,l
                olapcv = 0
                do j = 1,nindx(l,itype)
                  lm = lm + 1
                  olapcv = olapcv + olapcb(j)*cmt(lm,icent,:,:nk,ispin)
                  if(l_soc) olapcv = olapcv + olapcb(j)*cmt(lm,icent,:,:nk,2)
                enddo
                rdum                            = sum    ( abs(olapcv(:,:))**2 )
                rdum1                           = maxval ( abs(olapcv(:,:))    )
                iarr                            = maxloc ( abs(olapcv(:,:))    )
                olapcv_avg(  m,i,l,itype,ispin) = sqrt( rdum / nk / sum(nband) * size(nband,1) )
                olapcv_max(  m,i,l,itype,ispin) = rdum1
                olapcv_loc(:,m,i,l,itype,ispin) = iarr
              enddo
# endif
              write(6,*)

            enddo
          enddo
        enddo
      enddo

      if(lcore_soc) then
        Info('<core|val> overlap not implemented for CORESOC. Skipping...')
        goto 2
      endif

# ifndef LOAD

      write(6,'(/A)') 'Average overlap <core|val>'
      do ispin=1,nspin1
        if(nspin1==2.and.ispin==1) write(6,'(A)') 'Spin up'
        if(nspin1==2.and.ispin==2) write(6,'(A)') 'Spin down'
        do itype=1,ntype
          if(all(nindxc(:lcutc(itype),itype)==0)) cycle
          if(ntype>1) write(6,'(A,I3)') 'Atom type',itype
          do l=0,lcutc(itype)
            do i=1,nindxc(l,itype)
              write(6,'(I1,A,2X'NoA) i+l,lchar(l)
              write(6,'('//chr(2*l+1)//'F10.6)') olapcv_avg(-l:l,i,l,itype,ispin)
            enddo
          enddo
        enddo
      enddo

      write(6,'(/A)') 'Maximum overlap <core|val> at (band/kpoint)'
      do ispin=1,nspin1
        if(nspin==2.and.ispin==1) write(6,'(A)') 'Spin up'
        if(nspin==2.and.ispin==2) write(6,'(A)') 'Spin down'
        do itype=1,ntype
          if(all(nindxc(:lcutc(itype),itype)==0)) cycle
          if(ntype>1) write(6,'(A,I3)') 'Atom type',itype
          do l=0,lcutc(itype)
            do i=1,nindxc(l,itype)
              write(6,'(I1,A,2X'NoA) i+l,lchar(l)
              write(6,'('//chr(2*l+1)//'(F10.6,'' ('',A,''/'',A,'')''))')
     &                  (olapcv_max(  m,i,l,itype,ispin),
     &          trim(chr(olapcv_loc(1,m,i,l,itype,ispin),'I3.3')),
     &          trim(chr(olapcv_loc(2,m,i,l,itype,ispin),'I4.3')),m=-l,l)
            enddo
          enddo
        enddo
      enddo

# else
      write(6,'(/A)') 'LOAD: Average and maximum overlap not calculated.'
# endif

 2    deallocate (olapcb,olapcv,olapcv_avg,olapcv_max,olapcv_loc)

      write(6,'(/A)') 'Overlap <basis|basis>'
      allocate ( olapbb(maxindx,maxindx,0:maxlcut,ntype,nspin),carr(max(maxindx,maxgpt)) )
      do ispin=1,nspin
        if(nspin==2.and.ispin==1) write(6,'(A)') 'Spin up'
        if(nspin==2.and.ispin==2) write(6,'(A)') 'Spin down'
        do itype=1,ntype
          if(ntype>1) write(6,'(A,I3)') 'Atom type',itype
          do l=0,lcut(itype)
            do i=1,nindx(l,itype)
              select case(i)
                case(1);      write(6,'(''   u('',A,'')'''NoA) lchar(l)
                case(2);      write(6,'(''udot('',A,'')'''NoA) lchar(l)
                case default; write(6,'('' ulo('',A,'')'''NoA) lchar(l)
              end select
              do j=1,i
                olapbb(i,j,l,itype,ispin) = intgr(
     &            bas1(:,i,l,itype,ispin)*bas1(:,j,l,itype,ispin)+
     &            bas2(:,i,l,itype,ispin)*bas2(:,j,l,itype,ispin),
     &            grid(itype) )
                write(6,'(F10.6'NoA) olapbb(i,j,l,itype,ispin)
                olapbb(j,i,l,itype,ispin) = olapbb(i,j,l,itype,ispin)
              enddo
              write(6,*)
            enddo
          enddo
        enddo
      enddo

      call getkey(inp,'CHKOLAP', ldum, default=.false.)
      if(ldum) then
        olapww_avg = 0
        olapww_max = 0
        write(6,'(/A)') 'Timing:'
        write(6,'(A)')  '#kpt   pw     mt'
        allocate ( olappp(maxgpt,maxgpt) )
        do ispin = 1,nspin1
          if(lkptadd) then ; nk = nkpt * 2
          else             ; nk = nkpt
          endif
          do ikpt = 1,nk
            if(storeibz.and.kptp(ikpt)/=ikpt) cycle
# ifdef LOAD
            Allocate_ ( cmt,(maxlmindx,ncent,nband(ikpt,ispin),ikpt:ikpt,nspin3) )
            Allocate_ ( cpw,(maxgpt,         nband(ikpt,ispin),ikpt:ikpt,nspin3) )
            call read_wavef2([(i,i=1,nband(ikpt,ispin))],nband(ikpt,ispin),ikpt,ispin,cmt,cpw)
            kindx(ikpt) = ikpt
# endif
            call olap_gpt(olappp,maxgpt,ikpt,ikpt)
            ikindx = kindx(ikpt)
            allocate ( olapww(nband(ikpt,ispin),nband(ikpt,ispin)) )
            write(6,'(I4'NoA) ikpt
            call cpu_time(time2)
            n = ngpt(ikpt)
            do i=1,nband(ikpt,ispin)
              carr(1:n) = matmul ( olappp(1:n,1:n),cpw(1:n,i,ikindx,ispin) )
              do j=1,i
                olapww(j,i) = dot_product ( cpw(1:n,j,ikindx,ispin),carr(1:n) )
              enddo
              if(l_soc) then
                carr(1:n) = matmul ( olappp(1:n,1:n),cpw(1:n,i,ikindx,2) )
                do j=1,i
                  olapww(j,i) = olapww(j,i) + dot_product ( cpw(1:n,j,ikindx,2),carr(1:n) )
                enddo
              endif
            enddo
            call cpu_time(time1); write(6,'(F7.2'NoA) time1-time2
            icent  = 0
            do itype = 1,ntype
              do ineq = 1,neq(itype)
                icent = icent+1
                lm    = 0
                do l=0,lcut(itype)
                  n = nindx(l,itype)
                  do m=1,2*l+1
                    do i=1,nband(ikpt,ispin)
                      carr(1:n) = matmul ( olapbb(1:n,1:n,l,itype,ispin),
     &                                     cmt(lm+1:lm+n,icent,i,ikindx,ispin) )
                      do j=1,i
                        olapww(j,i) = olapww(j,i) + dot_product ( cmt(lm+1:lm+n,icent,j,ikindx,ispin),
     &                                                            carr(1:n) )
                      enddo
                      if(l_soc) then
                        carr(1:n) = matmul ( olapbb(1:n,1:n,l,itype,nspin),
     &                                       cmt(lm+1:lm+n,icent,i,ikindx,2) )
                        do j=1,i
                          olapww(j,i) = olapww(j,i) + dot_product ( cmt(lm+1:lm+n,icent,j,ikindx,2),
     &                                                              carr(1:n) )
                        enddo
                      endif
                    enddo
                    lm = lm + n
                  enddo
                enddo
              enddo
            enddo

            call cpu_time(time2); write(6,'(F7.2'NoA) time2-time1

            do i=1,nband(ikpt,ispin)
              olapww(i,i) = 1 - olapww(i,i)
              do j=1,i
                olapww(i,j) = conjg(olapww(j,i))
              enddo
            enddo
            rdum       = sqrt( sum ( abs(olapww(:,:))**2 ) ) / nband(ikpt,ispin) / nkpt
            olapww_max = max ( olapww_max,maxval(abs(olapww(:,:))) )
            olapww_avg = olapww_avg + rdum
            write(6,'(F24.20)') rdum
            deallocate ( olapww )

            Load( Deallocate_(cmt) )
            Load( Deallocate_(cpw) )

          enddo

        enddo
        deallocate (olappp)

        write(6,'(/A)')       'Overlap <val|val>'
        write(6,'(A,F24.20)') 'Average deviation from unity matrix:',olapww_avg
        write(6,'(A,F24.20)') 'Maximum deviation from unity matrix:',olapww_max
      endif

 1    deallocate (olapbb,carr)


      call getkey(inp,'CHKMISM', ldum, default=.false.)

      if(ldum) then
        Load ( write(6,'(/A)') 'LOAD: CHKMISM disabled.' ; goto 3 )

        write(6,'(/A)') 'Mismatch of wave functions at the MT-sphere boundaries'
        allocate (carr2(maxband,(maxlcut+1)**2))
        do ispin = 1,nspin2
          if(ispin==1) write(6,'(2X,A)') 'Spin up'
          if(ispin==2) write(6,'(2X,A)') 'Spin down'
          icent = 0
          do itype = 1,ntype
            do ineq = 1,neq(itype)
              icent = icent + 1
              if(ncent>1) write(6,'(2X,A,I3)') 'Atom',icent
              write(6,'(2X,A)') 'k point    average    (   maximum    )'

              if(lkptadd) then ; nk = nkpt * 2
              else             ; nk = nkpt
              endif
              do ikpt = 1,nk
                if(storeibz.and.kptp(ikpt)/=ikpt) cycle
                ikindx = kindx(ikpt)
                write(6,'(I6,2X'NoA) ikpt
                carr2 = 0

                ! PW
                do igpt = 1,ngpt(ikpt)
                  cexp  = exp(img * 2*pi * dot_product(kpt(:,ikpt)+gpt(:,pgpt(igpt,ikpt)),cent(:,icent)) )
                  q     = matmul(rlat,kpt(:,ikpt)+gpt(:,pgpt(igpt,ikpt)))
                  qnorm = sqrt(sum(q**2))
                  call sphbessel(sphbes,grid(itype)%radius*qnorm,lcut(itype))
                  call harmonicsr(y,q,lcut(itype))
                  y  = conjg(y)
                  lm = 0
                  do l = 0,lcut(itype)
                    cdum = 4*pi*img**l/sqrt(vol) * sphbes(l) * cexp
                    do m = -l,l
                      lm = lm + 1
                      do iband = 1,nband(ikpt,min(ispin,nspin1))
                        carr2(iband,lm) = carr2(iband,lm) + cdum * cpw(igpt,iband,ikindx,ispin) * y(lm)
                      enddo
                    enddo
                  enddo
                enddo
                ! MT
                lm  = 0
                lm1 = 0
                do l = 0,lcut(itype)
                  do m = -l,l
                    lm = lm + 1
                    do n = 1,nindx(l,itype)
                      lm1  = lm1 + 1
                      rdum = bas1(grid(itype)%number,n,l,itype,ispin) / grid(itype)%radius
                      do iband = 1,nband(ikpt,min(ispin,nspin1))
                        carr2(iband,lm) = carr2(iband,lm) - cmt(lm1,icent,iband,ikindx,ispin) * rdum
                      enddo
                    enddo
                  enddo
                enddo

                rarr = 0
                lm   = 0
                do l = 0,lcut(itype)
                  do m = -l,l
                    lm   = lm + 1
                    rarr = rarr + abs(carr2(:,lm))**2
                  enddo
                enddo
                rarr = sqrt ( rarr / (4*pi) )
                write(6,'(2X,F14.12,''  ('',F14.12,'')'')') sum(rarr**2/nband(ikpt,min(ispin,nspin1))),maxval(rarr)

              enddo
            enddo
          enddo

        enddo

        deallocate (carr2)

      endif

 3    allocate ( mtaccur(2) )
      call getkey(inp,'MTACCUR', mtaccur, status=i, allow_eV=.true.)
      if(i/=0) then
        if(i==1) mtaccur = [ minval(ene(:maxeband,:nkpt,:)) - 1 , maxval(ene(:maxeband,:nkpt,:) + 1) ]
        if(mtaccur(1)>mtaccur(2)) Error('First argument to MTACCUR must be smaller than second.')
        write(6,'(/,A'NoA) 'The MT representation error is written to files spex.mt....'
        call cpu_time(time1)
        emin = mtaccur(1)
        emax = mtaccur(2)
        do ispin = 1,nspin
          do itype = 1,ntype
            if(nspin==2) then ; iunit = fopen('spex.mt.'//trim(chr(ispin))//'.'//trim(chr(itype)),status='unknown')
            else              ; iunit = fopen('spex.mt.'//trim(chr(itype)),                       status='unknown')
            endif
            do l = 0,lcut(itype)
              n     = nindx(l,itype)
              write(iunit,'(A)')              '# MT accuracy'
              write(iunit,'(A,I3,A,I3,A,I3)') '# Atom',itype,',  l =',l,',  Spin =',ispin

# if 0
              do i = 1,n
                if(i/=2) then ; call dirac_hom(u1,u2,du,l,ebas(i,l,itype,ispin),itype,ispin)
                else          ; call dirac_inhom(u1,u2,du,l,ebas(i,l,itype,ispin),
     &                               -bas1(:,1,l,itype,ispin),-bas2(:,1,l,itype,ispin),itype,ispin)
                  rdum = intgrf(u1*bas1(:,1,l,itype,ispin)+u2*bas2(:,1,l,itype,ispin),itype)
                  u1 = u1 - rdum * bas1(:,1,l,itype,ispin)
                  u2 = u2 - rdum * bas2(:,1,l,itype,ispin)
                  rdum = sqrt(intgrf(u1**2+u2**2,itype))
                  u1 = u1 / rdum
                  u2 = u2 / rdum
                endif
                bas1(:,i,l,itype,ispin) = u1
                bas2(:,i,l,itype,ispin) = u2
              enddo
# endif

              do i = 1,n
                do j = 1,n
                  olap(i,j) = intgrf( bas1(:,i,l,itype,ispin) * bas1(:,j,l,itype,ispin)
     &                              + bas2(:,i,l,itype,ispin) * bas2(:,j,l,itype,ispin) , itype )
c                  if(l==0) write(*,*) i,j,olap(i,j)
                enddo
              enddo
              call diagonalize(evec(:n,:n),vec(:n),olap(:n,:n))
              do i = 1,n
                if(vec(i)>1d-8) then
                  evec(:,i) = evec(:,i) / sqrt(vec(i))
                else
                  evec(:,i) = 0
                  Warn('MT basis over-complete, eigenvalue ='//chr(vec(i),'ES9.1'))
                  if(nspin==2) then ; write(0,'(A,3I3)')   '            Spin, atom type, l:',ispin,itype,l
                  else              ; write(0,'(A,2I3)')   '            Atom type, l:',itype,l
                  endif
                endif
              enddo
              olap = matmul ( evec(:n,:n), transpose ( evec(:n,:n) ) )
              de   = ( emax - emin ) / 1000
              e    = emin
c              qnorm = 0
              do i = 0,1000
                call dirac_hom(u1,u2,du,l,e,itype,ispin)
c                call radsra(e,l,vmt(:,1,itype,ispin)*rgrid(:,itype),grid(itype)%first,grid(itype)%increment,grid(itype)%number,
c     &            maxgrid,137.0359895d0,du,du,j,u1,u2)
# if 0
                rdum1 = sqrt(intgrf(u1**2+u2**2,itype))
                u1 = u1 / rdum1
                u2 = u2 / rdum1
                qnorm = qnorm + intgrf((u1-v1)**2+(u2-v2)**2,itype)
                do j = 1,n
                  vec(j) = intgrf ( bas1(:,j,l,itype,ispin) * v1 + bas2(:,j,l,itype,ispin) * v2 , itype )
                enddo
                vec = matmul(olap(:n,:n),vec(:n))
                do j = 1,n
                  v1 = v1 - vec(j) * bas1(:,j,l,itype,ispin)
                  v2 = v2 - vec(j) * bas2(:,j,l,itype,ispin)
                enddo
                if(itype==1.and.l==0) write(100,'(2F22.15)') e,intgrf(v1**2+v2**2,itype)
# endif
                do j = 1,n
                  vec(j) = intgrf ( bas1(:,j,l,itype,ispin) * u1 + bas2(:,j,l,itype,ispin) * u2 , itype )
                enddo
                vec = matmul ( olap(:n,:n),vec(:n) )
                do j = 1,n
                  u1 = u1 - vec(j) * bas1(:,j,l,itype,ispin)
                  u2 = u2 - vec(j) * bas2(:,j,l,itype,ispin)
                enddo
c                if(itype==2.and.l==0) write(101,'(2F22.15)') e,intgrf(u1**2+u2**2,itype)
                write(iunit,'(2F22.15)') e,intgrf(u1**2+u2**2,itype)
                e = e + de
              enddo
c              write(*,*) qnorm
              write(iunit,*)
              write(iunit,*)
            enddo
            call fclose(iunit)
          enddo
        enddo
        call cpu_done(time1)
      endif
      deallocate ( mtaccur )

      endSingle

      end

