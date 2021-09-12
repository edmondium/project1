c--------------------------------------------------------------------------------
c Copyright (c) 2021 Peter Grünberg Institut, Forschungszentrum Jülich, Germany
c This file is part of SPEX and available as free software under the conditions
c of the MIT license as expressed in the LICENSE file in more detail.
c--------------------------------------------------------------------------------

# include "cppmacro.h"

c     Calculates the derivative
c                                  ikr    s       s'
c     dcprod(n',n,k,xyz) = d'   < e    phi   | phi       > / sqrt(vol)
c                           xyz           qn      q+k,n'
c
c                                  s             s'
c                             < phi  | d    | phi    >
c                                  nq   xyz      n'q
c                        = -i ------------------------ / sqrt(vol)  ,   s = ispin1   ,   s' = ispin2
c                                    s'     s                       ,   n = occ.     ,   n' = unocc.
c                                   e    - e                        ,   bandi1 <= n <= bandf1 , bandi2 <= n' <= bandf2
c                                    n'q    nq
c
c     with kp perturbation theory and
c
c             d    d    d               d     d     d
c     d    = -- , -- , --  ;   d'   =  ---,  ---,  ---  .
c      xyz   dx   dy   dz       xyz    dk    dk    dk
c                                        x     y     z
c
      subroutine dwavefproducts(dcprod,kpt1,nkpt1,ispin1,ispin2,bandi1,bandf1,bandi2,bandf2,lwrite,lacc)

      use global
      use wrapper
      use, intrinsic :: iso_fortran_env
      implicit none
      integer,    intent(in)    :: nkpt1,ispin1,ispin2,bandi1,bandf1,bandi2,bandf2,kpt1(nkpt1)
      logical,    intent(in)    :: lwrite,lacc
      MCOMPLEX_dp,intent(inout) :: dcprod(bandi2:bandf2,bandi1:bandf1,nkpt1,3)
      MCOMPLEX_dp               :: moment(bandi2:bandf2,bandi1:bandf1,nkpt1,3)
      real_dp                   :: rdum
      integer                   :: ikpt,ikpt1,iband1,iband2
      real                      :: time1,time2
      call cpu_time(time1)

      if(.not.lacc) dcprod = 0

c                                            __
c     Get momentum-matrix elements -i < uj | \/ | ui >
      call momentum_matrix(moment,kpt1,nkpt1,ispin1,ispin2,bandi1,bandf1,bandi2,bandf2,.false. MpiC(.false.) LoadC(.false.) )

      do ikpt1 = 1,nkpt1
        ikpt = kpt1(ikpt1)
c                                                  __
c       Calculate expansion coefficients -i < uj | \/ | ui > / ( ei - ej ) for periodic function ui
        do iband1 = bandi1,bandf1
          do iband2 = bandi2,bandf2
            rdum = ( ene(iband2,ikpt,min(ispin2,nspin1)) - ene(iband1,ikpt,min(ispin1,nspin1)) ) * sqrt(vol)
            if(abs(rdum)>1d-8) then
              dcprod(iband2,iband1,ikpt1,:) = dcprod(iband2,iband1,ikpt1,:) + moment(iband2,iband1,ikpt1,:) / rdum
            endif
          enddo
        enddo

      enddo ! end loop over k-points

      call cpu_time(time2);if(lwrite) write(6,'(F8.2'NoA) time2-time1;call cpu_time(time1)

      end

c     -----------------

c     Calculates the momentum matrix elements
c
c                                    s             s'
c     momentum(n',n,q,xyz) = -i < phi  | d    | phi    >
c                                    nq   xyz      n'q
c
c     s  = ispin1 , bandi1 <= n  <= bandf1
c     s' = ispin2 , bandi2 <= n' <= bandf2
c
c     lpar  = .true.  : parallel execution
c     lpar  = .false. : serial   execution
c     lload = .true.  : cmtq,cpwq
c     lload = .false. : cmt,cpw
c
c     Note that the momentum matrix is not completely Hermitian because of the slight mismatch at the MT boundaries (error ~ 1e-5).

# include "cppmacro.h"

      recursive subroutine momentum_matrix(momentum,kpt1,nkpt1,ispin1,ispin2,bandi1,bandf1,bandi2,bandf2,lacc
     &                                     MpiC(lpar) LoadC(lload) )

      use global
      use wrapper
      Mpi( use Mwrapper )
      use, intrinsic :: iso_fortran_env
      implicit none
      integer,     intent(in)    :: nkpt1,ispin1,ispin2,bandi1,bandf1,bandi2,bandf2,kpt1(nkpt1)
      logical,     intent(in)    :: lacc MpiC(lpar) LoadC(lload)
      MCOMPLEX_dp, intent(inout) :: momentum(bandi2:bandf2,bandi1:bandf1,nkpt1,3)
      complex_dp                 :: cvec1(maxlmindx),cvec2(maxlmindx),cvec3(maxlmindx),cvec(3)
      complex_dp                 :: cmt1(maxlmindx)!,cmt2(maxlmindx)
      complex_dp,  pointer_cnt   :: cmt2(:)
# ifdef DCHECK_ifort18bug
#   warning removed contiguous attribute
      MCOMPLEX_dp, pointer       :: vec2(:)
# else
      MCOMPLEX_dp, pointer_cnt   :: vec2(:)
# endif
      MCOMPLEX_dp, allocatable   :: vec1(:)
      real_dp,     allocatable   :: qg(:,:)
      real_dp                    :: kp(3),rot(3,3)
      real_dp                    :: qmat1(maxindx,maxindx,0:maxlcut,ntype),dbas1(maxgrid)
      real_dp                    :: qmat2(maxindx,maxindx,0:maxlcut,ntype),dbas2(maxgrid)
      real_dp                    :: fcoeff((maxlcut+1)**2,-1:1),gcoeff((maxlcut+1)**2,-1:1)
      integer                    :: itype,ieq,ic,l,m,lm,n1,n2,ikpt,iband1,iband2
      integer                    :: lm_0,lm_1,lm0,lm1,lm2,lm3,n0,nn,ikpt1,ikindx
      real_dp                    :: intgrf
      integer                    :: s1,s2 ! spin index for cmt and cpw

      s1 = ispin1 ; s2 = ispin2
# ifdef LOAD
      if(lacc) then ; s1 = 2 ; s2 = 2
      else          ; s1 = 1 ; s2 = 1
      endif
# endif

      if(any(kindx(kpt1)>size(cmt,4))) Bug('k point indices exceed cmt allocation.')

c
c     Define coefficients F and G
      lm = 0
      do l = 0,maxlcut
        do m = -l,l
          lm = lm + 1
          fcoeff(lm,-1) = - sqrt( 1d0 * (l+m+1)*(l+m+2) / (2*(2*l+1)*(2*l+3)) )
          fcoeff(lm, 0) =   sqrt( 1d0 * (l-m+1)*(l+m+1) /   ((2*l+1)*(2*l+3)) )
          fcoeff(lm, 1) = - sqrt( 1d0 * (l-m+1)*(l-m+2) / (2*(2*l+1)*(2*l+3)) )
          gcoeff(lm,-1) =   sqrt( 1d0 * (l-m)*(l-m-1)   / (2*(2*l-1)*(2*l+1)) )
          gcoeff(lm, 0) =   sqrt( 1d0 * (l-m)*(l+m)     /   ((2*l-1)*(2*l+1)) )
          gcoeff(lm, 1) =   sqrt( 1d0 * (l+m)*(l+m-1)   / (2*(2*l-1)*(2*l+1)) )
        enddo
      enddo

c
c     Calculate olap int r**2*u*u' + w * int r*u*u, w = -l,l+1 ( -> qmat1/2 )
      qmat1 = 0
      qmat2 = 0
      do itype = 1,ntype
        do l = 0,lcut(itype)
          do n2 = 1,nindx(l,itype)
            call derivative(dbas1,bas1(:,n2,l,itype,ispin2),itype) ; dbas1 = dbas1 - bas1(:,n2,l,itype,ispin2)/rgrid(:,itype)
            call derivative(dbas2,bas2(:,n2,l,itype,ispin2),itype) ; dbas2 = dbas2 - bas2(:,n2,l,itype,ispin2)/rgrid(:,itype)
            if(l/=0) then
              do n1 = 1,nindx(l-1,itype)
                qmat1(n1,n2,l,itype) = intgrf(   dbas1(:)                  * bas1(:,n1,l-1,itype,ispin1) +
     &                                           dbas2(:)                  * bas2(:,n1,l-1,itype,ispin1) , itype )
     &                       + (l+1) * intgrf( ( bas1(:,n2,l,itype,ispin2) * bas1(:,n1,l-1,itype,ispin1) +
     &                                           bas2(:,n2,l,itype,ispin2) * bas2(:,n1,l-1,itype,ispin1) ) / rgrid(:,itype) ,
     &                                           itype )
              enddo
            endif
            if(l/=lcut(itype)) then
              do n1 = 1,nindx(l+1,itype)
                qmat2(n1,n2,l,itype) = intgrf(   dbas1(:)                  * bas1(:,n1,l+1,itype,ispin1) +
     &                                           dbas2(:)                  * bas2(:,n1,l+1,itype,ispin1) , itype )
     &                           - l * intgrf( ( bas1(:,n2,l,itype,ispin2) * bas1(:,n1,l+1,itype,ispin1) +
     &                                           bas2(:,n2,l,itype,ispin2) * bas2(:,n1,l+1,itype,ispin1) ) / rgrid(:,itype) ,
     &                                           itype )
              enddo
            endif
          enddo
        enddo
      enddo

c                                                  __
c     Calculate momentum matrix elements -i < uj | \/ | ui > wrt wave functions u (->momentum)

      if(.not.lacc) momentum = 0d0

      do ikpt1 = 1,nkpt1
        ikpt   = kpt1(ikpt1) ; if(storeibz) ikpt = kptp(ikpt)
        ikindx = kindx(ikpt)
        kp     = kpt(:,ikpt)

c       MT contribution

        ic = 0
        do itype = 1,ntype
          if(lcut(itype)<0) cycle
          do ieq = 1,neq(itype)
            ic = ic + 1
            nn = sum( [ ((2*l+1)*nindx(l,itype),l=0,lcut(itype)) ] )
            do iband1 = bandi1,bandf1 ; Mpi( if(lpar.and.mod(iband1,Msize)/=Mrank) cycle )

              Load( if(lload) then ; cmt1  = cmtq(:,ic,iband1,s1) ; else )
              cmt1  = cmt(:,ic,iband1,ikindx,s1)
              Load( endif )
              cvec1 = 0 ; cvec2 = 0 ; cvec3 = 0
              ! build up left vector(s) ( -> cvec1/2/3 )
              lm_0 = 0              ! we start with s-functions (l=0)
              lm_1 = nindx(0,itype) ! we start with p-functions (l=0+1)
              lm   = 0
              do l = 0,lcut(itype)-1
                n0 = nindx(l,  itype)
                n1 = nindx(l+1,itype)
                do m = -l,l
                  lm  = lm + 1
                  lm0 = lm_0 + (m+l)    *n0
                  lm1 = lm_1 + (m+1+l+1)*n1
                  lm2 = lm_1 + (m  +l+1)*n1
                  lm3 = lm_1 + (m-1+l+1)*n1
                  cvec1(lm0+1:lm0+n0) = fcoeff(lm,-1) * matmul ( cmt1(lm1+1:lm1+n1) , qmat2(:n1,:n0,l,itype) )
                  cvec2(lm0+1:lm0+n0) = fcoeff(lm, 0) * matmul ( cmt1(lm2+1:lm2+n1) , qmat2(:n1,:n0,l,itype) )
                  cvec3(lm0+1:lm0+n0) = fcoeff(lm, 1) * matmul ( cmt1(lm3+1:lm3+n1) , qmat2(:n1,:n0,l,itype) )
                enddo
                lm_0 = lm_0 + (2*l+1)*n0
                lm_1 = lm_1 + (2*l+3)*n1
              enddo
              lm_0 = nindx(0,itype) ! we start with p-functions (l=1)
              lm_1 = 0              ! we start with s-functions (l=1-1)
              lm   = 1
              do l = 1,lcut(itype)
                n0 = nindx(l,  itype)
                n1 = nindx(l-1,itype)
                do m = -l,l
                  lm  = lm + 1
                  lm0 = lm_0 + (m+l)    *n0
                  lm1 = lm_1 + (m+1+l-1)*n1
                  lm2 = lm_1 + (m  +l-1)*n1
                  lm3 = lm_1 + (m-1+l-1)*n1
                  if(abs(m+1)<=l-1) cvec1(lm0+1:lm0+n0) = cvec1(lm0+1:lm0+n0) +
     &                              gcoeff(lm,-1) * matmul ( cmt1(lm1+1:lm1+n1) , qmat1(:n1,:n0,l,itype) )
                  if(abs(m)  <=l-1) cvec2(lm0+1:lm0+n0) = cvec2(lm0+1:lm0+n0) +
     &                              gcoeff(lm, 0) * matmul ( cmt1(lm2+1:lm2+n1) , qmat1(:n1,:n0,l,itype) )
                  if(abs(m-1)<=l-1) cvec3(lm0+1:lm0+n0) = cvec3(lm0+1:lm0+n0) +
     &                              gcoeff(lm, 1) * matmul ( cmt1(lm3+1:lm3+n1) , qmat1(:n1,:n0,l,itype) )
                enddo
                lm_0 = lm_0 + (2*l+1)*n0
                lm_1 = lm_1 + (2*l-1)*n1
              enddo
              ! multiply with right vector
              do iband2 = bandi2,bandf2
                Load( if(lload) then ; cmt2 => cmtq(:,ic,iband2,s2) ; else )
                cmt2    => cmt(:,ic,iband2,ikindx,s2)
                Load( endif )
                cvec(1) = dotprod(cvec1(:nn),cmt2(:nn)) / img
                cvec(2) = dotprod(cvec2(:nn),cmt2(:nn)) / img
                cvec(3) = dotprod(cvec3(:nn),cmt2(:nn)) / img
                ! add up and transform to cartesian coordinates
                momentum(iband2,iband1,ikpt1,1) = momentum(iband2,iband1,ikpt1,1) + (cvec(1)-cvec(3)) / sqrt(2d0)
                momentum(iband2,iband1,ikpt1,2) = momentum(iband2,iband1,ikpt1,2) - (cvec(1)+cvec(3)) / sqrt(2d0) * img
                momentum(iband2,iband1,ikpt1,3) = momentum(iband2,iband1,ikpt1,3) +  cvec(2)
              enddo
            enddo

          enddo
        enddo
        
c       plane-wave contribution

        allocate ( vec1(ngpt(ikpt)),qg(ngpt(ikpt),3) )

        do nn = 1,ngpt(ikpt)
          qg(nn,:) = matmul(rlat,kp+gpt(:,pgpt(nn,ikpt)))
        enddo

        nn = ngpt(ikpt)
        do iband1 = bandi1,bandf1 ; Mpi( if(lpar.and.mod(iband1,Msize)/=Mrank) cycle )
          Load( if(lload) then ; vec2 => cpwq(:nn,iband1,s1) ; else )
          vec2 => cpw(:nn,iband1,ikindx,s1)
          Load( endif )
          call olapvec(vec1,vec2,gpt(:,pgpt(:ngpt(ikpt),ikpt)),ngpt(ikpt))
          do iband2 = bandi2,bandf2
            Load( if(lload) then ; vec2 => cpwq(:nn,iband2,s2) ; else )
            vec2 => cpw(:nn,iband2,ikindx,s2)
            Load( endif )
            momentum(iband2,iband1,ikpt1,1) = momentum(iband2,iband1,ikpt1,1) + dotprod(vec1,vec2*qg(:,1))
            momentum(iband2,iband1,ikpt1,2) = momentum(iband2,iband1,ikpt1,2) + dotprod(vec1,vec2*qg(:,2))
            momentum(iband2,iband1,ikpt1,3) = momentum(iband2,iband1,ikpt1,3) + dotprod(vec1,vec2*qg(:,3))
          enddo
        enddo
        
        deallocate ( vec1,qg )
      enddo ! end loop over k-points

c
c     Spin loop for SOC

      if(l_soc) then
c        call teststop('SOC summation in momentum_matrix')
        if(lacc) then
          if(ispin1/=2.or.ispin2/=2) Bug('Error in spin indices (second SOC run).')
        else
          if(ispin1/=1.or.ispin2/=1) Bug('Error in spin indices (first SOC run).')
          call momentum_matrix(momentum,kpt1,nkpt1,2,2,bandi1,bandf1,bandi2,bandf2,.true. MpiC(lpar) LoadC(lload) )
          return ! trafo and distribution (see below) already done by previous call
        endif
      endif

c
c     kpt transformation (storeibz kptp->kpt)

      if(storeibz) then
        do ikpt1 = 1,nkpt1
          ikpt = kpt1(ikpt1)
          rot  = matmul(lat,matmul(sym(symkpt(ikpt))%rot,transpose(rlat)))/(2*pi)          
          do iband1 = bandi1,bandf1 ; Mpi( if(lpar.and.mod(iband1,Msize)/=Mrank) cycle )
            do iband2 = bandi2,bandf2
              momentum(iband2,iband1,ikpt1,:) = matmul ( rot , momentum(iband2,iband1,ikpt1,:) )
              NoInv( if(symkpt(ikpt)>nsymt) momentum(iband2,iband1,ikpt1,:) = -conjg( momentum(iband2,iband1,ikpt1,:) ) )
            enddo
          enddo
        enddo
      endif

c
c     Distribution among processes (if lpar=.true.)
      
# ifdef MPI
      if(lpar) then
        Mdistr(momentum(:,iband1,:,:),iband1,bandi1,bandf1)
      endif
# endif

      end

# if 0
      !test of use_sym=.false.
      subroutine check_mom
      use global
      implicit none
      MCOMPLEX_dp :: moment(3),moment1(3)
      real_dp :: rot(3,3),rdum
      integer i,j,ikpt
      do ikpt = 1,nkpt
        write(*,*) 'ikpt',ikpt,kptp(ikpt),kindx(ikpt)
        do i = 1,5
          do j = 1,5
            call momentum_matrix(moment,(/ikpt/),1,1,1,i,i,j,j,.false., MpiC(.false.) LoadC(.false.) )
            call momentum_matrix(moment1,(/kptp(ikpt)/),1,1,1,i,i,j,j,.false., MpiC(.false.) LoadC(.false.) )
            rot  = matmul(lat,matmul(sym(symkpt(ikpt))%rot,transpose(rlat)))/(2*pi)
            moment1 = matmul ( rot,moment1 )
            rdum = rdum + sum(abs(moment-moment1))
          enddo
        enddo
        write(*,*) rdum
        write(*,*)
      enddo
      stop
      end
# endif
      
c     -----------------      

      subroutine momentum_core(momentum,kpt1,nkpt1,ispin1,ispin2,ic,lc,mc,nc,bandi2,bandf2)

      use global
      use wrapper
      use, intrinsic :: iso_fortran_env
      implicit none
      integer,    intent(in)  :: nkpt1,ispin1,ispin2,ic,lc,mc,nc,bandi2,bandf2,kpt1(nkpt1)
      complex_dp, intent(out) :: momentum(bandi2:bandf2,nkpt1,3)
      complex_dp              :: hlp(3,3)
      real_dp                 :: qmat1(maxindx,ntype),dbas1(maxgrid)
      real_dp                 :: qmat2(maxindx,ntype),dbas2(maxgrid)
      real_dp                 :: fcoeff((maxlcut+1)**2,-1:1),gcoeff((maxlcut+1)**2,-1:1)
      integer                 :: itype,ieq,l,m,lm,n2,ikpt,iband2,mm
      integer                 :: lm1,n,ikpt1,ikindx
      real_dp                 :: intgrf

      Bug('Disabled.')

      if(any(kptp(kpt1)/=kpt1)) then
        if(storeibz) then
          Error('only k points in the IBZ are allowed when STOREIBZ is specified.')
        else if(any(kpt1>size(cmt,4))) then
          Bug('k point indices exceed cmt allocation.')
        endif
      endif

c
c     Define coefficients F and G
      lm = 0
      do l = 0,maxlcut
        do m = -l,l
          lm = lm + 1
          fcoeff(lm,-1) = - sqrt( 1d0 * (l+m+1)*(l+m+2) / (2*(2*l+1)*(2*l+3)) )
          fcoeff(lm, 0) =   sqrt( 1d0 * (l-m+1)*(l+m+1) /   ((2*l+1)*(2*l+3)) )
          fcoeff(lm, 1) = - sqrt( 1d0 * (l-m+1)*(l-m+2) / (2*(2*l+1)*(2*l+3)) )
          gcoeff(lm,-1) =   sqrt( 1d0 * (l-m)*(l-m-1)   / (2*(2*l-1)*(2*l+1)) )
          gcoeff(lm, 0) =   sqrt( 1d0 * (l-m)*(l+m)     /   ((2*l-1)*(2*l+1)) )
          gcoeff(lm, 1) =   sqrt( 1d0 * (l+m)*(l+m-1)   / (2*(2*l-1)*(2*l+1)) )
        enddo
      enddo

c
c     Calculate olap int r**2*u*u' + w * int r*u*u, w = -l,l+1 ( -> qmat1/2 )
      qmat1 = 0
      qmat2 = 0
      do itype = 1,ntype
        l = lc + 1
        if(l<=lcut(itype)) then
          do n2 = 1,nindx(l,itype)
            call derivative(dbas1,bas1(:,n2,l,itype,ispin2),itype) ; dbas1 = dbas1 - bas1(:,n2,l,itype,ispin2)/rgrid(:,itype)
            call derivative(dbas2,bas2(:,n2,l,itype,ispin2),itype) ; dbas2 = dbas2 - bas2(:,n2,l,itype,ispin2)/rgrid(:,itype)
            qmat1(n2,itype) = intgrf(   dbas1(:)                  * core1(:,nc,lc,itype,ispin1) +
     &                                  dbas2(:)                  * core2(:,nc,lc,itype,ispin1) , itype )
     &              + (l+1) * intgrf( ( bas1(:,n2,l,itype,ispin2) * core1(:,nc,lc,itype,ispin1) +
     &                                  bas2(:,n2,l,itype,ispin2) * core2(:,nc,lc,itype,ispin1) ) / rgrid(:,itype) , itype )
          enddo
        endif
        l = lc - 1
        if(l>=0) then
          do n2 = 1,nindx(l,itype)
            call derivative(dbas1,bas1(:,n2,l,itype,ispin2),itype) ; dbas1 = dbas1 - bas1(:,n2,l,itype,ispin2)/rgrid(:,itype)
            call derivative(dbas2,bas2(:,n2,l,itype,ispin2),itype) ; dbas2 = dbas2 - bas2(:,n2,l,itype,ispin2)/rgrid(:,itype)
            qmat2(n2,itype) = intgrf(   dbas1(:)                  * core1(:,nc,lc,itype,ispin1) +
     &                                  dbas2(:)                  * core2(:,nc,lc,itype,ispin1) , itype )
     &                  - l * intgrf( ( bas1(:,n2,l,itype,ispin2) * core1(:,nc,lc,itype,ispin1) +
     &                                  bas2(:,n2,l,itype,ispin2) * core2(:,nc,lc,itype,ispin1) ) / rgrid(:,itype) , itype )
          enddo
        endif
      enddo

c                                                  __
c     Calculate momentum matrix elements -i < uj | \/ | ui > wrt wave functions u (->momentum)

      momentum = 0d0

      n = 0
      do itype = 1,ntype
        do ieq = 1,neq(itype)
          n = n + 1 ; if(n==ic) goto 1
        enddo
      enddo
      Bug('atom index wrong.')
 1    continue

      do ikpt1 = 1,nkpt1
        ikpt   = kpt1(ikpt1)
        ikindx = kindx(ikpt)

        lm  = 0
        lm1 = 0
        do l = 0,lcut(itype)
          do m = -l,l
            lm = lm + 1
            do n = 1,nindx(l,itype)
              lm1 = lm1 + 1
              do mm = -1,1
                if(mc+mm==m) then
                  if     (l==lc+1) then
                    do iband2 = bandi2,bandf2
                      momentum(iband2,ikpt1,2+mm) = momentum(iband2,ikpt1,2+mm) +
     &                                              gcoeff(lm,mm) * cmt(lm1,ic,iband2,ikindx,ispin2) * qmat1(n,itype) / img
                    enddo
                  else if(l==lc-1) then
                    do iband2 = bandi2,bandf2
                      momentum(iband2,ikpt1,2+mm) = momentum(iband2,ikpt1,2+mm) +
     &                                              fcoeff(lm,mm) * cmt(lm1,ic,iband2,ikindx,ispin2) * qmat2(n,itype) / img
                    enddo
                  endif
                endif
              enddo
            enddo
          enddo
        enddo

        ! Transform to cartesian coordinates
        hlp        =  0
        hlp( 1, 1) =  1d0/sqrt(2d0)
        hlp( 1, 3) = -1d0/sqrt(2d0)
        hlp( 2, 1) = -img/sqrt(2d0)
        hlp( 2, 3) = -img/sqrt(2d0)
        hlp( 3, 2) =  1d0
        do iband2 = bandi2,bandf2
          momentum(iband2,ikpt1,:) = matmul(momentum(iband2,ikpt1,:),transpose(hlp))
        enddo

      enddo ! end loop over k-points

      end

c     -----------------

# if 0

c     Calculates contracted wave-function products:
c
c               s        s
c     < M    phi    | phi    M    >
c        k,I    q,n      q,n  k,J
c
c                         s        s'             s'          s
c     = SUM(n') < M    phi    | phi       >  < phi       | phi    M    >
c                  k,I    q,n      q+k,n'         q+k,n'      q,n  k,J
c
c     returned in contract_mt (MT part) and contract_pw (PW part).
c     The result is independent of k because the phase factor cancels out. There is only the formal k-dependence of the IPW set.
c     IPWs with G vectors gptm(:,:ngptm2) are used.
c
c     k1(:nk1)  : set of q points
c     s1        : spin s
c     ib1       : band n
c     mode = 1  : only MT part
c     mode = 2  : only PW part
c     mode = 3  : both
c     approx(1) : approximate MT part
c     approx(2) : approximate PW part
c
      subroutine contraction(contract_mt,contract_pw,gptmm,ngptmm,k1,nk1,s1,ib1,mode,approx)
      use global
      use wrapper
      use util
      use, intrinsic :: iso_fortran_env
      implicit none
      integer,     intent(in)          :: nk1,k1(nk1),s1,ib1,mode,ngptmm
      integer,     intent(in),  target :: gptmm(3,ngptmm)
      logical,     intent(in)          :: approx(2)
      complex_dp,  intent(out)         :: contract_mt(maxbasp*(maxbasp+1)/2,nk1,ncent)
      MCOMPLEX_dp, intent(out)         :: contract_pw(ngptmm,nk1)
      MCOMPLEX_dp, allocatable         :: cprod(:,:)
      complex_dp,  allocatable         :: cprod1(:,:)
      complex_dp                       :: hlp(nk1)
      real_dp                          :: integral(maxindxm,0:maxlcutm,maxindxm,0:maxlcutm,maxindxm,0:maxlcutm)
      real_dp,     allocatable, target :: basmm(:,:,:,:)
      integer,     allocatable, target :: pgptmm(:,:)
      integer,                  target :: nindxmm(0:2*maxlcutm,ntype)
      integer,     allocatable         :: pnt(:,:,:)
      real_dp,     pointer_cnt         :: basm0(:,:,:,:)
      integer,     pointer_cnt         :: nindxm0(:,:),gptm0(:,:),pgptm0(:,:)
      integer                          :: pindx(0:maxlcutm,maxindxm,0:maxlcutm,maxindxm,0:2*maxlcutm,ntype)
      integer                          :: nbasp0,lcutm0(ntype)
      integer                          :: g(3),ngptm0
      integer                          :: ik,ik1
      integer                          :: itype,ieq,ic
      integer                          :: i,j,l1,l2,l,n1,n2,n,ll,m1,m2,m,lmn1,lmn2,ibas,ibas0,igpt1,igpt2
      real                             :: time1,time2
      real_dp                          :: gaunt1,rdum
      real_dp                          :: gaunt,intgrf

      Error('deprecated.')

      if (mode==1.or.mode==3)                  contract_mt = 0
      if((mode==2.or.mode==3).and.ngptmm>0) contract_pw = 0

c     We construct a "product mixed basis" (suffix mm) temporarily, such that the subroutine wavefproducts2 can be used for the products.
c
c     MT
      ! define pointers and store values of the original mixed basis (suffix m0)
      nindxm0 => nindxm
      basm0   => basm
      nbasp0  =  nbasp
      lcutm0  =  lcutm
      if(mode==2) then
        lcutm = -1 ! switch off MTs
      else if(.not.approx(1)) then
        ! get dimensions for new MT functions
        nindxmm = 0
        nbasp   = 0
        do itype = 1,ntype
          do l1 = 0,lcutm(itype)
            do n1 = 1,nindxm(l1,itype)
              do l2 = 0,l1
                do n2 = 1,nindxm(l2,itype)
                  if(l1==l2.and.n2>n1) cycle
                  do l = abs(l1-l2),l1+l2,2 ! note the Gaunt condition l+l1+l2 = even
                    nindxmm(l,itype)           = nindxmm(l,itype) + 1
                    pindx(l1,n1,l2,n2,l,itype) = nindxmm(l,itype)
                    pindx(l2,n2,l1,n1,l,itype) = nindxmm(l,itype)
                    nbasp                      = nbasp + (2*l+1) * neq(itype)
                  enddo
                enddo
              enddo
            enddo
          enddo
        enddo
        ! define mixed-basis products -> basmm
        maxlcutm = 2*maxlcutm
        maxindxm = maxval(nindxmm)
        allocate ( basmm(maxgrid,maxindxm,0:maxlcutm,ntype) )
        nindxmm = 0
        do itype = 1,ntype
          do l1 = 0,lcutm(itype)
            do n1 = 1,nindxm(l1,itype)
              do l2 = 0,l1
                do n2 = 1,nindxm(l2,itype)
                  if(l1==l2.and.n2>n1) cycle
                  do l = abs(l1-l2),l1+l2,2
                    nindxmm(l,itype)   = nindxmm(l,itype) + 1
                    n                  = nindxmm(l,itype)
                    basmm(:,n,l,itype) = basm(:,n1,l1,itype) * basm(:,n2,l2,itype) / rgrid(:,itype)
                  enddo
                enddo
              enddo
            enddo
          enddo
          lcutm(itype) = 2*lcutm(itype)
        enddo
        ! point to new "product mixed basis" (MT)
        nindxm => nindxmm
        basm   => basmm
      endif
c     IPW
      ! define pointers and store values for original mixed basis
      ngptm0 =  ngptm(1)
      pgptm0 => pgptm
      gptm0  => gptm
      if(mode==1) then
        ngptm(1) = 0 ! switch off PW
      else if(.not.approx(2)) then
        allocate ( pgptmm(ngptmm,1) )
        ngptm(1)    = ngptmm
        maxgptm     = ngptmm
        pgptmm(:,1) = [ (i,i=1,ngptmm) ]
        ! point to new "product mixed basis" (IPW)
        gptm   => gptmm
        pgptm  => pgptmm
      endif
      allocate ( cprod(nk1,nbasp+ngptm(1)) )
      ! use subroutine wavefproducts2
      Error('contraction disabled!')
c      call wavefproducts2(cprod,k1,nk1,1,s1,s1,ib1,ib1,ib1,ib1,.false.,.false.)
c     Obtain final result from cprod
      ! MT
      if(mode/=2) then
        allocate ( cprod1(nk1,nbasp) )
        cprod1 = cprod(:,:nbasp)
# ifdef INV
        do ik = 1,nk1
          call desymmetrize(cprod1(ik,:),nbasp,1,1)
        enddo
# endif
        ibas0 = 0
        ic    = 0
        do itype = 1,ntype
          if(approx(1)) then
            do l1 = 0,lcutm(itype) ; do n1 = 1,nindxm(l1,itype)
              do l2 = 0,lcutm(itype) ; do n2 = 1,nindxm(l2,itype)
                do l = abs(l1-l2),min(l1+l2,lcutm(itype)),2 ; do n = 1,nindxm(l,itype)
                  integral(n,l,n2,l2,n1,l1) =
     &              intgrf ( basm(:,n,l,itype) * basm(:,n1,l1,itype) * basm(:,n2,l2,itype) / rgrid(:,itype), itype )
                enddo ; enddo
              enddo ; enddo
            enddo ; enddo
          endif
          do ieq = 1,neq(itype)
            ic   = ic + 1
            lmn1 = 0
            do l1 = 0,lcutm0(itype)
              do m1 = -l1,l1
                do n1 = 1,nindxm0(l1,itype)
                  lmn1 = lmn1 + 1
                  lmn2 = 0
                  do l2 = 0,lcutm0(itype)
                    do m2 = -l2,l2
                      do n2 = 1,nindxm0(l2,itype)
                        lmn2 = lmn2 + 1 ; if(lmn2>lmn1) cycle ! skip upper triangle
                        hlp  = 0
                        do l = abs(l1-l2),min(l1+l2,lcutm(itype)),2
                          m      = m1 - m2 ; if(abs(m)>l) cycle
                          gaunt1 = gauntarr(1,l,l1,l2,m,m2)
                          if(approx(1)) then
                            ibas = ibas0 + sum ( [ ((2*ll+1)*nindxm(ll,itype),ll=0,l-1)  ] )
     &                                   + nindxm(l,itype) * (l+m)
                            do n = 1,nindxm(l,itype)
                              ibas = ibas + 1
                              hlp  = hlp + gaunt1 * cprod1(:,ibas) * integral(n,l,n2,l2,n1,l1)
                            enddo
                          else
                            n      = pindx(l1,n1,l2,n2,l,itype)
                            ibas   = ibas0 + sum ( [ ((2*ll+1)*nindxm(ll,itype),ll=0,l-1)  ] )
     &                                     + nindxm(l,itype) * (l+m)
     &                                     + n
                            hlp    = hlp + gaunt1 * cprod1(:,ibas)
                          endif
                        enddo
                        contract_mt(lmn1*(lmn1-1)/2+lmn2,:,ic) = conjg(hlp) ! store the upper triangle
                      enddo
                    enddo
                  enddo
                enddo
              enddo
            enddo
            ibas0 = ibas0 + sum ( [ ((2*l+1)*nindxm(l,itype),l=0,lcutm(itype)) ] )
          enddo
        enddo
        deallocate ( cprod1 )
      endif
      ! IPW
      if(mode/=1) then
        if(approx(2)) then
          g = [ (maxval(abs(gptmm(i,:))),i=1,3) ]
          allocate ( pnt(-g(1):g(1),-g(2):g(2),-g(3):g(3)) )
          pnt = 0
          do n = 1,ngptmm
            g = gptmm(:,n)
            pnt(g(1),g(2),g(3)) = n
          enddo
          do n = 1,ngptm(1)
            g = gptm(:,pgptm(n,1))
            i = pnt(g(1),g(2),g(3))
            if(i/=0) contract_pw(i,:) = cprod(:,nbasp+n) / sqrt(vol)
          enddo
          deallocate ( pnt )
        else
          do n = 1,ngptm(1)
            contract_pw(n,:) = cprod(:,nbasp+n) / sqrt(vol)
          enddo
        endif
      endif
      deallocate ( cprod )
      ! undo modification of mixed basis
      lcutm = lcutm0
      if(.not.approx(1).and.mode/=2) then
        nbasp    =  nbasp0
        nindxm   => nindxm0
        basm     => basm0
        maxlcutm =  maxval(lcutm)
        maxindxm =  maxval(nindxm)
        deallocate ( basmm )
      endif
      ngptm(1) = ngptm0
      if(.not.approx(2).and.mode/=1) then
        gptm     => gptm0
        pgptm    => pgptm0
        maxgptm  =  maxval(ngptm)
        deallocate ( pgptmm )
      endif

      end

# endif
