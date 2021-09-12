c--------------------------------------------------------------------------------
c Copyright (c) 2021 Peter Grünberg Institut, Forschungszentrum Jülich, Germany
c This file is part of SPEX and available as free software under the conditions
c of the MIT license as expressed in the LICENSE file in more detail.
c--------------------------------------------------------------------------------

c     Calculates plane-wave overlap matrix OLAP defined by GPT(1:3,1:NGPT).
c     (Muffin-tin spheres are cut out.)

# include "cppmacro.h"

      subroutine olap_pw(olap,gpt,ngpt)

      use global, only: pi,ntype,neq,vol,grid,cent,img
      use, intrinsic :: iso_fortran_env

      implicit none
      integer,     intent(in)  :: ngpt,gpt(3,ngpt)
      MCOMPLEX_dp, intent(out) :: olap(ngpt,ngpt)
      integer                  :: dg(3),i,j,itype,icent,ieq
      real_dp                  :: g,r,fgr
      real_dp                  :: gptnorm

      do i=1,ngpt
        do j=1,i
          dg        = gpt(:,j)-gpt(:,i)
          g         = gptnorm(dg)
          olap(i,j) = 0
          if(g==0) then
            do itype=1,ntype
              r         = grid(itype)%radius
              olap(i,j) = olap(i,j) - neq(itype) * 4*pi*r**3/3 / vol
            enddo
          else
            icent = 0
            do itype=1,ntype
              r   = g * grid(itype)%radius
              fgr = 4*pi* ( sin(r) - r*cos(r) ) /g**3 / vol
              do ieq=1,neq(itype)
                icent     = icent+1
                olap(i,j) = olap(i,j) - fgr * exp( img * 2*pi*dot_product(dg,cent(:,icent)) )
              enddo
            enddo
          endif
          if(i==j) olap(i,j) = olap(i,j) + 1
          olap(j,i) = MCONJG(olap(i,j))
        enddo
      enddo

      end

c ---------------------

c     The same as olap_pw but returns a matrix in packed storage.

      subroutine olap_pwp(olap,gpt,ngpt)

      use global, only: pi,ntype,neq,vol,grid,cent,img
      use, intrinsic :: iso_fortran_env

      implicit none
      integer,     intent(in)  :: ngpt,gpt(3,ngpt)
      MCOMPLEX_dp, intent(out) :: olap(ngpt*(ngpt+1)/2)
      integer                  :: dg(3),i,j,k,itype,icent,ieq
      real_dp                  :: g,r,fgr
      real_dp                  :: gptnorm

      k = 0
      do i=1,ngpt
        do j=1,i
          k       = k + 1
          dg      = gpt(:,i)-gpt(:,j)
          g       = gptnorm(dg)
          olap(k) = 0
          if(g==0) then
            do itype=1,ntype
              r       = grid(itype)%radius
              olap(k) = olap(k) - neq(itype) * 4*pi*r**3/3 / vol
            enddo
          else
            icent = 0
            do itype=1,ntype
              r   = g * grid(itype)%radius
              fgr = 4*pi* ( sin(r) - r*cos(r) ) /g**3 / vol
              do ieq=1,neq(itype)
                icent   = icent+1
                olap(k) = olap(k) - fgr * exp( img * 2*pi*dot_product(dg,cent(:,icent)) )
              enddo
            enddo
          endif
          if(i==j) olap(k) = olap(k) + 1
        enddo
      enddo

      end

c ---------------------

c     Same as above but with a folding vector gfold
      subroutine olap_pw1(olap,gfold)

      use global, only: gpt,pi,ntype,neq,vol,grid,cent,img,ngptall
      use, intrinsic :: iso_fortran_env

      implicit none
      integer,     intent(in)  :: gfold(3)
      MCOMPLEX_dp, intent(out) :: olap(ngptall,ngptall)
      integer                  :: dg(3),i,j,itype,icent,ieq
      real_dp                  :: g,r,fgr
      real_dp                  :: gptnorm

      do i=1,ngptall
        do j=1,ngptall
          dg        = gpt(:,j) - gfold - gpt(:,i)
          g         = gptnorm(dg)
          olap(i,j) = 0
          if(g==0) then
            olap(i,j) = 1
            do itype=1,ntype
              r         = grid(itype)%radius
              olap(i,j) = olap(i,j) - neq(itype) * 4*pi*r**3/3 / vol
            enddo
          else
            icent = 0
            do itype=1,ntype
              r   = g * grid(itype)%radius
              fgr = 4*pi* ( sin(r) - r*cos(r) ) /g**3 / vol
              do ieq=1,neq(itype)
                icent     = icent+1
                olap(i,j) = olap(i,j) - fgr * exp( img * 2*pi*dot_product(dg,cent(:,icent)) )
              enddo
            enddo
          endif
        enddo
      enddo

      end

c ---------------------

      function stepfunction(g)
      use global
      use, intrinsic :: iso_fortran_env
      implicit none
      MCOMPLEX_dp         :: stepfunction
      integer, intent(in) :: g(3)
      real_dp             :: gnorm,r,fgr
      integer             :: itype,ieq,icent
      gnorm = sqrt(sum(matmul(rlat,g)**2))
c      if(gnorm==0) then ; stepfunction = 1
c      else              ; stepfunction = 0
c      endif
c      return
      if(gnorm==0) then
        stepfunction = 1
        do itype = 1,ntype
          r            = grid(itype)%radius
          stepfunction = stepfunction - neq(itype) * 4*pi*r**3/3 / vol
        enddo
      else
        stepfunction = 0
        icent        = 0
        do itype = 1,ntype
          r   = gnorm * grid(itype)%radius
          fgr = 4*pi* ( sin(r) - r*cos(r) ) / gnorm**3 / vol
          do ieq = 1,neq(itype)
            icent        = icent + 1
            stepfunction = stepfunction - fgr * exp( - img * 2*pi*dot_product(cent(:,icent),g) )
          enddo
        enddo
      endif
      end

c ---------------------

c     Returns wavefunction overlap < cmt1/cpw1 | cmt2/cpw2 >
c     Requires PW overlap olappw=<G|G'> and, if ikpt1/=ikpt2, MT overlap olapmt. (olapmt is not referenced if ikpt1/=ikpt2)
      function wfolap(cmt1,cpw1,cmt2,cpw2,ikpt1,ikpt2,ispin,olapmt,olappw)
      use global
      use wrapper
      use, intrinsic :: iso_fortran_env
      implicit none
      complex_dp             :: wfolap
      integer,    intent(in) :: ispin,ikpt1,ikpt2
      complex_dp, intent(in) :: cmt1(maxlmindx,ncent),cmt2(maxlmindx,ncent)
      MCOMPLEX_dp,intent(in) :: cpw1(ngpt(ikpt1)),cpw2(ngpt(ikpt2))
      MCOMPLEX_dp,intent(in) :: olappw(ngpt(ikpt1),ngpt(ikpt2))
      complex_dp, intent(in) :: olapmt(maxlmindx,maxlmindx,ncent,nspin)
      integer                :: itype,ieq,ic,l,lm,m,n      
      real_dp                :: olapmt0(maxindx,maxindx,0:maxlcut)
      real_dp                :: intgrf
      wfolap = 0
      ic     = 0
      do itype = 1,ntype
        if(ikpt1==ikpt2) then
          do l = 0,lcut(itype)
            do n = 1,nindx(l,itype)
              do m = 1,n
                olapmt0(m,n,l) = intgrf ( bas1(:,m,l,itype,ispin)*bas1(:,n,l,itype,ispin)+
     &                                    bas2(:,m,l,itype,ispin)*bas2(:,n,l,itype,ispin) , itype )
                olapmt0(n,m,l) = olapmt0(m,n,l)
              enddo
            enddo
          enddo
        endif
        do ieq = 1,neq(itype)
          ic = ic + 1
          if(ikpt1==ikpt2) then
            lm = 0
            do l = 0,lcut(itype)
              n = nindx(l,itype)
              do m = -l,l
                wfolap = wfolap + dot_product ( cmt1(lm+1:lm+n,ic) , matmul( olapmt0(:n,:n,l) , cmt2(lm+1:lm+n,ic) ) )
                lm     = lm + n
              enddo
            enddo
          else
            lm     = sum ( [ (nindx(l,itype)*(2*l+1),l=0,lcut(itype)) ] )
            wfolap = wfolap + dot_product ( cmt1(:lm,ic) , matmul(olapmt(:lm,:lm,ic,ispin),cmt2(:lm,ic)) )
          endif
        enddo
      enddo
      wfolap = wfolap + dotprod(cpw1,matvec(olappw,cpw2))
      end

c     ------------------

c     Returns PW overlap matrix in olappw=<G|G'>; kpoints can be different: ikpt1, ikpt2
c     
      subroutine olap_gpt(olappw,dim,ikpt1,ikpt2)
      use global, only: gpt,ngpt,pgpt,cstep
      use, intrinsic :: iso_fortran_env
      implicit none
      integer,     intent(in)  :: ikpt1,ikpt2,dim
      MCOMPLEX_dp, intent(out) :: olappw(dim,*)
      integer                  :: i,j,g(3)
      if(dim<ngpt(ikpt1)) Error('First dimension too small.')
      if(.not.allocated(cstep)) Bug('Array cstep not allocated.')
      do j = 1,ngpt(ikpt2)
        if(ikpt1==ikpt2) then
          do i = 1,j
            g           = gpt(:,pgpt(i,ikpt1)) - gpt(:,pgpt(j,ikpt2))
            olappw(i,j) = cstep(g(1),g(2),g(3))
            olappw(j,i) = MCONJG( olappw(i,j) )
          enddo
        else
          do i = 1,ngpt(ikpt1)
            g           = gpt(:,pgpt(i,ikpt1)) - gpt(:,pgpt(j,ikpt2))
            olappw(i,j) = cstep(g(1),g(2),g(3))
          enddo
        endif
      enddo
      end

      ! Same in packed storage (and identical kpoints)
      subroutine olapp_gpt(olappw,ikpt)
      use global, only: gpt,ngpt,pgpt,cstep
      use, intrinsic :: iso_fortran_env
      implicit none
      integer,     intent(in)  :: ikpt
      MCOMPLEX_dp, intent(out) :: olappw(*)
      integer                  :: i,j,k,g(3)
      call teststop('Routine olapp_gpt') 
      if(.not.allocated(cstep)) Bug('Array cstep not allocated.')
      k = 0
      do j = 1,ngpt(ikpt)
        do i = 1,j          
          g         = gpt(:,pgpt(i,ikpt)) - gpt(:,pgpt(j,ikpt))
          k         = k + 1
          olappw(k) = cstep(g(1),g(2),g(3))
        enddo
      enddo
      end

      ! Same for mixed-basis overlap.
      subroutine olap_gptm(olappw,dim,ikpt1,ikpt2)
      use global, only: gptm,ngptm,pgptm,cstep
      use, intrinsic :: iso_fortran_env
      implicit none
      integer,     intent(in)  :: ikpt1,ikpt2,dim
      MCOMPLEX_dp, intent(out) :: olappw(dim,*)
      integer                  :: i,j,g(3)
      if(dim<ngptm(ikpt1)) Error('First dimension too small.')
      if(.not.allocated(cstep)) Bug('Array cstep not allocated.')
      do j = 1,ngptm(ikpt2)
        do i = 1,ngptm(ikpt1)
          g           = gptm(:,pgptm(i,ikpt1)) - gptm(:,pgptm(j,ikpt2))
          olappw(i,j) = cstep(g(1),g(2),g(3))
        enddo
      enddo
      end

c     ------------------

c     Returns out = matmul(olap,vec), where olap is defined according to gpt(:ngpt).
c     (Does not require large pw overlap matrix.)
c
      subroutine olapvec(out,vec,gpt,ngpt)
      use wrapper
      use global, only: cstep
      use, intrinsic :: iso_fortran_env
      implicit none
      integer,     intent(in)  :: ngpt,gpt(3,ngpt)
      MCOMPLEX_dp, intent(in)  :: vec(ngpt)
      MCOMPLEX_dp, intent(out) :: out(ngpt)
      MCOMPLEX_dp              :: row(ngpt)
      integer                  :: i,j,g(3)
      if(.not.allocated(cstep)) Bug('Array cstep not allocated.')
      do i = 1,ngpt
        do j = 1,ngpt
          g      = gpt(:,j) - gpt(:,i)
          row(j) = cstep(g(1),g(2),g(3))
        enddo
        out(i) = dotprod(row,vec)
      enddo
      end

      ! "inplace" routine
      subroutine olapvec1(vec,gpt,ngpt)
      use wrapper
      use global, only: cstep
      use, intrinsic :: iso_fortran_env
      implicit none
      integer,     intent(in)    :: ngpt,gpt(3,ngpt)
      MCOMPLEX_dp, intent(inout) :: vec(ngpt)
      MCOMPLEX_dp                :: row(ngpt),out(ngpt)
      integer                    :: i,j,g(3)
      do i = 1,ngpt
        do j = 1,ngpt
          g      = gpt(:,j) - gpt(:,i)
          row(j) = cstep(g(1),g(2),g(3))
        enddo
        out(i) = dotprod(row,vec)
      enddo
      vec = out
      end

c     ------------------

c     Returns MT part of wavefunction overlap.
c       kzero=.true.  : kpoints are identical (olapmt not referenced)
c       kzero=.false. : kpoints are different (olapmt required)
c
      function wfolap_mt(cmt1,cmt2,ispin,olapmt,kzero)
      use global
      use, intrinsic :: iso_fortran_env
      implicit none
      complex_dp             :: wfolap_mt
      complex_dp, intent(in) :: cmt1(maxlmindx,ncent),cmt2(maxlmindx,ncent)
      complex_dp, intent(in) :: olapmt(maxlmindx,maxlmindx,ncent,nspin)
      logical,    intent(in) :: kzero
      integer,    intent(in) :: ispin
      integer                :: itype,ieq,ic,l,lm,m,n
      real_dp                :: olapmt0(maxindx,maxindx,0:maxlcut)
      real_dp                :: intgrf
      wfolap_mt = 0
      ic        = 0
      do itype = 1,ntype
        if(kzero) then
          do l = 0,lcut(itype)
            do n = 1,nindx(l,itype)
              do m = 1,n
                olapmt0(m,n,l) = intgrf ( bas1(:,m,l,itype,ispin)*bas1(:,n,l,itype,ispin)+
     &                                    bas2(:,m,l,itype,ispin)*bas2(:,n,l,itype,ispin) , itype )
                olapmt0(n,m,l) = olapmt0(m,n,l)
              enddo
            enddo
          enddo
        endif
        do ieq = 1,neq(itype)
          ic = ic + 1
          if(kzero) then
            lm = 0
            do l = 0,lcut(itype)
              n = nindx(l,itype)
              do m = -l,l
                wfolap_mt = wfolap_mt + dot_product ( cmt1(lm+1:lm+n,ic) , matmul( olapmt0(:n,:n,l) , cmt2(lm+1:lm+n,ic) ) )
                lm        = lm + n
              enddo
            enddo
          else
            lm        = sum ( [ (nindx(l,itype)*(2*l+1),l=0,lcut(itype)) ] )
            wfolap_mt = wfolap_mt + dot_product ( cmt1(:lm,ic) , matmul(olapmt(:lm,:lm,ic,ispin),cmt2(:lm,ic)) )
          endif
        enddo
      enddo
      end

c     ------------------

      subroutine wfolap_mt2(cmt1,olapmt,cmt2,ispin,kzero)
      use global
      use, intrinsic :: iso_fortran_env
      implicit none
      complex_dp, intent(out) :: cmt1(maxlmindx,ncent)
      complex_dp, intent(in)  :: cmt2(maxlmindx,ncent)
      complex_dp, intent(in)  :: olapmt(maxlmindx,maxlmindx,ncent,nspin)
      logical,    intent(in)  :: kzero
      integer,    intent(in)  :: ispin
      integer                 :: itype,ieq,ic,l,lm,m,n
      ic = 0
      do itype = 1,ntype
        do ieq = 1,neq(itype)
          ic = ic + 1
          if(kzero) then
            lm = 0
            do l = 0,lcut(itype)
              n = nindx(l,itype)
              do m = -l,l
                cmt1(lm+1:lm+n,ic) = matmul(olapmt(lm+1:lm+n,lm+1:lm+n,ic,ispin),cmt2(lm+1:lm+n,ic))
                lm                 = lm + n
              enddo
            enddo
          else
            lm           = sum ( [ (nindx(l,itype)*(2*l+1),l=0,lcut(itype)) ] )
            cmt1(:lm,ic) = matmul(olapmt(:lm,:lm,ic,ispin),cmt2(:lm,ic))
          endif
        enddo
      enddo
      end

c     ------------------

c     Initializes (MT) wavefunction overlap calculations with wfolap() or wfolap_mt()
c     Only required if ikpt1/=ikpt2  in wfolap()
c                  and kzero=.false. in wfolap_mt(),
c                  in which case kpoint is not the nullvector.      
c
      subroutine wfolap_init_mt(olapmt,kpoint)
      use global
      use, intrinsic :: iso_fortran_env
      implicit none
      complex_dp,  intent(out) :: olapmt(maxlmindx,maxlmindx,ncent,nspin)
      real_dp,     intent(in)  :: kpoint(3)
      integer                  :: ispin,itype,ieq,ic,l,l1,l2,m1,m2,nn,n1,n2,lm,lm1,lm2,m
      complex_dp               :: y((2*maxlcut+1)**2),cexp
      real_dp                  :: r,k,sphbes(maxgrid,0:2*maxlcut,ntype),gnt
      real_dp                  :: integral(maxindx,0:maxlcut,maxindx,0:maxlcut,0:2*maxlcut)
      real_dp                  :: intgrf,gaunt
      if(all(kpoint==0)) Info('wfolap_init_mt called with nullvector.')
      k = sqrt(sum(matmul(rlat,kpoint)**2))
      do itype = 1,ntype
        r = grid(itype)%first ; call sphbessel(sphbes(1,:,itype),k*r,2*lcut(itype))
        do nn = 2,grid(itype)%number
          r = r * exp(grid(itype)%increment) ; call sphbessel(sphbes(nn,:,itype),k*r,2*lcut(itype))
        enddo
      enddo
      call harmonicsr(y,matmul(rlat,kpoint),2*maxlcut)
      olapmt = 0
      do ispin = 1,nspin
        ic = 0
        integral = 0
        do itype = 1,ntype
          do l1 = 0,lcut(itype)
            do n1 = 1,nindx(l1,itype)
              do l2 = 0,lcut(itype)
                do n2 = 1,nindx(l2,itype)
                  do l = abs(l1-l2),l1+l2,2
                    integral(n1,l1,n2,l2,l) = intgrf ( sphbes(:,l,itype) * (
     &                                                 bas1(:,n1,l1,itype,ispin)*bas1(:,n2,l2,itype,ispin)+
     &                                                 bas2(:,n1,l1,itype,ispin)*bas2(:,n2,l2,itype,ispin) ) , itype )
                  enddo
                enddo
              enddo
            enddo
          enddo
          do ieq = 1,neq(itype)
            ic   = ic + 1
            cexp = exp(-img * 2*pi * dot_product(kpoint,cent(:,ic)))
            lm1  = 0
            do l1 = 0,lcut(itype)
              do m1 = -l1,l1
                do n1 = 1,nindx(l1,itype)
                  lm1 = lm1 + 1
                  lm2 = 0
                  do l2 = 0,lcut(itype)
                    do m2 = -l2,l2
                      do n2 = 1,nindx(l2,itype)
                        lm2 = lm2 + 1
                        do l = abs(l1-l2),l1+l2,2
                          m   = m2 - m1 ; if(abs(m)>l) cycle
                          gnt = gaunt(l1,l2,l,m1,m2,m) !; if(abs(gnt-gaunt(l,l2,l1,m,m2,m1))>1d-12) Error('gaunt')
                          lm  = l**2 + l + m + 1
                          olapmt(lm1,lm2,ic,ispin) = olapmt(lm1,lm2,ic,ispin) +
     &                                               4*pi*(-img)**l * y(lm) * cexp * gnt * integral(n1,l1,n2,l2,l)
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
      end

c     ------------------

# if 0

      subroutine wfolap_init(olappw,olapmt,kpoint)
      use global
      use, intrinsic :: iso_fortran_env
      implicit none
      MCOMPLEX_dp, intent(out) :: olappw(ngpt,ngpt)
      complex_dp,  intent(out) :: olapmt(maxlmindx,maxlmindx,ncent,nspin)
      real_dp,     intent(in)  :: kpoint(3)
      call wfolap_init_mt(olapmt,kpoint)
      call olap_pw(olappw,gpt,ngpt)
      end

c     ------------------

      subroutine cpwfold(cpw1,cpw0,gfold)
      use global, only: gpt,ngpt,grid,ntype,neq,cent,vol,img,pi
      use, intrinsic :: iso_fortran_env
      implicit none
      MCOMPLEX_dp, intent(out) :: cpw1(ngpt)
      MCOMPLEX_dp, intent(in)  :: cpw0(ngpt)
      integer,     intent(in)  :: gfold(3)
      MCOMPLEX_dp              :: olap
      integer                  :: i,j,dg(3),itype,ieq,icent
      real_dp                  :: g,r,fgr
      real_dp                  :: gptnorm
      cpw1 = 0
      do i = 1,ngpt
        do j = 1,ngpt
          dg   = gpt(:,j) - gfold - gpt(:,i)
          g    = gptnorm(dg)
          olap = 0
          if(g==0) then
            olap = 1
            do itype = 1,ntype
              r    = grid(itype)%radius
              olap = olap - neq(itype) * 4*pi*r**3/3 / vol
            enddo
          else
            icent = 0
            do itype = 1,ntype
              r   = g * grid(itype)%radius
              fgr = 4*pi* ( sin(r) - r*cos(r) ) /g**3 / vol
              do ieq = 1,neq(itype)
                icent = icent + 1
                olap  = olap - fgr * exp( img * 2*pi*dot_product(dg,cent(:,icent)) )
              enddo
            enddo
          endif
          cpw1(i) = cpw1(i) + olap * cpw0(j)
        enddo
      enddo
      end

c     ------------------

c     Returns the product
c
c     invert(olap1(:,:)) * olap2(:,:)
c
c     where olap1(:,:) = overlap of interstitial plane waves with cutoff   gcutm (mixed basis)
c           olap2(:,:) = overlap of interstitial plane waves with cutoff 2*gcut.
c
c     for the k point ikpt0 or all irreducible k points, if ikpt0 = 0.
c
      subroutine olap_pw_prod(ikpt0)
      use global
      use wrapper
      use, intrinsic :: iso_fortran_env
      implicit none
      integer,     intent(in)  :: ikpt0
      MCOMPLEX_dp, allocatable :: olap(:)
      integer,     allocatable :: ngpt1(:),pgpt1(:,:)
      integer                  :: icycle,igpt,igpt1,igpt2,ikpt,ikpt1,ikpt2,itype,ieq,icent
      integer                  :: g(3),dg(3)
      integer                  :: klow,kup,nkpt1,ik
      real                     :: time
      real_dp                  :: gnorm,r,fgr,kvec(3)
      real_dp                  :: gptnorm

      integer ik1,ik2,gpt2(3),indx(maxgptm)
      MCOMPLEX_dp, allocatable :: olap_prod1(:,:),olap1(:,:),olap0(:,:),olap_prod3(:,:)
      MCOMPLEX_dp :: olap_prod2(maxgptm)
      integer pnt_dim1(3,2),g2(3)
      integer, allocatable :: pnt_prod1(:,:,:),pnt_prod3(:,:,:)
      integer olap_dim1,i,j,kptsum,i1,i2,ig1,ig2,ix,iy,iz
      real time1,time2
      MCOMPLEX_dp cexp(maxgptm),cexp2

      if(cutzero) Warn('Singular-value decomposition perhaps needed but not implemented.')
      if(ikpt0==0) write(6,*)
      write(6,'(A'NoA) 'Calculate overlap-correction matrix... '
      call cpu_time(time)
      if(ikpt0==0) then ; klow = 1     ; kup = nkpti ; if(lkptadd.and.any(job(:)%type>2)) kup = kup + 1
      else                ; klow = ikpt0 ; kup = ikpt0
      endif
      ! Define ngpt1, pgpt1 locally to make the calculation faster
      if(lkptadd) then ; nkpt1 = nkpt*2
      else             ; nkpt1 = nkpt
      endif
      allocate ( ngpt1(nkpt1),pgpt1(ngpt,nkpt1) )
      do ikpt = 1,nkpt1
        if(ikpt>nkpt) then ; kvec = modulo1r(kpt(:,ikpt-nkpt)+kptadd)
        else                  ; kvec = kpt(:,ikpt)
        endif
        i = 0
        do igpt = 1,ngpt
          if(sum(matmul(rlat,kvec+gpt(:,igpt))**2)<=gcut**2) then
            i             = i + 1
            pgpt1(i,ikpt) = igpt
          endif
        enddo
        ngpt1(ikpt) = i
      enddo
      !
      allocate ( olap_dim(klow:kup) )
      g(1) = 2*maxval(abs(gpt(1,:))) + 1
      g(2) = 2*maxval(abs(gpt(2,:))) + 1
      g(3) = 2*maxval(abs(gpt(3,:))) + 1
      allocate ( pnt_prod(-g(1):g(1),-g(2):g(2),-g(3):g(3),klow:kup) )
      allocate ( pnt_prod3(-g(1):g(1),-g(2):g(2),-g(3):g(3)) )
      pnt_prod = 0
      icycle   = 1
      ! first  cycle: determine dimension olap_dim
      ! second cycle: calculate olap_prod and define the pointer pnt_prod
 1    do ik = klow,kup
        if(ik==nkpti+1) then ; ikpt = nkpt + 1
        else                 ; ikpt = ik
        endif
        do igpt = 1,ngptm(ikpt)
          g                           = gptm(:,pgptm(igpt,ikpt))
          pnt_prod(g(1),g(2),g(3),ik) = -1 ! cut out sphere with radius gcutm
        enddo
c       this allows to restrict to a sphere with radius smaller than 2*gcut
c        do iz = lbound(pnt_prod,3),ubound(pnt_prod,3)
c        do iy = lbound(pnt_prod,2),ubound(pnt_prod,2)
c        do ix = lbound(pnt_prod,1),ubound(pnt_prod,1)
c          g = [ix,iy,iz]
c          r = sum(matmul(rlat,kpt(:,ikpt)+g)**2)
c          if(r>(gcut+gcutm/2)**2) then
c            pnt_prod(g(1),g(2),g(3),ik) = -2
c          endif
c        enddo
c        enddo
c        enddo
        i = 0
        do ikpt1 = 1,nkpt
          ikpt2 = kptsum(ikpt1,ikpt)
          if(ikpt==nkpt+1) then ; gpt2 = nint ( modulo1r(kpt(:,ikpt2-nkpt)+kptadd) - (kpt(:,ikpt1) + kptadd) ) !
          else                  ; gpt2 = nint ( kpt(:,ikpt2) - ( kpt(:,ikpt1) + kpt(:,ikpt) ) )                  ! nint avoids rounding errors
          endif
          do ig1 = 1,ngpt1(ikpt1) ; igpt1 = pgpt1(ig1,ikpt1)
          do ig2 = 1,ngpt1(ikpt2) ; igpt2 = pgpt1(ig2,ikpt2)
            g = -gpt(:,igpt1) + gpt(:,igpt2) + gpt2
            if(pnt_prod(g(1),g(2),g(3),ik)==0) then
              i                           = i + 1
              pnt_prod(g(1),g(2),g(3),ik) = i
              if(icycle==2) then
                do igpt = 1,ngptm(ikpt)
                  dg                   = g  - gptm(:,pgptm(igpt,ikpt))
                  gnorm                = gptnorm(dg)
                  olap_prod(igpt,i,ik) = 0
                  icent                = 0
                  if(gnorm==0) Error('zero G vector.')
                  do itype = 1,ntype
                    r   = gnorm * grid(itype)%radius
                    fgr = 4*pi * ( sin(r) - r*cos(r) ) / gnorm**3 / vol
                    do ieq = 1,neq(itype)
                      icent                = icent + 1
                      olap_prod(igpt,i,ik) = olap_prod(igpt,i,ik) - fgr * exp( img * 2*pi*dot_product(dg,cent(:,icent)) )
                    enddo
                  enddo
                enddo
              endif
            endif
          enddo
          enddo
        enddo
        if(icycle==1) olap_dim(ik) = i
      enddo
      if(icycle==1) then
        maxolap_dim = maxval(olap_dim)
        call checkmem('olap_prod',MBYTES*maxgptm*maxolap_dim*(kup-klow+1))
        allocate ( olap_prod(maxgptm,maxolap_dim,klow:kup) )
        pnt_prod = 0
        icycle   = 2
        goto 1
      endif
      ! multiply with inverse of overlap
      do ik = klow,kup
        if(ik==nkpti+1) then ; ikpt = nkpt + 1
        else                 ; ikpt = ik
        endif
        allocate ( olap(ngptm(ikpt)*(ngptm(ikpt)+1)/2) )
        call olap_pwp ( olap,gptm(:,pgptm(:ngptm(ikpt),ikpt)),ngptm(ikpt) )
        call inverse ( olap )
        olap_prod(:ngptm(ikpt),:olap_dim(ik),ik) = matmat(olap,olap_prod(:ngptm(ikpt),:olap_dim(ik),ik) )
        deallocate ( olap )
      enddo
      call cpu_done(time)
      if(ikpt0==0) write(6,'(A,F6.1," MB")') 'Size of overlap-correction matrix:',MBYTES/megabyte*size(olap_prod)
      deallocate ( ngpt1,pgpt1 )
      end

c ---------------------

c     Transforms olap_prod (see above) at ikpt from olap_prod at the parent k point in the IBZ.
c     The result is written to olap_prod1, and the corresponding new pointer is returned in pnt_prod1.
c
      subroutine olap_pw_prod_trafo(olap_prod1,pnt_prod1,ikpt)
      use global
      use, intrinsic :: iso_fortran_env
      implicit none
      MCOMPLEX_dp, intent(out) :: olap_prod1(maxgptm,maxolap_dim)
      integer,     intent(out) :: pnt_prod1(lbound(pnt_prod,1):ubound(pnt_prod,1),
     &                                      lbound(pnt_prod,2):ubound(pnt_prod,2),
     &                                      lbound(pnt_prod,3):ubound(pnt_prod,3))
      integer,     intent(in)  :: ikpt
      MCOMPLEX_dp              :: cexp(maxgptm),cexp2
      logical                  :: done(0:maxolap_dim)
      integer                  :: indx(maxgptm),g(3),g1(3),gpt2(3)
      integer                  :: i,j,isym,isymi,igpt,igpt1,igpt2,ig1,ig2,ikpt1,ikpt2
      integer                  :: lb(3),ub(3),ix,iy,iz
      integer                  :: kptsum
      if(.not.allocated(olap_prod)) Error('olap_prod not allocated.')
      lb(1) = lbound(olap_prod,3)
      ub(1) = ubound(olap_prod,3)
      if(lb(1)>kptp(ikpt).or.ub(1)<kptp(ikpt))
     &                              Error('olap_prod not defined for parent k point.')
      ! copy olap_prod if ikpt lies in the IBZ
      if(ikpt==nkpt+1) then ; ikpt1 = nkpti + 1
      else                    ; ikpt1 = ikpt
      endif
      if(ikpt1>=lb(1).and.ikpt1<=ub(1)) then
        olap_prod1(:ngptm(ikpt),:olap_dim(ikpt1)) = olap_prod(:ngptm(ikpt),:olap_dim(ikpt1),ikpt1)
        pnt_prod1(:,:,:)                          = pnt_prod(:,:,:,ikpt1)
        return
      endif
      ! cut out sphere with radius gcutm
      pnt_prod1 = 0
      do igpt = 1,ngptm(ikpt)
        g                         = gptm(:,pgptm(igpt,ikpt))
        pnt_prod1(g(1),g(2),g(3)) = -1
      enddo
      ! perform trafo
      isym  = symkpt(ikpt) ; lsymmor = any(abs(sym(isym)%transl)>1d-10)
      isymi = sym(isym)%inv
      done  = .false.
      do igpt = 1,ngptm(ikpt)
        g          = matmul(sym(isymi)%rrot,gptm(:,pgptm(igpt,ikpt))-gkptsym(:,kptp(ikpt),isym))
        indx(igpt) = pntgptm(g(1),g(2),g(3),kptp(ikpt))
        cexp(igpt) = exp(-2*pi*img*dot_product(gptm(:,pgptm(igpt,ikpt)),sym(isym)%transl))
      enddo
      lb = lbound(pnt_prod1)
      ub = ubound(pnt_prod1)
      do iz = lb(3),ub(3)
        do iy = lb(2),ub(2)
          do ix = lb(1),ub(1)
            if(pnt_prod1(ix,iy,iz)==0) then
              g  = [ ix,iy,iz ]
              g1 = matmul(sym(isymi)%rrot,g-gkptsym(:,kptp(ikpt),isym))
              if(all(g1>=lb).and.all(g1<=ub)) then
                j = pnt_prod(g1(1),g1(2),g1(3),kptp(ikpt))
                if(j>0.and..not.done(j)) then
                  pnt_prod1(ix,iy,iz) = j
                  if(lsymmor) then
                    olap_prod1(:ngptm(ikpt),j) = olap_prod(indx(:ngptm(ikpt)),j,kptp(ikpt))
                  else
                    cexp2                      = exp(2*pi*img*dot_product(g,sym(isym)%transl))
                    olap_prod1(:ngptm(ikpt),j) = olap_prod(indx(:ngptm(ikpt)),j,kptp(ikpt)) * cexp2 * cexp(:ngptm(ikpt))
                  endif
                  done(j) = .true.
                endif
              endif
            endif
          enddo
        enddo
      enddo
      end

# endif
