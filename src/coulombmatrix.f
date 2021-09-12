c--------------------------------------------------------------------------------
c Copyright (c) 2021 Peter Grünberg Institut, Forschungszentrum Jülich, Germany
c This file is part of SPEX and available as free software under the conditions
c of the MIT license as expressed in the LICENSE file in more detail.
c--------------------------------------------------------------------------------

c     Calculates the Coulomb matrix
c
c     v      =  < M    | v | M    >
c      k,IJ        k,I        k,J
c
c     with the mixed-basis functions M (indices I and J).
c
c     Note that
c                 *
c     v      =  v     .
c      k,JI      k,IJ
c
c     In the code: coulomb1(IJ,k) = v     where only the upper triangle (I<=J) is stored.
c                                    k,IJ
c
c     The Coulomb matrix v(IJ,k) diverges at the Gamma-point. Here, we formally apply the decomposition
c
c              (0)        (1)   *        2-l              (0)*   (0)    (1)*        m  (1)
c     v     = v    + SUM v   * Y  (k) / k        with    v    = v   ,  v      = (-1)  v        .
c      k,IJ    IJ     lm  IJ    lm                        JI     IJ     JI,lm          IJ,l,-m
c
c
c     For the PW contribution we have to construct plane waves inside the MT spheres with the help
c     of spherical Bessel functions. The value lexp (LEXP in spex.inp) is the corresponding cutoff.
c
c     coulomb1 : Output Coulomb matrix in packed storage (:dim*(dim+1)/2,nkpt1)
c     dim      : Maximal dimension of Coulomb matrix
c     nkpt1    : Number of kpoints
c     kpt1     : Kpoint list (in indices of kpt)
c     lwrite   : Detailed output if true, otherwise short output.
c
c     Symmetry is used to accelerate calculation. Can be switched off here:
c # define switch_off_symmetry_coulomb
c
# include "cppmacro.h"

      subroutine coulombmatrix(Win(coulomb1),dim,nkpt1,kpt1,lwrite)
      use global
      use wrapper
      use util
      use file
      use key
      use, intrinsic :: iso_fortran_env
      use, intrinsic :: iso_c_binding
      Mpi( use Mwrapper )
      implicit none
      integer,      intent(in)    :: nkpt1,kpt1(nkpt1),dim MpiC(win_coulomb1)
      logical,      intent(in)    :: lwrite
      MCOMPLEX_dp,  intent(out)   :: coulomb1(dim*(dim+1)/2,nkpt1)
      MCOMPLEX_dp                 :: cint
      MCOMPLEX_dp,  allocatable   :: olapm(:),olaphlp(:,:)
      complex_dp,   pointer_cnt   :: structconst(:,:,:,:),coulmat(:,:)
      complex_dp,   allocatable   :: structconst1(:,:)
      complex_dp,   allocatable   :: y(:),y1(:),y2(:),imgl(:),work(:)
      complex_dp,   allocatable   :: carr(:),carr1(:),carr2(:,:),carr2a(:,:),carr2b(:,:),carr3(:,:,:),step(:),coulstep(:)
      complex_dp                  :: csum,csumf(9),cdum,cdum1,cexp,cexp1(ntype)
      real_dp,      allocatable   :: sphbesmoment(:,:,:),olap(:,:,:,:),integral(:,:,:,:),sphbes0(:,:,:),rarr(:),gmat(:,:)
      real_dp,      allocatable   :: qnrm(:),qnrmstep(:),olapstep(:,:,:,:),csphbes(:,:,:,:),error(:,:),rwork(:),eig(:)
      real_dp                     :: moment(maxindxm,0:maxlcutm,ntype),moment2(maxindxm,ntype)
      real_dp                     :: sphbes(maxgrid,0:maxlcutm),sphbesmoment1(maxgrid,0:maxlcutm)
      real_dp                     :: rarr1(0:maxlcutm),mat(maxindxm*(maxindxm+1)/2)
      real_dp                     :: primf1(maxgrid),primf2(maxgrid),integrand(maxgrid)
      real_dp                     :: svol,qnorm,qnorm1,qnorm2,q(3),q1(3),q2(3),rdum,rdum1,rdum2,gnorm
      integer,      allocatable   :: iarr(:),pointer(:,:,:),gptstep(:,:),pqnrm(:,:),pqnrmstep(:,:)
      integer                     :: g(3),pgptm1(maxgptm,nkpt1),ngptstep,nqnrmstep
      integer                     :: nqnrm,ng
      integer                     :: itype,itype1,itype2,ieq,ieq1,ieq2,ic,ic1,ic2,ikpt
      integer                     :: iqnrm,iqnrm1,iqnrm2,igpt,igpt0,igpt1,igpt2,igptp,igptp1,igptp2
      integer                     :: isym,isym1,nsym1(nkpt1),sym1(nsym,nkpt1),ngptm1(nkpt1)
      integer                     :: l,l1,l2,m,m1,m2,lm,lm1,lm2,n,n1,n2
      integer                     :: i,j,k,ix,iy,iy0,idum
      real                        :: time1,time2
      logical                     :: tstcoul,chkcoul
      real_dp                     :: gcutstep
      real_dp                     :: intgrf,gptnorm
# ifdef MPI
      type(c_ptr)                 :: ptr
      integer                     :: win_structconst,win_coulmat,Merr,mm
# endif

      ! Read spex.inp
      Rbegin
      call getkey(inp,'CHKCOUL',   chkcoul, section='COULOMB', default=.false.,       writeout=.false.)
      call getkey(inp,'TSTCOUL',   tstcoul, section='COULOMB', default=.false.,       writeout=.false.)
      call getkey(inp,'STEPRAD',  gcutstep, section='COULOMB', default=0d0, mini=0d0, writeout=.false.)
      if(tstcoul) chkcoul = .true.
      Rend
      Mpi( call Mcast(chkcoul) ; call Mcast(tstcoul) ; call Mcast(gcutstep) )

      allocate ( y((lexp+1)**2),y1((lexp+1)**2),y2((lexp+1)**2),imgl(0:lexp),rarr(0:max(maxlcutm,lexp)+1),
     &           gmat((lexp+1)**2,(lexp+1)**2) )

      svol = sqrt(vol)

      ! Get symmetry operations that leave kpt(:,ikpt) invariant -> sym1
      if(use_sym) then
        do k = 1,nkpt1 ; ikpt = kpt1(k)
          isym1 = 0
          do isym = 1,nsym
            if(all(abs(matmul(sym(isym)%rrot,kpt(:,ikpt))-kpt(:,ikpt))<1d-12)) then
              isym1         = isym1 + 1
              sym1(isym1,k) = isym
            endif
          enddo
          nsym1(k) = isym1
        enddo
      else
        nsym1 = 1
        sym1  = 1
      endif
#     ifdef switch_off_symmetry_coulomb
      nsym1 = 1
      sym1  = 1
#     endif

      ! Define reduced lists of G points -> pgptm1(:,ikpt), ikpt=1,..,nkpt
      allocate ( iarr(maxgptm),
     &           pointer(minval(gptm(1,:)):maxval(gptm(1,:)),
     &                   minval(gptm(2,:)):maxval(gptm(2,:)),
     &                   minval(gptm(3,:)):maxval(gptm(3,:))) )
      pgptm1 = 0
      do k = 1,nkpt1 ; ikpt = kpt1(k)
        pointer = 0
        do igpt = 1,ngptm(ikpt)
          g                       = gptm(:,pgptm(igpt,ikpt))
          pointer(g(1),g(2),g(3)) = igpt
        enddo
        iarr = 0
        j    = 0
        do igpt = ngptm(ikpt),1,-1
          if(iarr(igpt)==0) then
            j           = j + 1
            pgptm1(j,k) = igpt
            do isym1 = 1,nsym1(k)
              g = matmul ( sym(sym1(isym1,k))%rrot , gptm(:,pgptm(igpt,ikpt)) )
              i = pointer(g(1),g(2),g(3))
              if(i==0) Bug('zero pointer')
              iarr(i) = 1
            enddo
          endif
        enddo
        ngptm1(k) = j
      enddo
      deallocate ( iarr,pointer )

      ! Define imgl(l) = img**l
      imgl(0) = 1
      do i = 1,lexp
        imgl(i) = imgl(i-1) * img
      enddo

c
c     Calculate the structure constant
      Nallocate ( structconst,(S_ (2*lexp+1)**2,ncent,ncent,nkpt1 S_) )
      call structureconstant(Win(structconst),kpt(:,kpt1),nkpt1,lwrite)

      if(lwrite) then ; Rwrite(6,'(//A)') '### subroutine: coulombmatrix ###'
      else            ; Rwrite(6,'(A)') 'Calculate Coulomb matrix ... '
      endif

c
c     Calculate the Fourier coefficients (step) of the step function
      if(gcutstep>0) then
        Rwrite(6,'(/A)') 'Calculation of step function'
        call stepfunction(step,gptstep,ngptstep)
        i = maxval(abs(gptstep(1,:)))
        l = maxval(abs(gptstep(2,:)))
        m = maxval(abs(gptstep(3,:)))
        allocate ( pointer(-i:i,-l:l,-m:m) )
        do i=1,ngptstep
          pointer(gptstep(1,i),gptstep(2,i),gptstep(3,i)) = i
        enddo
      else
        ngptstep = 0
      endif

c
c     Testing of Coulomb-matrix calculation / redefine product basis => spherical Bessel functions
      if(tstcoul) then
        write(6,'(/A)') 'Product basis functions are replaced by orthonormalized sph. Bessel functions.'
        write(6,'(2X,A)') 'Reduction due to overlap:'
        call getnorm(kpt(:,2),1,gptm,ngptm(2),pgptm(:,2),size(pgptm,1),qnrm,nqnrm,pqnrm)
        deallocate (basm)
        allocate   (basm(maxgrid,nqnrm,0:maxlcutm,ntype))
        nindxm = nqnrm
        nbasp  = 0
        do itype=1,ntype
          if(ntype>1) write(6,'(4X,A,I3)') 'Atom type',itype
          do iqnrm=1,nqnrm
            do i=1,grid(itype)%number
              rdum = rgrid(i,itype)*qnrm(iqnrm)
              call sphbessel(rarr,rdum,lcutm(itype))
              basm(i,iqnrm,:,itype) = rarr(:maxlcutm) * rgrid(i,itype)
            enddo
          enddo
          do l=0,lcutm(itype)
            write(6,'(4X,A'NoA) lchar(l)//':'
            call orthonormalize(basm(:,:,l,itype),nindxm(l,itype),itype)
            nbasp = nbasp + neq(itype) * nindxm(l,itype) * (2*l+1)
          enddo
        enddo
        maxbasm = nbasp + maxgptm
        do k = 1,nkpt1 ; ikpt = kpt1(k)
          nbasm(ikpt) = nbasp + ngptm(ikpt)
        enddo
        deallocate(qnrm,pqnrm)
      endif

      Rif(lwrite) then
        write(6,'(/A,F7.1," MB")') 'Size of Coulomb matrix:',MBYTES/megabyte*size(coulomb)
        write(6,'(/A'NoA) 'Preparations...'
      endif
      call cpu_time(time1)

      ! Define gmat (symmetric)
      lm1 = 0
      do l1 = 0,lexp
        do m1 = -l1,l1
          lm1 = lm1 + 1
          lm2 = 0
      lp1:do l2 = 0,l1
            do m2 = -l2,l2
              lm2 = lm2 + 1
              if(lm2>lm1) exit lp1 ! Don't cross the diagonal!
              gmat(lm1,lm2) = sfac(l1+l2+m2-m1) * sfac(l1+l2+m1-m2)/
     &                        ( sfac(l1+m1) * sfac(l1-m1) * sfac(l2+m2) * sfac(l2-m2) ) /
     &                        sqrt( 1d0*(2*l1+1)*(2*l2+1)*(2*(l1+l2)+1) ) * (4*pi)**1.5d0
              gmat(lm2,lm1) = gmat(lm1,lm2)
            enddo
          enddo lp1
        enddo
      enddo

      ! Calculate moments of MT functions
      do itype = 1,ntype
        do l = 0,lcutm(itype)
          do i = 1,nindxm(l,itype)
            moment(i,l,itype) = intgrf(rgrid(:,itype)**(l+1)*basm(:,i,l,itype),itype) ! note that basm already contains the factor rgrid
          enddo
        enddo
        do i = 1,nindxm(0,itype)
          moment2(i,itype) = intgrf(rgrid(:,itype)**3*basm(:,i,0,itype),itype)
        enddo
      enddo

      ! Look for different qnorm = |k+G|, definition of qnrm and pqnrm.
      call getnorm(kpt(:,kpt1),nkpt1,gptm,ngptm(kpt1),pgptm(:,kpt1),size(pgptm,1),qnrm,nqnrm,pqnrm)

      ! Calculate moments of spherical Bessel functions (for (2) and (3))              (->sphbesmoment)
      ! Calculate overlap of spherical Bessel functions with basis functions (for (2)) (->olap)
      ! Calculate overlap of sphbesmoment1(r,l)         with basis functions (for (2)) (->integral)
      ! We use           sphbes(r,l) = j_l(qr)
      ! and       sphbesmoment1(r,l) = 1/r**(l-1) * INT(0..r) r'**(l+2) * j_l(qr') dr'
      !                                + r**(l+2) * INT(r..S) r'**(1-l) * j_l(qr') dr' .
      allocate ( sphbesmoment(0:lexp,ntype,nqnrm),
     &           olap(maxindxm,0:maxlcutm,ntype,nqnrm),
     &           integral(maxindxm,0:maxlcutm,ntype,nqnrm) )
      sphbes        = 0
      sphbesmoment  = 0
      sphbesmoment1 = 0

      do iqnrm = 1,nqnrm
        qnorm = qnrm(iqnrm)
        do itype = 1,ntype
          ng            = grid(itype)%number
          rdum          = grid(itype)%radius
          sphbes        = 0
          sphbesmoment1 = 0
          if(qnorm==0) then
            sphbesmoment(0,itype,iqnrm) = rdum**3 / 3
            do i = 1,ng
              sphbes(i,0)        = 1
              sphbesmoment1(i,0) = rgrid(i,itype)**2 / 3 + ( grid(itype)%radius**2 - rgrid(i,itype)**2 ) / 2
            enddo
          else
            call sphbessel(rarr,qnorm*rdum,lexp+1)
            do l = 0,lexp
              sphbesmoment(l,itype,iqnrm) = rdum**(l+2) * rarr(l+1) / qnorm
            enddo
            do i = ng,1,-1
              rdum = rgrid(i,itype)
              call sphbessel(rarr,qnorm*rdum,lcutm(itype)+1)
              do l = 0,lcutm(itype)
                sphbes(i,l) = rarr(l)
                if(l/=0) then ; rdum1 = -rdum**(1-l) * rarr(l-1)
                else          ; rdum1 = -cos(qnorm*rdum) / qnorm
                endif
                if(i==ng) rarr1(l) = rdum1
                sphbesmoment1(i,l) = ( rdum**(l+2) * rarr(l+1) / rdum**(l+1) + ( rarr1(l) - rdum1 ) * rdum**l ) / qnorm
              enddo
            enddo
          endif
          do l = 0,lcutm(itype)
            do n = 1,nindxm(l,itype)
              olap    (n,l,itype,iqnrm) = intgrf( rgrid(:,itype) * basm(:,n,l,itype) * sphbes(:,l),        itype ) ! note that basm already contains one factor rgrid
              integral(n,l,itype,iqnrm) = intgrf( rgrid(:,itype) * basm(:,n,l,itype) * sphbesmoment1(:,l), itype ) !
            enddo
          enddo
        enddo
      enddo

      call cpu_time(time2)
      if(lwrite) then ; Rwrite(6,'(A,18X,A,F8.2,A)') ' done','( Timing:',time2-time1,' )'
      else            ; Rwrite(6,'(''    Timings:'',F7.2'NoA) time2-time1
      endif
      call cpu_time(time1)

c 123  continue
      Nfence(coulomb1)
      ifO coulomb1 = 0
      Nfence(coulomb1)

c
c     (1) Case < MT | v | MT >

      Rif(lwrite) write(6,'(A'NoA) '< MT | v | MT > contribution...'

c     (1a) r,r' in same MT

      Mpi( mm = -1 )
      ix  = 0
      iy  = 0
      iy0 = 0
      do itype = 1,ntype
        do ieq = 1,neq(itype) ! Here the diagonal block matrices do not depend on ieq. In (1b) they do depend on ieq, though,
          ifMODP(mm)
          do l = 0,lcutm(itype)
            do n2 = 1,nindxm(l,itype)
              call primitivef(primf1,basm(:,n2,l,itype)*rgrid(:,itype)**(l+1), itype) ! note that basm already contains the factor rgrid
              call primitivef(primf2,basm(:,n2,l,itype)/rgrid(:,itype)**l    ,-itype) ! -itype is to enforce inward integration
              primf1 = primf1 / rgrid(:,itype)**l
              primf2 = primf2 * rgrid(:,itype)**(l+1)
              do n1 = 1,n2
                integrand           = basm(:,n1,l,itype) * (primf1 + primf2)
                mat(n2*(n2-1)/2+n1) = (4*pi)/(2*l+1) * intgrf(integrand,itype)
              enddo
            enddo
            ! distribute mat for m=-l,l on coulomb in block-matrix form
            do m = -l,l
              do n2 = 1,nindxm(l,itype)
                ix = ix + 1
                iy = iy0
                do n1 = 1,n2
                  iy            = iy + 1
                  i             = ix*(ix-1)/2 + iy
                  j             = n2*(n2-1)/2 + n1
                  coulomb1(i,:) = mat(j)
                enddo
              enddo
              iy0 = ix
            enddo
          enddo
          elseMOD( ix = ix + sum( [ ((2*l+1)*nindxm(l,itype),l=0,lcutm(itype)) ] ) ; iy0 = ix )
        enddo
      enddo
      Nfence(coulomb1)

c     (1b) r,r' in different MT

      Nallocate ( coulmat,(S_ nbasp,nbasp S_) )

      do k = 1,nkpt1 ; ikpt = kpt1(k)

      Nfence(coulmat)
      ifO coulmat = 0
      Nfence(coulmat)

      Mpi( mm = -1 )
      ix  = 0
      ic2 = 0
      do itype2 = 1,ntype
        !nn2 = sum( [ ((2*l+1)*nindxm(l,itype2),l=0,lcutm(itype2)) ] )
        do ieq2 = 1,neq(itype2)
          ic2 = ic2 + 1
          lm2 = 0
          do l2 = 0,lcutm(itype2)
            do m2 = -l2,l2
              lm2 = lm2 + 1
              do n2 = 1,nindxm(l2,itype2)
                ix  = ix + 1 !; Inv( ix_ = ix + (pcent(ic2,invsym)-ic2)*nn2 - 2*m2*nindxm(l2,itype2) )

                iy  = 0
                ic1 = 0
           lp2: do itype1 = 1,itype2
                  !nn1 = sum( [ ((2*l+1)*nindxm(l,itype1),l=0,lcutm(itype1)) ] )
                  do ieq1 = 1,neq(itype1)
                    ic1 = ic1 + 1
                    lm1 = 0
                    do l1 = 0,lcutm(itype1)
                      do m1 = -l1,l1
                        lm1   = lm1 + 1
                        do n1 = 1,nindxm(l1,itype1)
                          iy   = iy + 1 ; if(iy>ix) exit lp2 ; McycleP(mm)
                          l    = l1 + l2
                          lm   = l**2 + l + m1 - m2 + 1
                          rdum = (-1)**(l2+m2)*moment(n1,l1,itype1)*moment(n2,l2,itype2)*gmat(lm1,lm2)
                          cdum = exp(img* 2*pi *dot_product(kpt(:,ikpt),cent(:,ic2)-cent(:,ic1))) * rdum * structconst(lm,ic1,ic2,k)
                          coulmat(iy,ix) = cdum
                          coulmat(ix,iy) = conjg(cdum)
                        enddo
                      enddo
                    enddo
                  enddo
                enddo lp2

              enddo
            enddo
          enddo
        enddo
      enddo

      Nfence(coulmat)
      Inv( Ocall symmetrize(coulmat,nbasp,nbasp,3) )
      ifO coulomb1(:nbasp*(nbasp+1)/2,k) = coulomb1(:nbasp*(nbasp+1)/2,k) + packmat(coulmat)
      Nfence(coulomb1)

      enddo ! kpoint loop

      Ndeallocate ( coulmat )

      call cpu_time(time2)
      if(lwrite) then ; Rwrite(6,'(A,2X,A,F8.2,A)') ' done','( Timing:',time2-time1,' )'
      else            ; Rwrite(6,'(F7.2'NoA) time2-time1
      endif
      call cpu_time(time1)

c
c     (2) Case < MT | v | PW >

      Rif(lwrite) write(6,'(A'NoA) '< MT | v | PW > contribution...'

c     (2a) r in MT, r' everywhere
c     (2b) r,r' in same MT
c     (2c) r,r' in different MT

      do k = 1,nkpt1 ; ikpt = kpt1(k)

      ! start to loop over interstitial plane waves
      allocate ( carr(nbasp) )
      do igpt0 = 1,ngptm1(k) ; Mcycle(igpt0)
        carr  = 0
        igpt  = pgptm1(igpt0,k)
        igptp = pgptm(igpt,ikpt)
        ix    = nbasp + igpt
        q     = matmul ( rlat, kpt(:,ikpt) + gptm(:,igptp) )
        qnorm = sqrt(sum(q**2))
        iqnrm = pqnrm(igpt,k)
        if(abs(qnrm(iqnrm)-qnorm)>1d-12)
     &    Bug('qnorm does not equal corresponding element in qnrm')

        call harmonicsr(y,q,lexp)
        y  = conjg(y)
        iy = 0
        ic = 0
        do itype = 1,ntype
          do ieq = 1,neq(itype)
            ic = ic + 1
            lm = 0
            do l = 0,lcutm(itype)
              do m = -l,l
                lm = lm + 1

                ! calculate sum over lm and centers for (2c) -> csum, csumf
                csum  = 0
                csumf = 0
                ic1   = 0
                do itype1 = 1,ntype
                  do ieq1 = 1,neq(itype1)
                    ic1  = ic1 + 1
                    cexp = 4*pi * exp( img * 2*pi * ( dot_product( kpt(:,ikpt)+gptm(:,igptp), cent(:,ic1) ) -
     &                                                dot_product( kpt(:,ikpt)              , cent(:,ic) ) ) )

                    lm1 = 0
                    do l1 = 0,lexp
                      l2   = l + l1 ! for structconst
                      idum = 1
                      cdum = sphbesmoment(l1,itype1,iqnrm) * imgl(l1) * cexp
                      do m1 = -l1,l1
                        lm1  = lm1 + 1
                        m2   = m - m1              ! for structconst
                        lm2  = l2**2 + l2 + m2 + 1 !
                        csum = csum - idum * gmat(lm1,lm) * y(lm1) * cdum * structconst(lm2,ic,ic1,k)
                        idum = -idum ! factor (-1)*(l1+m1)
                      enddo
                    enddo

                    ! add contribution of (2c) to csum and csumf coming from linear and quadratic orders of Y_lm*(G) / G * j_(l+1)(GS)
                    if(ikpt==1.and.l<=2) then
                      cexp      = exp(img*2*pi*dot_product(gptm(:,igptp),cent(:,ic1))) * gmat(lm,1) * 4*pi/vol
                      csumf(lm) = csumf(lm) - cexp * sqrt(4*pi) * img**l * sphbesmoment(0,itype1,iqnrm) / facfac(l-1)
                      if(l==0) then
                        if(igpt/=1) then
                          csum = csum - cexp * ( sphbesmoment(0,itype1,iqnrm)*grid(itype1)%radius**2 -
     &                                           sphbesmoment(2,itype1,iqnrm)*2d0/3 ) / 10
                        else
                          csum = csum - cexp * grid(itype1)%radius**5/30!10 / 3
                        endif
                      else if(l==1) then
                        csum = csum + cexp * img * sqrt(4*pi) * sphbesmoment(1,itype1,iqnrm) * y(lm) / 3
                      endif
                    endif

                  enddo
                enddo

                ! add contribution of (2a) to csumf
                if(ikpt==1.and.igpt==1.and.l<=2) then
                  csumf(lm) = csumf(lm) + (4*pi)**2 * img**l / facfac(l)
                endif

                ! finally calculate coulomb elements -> carr
                cdum = (4*pi)**2 * imgl(l) * y(lm) * exp(img * 2*pi * dot_product(gptm(:,igptp),cent(:,ic)))
                do n = 1,nindxm(l,itype)
                  iy = iy + 1
                  if(ikpt==1.and.igpt==1) then
                    if(l==0)
     &                carr(iy) =   - cdum * moment2(n,itype) / 6 / svol         ! (2a)
                    carr(iy)   = carr(iy) +
     &                           ( - cdum / (2*l+1) * integral(n,l,itype,iqnrm) ! (2b)
     &                             + csum * moment(n,l,itype) ) / svol          ! (2c)
                  else
                    carr(iy)   = (   cdum * olap(n,l,itype,iqnrm) / qnorm**2    ! (2a)
     &                             - cdum / (2*l+1) * integral(n,l,itype,iqnrm) ! (2b)
     &                             + csum * moment(n,l,itype) ) / svol          ! (2c)
                  endif
                enddo

              enddo
            enddo
          enddo
        enddo
        Inv( call symmetrize(carr,nbasp,1,1) )
        idum                          = ix*(ix-1)/2
        coulomb1(idum+1:idum+nbasp,k) = carr
      enddo
      deallocate ( carr )

      enddo ! kpoint loop
      Nfence(coulomb1)

      call cpu_time(time2)
      if(lwrite) then ; Rwrite(6,'(A,2X,A,F8.2,A)') ' done','( Timing:',time2-time1,' )'
      else            ; Rwrite(6,'(F7.2'NoA) time2-time1
      endif
      call cpu_time(time1)

c     If gcutstep > 0, we test the MT-PW contribution with the step function (otherwise we proceed with case (3)).
c     For this we have to calculate all <MT|v|PW> for all plane waves defined by gcutstep.
      if(gcutstep>0d0) then
        allocate ( coulstep(nbasp*ngptstep) )
        coulstep = 0
        ! Look for different qnorm = |k+G|, definition of qnrmstep and pqnrmstep.
        call cpu_time(time1)
        nqnrmstep = 0
        call getnorm(kpt(:,kpt1),nkpt1,gptstep,[(ngptstep,i=1,nkpt1)],[((i,i=1,ngptstep),j=1,nkpt1)],ngptstep,
     &               qnrmstep,nqnrmstep,pqnrmstep)
        allocate ( olapstep(maxindxm,0:maxlcutm,ntype,nqnrmstep),sphbes0(maxgrid,0:maxlcutm,1) )
        do iqnrm = 1,nqnrmstep
          qnorm = qnrmstep(iqnrm)
          if(qnorm/=0) then
            do itype=1,ntype
              do i=1,grid(itype)%number
                rdum = rgrid(i,itype)*qnorm
                call sphbessel(sphbes0(i,:,1),rdum,lcutm(itype))
              enddo
              do l=0,lcutm(itype)
                do n=1,nindxm(l,itype)
                  olapstep(n,l,itype,iqnrm) = intgrf ( rgrid(:,itype) * basm(:,n,l,itype) * sphbes0(:,l,1), itype ) ! note that basm already contains one factor rgrid
                enddo
              enddo
            enddo
          endif
        enddo
        deallocate ( sphbes0 )
        rdum1 = 0
        do k = 1,nkpt1 ; ikpt = kpt1(k) ; if(ikpt==1) cycle
          ix = 0
          do igpt = 1,ngptstep
            ix    = ix + 1
            q     = matmul ( rlat, kpt(:,ikpt) + gptstep(:,igpt) )
            qnorm = sqrt(sum(q**2))
            iqnrm = pqnrmstep(igpt,k)
            if(abs(qnrmstep(iqnrm)-qnorm)>1d-12)
     &        Bug('qnorm does not equal corresponding element in qnrmstep')
            rdum  = (4*pi)**2 / qnorm**2
            call harmonicsr(y,q,lexp)
            y    = conjg(y)
            idum = (ix-1)*nbasp
            iy   = 0
            ic   = 0
            do itype = 1,ntype

              do ieq = 1,neq(itype)
                ic   = ic + 1
                cexp = exp(img * 2*pi * dot_product(gptstep(:,igpt),cent(:,ic)))
                lm   = 0
                do l = 0,lcutm(itype)
                  do m = -l,l
                    lm   = lm + 1

                    ! finally define coulstep
                    cdum = imgl(l) * y(lm) * cexp
                    do n=1,nindxm(l,itype)
                      iy                = iy + 1
                      coulstep(idum+iy) = cdum * rdum * olapstep(n,l,itype,iqnrm) / svol
                    enddo

                  enddo
                enddo
              enddo

            enddo
          enddo

          Inv( if(ngptstep>0) call symmetrize(coulstep,nbasp,ngptstep,1) )

          do igpt = 1,ngptm(ikpt)
          do n = 1,nbasp

c          call random_number(rdum)                !
c          igpt  = nint((ngptm(ikpt)-1)*rdum) + 1  ! We only calculate the deviation for random matrix elements.
c          igptp = pgptm(max(1,igpt),ikpt)         !
c          call random_number(rdum)                ! (Otherwise it would take too long.)
c          n     = max(1,nint((nbasp-1)*rdum) + 1) !

c          igpt = ngptm(ikpt)
c          n    = nbasp

          igptp = pgptm(max(1,igpt),ikpt)

          cdum = 0
          do i=1,ngptstep
            g = gptstep(:,i) + gptm(:,igptp)
            if(gptnorm(g)<=gcutstep) then
              m    = pointer(g(1),g(2),g(3))
              cdum = cdum + step(i) * coulstep((m-1)*nbasp+n)
            endif
          enddo
          rdum1 = rdum1 + abs(cdum-coulomb1((nbasp+igpt)*((nbasp+igpt)-1)/2+n,k))**2

          coulomb1((nbasp+igpt)*((nbasp+igpt)-1)/2+n,k) = cdum
          enddo
          enddo

        enddo

        if(any(kpt1/=1)) then
          call cpu_time(time2)
          Rwrite(6,'(A,F15.10'NoA) '  Step-function test:',sqrt(rdum1/count(kpt1/=1))
          Rwrite(6,'(A,F8.2,A )')  '  ( Timing:',time2-time1,' )'
          call cpu_time(time1)
        endif
      endif

c
c     (3) Case < PW | v | PW >

      Rif(lwrite) write(6,'(A'NoA) '< PW | v | PW > contribution...'

      allocate ( sphbes0(-1:lexp+2,ntype,nqnrm) )
      do iqnrm = 1,nqnrm
        do itype = 1,ntype
          rdum = qnrm(iqnrm) * grid(itype)%radius
          call sphbessel(sphbes0(0,itype,iqnrm),rdum,lexp+2)
          if(rdum/=0) sphbes0(-1,itype,iqnrm) = cos(rdum)/rdum
        enddo
      enddo

      do k = 1,nkpt1 ; ikpt = kpt1(k)

c     (3a) r,r' everywhere; r everywhere, r' in MT; r in MT, r' everywhere

      do igpt0 = 1,ngptm1(k) ; Mcycle(igpt0)
        igpt2  = pgptm1(igpt0,k)
        igptp2 = pgptm(igpt2,ikpt)
        ix     = nbasp + igpt2
        iy     = nbasp
        q2     = matmul ( rlat, kpt(:,ikpt) + gptm(:,igptp2) )
        rdum2  = sum(q2**2)
        if(rdum2/=0) rdum2 = 4*pi / rdum2
        do igpt1 = 1,igpt2
          igptp1 = pgptm(igpt1,ikpt)
          ! Calculate cint = sum(a) integral(MT(a)) exp[i(Gj-Gi)r] dr
          cint   = 0
          g      = gptm(:,igptp2)-gptm(:,igptp1)
          gnorm  = sqrt(sum(matmul(rlat,g)**2))
          if(gnorm==0) then
            do itype = 1,ntype
              cint = cint + neq(itype) * 4*pi*grid(itype)%radius**3/3
            enddo
          else
            ic = 0
            do itype = 1,ntype
              rdum = grid(itype)%radius * gnorm
              rdum = 4*pi * ( sin(rdum) - rdum * cos(rdum) ) / gnorm**3
              do ieq = 1,neq(itype)
                ic   = ic + 1
                cint = cint + rdum * exp( img * 2*pi * dot_product(cent(:,ic),g) )
              enddo
            enddo
          endif
          ! Update Coulomb matrix
          iy     = iy + 1
          q1     = matmul ( rlat, kpt(:,ikpt) + gptm(:,igptp1) )
          idum   = ix*(ix-1)/2+iy
          rdum1  = sum(q1**2)
          if(rdum1/=0) rdum1 = 4*pi / rdum1
          if(ikpt==1) then
            if(igpt1/=1) then
              coulomb1(idum,k) =                  - cint * rdum1 / vol
            endif
            if(igpt2/=1) then
              coulomb1(idum,k) = coulomb1(idum,k) - cint * rdum2 / vol
            endif
          else
            coulomb1(idum,k) =                    - cint * ( rdum1 + rdum2 ) / vol
          endif
        enddo
        if(ikpt/=1.or.igpt2/=1) then              !
          coulomb1(idum,k) = coulomb1(idum,k) + rdum2 ! diagonal term
        endif                                         !
      enddo
      Nfence(coulomb1)

c     (3b) r,r' in different MT

      ! group together quantities which depend only on l,m and igpt -> carr2y
      allocate ( carr2a((lexp+1)**2,maxgptm),carr2b(ncent,maxgptm) )
      do igpt=1,ngptm(ikpt)
        igptp = pgptm(igpt,ikpt)
        iqnrm = pqnrm(igpt,k)
        q     = matmul ( rlat, kpt(:,ikpt) + gptm(:,igptp) )
        call harmonicsr(y,q,lexp)
        y     = conjg(y)
        lm    = 0
        do l = 0,lexp
          do m = -l,l
            lm              = lm + 1
            carr2a(lm,igpt) = 4*pi * imgl(l) * y(lm)
          enddo
        enddo
        do ic = 1,ncent
          carr2b(ic,igpt) = exp ( -img * 2*pi * dot_product(kpt(:,ikpt)+gptm(:,igptp),cent(:,ic)) )
        enddo
      enddo

      ! finally we can loop over the plane waves (G: igpt1,igpt2)
      allocate ( carr2(ncent,(lexp+1)**2),structconst1(ncent,(2*lexp+1)**2) )
      do igpt0 = 1,ngptm1(k) ; Mcycle(igpt0)
        igpt2  = pgptm1(igpt0,k)
        ix     = nbasp + igpt2
        igptp2 = pgptm(igpt2,ikpt)
        iqnrm2 = pqnrm(igpt2,k)
        ic2    = 0
        carr2  = 0
        do itype2 = 1,ntype
          do ieq2 = 1,neq(itype2)
            ic2   = ic2 + 1
            cexp  = conjg ( carr2b(ic2,igpt2) )
            lm2   = 0
            do ic1 = 1,ncent
              structconst1(ic1,:) = structconst(:,ic1,ic2,k)
            enddo
            do l2 = 0,lexp
              idum = 1
              do m2 = -l2,l2
                lm2  = lm2 + 1
                cdum = idum * sphbesmoment(l2,itype2,iqnrm2) * cexp * carr2a(lm2,igpt2)
                if(cdum/=0) then
                  lm1 = 0
                  do l1 = 0,lexp
                    l  =  l1 + l2
                    m  = -l1 - m2 ! first loop of m1
                    lm = l**2 + l + m
                    do m1 = -l1,l1
                      lm1   = lm1 + 1
                      lm    = lm  + 1
                      cdum1 = cdum * gmat(lm1,lm2)
                      do ic1 = 1,ncent
                        carr2(ic1,lm1) = carr2(ic1,lm1) + cdum1 * structconst1(ic1,lm)
                      enddo
                    enddo
                  enddo
                endif
                idum = -idum ! factor (-1)**(l+m)
              enddo
            enddo
          enddo
        enddo
        iy = nbasp
        do igpt1 = 1,igpt2
          iy      = iy + 1
          igptp1  = pgptm(igpt1,ikpt)
          iqnrm1  = pqnrm(igpt1,k)
          csum    = 0
          ic      = 0
          do itype = 1,ntype
            do ieq = 1,neq(itype)
              ic   = ic  + 1
              cexp = carr2b(ic,igpt1)
              lm   = 0
              do l = 0,lexp
                cdum = cexp * sphbesmoment(l,itype,iqnrm1)
                do m = -l,l
                  lm   = lm + 1
                  csum = csum + cdum * carr2(ic,lm) * conjg ( carr2a(lm,igpt1) ) ! for coulomb
                enddo
              enddo
            enddo
          enddo
          idum             = ix*(ix-1)/2+iy
          coulomb1(idum,k) = coulomb1(idum,k) + csum / vol
        enddo
      enddo
      Nfence(coulomb1)

      deallocate ( carr2,carr2a,carr2b,structconst1 )

c     Add corrections from higher orders in (3b) to coulomb1(:,1)
      ! (1) igpt1 > 1 , igpt2 > 1  (finite G vectors)
      if(ikpt==1) then
        rdum = (4*pi)**(1.5d0)/vol**2 * gmat(1,1)
        do igpt0 = 1,ngptm1(k) ; Mcycle(igpt0)
          igpt2  = pgptm1(igpt0,k) ; if(igpt2==1) cycle
          ix     = nbasp + igpt2
          iqnrm2 = pqnrm(igpt2,k)
          igptp2 = pgptm(igpt2,1)
          q2     = matmul(rlat,gptm(:,igptp2))
          qnorm2 = sqrt(sum(q2**2))
          iy     = nbasp + 1
          do igpt1 = 2,igpt2
            iy     = iy + 1
            idum   = ix*(ix-1)/2+iy
            iqnrm1 = pqnrm(igpt1,k)
            igptp1 = pgptm(igpt1,1)
            q1     = matmul(rlat,gptm(:,igptp1))
            qnorm1 = sqrt(sum(q1**2))
            rdum1  = dot_product(q1,q2) / (qnorm1*qnorm2)
            ic1    = 0
            do itype1 = 1,ntype
              do ieq1 = 1,neq(itype1)
                ic1 = ic1 + 1
                ic2 = 0
                do itype2 = 1,ntype
                  do ieq2 = 1,neq(itype2)
                    ic2  = ic2 + 1
                    cdum = exp ( img * 2*pi * ( - dot_product(gptm(:,igptp1),cent(:,ic1)) +
     &                                            dot_product(gptm(:,igptp2),cent(:,ic2)) ) )
                    coulomb1(idum,k) = coulomb1(idum,k) + rdum * cdum * (
     &                - sphbesmoment(1,itype1,iqnrm1) * sphbesmoment(1,itype2,iqnrm2) * rdum1  / 3
     &                - sphbesmoment(0,itype1,iqnrm1) * sphbesmoment(2,itype2,iqnrm2)          / 6
     &                - sphbesmoment(2,itype1,iqnrm1) * sphbesmoment(0,itype2,iqnrm2)          / 6
     &                + sphbesmoment(0,itype1,iqnrm1) * sphbesmoment(1,itype2,iqnrm2) / qnorm2 / 2
     &                + sphbesmoment(1,itype1,iqnrm1) * sphbesmoment(0,itype2,iqnrm2) / qnorm1 / 2 )
                  enddo
                enddo
              enddo
            enddo
          enddo
        enddo
        Nfence(coulomb1)
        Rbegin
        ! (2) igpt1 = 1 , igpt2 > 1  (first G vector vanishes, second finite)
        iy = nbasp + 1
        do igpt0 = 1,ngptm1(k)
          igpt2  = pgptm1(igpt0,k) ; if(igpt2==1) cycle
          ix     = nbasp + igpt2
          iqnrm2 = pqnrm(igpt2,k)
          igptp2 = pgptm(igpt2,1)
          qnorm2 = qnrm(iqnrm2)
          idum   = ix*(ix-1)/2+iy
          do itype1 = 1,ntype
            do ieq1 = 1,neq(itype1)
              ic2 = 0
              do itype2 = 1,ntype
                do ieq2 = 1,neq(itype2)
                  ic2              = ic2 + 1
                  cdum             = exp ( img * 2*pi * dot_product(gptm(:,igptp2),cent(:,ic2)) )
                  coulomb1(idum,k) = coulomb1(idum,k) + rdum * cdum * grid(itype1)%radius**3 * (
     &              + sphbesmoment(0,itype2,iqnrm2) / 30 * grid(itype1)%radius**2
     &              - sphbesmoment(2,itype2,iqnrm2) / 18
     &              + sphbesmoment(1,itype2,iqnrm2) /  6 / qnorm2 )
                enddo
              enddo
            enddo
          enddo
        enddo
        ! (2) igpt1 = 1 , igpt2 = 1  (vanishing G vectors)
        iy   = nbasp + 1
        ix   = nbasp + 1
        idum = ix*(ix-1)/2+iy
        do itype1 = 1,ntype
          do ieq1 = 1,neq(itype1)
            do itype2 = 1,ntype
              do ieq2 = 1,neq(itype2)
                coulomb1(idum,k) = coulomb1(idum,k) + rdum * grid(itype1)%radius**3 * grid(itype2)%radius**3 *
     &            ( grid(itype1)%radius**2 + grid(itype2)%radius**2 ) / 90
              enddo
            enddo
          enddo
        enddo
        Rend
        Nfence(coulomb1)
      endif

c     (3c) r,r' in same MT

      allocate ( carr2((lexp+1)**2,maxgptm) )
      do igpt = 1,ngptm(ikpt)
        igptp = pgptm(igpt,ikpt)
        q     = matmul ( rlat, kpt(:,ikpt) + gptm(:,igptp) )
        call harmonicsr(carr2(:,igpt),q,lexp)
      enddo

      do igpt0 = 1,ngptm1(k) ; Mcycle(igpt0)
        igpt2          = pgptm1(igpt0,k)
        ix             = nbasp + igpt2
        igptp2         = pgptm(igpt2,ikpt)
        iqnrm2         = pqnrm(igpt2,k)
        q2             = matmul ( rlat, kpt(:,ikpt) + gptm(:,igptp2) )
        y2             = conjg ( carr2(:,igpt2) )
        iy             = nbasp
        do igpt1 = 1,igpt2
          iy     = iy + 1
          igptp1 = pgptm(igpt1,ikpt)
          iqnrm1 = pqnrm(igpt1,k)
          q1     = matmul ( rlat, kpt(:,ikpt) + gptm(:,igptp1) )
          y1     = carr2(:,igpt1)
          cexp1  = 0
          ic     = 0
          do itype = 1,ntype
            do ieq = 1,neq(itype)
              ic           = ic + 1
              cexp1(itype) = cexp1(itype) +
     &                       exp(img * 2*pi * dot_product((gptm(:,igptp2)-gptm(:,igptp1)),cent(:,ic)) )
            enddo
          enddo
          lm   = 0
          cdum = 0
          do l = 0,lexp
            cdum1 = 0
            do itype = 1,ntype
              cdum1 = cdum1 + cexp1(itype)*sphbessel_integral() / (2*l+1)
            enddo
            do m = -l,l
              lm   = lm + 1
              cdum = cdum + cdum1 * y1(lm) * y2(lm)
            enddo
          enddo
          idum             = ix*(ix-1)/2+iy
          coulomb1(idum,k) = coulomb1(idum,k) + (4*pi)**3 * cdum / vol
        enddo
      enddo
      Nfence(coulomb1)

      deallocate ( carr2 )

      enddo ! kpoint loop

      deallocate ( sphbesmoment,sphbes0 )
      Ndeallocate ( structconst )

      call cpu_time(time2)
      if(lwrite) then ; Rwrite(6,'(A,2X,A,F8.2,A)') ' done','( Timing:',time2-time1,' )'
      else            ; Rwrite(6,'(F8.2'NoA) time2-time1
      endif
      call cpu_time(time1)

      Obegin
      Mpi( call Msum(coulomb1,comm=Ocomm) )

c
c     Symmetry-equivalent G vectors
c       <M|v|G1> = <M|v|PG2>   * exp[i*(k+G1)a] | P=(A,a,C) rotates from G1 to G2 (with phase factor)
c       = [ C <P^-1 M|v| G2> ] * exp[i*(k+G1)a] | v is invariant wrt (A,a)
      Rif(lwrite) write(6,'(A'NoA) 'Symm.-equiv. matrix elements...'

      do k = 1,nkpt1 ; ikpt = kpt1(k)

      allocate ( carr2(nbasm(ikpt),2),iarr(maxgptm) )
      iarr                       = 0
      iarr(pgptm1(:ngptm1(k),k)) = 1
      do igpt2 = ngptm(ikpt),1,-1 ; if(iarr(igpt2)/=1) Bug('G column not defined.')
        j = (nbasp+igpt2-1) * (nbasp+igpt2) / 2
        do i = 1,nbasp+igpt2
          j          = j + 1
          carr2(i,2) = coulomb1(j,k)
        enddo
        do i = nbasp+igpt2+1,nbasm(ikpt)
          j          = j + i - 1
          carr2(i,2) = MCONJG(coulomb1(j,k))
        enddo
        do isym1 = 2,nsym1(k)
          isym  = sym1(isym1,k)
          g     = nint ( matmul(sym(isym)%rrot,kpt(:,ikpt))-kpt(:,ikpt) )
          g     = matmul(sym(isym)%rrot,gptm(:,pgptm(igpt2,ikpt))) + g
          igpt1 = pntgptm(g(1),g(2),g(3),ikpt) ! G1 = P*G2 * phase
          if(iarr(igpt1)==0) then
            call mtrafo(carr2(:,1),carr2(:,2),nbasm(ikpt),1,ikpt,isym,1,.false.) ! C <P^-1 M|v|G2>
            carr2(:,1) = carr2(:,1) * exp(img * 2*pi * dot_product(kpt(:,ikpt)+gptm(:,pgptm(igpt1,ikpt)),sym(isym)%transl)) ! phase
            l                             = (nbasp+igpt1-1) * (nbasp+igpt1) / 2
            coulomb1(l+1:l+nbasp+igpt1,k) = carr2(:nbasp+igpt1,1)
            iarr(igpt1)                   = 1
          endif
        enddo
      enddo
      deallocate ( carr2,iarr )

      enddo ! kpoint loop

      call cpu_time(time2)
      if(lwrite) then ; Rwrite(6,'(A,2X,A,F8.2,A)') ' done','( Timing:',time2-time1,' )'
      else            ; Rwrite(6,'(F7.2)') time2-time1
      endif
      call cpu_time(time1)

c     If gcutstep > 0, we test the PW-PW contribution with the step function.
      rdum1 = 0
      if(gcutstep>0) then
        do k = 1,nkpt1 ; ikpt = kpt1(k) ; if(ikpt==1) cycle

          do igpt2 = 1,ngptm(ikpt)
          do igpt1 = 1,igpt2

c          call random_number(rdum)             !
c          igpt1 = nint((ngptm(ikpt)-1)*rdum)+1 !
c          call random_number(rdum)             !
c          igpt2 = nint((ngptm(ikpt)-1)*rdum)+1 ! We only calculate the deviation for random matrix elements.
c          if(igpt1>igpt2) then                 ! (Otherwise it would take too long.)
c            idum  = igpt1                      !
c            igpt1 = igpt2                      !
c            igpt2 = idum                       !
c          endif                                !

c          igpt1 = ngptm(ikpt)
c          igpt2 = ngptm(ikpt)
          igptp1 = pgptm(igpt1,ikpt)
          igptp2 = pgptm(igpt2,ikpt)
          cdum   = 0
          do i=1,ngptstep
            g    = gptstep(:,i) - gptm(:,igptp1) ; if(gptnorm(g)>gcutstep) cycle
            m1   = pointer(g(1),g(2),g(3))       ; if(m1==0) cycle
            g    = gptstep(:,i) - gptm(:,igptp2) ; if(gptnorm(g)>gcutstep) cycle
            m2   = pointer(g(1),g(2),g(3))       ; if(m2==0) cycle
            q    = matmul(rlat,(kpt(:,ikpt)+gptstep(:,i)))
            cdum = cdum + conjg(step(m1)) * step(m2) / sum(q**2)
          enddo
          rdum1 = rdum1 + abs(4*pi*cdum-coulomb1((nbasp+igpt2)*(nbasp+igpt2-1)/2+nbasp+igpt1,k))**2

          coulomb1((nbasp+igpt2)*(nbasp+igpt2-1)/2+nbasp+igpt1,k) = 4*pi*cdum
          enddo
          enddo

        enddo

        if(any(kpt1/=1)) then
          call cpu_time(time2)
          Rwrite(6,'(A,F15.10'NoA) '  Step-function test:',sqrt(rdum1/count(kpt1/=1))
          Rwrite(6,'(A,F8.2,A )')  '  ( Timing:',time2-time1,' )'
          call cpu_time(time1)
        endif
      endif

      if(chkcoul) then

      if(maxgptm<=0)        then ; Warn('Coulomb check (CHKCOUL) impossible because maxgptm<=0.')
      else if(all(kpt1==1)) then ; Warn('Coulomb check (CHKCOUL) impossible for k=0.')
      else

c
c     FINAL TEST
c
c     Finally, we test the whole Coulomb matrix by representing plane waves |PW> in the mixed basis { |MB(i)> },
c     |PW'> = sum(i) c(i) |MB(i)>, and testing whether <PW'|v|PW'> = 4*pi/(k+G)**2

      write(6,'(/A)') 'Tests of Coulomb matrix'
      call cpu_time(time1)
      write(6,'(A)')  '  (1) The matrix <PW''|v|PW''> is calculated and compared to 4*pi/k**2 where',
     &                '      the |PW''> are plane waves represented in the (incomplete) mixed basis.'

      ! First, we have to represent the spherical Bessel functions in the MT product basis. => coefficients csphbes
      allocate ( csphbes(maxindxm,0:maxlcutm,ntype,nqnrm),error(0:maxlcutm,ntype),
     &           carr3(nbasp,maxgptm,nkpt1),carr1(nbasp),coulmat(nbasp,nbasp),
     &           work(2*maxbasm),rwork(3*maxbasm) )
      error = 0
      do iqnrm=1,nqnrm
        qnorm = qnrm(iqnrm)
        do itype=1,ntype
          ng = grid(itype)%number
          do i=1,ng
            rdum = rgrid(i,itype)*qnorm
            call sphbessel(sphbes(i,:),rdum,lcutm(itype))
          enddo
          do l=0,lcutm(itype)
            do n=1,nindxm(l,itype)
              csphbes(n,l,itype,iqnrm) = intgrf( rgrid(:ng,itype)*basm(:ng,n,l,itype)*sphbes(:ng,l), itype )
            enddo
            n1             = nindxm(l,itype)
            rdum           = intgrf((rgrid(:,itype)*sphbes(:,l))**2,itype)
            if(rdum==0) cycle
            integrand      = matmul(basm(:ng,:n1,l,itype),csphbes(:n1,l,itype,iqnrm)) - rgrid(:ng,itype)*sphbes(:ng,l)
            error(l,itype) = error(l,itype) + intgrf(integrand**2,itype)/rdum
          enddo
        enddo
      enddo
      write(6,'(2X,A)') 'Accuracy of representation of spherical Bessel functions:'
      do itype=1,ntype
        if(ntype>1) write(6,'(4X,A,I3)') 'Atom type',itype
        write(6,'(4X,2A,F15.8)') (lchar(l),':',sqrt(error(l,itype)/(nqnrm-min(1,l))),l=0,lcutm(itype))
        if(sqrt(error(0,itype)/nqnrm)>0.1) then
          write(6,'(4X,A)') 'The s-Bessel function is not well represented. Large GCUT?'
        endif
      enddo

      ! Calculate the coefficients c(i) for |PW'> => carr3
      do k = 1,nkpt1 ; ikpt = kpt1(k)
        do igpt=1,ngptm(ikpt)
          igptp = pgptm(igpt,ikpt)
          q     = matmul ( rlat, kpt(:,ikpt) + gptm(:,igptp) )
          qnorm = sqrt(sum(q**2))
          iqnrm = pqnrm(igpt,k)
          if(abs(qnrm(iqnrm)-qnorm)>1d-12)
     &      Bug('qnorm does not equal corresponding element in qnrm')
          call harmonicsr(y,q,maxlcutm)
          y  = conjg(y)
          i  = 0
          ic = 0
          do itype=1,ntype
            do ieq=1,neq(itype)
              ic   = ic + 1
              cexp = 4*pi * exp(img* 2*pi *dot_product(gptm(:,igptp),cent(:,ic)))
              lm   = 0
              do l=0,lcutm(itype)
                do m=-l,l
                  lm = lm + 1
                  do n=1,nindxm(l,itype)
                    i               = i + 1
                    carr3(i,igpt,k) = cexp * imgl(l) * y(lm) * csphbes(n,l,itype,iqnrm) / svol
                  enddo
                enddo
              enddo
            enddo
          enddo
        enddo
      enddo

# ifdef INV
      do k = 1,nkpt1 ; ikpt = kpt1(k) ; if(ikpt==1) cycle
        do igpt1 = 1,ngptm(ikpt)
          call symmetrize(carr3(1,igpt1,k),nbasp,1,1)
        enddo
      enddo
# endif

      ! Now we are ready to calculate <PW'|v|PW'>.
      rdum  = 0
      rdum1 = 0

      m = 0
      n = 0

      do k = 1,nkpt1 ; ikpt = kpt1(k) ; if(ikpt==1) cycle
        i = 0
        do ix=1,nbasp
          do iy=1,ix
            i = i + 1
            ! use symmetry property of coulomb: coulomb(ji) = conjg(coulomb(ij))
            coulmat(iy,ix) = coulomb1(i,k)
            coulmat(ix,iy) = MCONJG(coulomb1(i,k))
          enddo
        enddo
        do igpt1=1,ngptm(ikpt)
          igptp1 = pgptm(igpt1,ikpt)
          q1     = matmul ( rlat, kpt(:,ikpt) + gptm(:,igptp1) )
          carr1  = matmul ( conjg(carr3(:,igpt1,k)),coulmat )
          do igpt2=1,igpt1
            igptp2 = pgptm(igpt2,ikpt)
            q2     = matmul ( rlat, kpt(:,ikpt) + gptm(:,igptp2) )
            if(all(q1==0).and.all(q2==0)) cycle
            ! MT-MT
            cdum = sum(carr1*carr3(:,igpt2,k))
            ! PW-PW
            cdum = cdum + MCONJG(coulomb1((nbasp+igpt1)*(nbasp+igpt1-1)/2+nbasp+igpt2,k))
            ! MT-PW, PW-MT
            do j=1,nbasp
              cdum = cdum + conjg(carr3(j,igpt1,k)) *        coulomb1((nbasp+igpt2)*(nbasp+igpt2-1)/2+j,k) +
     &                            carr3(j,igpt2,k)  * MCONJG(coulomb1((nbasp+igpt1)*(nbasp+igpt1-1)/2+j,k))
            enddo
            if(igpt1==igpt2) then
              rdum  = rdum  + abs( 1 - cdum * sum(q1**2) / (4*pi) )**2              ; m = m + 1
            else
              rdum1 = rdum1 + abs( cdum * sqrt(sum(q1**2)*sum(q2**2)) / (4*pi) )**2 ; n = n + 1
            endif
          enddo
        enddo
      enddo

      write(6,'(2X,A)') 'Average deviation of'
      write(6,'(4X,A,F12.8)') 'diagonal elements:    ',sqrt(rdum /m)
      write(6,'(4X,A,F12.8)') 'off-diagonal elements:',sqrt(rdum1/n)

      rdum1 = 1d30 ! ikpt = index of shortest non-zero k-point
      do i = 1,nkpt
        rdum = sum(matmul(rlat,kpt(:,i))**2)
        if(rdum==0) cycle
        if(rdum<rdum1) then ; ikpt = i ; rdum1 = rdum ; endif
      enddo
      do k = 1,nkpt1
        if(kpt1(k)==ikpt) then
          write(6,'(A)') '  (2) One of the eigenvalues e of VC=eSC (V=coulomb matrix, S=overlap matrix)',
     &                   '      should approximate 4*pi/k**2.',
     &                   '      We check this with the shortest, nonzero k-vector.'
          rdum = (4*pi)/sum(matmul(rlat,kpt(:,ikpt))**2)
          n    = nbasm(ikpt)*(nbasm(ikpt)+1)/2
          allocate ( olapm(n),olaphlp(ngptm(ikpt),ngptm(ikpt)),eig(nbasm(ikpt)) )
          call olap_pw(olaphlp,gptm(:,pgptm(:ngptm(ikpt),ikpt)),ngptm(ikpt))
          olapm = packmat ( blockmat ( identity(nbasp)*MUNIT , olaphlp ) )
          call diagonalize(eig,coulomb1(:n,k),olapm)
          write(6,'(2X,A,F15.11)') 'Relative error of nearest eigenvalue:',minval(abs(eig(:)-rdum))/rdum
          deallocate ( olapm,olaphlp,eig )
        endif
      enddo

      call cpu_time(time2) ; write(6,'(2X,A,F7.2)') 'Timing:',time2-time1
      deallocate (coulmat,csphbes,carr1,carr3)

      endif
      endif ! end CHKCOUL

      do k = 1,nkpt1 ; ikpt = kpt1(k)
        if(ikpt==1) call coulomb_sphaverage(coulomb1(:,k),0) ! subtract Coulomb spherical average
      enddo

      Oend
      Nfence(coulomb1)

      if(gcutstep>0) deallocate (pointer,coulstep,gptstep,step,qnrmstep,pqnrmstep,olapstep)
      deallocate ( y,y1,y2,imgl,rarr,qnrm,pqnrm )

      if(tstcoul) call finish

      contains

c     -----------------

c     Returns a list of (k+G) vector lengths in qnrm(1:nqnrm) and the corresponding pointer pqnrm(1:ngpt(ikpt),ikpt)
      subroutine getnorm(kpt,nkpt,gpt,ngpt,pgpt,dim,qnrm,nqnrm,pqnrm)
      use util
      use, intrinsic :: iso_fortran_env
      implicit none
      integer, intent(in)  :: nkpt,ngpt(nkpt),gpt(3,*),dim,pgpt(dim,nkpt)
      real_dp, intent(in)  :: kpt(3,nkpt)
      real_dp, allocatable :: qnrm(:)
      integer, intent(out) :: nqnrm
      integer, allocatable :: pqnrm(:,:)
      integer              :: i,j,ikpt,igpt,igptp
      real_dp              :: q(3),qnorm
      allocate ( qnrm(maxval(ngpt)*nkpt),pqnrm(maxval(ngpt),nkpt) )
      qnrm = 0
      i    = 0
      do ikpt = 1,nkpt
        do igpt = 1,ngpt(ikpt)
          igptp = pgpt(igpt,ikpt) ; if(igptp==0) Bug('zero pointer')
          q     = matmul ( rlat, kpt(:,ikpt) + gpt(:,igptp) )
          qnorm = sqrt(sum(q**2))
          do j=1,i
            if(abs(qnrm(j)-qnorm)<1d-12) then
              pqnrm(igpt,ikpt) = j
              goto 1
            endif
          enddo
          i                = i + 1
          qnrm(i)          = qnorm
          pqnrm(igpt,ikpt) = i
 1      enddo
      enddo
      nqnrm = i
      call reallocate(qnrm,nqnrm)
      end subroutine getnorm

c     -----------------

      subroutine orthonormalize(f,n,itype)
      use global
      use, intrinsic :: iso_fortran_env
      implicit none
      integer, intent(in)    :: itype
      integer, intent(inout) :: n
      real_dp, intent(inout) :: f(maxgrid,n)
      real_dp                :: olap(n,n),work(3*n),eig(n)
      integer                :: i,j,nn,index(n)
      do i=1,n
        f(:,i) = f(:,i) / sqrt(intgrf(f(:,i)**2,itype))
        do j=1,i
          olap(i,j) = intgrf(f(:,i)*f(:,j),itype)
          olap(j,i) = olap(i,j)
        enddo
      enddo
      call dsyev('V','U',n,olap,n,eig,work,3*n,i)
      if(i/=0) Error('dsyev failed.')
      nn = 0
      do i=1,n
        if(eig(i)/n>tolerance) then
          nn        = nn + 1
          index(nn) = i
        endif
      enddo
      do i=1,grid(itype)%number
        f(i,:nn) = matmul(f(i,:),olap(:,index(:nn))) / sqrt(eig(index(:nn)))
      enddo
      write(6,'(I4,A,I4)') n,' ->',nn
      n = nn
      end subroutine orthonormalize

c     -----------------

      function sphbessel_integral()
      use, intrinsic :: iso_fortran_env
      implicit none
      real_dp :: sphbessel_integral
      real_dp :: q1,q2,dq,s,sb01,sb11,sb21,sb31,sb02,sb12,sb22,sb32,a1,a2,da,b1,b2,db,c1,c2,dc,r1,r2
      q1 = qnrm(iqnrm1)
      q2 = qnrm(iqnrm2)
      s  = grid(itype)%radius
      if(q1==0.and.q2==0) then
        if(l>0) then ; sphbessel_integral = 0
        else         ; sphbessel_integral = 2*s**5/15
        endif
      else if(q1==0.or.q2==0) then
        if(l>0)        then ; sphbessel_integral = 0
        else if(q1==0) then ; sphbessel_integral = s**3/(3*q2**2) * ( q2*s * sphbes0(1,itype,iqnrm2) + sphbes0(2,itype,iqnrm2) )
        else                ; sphbessel_integral = s**3/(3*q1**2) * ( q1*s * sphbes0(1,itype,iqnrm1) + sphbes0(2,itype,iqnrm1) )
        endif
      else if(abs(q1-q2)<1d-6) then
        dq                 = q2 - q1
        sphbessel_integral = s**3/(2*q1**2) * ( (2*l+3) * sphbes0(l+1,itype,iqnrm1)**2 -
     &                                          (2*l+1) * sphbes0(l,itype,iqnrm1) * sphbes0(l+2,itype,iqnrm1)
     &                                      + ( (2*l+3) * 5*sphbes0(l+1,itype,iqnrm1)*sphbes0(l-1,itype,iqnrm1) -
     &                                          ( sphbes0(l+1,itype,iqnrm1) + 5*sphbes0(l-1,itype,iqnrm1) ) *
     &                                          q1*s * sphbes0(l,itype,iqnrm1) ) * dq/(2*q1) )
      else ! We use either of two formulas that are stable for high and small q1/q2, respectively.
        sb01 = sphbes0(l-1,itype,iqnrm1)
        sb11 = sphbes0(l  ,itype,iqnrm1)
        sb21 = sphbes0(l+1,itype,iqnrm1)
        sb31 = sphbes0(l+2,itype,iqnrm1)
        sb02 = sphbes0(l-1,itype,iqnrm2)
        sb12 = sphbes0(l  ,itype,iqnrm2)
        sb22 = sphbes0(l+1,itype,iqnrm2)
        sb32 = sphbes0(l+2,itype,iqnrm2)
        dq   = q1**2 - q2**2
        a1   = q2/q1 * sb21 * sb02
        a2   = q1/q2 * sb22 * sb01
        da   = a1 - a2
        b1   = sb31 * sb12
        b2   = sb32 * sb11
        db   = b1 - b2
        c1   = sb21 * sb22 / ( q1*q2 )
        c2   = db / dq * (2*l+1)/(2*l+3)
        dc   = c1 + c2
        r1   = abs(da/a1)
        r2   = min ( abs(db/b1) , abs(dc/c1) )
        ! Ensure numerical stability. If both formulas are not sufficiently stable, the program issues a warning.
        if(r1>r2) then
          if(r1<1d-7) then
            Warn('Formula One possibly unstable. Ratios: '//Chf(r1,'E10.5')//'('//Chf(r2,'E10.5')//')')
            write(0,'(A,2F15.10,I3)')              '                    Current qnorms and atom type:',q1,q2,itype
          endif
          sphbessel_integral = s**3 / dq * da
        else
          if(r2<1d-7) then
            Warn('Formula Two possibly unstable. Ratios: '//Chf(r2,'E10.5')//'('//Chf(r1,'E10.5')//')')
            write(0,'(A,2F15.10,I3)')              '                    Current qnorms and atom type:',q1,q2,itype
          endif
          sphbessel_integral = s**3 * dc
        endif
      endif
      end function sphbessel_integral

c     -----------------

      subroutine stepfunction(step,gptstep,ngptstep)
      use global
      use, intrinsic :: iso_fortran_env
      implicit none
# include "interface/getshells.inc"
      complex_dp, allocatable :: step(:)
      integer,    allocatable :: gptstep(:,:)
      integer,    intent(out) :: ngptstep
      real_dp,    allocatable :: radsh(:)
      integer                 :: i,itype,ieq,icent
      real_dp                 :: q(3),r,fgr
      call getshells(gptstep,ngptstep,radsh,i,gcutstep,rlat,[0d0,0d0,0d0],1,.true.)
      allocate (step(ngptstep))
      step = 0
      do i=1,ngptstep
        q     = matmul ( rlat,gptstep(:,i) )
        qnorm = sqrt(sum(q**2))
        if(qnorm==0) then
          step(i) = 1
          do itype=1,ntype
            r       = grid(itype)%radius
            step(i) = step(i) - neq(itype) * 4*pi*r**3/3 / vol
          enddo
        else
          icent = 0
          do itype=1,ntype
            r   = grid(itype)%radius
            fgr = 4*pi* ( sin(qnorm*r) - (qnorm*r) * cos(qnorm*r) ) / qnorm**3 / vol
            do ieq=1,neq(itype)
              icent   = icent + 1
              step(i) = step(i) - fgr * exp( - img * 2*pi*dot_product(cent(:,icent),gptstep(:,i)) )
            enddo
          enddo
        endif
      enddo
      end subroutine stepfunction

c     -----------------

      end

c     -----------------

c     Calls routine coulombmatrix for all required kpoints.

# include "jobtype.h"

      subroutine coulombmatrix0
      use global
      use, intrinsic :: iso_fortran_env
      use, intrinsic :: iso_c_binding
      Mpi( use Mwrapper )
      implicit none
      integer                  :: i
# ifdef MPI
      type(c_ptr)              :: ptr
      integer                  :: Merr
      MCOMPLEX_dp, pointer_cnt :: Ninit_coulomb(:)
# endif      
      if(associated(coulomb)) Bug('Array coulomb already allocated.')
      if(lkptadd.and.any([(any(job(i)%type==[J_SUSR,J_DIEL,J_SCR]),i=1,njob)]) ) then
        Nallocate0 ( coulomb, (S_ maxbasm*(maxbasm+1)/2,nkpti+1 S_) )
        call coulombmatrix(Win(coulomb),maxbasm,nkpti+1,[ [(i,i=1,nkpti)] , nkpt+1 ],.true.)
      else
        Nallocate0 ( coulomb, (S_ maxbasm*(maxbasm+1)/2,nkpti S_) )
        call coulombmatrix(Win(coulomb),maxbasm,nkpti,    [(i,i=1,nkpti)]           ,.true.)
      endif
      end

c     -----------------

c     Calculate body of Coulomb matrix at Gamma point: v_IJ = SUM(G) c^*_IG c_JG 4*pi/G**2 .
c     For this we must subtract from coulomb(:,1) the spherical average of a term that comes
c     from the fact that MT functions have k-dependent Fourier coefficients (see script).
c     mode = 0 : default: subtracts (as explained), packed matrix, <M|v|M>
c     mode = 1 : partially undoes subtraction (adds term that would destroy multipole-sparsity)
c     mode = 2 : full Coulomb matrix [coulomb(:dim,:dim)]
c     mode = 4 : matrix is interpreted as <M~|v|M~> (M~ biorthogonal set) instead of <M|v|M>
c     All combinations of modes are possible.
      subroutine coulomb_sphaverage(coulomb1,mode)
      use global
      use wrapper
      use, intrinsic :: iso_fortran_env
      implicit none
      integer,     intent(in)    :: mode      
      MCOMPLEX_dp, intent(inout) :: coulomb1(*)
      complex_dp                 :: coeff(nbasm(1)),cderiv(nbasm(1),-1:1),claplace(nbasm(1))
      integer                    :: l,m,i,j,n,itype,ieq
      real_dp                    :: intgrf
      MCOMPLEX_dp                :: stepfunction
      if(mode<0.or.mode>7) Bug('mode out of range.')
      ! Define coefficients (coeff) and their derivatives (cderiv,claplace)
      n        = nbasm(1)
      coeff    = 0
      cderiv   = 0
      claplace = 0
      j        = 0
      do itype = 1,ntype
        do ieq = 1,neq(itype)
          do l = 0,lcutm(itype)
            do m = -l,l
              do i = 1,nindxm(l,itype)
                j = j + 1
                if(l==0) then ; if(i==1) coeff(j) = sqrt( 4*pi*grid(itype)%radius**3 / (3*vol) ) ! identical to ...
                  coeff(j)    =  sqrt(4*pi)       * intgrf( rgrid(:,itype)    * basm(:,i,0,itype) , itype ) / sqrt(vol)
                  claplace(j) = -sqrt(4*pi)       * intgrf( rgrid(:,itype)**3 * basm(:,i,0,itype) , itype ) / sqrt(vol)
                  !write(*,'(2I3,2F12.7)') l,i,real(coeff(j)),real(claplace(j))
                else if(l==1) then
                  cderiv(j,m) = -sqrt(4*pi/3)*img * intgrf( rgrid(:,itype)**2 * basm(:,i,1,itype) , itype ) / sqrt(vol)
                  !write(*,'(2I3,F12.7)') l,i,imag(cderiv(j,m))
                endif
              enddo
            enddo
          enddo
        enddo
      enddo
      if(iand(mode,4)/=0) then
        coeff(nbasp+1) = 1
      else
        do i = 1,ngptm(1)
          coeff(nbasp+i) = stepfunction(-gptm(:,pgptm(i,1)))
        enddo
      endif
# ifdef INV
      call symmetrize(coeff,       1,nbasm(1),2)
      call symmetrize(claplace,    1,nbasm(1),2)
      call symmetrize(cderiv(:,-1),1,nbasm(1),2)
      call symmetrize(cderiv(:, 0),1,nbasm(1),2)
      call symmetrize(cderiv(:, 1),1,nbasm(1),2)
# endif
      !if(.not.allocated(constfunc)) allocate ( constfunc(n) )
      ! Subtract (or partially add) head contributions from (to) coulomb(:,1) to obtain the body
      l = 0
      do j = 1,n
        do i = 1,n ; if(iand(mode,2)==0.and.i>j) exit
          l = l + 1
          if(iand(mode,1)==0) then ! standard: subtract spherical average
            coulomb1(l) = coulomb1(l) - 4*pi/3 * ( dot_product(cderiv(i,:),cderiv(j,:)) +
     &                                           ( conjg(coeff(i)) * claplace(j) + conjg(claplace(i)) * coeff(j) ) / 2 )
          else                     ! partially undo subtraction: add term that would destroy Coulomb multipole-sparsity
            coulomb1(l) = coulomb1(l) + 4*pi/3 * ( conjg(coeff(i)) * claplace(j) + conjg(claplace(i)) * coeff(j) ) / 2
          endif
        enddo
      enddo
# if 0
      coeff(nbasp+1)  = 1d0
      coeff(nbasp+2:) = 0d0
# ifdef INV
      call desymmetrize(coeff,1,nbasm(1),2)
      call   symmetrize(coeff,nbasm(1),1,1)
# endif
      ! Explicit normalization here in order to prevent failure of the diagonalization in diagonalize_coulomb
      ! due to inaccuracies in the overlap matrix (which can make it singular).
      constfunc = coeff / sqrt ( ( sum(abs(coeff(:nbasp))**2) + stepfunction([0,0,0]) ) )
c      constfunc = coeff / sqrt ( ( sum(abs(coeff(:nbasp))**2) + dotprod ( coeff(nbasp+1:), matmul(olap,coeff(nbasp+1:)) ) ) )
# endif      

      end

c     -----------------

c     Diagonalizes the Coulomb matrix at the k point ikpt and returns eigenvalues in coul and eigenvectors in ctrafo.
c     nbasm0 is the number of requested/found eigenvalues and eigenvectors.
c
c     At the Gamma point:
c
c     The eigenvector to the eigenvalue zero is already defined in subtract_sphaverage above (constfunc).
c     This can be viewed as a constraint for the eigenvalue problem. We take this into account by Gram-Schmidt
c     orthogonalizing the basis functions |M> to constfunc |c>:
c
c     |M'> = |M> - <c|M> |c> .
c
c     We leave out the interstitial constant function (index nbasp+1), so the new basis |M'> has one basis function less.
c     The new Coulomb and overlap matrices are then given by
c
c     <M'|v|M'> = <M|v|M> - h <M|c> <c|M>  ( with h = <c|v|c> \approx 0 )
c     <M'|M'>   = <M|M> - <M|c> <c|M> .
c
c     Afterwards the resulting eigenvectors are back-transformed to the basis |M>.
c
      subroutine diagonalize_coulomb(Win(ctrafo),coul,Win(olap),coulomb1,ikpt,nbasm0,x)
      use global
      use wrapper
      use util
      use key
      use, intrinsic :: iso_fortran_env
      use, intrinsic :: iso_c_binding
      Mpi ( use Mwrapper )
      implicit none
      integer,     intent(in)    :: ikpt MpiC(win_ctrafo) MpiC(win_olap)
      integer,     intent(inout) :: nbasm0
      real_dp,     intent(in)    :: x
      real_dp,     intent(out)   :: coul(nbasm(ikpt))
      MCOMPLEX_dp, intent(in)    :: coulomb1(nbasm(ikpt)*(nbasm(ikpt)+1)/2)
      MCOMPLEX_dp, intent(out)   :: ctrafo(nbasm(ikpt),nbasm(ikpt))
      MCOMPLEX_dp, intent(out)   :: olap(ngptm(ikpt),ngptm(ikpt))
      MCOMPLEX_dp                :: cdum
      MCOMPLEX_dp, pointer_cnt   :: eigv(:,:)
      MCOMPLEX_dp, pointer_cnt   :: proj(:,:)
      complex_dp                 :: cexp,cproj(nbasp),ylm((maxlcutm+1)**2)
      real_dp                    :: jl(0:maxlcutm+1),dj
      real_dp                    :: val1,val2,val3,dval1,dval2,dval3,wronsk,a,b
      real_dp                    :: g(3),kg(3),kgr(3),kgn,kgnr
      integer                    :: n1,n2
      integer                    :: itype,ieq,ic,l,m,n,lm,igpt,nbasm1,nlo
      real_dp,     allocatable   :: eig(:)
      real_dp                    :: rdum
      integer                    :: i,j,k,nn,z,ncut
      real                       :: cputime
      real_dp                    :: intgrf
      MCOMPLEX_dp                :: stepfunction
      MCOMPLEX_dp, pointer_cnt   :: olapmat(:,:),coulmat(:,:),helpmat(:,:),helpmat1(:,:)
# ifdef MPI
      integer                    :: win_olapmat,win_coulmat,win_proj,win_helpmat,win_helpmat1,win_eigv,Merr
      type(c_ptr)                :: ptr
# endif

      Rwrite(6,'(A'NoA) 'Diagonalize Coulomb matrix... '
      call cpu_time(cputime)

      coul = 0

c     Calculate overlap <M|M> -> overlap
      Nfence(olap) ; ifO olap = 0 ; Nfence(olap)
      k = -1
      do j = 1,ngptm(ikpt)
        do i = 1,j ; NcycleP(k)
          olap(i,j) = stepfunction( gptm(:,pgptm(i,ikpt)) - gptm(:,pgptm(j,ikpt)) )
          olap(j,i) = MCONJG(olap(i,j))
        enddo
      enddo
      Nfence(olap)

      z = 0 ; if(ikpt==1) z = 1 ! add G=0 term to proj matrix for k=0
      nbasm1 = nbasm(ikpt) - z  ! constant function will be taken out (for k=0)

c     Define constant function, first as <M|c>
      Nfence(ctrafo)
      Obegin
      ctrafo = 0
      if(ikpt==1) then
        ic          = 0
        i           = 1
        do itype = 1,ntype
          rdum = sqrt( 4*pi*grid(itype)%radius**3 / (3*vol) )
          do ieq = 1,neq(itype)
            ic = ic + 1
            Inv( if(pcent(ic,invsym)==ic) ) ctrafo(i,1) = rdum
            Inv( if(pcent(ic,invsym)> ic)   ctrafo(i,1) = rdum * sqrt(2d0) )
            i  = i + sum( [ ((2*l+1)*nindxm(l,itype),l=0,lcutm(itype)) ] )
          enddo
        enddo
        do i = 1,ngptm(1)
          ctrafo(nbasp+i,1) = stepfunction( gptm(:,pgptm(i,1)) )
        enddo
      endif
      Oend
      Nfence(ctrafo)

c
c     Construct projection matrix -> proj
c

c     (A) NOAPW: trivial proj (project out constant function if k=0)

      if(noapw) then
        Nallocate ( proj, (S_ nbasp+z,nbasm1 S_) )
        Obegin
        proj = 0
        k    = 0
        do i = 1,nbasm1 ; if(ikpt==1.and.i==nbasp+1) k = k + 1
          k = k + 1
          if(k<=size(proj,1)) proj(i,i) = 1
          if(ikpt==1.and.ctrafo(k,1)/=0) then
            do j = 1,nbasp
              proj(j,i) = proj(j,i) - ctrafo(j,1) * MCONJG(ctrafo(k,1))
            enddo
            proj(nbasp+1,i) = -MCONJG(ctrafo(k,1))
          endif
        enddo
        Oend
        Nfence(proj)
        nlo = nbasp

      else

c     (B) MPB-APW construction

        ! (1) Number of MPB-APWs
        do itype = 1,ntype
          do l = 0,lcutm(itype)
            if(nindxm(l,itype)==1) then
              val1 = basm(grid(itype)%number,1,l,itype) / grid(itype)%radius
              if(val1>1d-4) nbasm1 = nbasm1 - (2*l+1) * neq(itype)
            else if(nindxm(l,itype)>=2) then
              nbasm1 = nbasm1 - 2 * (2*l+1) * neq(itype)
            endif
          enddo
        enddo
        Rwrite(6,'(A'NoA) ' (MPB-APW/LOs: '//Chr(nbasm1)//') ...'
        Nallocate ( proj, (S_ nbasp+z,nbasm1 S_) )

        ! (2) MPB-LOs
        Obegin
        proj = 0
        i    = 0 ! i = MPB index
        j    = 0 ! j = APW-LO index
        ic   = 0
        do itype = 1,ntype
          do ieq = 1,neq(itype)
            ic = ic + 1
            do l = 0,lcutm(itype)

              if(nindxm(l,itype)==1) then
                val1 = basm(grid(itype)%number,1,l,itype) / grid(itype)%radius
                if(val1<=1d-4) then ! treat as LO (instead of in APWs)
                  do m = 1,2*l+1
                    proj(i+m,j+m) = 1
                  enddo
                  j = j + (2*l+1)
                endif
              else if(nindxm(l,itype)>=3) then
                call maxwronsk(n1,n2)
                nn = nindxm(l,itype)
                do n = 1,nn
                  if(n==n1.or.n==n2) cycle
                  val1   = basm(grid(itype)%number,n1,l,itype) / grid(itype)%radius
                  val2   = basm(grid(itype)%number,n2,l,itype) / grid(itype)%radius
                  val3   = basm(grid(itype)%number,n, l,itype) / grid(itype)%radius
                  dval1  = dbasm(n1,l,itype)
                  dval2  = dbasm(n2,l,itype)
                  dval3  = dbasm(n, l,itype)
                  wronsk = val1 * dval2 - dval1 * val2
                  a      = - ( val3 * dval2 - dval3 * val2 ) / wronsk
                  b      =   ( val3 * dval1 - dval3 * val1 ) / wronsk
                  rdum   = sqrt ( intgrf( ( a * basm(:,n1,l,itype) + b * basm(:,n2,l,itype) + basm(:,n,l,itype) )**2 , itype ) )
                  j      = j + 1
                  do m = 0,2*l
                    proj(i+m*nn+n1,j+m*(nn-2)) = a / rdum
                    proj(i+m*nn+n2,j+m*(nn-2)) = b / rdum
                    proj(i+m*nn+n, j+m*(nn-2)) = 1 / rdum
                  enddo
                enddo
                j = j + 2*l * ( nindxm(l,itype) - 2 )
              endif
              i = i + (2*l+1) * nindxm(l,itype)

            enddo
          enddo
        enddo
        nlo = j
        Oend
        Nfence(proj)
        Mpi( call Mcast(nlo) )

        ! (3) APWs
        do igpt = 1,ngptm(ikpt)
          if(ikpt==1.and.igpt==1) cycle ; Ncycle(igpt)
          g     =               gptm(:,pgptm(igpt,ikpt))
          kg    = kpt(:,ikpt) + gptm(:,pgptm(igpt,ikpt))
          kgr   = matmul ( rlat, kg )
          kgn   = sqrt(sum(kgr**2))
          call harmonicsr(ylm,kgr,maxlcutm)
          ylm   = conjg(ylm)
          cproj = 0
          i     = 0
          ic    = 0
          do itype = 1,ntype
            do ieq = 1,neq(itype)
              ic   = ic + 1
              cexp = exp(img * 2*pi * dot_product(g,cent(:,ic)) ) * 4*pi / sqrt(vol)
              kgnr = kgn * grid(itype)%radius
              lm   = 0
              call sphbessel(jl,kgnr,lcutm(itype)+1)
              do l = 0,lcutm(itype)
                dj = 0
                if(kgnr==0) then ; if(l==1) dj = 1d0/3
                else             ;          dj = jl(l)*l/kgnr - jl(l+1)
                endif
                dj = kgn * dj

                if(nindxm(l,itype)==1) then ! match value
                  n1   = 1
                  val1 = basm(grid(itype)%number,1,l,itype) / grid(itype)%radius
                  if(val1>1d-4) then ; a = jl(l) / val1
                  else               ; a = 0
                  endif
                else if(nindxm(l,itype)>=2) then ! match value and derivative
                  call maxwronsk(n1,n2)
                  val1   = basm(grid(itype)%number,n1,l,itype) / grid(itype)%radius
                  val2   = basm(grid(itype)%number,n2,l,itype) / grid(itype)%radius
                  dval1  = dbasm(n1,l,itype)
                  dval2  = dbasm(n2,l,itype)
                  wronsk = val1 * dval2 - dval1 * val2
                  a      =   ( jl(l) * dval2 - dj * val2 ) / wronsk
                  b      = - ( jl(l) * dval1 - dj * val1 ) / wronsk
                endif

                do m = -l,l
                  lm          = lm + 1
                  cproj(i+n1) = cexp * ylm(lm) * a
                  if(nindxm(l,itype)>=2)
     &            cproj(i+n2) = cexp * ylm(lm) * b
                  i = i + nindxm(l,itype)
                enddo

                cexp = cexp * img ! img**l
              enddo
            enddo
          enddo
          if(i/=nbasp) Bug('count error')
# ifdef INV
          call symmetrize(cproj,nbasp,1,1)
          if(any(abs(imag(cproj))>1d-10)) Bug('Projection coefficients not real.')
# endif
          if(ikpt==1) then ; proj(:nbasp,nlo+igpt-1) = cproj
          else             ; proj(:nbasp,nlo+igpt)   = cproj
          endif
        enddo

        Nfence(proj)

        ! (4) Project out constant function (k=0)
        if(ikpt==1) then
          do i = 1,nbasm1 ; Ncycle(i)
            cdum            = dotprod(ctrafo(:nbasp+1,1),proj(:,i)) ; if(i>nlo) cdum = cdum + MCONJG(ctrafo(nbasp+i-nlo+1,1)) ! <c|E>
            proj(:nbasp,i)  = proj(:nbasp,i)  - cdum * ctrafo(:nbasp,1)                                                       ! |E> - <c|E> |c>
            proj(nbasp+1,i) = proj(nbasp+1,i) - cdum                                                                          !
          enddo
          Nfence(proj)
        endif

      endif

c     Projection of overlap matrix
      MrangeDef1(n1,n2,nbasm1)
      ! Right multiplication
      Nallocate ( helpmat, (S_ nbasm(ikpt),nbasm1 S_) )
      Obegin
      helpmat                  = 0
      helpmat(:nbasp,:)        = proj(:nbasp,:)
      helpmat(nbasp+1:,nlo+1:) = olap(:,z+1:)
      if(ikpt==1) then
        do i = 1,nbasm1
          if(proj(nbasp+1,i)/=0) helpmat(nbasp+1:,i) = helpmat(nbasp+1:,i) + olap(:,1) * proj(nbasp+1,i)
        enddo
      endif
      Oend
      Nfence(helpmat)
      ! Left multiplication of overlap
      Nallocate ( olapmat,  (S_ nbasm1, nbasm1 S_) )
      Nallocate ( helpmat1, (S_ nbasp+z,nbasm1 S_) )
      ifO olapmat  = 0
      ifO helpmat1 = helpmat(:nbasp+z,:)
      Nfence(olapmat)
      Nfence(helpmat1)      
      CHKMEM0(0)
      olapmat(n1:n2,:) = macmat( proj(:,n1:n2) , helpmat1 )
      Nfence(olapmat)
      Nfence(helpmat1)
      Ndeallocate(helpmat1)
      Obegin ; Mpi( call Msum(olapmat,comm=Ocomm) )
      CHKMEM0(1)
      do j = 1,nbasm1
        do i = 1,nbasm1-nlo
          olapmat(nlo+i,j) = olapmat(nlo+i,j) + helpmat(nbasp+z+i,j)
        enddo
      enddo
      CHKMEM0(2)
      Oend
      Nfence(olapmat)

c     CUTZERO: Remove linear dependencies
      ncut = 0
      if(cutzero>0) then
        Ndeallocate( helpmat ) ! will be used temporarily and reallocated below
        Rbegin
        allocate ( eigv(nbasm1,nbasm1),eig(nbasm1) )
        call diagonalize ( eigv,eig,olapmat )
        ncut = count(eig<cutzero)
        if(ncut==nbasm1) Error('All eigenvalues smaller than CUTZERO.')
        Rend
        Mpi( call Mcast(ncut) )
        Rwrite(6,'('' (cut:'',I3,'') '''NoA) ncut
        if(ncut>0) then
          nbasm1 = nbasm1 - ncut
          allocate(helpmat(nbasm(ikpt),nbasm1))
          Rbegin
          helpmat(:nbasp+z,  :nbasm1) = matmat(proj,eigv(:,     ncut+1:))
          helpmat(nbasp+z+1:,:nbasm1) =             eigv(nlo+1:,ncut+1:)
          Rend
          Ndeallocate(proj)
          Ndeallocate(olapmat)
          Nallocate ( proj,(S_ nbasm(ikpt),nbasm1 S_) )
          Nallocate ( olapmat,(S_ nbasm1,nbasm1 S_) )
          Nallocate ( coulmat,(S_ nbasm1,nbasm1 S_) )
          Rbegin
          proj    = helpmat(:,:nbasm1)
          coulmat = macmat(proj,matmat(coulomb1,proj))
          olapmat = macmat(proj(:nbasp,:),proj(:nbasp,:))
     &            + macmat(proj(nbasp+1:,:),matmat(olap,proj(nbasp+1:,:)))
          deallocate(helpmat)
          Rend
          Mpi( Ocall Mcast(proj,   comm=Ocomm) )
          Mpi( Ocall Mcast(olapmat,comm=Ocomm) )
          Mpi( Ocall Mcast(coulmat,comm=Ocomm) )
        endif
        ifR deallocate ( eig,eigv )
        Nallocate( helpmat, (S_ nbasm(ikpt),nbasm1 S_) )
      endif

c     Projection of Coulomb matrix
      if(ncut==0) then
        ! Right multiplication of coulomb
        Nallocate ( coulmat, (S_ nbasm(ikpt),nbasm(ikpt) S_) )
        Nfence(helpmat)
        ifO helpmat = 0                        ; Nfence(helpmat)
        ifO call p_unpackmat(coulmat,coulomb1) ; Nfence(coulmat)        
        helpmat(:,n1:n2) = matmat( coulmat(:,:nbasp+z) , proj(:,n1:n2) )
        Nfence(helpmat)
        Obegin ; Mpi( call Msum(helpmat,comm=Ocomm) )
        do j = 1,nbasm1-nlo
          do i = 1,nbasm(ikpt)
            helpmat(i,nlo+j) = helpmat(i,nlo+j) + coulmat(i,nbasp+z+j)
          enddo
        enddo
        Oend
        Nfence(helpmat)
        Ndeallocate(coulmat)
        ! Left multiplication of coulomb
        Nallocate(coulmat, (S_ nbasm1, nbasm1 S_) )
        Nallocate(helpmat1,(S_ nbasp+z,nbasm1 S_) )        
        ifO coulmat  = 0
        ifO helpmat1 = helpmat(:nbasp+z,:)
        Nfence(coulmat)
        Nfence(helpmat1)
        coulmat(n1:n2,:) = macmat( proj(:,n1:n2) , helpmat1 )
        Nfence(coulmat)
        Nfence(helpmat1)
        Ndeallocate(helpmat1)
        Obegin ; Mpi( call Msum(coulmat,comm=Ocomm) )
        do j = 1,nbasm1
          do i = 1,nbasm1-nlo
            coulmat(nlo+i,j) = coulmat(nlo+i,j) + helpmat(nbasp+z+i,j)
          enddo
        enddo
        Oend
        Nfence(coulmat)
        Ndeallocate(helpmat)
      endif

c     Diagonalize
      call cpu_time(cputime)
      nbasm0 = nbasm0 - z
      allocate  ( eig(nbasm1) )
      Nallocate ( eigv,(S_ nbasm1,nbasm1 S_) )
      n1 = 1
      if(nbasm0>=nbasm1) then
        call Mdiagonalize ( Win(eigv),eig,coulmat,olapmat MpiC(lwrite=.true.) )
        nbasm0 = nbasm1
      else if(nbasm0>0) then
        call Mdiagonalize ( Win(eigv),eig,coulmat,olapmat,nbasm1-nbasm0,nbasm1,0d0,0d0 MpiC(lwrite=.true.) )
        nbasm0 = nbasm1
        do while(abs(eig(1)-eig(n1))<1d-8)
          n1 = n1 + 1 ; if(n1>nbasm0) Error('Eigenvalue index out of range.')
        enddo
      else
        call Mdiagonalize ( Win(eigv),eig,coulmat,olapmat,1,nbasm0,x,1d30 MpiC(lwrite=.true.) )
      endif
      Nfence(eigv)
      Mpi( Ocall Msum(eigv,comm=Ocomm) )
      Nfence(eigv)

      Rcall cpu_done(cputime)

      Rif(any(eig(n1:nbasm0)<0))
     &  Error('Negative eigenvalue. Try to increase LEXP or use CUTZERO.')

c     Back projection
      Nfence(ctrafo)
      coul = 0
      MnoR(ifO ctrafo(:,1) = 0 )
      Rif(ikpt==1) then
        ctrafo(nbasp+1, 1) = 1 ! ctrafo(:,1) = <M~|c>
        ctrafo(nbasp+2:,1) = 0 !
        rdum               = dotprod ( ctrafo(:,1) , matvec(coulomb1,ctrafo(:,1)) )
        write(6,'(A,F17.12)') 'Residual eigenvalue:     ',rdum
        if(abs(rdum)>1d-7) Error('Large residual eigenvalue. Increase LEXP.')
      endif
      m = z
      do i = nbasm0,n1,-1
        m       = m + 1
        coul(m) = eig(i) ; Mcycle(m)
        if(ncut==0) then
          ctrafo(:nbasp+z,m)   = matmul(proj,eigv(:,i))
          ctrafo(nbasp+z+1:,m) =             eigv(nlo+1:,i)
        else
          ctrafo(:,m)          = matmul(proj,eigv(:,i))
        endif
      enddo
      nbasm0 = m
      Nfence(ctrafo)
      Mpi( Ocall Msum(ctrafo,comm=Ocomm) )
      Nfence(ctrafo)
# if 0
      write(60+Mrank,*) Orank,Nrank
      do i = 1,nbasm0
        write(60+Mrank,*) i,sum(abs(ctrafo(:,i))**2)
      enddo
      Nfence(ctrafo)
      call finish
# endif
c     Rcall cpu_done(cputime)

      Ndeallocate ( coulmat )
      Ndeallocate ( olapmat )
      Ndeallocate ( proj )
      Ndeallocate ( eigv )

# if 0
      Rbegin
      rdum = 0
      do i = 1,nbasm0
        rdum = rdum + abs( coul(i) - dotprod(ctrafo(:,i),matvec(coulomb1,ctrafo(:,i))) )
c        write(*,'(F20.16)') abs( coul(i) - dotprod(ctrafo(:,i),matvec(coulomb1,ctrafo(:,i))) )
      enddo
      write(*,*) 'coulomb-error:',rdum
      allocate ( olapmat(nbasm(ikpt),nbasm(ikpt)) )
      olapmat                    = 0
      olapmat(:nbasp,:nbasp)     = identity(nbasp)
      olapmat(nbasp+1:,nbasp+1:) = olap
      rdum = 0
      do i = 1,nbasm0
        do j = 1,nbasm0
          if(i==j) then ; rdum = rdum + abs( 1 - dotprod(ctrafo(:,i),matvec(olapmat,ctrafo(:,j))) )
          else          ; rdum = rdum + abs(     dotprod(ctrafo(:,i),matvec(olapmat,ctrafo(:,j))) )
          endif
        enddo
      enddo
      write(*,*) 'olap-error',rdum
      deallocate ( olapmat )
      Rend
# endif

      Nfence(ctrafo)

      contains
      subroutine maxwronsk(n1,n2)
      use, intrinsic :: iso_fortran_env
      implicit none
      integer, intent(out) :: n1,n2
      integer              :: i1,i2
      real_dp              :: val1,val2,dval1,dval2,wronsk,mwronsk
      n1      = 0
      mwronsk = 0
      do i2 = 2,nindxm(l,itype)
        do i1 = 1,i2-1
          val1   = basm(grid(itype)%number,i1,l,itype) / grid(itype)%radius
          val2   = basm(grid(itype)%number,i2,l,itype) / grid(itype)%radius
          dval1  = dbasm(i1,l,itype)
          dval2  = dbasm(i2,l,itype)
          wronsk = abs( val1 * dval2 - dval1 * val2 )
          if(wronsk>mwronsk) then
            n1      = i1
            n2      = i2
            mwronsk = wronsk
          endif
        enddo
      enddo
      if(n1==0) Bug('Could not determine MT pair.')
      end subroutine maxwronsk

      end

c     -----------------

c     Calculates the structure constant
c                                                        1               *      ^
c     structconst(lm,ic1,ic2,k) = SUM exp(ikT) -----------------------  Y  ( T + R(ic) )
c                                  T           | T + R(ic1) - R(ic2) |   lm
c
c     with T = lattice vectors
c
c     An Ewald summation method devised by O.K. Andersen is used for l<8 (see e.g. H.L. Skriver, "The LMTO method", Springer 1984).
c     (The real-space function G can be calculated with gfunction.f)
c

# define CONVPARAM1 1d-10
# define CONVPARAM2 (CONVPARAM1*10)
# define CONVTYPE   2
c
c CONVPARAM1 is the convergence parameter for the Ewald sums (l<=7).
c CONVPARAM2 is the convergence parameter for the real-space sums (l>=8). (Harder to converge!)
c
c The shells are summed until a shell contributes with less than CONVPARAM or, in other words, until a T (or G) point contributes with less than
c CONVPARAM / Npoints(r)
c where  Npoints = 4*pi*r**2 / R**2  estimates the number of points on the shell with r=|T| (or |G|) and R=vol^(1/3) [or rvol^(1/3)].
c
c Type of convergence criterion: CONVPARAM refers to
c CONVTYPE = 0 : absolute values of scaled structure constant, i.e., Ewald scaling (scale) not taken into account (for testing)
c CONVTYPE = 1 : absolute values of structure constants
c CONVTYPE = 2 : relative to the estimated contribution of the first shell Npoints(R) / R**(l+1) = 4*pi / R**(l+1)
c

c Splitting parameter is set automatically by default (EWALD_SCALE=0). Can be set to an explicit value (>0) here.
# define EWALD_SCALE 0

# define HLP9 (1+a*(1+a/2*(1+a/3*(1+a/4*(1+a/5*(1+a/6*(1+a/7*(1+a/8*(1+a/9*

      subroutine structureconstant(Win(structconst),kpt1,nkpt1,lwrite)

      use global
      Mpi( use Mwrapper )

      use, intrinsic :: iso_fortran_env
      implicit none
# include "interface/getshells.inc"
      integer,    intent(in)  :: nkpt1 MpiC(win_structconst)
      logical,    intent(in)  :: lwrite
      real_dp,    intent(in)  :: kpt1(3,nkpt1)
      complex_dp, intent(out) :: structconst((2*lexp+1)**2,ncent,ncent,nkpt1)
      complex_dp              :: shlp((2*lexp+1)**2,nkpt1)
      real_dp                 :: ewaldscale
      real_dp,    allocatable :: radsh(:),shift(:,:)
      integer,    allocatable :: ptsh(:,:),pnt(:)
      logical                 :: lconv(0:2*lexp)
      complex_dp              :: y((2*lexp+1)**2)
      complex_dp              :: cexp,cdum
      real_dp                 :: pref,scale,convpar(0:lexp*2),rad,rrad,rc(3),ra(3),a,a1,aa,g(0:lexp*2),k(3),ki(3),ka(3)
      real_dp                 :: rexp,rdum,ew1,ew2,pts,latcon,factor
      integer                 :: nptsh,nshell,maxl,nshift
      integer                 :: ic1,ic2,ikpt,i,l,m,lm,l1,l2 MpiC(mm) MpiC(Merr)
      real                    :: time1,time2

      latcon = vol**(1d0/3) ! define "average lattice parameter"

      Rbegin
      if(lwrite) write(6,'(//A)') '### subroutine: structureconstant ###'
      ewaldscale = EWALD_SCALE
      if(ewaldscale==0) then
        pts = points(1d0)
        ew1 = 1 ; do ; ew1 = ew1/(ew1+1) ; if(points(ew1)>pts) exit ; enddo
        ew2 = 1 ; do ; ew2 = ew2 + 1     ; if(points(ew2)>pts) exit ; enddo
        ew2 = (ew2-ew1)/1000
        pts = points(ew1)
        do while(points(ew1+ew2)<pts) ; ew1 = ew1 + ew2 ; pts = points(ew1) ; enddo
        ewaldscale = ew1
      endif
      Rend
      Mpi( call Mcast(ewaldscale) )
      pts = points(ewaldscale)

c     Minimum radius (Convergence for l>7 sets a lower limit)
# if   CONVTYPE == 0
      rdum = (4*pi/(CONVPARAM2*latcon**2*scale**9))**(1d0/7) ! 4pi*r**2/R**2 / (r*scale)**(l+1) = CONVPARAM2
# elif CONVTYPE == 1
      rdum = (4*pi/(CONVPARAM2*latcon**2))**(1d0/7)          ! 4pi*r**2/R**2 / r**(l+1) = CONVPARAM2
# elif CONVTYPE == 2
      rdum = (1/CONVPARAM2)**(1d0/7) * latcon                ! 4pi*r**2/R**2 / r**(l+1) = CONVPARAM2 * 4pi / R**(l+1)
# else
#   error CONVTYPE unknown.
# endif

      Rbegin
      if(lwrite) then
        write(6,'(/A,F11.5)')   'Ewald scale: ',ewaldscale
        write(6,'(A,2F11.5'NoA) 'Cutoff radii:',max(rad,rdum),rrad
        if(rdum>rad) write(6,'(A,F10.5'NoA) '   increased from',rad
        write(6,'(//A)') 'Real-space sum'
      else
        write(6,'(A)')        'Calculate structure constant ... '
        write(6,'(A,F11.5)')  '    Ewald scale: ',ewaldscale
        write(6,'(A,2F11.5)') '    Cutoff radii:',max(rad,rdum),rrad
      endif
      Rend

      rad = max(rdum,rad)

c
c     Determine atomic shells
      nshift = ncent*(ncent-1)/2 + 1
      allocate ( shift(3,nshift) )
      i = 0
      do ic2 = 1,ncent
        do ic1 = 1,max(1,ic2-1)
          i          = i + 1
          shift(:,i) = matmul(lat, cent(:,ic2) - cent(:,ic1) )
        enddo
      enddo
      call getshells(ptsh,nptsh,radsh,nshell,rad,lat,shift,nshift,lwrite)

      allocate ( pnt(nptsh) )

      Nfence(structconst)
      ifO structconst = 0
      Nfence(structconst)
      Mpi( mm = -1 )
      factor  = 4*pi/scale**2/latcon**2
      call cpu_time(time1)
      do ic2 = 1,ncent
        do ic1 = 1,max(1,ic2-1) ; McycleP(mm)
          rc = matmul(lat,cent(:,ic2)-cent(:,ic1))
          do i = 1,nptsh
            ra       = matmul(lat,ptsh(:,i)) + rc
            radsh(i) = sqrt(sum(ra**2))
          enddo
          call rorderpf(pnt,radsh,nptsh)
          maxl  = 2*lexp
          lconv = .false.
          a1    = 0
          shlp  = 0
          loop_ptsh: do i = 1,nptsh
            ra = matmul(lat,ptsh(:,pnt(i))) + rc
            a  = scale * sqrt(sum(ra**2))
            if(a==0) cycle
            if(abs(a-a1)>1d-10) then
              a1   = a
              rexp = exp(-a)
              if(.not.lconv(0)) g(0) = rexp / a    * (1+a*11/16*(1+a*3/11*(1+a/9)))
              if(.not.lconv(1)) g(1) = rexp / a**2 * (1+a*(1+a/2*(1+a*7/24*(1+a/7))))
              if(.not.lconv(2)) g(2) = rexp / a**3 * (1+a*(1+a/2*(1+a/3*(1+a/4*(1+a*3/16*(1+a/9))))))
              if(.not.lconv(3)) g(3) = rexp / a**4 * (1+a*(1+a/2*(1+a/3*(1+a/4*(1+a/5*(1+a/6*(1+a/8)))))))
              if(.not.lconv(4)) g(4) = rexp / a**5 * (1+a*(1+a/2*(1+a/3*(1+a/4*(1+a/5*(1+a/6*(1+a/7*(1+a/8*(1+a/10)))))))))
              if(.not.lconv(5)) g(5) = rexp / a**6 * HLP9
     &                                               (1+a/10))))))))))
              if(.not.lconv(6)) g(6) = rexp / a**7 * HLP9
     &                                               (1+a/10*(1+a/11*(1+a/12))))))))))))
              if(.not.lconv(7)) g(7) = rexp / a**8 * HLP9
     &                                               (1+a/10*(1+a/11*(1+a/12*(1+a/13)))))))))))))
              do l = 8,2*lexp
                if(.not.lconv(l)) g(l) = a**(-(l+1))
              enddo
              lconv = lconv .or. g*factor*a**2<convpar
              do while(lconv(maxl)) ; maxl = maxl - 1 ; if(maxl<0) exit loop_ptsh ; enddo
            endif
            call harmonicsr(y,ra,maxl)
            do ikpt = 1,nkpt1
              cexp = exp( img * 2*pi * dot_product(kpt1(:,ikpt),ptsh(:,pnt(i))) )
              do l = 0,maxl ; if(lconv(l)) cycle
                l1               = l**2 + 1
                l2               = (l+1)**2
                shlp(l1:l2,ikpt) = shlp(l1:l2,ikpt) + g(l) * cexp * conjg(y(l1:l2))
              enddo
            enddo
          enddo loop_ptsh
          structconst(:,ic1,ic2,:) = shlp
        enddo
      enddo

      deallocate ( ptsh,radsh,pnt,shift )

      Rbegin
      call cpu_time(time2)
      if(lwrite) then ; write(6,'(A,F7.2)')       '  Timing: ',time2-time1
      else            ; write(6,'(A,F12.2'NoA) '    Timings: ',time2-time1
      endif
      call cpu_time(time1)
      if(lwrite) write(6,'(/A)')  'Fourier-space sum'
      Rend

c
c     Determine reciprocal shells
      nshift = nkpt1
      allocate ( shift(3,nshift) )
      do ikpt = 1,nkpt1
        shift(:,ikpt) = matmul(rlat, kpt1(:,ikpt) - nint(kpt1(:,ikpt)) )
      enddo
      call getshells(ptsh,nptsh,radsh,nshell,rrad,rlat,shift,nshift,lwrite)

      ! minimal nonzero reciprocal-shell radius (needed for algorithm 1 in routines exchange and selfenergy)
      radshmin = radsh(min(2,size(radsh)))

      call cpu_time(time1)

      allocate ( pnt(nptsh) )

c
c     Fourier-space sum
      latcon = rvol**(1d0/3)
      factor = 4*pi*scale**2/latcon**2
      pref   = 4*pi / ( scale**3 * vol )
      do ikpt = 1,nkpt1
        k = kpt1(:,ikpt) - nint(kpt1(:,ikpt)) ! -nint(...) transforms to the parallelepiped with the origin in the center ( i.e. -0.5 <= x,y,z < 0.5 )
        do i = 1,nptsh
          ki       = ptsh(:,i) + k
          ka       = matmul(rlat,ki)
          radsh(i) = sqrt(sum(ka**2))
        enddo
        call rorderp(pnt,radsh,nptsh)
        a1    = huge(a1)
        maxl  = min(7,lexp*2)
        lconv = .false.
        loop_rptsh: do i = 1,nptsh
          ki = ptsh(:,pnt(i)) + k
          ka = matmul(rlat,ki)
          a  = sqrt(sum(ka**2)) / scale
          if(abs(a-a1)>1d-10) then
            a1 = a
            aa = (1+a**2)**(-1)
            if(.not.lconv(2))   g(2) = pref * aa**5        / 3
            if(.not.lconv(3))   g(3) = pref * aa**5 * a    / 15
            if(.not.lconv(4))   g(4) = pref * aa**6 * a**2 / 105
            if(.not.lconv(5))   g(5) = pref * aa**6 * a**3 / 945
            if(.not.lconv(6))   g(6) = pref * aa**7 * a**4 / 10395
            if(.not.lconv(7))   g(7) = pref * aa**7 * a**5 / 135135
            if(a==0) then
              g(0) = pref * (-4)
              g(1) = 0
            else
              if(.not.lconv(0)) g(0) = pref * aa**4 / a**2
              if(.not.lconv(1)) g(1) = pref * aa**4 / a
              lconv(:7) = lconv(:7) .or. g(:7)*factor*a**2<convpar(:7)
            endif
            do while(lconv(maxl)) ; maxl = maxl - 1 ; if(maxl<0) exit loop_rptsh ; enddo
          endif
          call harmonicsr(y,ka,maxl)

          cdum = 1
          do l = 0,maxl
            if(.not.lconv(l)) then
              l1       = l**2 + 1
              l2       = (l+1)**2
              y(l1:l2) = conjg(y(l1:l2)) * g(l) * cdum
            endif
            cdum = cdum * img
          enddo

          Mpi( mm = -1 )
          do ic2 = 1,ncent
            do ic1 = 1,max(1,ic2-1) ; McycleP(mm)
              cexp = exp(img*2*pi*dot_product(ki,cent(:,ic1)-cent(:,ic2)))
              do l = 0,maxl ; if(lconv(l)) cycle
                l1                              = l**2 + 1
                l2                              = (l+1)**2
                structconst(l1:l2,ic1,ic2,ikpt) = structconst(l1:l2,ic1,ic2,ikpt) + cexp * y(l1:l2)
              enddo
            enddo
          enddo

        enddo loop_rptsh
      enddo
      Nfence(structconst)
      MpiO( call Msum(structconst,comm=Ocomm) )
      Nfence(structconst)

      deallocate ( ptsh,radsh,pnt,shift )

      Obegin
c     Add contribution for l=0 to diagonal elements
      structconst(1,1,1,:) = structconst(1,1,1,:) - 5d0/16/sqrt(4*pi)
c     Atom offdiagonals
      do ic2 = 1,ncent
        do ic1 = 1,ic2-1
          lm = 0
          do l = 0,2*lexp
            do m = -l,l
              lm                        = lm + 1
              structconst(lm,ic2,ic1,:) = (-1)**(l+m) * conjg(structconst(lm-2*m,ic1,ic2,:))
            enddo
          enddo
        enddo
      enddo
c     Atom diagonals
      do i = 2,ncent
        structconst(:,i,i,:) = structconst(:,1,1,:)
      enddo
c     Scale
      do l=0,lexp*2
        structconst(l**2+1:(l+1)**2,:,:,:) = structconst(l**2+1:(l+1)**2,:,:,:) * scale**(l+1)
      enddo
      Oend
      Nfence(structconst)

      Rbegin
      call cpu_time(time2)
      if(lwrite) then ; write(6,'(A,F7.2)') '  Timing: ',time2-time1
      else            ; write(6,'(F11.2)') time2-time1
      endif
c     Calculate accuracy of Gamma-decomposition
      if(nkpt1>1.and.all(kpt(:,1)==0).and.lwrite) then
        a = 1d30 ! ikpt = index of shortest non-zero k-point
        do i = 2,nkpt1
          rdum = sum(matmul(rlat,kpt(:,i))**2)
          if(rdum<a) then ; ikpt = i ; a = rdum ; endif
        enddo
        rdum = sqrt(sum(matmul(rlat,kpt(:,ikpt))**2))
        a    = 0
        do ic2=1,ncent
          do ic1=1,max(1,ic2-1)
            a = a + abs( structconst(1,ic1,ic2,ikpt) -
     &                 ( structconst(1,ic1,ic2,1) + sqrt(4*pi)/vol/rdum**2 *
     &                   exp(-img*2*pi*dot_product(kpt(:,ikpt),cent(:,ic2)-cent(:,ic1)))))**2
          enddo
        enddo
        a  = sqrt(a/ncent**2)
        aa = sqrt(sum(abs(structconst(1,:,:,ikpt))**2)/ncent**2)
        write(6,'(/A,F8.5,A,F8.5,A)') 'Accuracy of Gamma-decomposition (structureconstant):',a,' (abs)',a/aa,' (rel)'
      endif
      Rend

      contains

c     -----------------

c     Returns number of points in direct and reciprocal space for Ewald summation.
c     Also defines scale, convpar, rad, and rrad.
      function points(ew)
      use, intrinsic :: iso_fortran_env
      implicit none
      real_dp              :: points
      real_dp, intent(in)  :: ew
      real_dp              :: a,aa,rexp,pref,g(0:7),latcon,da,factor
      integer              :: l

c     Determine cutoff radii for real-space and Fourier-space summation

      ! (1) real space
      latcon = vol**(1d0/3)
      scale  = ew / latcon
# if   CONVTYPE == 0
      do l = 0,7      ; convpar(l) = CONVPARAM1                    ; enddo
      do l = 8,2*lexp ; convpar(l) = CONVPARAM2                    ; enddo
# elif CONVTYPE == 1
      do l = 0,7      ; convpar(l) = CONVPARAM1 / scale**(l+1)     ; enddo
      do l = 8,2*lexp ; convpar(l) = CONVPARAM2 / scale**(l+1)     ; enddo
# else
      do l = 0,7      ; convpar(l) = CONVPARAM1 / ew**(l+1) * 4*pi ; enddo
      do l = 8,2*lexp ; convpar(l) = CONVPARAM2 / ew**(l+1) * 4*pi ; enddo
# endif
      factor = 4*pi/scale**2/latcon**2
      a      = 1
      da     = 1d0
      do
        rexp = exp(-a)
        g(0) = rexp / a    * (1+a*11/16*(1+a*3/11*(1+a/9)))
        g(1) = rexp / a**2 * (1+a*(1+a/2*(1+a*7/24*(1+a/7))))
        g(2) = rexp / a**3 * (1+a*(1+a/2*(1+a/3*(1+a/4*(1+a*3/16*(1+a/9))))))
        g(3) = rexp / a**4 * (1+a*(1+a/2*(1+a/3*(1+a/4*(1+a/5*(1+a/6*(1+a/8)))))))
        g(4) = rexp / a**5 * (1+a*(1+a/2*(1+a/3*(1+a/4*(1+a/5*(1+a/6*(1+a/7*(1+a/8*(1+a/10)))))))))
        g(5) = rexp / a**6 * HLP9
     &                       (1+a/10))))))))))
        g(6) = rexp / a**7 * HLP9
     &                       (1+a/10*(1+a/11*(1+a/12))))))))))))
        g(7) = rexp / a**8 * HLP9
     &                       (1+a/10*(1+a/11*(1+a/12*(1+a/13)))))))))))))
        g    = g * factor*a**2 ; if(all(g<convpar(:7))) then ; if(da<5d-5) exit ; a = a - da ; da = da / 10 ; endif
        a    = a + da
      enddo
      rad    = a / scale
      points = 4*pi/3 * rad**3/vol

      ! (2) Fourier space
      latcon = rvol**(1d0/3)
      factor = 4*pi*scale**2/latcon**2
      pref   = 4*pi / ( scale**3 * vol )
      a      = 1
      da     = 1
      do
        aa   = (1+a**2)**(-1)
        g(0) = pref * aa**4 / a**2
        g(1) = pref * aa**4 / a
        g(2) = pref * aa**5        / 3
        g(3) = pref * aa**5 * a    / 15
        g(4) = pref * aa**6 * a**2 / 105
        g(5) = pref * aa**6 * a**3 / 945
        g(6) = pref * aa**7 * a**4 / 10395
        g(7) = pref * aa**7 * a**5 / 135135
        g    = g * factor*a**2 ; if(all(g<convpar(:7))) then ; if(da<5d-5) exit ; a = a - da ; da = da / 10 ; endif
        a    = a + da
      enddo
      rrad   = a*scale
      points = points + 4*pi/3 * rrad**3/rvol
      end function points

c     -----------------

      end

c     -----------------

c     Determines all shells of the crystal defined by lat and vol with radii smaller than rad.
c     The lattice points (number = nptsh) are stored in ptsh, their corresponding lengths (shell radii) in radsh.

c begin interface
      subroutine getshells(ptsh,nptsh,radsh,nshell,rad,lat,shift,nshift,lwrite)

      use util
# ifdef MPI      
      use global, only: Msize,Mrank,Mcomm,mpi_integer,mpi_sum
      use Mwrapper
# endif

      use, intrinsic :: iso_fortran_env !inc
      implicit none
      logical,              intent(in)  :: lwrite
      integer,              intent(in)  :: nshift
      integer,              intent(out) :: nptsh,nshell
      integer, allocatable, intent(out) :: ptsh(:,:)
      real_dp, allocatable, intent(out) :: radsh(:)
      real_dp,              intent(in)  :: rad,lat(3,3),shift(3,nshift)
c end interface
      real_dp                           :: r(3),rdum
      integer, allocatable              :: pnt(:)
      integer                           :: n,i,ix,iy,iz,ishift MpiC(Merr)
      logical                           :: found,define

      define = .false.
      i      = 0

 1    n      = ifMpi( Mrank , 0 )
      found  = .true.
      do while(found)
        found = .false.
        do ix = -n,n
          do iy = -(n-abs(ix)),n-abs(ix)
            iz   = n-abs(ix)-abs(iy)
 2          r    = ix*lat(:,1) + iy*lat(:,2) + iz*lat(:,3)
            rdum = minval ( [ (sum((r+shift(:,ishift))**2),ishift=1,nshift) ] )
            if(rdum<rad**2) then
              found = .true.
              i     = i + 1
              if(define) then
                ptsh(:,i) = [ ix,iy,iz ]
                radsh(i)  = sqrt(rdum)
              endif
            endif
            if(iz>0) then
              iz = -iz
              goto 2
            endif
          enddo
        enddo
        n = n + ifMpi( Msize , 1 )
      enddo

      if(.not.define) then
        nptsh = i
        i     = 0 ; Mpi( call mpi_scan(nptsh,i,1,mpi_integer,mpi_sum,Mcomm,Merr) ; i = i - nptsh ; call Msum(nptsh) )
        allocate ( radsh(nptsh),ptsh(3,nptsh) ) ; Mpi( radsh = 0 ; ptsh = 0 )
        define = .true.
        goto 1
      Mpi( else ; call Msum(radsh,0) ; call Msum(ptsh,0) )
      endif

      Rbegin
      allocate ( pnt(nptsh) )
      call rorderpf(pnt,radsh,nptsh)
      radsh  = radsh(pnt)
      ptsh   = ptsh(:,pnt)
      nshell = 1
      do i=2,nptsh
        if(radsh(i)-radsh(i-1)>1d-10) nshell = nshell + 1
      enddo
      if(lwrite)
     &  write(6,'(A,F10.5,A,I7,A,I6,A)') '  Sphere of radius',rad,' contains',
     &                                   nptsh,' lattice points and',nshell,' shells.'
      Rend

      Mpi( call Mcast(radsh) ; call Mcast(ptsh) ; call Mcast(nshell) )

      end

c     -----------------

c     Routines for multipole order of Coulomb matrix

c     -----------------      

c
c     Returns average deviation from zero of Coulomb matrix elements that should
c     vanish if MT functions are multipole-ordered (multipole-sparsity of Coulomb matrix).
c     dim      : rank of Coulomb matrix [usually nbasm(ikpt)]
c     mode = 0 : Coulomb matrix stored as upper triangle [coulomb(dim*(dim+1)/2)]
c     mode = 1 : Coulomb matrix stored as full matrix    [coulomb(dim,dim)]
      function error_coulomb_multipole(coulomb,dim,mode) result(err)
      use global, only: ntype,neq,lcutm,nindxm,nbasp
      use, intrinsic :: iso_fortran_env
      implicit none
      real_dp                 :: err
      integer,     intent(in) :: dim,mode
      MCOMPLEX_dp, intent(in) :: coulomb(*)
      integer                 :: itype,itype1,ieq,ieq1,l,l1,m,m1,n,n1,i,j,ij,ierr
      integer                 :: n_coulomb_multipole
      if(all(mode/=[0,1])) Bug('mode neither 0 nor 1')
      err  = 0
      ierr = 0
      ij   = 0
      j    = 0      
      do itype1 = 1,ntype
        do ieq1 = 1,neq(itype1)
          do l1 = 0,lcutm(itype1)
            do m1 = -l1,l1
              do n1 = 1,nindxm(l1,itype1)
                j = j + 1 ! coulumn index
                i = 0
                iloop: do itype = 1,ntype
                  do ieq = 1,neq(itype)
                    do l = 0,lcutm(itype)
                      do m = -l,l
                        do n = 1,nindxm(l,itype)
                          i  = i + 1  ! row index
                          ij = ij + 1 ! compound index
                          if((n/=1.or.n1/=1).and.(itype/=itype1.or.ieq/=ieq1.or.l/=l1.or.m/=m1)) then
                            err  = err + abs(coulomb(ij))**2
                            ierr = ierr + 1
                          endif
                          if(mode==0.and.i==j) exit iloop ! upper triangle criterion
                        enddo
                        if(mode==1) ij = ij + dim-nbasp
                      enddo
                    enddo
                  enddo
                enddo iloop
              enddo
            enddo
          enddo
        enddo
      enddo
      if(mode==0.and.ij/=nbasp*(nbasp+1)/2.or.mode==1.and.ij/=nbasp*dim) Bug('ij index count error.')
      do j = nbasp+1,dim ! column index
        i = 0
        do itype = 1,ntype
          do ieq = 1,neq(itype)
            do l = 0,lcutm(itype)
              do m = -l,l
                do n = 1,nindxm(l,itype)
                  i  = i + 1  ! row index
                  ij = ij + 1 ! compound index
                  if(n/=1) then
                    if(mode==0) then ; err = err + abs(coulomb(ij))**2                      
                    else             ; err = err + abs(coulomb((j-1)*dim+i))**2
                    endif
                    ierr = ierr + 1
                  endif
                enddo
              enddo
            enddo
          enddo
        enddo
        ij = ij + j-nbasp
      enddo
      err = sqrt(err/ierr)
      end

c     -----------------
      
c
c     mode=0 : Returns number of (generally) nonzero elements of upper triangle of Coulomb matrix if MT functions are multipole-ordered
c     mode=1 : Returns actual size of multipole-packed array as returned by pack_coulomb_multipole. Differs from number returned by mode=0
c              because diagonal elements of MT functions with nonzero multipoles are counted double (see reorder_coulomb_multipole)
c     mode=2 : Returns size of large block of multipole-packed array
c     mode=3 : Returns size of local blocks of multipole-packed array
c     dim    : Dimension (rank) of matrix
      function n_coulomb_multipole(dim,mode) result(n)
      use global, only: ntype,neq,lcutm,nindxm,nbasp
      use, intrinsic :: iso_fortran_env
      implicit none
      integer             :: n
      integer, intent(in) :: dim,mode
      integer             :: itype,l
      if(all(mode/=[0,1,2,3])) Bug('Unknown mode.')
      n = 0
      if(mode/=3) n = dim - nbasp + sum( [ (neq(itype) * (lcutm(itype)+1)**2, itype=1,ntype) ] )      
      n = n*(n+1)/2
      if(mode/=2) n = n + sum( [ ((neq(itype)*(2*l+1)*(nindxm(l,itype)*(nindxm(l,itype)+1)/2), l=0,lcutm(itype)),itype=1,ntype) ] )
      if(mode==0) n = n - sum( [ (neq(itype) * (lcutm(itype)+1)**2, itype=1,ntype) ] )
      end

c     -----------------

c
c     Front end routines for packing and unpacking Coulomb matrix according to multipole order (reorder_coulomb_multipole is called)
      subroutine pack_coulomb_multipole(coulomb,dim)
      use, intrinsic :: iso_fortran_env
      implicit none
      integer,     intent(in)    :: dim
      MCOMPLEX_dp, intent(inout) :: coulomb(dim*(dim+1)/2)
      call reorder_coulomb_multipole(coulomb,dim,1)
      end
      subroutine unpack_coulomb_multipole(coulomb,dim)
      use, intrinsic :: iso_fortran_env
      implicit none
      integer,     intent(in)    :: dim
      MCOMPLEX_dp, intent(inout) :: coulomb(dim*(dim+1)/2)
      call reorder_coulomb_multipole(coulomb,dim,0)
      end

c     -----------------

c
c     Multipole reordering
c
c     Back end routine for pack_coulomb_multipole and unpack_coulomb_multipole      
c     mode=0: Undoes reordering
c     mode=1: Reorders Coulomb matrix according to non-vanishing elements (multipole order):
c             coulomb(1:R)            : upper triangle matrix of large block: MT functions with nonzero multipoles (one per lm channel)
c                                       and interstital plane waves, R = n_coulomb_multipole(dim,2)
c             coulomb(R+1:R+N0)       : upper triangle matrix of l=0 local block of first atom, N0 = nindxm(0,1) * (nindxm(0,1)+1)/2
c             coulomb(R+N0+1:R+N0+N1) : ----- " -----    of l=1, m=-1    -----"-----          , N1 = nindxm(1,1) * (nindxm(1,1)+1)/2
c             etc.
c             coulomb(dim*(dim+1)/2) is set to huge(0d0) as flag for if a multipole-packed Coulomb matrix
c             Note that the diagonal elements of the MT functions with nonzero multipole are repeated and therefore set to zero
c             in the large block [coulomb(1:R)].      
      subroutine reorder_coulomb_multipole(coulomb,dim,mode)
      use global, only: ntype,neq,lcutm,nindxm,nbasp
      use, intrinsic :: iso_fortran_env
      implicit none
      integer,     intent(in)    :: dim,mode
      MCOMPLEX_dp, intent(inout) :: coulomb(dim*(dim+1)/2)
      MCOMPLEX_dp, allocatable   :: coulloc(:)
      integer                    :: dim_pack,dim_blk,dim_loc
      integer                    :: itype,itype1,ieq,ieq1,l,l1,m,m1,n,n1,i,j,ij,kl
      integer                    :: n_coulomb_multipole
      if(all(mode/=[0,1])) Bug('Unknown mode.')
      dim_pack = n_coulomb_multipole(dim,1) ! actual size of multipole-packed array
      dim_blk  = n_coulomb_multipole(dim,2) ! size of large block
      dim_loc  = dim_pack - dim_blk         ! size of local blocks
      if(dim_pack>=dim*(dim+1)/2) Error('Multipole-packed array exceeds array bounds.')
      allocate(coulloc(dim_loc))
      if(mode==0) then
        if(dble(coulomb(dim*(dim+1)/2))/=huge(0d0)) Bug('Matrix not multipole-ordered.') ! check multipole-order flag
        coulloc = coulomb(dim_blk+1:dim_pack)
        kl      = dim*(dim+1)/2 - dim_blk ! shift elements to the end to avoid conflict in copying
        do i = dim_blk,1,-1
          coulomb(kl+i) = coulomb(i)
        enddo
      else
        if(dble(coulomb(dim*(dim+1)/2))==huge(0d0)) Bug('Matrix multipole-ordered.')     ! check multipole-order flag
        kl = 0
        ij = 0        
        j  = 0
        do itype = 1,ntype
          do ieq = 1,neq(itype)
            do l = 0,lcutm(itype)
              do m = -l,l
                do n1 = 1,nindxm(l,itype)
                  ij = ij + j
                  do n = 1,n1
                    ij          = ij + 1
                    kl          = kl + 1
                    coulloc(kl) = coulomb(ij)
                  enddo
                enddo
                j = j + nindxm(l,itype)
              enddo
            enddo
          enddo
        enddo
        if(j/=nbasp)              Bug('Count error j/=nbasp.')
        if(ij/=nbasp*(nbasp+1)/2) Bug('Count error ij.')
        if(kl/=dim_loc)           Bug('Count error kl/=dim_loc.')
        kl = 0
      endif
      ! Large block: MT-MT with nonzero multipoles
      j  = 0
      ij = 0
      do itype1 = 1,ntype
        do ieq1 = 1,neq(itype1)
          do l1 = 0,lcutm(itype1)
            do m1 = -l1,l1
              do n1 = 1,nindxm(l1,itype1)
                j = j + 1 ! column index of unpacked matrix
                i = 0
                iloop: do itype = 1,ntype
                  do ieq = 1,neq(itype)
                    do l = 0,lcutm(itype)
                      do m = -l,l
                        do n = 1,nindxm(l,itype)
                          i  = i + 1  ! row index of unpacked matrix
                          ij = ij + 1 ! compound index for unpacked matrix
                          if(n==1.and.n1==1) then
                            kl = kl + 1 ! compound index for packed matrix
                            if(mode==0) then
                              coulomb(ij) = coulomb(kl)
                            else if(any([itype,ieq,l,m]/=[itype1,ieq1,l1,m1])) then
                              coulomb(kl) = coulomb(ij)
                            else
                              coulomb(kl) = 0
                            endif
                          else if(mode==0) then
                            coulomb(ij) = 0
                          endif
                          if(i==j) exit iloop                          
                        enddo
                      enddo
                    enddo
                  enddo
                enddo iloop
              enddo
            enddo
          enddo
        enddo
      enddo      
      if(any([i,j]/=nbasp))     Bug('Count error i,j.')
      if(ij/=nbasp*(nbasp+1)/2) Bug('Count error ij.')
      ! Large block: MT-PW matrix elements between IPWs and MT functions with nonzero multipoles
      do j = nbasp+1,dim ! column index of unpacked matrix
        do itype = 1,ntype
          do ieq = 1,neq(itype)
            do l = 0,lcutm(itype)
              do m = -l,l
                kl = kl + 1 ! compound index for packed matrix
                ij = ij + 1 ! compound index for unpacked matrix
                if(mode==0) then ; coulomb(ij) = coulomb(kl) ; coulomb(ij+1:ij+nindxm(l,itype)-1) = 0 ! unpack
                else             ; coulomb(kl) = coulomb(ij)                                          ! pack
                endif
                ij = ij + nindxm(l,itype) - 1 ! shift to last index of current lm channel
              enddo
            enddo
          enddo
        enddo
        ! PW-PW matrix elements
        if(mode==0) then ; coulomb(ij+1:ij+j-nbasp) = coulomb(kl+1:kl+j-nbasp)
        else             ; coulomb(kl+1:kl+j-nbasp) = coulomb(ij+1:ij+j-nbasp)
        endif
        kl = kl + j - nbasp
        ij = ij + j - nbasp
      enddo
      if(ij/=dim*(dim+1)/2)   Bug('Final compound index for unpacked matrix deviates from total number of elements.')
      if(mode==0) then
        if(kl/=dim*(dim+1)/2) Bug('Final compound index for packed matrix deviates from total number of elements.')
      else
        if(kl/=dim_blk)       Bug('Final compound index for packed matrix deviates from total number of block elements.')
      endif
      ! MT-MT local blocks
      if(mode==0) then
        kl = 0
        ij = 0        
        j  = 0
        do itype = 1,ntype
          do ieq = 1,neq(itype)
            do l = 0,lcutm(itype)
              do m = -l,l
                do n1 = 1,nindxm(l,itype)
                  ij = ij + j
                  do n = 1,n1
                    ij          = ij + 1
                    kl          = kl + 1
                    coulomb(ij) = coulloc(kl)
                  enddo
                enddo
                j = j + nindxm(l,itype)
              enddo
            enddo
          enddo
        enddo
        if(kl/=dim_loc)           Bug('Final compound index for packed matrix deviates from total number of local-block elements.')
        if(ij/=nbasp*(nbasp+1)/2) Bug('Final compound index for unpacked matrix deviates from total number of MT-MT elements.')        
      else
        coulomb(dim_blk+1:dim_pack) = coulloc
        coulomb(dim*(dim+1)/2)      = huge(0d0) ! set multipole-packed flag
      endif
      deallocate(coulloc)
      end

c     -----------------

c
c     Adds product mat*coulomb*mat to result(:)
c     result  : result matrix
c     coulomb : multipole-packed Coulomb matrix
c     mat     : input matrix (destroyed on exit!)
      subroutine matcoulombmat_multipole(result,coulomb,mat,dim1,dim2,fac)
      use global,  only: ntype,neq,lcutm,nindxm,nbasp
      use wrapper, only: matmat,macmat,unitarytrafo
      use, intrinsic :: iso_fortran_env
      implicit none
      integer,     intent(in)    :: dim1,dim2
      MCOMPLEX_dp, intent(inout) :: result(dim2,dim2)
      MCOMPLEX_dp, intent(inout) :: mat(dim1,dim2)
      MCOMPLEX_dp, intent(in)    :: coulomb(dim1*(dim1+1)/2)      
      real_dp,     intent(in)    :: fac
      MCOMPLEX_dp                :: help(dim2,dim2)
      integer                    :: itype,ieq,l,m,n,nn,i,k,kl,dim_blk
      integer                    :: n_coulomb_multipole
      if(coulomb(dim1*(dim1+1)/2)/=huge(0d0)) Bug('Coulomb matrix not multipole-packed.')
      dim_blk = n_coulomb_multipole(dim1,2)
      kl      = dim_blk
      k       = 0
      i       = 0
      do itype = 1,ntype
        do ieq = 1,neq(itype)
          do l = 0,lcutm(itype)
            n  = nindxm(l,itype)
            nn = n*(n+1)/2
            do m = -l,l
              result   = result + macmat ( mat(i+1:i+n,:) , matmat ( coulomb(kl+1:kl+nn) , mat(i+1:i+n,:) ) ) * fac
              k        = k + 1
              mat(k,:) = mat(i+1,:)
              i        = i + n
              kl       = kl + nn
            enddo
          enddo
        enddo
      enddo
      do i = nbasp+1,dim1
        k        = k + 1
        mat(k,:) = mat(i,:)
      enddo
      call unitarytrafo(help,coulomb(:dim_blk),mat(:k,:),1)
      result = result + help * fac  !macmat ( mat(:k,:) , matmat ( coulomb(:dim_blk) , mat(:k,:) ) ) * fac
      end

c     -----------------

c
c     Extracts submatrices from packed matrix coulomb(:) (assumed represented in multipole-ordered MT basis)
c     coulloc(*)   - local blocks in packed storage
c     coulmat(:,:) - large block in full storage
c     coulomb(:)   - Coulomb matrix in packed storage
c     dim          - number of mixed-basis functions
      subroutine coulomb_multipole(coulloc,coulmat,coulomb,dim)
      use global,  only: ntype,neq,lcutm,nindxm,nbasp,maxindxm
      use util,    only: chr
      use wrapper, only: unpackmat
      use, intrinsic :: iso_fortran_env
      implicit none
      integer,     intent(in)  :: dim
      MCOMPLEX_dp, intent(in)  :: coulomb(dim*(dim+1)/2)
      MCOMPLEX_dp, intent(out) :: coulloc(*)
      integer                  :: itype
      MCOMPLEX_dp, intent(out) :: coulmat( dim - nbasp + sum( [ (neq(itype) * (lcutm(itype)+1)**2, itype=1,ntype) ] ) ,
     &                                     dim - nbasp + sum( [ (neq(itype) * (lcutm(itype)+1)**2, itype=1,ntype) ] ) )
      integer                  :: ieq,l,m,n
      integer                  :: itype1,ieq1,l1,m1,n1
      integer                  :: i,i1,ii1,j,j1,jj1,kk1
      ! MT-MT
      i1  = 0
      j1  = 0
      ii1 = 0
      kk1 = 0
      do itype1 = 1,ntype
        do ieq1 = 1,neq(itype1)
          do l1 = 0,lcutm(itype1)
            do m1 = -l1,l1
              do n1 = 1,nindxm(l1,itype1)
                if(n1==1) j1 = j1 + 1 ! column index of coulmat(:,:)
                i1 = i1 + 1           ! column index of unpacked coulomb
                i  = 0
                j  = 0
                iloop: do itype = 1,ntype
                  do ieq = 1,neq(itype)
                    do l = 0,lcutm(itype)
                      do m = -l,l
                        do n = 1,nindxm(l,itype)
                          if(n==1) j = j + 1 ! row index of coulmat(:,:)
                          i   = i   + 1      ! row index of unpacked coulomb
                          ii1 = ii1 + 1      ! compound index of coulomb(:)
                          if(n==1.and.n1==1) then
                            coulmat(j,j1) =         coulomb(ii1)
                            coulmat(j1,j) = MCONJG( coulomb(ii1) )                            
                          endif
                          if(itype==itype1.and.ieq==ieq1.and.l==l1.and.m==m1) then
                            kk1 = kk1 + 1    ! compound index of coulloc
                            if(n==1.and.n1==1) then ; coulloc(kk1) = 0
                            else                    ; coulloc(kk1) = coulomb(ii1)
                            endif
                          endif
                          if(i==i1) exit iloop                          
                        enddo
                      enddo
                    enddo
                  enddo
                enddo iloop
              enddo
              n                       = nindxm(l,itype)
              kk1                     = kk1 - n*(n+1)/2
              coulloc(kk1+1:kk1+n**2) = reshape ( unpackmat( coulloc(kk1+1:kk1+n*(n+1)/2) ) , [n**2] )
              kk1                     = kk1 + n**2
            enddo
          enddo
        enddo
      enddo      
      if(any([i,i1]/=nbasp))     Bug('Count error i,i1.')
      if(ii1/=nbasp*(nbasp+1)/2) Bug('Count error ii1.')
      if(kk1/=sum( [ (( neq(itype)*(2*l+1) * nindxm(l,itype)**2, l=0,lcutm(itype)), itype=1,ntype) ] ) )
     &                           Bug('Count error kk1.')
      ! MT-PW
      do i1 = nbasp+1,dim ! column index of coulomb(:)
        j1 = j1 + 1       ! column index of coulmat
        j  = 0
        do itype = 1,ntype
          do ieq = 1,neq(itype)
            do l = 0,lcutm(itype)
              do m = -l,l
                j             = j   + 1 ! row index for coulmat
                ii1           = ii1 + 1 ! compound index for coulomb(:)
                coulmat(j,j1) =         coulomb(ii1)
                coulmat(j1,j) = MCONJG( coulomb(ii1) )
                ii1           = ii1 + nindxm(l,itype) - 1 ! shift to last index of current lm channel
              enddo
            enddo
          enddo
        enddo
        ! PW-PW matrix elements
        do i = 1,i1-nbasp
          j             = j   + 1
          ii1           = ii1 + 1
          coulmat(j,j1) =         coulomb(ii1)
          coulmat(j1,j) = MCONJG( coulomb(ii1) )
        enddo
      enddo
      if(ii1/=dim*(dim+1)/2) Bug('Count error ii1. (2)') 
      if(j/= dim - nbasp + sum( [ (neq(itype) * (lcutm(itype)+1)**2, itype=1,ntype) ] ) ) Bug('Count error j.')
      end
      
c     -----------------

c
c     Returns cprod1 = matmat(coulomb,cprod0)
c     Coulomb matrix is multipole-sparse and defined in coulloc and coulmat (as returned by coulomb_multipole)
c     - cprod0  :  input cprod (in standard order)
c     - cprod1  : output cprod (in standard order)
c     - coulloc : local blocks of multipole-sparse Coulomb matrix
c     - coulmat : large block  of multipole-sparse Coulomb matrix
c     - dim     : number of cprod vectors (second dimension of cprod0 and cprod1)
c     - ikpt    : kpoint index (if ikpt==1, the Coulomb matrix must not contain one term of the Coulomb spherical 
c                 average -- coeff*claplace + claplace*coeff -- that would break the multipole-sparsity of the
c                 Coulomb matrix and which is, therefore, taken into account explicitly)
      subroutine matmat_coulomb_multipole(cprod1,coulloc,coulmat,cprod0,ikpt,dim)
      use global,  only: ntype,neq,lcutm,nindxm,nbasp,fullpw,pi,nbasm,ngptm,grid,rgrid,vol,basm,gptm,pgptm
      use wrapper, only: matmat,matvec
      use timer_util
      use, intrinsic :: iso_fortran_env
      implicit none
      integer,     intent(in)  :: ikpt,dim
      MCOMPLEX_dp, intent(out) :: cprod1(nbasm(ikpt),dim)
      MCOMPLEX_dp, intent(in)  :: cprod0(nbasm(ikpt),dim),coulloc(*)
      integer                  :: dim1,itype,ieq,l,m,n,nn,i,j,k,il
      MCOMPLEX_dp, intent(in)  :: coulmat( ngptm(ikpt) + sum( [ (neq(itype)*(lcutm(itype)+1)**2,itype=1,ntype) ] ) ,
     &                                     ngptm(ikpt) + sum( [ (neq(itype)*(lcutm(itype)+1)**2,itype=1,ntype) ] ) )
      MCOMPLEX_dp, allocatable :: help(:,:)
      complex_dp,  allocatable :: coeff(:),claplace(:)
      MCOMPLEX_dp              :: mcdum,dot1,dot2
      MCOMPLEX_dp              :: stepfunction
      real_dp                  :: intgrf
      real_dp,     parameter   :: fac = 4*pi/6
c      call timer_start('exchange_mat')
c      call timer_start('exchange_help')
      dim1 = size(coulmat,1)
      allocate(help(dim1,dim))
      ! copy cprod for large block -> help
      do k = 1,dim
        i = 1 ! cprod index
        j = 0 ! help index
        do itype = 1,ntype
          do ieq = 1,neq(itype)
            do l = 0,lcutm(itype)
              do m = -l,l
                j         = j + 1
                help(j,k) = cprod0(i,k)
                i         = i + nindxm(l,itype)
              enddo
            enddo
          enddo
        enddo
        help(j+1:j+ngptm(ikpt),k) = cprod0(i:i+ngptm(ikpt)-1,k)
      enddo
c      call timer_stop('exchange_help')
      ! matrix multiplication of large block: coulmat*help -> cprod1
c      call timer_start('exchangem')
      cprod1(nbasm(ikpt)-dim1+1:,:) = matmat(coulmat,help)
c      call timer_stop('exchangem')
      ! redistribute cprod1 and add local blocks -> cprod1
      ! ikpt==1 : also add Coulomb spherical average fac*(coeff*claplace+claplace*coeff)
      deallocate(help)
c      call timer_start('exchange_redist')
      if(ikpt==1) then
        if(fullpw) then ; allocate(coeff(nbasp+1),claplace(nbasp+1))
        else            ; allocate(coeff(nbasm(1)),claplace(nbasm(1)))
        endif
        coeff    = 0
        claplace = 0
        i        = 0
        do itype = 1,ntype
          do ieq = 1,neq(itype)
            coeff(i+1)    =  sqrt(4*pi) * intgrf( rgrid(:,itype)    * basm(:,1,0,itype) , itype ) / sqrt(vol)
            claplace(i+1) = -sqrt(4*pi) * intgrf( rgrid(:,itype)**3 * basm(:,1,0,itype) , itype ) / sqrt(vol)
            !coeff(i+1)    =  sqrt( 4*pi * grid(itype)%radius**3 /  (3*vol)     ) ! alternatives
            !claplace(i+1) = -sqrt( 4*pi * grid(itype)%radius**7 / (25*vol) * 3 ) !
            do n = 2,nindxm(0,itype)
              claplace(i+n) = -sqrt(4*pi) * intgrf( rgrid(:,itype)**3 * basm(:,n,0,itype) , itype ) / sqrt(vol)
            enddo
            i = i + sum( [ ((2*l+1)*nindxm(l,itype),l=0,lcutm(itype)) ] )
          enddo
        enddo
        if(i/=nbasp) Bug('Count error: i.')
        if(fullpw) then ; coeff(nbasp+1) = 1d0
        else            ; do i = 1,ngptm(1) ; coeff(nbasp+i) = stepfunction(-gptm(:,pgptm(i,1))) ; enddo
        endif
        Inv( call symmetrize(coeff,   1,size(coeff),   2) )
        Inv( call symmetrize(claplace,1,size(claplace),2) )
      endif
      do k = 1,dim
        i  = 0 ! cprod index
        j  = 0 ! help index
        il = 0 ! coulloc index
        if(ikpt==1) then
          dot1 = 0
          dot2 = 0
        endif
        do itype = 1,ntype
          do ieq = 1,neq(itype)
            if(ikpt==1) then
              n    = nindxm(0,itype)
              dot1 = dot1 + sum( claplace(i+1:i+n) * cprod0(i+1:i+n,k) )
              dot2 = dot2 +         coeff(i+1)     * cprod0(i+1,k)
            endif
            do l = 0,lcutm(itype)
              n  = nindxm(l,itype)
              nn = n**2
              do m = -l,l
                j                 = j + 1                
                mcdum             = cprod1(nbasm(ikpt)-dim1+j,k)
                cprod1(i+1:i+n,k) = matmul( reshape(coulloc(il+1:il+nn),[n,n]) , cprod0(i+1:i+n,k) )
                cprod1(i+1,    k) = cprod1(i+1,k) + mcdum
                i                 = i  + n
                il                = il + nn
              enddo
            enddo
          enddo
        enddo
        if(ikpt==1) then
          dot1 = fac * dot1
          if(fullpw) then
            dot2               = fac * ( dot2 + cprod0(nbasp+1,k) )
            cprod1(:nbasp+1,k) = cprod1(:nbasp+1,k) - dot1 * MCONJG(coeff) - dot2 * MCONJG(claplace)
          else
            dot2               = fac * ( dot2 + sum(coeff(nbasp+1:)*cprod0(nbasp+1:,k)) )
            cprod1(:,k)        = cprod1(:,k)        - dot1 * MCONJG(coeff) - dot2 * MCONJG(claplace)
          endif
        endif
      enddo
      if(ikpt==1) deallocate(coeff,claplace)
c      call timer_stop('exchange_redist')
c      call timer_stop('exchange_mat')
      end
