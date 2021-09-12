c--------------------------------------------------------------------------------
c Copyright (c) 2021 Peter Grünberg Institut, Forschungszentrum Jülich, Germany
c This file is part of SPEX and available as free software under the conditions
c of the MIT license as expressed in the LICENSE file in more detail.
c--------------------------------------------------------------------------------

c     This routine is obsolete. The valence exchange is now calculated in the routine self-energy.
c
c     Calculates the valence contribution to the static (HF) exchange term:
c
c                                          s*       s           s*           s
c                                       phi  (r) phi     (r) phi     (r') phi  (r')
c                         occ.             nq       n'q+k       n'q+k        nq
c     exchange(n,q)  =  - SUM  INT INT  ------------------------------------------- dr dr'
c                         k,n                           | r - r' |
c
c                         occ                 s          s  ~        ~       s       s
c                    =  - SUM SUM  v     < phi      | phi   M    > < M    phi   | phi      >
c                         k,n I,J   k,IJ      n'q+k      nq  k,I      k,J    nq      n'q+k
c
c     At the Gamma point (k=0) v diverges. After diagonalization of v at k=0 the divergence is
c     restricted to the head element I=1. Furthermore, we expand <...> with kp perturbation theory.
c     As a result, the total I=1 element is given by a sum of a divergent 1/k**2-term and an
c     angular dependent term. The former is separated from the numerical k-summation and treated
c     analytically while the latter is spherically averaged and added to the k=0 contribution of
c     the numerical k-summation. (A better knowledge of the integrand's behavior at the BZ edges
c     might further improve the integration.)

c The divergence at the Gamma point is integrated with the following algorithm.
c     A periodic Function (same as Massidda's) with expo->0 is used:
c     The function  F(k) = lim(expo->0) SUM(G) exp(-expo*|k+G|**2) / |k+G|**2  is subtracted from
c     the BZ integral and integrated analytically. The contribution to the BZ integral including
c     the "tail" is
c     vol/(8*pi**3) INT F(k) d^3k - P SUM(k) F(k)  ( P = principal value ) .
c     For expo->0 the two terms diverge. Therefore a cutoff radius q0 is introduced and related to
c     expo by exp(-expo*q0**2)=delta  ( delta = small value, e.g., delta = 1d-10 ) .
c     The resulting formula
c     vol/(4*pi**1.5*sqrt(expo)) * erf(sqrt(a)*q0) - sum(q,0<q<q0) exp(-expo*q**2)/q**2
c     converges well with q0.

c MPI: Routine partly parallelized, but deactivated (getinput) because parallelization does not work yet.

# include "cppmacro.h"
# include "jobtype.h"
# include "restype.h"      

      subroutine exchange(job1)

      use global
      use arrays, only: selfx
      use util,   only: same_eigenspace,chr
      use wrapper
      use file
      use readwrite, only: read_vxc
      Mpi( use Mwrapper )

      use, intrinsic :: iso_fortran_env
      implicit none
      type(jobtype)            :: job1
      MCOMPLEX_dp, allocatable :: cprod(:,:,:,:),olap(:,:),cprod1(:),mat(:,:),help(:,:),moment(:,:,:)
      MCOMPLEX_dp              :: cprod0(maxbasm),mom(3)
      real_dp                  :: kvec(3),rdum
      real_dp                  :: print_pole,print_tail
      real_dp,     allocatable :: selfxc(:)
      complex_dp               :: hessian(3,3)
      complex_dp               :: cdum
      integer,     allocatable :: band(:),dimsub(:)
      integer                  :: ikpt2,iselfx
      integer                  :: l,n,nn,ib,deg0,deg1,iblock
      integer                  :: bandi,bandf,npack
      integer                  :: ispin,ibandq,i,j,k,ikpt,iband,ikptq,isub,jsub,di,dj,iband1,iband2
      integer                  :: nkpt1,kpt1(nkpt),nkpts(nkpt),nsym1,sym1(nsym),ikpt1,ndeg,isym
      logical                  :: occup(maxband)
      logical                  :: lkpt(nkpti+2)
      character(4)             :: cspin(2) = [ ' up ','down' ]
      logical                  :: ldum
      integer                  :: kptsum
      real_dp                  :: logdiv,error_coulomb_multipole
      real                     :: time1,time2,cputime,time_mat,time_cprod,time_trafo,time_gamma,time_vxc,time_tot

      if(all(job1%type/=[J_HF,J_GW,J_RPA,J_HFE,J_PBE0])) Bug('Wrong job type.')

# ifdef MPI
      RWarn('MPI run without NOSTORE: exchange self-energy parallelized only over k loop. Use NOSTORE for large systems.')
# endif

      time_cprod = 0
      time_mat   = 0
      time_gamma = 0
      time_trafo = 0
      time_vxc   = 0
      divergence = 0

c
c     Allocate array selfx
      if(job1%full) then
        i = 0
        do iblock = 1,size(block,2)
          j = count(block(:,iblock)/=0)
          i = i + j*(j+1)/2
        enddo
      else
        i = size(job1%band)
      endif
      Allocate_ ( selfx,(i) )
      selfx = 0

c
c     Read selfx from spex.sigx if requested (RESTART==2)
      Rbegin
      write(6,*)
      lkpt = .false.
      if     (iand(restart,R_sigx)/=0) then ; call read_dump(1,lkpt,nkpti+2,'spex.sigx',job1)
      else if(iand(restart,W_sigx)/=0) then ; call write_dump(1,[0],0,'spex.sigx',job1 MpiC(0) MpiC(0) ) ! truncate sigfile
      endif
      Rend
      Mpi( call Mcast(lkpt) )

c
c     Core contribution
      if(.not.lkpt(nkpti+2)) then
        call exchange_core(job1)
        Rif(iand(restart,W_sigx)/=0) call write_dump(1,[nkpti+2],1,'spex.sigx',job1 MpiC(0) MpiC(0) )
      endif

      if(all(lkpt(:nkpti))) then
        return
      else if(any(lkpt(:nkpti))) then
        Error('Cannot continue with calculation of exchange self-energy without NOSTORE. Remove spex.sigx or use NOSTORE.')
      endif

      Rbegin
      write(6,'(//A/)') '### subroutine: exchange ###'

      if(l_qsgw) write(6,'(A)') 'File qsgw exists and will be read.'
      Rend

      call gamma_divergence(.true.)

c      write(*,*) k0,'  Warning: k0 set to 0.403484954226802d0'
c      k0 = 0.403484954226802d0

      ! We skip the following section as it is not used.
# if 0
      ! Integrate 1/k**2 for 0 <= kx,ky,kz <= k0   -> integral  (NOT USED AT THE MOMENT)
      ! (1) Numerical part (integration with n*n*n k-points)
      n        = 200
      integral = 0
      a        = k0**2
      k(:3)    = rlat(:,3) / n
      rdum1 = 0
      do i = 1,n
        do j = 1,n
          kvec(:) = ( (i-0.5) * rlat(:,1) + (j-0.5) * rlat(:,2) + 0.5 * rlat(:,3) ) / n
          do l = 1,n
            k2 = kvec(1)**2 + kvec(2)**2 + kvec(3)**2
            if(k2<=a) then
              integral = integral + (3-2*sqrt(k2)/k0)/a
            else
              integral = integral + 1d0/k2
            endif
            rdum1 = rdum1 + 1d0/k2
            kvec(1) = kvec(1) + k(1)
            kvec(2) = kvec(2) + k(2)
            kvec(3) = kvec(3) + k(3)
          enddo
        enddo
      enddo
      integral = integral * rvol / (1d0*n)**3
      ! (2) Analytical part (use Neper's rules)
      a        = acos(dot_product(rlat(:,2),rlat(:,3))/sqrt(sum(rlat(:,2)**2)*sum(rlat(:,3)**2))) !
      b        = acos(dot_product(rlat(:,1),rlat(:,3))/sqrt(sum(rlat(:,1)**2)*sum(rlat(:,3)**2))) ! Angles between reciprocal lattice vectors
      c        = acos(dot_product(rlat(:,1),rlat(:,2))/sqrt(sum(rlat(:,1)**2)*sum(rlat(:,2)**2))) !
      aa       = acos((cos(a)-cos(b)*cos(c))/sin(b)/sin(c)) !
      bb       = acos((cos(b)-cos(c)*cos(a))/sin(c)/sin(a)) ! spherical angles
      cc       = acos((cos(c)-cos(a)*cos(b))/sin(a)/sin(b)) !
      rdum     = aa+bb+cc-pi ! spherical excess
      integral = integral + 2*pi*k0 * rdum/(4*pi)
# endif

c      write(*,*) 'SELFX SET TO ZERO!'
c      selfx = 0

c
c     Calculate the inverse of the overlap matrix for any k point
      if(fullpw) then
        Rwrite(6,'(A'NoA) 'Calculate inverse overlap matrices... '
        call cpu_time(cputime)
        call checkmem('olap',MBYTES*maxgptm*(maxgptm+1)/2*nkpti)
        allocate ( olap(maxgptm*(maxgptm+1)/2,nkpti) ) ; olap = 0
        do ikpt = Mrange1(nkpti)
          call olap_pwp(olap(1,ikpt),gptm(:,pgptm(:ngptm(ikpt),ikpt)),ngptm(ikpt))
          call inverse(olap(:ngptm(ikpt)*(ngptm(ikpt)+1)/2,ikpt))
        enddo
        Rcall cpu_done(cputime)
        Mpi( call Msum(olap) )
      endif

      call cpu_time(time_tot)

# ifdef MT_MULTIPOLE
      Obegin
      rdum = 0
      do ikpt = 2,nkpti
        n    = nbasm(ikpt)
        nn   = n*(n+1)/2
        rdum = max( rdum , error_coulomb_multipole(coulomb(:nn,ikpt),n) )
        call pack_coulomb_multipole(coulomb(:nn,ikpt),n)
      enddo
      write(6,'(A)') 'Average residual interaction between multipole-free MT functions: '//Chf(rdum,'ES8.1')
      if(rdum>1d-10) Error('Residual interaction too large: '//Chf(rdum,'ES8.1'))
      Oend
# endif

      Rbegin
      write(6,*)
      if(nspin==2) write(6,'(A'NoA) 'spin'
      write(6,'(A'NoA) '        k point       band      tail      pole     total     HF energy'
      if(job1%full) write(6,'(A'NoA) ' (diagonal)'
      write(6,*)
      Rend

      iselfx = 0

      do iblock = 1,size(block,2)

        ib = block(1,iblock) ! first band of block (or current band index of job1%band if only diagonal elements are calculated)

        if(.not.job1%full) then
          if(ib/=1) then ! if the current state is in the same eigenspace as the previous state, we can use the same exchange contribution
            if(ikptq==job1%kpt(ib).and.ispin==job1%spin(ib)) then
              i      = ibandq
              ibandq = job1%band(ib)
              if(same_eigenspace(ibandq,i,ikptq,ispin)) then
                iselfx        = iselfx + 1
                Rbegin
                allocate ( selfxc(1) )
                selfxc(1)     = selfx(iselfx)
                selfx(iselfx) = selfx(iselfx-1)
                Rend
                goto 3
              endif
            endif
          endif
        endif

        ibandq = job1%band(ib)
        ikptq  = job1%kpt(ib)
        ispin  = job1%spin(ib)

c
c       Get irreducible k points kpt1(:nkpt1) (wrt current k point ikptq)
        call getkpt1(kpt1,nkpt1,nkpts,sym1,nsym1,ikptq,0,.false.)
        if(nkpt1==nkpt) then
          call getkpt1_fullBZ(kpt1,nkpt1,nkpts,sym1,nsym1,0)
        endif

        n         = nband(ikptq,ispin)
        occup     = .false.
        occup(:n) = ene(:n,ikptq,ispin) <= efermi

c
c       Determine which wave-function products <M band|occ> must be calculated (->band)
        if(job1%full) then
          ! according to block (full self-energy calculations)
          i    = count(block(:,iblock)/=0) ; allocate ( band(i) )
          band = job1%band(block(:i,iblock))
        else
          ! Determine first and last state of degenerate subspace (only diagonal self-energy elements)
          deg0 = ibandq
          deg1 = deg(ibandq,ikptq,ispin)
          if(deg1<deg0) then
            deg0 = deg1
            deg1 = deg(deg0,ikptq,ispin)
          endif
          ndeg = deg1 - deg0 + 1
          ! double spin degeneracies (SOC+inv.sym.) are not due to spatial symmetry and don't have to be averaged over: leave out!
          if(l_soc.and.invsym/=0.and.ndeg==2) then
            ndeg = 1
            deg1 = deg0
          endif
          ! define array band accordingly
          if(ndeg*nkpt1>=nkpt) then ! in this case summing over all k points is faster
            allocate ( band(1) )
            band = ibandq
            deg0 = ibandq
            deg1 = ibandq
            call getkpt1_fullBZ(kpt1,nkpt1,nkpts,sym1,nsym1,0) ! sum over the full BZ
          else
            allocate ( band(ndeg) )
            band = [ (i,i=deg0,deg1) ]
          endif
        endif

        ! store core exchange terms for printing tail and pole contributions below
        Rbegin
        allocate ( selfxc(size(band)) )
        j = 0
        do i = 1,size(band)
          j         = j + i
          selfxc(i) = selfx(iselfx+j)
          if(.not.job1%full) exit
        enddo
        Rend

        allocate ( mat(size(band),size(band)) )
        mat = 0

c
c       Divide into packets of occupied states if memory (MEM) is exceeded and loop over packets
        ifR rdum = ( maxmem - mem ) Mpi(/Nsize) ; Mpi( call Mcast(rdum) )
        npack    = int(rdum/(1d0*MBYTES*maxbasm*nkpt1*size(band)))
        do bandi = 1,bando,npack
          bandf = min(bandi+npack-1,bando)

c
c         Calculate wave-function products
          call cpu_time(time1)
          allocate ( cprod(maxbasm,size(band),bandi:bandf,Mcol1(nkpt1)) )
          call wavefproducts1_mt(cprod,maxbasm,band,size(band),ikptq,ispin,kpt1(Mcol1(nkpt1)),Mnum1(nkpt1),bandi,bandf,1)
          call wavefproducts1_pw(cprod,maxbasm,band,size(band),ikptq,ispin,kpt1(Mcol1(nkpt1)),Mnum1(nkpt1),bandi,bandf,1)
          call cpu_time(time2) ; time_cprod = time_cprod + time2 - time1

c
c         Numerical contribution

          ! Loop over EIBZ
          do ikpt1 = Mrange1(nkpt1)
            ikpt  = kpt1(ikpt1)
            ikpt2 = kptsum(ikptq,ikpt)

            ! The Coulomb matrix is only defined in the IBZ. We have to transform either (a) the Coulomb matrix
            ! from kptp(ikpt) to ikpt or (b) the vector cprod from ikpt back to kptp(ikpt) which is faster.
            ! Comment out one of the following!
            call cpu_time(time1)
            ! (a) (disabled)
c            allocate ( coul1(nbasm(ikpt),nbasm(ikpt)),coul2(nbasm(ikpt),nbasm(ikpt)) )
c            coul1 = unpackmat(coulomb(:nbasm(ikpt)*(nbasm(ikpt)+1)/2,kptp(ikpt)))
c            call matrixtrafo(coul2,coul1,ikpt,sym(symkpt(ikpt))%inv,.true.)
            ! (b)
            ! Transforming the vectors (case fullpw and time-reversal symmetry P: IBZ->BZ, P MI=M'I):
            ! (summation over I,J,K,L is implicit)
            !   SUM(n=occ) <n|j M'I> <M'I|M'J>^(-1) <M'J|v|M'K> <M'K|M'L>^(-1) <M'L i|n>
            ! = SUM(n=occ) <n|j P MI> <P MI|P MJ>^(-1) <P MJ|v|P MK> <P MK|P ML>^(-1) <P ML i|n>
            ! = SUM(n=occ) <n|j P MI> conjg [ <MI|MJ>^(-1) <MJ|v|MK> <MK|ML>^(-1) ] <P ML i|n>
            ! = conjg { SUM(n=occ) conjg[<n|j P MI>] <MI|MJ>^(-1) <MJ|v|MK> <MK|ML>^(-1) conjg[<P ML i|n>] }
            ! Note that we have to take the transpose instead of conjg if the operator v is complex (e.g., the screened interaction W):
            ! = SUM(n=occ) <P ML i|n> <ML|MK>^(-1) <MK|W(w)|MJ> <MJ|MI>^(-1) <n|j P MI>
            ! = SUM(n=occ) conjg[<n|i P ML>] <ML|MK>^(-1) <MK|W(w)|MJ> <MJ|MI>^(-1) conjg[<P MI j|n>]
            n = ngptm(ikpt)*(ngptm(ikpt)+1)/2
            do i = bandi,bandf
              if(lomit(1)) then
                if(any(omit==i)) cprod(:,:,i,ikpt1) = 0
              endif
              do j = 1,size(band)
                call mtrafo Inv(_r) (cprod0,cprod(:,j,i,ikpt1),nbasm(ikpt),1,kptp(ikpt),-symkpt(ikpt),1,.false.)
                if(fullpw) cprod0(nbasp+1:nbasm(ikpt)) = matvec(olap(:n,kptp(ikpt)),cprod0(nbasp+1:nbasm(ikpt)))
                cprod(:,j,i,ikpt1) = cprod0
              enddo
            enddo
            call cpu_time(time2) ; time_trafo = time_trafo + time2 - time1

            n  = nbasp + ngptm(ikpt)
            nn = n*(n+1)/2

            ! Calculate <i|selfx|j> = SUM(n=occ) <n|j MI> <MI|v|MJ> <MJ i|n>
            allocate ( cprod1(n) )
            do iband = bandi,bandf
              if(wintgr(ikpt2,iband,ispin)<1d-10) cycle ! skip if weight is tiny
# ifdef MT_MULTIPOLE
              if(ikpt>1) then
                rdum = wintgr(ikpt2,iband,ispin) * nkpts(ikpt1) / nsym1
                call matcoulombmat_multipole(mat,coulomb(:nn,kptp(ikpt)),cprod(:n,:,iband,ikpt1),n,size(band),rdum)
              else
# endif             
                do j = 1,size(band)
                  ! matrix-vector product: cprod1 = coulomb * cprod(j)
                  cprod1 = matvec ( coulomb(:nn,kptp(ikpt)) , cprod(:n,j,iband,ikpt1) )
                  do i = 1,j
                    ! scalar products: conjg(cprod1) * cprod(i)
                    cdum = dotprod ( cprod1 , cprod(:n,i,iband,ikpt1) ) *
     &                     wintgr(ikpt2,iband,ispin) * nkpts(ikpt1) / nsym1
                    if(i/=j.and.symkpt(ikpt)>nsymt) cdum = conjg(cdum)
                    mat(i,j) = mat(i,j) + cdum
                    mat(j,i) = MCONJG(mat(i,j))
                  enddo
                enddo
# ifdef MT_MULTIPOLE
              endif
# endif              
            enddo
            deallocate ( cprod1 )
            call cpu_time(time1) ; time_mat = time_mat + time1 - time2

          enddo

          deallocate ( cprod )

        enddo ! end loop over packets

c
c       Gamma point:
c       (A) Zero-order terms [ O(k^2) / k^2 ] (diagonal and offdiagonal elements) (-> "tail" in the output)
c       (B) Analytic integration over 1/k**2  (diagonal elements)                 (-> "pole" in the output)
        ifMOD(iblock)
        call cpu_time(time1)
        ! Calculate momentum matrix <nq|-i\nabla|n'q>
        if(ozero) then
          allocate ( moment(nband(ikptq,ispin),3,size(band)) )
          do i = 1,size(band)
            ibandq = band(i)
            call momentum_matrix(moment(:,:,i),[ikptq],1,ispin,ispin,ibandq,ibandq,1,nband(ikptq,ispin),.false. MpiC(.false.) )
          enddo
        endif
        ! Loop over self-energy matrix
        do j = 1,size(band) ; iband2 = band(j) ; if(lomit(1)) then ; if(any(iband2==omit)) cycle ; endif
          if(ozero) then
            do i = 1,j      ; iband1 = band(i)
              ! Calculate Hessian
              hessian = 0
              do iband = 1,nband(ikptq,ispin)
                if(   occup(iband1) .eqv. occup(iband2) ) then ! occ/occ or unocc/unocc
                  if( occup(iband1) .eqv. occup(iband) ) cycle
                  rdum = ( ene(iband,ikptq,ispin) - ene(iband1,ikptq,ispin) ) *
     &                   ( ene(iband,ikptq,ispin) - ene(iband2,ikptq,ispin) ) * vol
                  if( occup(iband1) ) rdum = -rdum
                else                                           ! occ/unocc
                  rdum = ( ene(iband2,ikptq,ispin) - ene(iband1,ikptq,ispin) ) * vol
                  if( occup(iband) ) then ; rdum = rdum * ( ene(iband,ikptq,ispin) - ene(iband2,ikptq,ispin) )
                  else                    ; rdum = rdum * ( ene(iband,ikptq,ispin) - ene(iband1,ikptq,ispin) )
                  endif
                endif
                do k = 1,3
                  do l = 1,3
                    hessian(k,l) = hessian(k,l) + moment(iband,k,i) * MCONJG(moment(iband,l,j)) / rdum
                  enddo
                enddo
              enddo
              ! (A) Add spherical average [ trace(Hessian)/3 ] to array mat
              cdum     = 4*pi * (hessian(1,1)+hessian(2,2)+hessian(3,3))/3 / nkpt / nsym1
              mat(i,j) = mat(i,j) + cdum
              mat(j,i) = MCONJG(mat(i,j))
            enddo
          endif
          ! (B) 1/k**2 term
          if(ologdiv) then
            if(ozero) then ; mom = moment(iband2,:,iband2)
            else           ; call momentum_matrix(mom,[ikptq],1,ispin,ispin,iband2,iband2,iband2,iband2,.false. MpiC(.false.) )
            endif
            mat(j,j) = mat(j,j) + 4*pi/vol * divergence / nsym1 * logdiv(ene(iband2,ikptq,ispin)-efermi,dble(mom))
          else if(occup(iband2)) then
            mat(j,j) = mat(j,j) + 4*pi/vol * divergence / nsym1
          endif
        enddo
        if(allocated(moment)) deallocate ( moment )
        call cpu_time(time2) ; time_gamma = time_gamma + time2 - time1
        endMOD

c
c       Add equivalent k points to array mat
        if(use_sym.and.job1%full) then
          ! look for subspaces in current block (dimsub(i)=dimension of subspace i)
          allocate ( dimsub(nsub) )
          dimsub = 0
          do i = 1,size(band)
            isub         = psub(block(i,iblock))
            dimsub(isub) = dimsub(isub) + 1
          enddo
          ! multiply with irreps
          if(nkpt1/=nkpt) then
            allocate ( help(size(band),size(band)) )
            help = mat
            mat  = 0
            i    = 0
            do isub = 1,nsub
              di = dimsub(isub) ; if(di==0) cycle
              j  = 0
              do jsub = 1,nsub
                dj = dimsub(jsub) ; if(dj==0) cycle
                do isym = 1,nsym1
                  if(sym1(isym)>nsymt) then
                    mat(i+1:i+di,j+1:j+dj) = mat(i+1:i+di,j+1:j+dj) + conjg (
     &                                       matmul (           transpose(conjg(irrep_sub(:di,:di,isub,sym1(isym)))) ,
     &                                       matmul ( help(i+1:i+di,j+1:j+dj) , irrep_sub(:dj,:dj,jsub,sym1(isym)) ) ) )
                  else
                    mat(i+1:i+di,j+1:j+dj) = mat(i+1:i+di,j+1:j+dj) +
     &                                       matmul (           transpose(conjg(irrep_sub(:di,:di,isub,sym1(isym)))) ,
     &                                       matmul ( help(i+1:i+di,j+1:j+dj) , irrep_sub(:dj,:dj,jsub,sym1(isym)) ) )
                  endif
                enddo
                j = j + dimsub(jsub)
              enddo
              i = i + dimsub(isub)
            enddo
            deallocate ( help )
          endif
          deallocate ( dimsub )
        endif

        Mpi( call Msum(mat,0) )

c
c       Add matrix elements (mat) to array selfx
        Rbegin
        if(job1%full) then
          do j = 1,size(band)
            do i = 1,j
              iselfx        = iselfx + 1
              selfx(iselfx) = selfx(iselfx) - mat(i,j)
            enddo
          enddo
        else
          ! use the great orthogonality theorem for diagonal elements
          iselfx        = iselfx + 1
          selfx(iselfx) = selfx(iselfx) - sum( [ (mat(i,i),i=1,size(band)) ] ) / size(band) * nsym1
        endif
        Rend

c
c       Print diagonal elements
 3      Rbegin
        do i = 1,size(block,1)
          ibandq = block(i,iblock) ; if(ibandq==0) exit
          ibandq = job1%band(ibandq)

          ldum = lomit(1)
          if(ldum) ldum = ldum.and.any(omit==ibandq)
          if(occup(ibandq).and..not.ldum) then ; print_pole = -4*pi/vol * divergence * hartree
          else                                 ; print_pole = 0
          endif

          ! Index of diagonal element
          if(job1%full) then ; j = iselfx - size(band)*(size(band)+1)/2 + i*(i+1)/2
          else               ; j = iselfx
          endif

          print_tail = (selfx(j)-selfxc(i))*hartree - print_pole

          ! k-point vector
          kvec = kpt(:,ikptq)

          ! Output
          allocate(help(1,1)) ; help = 0
          call cpu_time(time1)
          if(ovxc==0) then ; call read_vxc(help,[ibandq],1,kptp(ikptq),ispin,l_qsgw)
          else             ; call calc_vxc(help,[ibandq],1,kptp(ikptq),ispin,l_qsgw)
          endif
          call cpu_time(time2) ; time_vxc = time_vxc + time2 - time1
          rdum = help(1,1)
          deallocate(help)
          if(any(job1%type==[J_RPA,J_HFE])) rdum = rdum + selfxc(i)/2 ! subtract core-val exchange once because it is doubled for energy calculations (see below)
          if(nspin==1) write(6,'(   ''  ('',F5.3,'','',F5.3,'','',F5.3,'')'',I4,1X,3F10.5,F14.7)')
     &                   kvec,ibandq,print_tail,print_pole,print_tail+print_pole,
     &      (ene(ibandq,ikptq,ispin)+dble(selfx(j))-rdum)*hartree
          if(nspin==2) write(6,'(A4,''  ('',F5.3,'','',F5.3,'','',F5.3,'')'',I4,1X,3F10.5,F14.7)')
     &      cspin(ispin),kvec,ibandq,print_tail,print_pole,print_tail+print_pole,
     &      (ene(ibandq,ikptq,ispin)+dble(selfx(j))-rdum)*hartree
        enddo

        if(job1%full.and.iblock/=size(block,2))
     &    write(6,'(A)') '------------------------------------------------------------------------'

        Rend

        if(allocated(band))   deallocate ( band )
        if(allocated(mat))    deallocate ( mat )
        if(allocated(selfxc)) deallocate ( selfxc )

      enddo

      if(fullpw) then
        call checkmem('olap',-MBYTES*maxgptm*(maxgptm+1)/2*nkpti)
        deallocate (olap)
      endif

# ifdef MT_MULTIPOLE
      Obegin
      do ikpt = 2,nkpti
        n  = nbasm(ikpt)
        nn = n*(n+1)/2
        call unpack_coulomb_multipole(coulomb(:nn,ikpt),n)
        if(error_coulomb_multipole(coulomb(:nn,ikpt),n)>0) Error('Coulomb error nonzero.')
      enddo          
      Oend
# endif

      Rbegin

c
c     Write exchange self-energy to file spex.sigx
      Rif(iand(restart,W_sigx)/=0) call write_dump(1,[(k,k=1,nkpti)],nkpti,'spex.sigx',job1 MpiC(0) MpiC(0) )

      call cpu_time(time2) ; time_tot = time2 - time_tot

      write(6,'(29X,A)')   'Timings'
      write(6,'(A,F13.5)') 'Wave-function products:',time_cprod
      write(6,'(A,F13.5)') 'Matrix-vector products:',time_mat
      write(6,'(A,F13.5)') 'Transformations:       ',time_trafo
      if(ozero)
     &write(6,'(A,F13.5)') 'Gamma-point correction:',time_gamma
      write(6,'(A,F13.5)') 'vxc potential:         ',time_vxc
# ifdef MPI
      write(6,'(A,F13.5)') 'MPI idle time:         ',time_tot-(time_cprod+time_mat+time_trafo+time_gamma)
# endif      
      write(6,'(A)')       '------------------------------------'
      write(6,'(A,F13.5)') 'Total:                 ',time_tot

      Rend

      end

c     ----------------------------

c     Calculates the core contribution to the static (HF) exchange term. The only convergence
c     parameter is lcutm. (Low cutoffs, e.g. 4, are enough for convergence.)
c     Otherwise the calculation is exact (i.e. all core-valence product functions are
c     calculated irrespective of linear dependencies).
c
c     About (a) SOC-averaged and (b) SOC-splitted core states:
c     - The summation over the degnerate core states leads to a factor (a) 2l+1 and (b) l+1 for j=l+1/2 and l for j=l-1/2
c       (b: in the m-independent spin-diagonal contribution).
c     - Case (a): the core exchange potential is spherical symmetric in real space and does not lift degeneracies.
c     - Case (b): the core exchange potential is spherical symmetric in real/spin space and does not lift degeneracies provided
c       that SOC is also used for the valence states.
c     - In the inconsistent case (case (b) and no SOC for valence), the m-dependent terms of the core exchange potential will
c       break the spatial spherical symmetry, though. Therefore, these terms are neglected in this case.
c
      subroutine exchange_core(job1)

      use global
      use arrays, only: selfx
      Mpi( use Mwrapper )
      Load( use readwrite )

      use, intrinsic :: iso_fortran_env
      implicit none
      type(jobtype)            :: job1
      real_dp                  :: fprim1(maxgrid),fprim2(maxgrid)
      real_dp,     allocatable :: fprod(:,:),integral(:,:)
      complex_dp,  allocatable :: c1(:),c2(:),c3(:),d1(:),d2(:),d3(:)
      MCOMPLEX_dp, allocatable :: selfx1(:)
      integer,     allocatable :: band(:)
      integer                  :: icent,ikpt,ikindx,iband,ispin,ispinc,itype,ineq,ibnd,jbnd
      integer                  :: l,m,l1,n1,l2,m2,n2,lm,mm,i,j,ib,iselfx,iselfx0
      integer                  :: ifac,iblock MpiC(iloop)
      real_dp                  :: kp(3),qfac,qfac0,qfac1,qfac2,pfac
      real_dp                  :: gaunt,intgrf
      character(4)             :: cspin(2) = [ ' up ','down' ]
      real                     :: time1,time2

      if(all(job1%type/=[J_HF,J_GW,J_RPA,J_HFE,J_SX,J_COSX,J_PBE0])) Bug('Wrong job type.')

      if(maxindxc==0) then
        Rwrite(6,'(/A)') 'No core states.'
        return
      endif

      Rbegin

      call cpu_time(time1)
      write(6,'(//A)') '### subroutine: exchange_core ###'
      if(nspin==1) write(6,'(/A)')     '        k point       band      exchange (core contribution)'
      if(nspin==2) write(6,'(/A)') 'spin        k point       band      exchange (core contribution)'

      Rend

      qfac1   = 0
      iselfx0 = 0

      do iblock = 1,size(block,2)

        ib    = block(1,iblock) ! first band
        iband = job1%band(ib)
        ikpt  = job1%kpt(ib)
        ispin = job1%spin(ib)

        if(l_soc.or.lcore_soc) then ; ispinc = 1
        else                        ; ispinc = ispin
        endif

        if(job1%full) then
          i = count(block(:,iblock)/=0)
          allocate ( band(i),selfx1(i*(i+1)/2) )
          band = job1%band(block(:i,iblock))
        else
          allocate ( band(1),selfx1(1) )
          band = iband
        endif
        selfx1 = 0

        kp     = kpt(:,ikpt)
        ikindx = kindx(ikpt)
        Load( allocate ( cmt(maxlmindx,ncent,size(band),1,nspin3) ) )
        Load( call read_wavef2(band,size(band),ikpt,ispin,cmt) ; ikindx = 1 )
        Mpi( iloop = -1 )

        icent = 0
        do itype = 1,ntype
          do ineq = 1,neq(itype)
            icent = icent + 1
            do l1 = 0,lcutc(itype)
              do n1 = 1,nindxc(l1,itype)

                McycleP(iloop)

                if(lcore_soc.and.l1/=0) then
                  if(mod(n1,2)==1) then ; ifac = l1 + 1
                  else                  ; ifac = l1
                  endif
                else
                  ifac = 2*l1+1
                endif

                do l2 = 0,lcut(itype)

                  do l = abs(l1-l2),l1+l2,2!min(l1+l2,lcutm(itype)),2

                    pfac = 4*pi / ((2*l1+1)*(2*l+1))

                    allocate ( fprod(maxgrid,nindx(l2,itype)) )
                    allocate ( integral(nindx(l2,itype),nindx(l2,itype)) )
                    allocate ( c1(nindx(l2,itype)) )
                    allocate ( d1(nindx(l2,itype)) )
                    if(l_soc) then
                      allocate ( c2(nindx(l2,itype)) )
                      allocate ( d2(nindx(l2,itype)) )
                      if(lcore_soc.and.l2/=0.and.l1/=0) then
                        allocate ( c3(nindx(l2,itype)) )
                        allocate ( d3(nindx(l2,itype)) )
                      endif
                    endif

                    ! Define core-valence product functions
                    do n2 = 1,nindx(l2,itype)
                      fprod(:,n2) = ( core1(:,n1,l1,itype,ispinc)*bas1(:,n2,l2,itype,ispin) +
     &                                core2(:,n1,l1,itype,ispinc)*bas2(:,n2,l2,itype,ispin) ) / rgrid(:,itype)
                    enddo

                    ! Evaluate radial integrals (special part of Coulomb matrix : contribution from single MT)
                    do i = 1,nindx(l2,itype)
                      call primitivef(fprim1, fprod(:,i) * rgrid(:,itype)**(l+1), itype)
                      call primitivef(fprim2, fprod(:,i) / rgrid(:,itype)**l    ,-itype)
                      fprim1 = fprim1 / rgrid(:,itype)**l
                      fprim2 = fprim2 * rgrid(:,itype)**(l+1)
                      do j = 1,nindx(l2,itype)
                        integral(i,j) = intgrf( fprod(:,j) * (fprim1 + fprim2), itype)
                      enddo
                    enddo

                    qfac0 = pfac * gaunt(l1,l2,l,0,0,0) * sqrt((2*l1+1)*(2*l+1)/(4*pi*(2*l2+1))) * ifac
c                    qfac0 = pfac * sum([ ( gaunt(l1,l2,l,m,0,-m)**2, m=-l1,l1) ]) * ifac
                    if(l_soc.and.lcore_soc) then
                      qfac1 = pfac * sum([ ( gaunt(l1,l2,l,m,1,1-m)**2*m, m=-l1,l1) ])
                      qfac2 = pfac * sum([ ( gaunt(l1,l2,l,m,0,-m)*gaunt(l1,l2,l,m+1,1,-m)*sqrt((l1+1d0+m)*(l1-m)), m=-l1,l1-1) ])
                      if(ifac==l1) then
                        qfac1 = -qfac1
                        qfac2 = -qfac2
                      endif
                    endif

                    ! Add everything up
                    lm = sum( [ ((2*i+1)*nindx(i,itype),i=0,l2-1) ] )
                    do m2 = -l2,l2
                      iselfx = 0
                      do j = 1,size(band) ; jbnd = ifLoad(j,band(j))
                        c1 = cmt(lm+1:lm+nindx(l2,itype),icent,jbnd,ikindx, ifLoad(1,ispin) )
                        d1 = matmul(integral,c1)
                        if(l_soc) then
                          c2 = cmt(lm+1:lm+nindx(l2,itype),icent,jbnd,ikindx,2)
                          d2 = matmul(integral,c2)
                          if(lcore_soc.and.m2/=l2.and.l1/=0) then
c                            if(size(band)>1) call fatal('never tested: CORESOC+FULL')
                            c3 = cmt(lm+nindx(l2,itype)+1:lm+2*nindx(l2,itype),icent,jbnd,ikindx,2)
                            d3 = matmul(integral,c3)
                          endif
                        endif
                        do i = 1,j ; ibnd = ifLoad(i,band(i))
                          iselfx        = iselfx + 1
                          c1            = cmt(lm+1:lm+nindx(l2,itype),icent,ibnd,ikindx, ifLoad(1,ispin) )
                          selfx1(iselfx) = selfx1(iselfx) - dot_product(c1,d1) * ( qfac0 + m2*qfac1 ) ! SOC: up/up
                          if(l_soc) then
                            c2            = cmt(lm+1:lm+nindx(l2,itype),icent,ibnd,ikindx,2)
                            selfx1(iselfx) = selfx1(iselfx) - dot_product(c2,d2) * ( qfac0 - m2*qfac1 ) ! SOC: down/down
                            if(lcore_soc.and.m2/=l2.and.l1/=0) then
                              mm            = ( abs(2*m2+1) + 1 ) / 2
                              qfac          = sqrt((l2+mm)*(l2+1d0-mm)/(l2*(l2+1))) * qfac2
                              c3            = cmt(lm+nindx(l2,itype)+1:lm+2*nindx(l2,itype),icent,ibnd,ikindx,2)
                              selfx1(iselfx) = selfx1(iselfx) - qfac * ( dot_product(c1,d3) + dot_product(c3,d1) ) ! SOC: up/down + down/up
                            endif
                          endif
                        enddo
                      enddo
                      lm = lm + nindx(l2,itype)
                    enddo

                    deallocate(integral,c1,d1,fprod)
                    if(allocated(c2)) deallocate ( c2 )
                    if(allocated(c3)) deallocate ( c3 )
                    if(allocated(d2)) deallocate ( d2 )
                    if(allocated(d3)) deallocate ( d3 )

                  enddo

                enddo
              enddo
            enddo
          enddo
        enddo

        if(any(job1%type==[J_RPA,J_HFE])) selfx1 = selfx1 * 2 ! add cor*val*val*cor, which is identical to val*cor*cor*val calculated above

        Mpi ( call Msum(selfx1) )
        Rbegin
        selfx(iselfx0+1:iselfx0+size(band)*(size(band)+1)/2) = selfx(iselfx0+1:iselfx0+size(band)*(size(band)+1)/2) + selfx1
        do i = 1,size(band)
          iband = band(i)
          if(job1%full) then ; j = iselfx - size(band)*(size(band)+1)/2 + i*(i+1)/2
          else               ; j = iselfx
          endif
          if(nspin==1) write(6,'(   ''  ('',F5.3,'','',F5.3,'','',F5.3,'')'',I4,1X,F14.7)')
     &                                kp,iband,dble(selfx1(j)) * hartree
          if(nspin==2) write(6,'(A4,''  ('',F5.3,'','',F5.3,'','',F5.3,'')'',I4,1X,F14.7)')
     &                   cspin(ispin),kp,iband,dble(selfx1(j)) * hartree
        enddo
        if(job1%full.and.iblock/=size(block,2)) write(6,'(A)') '------------------------------------------'
        Rend

        iselfx0 = iselfx0 + size(band)*(size(band)+1)/2

        deallocate ( band,selfx1 )
        Load( deallocate ( cmt ) )

      enddo

      call cpu_time(time2)
      Rwrite(6,'(A,F7.2)') 'Timing:',time2-time1

      end

c     ----------------------------

      function exchange_energy_core()
      use global
      use, intrinsic :: iso_fortran_env
      implicit none
      real_dp :: exchange_energy_core
      real_dp :: fprod(maxgrid),fprim1(maxgrid),fprim2(maxgrid)
      real_dp :: result,integral
      integer :: ispin,itype
      integer :: l1,l2,l,m1,m2,m,n1,n2
      real_dp :: gaunt,intgrf  ,rdum

      exchange_energy_core = 0
      if(maxindxc==0) return
      if(lcore_soc) Error('CORESOC not implemented yet.')

      result = 0

      do ispin = 1,nspin
        do itype = 1,ntype

          do l1 = 0,lcutc(itype)   ; if(nindxc(l1,itype)==0) cycle
            do l2 = 0,lcutc(itype) ; if(nindxc(l2,itype)==0) cycle
              do l = abs(l1-l2),l1+l2,2

                ! Evaluate radial integrals and sum (-> integral)
                integral = 0
                do n1 = 1,nindxc(l1,itype)
                  do n2 = 1,nindxc(l2,itype)
                    fprod = ( core1(:,n1,l1,itype,ispin)*core1(:,n2,l2,itype,ispin) +
     &                        core2(:,n1,l1,itype,ispin)*core2(:,n2,l2,itype,ispin) ) / rgrid(:,itype)
                    call primitivef(fprim1, fprod * rgrid(:,itype)**(l+1), itype)
                    call primitivef(fprim2, fprod / rgrid(:,itype)**l    ,-itype)
                    fprim1   = fprim1 / rgrid(:,itype)**l
                    fprim2   = fprim2 * rgrid(:,itype)**(l+1)
                    integral = integral + intgrf( fprod * (fprim1 + fprim2), itype)
                  enddo
                enddo

                ! Check if identity used below is correct (can be removed later)
                rdum = 0
                do m1 = -l1,l1
                  do m2 = -l2,l2
                    m    = m1 - m2 ; if(abs(m)>l) cycle
                    rdum = rdum + gaunt(l,l1,l2,m,m1,m2)**2 * (4*pi)/(2*l+1)
                  enddo
                enddo
                if(abs(rdum-sqrt(4*pi*(2*l1+1)*(2*l2+1)/(2*l+1)) * gaunt(l,l1,l2,0,0,0))>1d-14) Error('!')

                ! Sum up HF energy (-> result)
                result = result + sqrt(4*pi*(2*l1+1)*(2*l2+1)/(2*l+1)) * gaunt(l,l1,l2,0,0,0) * integral * neq(itype)

              enddo
            enddo
          enddo

        enddo
      enddo

      exchange_energy_core = - result / nspin2

      end

