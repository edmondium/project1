c--------------------------------------------------------------------------------
c Copyright (c) 2021 Peter Grünberg Institut, Forschungszentrum Jülich, Germany
c This file is part of SPEX and available as free software under the conditions
c of the MIT license as expressed in the LICENSE file in more detail.
c--------------------------------------------------------------------------------

c Prepares calculation of offdiagonal self-energy matrix elements.
c
c The list of band indices is decomposed into subspaces of degenerate states.
c For each subspace the corresponding irrep is calculated.
c The array block contains the band indices of equivalent subspace that form the blocks of the self-energy matrix.
c
c The following constants/arrays are defined:
c nsub                     - number of subspaces
c psub(i)                  - index of the subspace that contains band i (in job1%band)
c irrep_sub(:,:,isub,isym) - irreducible representation for subspace isub and symmetry operation isym
c block(:,i)               - band indices (from job1%band) for block i
c
c Definition of irrep_sub:
c irrep_sub(i,j,isub,isym) = < i | P(isym) | j >
c
c Note: All irreps are unitary, but
c       matmul ( irrep_sub(:,:,isub,isym), irrep(:,:,isub,sym(isym)%inv) ) = identity
c       is not always valid:
c       (1) P(isym) involves time reversal. Then, use conjg(transpose(irrep(:,:,isub,sym(isym)%inv))) instead.
c       (2) Spin-orbit coupling. Rotation about 360Â° yields a factor -1 so that the right-hand side should be -identity.
c
c Uncomment the following to avoid automatic increase of irrep_thr
c # define NO_AUTO_INCREASE

# include "cppmacro.h"

      subroutine prepare_offdiag(job1)

      use global
      use util
      Mpi( use Mwrapper )
      Load( use readwrite )

      use, intrinsic :: iso_fortran_env
      implicit none

      type(jobtype), intent(inout) :: job1
      MCOMPLEX_dp,   allocatable   :: olappw(:,:)
      MCOMPLEX_dp,   allocatable   :: cpwhlp(:,:,:),cpwhlp1(:,:,:)
      complex_dp,    allocatable   :: cmthlp(:,:,:,:),cmthlp1(:,:,:,:),cmat(:,:)
      complex_dp                   :: cdum,ch(nsym),cphase
      real_dp,       allocatable   :: mat(:,:)
      real_dp,       parameter     :: irrep_threshold = 1d-8 ! choose 0 for k-point blocks instead of irrep blocks
      real_dp                      :: rdum,maxerr,off(nsym),irrep_thr
      logical,       allocatable   :: help(:),lmat(:,:),todo(:)
      integer,       allocatable   :: band(:),kpoint(:),spins(:),dsub(:),chsame(:),chdiff(:),psub1(:)
      integer                      :: nsym1
      integer                      :: ib1,ib2,nb,ikpt,ispin,isub,jsub,iloop,maxb,i,j,k,n,n1,ib,isym
      complex_dp                   :: wfolap,wfolap_mt
      real                         :: cputime
      logical                      :: isnan
      Load( logical                :: storeibz1 )

      irrep_thr = irrep_threshold

      if(irrep_thr<0) then
        Bug('Parameter irrep_thr negative.')
      else if(irrep_thr==0) then
        Info('irrep_thr is zero. K-point blocks used instead of irrep blocks.')
      endif

      if(trsoff) then ; nsym1 = nsym
      else            ; nsym1 = nsym/2
      endif

c     Dummy definition if only diagonal elements are calculated (define a block for each band index in job1%band containing only that band index)
      if(.not.job1%full) then
        nblock     = size(job1%band) ; allocate ( block(1,nblock),sizeblock(nblock) )
        block(1,:) = [ (i,i=1,nblock) ]
        sizeblock  = 1
        return
      endif

c     Dummy definition for use_sym=.false. (define blocks for each kpoint)
      if(.not.use_sym) then
        do
          ib    = 0
          n1    = 0
          ikpt  = 0
          ispin = 0
          do i = 1,size(job1%band)
            if(job1%kpt(i)/=ikpt.or.job1%spin(i)/=ispin) then
              ikpt  = job1%kpt(i)
              ispin = job1%spin(i)
              ib    = ib + 1
              n     = 0
            endif
            n  = n + 1
            n1 = max(n,n1)
            if(allocated(block)) then ; block(n,ib) = i ; sizeblock(ib) = n ; endif
          enddo
          if(allocated(block)) then ; exit
          else                      ; nblock = ib ; allocate ( block(n1,nblock), sizeblock(nblock) )
          endif
        enddo
        return
      endif

      Load( storeibz1 = storeibz ; storeibz  = .false. )

c     Make sure that the bands in job1%band include all degenerate states (if not, add the corresponding states)
c     (iloop = 1: count needed JOB states; iloop = 2: define states; then, redefine job1 if necessary)
      iloop = 1
      n     = size(job1%band)
 1    n1    = 0
      i     = 1
      do while(i<=n)
        ikpt  = job1%kpt(i)
        ispin = job1%spin(i)
        ib1   = job1%band(i)
        ib2   = deg(ib1,ikpt,ispin)
        if(ib2<ib1) then
          ib1 = ib2
          ib2 = deg(ib1,ikpt,ispin)
        endif
        do ib = ib1,ib2
          n1 = n1 + 1
          if(iloop==2) then
            band(n1)   = ib
            kpoint(n1) = ikpt
            spins(n1)  = ispin
            do k = 1,n1-1
              if(all([band(k),kpoint(k),spins(k)]==[ib,ikpt,ispin]))
     &          Error('Band definition double or not ordered.')
            enddo
          endif
        enddo
        do while(job1%kpt(i)==ikpt.and.job1%spin(i)==ispin.and.job1%band(i)>=ib1.and.job1%band(i)<=ib2) ! fast forward the subsequent degenerate states
          i = i + 1 ; if(i>n) exit
        enddo
      enddo
      if(n1>n) then
        if(iloop==1) then
          allocate ( band(n1),kpoint(n1),spins(n1) )
          iloop = 2
          goto 1
        endif
        Rwrite(6,'(/A)') 'Degenerate subspaces not complete. Additional states will be included.'
        deallocate ( job1%band,job1%kpt,job1%spin )
        allocate ( job1%band(n1),job1%kpt(n1),job1%spin(n1) )
        job1%band = band
        job1%kpt  = kpoint
        job1%spin = spins
        deallocate ( band,kpoint,spins )
      else
        Rwrite(6,'(/A)') 'Degenerate subspaces complete.'
      endif

c     Calculate irreps (iloop=1: count irreps, iloop=2: calculate irreps)
      Rwrite(6,'(A'NoA) 'Calculate irreps... '
      call cpu_time(cputime)
      allocate ( psub(size(job1%band)) )
      maxerr = 0
      iloop  = 1
      maxb   = 0
      psub   = 0
 2    isub   = 0
      do ispin = 1,nspin1
        do ikpt = 1,nkpt2
          n = count(job1%kpt==ikpt.and.job1%spin==ispin)
          if(n/=0) then
            allocate ( band(n),psub1(n) )
            band = pack(job1%band,job1%kpt==ikpt.and.job1%spin==ispin)
# ifdef LOAD
            if(iloop==2) then
              allocate ( cmt(maxlmindx,ncent,n,ikpt:ikpt,nspin3) )
              allocate ( cpw(maxgpt,         n,ikpt:ikpt,nspin3) )
              call read_wavef2(band,size(band),ikpt,ispin,cmt,cpw)
              kindx(ikpt) = ikpt
            endif
# endif
            i = 1
            do while(i<=n)
              ib1  = band(i)
              ib2  = deg(ib1,ikpt,ispin)
              nb   = ib2 - ib1 + 1
              isub = isub + 1
              maxb = max(maxb,nb)  ! max. dimension of irrep
              psub1(i:i+nb-1) = isub
              ifMOD(isub)
              if(iloop==2) then
                ! Calculate irreps
                allocate ( cmthlp(maxlmindx,ncent,nb,nspin3),cpwhlp(maxgpt,nb,nspin3) )
                allocate ( cmthlp1(maxlmindx,ncent,nb,nspin3),cpwhlp1(maxgpt,nb,nspin3) )
                allocate ( olappw(ngpt(ikpt),ngpt(ikpt)) )
                call olap_gpt(olappw,ngpt(ikpt),ikpt,ikpt)
                do isym = 1,nsym
                  if(kptsym(ikpt,isym)/=ikpt) cycle
                  Inv ( cphase = exp(-img * 2*pi * dot_product(kpt(:,ikpt),sym(isym)%transl) ) )
                  do j = 1,nb
# ifndef old_trafo
                    call wavefunction(cmthlp1(:,:,j,:),cpwhlp1(:,j,:),     j-1+ ifLoad(i,ib1) ,ikpt, ifLoad(1,ispin) )
                    call waveftrafo  (cmthlp (:,:,j,:),cpwhlp (:,j,:),isym,j-1+ ifLoad(i,ib1) ,ikpt, ifLoad(1,ispin) )
# else
                    call wavefunction(cmthlp1(:,:,j,1),cpwhlp1(:,j,1),          j-1+ ifLoad(i,ib1) ,ikpt, ifLoad(1,ispin) )
                    call waveftrafo  (cmthlp (:,:,j,1),cpwhlp (:,j,1),ikpt,isym,j-1+ ifLoad(i,ib1) ,      ifLoad(1,ispin) )
                    if(l_soc) then
                      call wavefunction(cmthlp1(:,:,j,2),cpwhlp1(:,j,2),          j-1+ ifLoad(i,ib1) ,ikpt, 2)
                      call waveftrafo  (cmthlp (:,:,j,2),cpwhlp (:,j,2),ikpt,isym,j-1+ ifLoad(i,ib1) ,      2)
                      call waveftrafo_soc(cmthlp(:,:,j,:),maxlmindx*ncent,isym)
                      call waveftrafo_soc(cpwhlp(:,j,:),maxgpt,isym)
                    endif
# endif
                  enddo
                  do j = 1,nb
                    do k = 1,nb
                      irrep_sub(k,j,isub,isym) = wfolap( cmthlp1(:,:,k,1),cpwhlp1(:,k,1),
     &                                                   cmthlp(:,:,j,1),cpwhlp(:,j,1),ikpt,ikpt,ispin,
     &                                                   cdum,olappw ) Inv( * cphase )
                      if(l_soc) then
                        irrep_sub(k,j,isub,isym) = irrep_sub(k,j,isub,isym) +
     &                                             wfolap( cmthlp1(:,:,k,2),cpwhlp1(:,k,2),
     &                                                     cmthlp(:,:,j,2),cpwhlp(:,j,2),ikpt,ikpt,nspin,
     &                                                     cdum,olappw ) Inv( * cphase )
                      endif
                    enddo
                  enddo
                enddo
                deallocate ( cmthlp,cpwhlp,cmthlp1,cpwhlp1,olappw )
              endif
              endMOD
              i = i + nb
            enddo
            if(iloop==2) then
              ! Determine whether matrix elements (e.g., for the self-energy) can possibly nonzero
              allocate ( cmthlp(maxlmindx,ncent,n,nspin3),cpwhlp(maxgpt,1,1) ) ! cpwhlp dummy array
              do i = 1,n
# ifndef old_trafo
                call wavefunction_mt(cmthlp(:,:,i,:),maxlmindx,0, ifLoad(i,band(i)) ,ikpt, ifLoad(1,ispin) )
# else
                if(storeibz) then
                  call wavefunction(cmthlp(:,:,i,1),cpwhlp,band(i),ikpt,ispin)
                  if(l_soc)
     &            call wavefunction(cmthlp(:,:,i,2),cpwhlp,band(i),ikpt,2)
                else
                  cmthlp(:,:,i,1) = cmt(:,:, ifLoad(i,band(i)) ,ikpt, ifLoad(1,ispin) )
                  if(l_soc)
     &            cmthlp(:,:,i,2) = cmt(:,:, ifLoad(i,band(i)) ,ikpt,2)
                endif
# endif
              enddo
              do i = 1,n ; Mcycle(i)
                do j = 1,n
                  cdum =        wfolap_mt(cmthlp(:,:,j,1),cmthlp(:,:,i,1),ispin,cdum,.true.)
                  if(l_soc)
     &            cdum = cdum + wfolap_mt(cmthlp(:,:,j,2),cmthlp(:,:,i,2),ispin,cdum,.true.)
                  mat(psub1(i),psub1(j)) = mat(psub1(i),psub1(j)) + abs(cdum)
                enddo
              enddo
              deallocate ( cmthlp,cpwhlp )
              Load( deallocate(cmt,cpw) )
            endif
            if(iloop==1) psub = unpack(psub1,job1%kpt==ikpt.and.job1%spin==ispin,psub)
            deallocate ( band,psub1 )
          endif
        enddo
      enddo
      if(iloop==1) then
        iloop = 2
        nsub  = isub
        allocate ( irrep_sub(maxb,maxb,nsub,nsym) )
        allocate ( mat(nsub,nsub) )
        allocate ( dsub(nsub) )
        do i = 1,nsub
          dsub(i) = count(psub==i)
        enddo
        irrep_sub = 0
        mat       = 0
        goto 2
      endif
      Mpi( call Msum(irrep_sub); call Msum(mat) )

c     Calculate maximal deviation from unitarity
      do isym = 1,nsym
        do isub = 1,nsub ; Mcycle(isub+(isym-1)*nsym)
          if(any(irrep_sub(:,:,isub,isym)/=0)) then
            nb   = dsub(isub) ; allocate ( cmat(nb,nb) )
            cmat = matmul ( irrep_sub(:nb,:nb,isub,isym) , conjg(transpose(irrep_sub(:nb,:nb,isub,isym))) )
            do j = 1,nb
              cmat(j,j) = cmat(j,j) - 1
            enddo
            rdum   = sum(abs(cmat)) ; deallocate ( cmat )
            maxerr = max(maxerr,rdum)
            if(isnan(rdum)) Bug('Detected NaN in irreps.')
            if(rdum>1d-10*nb**2*ncent) Error('Broken unitarity '//Chf(rdum,'ES8.1')//' for sub/sym: '//Chr(isub)//'/'//Chr(isym))
          endif
        enddo
      enddo

c     Calculate non-symmetry of overlap
      rdum = 0
      do j = 1,nsub
        do i = 1,j
          rdum     = rdum + abs(mat(i,j)-mat(j,i))
          mat(j,i) = mat(i,j)
        enddo
      enddo

      Rbegin
      call cpu_done(cputime)
      write(6,'(A,ES8.1)') 'Maximal deviation from unitarity: ',maxerr
      write(6,'(A,ES8.1)') 'Non-symmetry of overlap:          ',rdum

c     Prepare lists of band indices that define the nonzero matrix blocks of the self-energy (-> block)
c     ( iloop = 1: define dimension ; iloop = 2: define array block )
      allocate ( help(nsub),todo(nsub),lmat(nsub,nsub) )
      iloop = 1
 3    lmat  = mat>irrep_thr
      todo  = .true.
      n1    = 0
      n     = 0
      do i = 1,nsub
        help = lmat(:,i)
        do
          nb = count(help)
          do j = 1,nsub
            if(help(j)) help = help .or. lmat(:,j)
          enddo
          if(nb==count(help)) exit
        enddo
        help = help .and. todo
        if(count(help)>0) then
          n = n + 1
          k = 0
          do j = 1,nsub
            if(help(j)) then
              todo(j) = .false.
              k       = k + 1 ; if(iloop==2) block(k,n) = j
            endif
          enddo
          n1 = max(n1,k)
        endif
      enddo
      if(iloop==1) then
        allocate ( block(n1,n) )
        block = 0
        iloop = 2
        goto 3
      endif
      if(any(todo)) Error('Did not find all irreps from MT overlap.')
      nb = 0
      do i = 1,n
        nb = max( nb , sum ( [ (count(psub==block(j,i)),j=1,n1) ] ) )
      enddo
      if(nb<n1) Bug('nb < n1.')
      call reallocate(block,nb,n)
      allocate(band(nb))
      do i = 1,n
        band = 0
        ib   = 0
        do j = 1,n1
          do k = 1,size(job1%band)
            if(psub(k)==block(j,i)) then
              ib       = ib + 1
              band(ib) = k
            endif
          enddo
        enddo
        block(:,i) = band
      enddo
      deallocate(band)

      if(irrep_threshold==0) goto 5 ! Skip characters if k-point blocks (instead of irrep blocks) are used

c     Check characters
      ! loop over blocks
      do i = 1,size(block,2)
        isub = psub(block(1,i))
        nb   = dsub(isub)
        k    = job1%kpt(block(1,i))
        ! define characters
        do isym = 1,nsym1
          if(kptsym(k,isym)/=k) cycle
          ch(isym) = sum ( [ (irrep_sub(j,j,isub,isym),j=1,nb) ] )
        enddo
        ! check if characters are identical for the whole block.
        allocate ( chsame(size(block,1)),chdiff(size(block,1)) )
        n  = 1 ; chsame(1) = isub
        n1 = 0
        do j = nb+1,size(block,1),nb
          if(block(j,i)==0) exit
          jsub   = psub(block(j,i))
          maxerr = 0
          if(dsub(jsub)/=nb) Warn('Subgroups in current block have different dimensions. This will fail...')
          do isym = 1,nsym1
            if(kptsym(k,isym)/=k) cycle
            rdum   = abs(ch(isym)-sum([(irrep_sub(j,j,jsub,isym),j=1,nb)]))
            maxerr = max(maxerr,rdum)
          enddo
          if(maxerr<1d-8) then
            n          = n + 1
            chsame(n)  = jsub
          else
            n1         = n1 + 1
            chdiff(n1) = jsub
            if(rdum>1d-8) Warn('Character for symmetry '//Chr(isym)//' of subgroup '//Chr(jsub)//)
     &                         ' deviates from that of parent subgroup.'
          endif
        enddo
        ! If not, increase irrep_thr and go back to redefine blocks.
        if(n1/=0) then
          write(0,'(/A)') 'Deviation of symmetry characters detected for block '//Chr(i)
          write(0,'(A)')  'You might have picked up an accidental degeneracy. (Unlikely but possible.)'
          write(0,'(A)')  'Spex assumes all degeneracies to be caused by symmetry. Currently, there is no other solution '//
     &                    'than changing the kpoint in the JOB definition.'
          if(k==nkpt+1) write(0,'(A)') 'The problem occurred for the additional (+) kpoint. '//
     &                                 'Changing this kpoint slightly should solve the problem.'
          write(0,*)
          rdum = maxval ( [ ((mat(chsame(j),chdiff(k)),j=1,n),k=1,n1) ] )
# ifdef NO_AUTO_INCREASE
          Warn('Deviation in characters (block '//Chr(i)//');'//)
     &                              'irrep threshold should be increased to '//Chr(rdum)
          write(0,'(A)')            '        But NO_AUTO_INCREASE is set. Calculation is continued.'
# else
          irrep_thr = rdum
          if(irrep_thr>1d-4) then ; Error('Deviation in characters; too large irrep threshold increase: '//Chf(irrep_thr,'ES8.1'))
          else                    ; Warn('Deviation in characters; irrep threshold increased to '//Chf(irrep_thr,'ES8.1'))
          endif
          deallocate ( chsame,chdiff,block )
          iloop = 1
          goto 3
# endif
        endif
        deallocate ( chsame,chdiff )
      enddo

      deallocate ( help,lmat,todo )

      if(irrep_thr/=irrep_threshold)
     &  write(6,'(A,ES8.1)') 'Irrep threshold increased to:     ',irrep_thr
      write(6,'(A,I4)')      'Number of self-energy blocks:     ',size(block,2)
      write(6,'(A,I4)')      'Maximum number of bands per block:',size(block,1)

c     Write character table
      write(6,'(//A)') '- Character table -'
      allocate ( help(48) )
      do i = 1,size(block,2)
        isub = psub(block(1,i))
        nb   = dsub(isub)
        k    = job1%kpt(block(1,i))
        help = .false.
        write(6,'(/A,I4,A)') 'Block',i,' ('//Chr(k)//')'
        do isym = 1,nsym1
          if(kptsym(k,isym)/=k) cycle
          ch(isym)   = sum ( [ (irrep_sub(j,j,isub,isym),j=1,nb) ] )
          off(isym)  = sum(abs(irrep_sub(:nb,:nb,isub,isym))) - sum([(abs(irrep_sub(j,j,isub,isym)),j=1,nb)])
          help(isym) = .true.
        enddo
        n = 1
 4      if(off(n)>1d-10) then ; write(6,'(''*'''NoA)
        else                  ; write(6,'('' '''NoA)
        endif
        write(6,'(2F9.5,'' : '''NoA) ch(n)
        j = 0
        do isym = n,nsym1
          if(help(isym).and.abs(ch(isym)-ch(n))<1d-8) then
            write(6,'(I3'NoA) isym
            help(isym) = .false.
            j          = j + 1 ; if(j==16) exit
          endif
        enddo
        write(6,*)
        do while(.not.help(n))
          n = n + 1 ; if(n==nsym1+1) exit
        enddo
        if(n<=nsym1) goto 4
      enddo

      deallocate ( help,mat,dsub )

# if 0
      ! Topological invariants (not implemented yet)
      if(l_soc.and..not.trsoff.and.all(mod(nkpt3,2)==0)) then
        allocate ( dsub(8) ) ! \delta=Pf(irrep)/sqrt(det(irrep)) for the eight TRIMS
        n = 0
        do i = 0,1 ; do j = 0,1 ; do k = 0,1
          n    = n + 1
          ikpt = pkpt(i*nkpt3(1)/2,j*nkpt3(2)/2,k*nkpt3(3)/2) ; if(ikpt>nkpti) cycle
          ...
        enddo
      endif
# endif

      write(6,'(/A)') 'Band assignment to irreps (blocks)'
      write(6,'(A)')  '  Blk kpt spn bnd'
      do i = 1,size(block,2)
        ikpt  = job1%kpt(block(1,i))
        ispin = job1%spin(block(1,i))
        write(6,'(A'NoA) trim(chr(i,'I5'))//trim(chr(ikpt,'I4'))//trim(chr(ispin,'I4'))
        do j = 1,size(block,1)          
          if(block(j,i)/=0) then
            ib = job1%band(block(j,i))
            write(6,'(A'NoA) trim(chr(ib,'I4'))
          endif
        enddo
        write(6,*)
      enddo

 5    continue      

      Rend
      Mpi ( call Mcastl(block) )

      nblock = size(block,2) ; allocate(sizeblock(nblock))
      do i = 1,nblock
        sizeblock(i) = count(block(:,i)/=0)
      enddo

      Load( storeibz = storeibz1 )

c     ----------------

# ifdef MPI
      contains

      subroutine Mbroadcast
      use Mwrapper
      use, intrinsic :: iso_fortran_env
      implicit none
      call Mcastl(irrep_sub)
      call Mcastl(psub)
      call Mcastl(block)
      call Mcast(nsub)
      end subroutine Mbroadcast
# endif

c     ----------------

      end

c     ----------------

c     Returns irreducible representation in irrep(:dim,*) for symop isym of subspace b1:b2,ikpt,ispin
c     olappw        : plane-wave overlap <G|G'> at ikpt1=kptsym(ikpt,isym).
c     lcheck = true : Check unitarity.
c
c     Note that isym can rotate ikpt to a different kpoint ikpt1=kptsym(ikpt,isym).
c     Strictly speaking, irrep then corresponds to the subspace at b1:b2,ikpt1,ispin.
c
c     LOAD: If cmt is not allocated, allocate cmt/cpw and read the required wave functions from file. (Deallocate cmt/cpw at the end.)
c           If cmt is allocated, the arrays are used.
c           Structure: cmt( :maxlmindx , :ncent , ..bands.. , ..kpoints.. , :nspin3 )
c                      cpw( :maxgpt    ,          ..bands.. , ..kpoints.. , :nspin3 )
c           ..bands..   : Must comprise b1:b2
c           ..kpoints.. : Must comprise ikpt and ikpt1[=kptsym(ikpt,ispin)], kindx defined suitably
      subroutine get_irrep(irrep,dim,b1,b2,ikpt,ispin,isym,olappw,lcheck)
      use global
      use wrapper, only: identity
# ifdef LOAD
      use readwrite, only: read_wavef0
# endif
      implicit none
      integer,     intent(in)  :: b1,b2,ikpt,ispin,isym,dim
      logical,     intent(in)  :: lcheck
      MCOMPLEX_dp, intent(in)  :: olappw(*)
      complex_dp,  intent(out) :: irrep(dim,*)      
      complex_dp,  allocatable :: cmthlp(:,:,:,:),cmthlp1(:,:,:,:)
      MCOMPLEX_dp, allocatable :: cpwhlp(:,:,:),cpwhlp1(:,:,:)
      complex_dp               :: cphase,cdum
      complex_dp               :: wfolap
      real_dp                  :: err
      integer                  :: nb,ikpt1
      integer                  :: ib,j,k
# ifdef LOAD
      logical                  :: ldealloc
# else
      if(deg(b1,ikpt,ispin)/=b2) Bug('Not a degenerate subspace: b1:b2.')
# endif
      ikpt1 = kptsym(ikpt,isym)
      nb    = b2 - b1 + 1
      if(b2<b1)  Bug('Input band index is not first of degenerate subspace.')
      if(nb>dim) Bug('Leading dimension too small.')
# ifdef LOAD
      ldealloc = .false.
      if(.not.associated(cmt)) then
        ldealloc = .true.
        k        = 1 ; if(ikpt/=ikpt1) k = 2
        allocate ( cmt(maxlmindx,ncent,b1:b2,k,nspin3) )
        allocate ( cpw(maxgpt,         b1:b2,k,nspin3) )
        if(k==1) then ; call read_wavef0([(ib,ib=b1,b2)],[ikpt],      ispin,cmt,cpw)          
        else          ; call read_wavef0([(ib,ib=b1,b2)],[ikpt,ikpt1],ispin,cmt,cpw)
        endif
        kindx(ikpt)  = 1
        kindx(ikpt1) = k
      endif
# endif
      allocate ( cmthlp(maxlmindx,ncent,nb,nspin3),cpwhlp(maxgpt,nb,nspin3) )
      allocate ( cmthlp1(maxlmindx,ncent,nb,nspin3),cpwhlp1(maxgpt,nb,nspin3) )
      Inv ( cphase = exp(-img * 2*pi * dot_product(kpt(:,ikpt1),sym(isym)%transl) ) )
      do j = 1,nb
        ib = b1 + j - 1
# ifndef old_trafo
        call wavefunction(cmthlp (:,:,j,:),cpwhlp (:,j,:),     ib,ikpt1, ifLoad(1,ispin) )
        call waveftrafo  (cmthlp1(:,:,j,:),cpwhlp1(:,j,:),isym,ib,ikpt,  ifLoad(1,ispin) )
# else
        call wavefunction(cmthlp (:,:,j,1),cpwhlp (:,j,1),          ib,ikpt1, ifLoad(1,ispin) )
        call waveftrafo  (cmthlp1(:,:,j,1),cpwhlp1(:,j,1),ikpt,isym,ib,       ifLoad(1,ispin) )
        if(l_soc) then
          call wavefunction(cmthlp (:,:,j,2),cpwhlp (:,j,2),           ib ,ikpt1, 2)
          call waveftrafo  (cmthlp1(:,:,j,2),cpwhlp1(:,j,2),ikpt,isym, ib ,       2)
          call waveftrafo_soc(cmthlp1(:,:,j,:),maxlmindx*ncent,isym)
          call waveftrafo_soc(cpwhlp1(:,j,:),maxgpt,isym)
        endif
# endif
      enddo
      do j = 1,nb
        do k = 1,nb
          irrep(k,j) = wfolap( cmthlp(:,:,k,1),cpwhlp(:,k,1),
     &                         cmthlp1(:,:,j,1),cpwhlp1(:,j,1),ikpt1,ikpt1,ispin,
     &                         cdum,olappw ) Inv( * cphase )
          if(l_soc) then
            irrep(k,j) = irrep(k,j) + wfolap( cmthlp(:,:,k,2),cpwhlp(:,k,2),
     &                                        cmthlp1(:,:,j,2),cpwhlp1(:,j,2),ikpt1,ikpt1,nspin,
     &                                        cdum,olappw ) Inv( * cphase )
          endif
        enddo
      enddo
      deallocate (cmthlp,cpwhlp,cmthlp1,cpwhlp1)
      Load( if(ldealloc) deallocate ( cmt,cpw ) )
      if(lcheck) then
        err = sum(abs( matmul ( irrep(:nb,:nb) , conjg(transpose(irrep(:nb,:nb))) ) - identity(nb) ))
        if( err > 1d-10*nb**2*ncent ) then
          write(0,'(A,I3,A,I2,A,F14.10)') 'Irrep not unitary. Symop:',isym,', dimension:',nb,', error:',err
          write(0,'(2F20.10)') irrep(:nb,:nb)
          Warn('Irrep is not unitary.')
        endif
      endif
      end

c     ----------------            

c
c Calculates the representations of the symmetry operations in the Wannier basis. -> irrep_wan
c (These are not necessarily irreducible. So, "irrep" is a bit misleading.)
c
c For a given operation P, they are defined by
c
c irrep_wan(n',n)_R = < w_n'0 | P | w_nR >
c
c and are identical to the Fourier transform of the irreps of the Wannier Bloch functions
c
c irrep_wan(n',n)_R = 1/nkpt SUM(k) exp(-ikR) < w_n'k' | P | w_nk > .
c
c They vanish for all but one R if the Wannier centers coincide with the atomic centers:
c
c irrep_wan(n',n) = < w_n'k=0 | P | w_nk=0 >.
c
c This is the case we support here. (Otherwise, the program stops with an error.)
c
c The following can be uncommented to use a safer (but very slow) algorithm.
c# define TEST_SAFE
      subroutine irreps_wannier

      use wrapper, only: identity
      use global
      use util, only: chr
      use, intrinsic :: iso_fortran_env
      implicit none

      MCOMPLEX_dp,   allocatable :: cpwhlp(:,:),olappw(:,:)
      complex_dp,    allocatable :: cmthlp(:,:,:)
      complex_dp,    allocatable :: irrep(:,:),irrep_bwan(:,:,:,:)
      complex_dp                 :: cphase,cdum
      real_dp                    :: transl(3),rvec(3),err,rdum
      integer                    :: tst(nkpt,nsym)
      integer                    :: ib1,ib2,ikpt,k0,ikpt0,ikpt1,ispin,i,j,k,isym,isym0,isymk,isymk1,isymi
      complex_dp                 :: wfolap
      real                       :: time

      if(l_soc) Error('irreps_wannier not implemented for SOC.')
# ifdef LOAD      
      Error('IRREP (WANNIER) not implemented for -DLOAD.')
# endif      

      if(nspin1==1) then ; write(6,'(/A'NoA) 'Calculate irreducible representations for Wannier basis ...'
      else               ; write(6,'(/A)')   'Calculate irreducible representations for Wannier basis'
      endif
      call cpu_time(time)

      allocate ( irrep_wan(nwan,nwan,nsym,nspin) )
      allocate ( irrep_bwan(nwan,nwan,nkpt,nsym) )
      if(lkptadd) then
        allocate ( irrep_wan1(nwan,nwan,nsym,nspin) )
        irrep_wan1 = 0
      endif

c
c     Calculate irreps
c
      do ispin = 1,nspin
        if(nspin1>1) then
          if(ispin==1) write(6,'(A'NoA) '  Spin up... '
          if(ispin==2) write(6,'(A'NoA) '  Spin down... '
        endif
        irrep_bwan = 0
        tst        = 0

# ifdef TEST_SAFE
        if(storeibz) Error('TEST_SAFE not implemented without STOREBZ.')
        if(lkptadd)  Error('TEST_SAFE not implemented for additional k point.')
        Error('New trafo routines not implemented for TEST_SAFE.')
        do ikpt = 1,nkpt
          allocate ( olappw(ngpt(ikpt),ngpt(ikpt)) )
          allocate ( cmthlp(maxlmindx,ncent,wanbandi:wanbandf),cpwhlp(maxgpt,wanbandi:wanbandf) )
          allocate ( irrep(wanbandi:wanbandf,wanbandi:wanbandf) )
          do isym = 1,nsym
            ikpt1  = kptsym(ikpt,isym)
            cphase = exp(-img * 2*pi * dot_product(kpt(:,ikpt1),sym(isym)%transl) )
            call olap_gpt(olappw,ngpt(ikpt1),ikpt1,ikpt1)
            do ib2 = wanbandi,wanbandf
              call waveftrafo(cmthlp(:,:,ib2),cpwhlp(:,ib2),ikpt,isym,ib2,ispin)
            enddo
            do ib1 = wanbandi,wanbandf
              do ib2 = wanbandi,wanbandf
                irrep(ib1,ib2) = wfolap( cmt(:,:,ib1,ikpt1,ispin),cpw(:,ib1,ikpt1,ispin),
     &                                   cmthlp(:,:,ib2),cpwhlp(:,ib2),ikpt1,ikpt1,ispin,cdum,olappw ) * cphase
              enddo
            enddo
            irrep_bwan(:,:,ikpt,isym) = matmul ( conjg(transpose(uwan(wanbandi:wanbandf,:,ikpt1,ispin))) ,
     &                                  matmul ( irrep ,         uwan(wanbandi:wanbandf,:,ikpt ,ispin) ) )
            tst(ikpt,isym) = 1
          enddo
          deallocate(olappw,cmthlp,cpwhlp,irrep)
        enddo
# else

        ! loop over IBZ (k0)
        do k0 = 1,nkpti+1
          ikpt0 = k0
          if(k0==nkpti+1) then
            if(lkptadd) then ; ikpt0 = nkpt + 1
            else             ; exit
            endif
          endif
          allocate ( olappw(ngpt(ikpt0),ngpt(ikpt0)) )
          call olap_gpt(olappw,ngpt(ikpt0),ikpt0,ikpt0)

          ! loop over degenerate subspaces of k0
          i = wanbandi
          do while(i<=min(wanbandf,nband(ikpt0,ispin)))
            call getdeg(ib1,ib2,i,ikpt0,ispin)
            if(ib2>wanbandf) exit
            if(ib1<wanbandi) then ; i = ib2 + 1 ; cycle ; endif
            allocate ( cmthlp(maxlmindx,ncent,ib1:ib2),cpwhlp(maxgpt,ib1:ib2) )
            allocate ( irrep(ib1:ib2,ib1:ib2) )

            do isym0 = 1,nsym
              if(kptsym(ikpt0,isym0)==ikpt0) then

                ! calculate irrep of current subspace < phi_k0n' | P0 | phi_k0n > -> irrep
                cphase = exp(-img * 2*pi * dot_product(kpt(:,ikpt0),sym(isym0)%transl) )
                do j = ib1,ib2
# ifdef old_trafo                  
                  call waveftrafo(cmthlp(:,:,j),cpwhlp(:,j),ikpt0,isym0,j,ispin)
# else
                  call waveftrafo(cmthlp(:,:,j),cpwhlp(:,j),isym0,j,ikpt0,ispin)
# endif
                enddo
                do j = ib1,ib2
                  do k = ib1,ib2
                    irrep(k,j) = wfolap( cmt(:,:,k,k0,ispin),cpw(:,k,k0,ispin),
     &                                   cmthlp(:,:,j),cpwhlp(:,j),ikpt0,ikpt0,ispin,cdum,olappw ) * cphase
                  enddo
                enddo
                if(l_soc) then
                  Error('SOC not implemented in this routine.')
                  do j = ib1,ib2
                    call waveftrafo(cmthlp(:,:,j),cpwhlp(:,j),ikpt0,isym0,j,2)
                  enddo
                  do j = ib1,ib2
                    do k = ib1,ib2
                      irrep(k,j) = irrep(k,j) +
     &                             wfolap( cmt(:,:,k,k0,2),cpw(:,k,k0,2),
     &                                     cmthlp(:,:,j),cpwhlp(:,j),ikpt0,ikpt0,nspin,cdum,olappw ) * cphase
                    enddo
                  enddo
                endif

                ! loop over k and k' to calculate < w_k'n' | P | w_kn > -> irrep_bwan
                do ikpt = 1,nkpt    ; if(kptp(ikpt )/=ikpt0) cycle
                  do ikpt1 = 1,nkpt ; if(kptp(ikpt1)/=ikpt0) cycle

                    isymk  = symkpt(ikpt)
                    isymk1 = symkpt(ikpt1)
                    isymi  = sym(isymk)%inv
                    isym   = symtab(isymk1,isym0)
                    isym   = symtab(isym  ,isymi) ! the current isym: P = Pk' * P0 * Pk^(-1)
                    isymi  = sym(isymk1)%inv
                    transl = matmul( sym(isym) %rot , sym(isymk)%transl ) + sym(isym)  %transl              ! The translation in irrep is the one of P0, but
                    transl = matmul( sym(isymi)%rot , transl              - sym(isymk1)%transl )            ! Pk'^(-1)*P*Pk requires 'transl' defined here,
                    cphase = exp(-img * 2*pi * dot_product ( kpt(:,ikpt0) , -sym(isym0)%transl + transl ) ) ! leading to the correcting phase factor 'cphase'.
                    cphase = cphase / (conjg(phase(ikpt1))*phase(ikpt))                                     ! INV: also consider phases
                    if(sum(abs(transl-sym(isym0)%transl-nint(transl-sym(isym0)%transl)))>1d-12)
     &                                           Bug('transl-transl0 is not a lattice vector.')
                    if(ikpt1/=kptsym(ikpt,isym)) Bug('k-point symmetry error.')

                    if(isym>nsymt) then
                      call teststop('Time-reversal symmetry in irreps_wannier')
                      irrep_bwan(:,:,ikpt,isym) = irrep_bwan(:,:,ikpt,isym) +
     &                                            matmul ( conjg(transpose(uwan(ib1:ib2,:,ikpt1,ispin))) ,
     &                                            matmul ( irrep ,   conjg(uwan(ib1:ib2,:,ikpt ,ispin))) ) * cphase
                    else
                      irrep_bwan(:,:,ikpt,isym) = irrep_bwan(:,:,ikpt,isym) +
     &                                            matmul ( conjg(transpose(uwan(ib1:ib2,:,ikpt1,ispin))) ,
     &                                            matmul ( irrep ,         uwan(ib1:ib2,:,ikpt ,ispin) ) ) * cphase
     &
                    endif
                    tst(ikpt,isym) = 1
                  enddo
                enddo

                ! Case additional k point -> irrep_wan1
                if(ikpt0==nkpt+1) then
                  if(isym0>nsymt) then ; call teststop('Time-reversal symmetry in irreps_wannier')
                    irrep_wan1(:,:,isym0,ispin) = irrep_wan1(:,:,isym0,ispin) +
     &                                            matmul ( conjg(transpose(uwan(ib1:ib2,:,ikpt0,ispin))) ,
     &                                            matmul ( irrep ,   conjg(uwan(ib1:ib2,:,ikpt0,ispin))) )
                  else
                    irrep_wan1(:,:,isym0,ispin) = irrep_wan1(:,:,isym0,ispin) +
     &                                            matmul ( conjg(transpose(uwan(ib1:ib2,:,ikpt0,ispin))) ,
     &                                            matmul ( irrep ,         uwan(ib1:ib2,:,ikpt0,ispin) ) )
                  endif
                endif

              endif

            enddo
            deallocate ( cmthlp,cpwhlp,irrep )
            i = ib2 + 1
          enddo
          deallocate(olappw)
        enddo
# endif

        ! Check irrep_bwan and define irrep_wan
        call cpu_done(time)
        if(any(tst/=1)) Bug('Irreps incomplete.')
        err = 0
        do isym = 1,nsym
          do j = 1,nwan
            ! Determine R from which isym rotates Wannier center j into unit cell at origin
            rvec = matmul(sym(isym)%rot,wancent(:,j,ispin)) + sym(isym)%transl
            i    = minloc( [( sum(abs(rvec-wancent(:,i,ispin)-nint(rvec-wancent(:,i,ispin)))) ,i=1,nwan)] ,1)
            rvec = -matmul(sym(sym(isym)%inv)%rot,rvec - wancent(:,i,ispin))
            if(sum(abs(matmul(sym(isym)%rot,rvec+wancent(:,j,ispin))+sym(isym)%transl-wancent(:,i,ispin)))>1d-10)
     &        Error('No unique R for P(R+w)=w''. This case is not supported.')
            ! Check if irrep_bwan shows correct exp(ikR) periodicity
            do i = 1,nwan
              do ikpt = 2,nkpt
                rdum = abs( irrep_bwan(i,j,ikpt,isym) - irrep_bwan(i,j,1,isym) * exp(img*2*pi*dot_product(kpt(:,ikpt),rvec)) )
                err  = err + rdum
                if(rdum>1d-6) Error('Large deviation of irrep_bwan from exp(ikR) periodicity: '//Chf(rdum,'F13.10'))
              enddo
              ! Define irrep_wan
              irrep_wan(i,j,isym,ispin) = irrep_bwan(i,j,1,isym)
            enddo
          enddo
          ! Check if irrep_wan is unitary
          if(sum(abs( matmul(irrep_wan(:,:,isym,ispin),conjg(transpose(irrep_wan(:,:,isym,ispin))))
     &               - identity(nwan)))>1d-10) Error('Broken unitarity.')
        enddo
        write(6,'(A,F15.12)') '  Deviation from exp(ikR) periodicity:',err/(nwan**2*nsym*(nkpt-1))

      enddo ! spin loop

      deallocate(irrep_bwan)

      if(nspin==2) then
        if(sum(abs(irrep_wan(:,:,:,1)-irrep_wan(:,:,:,2)))<1d-8) then
          write(6,'(A)') 'Spin-up and -down transformation matrices are identical.'
        endif
      endif

      end























c backup
# if 0

      subroutine irreps_wannier0

      use wrapper, only: identity
      use util,    only: chr
      use global
      use, intrinsic :: iso_fortran_env
      implicit none

      real_dp                    :: olap_mt(maxindx,maxindx,0:maxlcut,ntype)
      MCOMPLEX_dp,   allocatable :: cpwhlp(:,:),olappw(:,:)
      complex_dp,    allocatable :: cmthlp(:,:,:)
      complex_dp,    allocatable :: irrep(:,:)
      complex_dp                 :: olapmt(maxlmindx,maxlmindx,ncent,nspin),cexp,cphase
      integer,       allocatable :: band(:)
      integer                    :: tst(nkpt,nsym,nspin)
      integer                    :: ib1,ib2,ikpt,ikpt0,ikpt1,ispin,i,j,k,isym,isym0
      complex_dp                 :: wfolap
      real                       :: time

      write(6,'(/A'NoA) 'Calculate irreducible representations for Wannier bands ...'
      call cpu_time(time)

      allocate ( irrep_wan(nwan,nwan,nsym,nspin) )
      irrep_wan = 0
      tst       = 0

c     Calculate irreps
      call wfolap_init_mt(olapmt,[0d0,0d0,0d0])
 1    do ispin = 1,nspin
        do ikpt0 = 1,nkpti
          allocate ( olappw(ngpt(ikpt0),ngpt(ikpt0)) )
          call olap_gpt(olappw,ngpt(ikpt0),ikpt0,ikpt0)
          i = wanbandi
          do while(i<=wanbandf)
            call getdeg(j,k,i,ikpt0,ispin) !!!!!!!!test
            ib1 = i
            ib2 = deg(ib1,ikpt0,ispin)
            if(ib2<ib1) then
              ib1 = ib2
              ib2 = deg(ib1,ikpt0,ispin)
            endif
            if(j/=ib1.or.k/=ib2) Error('getdeg failed.') !!!!!!!!!test
            if(ib2>wanbandf) exit
            if(ib1<wanbandi) then ; i = ib2 + 1 ; cycle ; endif
            allocate ( cmthlp(maxlmindx,ncent,ib1:ib2),cpwhlp(maxgpt,ib1:ib2) )
            allocate ( irrep(ib1:ib2,ib1:ib2) )
            do isym0 = 1,nsym

              if(kptsym(ikpt0,isym0)==ikpt0) then
                cphase = exp(-img * 2*pi * dot_product(kpt(:,ikpt0),sym(isym0)%transl) )
                do j = ib1,ib2
                  call waveftrafo(cmthlp(:,:,j),cpwhlp(:,j),ikpt0,isym0,j,ispin)
                enddo
                do j = ib1,ib2
                  do k = ib1,ib2
                    irrep(k,j) = wfolap( cmt(:,:,k,ikpt0,ispin),cpw(:,k,ikpt0,ispin),
     &                                   cmthlp(:,:,j),cpwhlp(:,j),ikpt0,ikpt0,ispin,olapmt,olappw,.true. ) * cphase
                  enddo
                enddo
                if(l_soc) then
                  do j = ib1,ib2
                    call waveftrafo(cmthlp(:,:,j),cpwhlp(:,j),ikpt0,isym0,j,2)
                  enddo
                  do j = ib1,ib2
                    do k = ib1,ib2
                      irrep(k,j) = irrep(k,j) +
     &                             wfolap( cmt(:,:,k,ikpt0,ispin),cpw(:,k,ikpt0,ispin),
     &                                     cmthlp(:,:,j),cpwhlp(:,j),ikpt0,ikpt0,nspin,olapmt,olappw,.true. ) * cphase
                    enddo
                  enddo
                endif
c                write(*,*) sum(abs(matmul(irrep,conjg(transpose(irrep)))-identity(ib2-ib1+1))); read(*,*)
                do ikpt = 1,nkpt    ; if(kptp(ikpt )/=ikpt0) cycle
                  do ikpt1 = 1,nkpt ; if(kptp(ikpt1)/=ikpt0) cycle
                    isym = symtab(symkpt(ikpt1),isym0)
                    isym = symtab(isym,sym(symkpt(ikpt))%inv)
                    cexp = exp(img * 2*pi * dot_product(kpt(:,ikpt1)+gkptsym(:,ikpt1,isym),sym(isym)%transl))
                    if(abs(1-cexp)>1d-10) Error('Non-symmorphic operations not supported.')
                    if(ikpt1/=kptsym(ikpt,isym)) Error('error.')
                    if(ib1==wanbandi) tst(ikpt1,isym,ispin) = tst(ikpt1,isym,ispin) + 1
                    if(isym>nsymt) then
                      call teststop('Time-reversal symmetry in irreps_wannier')
                      irrep_wan(:,:,isym,ispin) = irrep_wan(:,:,isym,ispin) + cexp / nkpt *
     &                                            matmul ( conjg(transpose(uwan(ib1:ib2,:,ikpt1,ispin))) ,
     &                                            matmul ( irrep ,   conjg(uwan(ib1:ib2,:,ikpt ,ispin))) )
                    else
                      irrep_wan(:,:,isym,ispin) = irrep_wan(:,:,isym,ispin) + cexp / nkpt *
     &                                            matmul ( conjg(transpose(uwan(ib1:ib2,:,ikpt1,ispin))) ,
     &                                            matmul ( irrep ,         uwan(ib1:ib2,:,ikpt ,ispin) ) )
                    endif
                  enddo
                enddo
              endif

            enddo
            deallocate ( cmthlp,cpwhlp,irrep )
            i = ib2 + 1
          enddo
          deallocate(olappw)
        enddo
      enddo

      call cpu_done(time)

      ! check if all terms are considered (once)
      if(any(tst/=1)) Error('error.')

      ! check if the irreps are unitary
      allocate ( irrep(nwan,nwan) )
      do ispin = 1,nspin
        do isym = 1,nsym
c          if(ispin==1) then
c            write(6,'(A,I3)') 'Symmetry',isym
c            write(6,'(5F12.8)') real(transpose(irrep_wan(:,:,isym,ispin)))
c          endif
          irrep = matmul(irrep_wan(:,:,isym,ispin),conjg(transpose(irrep_wan(:,:,isym,ispin))))
          do i = 1,nwan
            irrep(i,i) = irrep(i,i) - 1
          enddo
          if(sum(abs(irrep)**2)>1d-20) Error('Wannier symmetry broken.')
        enddo
      enddo

      do ispin = 1,nspin
        do isym = 1,nsym
          do isym0 = 1,nsym
            if(sum(abs(matmul(irrep_wan(:,:,isym,ispin),irrep_wan(:,:,isym0,ispin))-
     &        irrep_wan(:,:,symtab(isym,isym0),ispin)))>1d-12) Error('Wannier symmetry broken.')
          enddo
        enddo
      enddo
      deallocate ( irrep )

      end

c     ----------------

      subroutine irreps_wannier

      use global
      use util, only: chr
      use, intrinsic :: iso_fortran_env
      implicit none

      real_dp                    :: olap_mt(maxindx,maxindx,0:maxlcut,ntype)
      MCOMPLEX_dp                :: olappw(ngpt,ngpt)
      MCOMPLEX_dp,   allocatable :: cpwhlp(:,:)
      complex_dp,    allocatable :: cmthlp(:,:,:)
      complex_dp                 :: olapmt(maxlmindx,maxlmindx,ncent,nspin),cfac
      real_dp                    :: kvec(3)
      integer,       allocatable :: band(:)
      integer                    :: nband1,nband2,nkpt1
      integer                    :: ib1,ib2,nb,ikpt,ikpt1,ispin,iloop,i,j,k,isym,nelem,ielem,maxelem
      complex_dp                 :: wfolap
      real                       :: time

# ifdef INV
      Error('Not tested with inversion symmetry yet!')
# endif

      write(6,'(/A'NoA) 'Calculate irreducible representations for Wannier bands ...'
      call cpu_time(time)

c     Calculate irreps
      call wfolap_init(olappw,olapmt,[0d0,0d0,0d0])
      maxelem = 0
 1    do ispin = 1,nspin
        do ikpt = 1,nkpti
          nelem = 0
          ielem = 0
          i     = wanbandi
          do while(i<=wanbandf)
            ib1 = i
            ib2 = deg(ib1,ikpt,ispin)
            if(ib2<ib1) then
              ib1 = ib2
              ib2 = deg(ib1,ikpt,ispin)
            endif
            if(ib2>wanbandf) exit
            nb    = ib2 - ib1 + 1
            nelem = nelem + nb**2
            if(iloop==2) then
              allocate ( cmthlp(maxlmindx,ncent,nb),cpwhlp(ngpt,nb) )
              do isym = 1,nsym
                kvec  = matmul(sym(isym)%rrot,kpt(:,ikpt))!-kpt(:,ikpt)
                kvec  = modulo1r(kvec) * nkpt3
                ikpt1 = pkpt(nint(kvec(1)),nint(kvec(2)),nint(kvec(3)))

                if(all(kvec-nint(kvec)<1d-10)) then
                  do j = ib1,ib2
                    call waveftrafo(cmthlp(:,:,j-ib1+1),cpwhlp(:,j-ib1+1),cfac,ikpt,isym,j,ispin)
                  enddo
                  do j = 1,nb
                    do k = 1,nb
                      ielem = ielem + 1
                      irrep_wan(,isym,ikpt,ispin) = wfolap( cmt(:,:,k+ib1-1,ikpt,ispin),cpw(:,k+ib1-1,ikpt,ispin),
     &                                                           cmthlp(:,:,j),cpwhlp(:,j),ispin,
     &                                                           olappw,olapmt,.true. ) * cfac * conjg(phase(ikpt1))
                    enddo
                  enddo
                endif

              enddo
              deallocate ( cmthlp,cpwhlp )
            endif
            i = i + nb
          enddo
          maxelem = max(maxelem,nelem)
        enddo
      enddo
      if(iloop==1) then
        iloop = 2
        allocate ( irrep_wan(maxelem,nsym,nkpti,nspin) )
        goto 1
      endif

      call cpu_done(time)

      end

      subroutine try_irreps
      use global
      use, intrinsic :: iso_fortran_env
      implicit none
      real_dp :: kvec(3)
      integer :: sym2(3),i,isym
      kvec = 0
      do isym = 2,nsym
        sym2(1) = isym
        i       = 1
        do while(symtab(sym2(i),isym)/=1)
          sym2(i+1) = symtab(sym2(i),isym)
          i         = i + 1 ; if(i>2) exit
        enddo
        if(i<=2) write(*,*) i,sym2(:i)
      enddo
      stop 'done'
      end

# endif


