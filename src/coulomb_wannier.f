c--------------------------------------------------------------------------------
c Copyright (c) 2021 Peter Grünberg Institut, Forschungszentrum Jülich, Germany
c This file is part of SPEX and available as free software under the conditions
c of the MIT license as expressed in the LICENSE file in more detail.
c--------------------------------------------------------------------------------

c Updates Coulomb matrix in Wannier basis for current k point ikpt0.
c
c     barew(12,34) = <12|v|34> = INT 1*(r) 3(r) v(r,r') 2*(r') 4(r') dV dV'
c   screenw(12,34) = <12|W|34> = INT 1*(r) 3(r) W(r,r') 2*(r') 4(r') dV dV'
c
c In the routine, the indexing is temporarily transformed to "summation order" (13,42).
c
c Enforce old algorithm:
c # define oldalgo
c
# ifdef CMP_0405adj
#   define oldalgo
# endif
# ifdef oldalgo
#   warning Old algorithm selected
# endif

# include "cppmacro.h"
      subroutine coulomb_wannier(screenw,barew,screen,head,ctrafo,n,ikpt0,nfrq,rsite,nsite,spin)
      use global
      use wrapper
      use util
      Mpi( use Mwrapper )
      use, intrinsic :: iso_fortran_env
      use, intrinsic :: iso_c_binding
      implicit none
      integer,     intent(in)    :: n,ikpt0,nfrq,spin,nsite,rsite(3,nsite)
      complex_dp,  intent(inout) :: screenw(nwan*nwan,nwan*nwan,nfrq,nsite)
      complex_dp,  intent(inout) :: barew(nwan*nwan,nwan*nwan,nsite)
      complex_dp,  intent(in)    :: screen(n,n,nfrq),head(3,3,nfrq)
      MCOMPLEX_dp, intent(in)    :: ctrafo(nbasm(ikpt0),n)
      MCOMPLEX_dp                :: olap(ngptm(ikpt0)*(ngptm(ikpt0)+1)/2)
      MCOMPLEX_dp, allocatable   :: cprod(:,:,:,:,:)
      MCOMPLEX_dp, pointer_cnt   :: olappw(:,:)
      complex_dp,  allocatable   :: cprodw(:,:,:),cprod0(:,:,:),cprod1(:,:,:),weight(:,:)
      complex_dp                 :: carr(nbasm(ikpt0)),chlp(nwan*nwan,nwan*nwan)
      complex_dp                 :: carr2(nwanband,nwanband,nbasm(ikpt0)),carr3(nwan**2,nbasm(ikpt0))
      logical                    :: needed(nkpt,nsym)
      integer                    :: b(4,3)
      integer                    :: kpt1(nkpt),nkpts(nkpt),sym1(nsym),nkpt1,nsym1,kpt1p(nkpt),symkpt1(nkpt)
      integer                    :: nwan2
      integer                    :: ispin1,ispin2,ispin,minspin,maxspin
      integer                    :: ikpt0_,ikpt1,ikpt2,ikpt,ifreq,iwan,isite,isym
      integer                    :: i,j,k,m,nk,k1,k2
      integer                    :: kptsum
      real_dp,     external      :: scaling_linear
      complex_dp,  allocatable, save :: irrep_bwan(:,:)
      integer,     allocatable, save :: pirrep_bwan(:,:,:)
# ifdef oldalgo
      complex_dp                 :: cexp
# endif
# ifdef MPI
      integer                    :: nkpt0
      logical                    :: lpar
      type(c_ptr)                :: ptr
      integer                    :: win_olappw,Merr
      MCOMPLEX_dp, pointer_cnt   :: Ninit_olappw(:)
# endif
      real time1,time2

      if(l_soc) Error('SOC not implemented for this routine.')

      if     (spin==0) then ; ispin1 = 1 ; ispin2 = nspin ! average over both spins
      else if(spin==1) then ; ispin1 = 1 ; ispin2 = 1
      else if(spin==2) then ; ispin1 = 2 ; ispin2 = 2
      else if(spin==3) then ; ispin1 = 1 ; ispin2 = 2
      else if(spin==4) then ; ispin1 = 2 ; ispin2 = 1
      else                  ; Error('unknown spin index.')
      endif
      minspin = min(ispin1,ispin2)
      maxspin = max(ispin1,ispin2)

      nwan2 = nwan**2

      ! Transform to summation order
      do isite = 1,nsite
        call transpose2(barew(:,:,isite),1)
        do ifreq = 1,nfrq
          call transpose2(screenw(:,:,ifreq,isite),1)
        enddo
      enddo

      ! Weights and Gamma point
      allocate ( weight(nsite,count(kptp(:nkpt)==ikpt0)) )
      call cpu_time(time1)
      call weight_and_gamma
      Timing( Rwrite(*,'(A,$)') 'weight_and_gamma';Rcall cpu_done(time1) )

      ! Inverse IPW overlap <M~|M~>
      if(fullpw) then
        call olap_pwp(olap,gptm(:,pgptm(:ngptm(ikpt0),ikpt0)),ngptm(ikpt0))
        call Minverse(olap)
        Mpi( call Msum(olap) )
        Timing( Rwrite(*,'(A,$)') 'olap';Rcall cpu_done(time1) )
      endif

# ifdef LOAD
      if(wbloch) then ; call prepare_wbloch
      else            ; Error('New algorithm not implemented for LOAD. Please use old algorithm: WBLOCH.')
      endif
# endif

      ! Precalculate irrep_bwan (done here for all kpoints)
      if(.not.wbloch.and.use_sym) then
        if(.not.allocated(irrep_bwan) ) then
          call cpu_time(time1)
          do
            k = 1 ; Mpi( k1 = 0 )
            do ispin = 1,nspin1
              do ikpt = 1,nkpt
                if(allocated(irrep_bwan)) then
                  Nallocate0 ( olappw,(S_ ngpt(ikpt),ngpt(ikpt) S_) ) ! allocating olappw in shared memory improves the scalability!
                  Ocall olap_gpt(olappw,ngpt(ikpt),ikpt,ikpt)                  
                  Nfence(olappw)
                endif
                i = wanbandi
                do while(i<=wanbandf)
                  j = deg(i,ikpt,ispin)
                  if(j<wanbandi) then ; i = i + 1 ; cycle ; endif
                  if(j>wanbandf) exit
                  m = j - i + 1
                  if(allocated(irrep_bwan) lMOD(k1) ) then
                    pirrep_bwan(i,ikpt,ispin) = k ! note: ikpt is kpoint to which isym rotates!
                    do isym = 2,nsym 
                      ikpt0_ = kptsym(ikpt,sym(isym)%inv) ; if(.not.needed(ikpt0_,isym)) cycle
                      call get_irrep(irrep_bwan(k,isym),m,i,j,ikpt0_,ispin,isym,olappw,.true.) ! define -> irrep_bwan
                    enddo
                  endif
                  k = k + m**2 ; Mpi( k1 = k1 + 1 )
                  i = j + 1
                enddo
                if(allocated(irrep_bwan)) tNdeallocate ( olappw )
              enddo
            enddo
            if(allocated(irrep_bwan)) exit
            allocate ( irrep_bwan(k,2:nsym),pirrep_bwan(wanbandi:wanbandf,nkpt,nspin1) )
            ! check which ikpt/isym combinations are needed -> needed(:,:)
            irrep_bwan  = 0
            pirrep_bwan = 0
            needed      = .false.
            do k = 1,nkpt
              ikpt0_ = kptp(k)
              call getkpt1(kpt1,nkpt1,nkpts,sym1,nsym1,ikpt0_,0,.false.)
              call getkpt1p(kpt1p,kpt1,nkpt1,sym1,nsym1,symkpt1,.true.)
              do ikpt = 1,nkpt
                isym  = symtab( symkpt(k) , symkpt1(ikpt) )
                ikpt1 = kpt1(kpt1p(ikpt))    ; needed(ikpt1,isym) = .true.
                ikpt2 = kptsum(ikpt0_,ikpt1) ; needed(ikpt2,isym) = .true.
              enddo
            enddo
          enddo
          call cpu_time(time2) ; time2 = time2 - time1
          Mpi( call Msum(irrep_bwan) ; call Msum(pirrep_bwan) )
          Timing( Rwrite(*,'(A,$)') 'irrep_bwan';Rcall cpu_done(time1) ; Rwrite(*,*) ([size(irrep_bwan),size(pirrep_bwan)]) )
        endif
        ! Define EIBZ
        call getkpt1(kpt1,nkpt1,nkpts,sym1,nsym1,ikpt0,0,.false.)
        call getkpt1p(kpt1p,kpt1,nkpt1,sym1,nsym1,symkpt1,.true.)
      endif

# ifdef MPI      
      nkpt0 = count(kptp(:nkpt)==ikpt0)
# endif

      Allocate_ ( cprod0,(nwan2,nbasm(ikpt0),minspin:maxspin) ) ! mixed basis
      Allocate_ ( cprod1,(nwan2,n,minspin:maxspin) )            ! eigenfunctions of coulomb

      if(.not.wbloch.and.use_sym) then
# ifdef MPI
        call Mblocks(b,[nwanband,nkpt1,maxspin-minspin+1],3,scaling_linear)
        b(:2,3) = b(:2,3) - 1 + minspin
#   define BANDS  b(1,1):b(2,1)
#   define BANDS1 b(1,1)+wanbandi-1,b(2,1)+wanbandi-1
#   define KPTS   b(1,2):b(2,2)
#   define NKPTS  b(2,2)-b(1,2)+1
# else
#   define BANDS  :
#   define BANDS1 wanbandi,wanbandf
#   define KPTS   :
#   define NKPTS  nkpt1
# endif
        Allocate_ ( cprod,(nwanband,nwanband,nkpt1,nbasm(ikpt0),minspin:maxspin) )
        Mpi( cprod = 0 )
        do ispin = minspin,maxspin ; Mpi( if(ispin<b(1,3).or.ispin>b(2,3)) cycle )
          call cpu_time(time1)
          beginSplit(Mrank)
          call wavefproducts2_mt(cprod(:,BANDS,KPTS,:nbasp  ,ispin),kpt1(KPTS),NKPTS,ikpt0,ispin,ispin,BANDS1,wanbandi,wanbandf)          
          call wavefproducts2_pw(cprod(:,BANDS,KPTS,nbasp+1:,ispin),kpt1(KPTS),NKPTS,ikpt0,ispin,ispin,BANDS1,wanbandi,wanbandf)
          endSplit
          Timing( write(*,'(A,$)') 'wavefprod' Mpi(//Chr(Mrank)) ; call cpu_done(time1) )
        enddo
        Mpi( call Msum(cprod) )
      endif

      k1 = 1
      k2 = nkpt
      nk = nkpt
      
      k = 0
      do ikpt0_ = 1,nkpt
        if(kptp(ikpt0_)/=ikpt0) cycle ! loop over the symmetry-equivalent k points
        
        k = k + 1
# ifdef MPI
        if(k>nkpt0/Msize*Msize) then ; lpar = .true.  ; MrangeDef1(k1,k2,nkpt) ; nk = k2 - k1 + 1 ! calculate current k point in parallel
        else                         ; lpar = .false. ; Mcycle(k)                                 ! calculate current k point by one rank only
        endif
# endif
        if(wbloch)            then ; Allocate_ ( cprodw,(nwan2,k1:k2,nbasm(ikpt0)) )
        else if(.not.use_sym) then ; Allocate_ ( cprod,(nwanband,nwanband,k1:k2,nbasm(ikpt0),1) )
        endif 

        cprod0 = 0
        cprod1 = 0

        do ispin = minspin,maxspin

          call cpu_time(time1)
          if(wbloch) then
            call wavefproducts3_mt(cprodw,              [(ikpt,ikpt=k1,k2)],nk,ikpt0_,ispin,ispin,1,nwan,1,nwan)
            call wavefproducts3_pw(cprodw(:,:,nbasp+1:),[(ikpt,ikpt=k1,k2)],nk,ikpt0_,ispin,ispin,1,nwan,1,nwan)
            Inv( call symmetrize(cprodw,nwan2*nk,nbasm(ikpt0),-1) )
          else if(.not.use_sym) then
            beginSplit(Mrank)
            call wavefproducts2_mt(cprod,                  [(ikpt,ikpt=k1,k2)],nk,ikpt0_,ispin,ispin,wanbandi,wanbandf,
     &                                                                                               wanbandi,wanbandf)
            call wavefproducts2_pw(cprod(:,:,:,nbasp+1:,1),[(ikpt,ikpt=k1,k2)],nk,ikpt0_,ispin,ispin,wanbandi,wanbandf,
     &                                                                                               wanbandi,wanbandf)
            endSplit
          endif
          Timing( Rwrite(*,'(A,$)') 'wavefprod' ; Rcall cpu_done(time1) )

          ! Wannier BZ sum
          do ikpt = k1,k2
            if(wbloch) then
              carr3 = cprodw(:,ikpt,:)
            else if(use_sym) then
              carr2 = cprod(:,:,kpt1p(ikpt),:,ispin)
              isym  = symtab( symkpt(ikpt0_) , symkpt1(ikpt) )
              call mtrafo1(carr2,nwanband**2,nbasm(ikpt0),ikpt0,isym,-1)
              call uwan_cprod(carr3,carr2,ikpt0,kpt1(kpt1p(ikpt)),isym,ispin)
            else
              carr2 = cprod(:,:,ikpt,:,1)
              call uwan_cprod(carr3,carr2,ikpt0_,ikpt,1,ispin)              
            endif
            cprod0(:,:,ispin) = cprod0(:,:,ispin) + carr3 / nkpt ! gives a factor 1/nkpt**2 in the formula
          enddo
          Timing( Rwrite(*,'(A,$)') 'Wan-BZsum' ; Rcall cpu_done(time1) )
          
          Mpi( if(lpar) call Msum(cprod0(:,:,ispin)) )

          ! transform vectors back to ikpt0 (as coulomb and screen are defined there)
          if(ikpt0_/=ikpt0) then
            call mtrafo1(cprod0(:,:,ispin),nwan2,nbasm(ikpt0),ikpt0,-symkpt(ikpt0_),-1)
          endif
          ! transform to eigenvector basis
# ifdef MPI
          if(lpar) then
            cprod1(:,Mcol1(n),ispin) = matmul ( cprod0(:,:,ispin),MCONJG(ctrafo(:,Mcol1(n))) )
          else
# endif
          cprod1(:,:,         ispin) = matmul ( cprod0(:,:,ispin),MCONJG(ctrafo(:,:)) ) ; Mpi( endif )

          do iwan = 1,nwan2
            Mpi( if(.not.lpar.or.lpar.and.mod(iwan,Msize)==Mrank) then )
            if(fullpw) cprod0(iwan,nbasp+1:,ispin) = matvec ( olap , cprod0(iwan,nbasp+1:,ispin) )
            Mpi( else if(lpar) then ; cprod0(iwan,:,ispin) = 0 ; endif )
          enddo

          Timing( Rwrite(*,'(A,$)') 'Trafos' ; Rcall cpu_done(time1) )

        enddo

        Mpi( if(lpar) call Msum(cprod0) )
        Mpi( if(lpar) call Msum(cprod1) )

        if(spin==0) then
          ispin1 = 1
          ispin2 = 1
        endif

 1      continue

        chlp = 0
        do iwan = 1,nwan2 ; Mpi( if(lpar) then ; Mcycle(iwan) ; endif )
          chlp(:,iwan) = matvec ( conjg(cprod0(:,:,ispin1)) ,
     &                   matvec ( coulomb0(:nbasm(ikpt0)*(nbasm(ikpt0)+1)/2),cprod0(iwan,:,ispin2) ) )
        enddo
        if(spin==0.and.nspin==2) chlp = chlp / 2

        Timing( Rwrite(*,'(A,$)') 'coulomb0' ; Rcall cpu_done(time1) )

        do isite = 1,nsite
# ifdef oldalgo
          cexp             = exp( img * 2*pi * dot_product(kpt(:,ikpt0_) , rsite(:,isite)) )
          barew(:,:,isite) = barew(:,:,isite) + chlp * cexp / nkpt
# else
          barew(:,:,isite) = barew(:,:,isite) + chlp * weight(isite,k)
# endif
        enddo

        do ifreq = 1,nfrq
          do iwan = 1,nwan2 ; Mpi( if(lpar) then ; Mcycle(iwan) ; endif )
            chlp(:,iwan) = matvec ( conjg(cprod1(:,:,ispin1)),
     &                     matvec ( screen(:,:,ifreq),cprod1(iwan,:,ispin2)) )
          enddo
          if(spin==0.and.nspin==2) chlp = chlp / 2
          do isite = 1,nsite
# ifdef oldalgo
            cexp                     = exp( img * 2*pi * dot_product(kpt(:,ikpt0_) , rsite(:,isite)) )
            screenw(:,:,ifreq,isite) = screenw(:,:,ifreq,isite) + chlp * cexp / nkpt
# else
            screenw(:,:,ifreq,isite) = screenw(:,:,ifreq,isite) + chlp * weight(isite,k)
# endif
          enddo
        enddo

        Timing( Rwrite(*,'(A,$)') 'screenw' ; Rcall cpu_done(time1) )

        if(nspin==2) then
          if(spin==0) then
            if(ispin1==1) then
              ispin1 = 2
              ispin2 = 2
              goto 1
            endif
          endif
        endif

        if(wbloch)            then ; Deallocate_(cprodw)
        else if(.not.use_sym) then ; Deallocate_(cprod)
        endif

      enddo ! end loop over symmetry-equivalent k points

      Mpi( call Msum(barew,0); call Msum(screenw,0) )
      MnoR( barew = 0 ; screenw = 0 )
      Load( if(wbloch) tDeallocate(cmtu) )
      Load( if(wbloch) tDeallocate(cpwu) )
      if(.not.wbloch.and.use_sym) tDeallocate(cprod)
      Deallocate_(cprod0)
      Deallocate_(cprod1)
      deallocate(weight)
      if(ikpt0==nkpti.and.allocated(irrep_bwan)) deallocate(irrep_bwan,pirrep_bwan)

      ! Transform to normal order
      do isite = 1,nsite
        call transpose2(barew(:,:,isite),0)
        do ifreq = 1,nfrq
          call transpose2(screenw(:,:,ifreq,isite),0)
        enddo
      enddo

      contains

      subroutine uwan_cprod(cprod_out,cprod_in,ikpt0,ikpt1,isym,ispin)
      use global
      use wrapper
      use, intrinsic :: iso_fortran_env
      implicit none
      integer,    intent(in)  :: ikpt0,ikpt1,ispin,isym
      complex_dp, intent(in)  :: cprod_in(nwanband,nwanband,nbasm(ikpt0))
      complex_dp, intent(out) :: cprod_out(nwan,nwan,nbasm(ikpt0))
      complex_dp              :: uwan1(wanbandi:wanbandf,nwan),uwan2(wanbandi:wanbandf,nwan)
      integer                 :: ikpt2,i
      integer                 :: kptsum      
      if(any([ikpt0,ikpt1]<1.or.[ikpt0,ikpt1]>nkpt)) Bug('Parameter ikpt out of range.')
      ikpt2 = kptsum(ikpt0,ikpt1)
      if(isym==1) then
        uwan1 = conjg(uwan(:,:,ikpt1,ispin))
        uwan2 = conjg(uwan(:,:,ikpt2,ispin))
      else
        call get_uwan(uwan1,ikpt1,ispin,isym) ; uwan1 = conjg(uwan1)
        call get_uwan(uwan2,ikpt2,ispin,isym) ; uwan2 = conjg(uwan2)
      endif
      do i = 1,nbasm(ikpt0)
        cprod_out(:,:,i) = macmat ( uwan2 , matmat ( cprod_in(:,:,i) , uwan1 ) )
      enddo
      end subroutine uwan_cprod

      subroutine get_uwan(uwan1,ikpt,ispin,isym)
      implicit none
      integer,    intent(in)  :: ikpt,ispin,isym
      complex_dp, intent(out) :: uwan1(wanbandi:wanbandf,nwan)
      complex_dp, allocatable :: irrep(:,:)
      integer                 :: ikpt1,i1,i2,nb,irr
      uwan1 = 0
      ikpt1 = kptsym(ikpt,isym)
      i1    = wanbandi
      do while(i1<=wanbandf)
        i2 = deg(i1,ikpt,ispin)
        if(i2<wanbandi) then ; i1 = i1 + 1 ; cycle ; endif
        if(i2>wanbandf) exit
        nb             = i2 - i1 + 1 ; allocate(irrep(nb,nb))
        irr            = pirrep_bwan(i1,ikpt1,ispin)                            ; if(irr==0) Bug('Irrep pointer is zero.')
        irrep          = reshape ( irrep_bwan(irr:irr+nb**2-1,isym) , [nb,nb] ) ; if(all(irrep==0)) Bug('Found zero irrep.')
        uwan1(i1:i2,:) = macmat ( irrep , uwan(i1:i2,:,ikpt1,ispin) )
        i1             = i2 + 1
        deallocate(irrep)
      enddo
      end subroutine get_uwan

c -----------------------------

      ! Gamma contribution
      subroutine weight_and_gamma
      use, intrinsic :: iso_fortran_env
      implicit none
      complex_dp              :: mat(nwan,nwan,nspin)
      complex_dp              :: cdum
      complex_dp              :: div(nfrq,nsite),div0(nsite)
      integer                 :: i,j,k,l,m,n,isite,ikpt,ikpt1,ispin
      complex_dp, allocatable :: wtetra(:,:)

# ifndef oldalgo
      if(ikpt0==1) then ; Allocate_(wtetra,(nsite,nkpt))
      else              ; Allocate_(wtetra,(1,nkpt))
      endif
      do isite = 1,nsite ; j = 1 ; if(ikpt0==1) j = isite
        Mcycle(isite)
        call tetrahedron6_init(wtetra(j,:),rsite(:,isite))
c        cdum = 0
c        do i = 1,nkpt
c          cdum = cdum + wtetra(j,i) / exp( img * 2*pi * dot_product(kpt(:,i),rsite(:,isite)) ) / nkpt
c        enddo
c        carr(1) = 0
c        do i = 1,nkpt
c          carr(1) = carr(1) + abs(wtetra(j,i)/ exp( img * 2*pi * dot_product(kpt(:,i),rsite(:,isite)) )-cdum)
c        enddo
c        write(*,'(I3,2F20.16)') isite,carr(1)
        k = 0
        do ikpt1 = 1,nkpt ; if(kptp(ikpt1)/=ikpt0) cycle
          k               = k + 1
          weight(isite,k) = wtetra(j,ikpt1)
        enddo
      enddo
      Mdistr(weight(isite,:),isite,1,nsite)

      if(ikpt0==1) then
        Mdistr(wtetra(isite,:),isite,1,nsite)
        call gamma_divergence2(div0,wtetra,(1d0,0d0)*identity(3))
        div = 0
        do i = 1,nfrq
          cdum = ( head(1,1,i) + head(2,2,i) + head(3,3,i) ) / 3 ! trace
          if(abs(real(cdum))<1d60) then
            if(all(abs( head(:,:,i) - transpose(conjg(head(:,:,i))) )<1d-10)) then
              call gamma_divergence2(div(i,:),wtetra,head(:,:,i)) ! head is Hermitian
            else
              call gamma_divergence3(div(i,:),wtetra,head(:,:,i)) ! head is complex
            endif
          endif
        enddo
      endif
      Deallocate_(wtetra)
# else
      if(ikpt0==1) then
        if(divergence==0) call gamma_divergence(.false.)
        div0 = divergence
        div  = 0
        do i = 1,nfrq
          cdum = ( head(1,1,i) + head(2,2,i) + head(3,3,i) ) / 3 ! trace
          if(abs(real(cdum))<1d60) div(i,:) = divergence / cdum
        enddo
      endif
# endif

      if(ikpt0/=1) return
      MnoR(return)

      mat = 0
      do ispin = 1,nspin
        do ikpt = 1,nkpt
          ! Prepare matrix for Gamma contribution -> mat
          mat(:,:,ispin) = mat(:,:,ispin) + matmul ( transpose(uwan(:,:,ikpt,ispin)) , conjg(uwan(:,:,ikpt,ispin)) ) / nkpt
        enddo
      enddo
      m = 0
      do j = 1,nwan
        do i = 1,nwan
          m = m + 1
          n = 0
          do l = 1,nwan
            do k = 1,nwan
              n                = n + 1
              cdum             = conjg ( mat(i,j,ispin1) ) * mat(k,l,ispin2)
              barew(m,n,:)     = barew(m,n,:)     + 4*pi * cdum * div0 / vol
              screenw(m,n,:,:) = screenw(m,n,:,:) + 4*pi * cdum * div  / vol
            enddo
          enddo
        enddo
      enddo
      end subroutine weight_and_gamma

c -----------------------------

c     Similar to gamma_divergence but (a) integrand contains an exponential function and (b) a dielectric tensor (H=head) is considered.
c
c     div_out(R) = vol/(8pi^3) INT exp(ikR)/(kHk) d^3k - 1/nkpt SUM(k/=0) exp(ikR)/(kHk)
c
c     wtetra contains tetrahedron weights from tetrahedron6_init and, thus, includes the factor exp(ikR).
c
c     gamma_divergence2: head must be Hermitian (frequency zero or purely imaginary),
c     gamma_divergence3: general case (frequency may be real), slow.
c
      subroutine gamma_divergence2(div_out,wtetra,head)
      use global
# ifdef MPI
      use Mwrapper, only: Msum
# endif
      use, intrinsic :: iso_fortran_env
      implicit none
      complex_dp, intent(out) :: div_out(nsite)
      complex_dp, intent(in)  :: wtetra(nsite,nkpt),head(3,3)
      integer                 :: ix,iy,iz,n,ikpt,isite
      logical                 :: found
      real_dp,    parameter   :: expo = 5d-3
      integer,    parameter   :: nrad = 50000
      real_dp                 :: rrad,drad,kv1(3),kv2(3),kv3(3),k(3),rvec(3),k2,k1,r,rfac,si
      complex_dp              :: div(nsite),cfac
      real_dp                 :: rnorm(nsite)
      complex_dp              :: cdeterminant
# ifdef MPI
      logical                 :: ldum
      integer                 :: mm,Merr
# endif
      if(any(abs( head - transpose(conjg(head)) )>=1d-10)) Bug('Input matrix (head) must be Hermitian.')
      div_out = 0
      rrad    = sqrt(-log(expo)/expo)
      drad    = (rvol/nkpt)**(1d0/3)
      kv1     = rlat(:,1)/nkpt3(1)
      kv2     = rlat(:,2)/nkpt3(2)
      kv3     = rlat(:,3)/nkpt3(3)
      n       = 1 Mpi(+Mrank)
      found   = .true.
      do while(found)
        found = .false.
        div   = 0
        do ix = -n,n
          do iy = -(n-abs(ix)),n-abs(ix)
            iz   = n - abs(ix) - abs(iy)
 1          k(1) = ix*kv1(1) + iy*kv2(1) + iz*kv3(1)
            k(2) = ix*kv1(2) + iy*kv2(2) + iz*kv3(2)
            k(3) = ix*kv1(3) + iy*kv2(3) + iz*kv3(3)
            ikpt = pkpt(modulo(ix,nkpt3(1)),modulo(iy,nkpt3(2)),modulo(iz,nkpt3(3)))
            k2   = dot_product ( k , matmul(head,k) )
            k1   = sqrt(k2)
            if(k1<rrad+drad) then
              k1    = k1 - rrad
              rfac  = exp(-expo*k2)/k2 ; if(abs(k1)<drad) rfac = rfac * ( (k1/drad)**3 - 3*(k1/drad) + 2 ) / 4
              found = .true.
              div   = div - rfac * wtetra(:,ikpt)
            endif
            if(iz>0) then
              iz = -iz
              goto 1
            endif
          enddo
        enddo
        div_out = div_out + div
        n       = n + ifMpi(Msize,1)
      enddo
      Mpi( call Msum(div_out,0) )
      MnoR( div_out = 0 )
      do isite = 1,nsite
        rvec         = matmul(lat,rsite(:,isite))
        rnorm(isite) = sqrt( dot_product( rvec , matmul(invert(head),rvec) ) )
      enddo
      cfac = vol/(2*pi**2) * (rrad+drad)/nrad / sqrt(cdeterminant(head))
      do ix = 1,nrad ; Mcycle(ix)
        k1   = (rrad+drad)*(ix-0.5d0)/nrad
        k2   = k1 - rrad
        rfac = 1 ; if(abs(k2)<drad) rfac = ( (k2/drad)**3 - 3*(k2/drad) + 2 ) / 4
        k2   = k1*k1
        rfac = rfac * exp(-expo*k2)
        do isite = 1,nsite
          r              = rnorm(isite)
          si             = 1 ; if(r/=0) si = sin(k1*r) / (k1*r)
          div_out(isite) = div_out(isite) + si * rfac * cfac
        enddo
      enddo
      Mpi( call Msum(div_out) )
      end subroutine gamma_divergence2

      subroutine gamma_divergence3(div_out,wtetra,head)
      use global
# ifdef MPI
      use Mwrapper, only: Msum
# endif
      use, intrinsic :: iso_fortran_env
      implicit none
      complex_dp, intent(out) :: div_out(nsite)
      complex_dp, intent(in)  :: wtetra(nsite,nkpt),head(3,3)
      integer                 :: ix,iy,iz,n,ikpt,ll,l,opt,isite
      logical                 :: found
      real_dp,    parameter   :: expo = 5d-3
      integer,    parameter   :: nrad = 50000
      real_dp                 :: rrad,rrad2,drad,kv1(3),kv2(3),kv3(3),k(3),rvec(3),k2,k1,r,rfac,pre
      complex_dp              :: div(nsite),khk
      complex_dp              :: hharm(9)
      complex_dp, allocatable :: hinv(:),ylm(:)
      real_dp,    allocatable :: sphbes(:),integ(:)
      complex_dp              :: cerf
# ifdef MPI
      logical                 :: ldum
      integer                 :: Merr
# endif

      call matrix2harmonics(hharm(1),hharm(5),head) ; hharm(2:4) = 0
      call conjgharmonics(hharm,2)
      call invert_harm(hinv,ll,hharm,head)

      if(abs(head(1,1))+abs(head(2,2))+abs(head(3,3))>huge(0d0)) then
        if(hinv(1)/=0) Bug('hinv expected to be zero.')
        div_out = 0
        return
      endif

      if(all(rsite==0).and.sum(abs(head-head(1,1)*identity(3)))<1d-10) then ! scalar
        if(divergence==0) call gamma_divergence(.false.)
        div_out = divergence * hinv(1)/sqrt(4*pi)
        return
      else if(sum(abs(head))-abs(head(1,1))-abs(head(2,2))-abs(head(3,3))<1d-10) then ; opt = 1 ! diagonal matrix
      else                                                                            ; opt = 0 ! general matrix
      endif

      div_out = 0
      rrad    = sqrt(-log(expo)/expo)
      drad    = (rvol/nkpt)**(1d0/3)
      rrad2   = (rrad+drad)**2
      kv1     = rlat(:,1)/nkpt3(1)
      kv2     = rlat(:,2)/nkpt3(2)
      kv3     = rlat(:,3)/nkpt3(3)
      n       = 1 Mpi(+Mrank)
      found   = .true.
      do while(found)
        found = .false.
        div   = 0
        do ix = -n,n
          do iy = -(n-abs(ix)),n-abs(ix)
            iz   = n - abs(ix) - abs(iy)
 1          k(1) = ix*kv1(1) + iy*kv2(1) + iz*kv3(1)
            k(2) = ix*kv1(2) + iy*kv2(2) + iz*kv3(2)
            k(3) = ix*kv1(3) + iy*kv2(3) + iz*kv3(3)
            k2   = k(1)**2 + k(2)**2 + k(3)**2
            if(k2<rrad2) then
              if(opt==0) then ; khk = dot_product ( k , matmul(head,k) )
              else            ; khk = head(1,1)*k(1)**2 + head(2,2)*k(2)**2 + head(3,3)*k(3)**2
              endif
              ikpt  = pkpt(modulo(ix,nkpt3(1)),modulo(iy,nkpt3(2)),modulo(iz,nkpt3(3)))
              k1    = sqrt(k2)
              k1    = k1 - rrad
              rfac  = exp(-expo*k2) ; if(abs(k1)<drad) rfac = rfac * ( (k1/drad)**3 - 3*(k1/drad) + 2 ) / 4
              found = .true.
              div   = div - rfac * wtetra(:,ikpt) / khk
            endif
            if(iz>0) then
              iz = -iz
              goto 1
            endif
          enddo
        enddo
        div_out = div_out + div
        n       = n + ifMpi(Msize,1)
      enddo
      Mpi( call Msum(div_out,0) )
      MnoR( div_out = 0 )
      allocate ( sphbes(0:ll),ylm((ll+1)**2),integ(0:ll) )
      pre = vol/(2*pi**2) * (rrad+drad)/nrad
      do isite = 1,nsite
        rvec = matmul(lat,rsite(:,isite))
        r    = sqrt( sum(rvec**2) )
        if(r<1d-10) then
          ifR div_out(isite) = div_out(isite) + vol / (8*pi**2) / sqrt(expo) * cerf((1d0,0d0)*sqrt(expo)*rrad) * hinv(1)
        else
          call harmonicsr(ylm,rvec,ll)
          integ = 0
          do ix = Mrange1(nrad)
            k1    = (rrad+drad)*(ix-0.5d0)/nrad ; call sphbessel(sphbes,k1*r,ll)
            k2    = k1**2
            k1    = k1 - rrad
            rfac  = exp(-expo*k2) ; if(abs(k1)<drad) rfac = rfac * ( (k1/drad)**3 - 3*(k1/drad) + 2 ) / 4
            integ = integ + rfac * sphbes
          enddo
          do l = 0,ll,2
            div_out(isite) = div_out(isite) + pre*img**l * integ(l) * sum( hinv(l**2+1:(l+1)**2) * ylm(l**2+1:(l+1)**2) )
          enddo
        endif
      enddo

      deallocate(hinv,sphbes,ylm,integ)
      Mpi( call Msum(div_out) )
      end subroutine gamma_divergence3


      subroutine invert_harm(hinv,ll,hharm,head)
      use global,  only: pi
      use util,    only: chr
      use wrapper, only: trace
      implicit none
      integer,    intent(out)              :: ll
      complex_dp, intent(out), allocatable :: hinv(:)
      complex_dp, intent(in)               :: head(3,3),hharm(9)
      complex_dp                           :: h0,h1
      ll = 0
      h0 = sqrt(4*pi) * 3 / trace(head)
      do
        ll = ll + 2
        allocate ( hinv((ll+1)**2) )
        call expandfrac(hinv,[(1d0,0d0)*sqrt(4*pi)],hharm,ll,0,2)
        h1 = hinv(1)
        if(abs(h1-h0)<1d-6) then
          if(ll==2.and.abs(h1-h0)<1d-10) then ; ll = 0 ; deallocate( hinv ) ; allocate ( hinv(1) ) ; hinv = h0 ; endif
          return
        else if(ll==20) then ; Warn('Inversion not converged with l=20.')
        else if(ll==40) then ; Warn('Convergence not reached with l=40.') ; goto 1
        endif
        deallocate ( hinv )
        h0 = h1
      enddo
      if(ll>20) Warn('Inversion converged with l='//chr(ll))
 1    continue
      end subroutine invert_harm

      end

c -----------------------------

c     mode=0 : Transform from summation to normal order (13,42) -> (12,34)
c     mode=1 : Transform from normal to summation order (12,34) -> (13,42)
c     mode=2 : Transform to flip spins                  (12,34) -> (43,21)
      subroutine transpose2(mat,mode)
      use global, only: nwan
      use, intrinsic :: iso_fortran_env
      implicit none
      complex_dp, intent(inout) :: mat(nwan,nwan,nwan,nwan)
      integer,    intent(in)    :: mode
      integer                   :: iwan1,iwan2
      do iwan1 = 1,nwan
        do iwan2 = 1,nwan
          if(mode==0) then ; mat(iwan1,iwan2,:,:) = transpose(mat(iwan1,iwan2,:,:))
          else             ; mat(iwan1,:,:,iwan2) = transpose(mat(iwan1,:,:,iwan2))
          endif
        enddo
      enddo
      do iwan1 = 1,nwan
        do iwan2 = 1,nwan
          if     (mode==0) then ; mat(iwan1,:,:,iwan2) = transpose(mat(iwan1,:,:,iwan2))
          else if(mode==1) then ; mat(iwan1,iwan2,:,:) = transpose(mat(iwan1,iwan2,:,:))
          else if(mode==2) then ; mat(:,iwan1,iwan2,:) = transpose(mat(:,iwan1,iwan2,:))
          endif
        enddo
      enddo
      end


c BACKUP: Old gamma-divergence
# if 0

c     As gamma_divergence but integrand contains an exponential function
c     div_out(R) = 1/(8pi^3) INT exp(ikR)*Ylm(k)/k**2 d^3k - 1/nkpt SUM(k/=0) exp(ikR)*Ylm(k)/k**2
c     wtetra contains tetrahedron weights from tetrahedron6_init.
      subroutine gamma_divergence2(div_out,wtetra,lm)
      use global
# ifdef MPI
      use Mwrapper, only: Msum
# endif
      use, intrinsic :: iso_fortran_env
      implicit none
      Mpi( include 'mpif.h' )
      complex_dp, intent(out) :: div_out(nsite)
      complex_dp, intent(in)  :: wtetra(nsite,nkpt)
      integer,    intent(in)  :: lm
      complex_dp, allocatable :: ylm(:)
      real_dp,    allocatable :: sphbes(:)
      integer                 :: ix,iy,iz,n,ikpt,l
      logical                 :: found
      real_dp,    parameter   :: expo = 5d-3
      integer,    parameter   :: nrad = 50000
      real_dp                 :: rrad,drad,kv1(3),kv2(3),kv3(3),k(3),k2,k1,fsmooth,r
      complex_dp              :: div(nsite),fac(nsite)
      real_dp                 :: rnorm(nsite)
# ifdef MPI
      logical                 :: ldum
      integer                 :: m,mm,m1,m2,Merr
# endif
      l = 0 ; do while((l+1)**2<lm) ; l = l + 1 ; enddo ; allocate ( ylm((l+1)**2),sphbes(0:l) ) ; ylm(1) = 1/sqrt(4*pi)
      div_out = 0
      rrad    = sqrt(-log(expo)/expo)
      drad    = (rvol/nkpt)**(1d0/3)
      kv1     = rlat(:,1)/nkpt3(1)
      kv2     = rlat(:,2)/nkpt3(2)
      kv3     = rlat(:,3)/nkpt3(3)
      n       = 1
      found   = .true.
      Mpi( mm = -1 )
      do while(found)
        found = .false.
        div   = 0
        do ix = -n,n
          do iy = -(n-abs(ix)),n-abs(ix)
            McycleP(mm)
            iz   = n - abs(ix) - abs(iy)
 1          k(1) = ix*kv1(1) + iy*kv2(1) + iz*kv3(1)
            k(2) = ix*kv1(2) + iy*kv2(2) + iz*kv3(2)
            k(3) = ix*kv1(3) + iy*kv2(3) + iz*kv3(3)
            k2   = k(1)**2   + k(2)**2   + k(3)**2
            k1   = sqrt(k2)
            ikpt = pkpt(modulo(ix,nkpt3(1)),modulo(iy,nkpt3(2)),modulo(iz,nkpt3(3)))
            if(k1<rrad+drad) then
              if(l>0) call harmonicsr(ylm,k,l)
              k1 = k1 - rrad
              if(abs(k1)<drad) then ; fsmooth = ( (k1/drad)**3 - 3*(k1/drad) + 2 ) / 4
              else                  ; fsmooth = 1
              endif
              found = .true.
              div   = div - exp(-expo*k2)/k2 * fsmooth * ylm(lm) * wtetra(:,ikpt)
            endif
            if(iz>0) then
              iz = -iz
              goto 1
            endif
          enddo
        enddo
        Mpi( call mpi_allreduce(found,ldum,1,mpi_logical,mpi_lor,Mcomm,Merr) ; found = ldum )
        div_out = div_out + div
        n       = n + 1 !; Rif(l==4) write(*,*) 'n',n,found
      enddo
      Mpi( call Msum(div_out,0) )
      MnoR( div_out = 0 )
      Rwrite(*,*) 'bef',div_out
      Rwrite(*,*) '...'
      do isite = 1,nsite
        k            = matmul(lat,rsite(:,isite)) ; call harmonicsr(ylm,k,l)
        rnorm(isite) = sqrt(sum(k**2))
        fac(isite)   = vol/(8*pi**3) * (rrad+drad)/nrad * 4*pi*img**l * ylm(lm)
      enddo
      do ix = 1,nrad ; Mcycle(ix)
        k1 = (rrad+drad)*(ix-0.5d0)/nrad
        k2 = k1 - rrad
        if(abs(k2)<drad) then ; fsmooth = ( (k2/drad)**3 - 3*(k2/drad) + 2 ) / 4
        else                  ; fsmooth = 1
        endif
        k2 = k1*k1
        do isite = 1,nsite
          r              = rnorm(isite) ; call sphbessel(sphbes,k1*r,l)
          div_out(isite) = div_out(isite) + exp(-expo*k2) * sphbes(l) * fac(isite) * fsmooth
        enddo
      enddo
      Mpi( call Msum(div_out) )
      Rwrite(*,*) div_out
      deallocate ( ylm,sphbes )

      end subroutine gamma_divergence2


c     As gamma_divergence but integrand contains an exponential function
c     div_out(R) = 1/(8pi^3) INT exp(ikR) 1/k**2 d^3k - 1/nkpt SUM(k/=0) exp(ikR) 1/k**2
c     wtetra contains tetrahedron weights from tetrahedron6_init.
      subroutine gamma_divergence2(div_out,wtetra)
      use global
# ifdef MPI
      use Mwrapper, only: Msum
# endif
      use, intrinsic :: iso_fortran_env
      implicit none
      Mpi( include 'mpif.h' )
      real_dp,    intent(out) :: div_out(nsite)
      complex_dp, intent(in)  :: wtetra(nsite,nkpt)
      integer                 :: ix,iy,iz,n,ikpt
      logical                 :: found
      real_dp,    parameter   :: expo = 5d-3
      real_dp                 :: div(nsite),rrad,drad,kv1(3),kv2(3),kv3(3),k(3),k2,k1,fsmooth,rfac,si,r
# ifdef MPI
      logical                 :: ldum
      integer                 :: m,mm=0,m1,m2,Merr
# endif
      div_out = 0
      rrad    = sqrt(-log(expo)/expo)
      drad    = (rvol/nkpt)**(1d0/3) !; drad = 0 ; write(*,*) 'drad set to zero!'
      kv1     = rlat(:,1)/nkpt3(1)
      kv2     = rlat(:,2)/nkpt3(2)
      kv3     = rlat(:,3)/nkpt3(3)
      n       = 1
      found   = .true.
      do while(found)
        found = .false.
        div   = 0
        do ix = -n,n
          do iy = -(n-abs(ix)),n-abs(ix)
            McycleP(mm)
            iz   = n - abs(ix) - abs(iy)
 1          k(1) = ix*kv1(1) + iy*kv2(1) + iz*kv3(1)
            k(2) = ix*kv1(2) + iy*kv2(2) + iz*kv3(2)
            k(3) = ix*kv1(3) + iy*kv2(3) + iz*kv3(3)
            k2   = k(1)**2   + k(2)**2   + k(3)**2
            k1   = sqrt(k2)
            ikpt = pkpt(modulo(ix,nkpt3(1)),modulo(iy,nkpt3(2)),modulo(iz,nkpt3(3)))
            if(k1<rrad+drad) then
              k1 = k1 - rrad
              if(abs(k1)<drad) then ; fsmooth = ( (k1/drad)**3 - 3*(k1/drad) + 2 ) / 4
              else                     ; fsmooth = 1
              endif
              found = .true.
              div   = div - exp(-expo*k2)/k2 * fsmooth * wtetra(:,ikpt)
            endif
            if(iz>0) then
              iz = -iz
              goto 1
            endif
          enddo
        enddo
        Mpi( call mpi_allreduce(found,ldum,1,mpi_logical,mpi_lor,Mcomm,Merr) ; found = ldum )
        div_out = div_out + div
        n       = n + 1
      enddo
      Mpi( call Msum(div_out) )
      do isite = 1,nsite
        div(isite) = sqrt(sum(matmul(lat,rsite(:,isite))**2)) ! norm of rsite
      enddo
      rfac = vol/(2*pi**2) * rrad/10000
      do ix = 1,10000
        k1 = rrad*(ix-0.5d0)/10000
        k2 = k1*k1
        do isite = 1,nsite
          r = div(isite)
          if(r==0) then ; si = 1
          else          ; si = sin(k1*r) / (k1*r)
          endif
          div_out(isite) = div_out(isite) + exp(-expo*k2) * si * rfac
        enddo
      enddo
      end subroutine gamma_divergence2
# endif

c -----------------------------

