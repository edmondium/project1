c--------------------------------------------------------------------------------
c Copyright (c) 2021 Peter Grünberg Institut, Forschungszentrum Jülich, Germany
c This file is part of SPEX and available as free software under the conditions
c of the MIT license as expressed in the LICENSE file in more detail.
c--------------------------------------------------------------------------------

c Calculates the susceptibility matrices (suscepc) at k point ikpt and frequencies frq(:nfrq)
c projected onto the Wannier functions:
c
c  ss'          1     occ uno  / U*(q+km',s,n1) U(q+km',s,n3) U(km,s',n2) U*(km,s',n4)
c K   (q,w) = - - SUM SUM SUM <  -----------------------------------------------------       <--- iterm = 1, ispin1 = s', ispin2 = s
c  n1n2,n3n4    N  k   m   m'  \            w + e(km,s') - e(q+km',s) + i\eta
c
c                                U*(q+km,s,n1) U(q+km,s,n3) U(km',s',n2) U*(km',s',n4) \
c                             -  -----------------------------------------------------  >    <--- iterm = 2, ispin1 = s, ispin2 = s'
c                                           w + e(km',s') - e(q+km,s) - i\eta          /
c                s      s'
c           = i G (13) G (42+)                                                               <--- non-interacting two-particle propagator
c
c ikpt = q
c frq  = w
c spin = 0 =>                "uu+dd"
c spin = 1 => s = 1, s' = 1  "uu"
c spin = 2 => s = 2, s' = 2  "dd"
c spin = 3 => s = 1, s' = 2  "ud"
c spin = 4 => s = 2, s' = 1  "du"
c

# include "cppmacro.h"

      subroutine susceptibility_wannier ( ikpt,spin,symon,nfrq,frq,suscepc )

      use global
      use wrapper
      use util
      use file
      use key
      use, intrinsic :: iso_fortran_env
      Mpi ( use Mwrapper )
      Mpi2( use, intrinsic :: iso_c_binding )

      implicit none
      integer,    intent(in)  :: ikpt,spin
      logical,    intent(in)  :: symon

      integer,    intent(in)  :: nfrq
      complex_dp, intent(in)  :: frq(nfrq)
      complex_dp, intent(out) :: suscepc(nwan*nwan,nwan*nwan,nfrq)
      complex_dp, allocatable :: help0(:,:),help1(:,:)
      real_dp,    allocatable :: wghtSr(:,:),wghtSi(:,:)
      real_dp,    pointer_cnt :: wtetra(:,:)
      real_dp,    allocatable :: re(:),im(:)
      real_dp                 :: rdum
      integer,    pointer_cnt :: index(:,:)
      integer,    allocatable :: nindex(:)
      integer                 :: ispin1,ispin2,ikpt0,ikpt1,ikpt2,ifrq,ifrqs,bandi,bandf,iunit,iterm
      integer                 :: i,j,k,nu,nu1,n,isym,m1,m2
      integer                 :: n1,n2,n3,n4
      integer                 :: kpt1(nkpt),kpt1p(nkpt),nkpts(nkpt),sym1(nsym),nsym1,nkpt1,kptsum,bandow,wanbando,wanbandu
      real                    :: time1,time2
      complex_dp              :: cdum,cdum1

      integer                 :: nfrqs
      real_dp,    allocatable :: frqs(:)
      complex_dp, allocatable :: sf(:)

      real_dp, save           :: disorder,wghtthr,fspec(2)
      logical, save           :: firstcall = .true.

# ifdef MPI
      type(c_ptr)             :: ptr
      integer                 :: win_wtetra,win_index,Merr
# endif

# include "interface/getfrqs.inc"

c      if(ikpt==nkpt+1.and.symon) Error('WANNIER->IRREP not implemented for additional k point.')

c     Input from spex.inp
      if(firstcall) then
        firstcall = .false.
        Rbegin
        call getkey(inp,'WGHTTHR',  wghtthr,  section='SUSCEP',  default=1d-10, mini=0d0)
        call getkey(inp,'DISORDER', disorder, section='SUSCEP',  default=0d0,   mine=0d0)
        call getfspec(fspec,0,nfrq,frq,.true.)
        if(fspec(1)==0) Error('Hilbert mesh needed (HILBERT).')
        if(disorder/=0) disorder = 1/(2*disorder) ! eta = 1/(2*tau)
        Rend
        Mpi( call Mcast(wghtthr) ; call Mcast(disorder) )
      endif

      if(nfrq<=0) Error('no frequencies.')
      suscepc = 0

c      if(symon)          Error('symmetry disabled (symon).')
      if(wanbandi>bando) Error('no occupied states included in Wannier construction.')
      if(wanbandf<bandu) Error('no unoccupied states included in Wannier construction.')
      Rif(any(real(frq)<0)) Info('frequency has negative real part, which we accept ...')        

      Rwrite(6,'(/A)') '---------'

      Rwrite(6,'(/A)')           'Calculation of susceptibility'
      Rwrite(6,'(/A,I7)')        'Order of array suscep:',nwan
      Rwrite(6,'(A,F7.1," MB")') 'Size of array suscep:   ',(nwan**4*nfrq*16d0)/megabyte

c
c     Get irreducible k points kpt1(:nkpt1) (wrt current k point ikpt)
      if(symon) then ; call getkpt1       (kpt1,nkpt1,nkpts,sym1,nsym1,ikpt,0, ifMpi(Mrank==0,.true.) )
      else           ; call getkpt1_fullBZ(kpt1,nkpt1,nkpts,sym1,nsym1,0)
      endif
      ! Construct kpt1p (parent k points of EIBZ)
      do i = 1,nkpt1
        do isym = 1,nsym1
          kpt1p(kptsym(kpt1(i),sym1(isym))) = i
        enddo
      enddo

c
c     Calculate packet size (leading dimension) of cprod,wtetra,... according to memprod (->nu)
      rdum = ( maxmem - mem ) Mpi(/Nsize) ; Rwrite(6,'(A,F13.1,A)') 'Auxiliary storage:',rdum/megabyte,' MB'

      wanbando = min(wanbandf,bando)
      wanbandu = max(wanbandi,bandu)
      bandow   = min(wanbandf,bando) - wanbandi + 1

c
c                        iterm=1             iterm=2
c     ispin1 :   occupied states   unoccupied states
c     ispin2 : unoccupied states     occupied states
      if     (spin==0) then ; ispin1 = 1 ; ispin2 = 1
      else if(spin==1) then ; ispin1 = 1 ; ispin2 = 1 ! uu
      else if(spin==2) then ; ispin1 = 2 ; ispin2 = 2 ! dd
      else if(spin==3) then ; ispin1 = 2 ; ispin2 = 1 ! ud
      else if(spin==4) then ; ispin1 = 1 ; ispin2 = 2 ! du
      else                  ; Error('unknown spin index.')
      endif

      iterm = 0
 1    iterm = iterm + 1

c
c     Define Hilbert mesh
      Rbegin
      ! Get frequency range -> re(1)..re(2)
      allocate ( re(2) )
      if(gauss(1)/=0) then
        re(1) = 0 ; re(2) = maxval(ene(wanbandf,:,:)) - minval(ene(wanbandi,:,:))
      else
        call tetrahedron5_init(re,rdum,0,ikpt,kpt1,kpt1p,nkpt1,iterm,ispin1,ispin2,wanbandi,wanbando,wanbandu,wanbandf,1,2,1,1)
      endif
c      re(2) = maxval(ene(:maxeband,:nkpt,:)) - minval(ene(:maxeband,:nkpt,:)) ; write(*,*) 're(2) replaced'
      ! Define
      call getfrqs(frqs,nfrqs,fspec,re(1),re(2),.true.)
c      write(*,*) re(2),maxval(ene(wanbandf,:,:)) - minval(ene(wanbandi,:,:))
c      write(*,*) nfrqs,frqs ; stop
c      write(*,*) re(2)
c      write(*,*) frqs(nfrqs) + (frqs(nfrqs-1)-frqs(nfrqs-2))**2/(frqs(nfrqs)-frqs(nfrqs-2))

      deallocate ( re )
      if(iterm==1) then
        write(6,'(A,I7'NoA) 'Number of frequencies:',nfrq
        if(nfrqs/=0) write(6,'(I6'NoA) nfrqs
        write(6,*)
      endif
      ! Check
      if(all(imag(frq)==0)) then
        i = 0
        do ifrqs = 1,nfrqs-1
          i = max(i,count(real(frq)>frqs(ifrqs).and.real(frq)<frqs(ifrqs+1)))
        enddo
        if(i>10) Warn('Coarse Hilbert mesh. Line segments in imaginary part of length: '//trim(chr(i)))
      endif
      Rend
      Mpi( call Mcast(nfrqs) ; call Mcastl(frqs) )

      if(nfrqs/=0) then
        allocate ( wghtSr(nfrqs,nfrq),wghtSi(nfrqs,nfrq) )
        call getwghtS(wghtSr,wghtSi,frq,nfrq,frqs,nfrqs,disorder,0,iterm)
        wghtSr = -wghtSr ! K defined negative (see above)
        wghtSi = -wghtSi !
        if(wrtinfo) then ; allocate(sf(nfrqs)) ; sf = 0 ; endif
      endif

c
c     Define current packet
      nu = int(rdum/((16+(8d0*nfrqs+4*nfrqs) Mpi(/Nsize) )*nkpt1*bandow))
      if(nu==0) Error('Maximal allowed storage exceeded. Increase MEM or reduce number of frequencies!')

      if(nspin/=1) then
        Rif(ispin1==1.and.ispin2==1) write(6,'(/A'NoA) 'up/up    '
        Rif(ispin1==2.and.ispin2==2) write(6,'(/A'NoA) 'down/down'
        Rif(ispin1==1.and.ispin2==2) write(6,'(/A'NoA) 'up/down  '
        Rif(ispin1==2.and.ispin2==1) write(6,'(/A'NoA) 'down/up  '
      else
        Rwrite(6,'(/A'NoA) '         '
      endif
      if(nfrqs==0) then
        Rwrite(6,'(A)') 'skipped because Hilbert mesh has no points.'
        goto 2
      endif
      Rwrite(6,'(A)')          '        Bands | Timings'
      Rwrite(6,'(A)') '                       |'Mpi(//'    sync')//'   wghts'Mpi(//'    sync')//'  suscep'

c
c     Loop over unoccupied states in packets of nu bands

      do bandi = wanbandu,wanbandf,nu
        bandf = min(bandi+nu-1,wanbandf)
        nu1   = bandf-bandi+1 ! last packet can be smaller than nu

        Rwrite(6,'(I6,'' ->'',I4,''  ('',I4,'' )  '''NoA) bandi,bandf,nu1

        Nallocate( wtetra,(S_ nu1*bandow*nkpt1,nfrqs S_) )
        Nallocate( index, (S_ nu1*bandow*nkpt1,nfrqs S_) )
        allocate ( re     (   nu1*bandow*nkpt1)       )
        allocate ( im     (   nu1*bandow*nkpt1)       )
        allocate ( nindex (                    nfrqs) )

c
c       Loop over spectral-function frequencies and calculate integration weights
        call cpu_time(time1)
        Mpi( call mpi_barrier(Mcomm,i);call cpu_time(time2);Rwrite(6,'(F8.2'NoA) time2-time1;call cpu_time(time1) )
        Nfence(wtetra)
        ifO wtetra = 0
        Nfence(wtetra)
        if(gauss(1)==0) then
          n1=wanbandi;n2=wanbando;m1=bandi;m2=bandf
          Mpi( call distribute_bands(n1,n2,m1,m2) )
          call tetrahedron5_init(wtetra,frqs,nfrqs,ikpt,kpt1,kpt1p,nkpt1,iterm,ispin1,ispin2,n1,n2,m1,m2,
     &                           wanbandi,wanbando,bandi,bandf)
        else
          Ocall gauss4_init     (wtetra,frqs,nfrqs,ikpt,kpt1,kpt1p,nkpt1,iterm,ispin1,ispin2,wanbandi,wanbando,bandi,bandf)
        endif
        call cpu_time(time2);Rwrite(6,'(F8.2'NoA) time2-time1;call cpu_time(time1)
        Nfence(wtetra)
        Mpi( call cpu_time(time2);Rwrite(6,'(F8.2'NoA) time2-time1;call cpu_time(time1) )
c        Mpi( if(gauss(1)==0) call distribute_wtetra(wtetra) )
        ! Remove zero weights and reorder.
        nindex = 0
        do ifrqs = Nrange(1,nfrqs)
          n = 0
          do i = 1,size(wtetra,1)
            if(wtetra(i,ifrqs)>wghtthr) then
              n = n + 1 ; index(n,ifrqs) = i ; wtetra(n,ifrqs) = wtetra(i,ifrqs) / nsym1
            endif
          enddo
          nindex(ifrqs) = n
        enddo
        Nfence(wtetra)
        Nfence(index)
        Mpi( call Msum(nindex,comm=Ncomm) )

c
c       Finally calculate contribution of current packet
        ! loop over elements of suscepc
        n = size(wtetra,1)
        do m1 = 1,nwan*nwan         ; n1 = mod(m1-1,nwan)+1 ; n2 = (m1-n1)/nwan+1
          do m2 = Mrange(1,n2*nwan) ; n3 = mod(m2-1,nwan)+1 ; n4 = (m2-n3)/nwan+1
            k = 0
            do ikpt0 = 1,nkpt1
              ikpt1 = kpt1(ikpt0)
              ikpt2 = kptsum(ikpt1,ikpt)
              if(iterm==2) then
                do i = wanbandi,wanbando
                  cdum1 = uwan(i,n3,ikpt2,ispin1) * conjg ( uwan(i,n1,ikpt2,ispin1) )
                  do j = bandi,bandf
                    k     = k + 1
                    cdum  = uwan(j,n2,ikpt1,ispin2) * conjg( uwan(j,n4,ikpt1,ispin2) ) * cdum1
                    re(k) = real(cdum)
                    im(k) = imag(cdum)
                  enddo
                enddo
              else
                do i = wanbandi,wanbando
                  cdum1 = uwan(i,n2,ikpt1,ispin1) * conjg ( uwan(i,n4,ikpt1,ispin1) )
                  do j = bandi,bandf
                    k     = k + 1
                    cdum  = uwan(j,n3,ikpt2,ispin2) * conjg( uwan(j,n1,ikpt2,ispin2) ) * cdum1
                    re(k) = real(cdum)
                    im(k) = imag(cdum)
                  enddo
                enddo
              endif
            enddo
            call integrate_spectralf(n2/=n4)

          enddo
        enddo
        deallocate ( re,im,nindex )
        Ndeallocate( wtetra )
        Ndeallocate( index )
        call cpu_time(time2);Rwrite(6,'(F8.2)') time2-time1;call cpu_time(time1)

      enddo ! End loop over packets

      Rif(wrtinfo.and.nfrqs/=0) then
        iunit = fopen('spex.hil',status='unknown')
        write(iunit,'(A)') '# Hilbert frequency mesh'
        do ifrqs = 1,nfrqs
          write(iunit,'(F20.10)') frqs(ifrqs)
        enddo
        call fclose(iunit)
      endif

 2    deallocate ( frqs )
      if(allocated(wghtSr)) deallocate ( wghtSr,wghtSi )
      if(allocated(sf))     deallocate ( sf )

c
c     Loop over spin down/down if necessary
      if(spin==0) then
        Error('spin=0 never checked.')
        if(nspin==2) then
          if(ispin1==1) then ; ispin1 = 2 ; ispin2 = 2 ; goto 1 ; endif
        else
          if(nfrq>0) suscepc = suscepc * 2
        endif
      endif

c
c     Loop over spin down/up (up/down) if necessary
      if(spin==3) then
        if(ispin1==2) then ; ispin1 = 1 ; ispin2 = 2 ; goto 1 ; endif
      else if(spin==4) then
        if(ispin1==1) then ; ispin1 = 2 ; ispin2 = 1 ; goto 1 ; endif
      endif

      Mpi( call Msum(suscepc) )

c
c     Symmetrize if symmetry is used (note that ispin1 and ispin2 are switched after the second loop)
      if(symon) then
        if(ispin1==ispin2) Error('Symmetrization not implemented for spin=0-2.')
        Rwrite(6,'(/A'NoA) 'Add contribution of symmetry-equivalent k points... '
        call cpu_time(time1)
        allocate ( help0(nwan*nwan,nwan*nwan) )
        allocate ( help1(nwan*nwan,nwan*nwan) )
        do ifrq = Mrange1(nfrq)
          help0 = suscepc(:,:,ifrq)
          do isym = 2,nsym1
            help1 = help0
            if(ikpt==nkpt+1) then
              call convolute(help1,      irrep_wan1(:,:,sym1(isym),ispin1) ,nwan,1)
              call convolute(help1,conjg(irrep_wan1(:,:,sym1(isym),ispin1)),nwan,3)
            else
              call convolute(help1,      irrep_wan(:,:,sym1(isym),ispin1) ,nwan,1)
              call convolute(help1,conjg(irrep_wan(:,:,sym1(isym),ispin1)),nwan,3)
            endif
            call convolute(help1,conjg(irrep_wan(:,:,sym1(isym),ispin2)),nwan,2)
            call convolute(help1,      irrep_wan(:,:,sym1(isym),ispin2) ,nwan,4)
            suscepc(:,:,ifrq) = suscepc(:,:,ifrq) + help1
          enddo
        enddo
        deallocate ( help0,help1 )
        Rcall cpu_done(time1)
        MrangeDistr(suscepc(:,:,McolD1(nfrq,i)),i)
      endif

      Rwrite(6,'(/A/)') '---------'

c     ---------------

      contains

c     Calculate spectral function at frequencies frqs (->reS+img*imS) and integrate over spectral function (-> out).
      subroutine integrate_spectralf(offdiag)
      implicit none
      logical, intent(in) :: offdiag
      integer             :: ifrqs,k,ifrq,i1
      real_dp             :: reS(nfrqs),imS(nfrqs),rdum1,rdum2
      real_dp             :: ddot
      do ifrqs = 1,nfrqs
        rdum1 = 0d0
        rdum2 = 0d0
        do k = 1,nindex(ifrqs)
          i1    = index(k,ifrqs)
          rdum1 = rdum1 + re(i1) * wtetra(k,ifrqs)
          rdum2 = rdum2 + im(i1) * wtetra(k,ifrqs) !; if(abs(im(i1))>1d-6) write(*,*) im(i1),wtetra(k,ifrqs)
        enddo
        reS(ifrqs) = rdum1
        imS(ifrqs) = rdum2 !; write(*,*) rdum2 !if(abs(rdum2)>1d-9) then ; write(*,*) rdum2 ; read(*,*) ; endif
      enddo

      do ifrq = 1,nfrq
        rdum1     = ddot(nfrqs,reS,1,wghtSr(1,ifrq),1)
        rdum2     = ddot(nfrqs,imS,1,wghtSr(1,ifrq),1)
                    suscepc(m1,m2,ifrq) = suscepc(m1,m2,ifrq) + rdum1 + img*rdum2
        if(offdiag) suscepc(m2,m1,ifrq) = suscepc(m2,m1,ifrq) + rdum1 - img*rdum2
        rdum1     = ddot(nfrqs,reS,1,wghtSi(1,ifrq),1)
        rdum2     = ddot(nfrqs,imS,1,wghtSi(1,ifrq),1)
                    suscepc(m1,m2,ifrq) = suscepc(m1,m2,ifrq) + ( rdum1 + img*rdum2 ) * img
        if(offdiag) suscepc(m2,m1,ifrq) = suscepc(m2,m1,ifrq) + ( rdum1 - img*rdum2 ) * img
      enddo

      end subroutine integrate_spectralf

c     ---------------

# ifdef MPI
      subroutine distribute_bands(n1,n2,m1,m2) ! distribute over processes in local node
      implicit none
      integer, intent(inout) :: n1,n2,m1,m2
      integer                :: prim(ceiling(log(1d0*Nsize)/log(2d0)))
      integer                :: size1,size2,rank1,rank2
      integer                :: i,j,n
c      if(option==2.or.Nsize==1) return
      if(Nsize==1) return
      call primfac(prim,n,Nsize)
      if(product(prim(:n))/=Nsize) Bug('Prime factorization failed.')
      j = 1
      do i = 1,n
        if(product(prim(:i))>product(prim(i+1:n))) then ; j = i ; exit ; endif
      enddo
      if(j/=1) then
        if(abs(product(prim(:j))  -product(prim(j+1:n)))>
     &     abs(product(prim(:j-1))-product(prim(  j:n)))) j = j-1
      endif
      size1 = product(prim(:j))
      size2 = product(prim(j+1:n))
      rank1 = mod ( Nrank , size1 )
      rank2 = Nrank / size1
      n     = n2-n1+1
      n2    = n1+(1+rank1)*n/size1-1
      n1    = n1+   rank1 *n/size1
      n     = m2-m1+1
      m2    = m1+(1+rank2)*n/size2-1
      m1    = m1+   rank2 *n/size2
      end subroutine distribute_bands

#   if 0
      subroutine distribute_wtetra(wtetra_in)
      implicit none
      real_dp :: wtetra_in(*)
c      if(option==2) then ; call distribute_wtetra_opt2(wtetra_in)
c      else               ; call distribute_wtetra0    (wtetra_in)
c      endif
      call distribute_wtetra0    (wtetra_in)
      end subroutine distribute_wtetra

      subroutine distribute_wtetra_opt2(wtetra_in)
      implicit none
      real_dp :: wtetra_in(nu1*bandow*Mnum1(nkpt1),nfrqs)
      real_dp :: hlp      (nu1*bandow*Mnum1(nkpt1),nfrqs)
      integer :: i
      hlp                                          = wtetra_in
      wtetra(nu1*bandow*Mcol1(nkpt1)*nu1*bandow,:) = hlp
      MrangeDistr(wtetra(nu1*bandow*McolD1(nkpt1,i)*nu1*bandow,:),i)
      end subroutine distribute_wtetra_opt2

      subroutine distribute_wtetra0(wtetra_in)
      real_dp, intent(in) :: wtetra_in(m1:m2,n1:n2,nkpt1,nfrqs)
      real_dp             :: hlp      (m1:m2,n1:n2,nkpt1,nfrqs)
      hlp = wtetra_in
      call distribute_wtetra1(wtetra,hlp)
      end subroutine distribute_wtetra0

      subroutine distribute_wtetra1(wtetra,hlp)
      implicit none
      real_dp, intent(out) :: wtetra(bandi:bandf,wanbandi:wanbando,nkpt1,nfrqs)
      real_dp, intent(in)  :: hlp   (m1:m2,n1:n2,nkpt1,nfrqs)
      real_dp, allocatable :: hlp1  (:,:,:,:)
      integer              :: mm1,mm2,nn1,nn2,rank
      do rank = 0,Msize-1
        if(rank==Mrank) then
          call Mcast(n1,Mrank)  ; call Mcast(n2,Mrank)
          call Mcast(m1,Mrank)  ; call Mcast(m2,Mrank)
          call Mcast(hlp,Mrank) ; wtetra(m1:m2,n1:n2,:,:) = hlp
        else
          call Mcast(nn1,rank)  ; call Mcast(nn2,rank)
          call Mcast(mm1,rank)  ; call Mcast(mm2,rank) ; allocate(hlp1(mm1:mm2,nn1:nn2,nkpt1,nfrqs))
          call Mcast(hlp1,rank) ; wtetra(mm1:mm2,nn1:nn2,:,:) = hlp1 ; deallocate(hlp1)
        endif
      enddo
      end subroutine distribute_wtetra1
#   endif

# endif

c     ---------------

      end

c     ---------------

      subroutine convolute(a,p,n,i)
      use, intrinsic :: iso_fortran_env
      implicit none
      integer,    intent(in)    :: n,i
      complex_dp, intent(in)    :: p(n,n)
      complex_dp, intent(inout) :: a(n,n,n,n)
      complex_dp                :: b(n,n,n,n)
      integer                   :: n1,n2,n3,n4
      b = a
      if(i==1) then
        do n4 = 1,n ; do n3 = 1,n ; do n2 = 1,n ; do n1 = 1,n
          a(n1,n2,n3,n4) = sum(p(n1,:)*b(:,n2,n3,n4))
        enddo ; enddo ; enddo ; enddo
      else if(i==2) then
        do n4 = 1,n ; do n3 = 1,n ; do n2 = 1,n ; do n1 = 1,n
          a(n1,n2,n3,n4) = sum(p(n2,:)*b(n1,:,n3,n4))
        enddo ; enddo ; enddo ; enddo
      else if(i==3) then
        do n4 = 1,n ; do n3 = 1,n ; do n2 = 1,n ; do n1 = 1,n
          a(n1,n2,n3,n4) = sum(p(n3,:)*b(n1,n2,:,n4))
        enddo ; enddo ; enddo ; enddo
      else if(i==4) then
        do n4 = 1,n ; do n3 = 1,n ; do n2 = 1,n ; do n1 = 1,n
          a(n1,n2,n3,n4) = sum(p(n4,:)*b(n1,n2,n3,:))
        enddo ; enddo ; enddo ; enddo
      else
        Bug('wrong input i.')
      endif
      end

