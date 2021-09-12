c--------------------------------------------------------------------------------
c Copyright (c) 2021 Peter Grünberg Institut, Forschungszentrum Jülich, Germany
c This file is part of SPEX and available as free software under the conditions
c of the MIT license as expressed in the LICENSE file in more detail.
c--------------------------------------------------------------------------------

# include "cppmacro.h"
# include "jobtype.h"

c# define additional

c Solves quasiparticle equation
c
c [ hks + sigmax + sigmac(eqp) - vxc ] PHI = eqp * PHI
c
c
c (1) Linearized one-shot solution:
c
c     The quasiparticle equation is solved perturbatively, i.e., using the KS states as approximate quasiparticle amplitudes PHI:
c     eqp = eks + sigmax + sigmac(eqp) - vxc    ! approximation: sigmac(eqp) = sigmac(eks) + dsigmac(eqp) (eqp-eks)
c         = eks + [1-dsigmac(eks)]**(-1) * [ sigmax + sigmac(eks) - vxc ]
c
c (2) Direct one-shot solution:
c
c     The nonlinear equation is solved by iteration. Let eqp1 be the old and eqp2 the new energy.
c     eqp2 = eks + sigmax + sigmac(eqp2) - vxc  ! approximation: sigmac(eqp2) = sigmac(eqp1) + dsigmac(eqp1) (eqp2-eqp1)
c          = [1-dsigmac(eqp1)]**(-1) * [ eks + sigmax + sigmac(eqp1) - vxc - dsigmac(eqp1) eqp1 ]
c     The start value is eqp1=eks so that the linearized solution (1) corresponds to the first iteration.
c
c     Notes about CONTOUR calculations:
c     - The direct solution requires a real mesh to be defined ({...} or [{...}]).
c     - Since we employ a real mesh for sigmac, we can interpolate sigmac(eqp) only for real eqp. So, to start
c       a new iteration we set eqp1=real(eqp2). We could also define sigmac(eqp) for complex eqp by:
c       sigmac(eqp) = sigmac(real(eqp)) + dsigmac(real(eqp)) * (eqp-real(eqp))
c       dsigmac(eqp) = dsigmac(real(eqp)),
c       but it is easy to see that this would lead to exactly the same eqp2 with the above formula.
c
c (3) Full quasiparticle solution:
c
c     The quasiparticle Hamiltonian   hqp = diagonal(hks) + sigmax + sigmac(eqp) - vxc   is diagonalized in
c     the basis of KS states (iteratively for each state until self-consistency in the energy eqp is reached).
c     Only a direct solution is calculated. Starting from eqp1=eks, we get a new energy eqp2 by
c     eqp2 = [1-dsigmac(eqp1)]**(-1) * [ eval - dsigmac(eqp1)*eqp1 ],
c     where eval is one of the eigenvalues of hqp.
c     CONTOUR calculations are only possible with [{...}]. Since the diagonalization changes the eigenstates from
c     iteration to iteration, it does make a difference whether we set eqp1=real(eqp2) or use the sigmac(eqp) for
c     complex eqp [see note in (2)]. The difference is very small, though, and we use eqp1=real(eqp2) for consistency.
c
c     Degeneracies should be preserved. However, hqp is complex. Then, the following applies
c
c     P = symop w/o  time reversal => If c is a right-eigenvector (hqp*c=e*c), then Pc is a right-eigenvector:       hqp * Pc  = e * Pc
c     P = symop with time reversal => If c is a right-eigenvector (hqp*c=e*c), then Pc is a  left-eigenvector: conjg(Pc) * hqp = e * conjg(Pc)
c     (Here, conjg includes transposition. For a Hermitian matrix, right- and left-eigenvectors are identical.)
c
c     In SOC calculations with inversion symmetry, all states are (at least) doubly degenerate, and the pairs of
c     KS states are related by trs*inv, i.e., the combination of time reversal and inversion. This is not so for
c     hqp, although there are still double degeneracies. So, while the quasiparticle energies are doubly degenerate,
c     the corresponding vxc values (etc., note that vxc is Hermitian) can differ (because the two eigenvectors are
c     not related by trs*inv).
c
c     ALIGNVXC: Instead of E = e + SIGMA(E-ef) - vxc, we solve E = e + SIGMA(E-D-ef) - vxc. Interpreting D as a shift of vxc, we have
c               vxc'=vxc+D and consequently e'=e+D and ef'=ef+D. Then, the equation E = e' + SIGMA(E-ef') - vxc' is equivalent to the one above.
c      
      subroutine quasiparticle(job1)
      use global
      use arrays
      use util, only: chr
      use wrapper
      use file
      use key
      use readwrite
      use, intrinsic :: iso_fortran_env
      Mpi( use Mwrapper )
      implicit none
      type(jobtype)                       :: job1 ! intent(in) omitted because routine set_stage temporarily modifies job1%full (because of call to read_dump)
      complex_dp,     allocatable         :: sigmac(:,:),dsigmac(:,:)
      complex_dp,     allocatable         :: evec(:,:),evec0(:,:),eval(:),evecc(:,:),evalc(:),sigf(:),hqp(:,:)
      complex_dp,     allocatable, target :: selft(:,:)
      MCOMPLEX_dp,    allocatable         :: sigmax(:,:),vxc(:,:),evecx(:,:),ihqp(:,:),eigv(:,:)
      real_dp,        allocatable         :: evalx(:),freq1(:),hks(:),olap(:),eig(:),freqt(:),overlap(:,:)
      integer,        allocatable         :: pnt(:),blockt(:,:)
      logical,        allocatable         :: found(:),fail(:)
      integer,        parameter           :: minstep = 10
      integer                             :: npole,np,mx(2)
      logical                             :: fitconst,fitconstr,allowup,newkpt
      complex_dp                          :: eqp,eqp0,eqp1,sigc,sigc0,dsigc,dsigc0,cdum
      real_dp                             :: vxc0,eks,sigx,rdum,w,fmin=0,fmax=0
      real_dp,        parameter           :: spec_warn=0.02,spec_pad=1.5
      real_dp,        allocatable         :: spec(:)
      real_dp                             :: wspec(4),alignvxc
      integer                             :: nspec
      integer                             :: iselfx0,iselfx,iselfc0,iselfc
      integer                             :: iblock,ib,ibandq,ikptq,ikptq0,ikindx,ispin,spin,s1,s2,n,ndeg
      integer                             :: istep,nfreqt,oselft
      integer                             :: i,j,k,ii,jb,imax
      integer                             :: iunit,iunits,iunitx
      integer                             :: itype,l
      logical                             :: ldum,lwrite,lfail,fullt
      character(11)                       :: cspin
      character(52),  parameter           :: hash = '####################################################'
      character(:),   allocatable         :: arg(:)
      character(20)                       :: parse_out(3)
      character(30)                       :: chr_band
      real                                :: cputime
      real_dp                             :: ch2r
      integer                             :: ch2i
      type padetype
        integer                           :: n
        complex_dp,   allocatable         :: a(:),b(:)
        complex_dp                        :: c
      end type padetype
      type(padetype), allocatable         :: pade(:,:),pade_t(:,:)
      type(statetype),allocatable         :: state(:)
      complex_dp,     allocatable         :: enec(:,:,:),enec_hlp(:,:,:)
      real_dp,        allocatable         :: wintgr1(:,:,:),ene1(:,:,:)
      type(padetype), allocatable         :: pade_w(:,:,:,:)
      logical                             :: l_inter
# ifdef MPI
      integer                             :: Merr,m
# endif

      call cpu_time(cputime)

      alignvxc = 0

      Rwrite(6,'(////A)') '### subroutine: quasiparticle ###'

# ifdef CHECK_SENERGY
      Warn('CHECK_SENERGY is enabled. Results below are meaningless.')
# endif      

      if(all(job1%type/=[J_HF,J_GW,J_SX,J_COSX,J_PBE0,J_GT,J_KS])) Bug('Wrong job number.')
      if(job1%type==J_KS.and.job1%full) Bug('FULL calculations not possible with JOB KS.')
      if(job1%type==J_PBE0) then
        beginSingle
        call read_pot(2)
        if(ogrid==2) then
          do ispin = 1,nspin
            do itype = 1,ntype
              do l = 1,nlh(itype)
                call interpolate_radial(vmt_xc(:,l,itype,ispin),rgrid0(:,itype),0d0,3,grid(itype))
              enddo
              n                        = grid(itype)%number
              vmt_xc(:n,1,itype,ispin) = vmt_xc(:n,1,itype,ispin) / rgrid(:n,itype)
            enddo
          enddo
        endif
        vmt_xc = vmt_xc / 4
        vpw_xc = vpw_xc / 4        
        endSingle
        MpiO( call Mcast(vmt_xc,comm=Ocomm) ; call Mcast(vpw_xc,comm=Ocomm) )
        Nfence(vmt) ; Nfence(vmt_xc) ; Nfence(vpw) ; Nfence(vpw_xc)
      endif

c     We redefine freq as freqr if oselfc>1 because the imag. frequencies (freq) are not needed anymore.
      if(any(job1%type==[J_GW,J_GT]).and.oselfc/=1) then
        nfreq = nfreqr ; deallocate(freq) ; allocate(freq(nfreq))
        freq  = freqr  ; deallocate(freqr)
      endif

c     Read parameters from spex.inp
      Rbegin
      if(any(job1%type==[J_GW,J_GT]).and.oselfc==1) then
        call getkey(inp,'CONTINUE', arg, section='SENERGY', status=i)
        if(i<=1) then
          npole = 0
          if(mod(nfreq,2)==0) then ; write(6,'(/A)') 'Number of mesh points even: No constraint will be used.'
          else                     ; write(6,'(/A)') 'Number of mesh points odd: Constraint Im dSIG/de = 0 will be used.'
          endif
        else
          fitconst  = .false.
          fitconstr = .false.
          allowup   = .false.
 1        i         = len_trim(arg(1))
          if     (arg(1)(i:i)=='+') then ; fitconst  = .true. ; arg(1)(i:i) = ' ' ; goto 1
          else if(arg(1)(i:i)=='c') then ; fitconstr = .true. ; arg(1)(i:i) = ' ' ; goto 1
          else if(arg(1)(i:i)=='*') then ; allowup   = .true. ; arg(1)(i:i) = ' ' ; goto 1
          endif
          npole = ch2i(arg(1),'CONTINUE')
          if(npole<=0) Error('Number of poles after CONTINUE must be positive.')
          i = 4*npole ; if(fitconst)  i = i + 2
          j = 2*nfreq ; if(fitconstr) j = j + 2
          if(i>j) Error('Fit is under-determined (more parameters than conditions).')
        endif
        if(allocated(arg)) deallocate ( arg )
      endif
      if(any(job1%type==[J_GW,J_GT])) then
        call getkey(inp,'ALIGNVXC', alignvxc, section='SENERGY', status=i, allow_eV=.true.)
        if(i==1) alignvxc = huge(0d0)
      endif
      call getkey(inp,'INTERPOL',l_inter, section='WANNIER', default=.false., status=i)
      if(l_inter) then
        if(.not.allocated(kpt_path).and.i==1) Warn('KPTPATH or file containing q points required for WANNIER/INTERPOL')
        if(nwan==0)                           Warn('Wannier definition required for WANNIER/INTERPOL')
        l_inter = (allocated(kpt_path).or.i==2) .and. nwan>0
      endif
      if(l_qsgw) write(6,'(A)') 'File qsgw exists and will be read.'

      Rend

c     Read sigt if required
      if(job1%label=='GWT') call read_sigt

      Mpi( call Mbroadcast_inp )

      if(oselfc> 1) allocate ( freq1(nfreq) )
      if(oselfc==4) freq1 = freq

c     Define SPECTRAL frequency range
      wspec = [ 0d0,0d0,0d0,-huge(0d0) ]
      nspec = 0
      Rif(any(job1%type==[J_GW,J_GT])) then
        call getkey(inp,'SPECTRAL',arg,section='SENERGY',status=i)
        if(i/=0) then
          if(oselfc==2) then
            Warn('No SPECTRAL plot possible with this CONTOUR definition.')
          else
            if(i==1) then
              if(oselfc==1) then
                wspec(:3) = [ huge(0d0),tiny(0d0),0.002d0 ]
                do i = 1,size(job1%band)
                  eks      = ene(job1%band(i),job1%kpt(i),job1%spin(i)) - efermi
                  wspec(1) = min(wspec(1),eks-spec_pad)
                  wspec(2) = max(wspec(2),eks+spec_pad)
                enddo
              else if(oselfc==3) then
                n = size(job1%band) ; allocate(spec(2*n))
                do i = 1,size(job1%band)
                  eks       = ene(job1%band(i),job1%kpt(i),job1%spin(i)) - efermi
                  spec(i)   = minval ( eks + sign(1d0,eks)*freq )
                  spec(i+n) = maxval ( eks + sign(1d0,eks)*freq )
                enddo
                wspec(1) = minval(spec(:n))
                wspec(2) = maxval(spec(n+1:))
                wspec(3) = (freq(2)-freq(1))/5
                if(any(spec(:n)>wspec(1)-spec_warn).or.any(spec(n+1:)<wspec(2)+spec_warn))
     &            Warn('SPECTRAL frequ. range goes beyond calculated range. You will get warnings...')
                deallocate(spec)
              else if(oselfc==4) then
                wspec(:3) = [ freq1(1),freq1(nfreq),(freq1(2)-freq1(1))/5 ]
              endif
            else if(i==2) then
              if(index(arg(1),'..')/=0) then
                call parser(parse_out,arg(1),'{1..2,3}','SPECTRAL')
                Info('".." obsolete; use ":" instead.')
              else
                call parser(parse_out,arg(1),'{1:2,3}','SPECTRAL')
              endif
              do i = 1,3 ; wspec(i) = ch2r(parse_out(i),'*SPECTRAL') ; enddo
              if(size(arg)==2) wspec(4) = ch2r(arg(2),'*SPECTRAL')
              if(wspec(1)>wspec(2).or.wspec(3)<0) Error('SPECTRAL range incorrect. {a..b,c}, a<=b, c>0.')
              deallocate(arg)
            endif
            if(njob>1) then
              iunits = fopen('spectral'//Chr(job1%indx),status='unknown')
              if(wrtext) iunitx = fopen('senergy'//Chr(job1%indx), status='unknown')
            else
              iunits = fopen('spectral',status='unknown')
              if(wrtext) iunitx = fopen('senergy', status='unknown')
            endif
            nspec = nint((wspec(2)-wspec(1))/wspec(3))+1
          endif
        endif
      endif
      Mpi( call Mcast(wspec); call Mcast(nspec) )

c     Prepare pade_t for GWT in case of combination "GW FULL" and "GT"(diag) (->pade_t)
      if(job1%label=='GWT'.and.oselft==1.and.job1%full.and..not.fullt) then
        n       = size(job1%band) ; allocate(pade_t(n,1),state(n))
        state%b = job1%band(block(:n,iblock))
        state%k = job1%kpt(block(:n,iblock))
        state%s = job1%spin(block(:n,iblock))
        call set_stage(2) ; if(job1%full) Bug('job1full /= fullt.')        
        call prepare_pade(pade_t,n,1,0,state,.false.)
        call set_stage(1)
        deallocate(state)
      endif

c     Prepare pade_w for Wannier interpolation
      if(l_inter) then
        if(wspec(3)/=0) then
          n = 1 ; if(job1%full) n = nwanband
          if(oselfc==1) then
            if(npole==0) then ; np = (nfreq+1)/2
            else              ; np = npole ; if(fitconst) Warn('SPECTRAL/CONTINUE-c ignored.')
            endif
            allocate ( pade_w(nwanband,n,nkpti,nspin1) )
            do ispin = 1,nspin1 ; do ikptq = 1,nkpti ; do i = 1,n ; do j = 1,nwanband
              allocate( pade_w(j,i,ikptq,ispin)%a(np) )
              allocate( pade_w(j,i,ikptq,ispin)%b(np) )
            enddo ; enddo ; enddo ; enddo
          endif
        endif
        ifR allocate ( enec( min(minval(job1%band),wanbandi) : max(maxval(job1%band),bando+1,wanbandf) ,nkpti,
     &                                     minval(job1%spin) : maxval(job1%spin) ) )
      else
        ifR allocate ( enec( minval(job1%band) : max(maxval(job1%band),bando+1) ,nkpti, minval(job1%spin) : maxval(job1%spin) ) )
      endif
      ifR enec = huge(0d0)

      beginSingle

c     Write QSGW self-energy
      if(job1%full) call write_qsgw

c     Open file to write quasiparticle states
      if(job1%full) then
        iunit = fopen('spex.qp',form='unformatted',status='unknown')
        write(iunit) job1%type,nblock
      endif

c     ALIGNVXC: Set vxc correction from VBM
      if(alignvxc==huge(0d0)) then
        rdum = maxval(ene(:maxeband,:nkpt,:),ene(:maxeband,:nkpt,:)<efermi+1d-8)
        i    = minloc( [ (abs(ene(job1%band(j),job1%kpt(j),job1%spin(j))-rdum),j=1,size(job1%band)) ] ,1)
        if(abs(ene(job1%band(i),job1%kpt(i),job1%spin(i))-rdum)>1d-6) then
          Warn('State closest below Fermi energy not defined in JOB.')
          write(0,'(A)')         '               State used instead: '//chr_band(job1%band(i),job1%kpt(i),job1%spin(i))
          write(0,'(A,F10.5,A)') '               with energy',ene(job1%band(i),job1%kpt(i),job1%spin(i))*hartree,' eV.'
        else
          write(6,'(/A)')        'State closest below Fermi energy: '//chr_band(job1%band(i),job1%kpt(i),job1%spin(i))
          write(6,'(A,F10.5,A)') 'with energy',ene(job1%band(i),job1%kpt(i),job1%spin(i))*hartree,' eV.'
        endif
        if(job1%type==J_GT) then
          vxc0 = 0
          sigx = 0
        else
          allocate(vxc(1,1))
          if(ovxc==0) then ; call read_vxc(vxc,[job1%band(i)],1,kptp(job1%kpt(i)),job1%spin(i),l_qsgw)
          else             ; call calc_vxc(vxc,[job1%band(i)],1,kptp(job1%kpt(i)),job1%spin(i),l_qsgw)
          endif
          vxc0 = vxc(1,1)
          sigx = get_self(i,1)
          deallocate(vxc)
        endif
        eks      = ene(job1%band(i),job1%kpt(i),job1%spin(i)) - efermi ; if(oselfc==2.or.oselfc==3) freq1 = eks + sign(1d0,eks)*freq
        sigc     = get_self(i,2)
        alignvxc = sigx + sigc - vxc0
        write(6,'(A)')         'Chemical potential (xc potential) will be aligned to that of the self-energy.'
        write(6,'(A,F10.5,A)') 'Quasiparticle correction (Z=1)',alignvxc*hartree,' eV is set to zero.'
      endif

      endSingle

      Mpi( call Mcast(alignvxc) )

      if(alignvxc/=0) then
        ene    = ene    + alignvxc
        efermi = efermi + alignvxc
        Rwrite(6,'(/A,F10.5,A)') 'Fermi energy shifted to ',efermi*hartree,' eV according to ALIGNVXC.'
      endif

c
c     Write list of k points
      Rbegin
      write(6,'(/A)') 'List of k points'
      allocate(found(nkpt2)) ; found = .false.
      do iblock = 1,nblock
        ib     = block(1,iblock)
        ikptq  = job1%kpt(ib)
        ikptq0 = get_ikptq0(ikptq)
        if(.not.found(ikptq)) then
          found(ikptq) = .true.
          write(6,'(I5,3F10.5'NoA) abs(ikptq0),kpt(:,ikptq) ; if(ikptq0<=0) write(6,'(A'NoA) '  (additional)'
          write(6,*)
        endif
      enddo
      deallocate(found)
      Rend

      iselfx0 = 0
      iselfc0 = 0

c
c     Loop over blocks

      do iblock = 1,nblock

        ib     = block(1,iblock)
        ibandq = job1%band(ib)
        ikptq  = job1%kpt(ib) ; ikindx = kindx(ikptq)
        ispin  = job1%spin(ib)
        n      = sizeblock(iblock)
        ndeg   = deg(ibandq,ikptq,ispin) ; if(ndeg<ibandq) ndeg = deg(ndeg,ikptq,ispin)
        ndeg   = ndeg - ibandq + 1 ! dimension of current irrep

        allocate ( state(n) )
        state%b = job1%band(block(:n,iblock))
        state%k = job1%kpt(block(:n,iblock))
        state%s = job1%spin(block(:n,iblock))

        if(nspin==2.and.ispin==1) cspin = ', spin up'
        if(nspin==2.and.ispin==2) cspin = ', spin down'
        if(nspin==1)              cspin = ' '

        lwrite = .false.
        newkpt = ib==1.or.job1%kpt(max(1,ib-1))/=ikptq.or.job1%spin(max(1,ib-1))/=ispin

        ikptq0 = abs(get_ikptq0(ikptq))
        if(job1%full) then
          Rwrite(6,'(///A)') hash(:40+len_trim(cspin))
          Rwrite(6,'(2(A,I5),A,I2,A)') '### BLOCK:',iblock,'  (kpt:',ikptq0,', dim:',ndeg,trim(cspin)//') ###'
          Rwrite(6,'(A)') hash(:40+len_trim(cspin))
        else if(newkpt) then
          lwrite = .true.
          Rwrite(6,'(///A)') hash(:21+len_trim(cspin))
          Rwrite(6,'(A,I5,A)')    '### K POINT:',ikptq0,trim(cspin)//' ###'
          Rwrite(6,'(A)') hash(:21+len_trim(cspin))
        endif

        if(wspec(3)/=0.and.newkpt) then
          if(allocated(spec)) then ; Mpi( call Msum(spec) ; ) Rcall write_spectral(job1%kpt(ib-1),job1%spin(ib-1))
          else                     ; allocate(spec(nspec))
          endif
          spec = 0
        endif

        if(use_sym.and.job1%full.and.job1%type==J_GW) then
          Nfence(selfc)
          Obegin
          rdum = 0
          do i = 1,size(selfc,2)
            call symmetrize_matrix(selfc(iselfc0+1:iselfc0+n**2,i),n,iblock,rdum)
          enddo
          Rwrite(6,'(/A,ES8.1)') 'Symmetry fulfilled to within',rdum/size(selfc)
          Oend
          Nfence(selfc)
        endif

        ! Analytic continuation
        if(any(job1%type==[J_GW,J_GT]).and.oselfc==1) then
          Rif(job1%full) write(6,'(//A/)') '--- ANALYTIC CONTINUATION (diagonal elements) ---'
          allocate(pade(n,n))
          call prepare_pade(pade,n,n,iselfc0,state,.true.)
          if(job1%label=='GWT'.and.(job1%full.eqv.fullt)) then
            allocate(pade_t(n,n))
            call set_stage(2)
            call prepare_pade(pade_t,n,n,iselfc0,state,.false.)
            call set_stage(1)
          endif
          if(allocated(pade_w).and.ikptq<=nkpti) then
            if(job1%label=='GWT') Error('Not implemented: Wannier interpolation for GWT.')
            do i = 1,n   ; ib = job1%band(block(i,iblock)) - wanbandi + 1
              do j = 1,n ; jb = job1%band(block(j,iblock)) - wanbandi + 1
                if(min(ib,jb)>=1.and.max(ib,jb)<=nwanband) then
                  if(.not.job1%full) jb = 1
                  pade_w(ib,jb,ikptq,ispin) = pade(i,j)
                endif
              enddo
            enddo
          endif
        endif

c
c       DIAGONAL ELEMENTS
c
        if(n>ndeg.or.lwrite) then ; Rwrite(6,'(//A/)') '--- DIAGONAL ELEMENTS [eV] ---'
        else if(job1%full)   then ; Rwrite(6,'(//A/)') '--- DIAGONAL ELEMENTS [eV] (identical to full solution) ---'
        endif

        Rif(job1%full.or.lwrite) then
          select case(job1%type)
            case(J_HF)   ; write(6,'(A)')   ' Bd       vxc    sigmax        KS        HF'
            case(J_PBE0) ; write(6,'(A)')   ' Bd       vxc    sigmax        KS      PBE0'
            case(J_SX)   ; write(6,'(A)')   ' Bd       vxc    sigmax        KS        SX'
            case(J_COSX) ; write(6,'(A)')   ' Bd       vxc    sigmax        KS      COSX'
            case(J_GT)   ; write(6,'(A)')   ' Bd    sigmac         Z        KS        GT'
            case(J_KS)   ; write(6,'(A)')   ' Bd       vxc        KS'
            case default ; write(6,'(A'NoA) ' Bd       vxc    sigmax    sigmac         Z        KS        HF'
# ifndef additional
            write(6,'(A'NoA) '    '
# endif
            if(job1%label=='GW') write(6,'(A'NoA) ' '
            if(oselfc==1.or.nfreq>2) then ; write(6,'(A'NoA) '     '//trim(job1%label)//' lin/dir'
            else                          ; write(6,'(A'NoA) '         '//trim(job1%label)//' lin'
            endif
# ifdef additional
            write(6,'(A'NoA) '/Z=1       '//trim(job1%label)//' correction'
# endif
            write(6,*)
          end select
        endif

        beginSingle
        do i = 1,n,ndeg
          ib     = block(i,iblock)
          ibandq = job1%band(ib)

          iselfx = iselfx0 + i*(i+1)/2
          iselfc = iselfc0 + (i-1)*n+i

          if(job1%type==J_GT) then
            vxc0 = alignvxc ! vxc assumed not to contain GT, so no subtraction of vxc necessary; deltaex considered correction of vxc, so, again, no correction
            sigx = 0
          else
            allocate(vxc(1,1))
            if(ovxc==0) then ; call read_vxc(vxc,[ibandq],1,kptp(ikptq),ispin,l_qsgw)
            else             ; call calc_vxc(vxc,[ibandq],1,kptp(ikptq),ispin,l_qsgw)
            endif
            vxc0 = vxc(1,1) + alignvxc + deltaex * (2*ispin-3)/2 ! deltaex considered correction of vxc; since we subtract vxc, we also have to correct for deltaex
            deallocate(vxc)
            if(job1%type/=J_KS)   sigx = selfx(iselfx)
            if(job1%type==J_PBE0) sigx = sigx / 4
          endif
          eks = ene(ibandq,ikptq,ispin) - efermi


          if(all(job1%type/=[J_GW,J_GT])) goto 100

          if(oselfc==3) freq1 = eks + sign(1d0,eks)*freq

          ! SPECTRAL function
          Rif(wspec(3)/=0) then
            ldum = .not.job1%full.or.n<=ndeg
            if(wrtext) call write_spectral_head(iunitx,ikptq,ispin,'Self-energy',.true.)
            if(wrtext.or.ldum) then
              if(oselfc>=3) then
                fmin = minval(freq1) - 1d-10
                fmax = maxval(freq1) + 1d-10
                if(wspec(1)<fmin-spec_warn) 
     &            Warn('SPECTRAL frequ. below calculated range by more than '//Chr(spec_warn)//'htr: '//chr(fmin-wspec(1)))
                if(wspec(2)>fmax+spec_warn)
     &            Warn('SPECTRAL frequ. above calculated range by more than '//Chr(spec_warn)//'htr: '//chr(wspec(2)-fmax))
              endif
              w = wspec(1)
              do j = 1,nspec
                call selfenergy_interpolate(sigc,dsigc,w*(1d0,0d0),1,i,i,iblock)
                if(wrtext) then
                  if(oselfc>=3.and.(w<fmin.or.w>fmax)) then ; write(iunitx,'(F12.7,2X,2F20.10,A)') w*escale,sigx+sigc,' *'
                  else                                      ; write(iunitx,'(F12.7,2X,2F20.10)')   w*escale,sigx+sigc
                  endif
                endif
                if(ldum) then
                  sigc    = real(sigc) + img * max( imag(sigc)*sign(1d0,-w) , wspec(4) )
                  cdum    = w - eks - sigx - sigc + vxc0
                  spec(j) = spec(j) + imag(cdum**(-1)) / pi * min(n,ndeg)
                endif
                w = w + wspec(3)
              enddo
            endif
            if(wrtext) write(iunitx,'(/)')
          endif

c          if(metal) write(6,'(A,F8.4)') 'Bare Drude pole:',eks+sign(1d0,eks)*plasma

          if(oselfc==2) then
            ! (a) Linearized solution
            sigc   = ( selfc(iselfc,2) + selfc(iselfc,1) ) / 2
            dsigc  = ( selfc(iselfc,2) - selfc(iselfc,1) ) / (freq(2)-freq(1)) * sign(1d0,eks)
            eqp    = eks + (1-dsigc)**(-1) * ( sigx + sigc - vxc0 )
            sigc0  = sigc ; dsigc0 = dsigc ; eqp0 = eqp ! values at KS energy
          else
            ! (b) Nonlinear solution
            ! solve quasiparticle equation (diagonal elements) with Newton method
            eqp   = eks ! start value
            istep = 0
            lfail = .false.
            do
              eqp1 = eqp
              if(oselfc/=1) eqp = real(eqp)
              call selfenergy_interpolate(sigc,dsigc,eqp,1,i,i,iblock)
              istep = istep + 1
              eqp   = (1-dsigc)**(-1) * ( eks + sigx + sigc - vxc0 - eqp*dsigc )
              eqp   = eqp1 + (eqp-eqp1) / (1+istep/500)
              if(istep==1) then ; sigc0 = sigc ; dsigc0 = dsigc ; eqp0 = eqp ; endif ! values at KS energy
              if(abs(eqp-eqp1)<1d-10) exit
              if(istep>5000) then
                lfail = .true.
                exit
              endif
c                if(oselfc/=1) write(6,'(3F15.10)') (freq1(j),selfc(iselfc,j),j=1,nfreq)
c                write(6,'(A,ES8.1)')  'Last residual error:',abs(eqp-eqp1)
c                write(6,'(A,2F10.5)') 'Last energy:',eqp
c                Error('Solution of quasiparticle equation did not converge after 5000 iterations.')
c              endif
            enddo
          endif

          ! Output
 100      continue
          do j = i,min(n,i+ndeg-1)
            ib     = block(j,iblock)
            ibandq = job1%band(ib)
            ldum   = oselfc==1.or.nfreq>2
            if(job1%type==J_GT) then
              write(6,'(I3,5F10.5)')     ibandq, hartree*real(sigc0),                 real((1-dsigc0)**(-1)),
     &                                           hartree*(eks+efermi),                hartree*(eqp0+efermi)
              write(6,'(3X,2F10.5'NoA)           hartree*imag(sigc0),                 imag((1-dsigc0)**(-1))
              if(ldum) write(6,'(10X,2F10.5'NoA) hartree*(eqp+efermi)
              if(ldum.and.lfail) write(6,'(A'NoA) ' ?'
              write(6,*)
              if(ldum) then ; cdum = eqp  + efermi
              else          ; cdum = eqp0 + efermi
              endif
            else if(job1%type==J_KS) then
              write(6,'(I3,2F10.5)')    ibandq,  hartree*vxc0,                        hartree*(eks+efermi)
              cdum = eks+efermi
            else if(job1%type/=J_GW) then
              write(6,'(I3,7F10.5)')    ibandq,  hartree*vxc0,                        hartree*sigx,
     &                                           hartree*(eks+efermi),                hartree*(eks+efermi+sigx-vxc0)
              cdum = eks+efermi+sigx-vxc0
            else
              write(6,'(I3,8F10.5'NoA) ibandq,   hartree*vxc0,                        hartree*sigx,
     &                                           hartree*real(sigc0),                 real((1-dsigc0)**(-1)),
     &                                           hartree*(eks+efermi),                hartree*(eks+efermi+sigx-vxc0),
     &                                           hartree*(eqp0+efermi)
# ifdef additional
              write(6,'(2F10.5'NoA)              hartree*(eqp0-eks)
# endif
              if(oselfc>=3) then
                rdum = real(eqp0)
                if(rdum<minval(freq1).or.rdum>maxval(freq1)) write(6,'(A'NoA) ' *'
              endif
              write(6,*)
              write(6,'(3X,20X,2F10.5'NoA)       hartree*imag(sigc0),                 imag((1-dsigc0)**(-1))
# ifdef additional
              if(ldum) write(6,'(20X,4F10.5)')   hartree*(eqp+efermi),                hartree*(eqp-eks)
              write(6,'(63X,4F10.5'NoA)          hartree*(eks+efermi+sigx+sigc0-vxc0),hartree*(sigx+sigc0-vxc0)
# else
              if(ldum) write(6,'(20X,2F10.5'NoA) hartree*(eqp+efermi)
# endif
              if(ldum.and.lfail) write(6,'(A'NoA) ' ?'
              if(oselfc==4) then
                rdum = real(eqp)
                if(rdum<minval(freq1).or.rdum>maxval(freq1)) write(6,'(A'NoA) ' *'
              endif
              write(6,*)
              if(ldum) then ; cdum = eqp  + efermi
              else          ; cdum = eqp0 + efermi
              endif
            endif

            if((.not.job1%full.or.n<=ndeg).and.ikptq<=nkpt) enec(ibandq,kptp(ikptq),ispin) = cdum
          enddo

        enddo

        endSingle

c
c       FULL MATRIX SOLUTION
c
        if(job1%full.and.n>ndeg) then
          Rwrite(6,'(//A)') '--- FULL SOLUTION [eV] ---'

          ! define HKS , VXC , SIGMAX
          allocate ( sigmax(n,n),vxc(n,n),evecx(n,n),evalx(n),hks(n),evec(n,n) )
          iselfx = iselfx0
          do j = 1,n   ; jb = job1%band(block(j,iblock))
            do i = 1,j ; ib = job1%band(block(i,iblock))
              iselfx      = iselfx + 1
              sigmax(i,j) = selfx(iselfx)
              sigmax(j,i) = MCONJG(sigmax(i,j))
            enddo
            hks(j) = ene(jb,ikptq,ispin) - efermi
          enddo
          beginSingle
          if(ovxc==0) then ; call read_vxc(vxc,job1%band(block(:n,iblock)),n,kptp(ikptq),ispin,l_qsgw)
          else             ; call calc_vxc(vxc,job1%band(block(:n,iblock)),n,kptp(ikptq),ispin,l_qsgw)
          endif
          rdum = alignvxc + deltaex * (2*ispin-3)/2
          evec = vxc    ; if(use_sym) call symmetrize_matrix(evec,n,iblock) ; vxc    = evec + rdum*identity(n)
          evec = sigmax ; if(use_sym) call symmetrize_matrix(evec,n,iblock) ; sigmax = evec
          if(job1%type==J_PBE0) sigmax = sigmax / 4
          endSingle
          Mpi( call Mcast(vxc) ; call Mcast(sigmax) )

          ! HF: diagonalize HF Hamiltonian     HamHF evecx = evalx evecx
          !                           with     HamHF = HKS - VXC + SIGMAX
          Rcall diagonalize(evecx,evalx,diagonal(hks)+sigmax-vxc)
          if(job1%type==J_GW) then

            ! GW: solve quasiparticle equation   HamQP(eval) evec = eval evec
            !                             with   HamQP(eval) = HKS - VXC + SIGMAX + SIGMAC(eval)
            allocate ( sigmac(n,n),dsigmac(n,n),eval(n),evec0(n,n) )
            allocate ( evecc(n,n),evalc(n),sigf(n),olap(n),hqp(n,n),pnt(n),fail(n) )
            Mpi( evecc = 0 ; evalc = 0 ; sigf  = 0 ; olap  = 0 ; fail = .false. )

            !     SPECTRAL function
            if(wspec(3)/=0) then
              allocate ( ihqp(n,n),eig(n),eigv(n,n) )
              do j = Mrange1(nspec)
                w       = wspec(1) + (j-1) * wspec(3)               ; call qp_hamiltonian(hqp,sigmac,dsigmac,w*(1d0,0d0))
                ihqp    = ( hqp - conjg(transpose(hqp)) ) / (2*img) ; call diagonalize(eigv,eig,ihqp)
                hqp     = ( hqp + conjg(transpose(hqp)) ) / 2
                eig     = max( eig*sign(1d0,-w) , wspec(4) )
                ihqp    = matmul(eigv,matmul(diagonal(eig),transpose(MCONJG(eigv))))
                hqp     = hqp + img * ihqp
                hqp     = w*identity(n) - hqp                       ; call inverse(hqp)
                spec(j) = spec(j) + sum ( [ (imag(hqp(ii,ii)),ii=1,n) ] ) / pi
              enddo
              deallocate ( ihqp,eig,eigv )
            endif

            !     solve nonlinear equation starting from KS energies, then adiabatically switch on offdiagonal elements
            lfail = .false.
            i     = 1
            Mpi(m =-1)
            do while(i<=n)
              ifMODP(m)
              istep = 0
              eqp   = hks(i)
              imax  = i
              evec  = identity(n)
              do
                evec0 = conjg(transpose(invert(evec)))
                eqp1  = eqp
                if(oselfc/=1) eqp = real(eqp)
                call qp_hamiltonian(hqp,sigmac,dsigmac,eqp)
                call qp_adiabatic(hqp,i,istep)
                call diagonalize_gen(evec,eval,hqp)
                call rorderp(pnt,real(eval),n)
                eval  = eval(pnt)
                evec  = evec(:,pnt)
                call maxolap(imax)                                   ! largest overlap with previous eigenvectors
                dsigc = dot_product(evec(:,imax),matmul(dsigmac,evec(:,imax)))
                eqp   = (1-dsigc)**(-1) * ( eval(imax) - dsigc*eqp ) ! extrapolated energy
                istep = istep + 1
                if(istep>1000) then
                  lfail = .true.
                  exit
                else
                  eqp = eqp1 + (eqp-eqp1) / (1+istep/50)
                  if(abs(eqp-eqp1)<1d-6.and.istep>minstep) exit
                endif
              enddo
              ii   = (imax-1) / ndeg * ndeg
              rdum = 0
              do j = ii+1,ii+ndeg
                rdum       = rdum + abs(eval(j)-eval(imax))
                fail(i)    = lfail
                evalc(i)   = eqp
                evecc(:,i) = evec(:,j)
                olap(i)    = sqrt(sum(abs(evec(ii+1:ii+ndeg,j))**2))
                sigf(i)    = hartree * dot_product(evecc(:,i),matmul(sigmac,evecc(:,i)))
                i          = i + 1
              enddo
              if(rdum>1d-8) then
                Warn('Bad degeneracy in subspace of eigenvalue '//Chr(real(eval(imax)))//': '//Chf(rdum,'F15.10'))
                if(rdum>1d-4) then
                  write(0,'(A)') '     123456789012345     123456789012345'
                  write(0,'(3F20.15)') (abs(eval(min(n,ii+1))-eval(ii)),eval(ii),ii=1,n)
                  Error('Symmetry break.')
                endif
              endif
              if(ndeg>1) then ! orthonormalize degenerate states
                evec = evecc
                do ib = 1,ndeg-1
                  do jb = 1,ib-1
                    evecc(:,ii+ib) = evecc(:,ii+ib) - evecc(:,ii+jb) * dot_product(evecc(:,ii+jb),evec(:,ii+ib))
                  enddo
                  evecc(:,ii+ib) = evecc(:,ii+ib) / sqrt(sum(abs(evecc(:,ii+ib))**2))
                enddo
              endif
              elseMOD( i = i + ndeg )
            enddo
            Mpi( call Msum(evecc) ; call Msum(evalc) ; call Msum(sigf) ; call Msum(olap) )
            MpiR( call mpi_reduce(mpi_in_place,fail,n,mpi_logical,mpi_lor,0,Mcomm,Merr) )
            MnoR( call mpi_reduce(fail,       lfail,n,mpi_logical,mpi_lor,0,Mcomm,Merr) )
            deallocate ( sigmac,dsigmac,evec0,eval,hqp,pnt )
          endif

          Rbegin

          ! assign HF(x)->KS states
          allocate ( overlap(n/ndeg,n/ndeg),pnt(n/ndeg) )
          do i = 1,n/ndeg
            do j = 1,n/ndeg
              overlap(j,i) = sum(abs(evecx((j-1)*ndeg+1:j*ndeg,(i-1)*ndeg+1:i*ndeg))**2) / ndeg
            enddo
          enddo
          do i = 1,n/ndeg
            mx               = maxloc(overlap) ; if(overlap(mx(1),mx(2))<0) Bug('Negative overlap.')
            pnt(mx(1))       = (mx(2)-1)*ndeg + 1
            overlap(mx(1),:) = -1
            overlap(:,mx(2)) = -1
          enddo
          deallocate(overlap)

          ! Output
          write(iunit) iblock,ikptq,ispin,ndeg,n,job1%band(block(:n,iblock))
          if(job1%type==J_GW) then
            write(iunit) evalc
            write(iunit) evecc ; evec = evecc
          else
            write(iunit) evalx
            write(iunit) evecx ; evec = evecx
          endif

          select case(job1%type)
            case(J_HF)   ; write(6,'(/A)') ' Bd        KS      olap        HF'
            case(J_PBE0) ; write(6,'(/A)') ' Bd        KS      olap      PBE0'
            case(J_SX)   ; write(6,'(/A)') ' Bd        KS      olap        SX'
            case(J_COSX) ; write(6,'(/A)') ' Bd        KS      olap      COSX'
            case(J_GW)   ; write(6,'(/A)') ' Bd        KS      olap        HF      olap        GW'
          end select

          do i = 1,n,ndeg ; k = pnt((i-1)/ndeg+1)
            do ii = 0,ndeg-1
              rdum = sqrt( sum ( abs(evecx(i:i+ndeg-1,k+ii))**2 ) )
              if(job1%type/=J_GW) then
                write(6,'(I3,3F10.5)')   job1%band(block(i+ii,iblock)), hartree*(  hks(i+ii)+efermi),
     &                                                            rdum, hartree*(evalx(k+ii)+efermi)
                cdum = evalx(k+ii) + efermi
              else
                write(6,'(I3,6F10.5'NoA) job1%band(block(i+ii,iblock)), hartree*(  hks(i+ii)+efermi),
     &                                                            rdum, hartree*(evalx(k+ii)+efermi),
     &                                                      olap(i+ii), hartree*(evalc(i+ii)+efermi)
                if(fail(i+ii)) write(6,'(A'NoA) ' ?'
                if(oselfc==4) then
                  rdum = real(evalc(i+ii))
                  if(rdum<minval(freq1).or.rdum>maxval(freq1)) write(6,'(A'NoA) ' *'
                endif
                write(6,*)
                cdum = evalc(i+ii) + efermi
              endif
              if(ikptq<=nkpt) enec(job1%band(block(i+ii,iblock)),kptp(ikptq),ispin) = cdum
            enddo
          enddo
          deallocate ( pnt )
          if(bandinfo) call write_bandinfo

          Rend

          if(job1%type==J_GW) deallocate ( fail )

c          ! Write natural orbitals and corresponding occupation numbers to fort.999 (spex.nat)
c          if(job1%type==2) then
c            if(any(job1%band(block(:n,iblock))<=bando)) then
c              allocate ( densmat(n,n) )
c              densmat = 0
c              do ib = 1,n
c                ibandq = job1%band(block(ib,iblock))
c                if(ibandq<=bando) then
c                  do i = 1,n ; do j = 1,n
c                    densmat(i,j) = densmat(i,j) + wintgr(ikptq,ibandq,ispin)*nkpt * evec(i,ib) * conjg(evec(j,ib))
c                  enddo ; enddo
c                endif
c              enddo
c              if(any(densmat/=0)) then
c                call diagonalize ( evecx,evalx,densmat )
c                do i = n,1,-1
c                  write(999,'(2I3,F15.10'NoA) ikptq,i,evalx(i)
c                  write(999,'(    F15.10'NoA) evecx(:,i)
c                  write(999,*)
c                enddo
c              endif
c              deallocate ( densmat )
c            endif
c          endif
          deallocate ( sigmax,evalx,evecx,hks,vxc,evec )
          if(job1%type==J_GW) deallocate ( evalc,evecc,sigf,olap )
        endif

        iselfx0 = iselfx0 + n*(n+1)/2
        iselfc0 = iselfc0 + n**2

        if(allocated(pade))                             deallocate ( pade )
        if(allocated(pade_t).and.(job1%full.eqv.fullt)) deallocate ( pade_t )
        deallocate ( state )

      enddo

      if(allocated(pade_t).and.(job1%full.neqv.fullt)) deallocate ( pade_t )

      Rif(job1%type==J_PBE0) call read_pot(1)

c     Close spex.qp
      Rif(job1%full.and.n>ndeg) call fclose(iunit)

c     Write last spectral function and close "spectral"
      if(allocated(spec)) then
        Mpi( call Msum(spec) )
        Rcall write_spectral(ikptq,ispin)
        Rcall fclose(iunits) ; Rif(wrtext) call fclose(iunitx)
        deallocate(spec)
      endif

      rdum = cputime
      call cpu_time(cputime)
      Rwrite(6,'(/A,F12.5)') 'Timing (quasiparticle equation):',cputime-rdum

      ! Determine spin:
      ! spin = 0 : Both spins have complete IBZ
      ! spin = s : Only spin s has complete IBZ
      ! spin = 3 : Both spins have incomplete IBZ
      spin = 3
      do ispin = minval(job1%spin),maxval(job1%spin)
        spin = spin - (3-ispin)
        do k = 1,nkpti
          if( all( pack(kptp(job1%kpt),job1%spin==ispin) /= k ) ) then ; spin = spin + (3-ispin) ; exit ; endif
        enddo
      enddo
      if(spin==0) then ; s1 = 1    ; s2 = nspin1
      else             ; s1 = spin ; s2 = spin
      endif

c     Wannier interpolation
      if(l_inter.and.spin<3) then
        Rbegin
        write(6,*)
        allocate(enec_hlp,mold=enec)
        enec_hlp = enec
        call complete_ene(enec_hlp(wanbandi:wanbandf,:,s1:s2),wanbandi,wanbandf,spin,2,i,j,k,0,.true.)
        if     (k>0) then; Bug('Energies at '//Chr(k)//' kpoints missing.')
        else if(j>0) then; Warn('Energies of '//Chr(j)//' bands missing, now adjusted to KS. Interpolation might be inaccurate.')
        else if(i>0) then; Warn('Energies of '//Chr(i)//' states missing, now adjusted to KS. Check interpolation.')
        endif
        write(6,*)
        call wannier_interpolation('1',enec_hlp(wanbandi:wanbandf,:,s1:s2),spin)
        deallocate(enec_hlp)
        Rend
        if(wspec(3)/=0) then
          call spectral_wannier_interpolation(spin)
          if(oselfc==1) deallocate(pade_w)
        endif
      else if(l_inter) then
        Warn('Energies not complete in (irreducible) Brillouin zone. No Wannier interpolation.')
      endif

c     New Fermi energy
      Rbegin
      if(spin==0.or.nspin1==1.and.spin==1) then
        if(nelec-lbound(enec,1)+1<=0) then
          Warn('Occupied states missing. Cannot calculate new Fermi energy.')
        else
          write(6,*)
          k = 0
          call complete_ene(enec,lbound(enec,1),ubound(enec,1),spin,2,i,j,k,1,.true.)
          if     (k>0) then; Bug('Energies at '//Chr(k)//' kpoints missing.')
          else if(j>0) then; Warn('Energies of '//Chr(j)//' bands missing, now adjusted to KS. New Fermi energy may be inaccurate.')
          else if(i>0) then; Warn('Energies of '//Chr(i)//' states missing, now adjusted to KS. Check new Fermi energy.')
          endif
          allocate(wintgr1(nkpt,maxeband,nspin1),ene1(lbound(enec,1):ubound(enec,1),nkpt,nspin1))
          ene1                                        = ene(lbound(enec,1):ubound(enec,1),:nkpt,:)
          ene(lbound(enec,1):ubound(enec,1),:nkpti,:) = enec(:,:nkpti,:)
          do k = nkpti+1,nkpt
            ene(lbound(enec,1):ubound(enec,1),k,:) = enec(:,kptp(k),:)
          enddo
          call fermi(rdum,cdum,ldum,wintgr1,.true.)
          ene(lbound(enec,1):ubound(enec,1),:nkpt,:) = ene1
          deallocate(wintgr1,ene1)
        endif
      endif
      deallocate ( enec )
      Rend

      if(oselfc>1) deallocate ( freq1 )

      if(alignvxc/=0) then
        ene    = ene    - alignvxc
        efermi = efermi - alignvxc
      endif

      contains

c     ------------

c     Interpolates selfenergy at e.
c
c     iblock present : Get iselfc from index_selfc according to iblock,i,j (normal case)
c     else           : Set iselfc:=i
c
c     oselfc=1 : p=0 - use Padé approximant
c                p=1 - use precalculated pade(i,j)
c                p=2 - use precalculated pade_t(i,j)
      recursive subroutine selfenergy_interpolate(self,dself,e,p,i,j0,iblock)
      use global
      implicit none
      integer,    intent(in)           :: i,p
      integer,    intent(in), optional :: j0,iblock
      complex_dp, intent(in)           :: e
      complex_dp, intent(out)          :: self,dself
      complex_dp                       :: c,dc
      integer                          :: j,np,iselfc
      type(padetype)                   :: pd0(1,1),pd
      if(present(iblock)) then ; j = j0 ; iselfc = index_selfc(iblock,i,j)
      else                     ; j = 1  ; iselfc = i
      endif
      if(oselfc==1) then
        if     (p==0) then ; call prepare_pade(pd0,1,1,iselfc-1)
                             np = pd0   (1,1)%n ; allocate(pd%a(np),pd%b(np)) ; pd = pd0   (1,1)
        else if(p==1) then ; np = pade  (i,j)%n ; allocate(pd%a(np),pd%b(np)) ; pd = pade  (i,j)
        else if(p==2) then ; np = pade_t(i,j)%n ; allocate(pd%a(np),pd%b(np)) ; pd = pade_t(i,j)
        else               ; Bug('Wrong p argument.')
        endif
        if(real(e)<0) then
          pd%b = conjg(pd%b)
          pd%a = conjg(pd%a)
          dc   = selfc(iselfc,1) ; if(dc/=0) pd%a = pd%a * dc/conjg(dc)
        endif
        self  = pd%c + sum ( pd%a / ( e - pd%b )    )
        dself =      - sum ( pd%a / ( e - pd%b )**2 )
      else
        call interpolate(c,dc,real(e),selfc(iselfc,:nfreq),freq1,nfreq)
        self  = c + dc*imag(e)
        dself = dc
      endif
      if(job1%label=='GWT'.and.p/=2) then
        if(.not.present(iblock)) Bug('Need iblock but is not present.')
        call add_sigt(self,dself,e,i,j,iblock)
      endif
      end subroutine selfenergy_interpolate

c     ------------

c
c     Read spex.sigt and define
c     oselft           - GT self-energy type
c     freqt(:nfreqt)   - frequency mesh for GT self-energy
c     selft(:,:nfreqt) - GT self-energy
c     blockt(:,:)      - GT self-energy blocks
c     fullt            - true (GT FULL), false (GT)      
      subroutine read_sigt
      use util, only: same_eigenspace
      implicit none
      logical :: full,lkpt(nkpt),done(size(job1%band)),lsym
      integer :: i,j,k,jb,kb,iband,ispin,ikpt
      ! store current stage in (pnt,oself,frq,nfrq,full,blk)
      full = job1%full
      call set_stage(0)
      Rbegin
      ! read GT stage (selft,oselft,freqt,nfreqt,fullt) from file
      lkpt = .false.
      deallocate(freq)
      call read_dump(0,lkpt,nkpt,'spex.sigt',job1)
      fullt = job1%full
      if(full.eqv.fullt) then ; allocate(blockt(size(block,1),nblock)) ; blockt      = block
      else                    ; allocate(blockt(1,size(job1%band)))    ; blockt(1,:) = (/(i,i=1,size(job1%band))/)
      endif
      if(oselfc>1) then
        nfreq = nfreqr ; allocate(freq(nfreq))
        freq  = freqr  ; deallocate(freqr)
      endif
      allocate(selft(size(job1%band),nfreq))
      oselft =  oselfc
      selfc  => selft
      nfreqt =  nfreq ; allocate(freqt(nfreqt))
      freqt  =  freq  ; deallocate(block) ; allocate(block(size(blockt,1),size(blockt,2)))
      block  =  blockt
      call read_dump(2,lkpt,nkpt,'spex.sigt',job1)
      if(all(.not.lkpt)) Error('JOB GWT requires the GT self-energy stored in spex.sigt.')
      if(any(.not.lkpt)) then
        if(all(lkpt(:nkpti)).and.all(.not.(lkpt(nkpti+1:)))) then
          Info('GT data obviously from IRREP calculation.')
        else
          Error('GT calculation not complete according to spex.sigt.')
        endif
      endif
      ! Symmetrize selft (to avoid warnings about bad degeneracies)
      if(full) then
        if(fullt) Bug('GT FULL not implemented.')        
        lsym = .false.
        rdum = 0
        do i = 1,size(selft,2)
          done = .false.
          do j = 1,size(selft,1)
            if(done(j)) cycle
            jb    = block(1,j)
            iband = job1%band(jb)
            ikpt  = job1%kpt(jb)
            ispin = job1%spin(jb)
            cdum  = selft(j,i)
            n     = 1
            do k = j+1,size(selft,1)
              kb = block(1,k)
              if(all([job1%kpt(kb),job1%spin(kb)]==[ikpt,ispin]).and.same_eigenspace(job1%band(kb),iband,ikpt,ispin)) then
                cdum = cdum + selft(k,i)
                n    = n + 1
              endif
            enddo
            if(n>1) then
              lsym = .true.
              cdum = cdum / n
              do k = j,size(selft,1)
                kb = block(1,k)
                if(all([job1%kpt(kb),job1%spin(kb)]==[ikpt,ispin]).and.same_eigenspace(job1%band(kb),iband,ikpt,ispin)) then
                  rdum       = rdum + abs(cdum-selft(k,i))
                  selft(k,i) = cdum
                  done(k)    = .true.
                endif
              enddo
            endif
          enddo
        enddo
        Rif(lsym) write(6,'(/A,ES8.1)') 'GT symmetry fulfilled to within',rdum/size(selft)
      endif
      ! restore
      Rend      
      call set_stage(1)
# ifdef MPI      
      call Mcast(oselft) ; call Mcastl(selft)
      call Mcastl(freqt) ; call Mcast(nfreqt)
      call Mcast(fullt)  ; call Mcastl(blockt)
# endif
      end subroutine read_sigt

c     ------------

c
c     Adds GT self-energy to selfc
c     We have to restage (selfc,oselfc,freq,nfreq,job1%full,block) by (selft,oselft,freqt,nfreqt,full,blockt) and restore afterwards.
      subroutine add_sigt(self,dself,e,i,j,iblock)
      use global
      implicit none
      integer,    intent(in)    :: i,j,iblock
      complex_dp, intent(inout) :: self,dself
      complex_dp, intent(in)    :: e
      complex_dp                :: cdum,dcdum
      integer                   :: i1
      logical                   :: full
      if(.not.allocated(selft)) Bug('Array selft not allocated.')
      full = job1%full
      i1   = block(i,iblock)
      ! restage
      call set_stage(2)
      ! add sigt
      cdum  = 0
      dcdum = 0
      if(fullt.eqv.full)        then ; call selfenergy_interpolate(cdum,dcdum,e,2,i,j,iblock) ! GT(diag) and GW(diag)   (GT FULL and GW FULL not impl.)
      else
        if(.not.fullt.and.i==j) then ; call selfenergy_interpolate(cdum,dcdum,e,2,i1)         ! GT(diag) and GW FULL
        else if(fullt)          then ; Bug('Not implemented: GT FULL.')                       ! GT(FULL) and GT(diag)   (GT FULL not impl.)
        endif
      endif
      self  = self  + cdum        
      dself = dself + dcdum
      ! restore
      call set_stage(1)
      end subroutine add_sigt

c     ------------
      
c
c     Sets stage for GT definition (mode=0), GT calculation (mode=2), or resets stage for GW (mode=1)
      subroutine set_stage(mode)
      use file
      implicit none
      integer,    intent(in)        :: mode
      complex_dp, pointer_cnt, save :: pnt(:,:)
      real_dp,    allocatable, save :: frq(:)
      integer,    allocatable, save :: blk(:,:)
      integer,                 save :: nfrq,oself
      logical,                 save :: full
      if(all(mode/=[0,1,2])) Bug('Wrong mode.')
      if(mode==1) then ; if(.not.associated(pnt)) Bug('pnt not assigned.')
        selfc     => pnt  ; nullify(pnt)
        oselfc    =  oself
        job1%full =  full
        nfreq     =  nfrq ; deallocate(freq) ; allocate(freq(nfreq))
        freq      =  frq  ; deallocate(frq)  ; deallocate(block) ; allocate(block(size(blk,1),size(blk,2)))
        block     =  blk  ; deallocate(blk)
      else
        pnt       => selfc
        oself     =  oselfc
        full      =  job1%full
        nfrq      =  nfreq ; allocate(frq(nfrq))
        frq       =  freq  ; allocate(blk(size(block,1),nblock))
        blk       =  block
        if(mode==2) then   ; deallocate(freq,block)
          selfc     => selft
          oselfc    =  oselft
          job1%full =  fullt
          nfreq     =  nfreqt ; allocate(freq(nfreq))
          freq      =  freqt  ; allocate(block(size(blockt,1),size(blockt,2)))
          block     =  blockt
        endif
      endif
      if(mode>0.and.allocated(freq1)) then
        deallocate(freq1)
        allocate(freq1(nfreq))
        if(oselfc==2.or.oselfc==3) then ; freq1 = eks + sign(1d0,eks)*freq
        else if(oselfc==4)         then ; freq1 = freq
        endif
      endif
      end subroutine set_stage

c     ------------

c     Construct Hamiltonian HQP = HKS + SIGX + SIGC - VXC
c     Also returns sigmac (SIGC) and dsigmac.
      subroutine qp_hamiltonian(hqp,sigmac,dsigmac,eqp)
      use wrapper, only: diagonal
      implicit none
      complex_dp, intent(out) :: hqp(n,n),sigmac(n,n),dsigmac(n,n)
      complex_dp, intent(in)  :: eqp
      integer                 :: j,k
      do k = 1,n
        do j = 1,n
          call selfenergy_interpolate(sigmac(j,k),dsigmac(j,k),eqp,1,j,k,iblock)
        enddo
      enddo
      if(use_sym.and.oselfc==1) then
        call symmetrize_matrix( sigmac,n,iblock)
        call symmetrize_matrix(dsigmac,n,iblock)
      endif
      hqp = diagonal(hks) + sigmax + sigmac - vxc
      end subroutine qp_hamiltonian

c     ------------

c     Adiabatically switch on diagonal elements and apply scissor shift to other subspaces (not containing i) in hqp if istep < minstep.
      subroutine qp_adiabatic(hqp,i,istep)
      implicit none
      complex_dp, intent(inout) :: hqp(n,n)
      integer,    intent(in)    :: i,istep
      integer                   :: j,k
      if(istep<minstep) then
        if(istep<0) Bug('istep<0.')
        do k = 1,n ; do j = 1,n
          if(j==k) then
            if(j<i)       hqp(j,j) = hqp(j,j) - 0.1d0 * (minstep-istep)/minstep
            if(j>=i+ndeg) hqp(j,j) = hqp(j,j) + 0.1d0 * (minstep-istep)/minstep
          else
            hqp(j,k) = hqp(j,k) * istep / minstep
          endif
        enddo ; enddo
      endif
      end subroutine qp_adiabatic

c     ------------

c     Return eigenvector index that has the largest overlap with previous eigenvector
      subroutine maxolap(j)
      implicit none
      integer, intent(inout) :: j
      integer                :: k,l
      real_dp                :: olap(n/ndeg)
      do k = ndeg,n,ndeg
        olap(k/ndeg) = sum ( [ (abs(dot_product(evec0(:,j),evec(:,l)))**2,l=k-ndeg+1,k) ] )
      enddo
      j = maxloc(olap,1) * ndeg
      end subroutine maxolap

c     ------------

      subroutine write_bandinfo
      implicit none
      character(:),  allocatable :: lines(:)
      complex_dp,    allocatable :: cmt1(:,:,:,:)
      integer                    :: i,s,ic
# include "interface/band_info.inc"
      call getkey(inp,'PROJECT',lines,status=i)
      if(i==0)      then ; Bug('PROJECT not found.')
      else if(i==1) then ; allocate(character(80) :: lines(ntype)) ; lines = '[s,p,d,f,g]'
      endif
      allocate( cmt1(maxlmindx,ncent,n,nspin3) )
      do i = 1,n
        call wavefunction_mt(cmt1(:,:,i,:),maxlmindx,0,job1%band(block(i,iblock)),ikptq,ispin)
      enddo
      do s = 1,nspin3
        do ic = 1,ncent
          cmt1(:,ic,:,s) = matmul ( cmt1(:,ic,:,s) , evec )
        enddo
      enddo
      call band_info(lines,cmt1,spin=ispin,band=job1%band(block(:n,iblock)),title=
     &  trim(job1%label)//' kpt: '//Chr(ikptq)//', block: '//Chr(iblock))
      deallocate(cmt1)
      end subroutine write_bandinfo      

      subroutine write_bandinfo_old
      use util, only: chr
      implicit none
      complex_dp, allocatable :: proj(:,:,:,:),a(:)
      real_dp,    allocatable :: olap_mt(:,:)
      integer                 :: itype,ieq,ic
      integer                 :: i,j,l,m,ib,jb,nn,ii
      real_dp                 :: intgrf
      write(6,'(/A)') 'Band info (s,p,d,f,g character for each atom type):'
      allocate ( proj(n,n,0:4,ntype) )
      proj = 0
      ic   = 0
      do itype = 1,ntype
        do ieq = 1,neq(itype)
          ic = ic + 1
          ii = 0
          do l = 0,lcut(itype)
            nn = nindx(l,itype)
            allocate ( olap_mt(nn,nn),a(nn) )
            do i = 1,nn
              do j = 1,nn
                olap_mt(i,j) = intgrf( bas1(:,i,l,itype,ispin) * bas1(:,j,l,itype,ispin) +
     &                                 bas2(:,i,l,itype,ispin) * bas2(:,j,l,itype,ispin) , itype )
              enddo
            enddo
            do m = -l,l
              if(l<=4) then
                do i = 1,n   ; ib = job1%band(block(i,iblock))
                  a = matmul ( olap_mt , cmt(ii+1:ii+nn,ic,ib,ikindx,ispin) )
                  do j = 1,n ; jb = job1%band(block(j,iblock))
                    proj(j,i,l,itype) = proj(j,i,l,itype) + dot_product ( cmt(ii+1:ii+nn,ic,jb,ikindx,ispin) , a )
                  enddo
                enddo
                if(l_soc) then
                  do i = 1,n   ; ib = job1%band(block(i,iblock))
                    a = matmul ( olap_mt , cmt(ii+1:ii+nn,ic,ib,ikindx,2) )
                    do j = 1,n ; jb = job1%band(block(j,iblock))
                      proj(j,i,l,itype) = proj(j,i,l,itype) + dot_product ( cmt(ii+1:ii+nn,ic,jb,ikindx,2) , a )
                    enddo
                  enddo
                endif
              endif
              ii = ii + nn
            enddo
            deallocate ( olap_mt,a )
          enddo
        enddo
      enddo
      do i = 1,n
        ! diagonal
        write(6,'(I3'NoA) job1%band(block(i,iblock))
        do itype = 1,ntype
          m = min(4,lcut(itype))
          write(6,'('//chr(m+1)//'F8.4'NoA) (real(proj(i,i,l,itype)),l=0,m)
        enddo
        write(6,'(A)') '  (KS)'
        ! full
        write(6,'(I3'NoA) job1%band(block(i,iblock))
        do itype = 1,ntype
          m = min(4,lcut(itype))
          write(6,'('//chr(m+1)//'F8.4'NoA) (real(dot_product(evec(:,i),matmul(proj(:,:,l,itype),evec(:,i)))),l=0,m)
        enddo
        if(job1%type==J_HF) then ; write(6,'(A)') '  (HF)'
        else                     ; write(6,'(A)') '  (GW)'
        endif
      enddo
      deallocate ( proj )
      end subroutine write_bandinfo_old

c     ------------

      subroutine symmetrize_matrix(matrix,n,iblock,error)
      implicit none
      integer,    intent(in)              :: iblock,n
      real_dp,    intent(inout), optional :: error
      complex_dp, intent(inout)           :: matrix(n,n)
      complex_dp                          :: matsum(n,n)
      integer                             :: isub,jsub
      integer                             :: ikpt1,ikptq,ispin
      integer                             :: isym,nsym1
      integer                             :: ndeg
      integer                             :: i,j,i1,j1
      if(n/=sizeblock(iblock)) Error('Count error.')
      i     = block(1,iblock)
      ikptq = job1%kpt(i)
      ispin = job1%spin(i)
      i     = job1%band(i)
      ndeg  = deg(i,ikptq,ispin) ; if(ndeg<i) ndeg = deg(ndeg,ikptq,ispin)
      ndeg  = ndeg - i + 1
      ! rotate and add
      matsum = matrix
      nsym1  = 1
      do isym = 2,nsym
        ikpt1 = kptsym(ikptq,isym) ; if(ikpt1/=ikptq) cycle
        nsym1 = nsym1 + 1
        do i = 1,n,ndeg   ; isub = psub(block(i,iblock))
          do j = 1,n,ndeg ; jsub = psub(block(j,iblock))
            i1 = i + ndeg - 1
            j1 = j + ndeg - 1
            if(isym>nsymt) then
              matsum(i:i1,j:j1) = matsum(i:i1,j:j1) + transpose (
     &                            matmul(  transpose(conjg(irrep_sub(:ndeg,:ndeg,jsub,isym))),
     &                            matmul(matrix(j:j1,i:i1),irrep_sub(:ndeg,:ndeg,isub,isym))) )
            else
              matsum(i:i1,j:j1) = matsum(i:i1,j:j1) +
     &                            matmul(  transpose(conjg(irrep_sub(:ndeg,:ndeg,isub,isym))),
     &                            matmul(matrix(i:i1,j:j1),irrep_sub(:ndeg,:ndeg,jsub,isym)))
            endif
          enddo
        enddo
      enddo
      if(present(error)) error = error + sum(abs(matrix-matsum/nsym1))
      matrix = matsum / nsym1
      end subroutine symmetrize_matrix

c     ------------

      subroutine write_qsgw
      use wrapper
      use file
      implicit none
      MCOMPLEX_dp, allocatable :: self(:,:),vxc(:,:),sxc(:,:)
      complex_dp               :: s1,s2,ds1,ds2,ei,ej
      real_dp                  :: errsym
      integer                  :: iblock,iunit,iselfc0,iselfx0
      integer                  :: i,k,s,j,i1,j1,n
# ifdef INV
      complex_dp,  allocatable :: cself(:,:)
# endif
      if(oselfc==1.and.npole/=0)
     &  write(6,'(A)') 'Calculate Hermitian self-energy for QSGW with numerical fits. This may take some time.'
      errsym  = 0
      iunit   = fopen('spex.qsgw',status='unknown')
      iselfc0 = 0
      iselfx0 = 0
      if(job1%type==J_HF) then
        write(iunit,'(A)') '# Energy-independent and Hermitian xc operator for scHF: \Sigmax - Vxc'
      else if(job1%type==J_PBE0) then
        write(iunit,'(A)') '# Energy-independent and Hermitian xc operator for PBE0: (\Sigmax-Vx) / 4'
      else
        write(iunit,'(A)') '# Energy-independent and Hermitian xc operator for QSGW: \Sigma - Vxc'
      endif
      write(iunit,'(A)') '#  s   k  b1  b2'
      do iblock = 1,nblock
        n    = sizeblock(iblock) ; allocate ( self(n,n),vxc(n,n),sxc(n,n) )
        self = unpackmat ( selfx(iselfx0+1:iselfx0+n*(n+1)/2) )
        k    = job1%kpt(block(1,iblock))
        s    = job1%spin(block(1,iblock))
        if(ovxc==0) then ; call read_vxc(vxc,job1%band(block(:n,iblock)),n,kptp(k),s,.false.)
        else             ; call calc_vxc(vxc,job1%band(block(:n,iblock)),n,kptp(k),s,.false.)
        endif
        if(ovxc==0) then ; call read_vxc(sxc,job1%band(block(:n,iblock)),n,kptp(k),s,l_qsgw)
        else             ; call calc_vxc(sxc,job1%band(block(:n,iblock)),n,kptp(k),s,l_qsgw)
        endif
        if(job1%type==J_PBE0) self = self / 4
        sxc  = sxc  - vxc ! Hermitian self-energy
        self = self - vxc
        if(job1%type==J_GW) then
          do j = 1,n
            do i = 1,j
              i1 = job1%band(block(i,iblock))
              j1 = job1%band(block(j,iblock))
              ei = ene(i1,k,s) - efermi
              ej = ene(j1,k,s) - efermi
              if(i/=j) then
                call selfenergy_interpolate(s1,ds1,ei,0,i,j,iblock)
                call selfenergy_interpolate(s2,ds2,ei,0,j,i,iblock)
                self(i,j) = self(i,j) + ( s1 + conjg(s2) ) / 4
                self(j,i) = self(j,i) + ( s2 + conjg(s1) ) / 4
                call selfenergy_interpolate(s1,ds1,ej,0,i,j,iblock)
                call selfenergy_interpolate(s2,ds2,ej,0,j,i,iblock)
                self(i,j) = self(i,j) + ( s1 + conjg(s2) ) / 4
                self(j,i) = self(j,i) + ( s2 + conjg(s1) ) / 4
              else
                call selfenergy_interpolate(s1,ds1,ei,0,i,j,iblock)
                s1        = s1 + ds1 * ( (self(i,i)-sxc(i,i)+s1) / (1-ds1) ) ! take self-energy at quasiparticle energy (Taylor 1st order); comment this line to take the self-energy at the KS energy (
                self(i,i) = self(i,i) + ( s1 + conjg(s1) ) / 2
              endif
            enddo
          enddo
        endif
        iselfc0 = iselfc0 + n**2
        iselfx0 = iselfx0 + n*(n+1)/2
# ifdef INV
        if(use_sym) then
          allocate ( cself(n,n) )                       ; cself =  self
          call symmetrize_matrix(cself,n,iblock,errsym) ;  self = cself
          deallocate ( cself )
        endif
        write(iunit,'(4I4, F18.13)') ((s,k,job1%band(block(i,iblock)),job1%band(block(j,iblock)),self(i,j),i=1,j),j=1,n)
# else
        if(use_sym) call symmetrize_matrix(self,n,iblock,errsym)
        write(iunit,'(4I4,2F18.13)') ((s,k,job1%band(block(i,iblock)),job1%band(block(j,iblock)),self(i,j),i=1,j),j=1,n)
# endif
        deallocate ( self,vxc,sxc )
      enddo
      call fclose(iunit)
      write(6,'(/A,ES8.1,A)') 'New QSGW self-energy written to spex.qsgw.  (',errsym,' )'
      if(l_soc) then
        iunit = fopen('spex.su2',status='unknown')
        do i = 1,nsym
          write(iunit,'(I3)') i
          write(iunit,'(2(F24.16,F21.16))') sym(i)%esoc
        enddo
        call fclose(iunit)
      endif
      end subroutine write_qsgw

c     ------------

      subroutine write_spectral(ikpt,ispin)
      use file
      implicit none
      integer, intent(in) :: ikpt,ispin
      integer             :: i
      real_dp             :: w
      call write_spectral_head(iunits,ikpt,ispin,'Spectral function',.false.)
      w = wspec(1)
      do i = 1,nspec
        write(iunits,'(F12.7,F28.10)') w*escale,spec(i)
        w = w + wspec(3)
      enddo
      write(iunits,'(/)')
      end subroutine write_spectral

      subroutine write_spectral_head(iunit,ikpt,ispin,title,lband)
      use global, only: kpt,escale
      use file
      implicit none
      integer, intent(in) :: iunit,ikpt,ispin
      logical, intent(in) :: lband
      character(*)        :: title
      character(4)        :: eunit
      integer             :: ikpt0
      if(escale==hartree) then ; eunit = '(eV)'
      else if(escale==1)  then ; eunit = '(Ha)'
      else                     ; eunit = '    '
      endif
      ikpt0 = get_ikptq0(ikpt)
      write(iunit,'(A)')        '# '//trim(title)
      write(iunit,'(A,9F10.5)') '# lattvec:',lat
      write(iunit,'(A,3F10.5)') '# k point:',kpt(:,ikpt)
      write(iunit,'(A,I6'NoA)   '# k index:',abs(ikpt0) ; if(ikpt0<=0) write(iunit,'(A'NoA) '  (additional)'
      write(iunit,*)
      if(nspin==2) write(iunit,'(A,I2)') '# spin:',ispin
      if(lband)    write(iunit,'(A,I5)') '# band:',ibandq
      write(iunit,'(A,11X,A)') '# Freq. '//eunit,trim(title)
      end subroutine write_spectral_head

c     ------------

c     Writes interpolated spectral function to spectralw (e.g., gnuplot: cbrange [0:?]; plot "spectralw" u 1:2:3 lw 1 palette w l)
      subroutine spectral_wannier_interpolation(spin)
      use file
      use key
      use wrapper
      Mpi2( use util, only: chr )
      implicit none
      Mpi( include 'mpif.h' )
      integer,     intent(in)  :: spin
      integer,     allocatable :: lpt(:,:),ldeg(:,:)
      complex_dp,  allocatable :: hamr(:,:,:),a(:),b(:)
      complex_dp               :: ham (nwanband,nwanband,nkpti),hamk(nwan,nwan)
      MCOMPLEX_dp              :: ham0(nwanband,nwanband,nkpti),iham(nwan,nwan),eigv(nwan,nwan)
      integer                  :: self(nwanband,nwanband,nkpti)
      real_dp                  :: w,spec,qvec0(3),qpath,qstep,eig(nwan)
      real_dp,     allocatable :: qvec(:,:)
      integer                  :: iunit,iblock
      integer                  :: ikpt,ispin,iselfc,iselfx,iselfc0,iselfx0
      integer                  :: i,j,iw,ib,jb,np,nlpt,ndiff,nw,w1,w2,nq,iq
      character(16)            :: filename
      character(256)           :: qfile
      real                     :: cputime
# ifdef MPI
      character(80)            :: line
      integer                  :: Merr,Mstat(mpi_status_size),iunit1
# endif
# include "interface/getqlist.inc"
# include "interface/wannier_spatial.inc"

      qfile = ' '
      Rcall getkey(inp,'INTERPOL', qfile, section='WANNIER', status=i)
      Rcall getqlist(qvec,nq,100,qfile)
      Rif(nq==0) then
        Warn('No list of q vectors provided (KPTPATH or file). Wannier interpolation skipped.')
        deallocate(qvec)
        return
      endif
      Mpi( call Mcast(nq) ; call Mcastl(qvec) )

      Rwrite(6,'(A'NoA) 'Interpolating spectral function'
      call cpu_time(cputime)

      nw = 0
      do ispin = 1,nspin1 ; if(spin/=0.and.ispin/=spin) cycle
        nw = nw + nspec Mpi(/Msize) ! for progress bar
      enddo

      do ispin = 1,nspin1 ; if(spin/=0.and.ispin/=spin) cycle

        Rbegin

        ! HF Hamiltonian -> ham0
        ! (a) hks - vxc
        ham0 = 0
        do ikpt = 1,nkpti
          if(job1%type/=J_GT) then
            if(ovxc==0) then ; call read_vxc(ham0(:,:,ikpt),[(i,i=wanbandi,wanbandf)],nwanband,ikpt,ispin,l_qsgw)
            else             ; call calc_vxc(ham0(:,:,ikpt),[(i,i=wanbandi,wanbandf)],nwanband,ikpt,ispin,l_qsgw)
            endif
            ham0(:,:,ikpt) = ham0(:,:,ikpt) + deltaex * (2*ispin-3)/2 * identity(nwanband)
          endif
          ham0(:,:,ikpt) = ham0(:,:,ikpt) + alignvxc * identity(nwanband)
          ham0(:,:,ikpt) = diagonal(ene(wanbandi:wanbandf,ikpt,ispin)-efermi) - ham0(:,:,ikpt)
        enddo
        ! (b) sigmax / sigmac
        self    = 0
        iselfx0 = 0
        iselfc0 = 0
        do iblock = 1,nblock
          ib   = block(1,iblock)
          ikpt = job1%kpt(ib)
          n    = sizeblock(iblock)
          if(job1%spin(ib)==ispin.and.ikpt<=nkpti) then
            iselfx = iselfx0
            iselfc = iselfc0
            do j = 1,n   ; jb = job1%band(block(j,iblock)) - wanbandi + 1
              do i = 1,n ; ib = job1%band(block(i,iblock)) - wanbandi + 1
                iselfc = iselfc + 1
                if(i<=j) iselfx = iselfx + 1
                if(min(ib,jb)>=1.and.max(ib,jb)<=nwanband) then
                  if(job1%type/=J_GT) then
                    if(i<=j) ham0(ib,jb,ikpt) = ham0(ib,jb,ikpt) +        selfx(iselfx)
                    if(i< j) ham0(jb,ib,ikpt) = ham0(jb,ib,ikpt) + MCONJG(selfx(iselfx))
                  endif
                  self(ib,jb,ikpt) = iselfc
                endif
              enddo
            enddo
          endif
          iselfx0 = iselfx0 + n*(n+1)/2
          iselfc0 = iselfc0 + n**2
        enddo
        
        Rend
        Mpi( call Mcast(ham0); call Mcast(self) )

        filename = 'spectralw'
        if(nspin1==2) write(filename(10:),'(I1)') ispin
        MnoR( filename(len_trim(filename)+1:) = '.'//adjustl(chr(Mrank)) )
        iunit = fopen(filename,status='unknown')
        
        Rif(nspin1==2) then
        if(ispin==1) write(iunit,'(A)') '# Spin up'
        if(ispin==2) write(iunit,'(A)') '# Spin down'
      endif
      
      MrangeDef1(w1,w2,nspec)
      w = wspec(1) + wspec(3) * (w1-1)
      do iw = w1,w2
        
        ! QP Hamiltonian: add sigmac to ham0 -> ham
        ham = 0
        do ikpt = 1,nkpti
          do i = 1,nwanband
            do j = 1,nwanband ; jb = j ; if(.not.job1%full) then ; if(i/=j) cycle ; jb = 1 ; endif
              iselfc = self(i,j,ikpt)  ; if(iselfc==0) cycle
              if(oselfc==1) then
                np = pade_w(i,jb,ikpt,ispin)%n ; if(np==0) cycle ; allocate(a(np),b(np))
                a  = pade_w(i,jb,ikpt,ispin)%a(:np)
                b  = pade_w(i,jb,ikpt,ispin)%b(:np)
                if(w<0) then
                  a(:np) = conjg(a(:np))
                  b(:np) = conjg(b(:np))
                  sigc   = selfc(iselfc,1) ; if(sigc/=0) a(:np) = a(:np) * sigc/conjg(sigc)
                endif
                sigc = sum ( a(:np) / ( w - b(:np) ) )
                deallocate(a,b)
              else
                call interpolate(sigc,dsigc,w,selfc(iselfc,:),freq1,nfreq)
              endif
              ham(i,j,ikpt) = ham0(i,j,ikpt) + sigc
            enddo
          enddo
        enddo
        
        ! Fourier transform ham to R -> hamr
        if(job1%full) then
          call wannier_spatial(hamr,ham,lpt,nlpt,ldeg,ndiff,ispin,1)
        else
          call wannier_spatial(hamr,[((ham(i,i,ikpt),i=1,nwanband),ikpt=1,nkpti)],lpt,nlpt,ldeg,ndiff,ispin,0)
        endif
        
        ! Interpolate along kpt_path
        qpath = 0
        qvec0 = qvec(:,1)
        do iq = 1,nq
          qstep = sqrt(sum(matmul(rlat,qvec(:,iq)-qvec0)**2))
          qpath = qpath + qstep                               ; call wannier_interpolate(hamk,hamr,nwan**2,lpt,nlpt,qvec(:,iq))
          iham  = ( hamk - conjg(transpose(hamk)) ) / (2*img) ; call diagonalize(eigv,eig,iham)
          hamk  = ( hamk + conjg(transpose(hamk)) ) / 2
          eig   = max( eig*sign(1d0,-w) , wspec(4) )
          iham  = matmul(eigv,matmul(diagonal(eig),MCONJG(transpose(eigv))))
          hamk  = hamk + img * iham
          hamk  = w*identity(nwan) - hamk                     ; call inverse(hamk)
          spec  = sum ( [ (imag(hamk(j,j)),j=1,nwan) ] ) / pi * sign(1d0,-w)
          qvec0 = qvec(:,iq)
          write(iunit,'(2F12.7,F30.10)') qpath,w*escale,spec
        enddo
        do i = 23*(iw-1)/nw+1,23*iw/nw ; Rwrite(6,'(A'NoA) '.' ; enddo ! Progress bar
          write(iunit,*)
          w = w + wspec(3)
        enddo
        
        deallocate(lpt,hamr,ldeg)
        
# ifdef MPI
        if(Mrank==0) then
          do i = 1,Msize-1
            call mpi_recv(filename,16,mpi_character,i,0,Mcomm,Mstat,Merr)
            iunit1 = fopen(filename,status='old')
            do
              read(iunit1,'(A)',iostat=j) line ; if(j/=0) exit
              write(iunit,'(A)') trim(line)
            enddo
            call fclose(iunit1,status='delete')
          enddo
        else
          call fclose(iunit)
          call mpi_send(filename,16,mpi_character,0,0,Mcomm,Merr)
        endif
# endif
        Rcall fclose(iunit)
        
      enddo
      
      deallocate(qvec)

      Rcall cpu_done(cputime)
      end subroutine spectral_wannier_interpolation

c     ------------

c     Return diagonal element (indx) of self-energy. (mode=1: selfx, mode=2: selfc)
      function get_self(indx,mode)
      implicit none
      complex_dp          :: get_self
      integer, intent(in) :: indx,mode
      complex_dp          :: cdum
      real_dp             :: eks
      integer             :: iblock,iself,iband,ikpt,ispin,iband1,ikpt1,ispin1
      integer             :: i,j,n
      iband = job1%band(indx)
      ikpt  = job1%kpt (indx)
      ispin = job1%spin(indx)
      eks   = ene(iband,ikpt,ispin) - efermi
      iself = 0
      do iblock = 1,nblock
        ikpt1  = job1%kpt (block(1,iblock))
        ispin1 = job1%spin(block(1,iblock))
        n      = sizeblock(iblock)
        do j = 1,n   ; iband1 = job1%band(block(j,iblock))
          do i = 1,n ; if(mode==1.and.i>j) exit
            iself = iself + 1
            if(all([i,iband1,ikpt1,ispin1]==[j,iband,ikpt,ispin])) then
              if(mode==1) then ; get_self = selfx(iself)
              else             ; call selfenergy_interpolate(get_self,cdum,eks*(1d0,0d0),0,i,j,iblock)
              endif
              return
            endif
          enddo
        enddo
      enddo
      Bug('State index not found in JOB.')
      end function get_self

c     ------------

c
c     Returns index of the array selfc for the self-energy block blk and matrix element i,j.
      function index_selfc(blk,i,j)
      implicit none
      integer             :: index_selfc
      integer, intent(in) :: blk,i,j
      integer             :: iblock,n
      index_selfc = 0
      do iblock = 1,nblock
        n = sizeblock(iblock)
        if(iblock<blk) then ; index_selfc = index_selfc + n**2
        else                ; index_selfc = index_selfc + (j-1)*n+i ; return
        endif
      enddo
      Bug('Did not find block index for: '//Chr(blk)//' '//Chr(i)//','//Chr(j))
      end function index_selfc

c     ------------

c
c     Returns k-point index:
c     = ikptq  if ikptq <= nkpt
c     = 0      if ikptq  = nkpt+1 and nkpt+1st k point is different from all others
c     = -ikpt  if ikptq  = nkpt+1 and nkpt+1st k point is identical to kpt(:,ikpt)
      function get_ikptq0(ikptq)
      use global, only: kpt,nkpt
      implicit none
      integer             :: get_ikptq0
      integer, intent(in) :: ikptq
      integer             :: k
      if(ikptq<=nkpt) then ; get_ikptq0 = ikptq
      else                 ; get_ikptq0 = 0 ; do k = 1,nkpt ; if(sum(abs(kpt(:,ikptq)-kpt(:,k)))<1d-8) get_ikptq0 = -k ; enddo
      endif
      end function get_ikptq0

c     ------------

c
c     Returns pade approximant in pade(:n,:m) of the self-energy block with offset iselfc0.
c     (normally, n=m, unless m=1. Then pade(:n,1) is essentially one-dimensional.)
      subroutine prepare_pade(pade,n,m,iselfc0,state,lwrite0)
      implicit none
      integer,         intent(in)           :: n,m,iselfc0
      type(padetype),  intent(out)          :: pade(:,:)
      type(statetype), intent(in), optional :: state(:)
      logical,         intent(in), optional :: lwrite0
      complex_dp,      allocatable          :: a(:),b(:)
      complex_dp                            :: cpade(nfreq+1)
      logical                               :: lwrite
      integer                               :: np,i,j,j1
      if(npole==0) then ; np = (nfreq+1)/2 * smooth(1)
      else              ; np = npole
      endif
      allocate( a(np),b(np) )
      if(m/=1.and.n/=m.or.n<=0) Bug('Pade dimensions incorrect.')
      lwrite = .false. ; if(present(lwrite0)) lwrite = lwrite0
      pade%c = 0
      pade%n = 0
      iselfc = iselfc0
      do j = 1,m
        do i = 1,n
          iselfc = iselfc + 1
          ldum   = lwrite.and.job1%full.and.i==j andR
          if(ldum.and.ndeg>0) ldum = ldum.and.mod(i,ndeg)==0
          if(npole==0) then ! Pade
            call pade_init(cpade,freq(nfreq:1:-1)*img,selfc(iselfc,nfreq:1:-1),nfreq,2)
            call pade_poles(b,a,np,freq(nfreq:1:-1)*img,selfc(iselfc,nfreq:1:-1),cpade,nfreq,smooth(1),2,ldum)
          else              ! Fit
            pade(i,j)%n = npole
            call continuation(freq*img,selfc(iselfc,:),nfreq,a,b,pade(i,j)%c,np,0,allowup,fitconst,fitconstr,ldum)
          endif
          allocate(pade(i,j)%a(np),pade(i,j)%b(np))
          pade(i,j)%n = np
          pade(i,j)%a = a(:np)
          pade(i,j)%b = b(:np)
          Rif(lwrite.and.wrtinfo) then
            if(.not.present(state)) Bug('Argument state missing.')
            j1 = j ; if(m/=n) j1 = i
            call write_seinfo(freq*img,selfc(iselfc,:),nfreq,a,b,pade(i,i)%c,np,state(i),state(j1))
          endif
          Rif(ldum) then
            write(6,'(A,F15.10)') '  Integrability of self-energy      :',imag ( sum(a(:np))      )
            if(npole/=0)
     &      write(6,'(A,F15.10)') '  Discontinuity at w=0 (self-energy):',imag ( sum(a(:np)/b(:np))    )
            write(6,'(A,F15.10/)')'  Discontinuity at w=0 (derivative) :',imag ( sum(a(:np)/b(:np)**2) )
          endif
        enddo
      enddo
      deallocate(a,b)
      end subroutine prepare_pade

c     ------------

# ifdef MPI
#   define CB call Mcast
#   define CA call Mcastl
      subroutine Mbroadcast_inp
      implicit none
      CB(fitconst) ; CB(fitconstr) ; CB(allowup) ; CB(npole) ; CB(bandinfo) ; CB(l_qsgw) ; CB(alignvxc)
      CB(l_inter)
      end subroutine Mbroadcast_inp
#   undef CB
#   undef CA
# endif

c     ------------

      end
