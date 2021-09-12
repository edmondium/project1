c--------------------------------------------------------------------------------
c Copyright (c) 2021 Peter Grünberg Institut, Forschungszentrum Jülich, Germany
c This file is part of SPEX and available as free software under the conditions
c of the MIT license as expressed in the LICENSE file in more detail.
c--------------------------------------------------------------------------------

c # define gammatest
c # define test_mtrafo

c Main routine for calculating
c - the self-energy (correlation part),
c - the spectra.
c
c Calculates exchange contribution, too, if NOSTORE is specified (default for MPI).
c
c The storage of the quantities susceptibility (suscep(c)), dielectric matrix (dielec(c)) and screened potential
c (screen(c)) needs a lot of memory. Therefore, they are calculated one after another and occupy the same memory
c space (matrix(c)).
c
c If iand(restart,R_cor)/=0, the dielectric matrix is read from spex.cor, if it exists.
c If iand(restart,W_cor)/=0, the dielectric matrix is written to spex.cor.
c If iand(restart,R_sig?)/=0, the exchange (?=x) and correlation (?=c) self-energies are read from spex.sigx and spex.sigc, respectively, if they exist.
c If iand(restart,W_sig?)/=0, the self-energies are written to spex.sigx(c).
c etc.
c
c Logical table:
c
c      HF HFE PBE0 SX COSX GW RPA SCRW SPEC Line (approx)
c sigx 1  1   1    1  1    2  2   0    0    166  read sigx file
c sigc 0  0   0    0  0    1  0   0    0    166  read sigc file
c cor  0  0   0    1  1    1  1   1    0    260  read cor file
c bas  0  0   0    1  1    1  1   1    1    276  basis needed
c X0   1  1   1    0  0    2  2   0    0    281  call x selfenergy
c skip 1  1   1    0  0    2  0   0    0    288  skip c (suscep ...)
c X    0  0   0    1  1    2  0   0    0    790  eval(1): x needed in c selfenergy
c C    0  0   0    0  0    1  0   0    0    791  eval(2): c needed in c selfenergy
c
c 0 - no
c 1 - yes
c 2 - perhaps
c
c
c PLASMA METAL: Scales H(w)=Whead(w) (or rather H(w)-4pi) such that metallic
c screening is in effect at w=0. "plasma" is set to zero, so H(0)=H>0 and
c H(w)=4pi for w->infty. We transform
c H(w) := 4pi-4pi[H(w)-4pi]/(H-4pi) = 4pi[H-H(w)]/(H-4pi). This gives
c H(0)=0 and H(w)=4pi for w->infty. However, for "divergence_h/c" we need
c the macro dielectric function D(w)=4pi/H(w). The corresponding trafo is thus
c D(w) := D(w)(D-1)/[D-D(w)] with D=4pi/H. This gives
c D(0)=infty and D(w)=1 for w->infty.
c In the more general case where D(w) is a tensor, we interpret the above trafo
c as a matrix equation and symmetrize afterwards. (Only the symmetrized part
c [D(w)+D(w)^T]/2 is relevant in kD(w)k.) If D(w) is diagonal, the equation
c separates into three scalar equations. If it is not (in particular if the
c eigenvectors of D(w) differ), one might still think about the symmetrization
c again. (Should it be done with sqrtmat?)
c
c
c Note: returns immediately if HF,HFE,PBE0+coulstore

c MPI: Define kind of ordering of output (MPI_order 0, 1 for ii==0 and ii>0, respectively, see Mwrite_n)
# define MPI_order 0

# if MPI_order == 0
#   define IKPT 0
# elif MPI_order != 1
#   error unknown MPI_order
# endif

# include "cppmacro.h"
# include "jobtype.h"
# include "restype.h"

      subroutine correlation(job1)

      use global
      use arrays
      use wrapper
      use key
      use file
      use util
      use, intrinsic :: iso_fortran_env
      Mpi ( use Mwrapper )
      Mpi2( use, intrinsic :: iso_c_binding )

      implicit none
# include "interface/selfenergy.inc"
# include "interface/define_freq.inc"
      type(jobtype), intent(in)  :: job1
      MCOMPLEX_dp,   pointer_cnt :: matrix(:,:),suscep(:,:),dielec(:,:),screen(:,:)             ! Hermitian matrices for imaginary frequencies
      MCOMPLEX_dp,   allocatable :: wing(:,:,:),head(:,:)                                       !
      complex_dp,    pointer_cnt :: matrixc(:,:,:),suscepc(:,:,:),dielecc(:,:,:),screenc(:,:,:) ! Complex matrices for real frequencies (for W and spectra)
      complex_dp,    allocatable :: wingc(:,:,:,:),headc(:,:,:)                                 !
      integer,       allocatable :: gpw(:),iarr(:)
      MCOMPLEX_dp,   pointer_cnt :: ctrafo(:,:)
      MCOMPLEX_dp,   allocatable :: wing2(:,:)
      MCOMPLEX_dp                :: head0(6),head1(3,3)
      character(256)             :: line
      character(6)               :: cjob
      logical                    :: symon,ldum,trafo,eval(3)
      logical,       allocatable :: lkpt(:,:)
      real_dp,       allocatable :: coul(:),sqcoul(:),eig(:)
      real_dp                    :: plasma,vec(3),contour_par(4)
      real                       :: cputime
      real_dp                    :: rdum
      complex_dp                 :: mat(3,3),h0(3,3)
      integer,       allocatable :: kpt_list(:)
      integer                    :: nbas
      integer                    :: iunit
      integer                    :: ikpt,ikpt1,kpt_up
      integer                    :: iblock,spin,typ
      integer                    :: i,j,k,nk,m
      integer                    :: g(3),g1(3),g2(3)
      real_dp                    :: ch2r
      complex_dp                 :: cdum,cdum1,cdum2
# ifdef MPI
      type(c_ptr)                :: ptr
      integer                    :: win_matrix,win_matrixc,win_ctrafo,Merr,Mcolor,Maxcolor,Mnkpt,Mfile
      integer,     allocatable   :: Msub(:)
      complex_dp,  pointer_cnt   :: Ninit_selfc(:),Ninit_matrixc(:)
      MCOMPLEX_dp, pointer_cnt   :: Ninit_coulomb0(:),Ninit_matrix(:),Ninit_ctrafo(:)
      real                       :: time_idle,time_maxidle
      integer                    :: step,irank_,color_
      logical                    :: pcycle_skip
# endif

      if(any(job1%type==[J_HF,J_HFE,J_PBE0]).and.coulstore) return
      if(any(job1%type==[J_GT,J_GOLD]))                     return
      if(    job1%type== J_SUSR .and. job1%kern=='BSE')     return

      plasma = 0

      spin = 0
      if(any(job1%type==J_SPEC)) then ; spin  = job1%spin(1)
      else if(job1%type==J_SCRW) then ; spinw = job1%spin(1)
      endif

      if(any(job1%type==J_SPEC)) then ; nk = size(job1%kpt) ; allocate(kpt_list(nk)) ; kpt_list = job1%kpt
      else                            ; nk = nkpti          ; allocate(kpt_list(nk)) ; kpt_list = [(i,i=1,nkpti)]
      endif
      if(any(kpt_list<=0)) Bug('Found non-positive kpoint index.')

c
c     eval(1)   : selfx   to be calculated
c     eval(2)   : selfc   to be calculated
c     eval(3)   : screenk is required (contracted W)
c     lkpt(k,i) : kpoint k contained in dump file i(=1,2,3)
      if     (any(job1%type==[J_HF,J_HFE,J_PBE0])) then ; eval = [ .not.coulstore , .false. , .false. ]
      else if(    job1%type== J_SX               ) then ; eval = [ .true.         , .false. , .false. ]
      else if(    job1%type== J_COSX             ) then ; eval = [ .true.         , .false. , .true.  ]
      else if(    job1%type== J_GW               ) then ; eval = [ .not.coulstore , .true.  , oibc/=0 .or. any(cores) ]
      else if(    job1%type== J_RPA              ) then ; eval = [ .not.coulstore , .true.  , .false. ]
      else                                              ; eval = [ .false.        , .true.  , .false. ]
      endif
      allocate(lkpt(size(kpt_list)+2,3))                               ! lkpt map : true (false) = already calculated (not calculated yet)
      lkpt                      = .false.                              ! nkpti+1 is index for contracted W
      lkpt(size(kpt_list)+2, 1) = any(job1%type==J_EXCH).and.coulstore ! nkpti+2 is index for core exchange: set to true if already calculated

      Rbegin

      write(6,'(//A)') '### subroutine: correlation ###'

      select case(job1%type)
        case(J_HF)       ; write(6,'(/A)') '< Hartree-Fock >'
        case(J_HFE)      ; write(6,'(/A)') '< Hartree-Fock energy >'
        case(J_GW)       ; write(6,'(/A)') '< GW self-energy >'
        case(J_SX)       ; write(6,'(/A)') '< SX self-energy >'
        case(J_COSX)     ; write(6,'(/A)') '< COSX self-energy >'
        case(J_PBE0)     ; write(6,'(/A)') '< PBE0 potential >'
        case(J_RPA)      ; write(6,'(/A)') '< RPA xc energy >'
        case(J_SUS)      ; write(6,'(/A)') '< Susceptibility >'
        case(J_SUSR)     ; write(6,'(/A)') '< Renormalized susceptibility '//trim(job1%kern)//' >'
        case(J_DIEL)     ; write(6,'(/A)') '< Dielectric matrix and its inverse >'
        case(J_SCR)      ; write(6,'(/A)') '< Screened potential >'
        case(J_SCRW)     ; write(6,'(/A)') '< Coulomb matrix in Wannier basis >'
        case default     ; Error('job number unknown.')
      end select

c
c     Define frequencies
      call define_freq(job1,contour_par,0,.true.,.true.)
      if(iand(restart,R_sigc)/=0.and.eval(2).and.job1%type==J_GW) call read_dump(0,ldum,0,'spex.sigc',job1)

      Rend
      Mpi( call Mbroadcast0 )

c
c     Core contribution to susceptibility
      if(any(cores).and.all(job1%type/=[J_HF,J_HFE,J_PBE0])) then
        Rcall susceptibility_core(spin,[img*freq,freqc],nfreq+nfreqc)
        Mpi( call mpi_barrier(Mcomm,Merr) )
      endif

c
c     Read RSITE from spex.inp and allocate (unless the Hubbard file is written in this run, in which case arrays are not needed)
      if(job1%type==J_SCRW) then
        call read_rsite
        Rcall getkey(inp,'HUBBARD',line,status=i) ; ldum=.false.
        Rif(i==2) inquire(file=line,exist=ldum)   ; Mpi( call Mcast(ldum) )        
        if(.not.ldum) then
          Allocate_ (   barew,(nwan,nwan,nwan,nwan,nsite)        ) ; barew   = 0
          Allocate_ ( screenw,(nwan,nwan,nwan,nwan,nfreqc,nsite) ) ; screenw = 0
        endif
      endif

c
c     Allocate arrays selfc/x
      if(job1%type==J_GW) then ! selfc
        i = size(job1%band)
        if(job1%full) then
          i = 0
          do iblock = 1,nblock
            j = sizeblock(iblock)
            if(i+dble(j)**2>huge(0)) Error('Integer overflow for selfc. Please report.')
            i = i + j**2
          enddo
        endif
        if(oselfc==1) then ; Nallocate0 ( selfc,(S_ i,nfreq  S_) )
        else               ; Nallocate0 ( selfc,(S_ i,nfreqr S_) )
        endif
      endif
      if(eval(1)) then ! selfx
        if(job1%full) then
          i = 0
          do iblock = 1,nblock
            j = sizeblock(iblock)
            i = i + j*(j+1)/2
          enddo
        else
          i = size(job1%band)
        endif
        Allocate_ ( selfx,(i) )
        selfx = 0
      endif

c
c     Allocate screenk_mt/pw
      if(eval(3)) then
        allocate ( screenk_mt(maxlmindxm,maxlmindxm,ncent,nfreq) ) ; screenk_mt = 0
        if(job1%type==J_COSX) then
          g1 =  huge(0)
          g2 = -huge(0)
          do k = 1,nkpt
            do j = 1,ngptm(k)
              do i = 1,ngptm(k)
                g  = gptm(:,pgptm(i,k)) - gptm(:,pgptm(j,k))
                g1 = min(g1,g)
                g2 = max(g2,g)
              enddo
            enddo
          enddo
          allocate ( screenk_pw(g1(1):g2(1),g1(2):g2(2),g1(3):g2(3)) )
          screenk_pw = 0
        endif
      endif

c
c     Read dump files if requested
      Rif(iand(restart,R_sigx)/=0.and.eval(1))                     call read_dump(1,lkpt(:,1),nkpti+2,'spex.sigx',job1) 
      Rif(iand(restart,R_sigc)/=0.and.eval(2).and.job1%type==J_GW) call read_dump(2,lkpt(:,2),nkpti+1,'spex.sigc',job1)
      Rif(iand(restart,R_wcou)/=0.and.allocated(barew))            call read_dump(3,lkpt(:,2),nkpti,  'spex.wcou',job1)
      Rif(iand(restart,R_ccou)/=0.and.eval(3)) then
        if(job1%type==J_COSX) then ; i = 1 ; else ; i = 2 ; endif
        if(lkpt(nkpti+1,i))   then ; lkpt(:nkpti,3) = .true. ; eval(3) = .false. ! If contracted part already obtained from dump file, we don't need contracted W.
        else                       ; call read_dump(4,lkpt(:,3),nkpti,'spex.ccou',job1)
        endif
      endif

c
c     Modify lkpt map : true (false) = already calculated or not needed (not calculated yet) [we later set eval=.not.lkpt(i,:)]
      do i = 1,size(lkpt,1)
        lkpt(i,:) = lkpt(i,:) .or. .not.eval
      enddo

      Mpi( call Mbroadcast1 ; if(job1%type==J_GW) then ; Nfence(selfc) ; endif )

c
c     Remove kpoints from kpt_list that have already been calculated
      nk = size(kpt_list)
      m  = 0
      do i = 1,nk
        if(.not.all(lkpt(i,:))) then
          m = m + 1
          kpt_list(m) = kpt_list(i)
        endif
      enddo
      call reallocate(kpt_list,m)

c
c     Gamma divergence
      if(any(job1%type==[J_EXCH,J_SX,J_COSX]).and.any(.not.lkpt(:2,1)).and.divergence==0) call gamma_divergence(.true.)

c
c     MPIKPT
# ifdef MPI

# ifdef old_mpikpt
      kpt_up = size(kpt_list)
      
      allocate(Msub(kpt_up))
      Msub   = 0
      Mcolor = 0
      nk     = count ( [ (any(.not.lkpt(ikpt,:)),ikpt=1,kpt_up) ] )
      if(all(job1%type/=J_SPEC).and.Osize>1.and.nk>1) then
        Rcall getkey(inp,'MPIKPT',ldum,default=.false.)
        call Mcast(ldum)
        if(ldum) call Mdistribute_kpt
        if(Mcolor/=0) then
          Rwrite(6,'(/A,'//chr(nkpti)//'(1X,A))') 'Subcommunicators for k loop:',(trim(chr(Msub(i))),i=1,nkpti)
          call begin_split_nodes(Mcolor)
          call cpu_time(time_idle)
        endif
      endif
      Maxcolor = maxval(Msub)

# else

      Mcolor   = 0
      Maxcolor = 0
      ldum     = .false.
      if(size(kpt_list)>1.and.Osize>1) then
        Rcall getkey(inp,'MPIKPT',ldum,default=.false.)
        call Mcast(ldum)
      endif
      if(ldum) then
        call Mdistribute_kpt(Mcolor,Maxcolor,kpt_list)
        if(Maxcolor>1) then
          Rwrite(6,'(A)') 'Number of MPIKPT subcommunicators: '//Chr(Maxcolor)
          Oif(Maxcolor<size(kpt_list)) call Mtask_init(size(kpt_list),Maxcolor,Mfile,Orank,Ocomm) ! we need dynamical work sharing
          call begin_split_nodes(Mcolor)
          call cpu_time(time_idle)
        else
          Rwrite(6,'(A)') 'No MPIKPT subcommunicators will be used.'
        endif
      endif

# endif

# endif

c
c     Prepare dump files:
c     wcoufiles/sigfiles: Truncate dump files (n=0); if iand(restart,W_sig?)/=0, write global dump file (Mcolor=0) or dump file of first subcommunicator (Mcolor=1)
c     ccoufiles:                 --"--             ; if iand(restart,W_cor)/=0,                                    -- " --
# define Ifwrt if(m>0 Mpi(.and.Mcolor<=1) )
      j = size(lkpt,1)
      m = count(lkpt(:j,1))
      Rif(iand(restart,W_sigx)/=0 .and. (m<j.and.eval(1)) ) then
              call write_dump(1,[0],                         0,'spex.sigx',job1 MpiC(Mcolor) MpiC(Maxcolor) )
        Ifwrt call write_dump(1,pack([(k,k=1,j)],lkpt(:j,1)),m,'spex.sigx',job1 MpiC(Mcolor) MpiC(Maxcolor) )
      endif
      m = count(lkpt(:j,2))
      Rif(iand(restart,W_sigc)/=0 .and. (m<j.and.eval(2).and.job1%type==J_GW) ) then
              call write_dump(2,[0],                         0,'spex.sigc',job1 MpiC(Mcolor) MpiC(Maxcolor) )
        Ifwrt call write_dump(2,pack([(k,k=1,j)],lkpt(:j,2)),m,'spex.sigc',job1 MpiC(Mcolor) MpiC(Maxcolor) )
      endif
      j = size(kpt_list)
      Rif(iand(restart,W_wcou)/=0 .and. (m<j.and.allocated(barew)) ) then
              call write_dump(3,[0],                         0,'spex.wcou',job1 MpiC(Mcolor) MpiC(Maxcolor) )
        Ifwrt call write_dump(3,pack([(k,k=1,j)],lkpt(:j,2)),m,'spex.wcou',job1 MpiC(Mcolor) MpiC(Maxcolor) )
      endif
      m = count(lkpt(:j,3))
      Rif(iand(restart,W_ccou)/=0 .and. (m<j.and.eval(3)) ) then
              call write_dump(4,[0],                         0,'spex.ccou',job1 MpiC(Mcolor) MpiC(Maxcolor) )
        Ifwrt call write_dump(4,pack([(k,k=1,j)],lkpt(:j,3)),m,'spex.ccou',job1 MpiC(Mcolor) MpiC(Maxcolor) )
      endif
# undef Ifwrt

c
c     Loop over k points

# ifdef old_mpikpt

      Mpi( Mnkpt = kpt_up )

      do ikpt1 = 1,kpt_up
        eval = .not.lkpt(ikpt1,:)

        if(any(job1%type==J_SPEC)) then ; ikpt = job1%kpt(ikpt1) ; if(ikpt==0) cycle
        else                            ; ikpt = ikpt1
        endif

        if(.not.any(eval) Mpi(.or.Mcolor/=Msub(ikpt1)) ) cycle ! cycle if not required (.not.eval) or done by other subcommunicator

# else

# ifdef MPI
      Mnkpt = size(kpt_list) ! needed for Mwrite_n
      ikpt1 = 0
      do
        if(Mcolor==0) then
          ikpt1 = ikpt1 + 1
        else
          if(Maxcolor<size(kpt_list)) then
            Rcall Mtask_indx(ikpt1,Mcolor,Maxcolor,Mfile) ! get next task from file
            call Mcast(ikpt1)
          else
            if(ikpt1/=0) exit ! only one loop
            ikpt1 = Mcolor
          endif
        endif
        if(ikpt1>size(kpt_list)) exit
# else
      do ikpt1 = 1,size(kpt_list)
# endif          
        
        ikpt = kpt_list(ikpt1)

# endif

        Mpi( Rif(Mcolor/=0) call Mwrite_n(0,IKPT,Mnkpt) )

        Rwrite(6,'(/A)')              '========='
        Rwrite(6,'(/A,I4,A,3F10.5/)') 'K point (',ikpt,' ):',kpt(:,ikpt)

        ! Define suffix for output files
        cjob = ' '
        if(maxval(job%indx)>1) cjob = trim(chr(job1%indx))
        if(size(kpt_list)>1)   cjob = trim(cjob)//'.'//chr(ikpt1,'I3.3')

        ! J_SPEC: Cycle if output file exists
        if(iand(restart,R_spec)/=0.and.any(job1%type==J_SPEC)) then
          Rbegin
          select case(job1%type)
            case(J_SUS)  ; line = 'suscep' //cjob
            case(J_SUSR) ; line = 'suscepR'//cjob
            case(J_DIEL) ; line = 'dielecR'//cjob
            case(J_SCR)  ; line = 'screen' //cjob
            case default ; Bug('Unkown job type.')
          end select
          inquire(file=line,exist=ldum)
          Rend
          Mpi( call Mcast(ldum) )
          if(ldum) then
            Rwrite(6,'(A)') 'Output file '//trim(line)//' exists. Skipping...'
            cycle
          endif
        endif

        ! RESTART: Read restart file if present ( typ = -1:none 0:coulomb0, 1:suscep, 2:dielec, 3:inverse(dielec) )
        typ = -1
        if(iand(restart,R_cor)/=0) then
          if(any(job1%type==[J_HF,J_HFE,J_PBE0])) then ; typ = 0
          else if(job1%type==J_SUS)               then ; typ = 1
          else if(job1%type==J_RPA)               then ; typ = 2
          else if(job1%type==J_GW)                then ; typ = 3 ; if(.not.eval(2)) typ = 0
          else                                         ; typ = 3
          endif
          call readcor(typ)
        endif

        ! Set current Coulomb matrix: calculate (NOSTORE) or copy -> coulomb0
        if(typ==-1) then
          Nallocate0 ( coulomb0,(S_ nbasm(ikpt)*(nbasm(ikpt)+1)/2 S_) )
          if(associated(coulomb)) then
            coulomb0 = coulomb(:size(coulomb0,1),ikpt)
          else
            call coulombmatrix(Win(coulomb0),nbasm(ikpt),1,ikpt,.false.)
          endif
          if(iand(restart,W_cor)/=0) call writecor(0)
        endif

        ! ldum = .true. if correlation not needed (then Coulomb eigenbasis not needed)
        ldum = any(job1%type==[J_HF,J_HFE,J_PBE0]).or.(job1%type==J_GW.and..not.eval(2))

        ! Determine basis set - eigenvectors of Coulomb matrix or plane waves - and whether we reduce its size.
        if(typ<=0.and..not.ldum) call get_basis

        ! Evaluate exchange self-energy for ikpt if necessary (and cycle)
        if(eval(1) .and. (ldum.or.job1%type==J_RPA)) then
          call  selfenergy(job1,ikpt,[.true.,.false.],0)
          Rif(iand(restart,W_sigx)/=0) call write_dump(1,[ikpt],1,'spex.sigx',job1 MpiC(Mcolor) MpiC(Maxcolor) )
#         ifdef CHECK_SENERGY
          Rwrite(*,*) selfx  ; Rread(*,*) ; Mpi( call mpi_barrier(Mcomm,i); ) selfx = 0
#         endif
          if(job1%type/=J_RPA) then
            Ndeallocate ( coulomb0 )
            goto 10
          endif
        endif

        ! RESTART: Matrix already calculated if typ>0; skip calculations
        select case(typ)
          case(1) ; goto 1
          case(2) ; goto 2
          case(3) ; goto 3
        end select

c
c       Matrix allocation
# ifdef CHECK
#   ifndef noWARN
#     warning(set max(nfreq[c],1) for CHECK!)
#   endif
#   define NFREQ  max(nfreq,1)
#   define NFREQC max(nfreqc,1)
# endif
        Rcall checkmem('matrix/matrixc', MBYTES*NFREQ*nbas*(nbas+1)/2 + 16d0*NFREQC*nbas**2)
        Rcall checkmem('matrix/matrixc',-MBYTES*NFREQ*nbas*(nbas+1)/2 - 16d0*NFREQC*nbas**2)
        Nallocate0 ( matrix, (S_ nbas*(nbas+1)/2,NFREQ  S_))
        Nallocate0 ( matrixc,(S_ nbas,nbas,      NFREQC S_))
# ifdef gammatest
        if(ikpt==1) then
# endif
        allocate  ( wing(3,nbas,nfreq),head(6,nfreq) )         ! head and wing not used if ikpt/=1 or spin>=3
        allocate  ( wingc(3,nbas,2,nfreqc),headc(3,3,nfreqc) ) ! head and wing not used if ikpt/=1 or spin>=3
# ifdef gammatest
        endif
# endif

        suscep  => matrix
        suscepc => matrixc

        call susceptibility ( ikpt,spin,nbas,ctrafo,coul,trafo,symon,plasma,nfreqc,freqc,   Win(matrixc),headc,wingc ,
     &                                                                      nfreq, freq*img,Win(matrix), head, wing  )

c                                                               _
c       APPROXPW: from now on we need < M | eigv > instead of < M | eigv >
        if(.not.fullpw) then
          Nfence(ctrafo)
          do i = 1,nbas ; Ncycle(i-1)
            call olapvec1(ctrafo(nbasp+1:,i),gptm(:,pgptm(:ngptm(ikpt),ikpt)),ngptm(ikpt))
          enddo
          Nfence(ctrafo)
        endif

c
c       Transform to Coulomb eigenbasis
        if(.not.trafo) then
          call cpu_time(cputime)
          Rwrite(6,'(A'NoA) 'Transform susceptibility matrices... '
          do i = 1,nfreq
            if(ikpt==1) then
              do j = 1,3
                wing(j,:,i) = matmul ( wing(j,:,i), MCONJG(ctrafo) )
              enddo
            endif
            Ncycle(i-1)
            call unitarytrafo( suscep(:,i) , ctrafo , 1 )
          enddo
          do i = 1,nfreqc
            if(ikpt==1) then
              do j = 1,3
                wingc(j,:,1,i) = matmul ( wingc(j,:,1,i), MCONJG(ctrafo) )
                wingc(j,:,2,i) = matmul ( wingc(j,:,2,i),        ctrafo  )
              enddo
            endif
            Ncycle(i-1)
            call unitarytrafo( suscepc(:,:,i) , ctrafo , 1 )
          enddo
          Rcall cpu_done(cputime)
        endif

        if(job1%type==J_SUS.and.iand(restart,W_cor)/=0.and.typ<=0) call writecor(1)

 1      if(wrtsuscep.or.any(job1%type==[J_SUS,J_SUSR])) then
          Nfence(matrix)
          Nfence(matrixc)
          Rif(ikpt==1.and.metal.and.spin<=2) then
            write(6,'(A)') 'correlation: Note! Head and wings must vanish. With the present implementation'
            write(6,'(A)') '             this is not the case for metals.'
          endif
          Rcall writeout(suscep,suscepc,nbas,job1%out,'suscep','Susceptibility')
          Nfence(matrix)
          Nfence(matrixc)
          if(job1%type==J_SUS) goto 4
        endif

c
c       If we treat the Gamma point, enforce suscep to be the body, i.e. set suscep(1,:), suscep(:,1), and wing(:,1) to zero
        if(ikpt==1) then
          Nfence(matrix)
          Nfence(matrixc)
          Rbegin
          if(nfreq/=0) then
            rdum = sum( [ (sum(abs(suscep(1+i*(i-1)/2,:))**2),i=1,nbas) ] )
            write(6,'(A,F20.15)')   'Residual head: ',sqrt(sum(abs(wing(:,1,:))**2)/nfreq)
            write(6,'(A,F20.15'NoA) 'Residual wings:',sqrt(rdum/(nbas*nfreq))
            if     (metal)                           then ; write(6,'(A)') '  (not small for a metallic system)'
            else if(.not.fullpw.and.any(lcutp<lcut)) then ; write(6,'(A)') '  (small for large GCUT and LCUT in sec. WFPROD)'
            else if(.not.fullpw)                     then ; write(6,'(A)') '  (small for large GCUT)'
            else if(any(lcutp<lcut))                 then ; write(6,'(A)') '  (small for large LCUT in sec. WFPROD)'
            else                                          ; write(6,*)
            endif
          endif
          if(nfreqc/=0) then
            write(6,'(A,F20.15)')   'Residual head: ',sqrt(sum(abs(wingc(:,1,:,:))**2)/nfreqc)
            write(6,'(A,F20.15'NoA) 'Residual wings:',sqrt((sum(abs(suscepc(:,1,:))**2)+sum(abs(suscepc(1,:,:))**2))/(nbas*nfreqc))
            if     (metal)                           then ; write(6,'(A)') '  (not small for a metallic system)'
            else if(.not.fullpw.and.any(lcutp<lcut)) then ; write(6,'(A)') '  (small for large GCUT and LCUT in sec. WFPROD)'
            else if(.not.fullpw)                     then ; write(6,'(A)') '  (small for large GCUT)'
            else if(any(lcutp<lcut))                 then ; write(6,'(A)') '  (small for large LCUT in sec. WFPROD)'
            else                                          ; write(6,*)
            endif
          endif
          Rend
          Nfence(matrix)
          Nfence(matrixc)
          j = 1
          do i = 1,nbas
            ifO suscep(j,:) = 0
            j = j + i
          enddo
          ifO suscepc(:,1,:) = 0
          ifO suscepc(1,:,:) = 0
              wing(:,1,:)    = 0
              wingc(:,1,:,:) = 0
          Nfence(matrix)
          Nfence(matrixc)
        endif

c
c       Calculate dielectric matrices (Drude term is added later)
        allocate ( sqcoul(nbas) )
        sqcoul =  sqrt(coul(:nbas))
        dielec => matrix
        do i = 1,nfreq ; Ncycle(i-1)
          suscep(:,i) = -suscep(:,i)                !   -       P
          call p_diagmatdiag( sqcoul, suscep(:,i) ) !   - v^1/2 P v^1/2
          call p_plusdiag( 1d0, suscep(:,i) )       ! 1 - v^1/2 P v^1/2 ! if ikpt=1, the head is set to 1 to enable inversion
        enddo
        nullify ( suscep )
        if(ikpt==1) then
          do i = 1,nfreq
            head(:,i) = identityp(3) - 4*pi * head(:,i) ! interband term
          enddo
          do i = 1,nbas
            wing(:,i,:) = -sqrt(4*pi) * sqcoul(i) * wing(:,i,:)
          enddo
        endif

        dielecc => matrixc
        do i = 1,nfreqc ; Ncycle(i-1)
          call p_diagmatdiag( -sqcoul, suscepc(:,:,i), sqcoul ) !   - v^1/2 P v^1/2
          call p_plusdiag   ( 1d0, suscepc(:,:,i) )             ! 1 - v^1/2 P v^1/2 ! if ikpt=1, the head is set to 1 to enable inversion
        enddo
        if(ikpt==1) then
          do i = 1,nfreqc
            headc(:,:,i) = identity(3) - 4*pi * headc(:,:,i) ! interband term
          enddo
          if(nfreqc/=0) then
            do i = 1,nbas
              wingc(:,i,:,:) = -sqrt(4*pi) * sqcoul(i) * wingc(:,i,:,:)
            enddo
          endif
        endif
        nullify ( suscepc )

        if(job1%type==J_RPA.and.iand(restart,W_cor)/=0.and.typ<=0) call writecor(2)

 2      if(wrtdielec.or.job1%type==J_DIEL) then
          Nfence(matrix)
          Nfence(matrixc)
          Rbegin
          if(ikpt==1) then
            if(lkptadd) then
              write(6,'(A)') 'Additional k vector present: Head is projected on corresponding direction.'
              vec = matmul(rlat,kptadd)
              vec = vec / sqrt(sum(vec**2))
              do i = 1,nfreq  ; dielec (1,i)   = dotprod(vec,matvec(head   (:,i),vec)) - drude (freq (i)) ; enddo
              do i = 1,nfreqc ; dielecc(1,1,i) = dotprod(vec,matmul(headc(:,:,i),vec)) - drudec(freqc(i)) ; enddo
            else
              write(6,'(A)') 'No additional k vector present: Head is angularly averaged.'
              do i = 1,nfreq  ; dielec (1,i)   = trace(head   (:,i))/3 - drude (freq (i)) ; enddo
              do i = 1,nfreqc ; dielecc(1,1,i) = trace(headc(:,:,i))/3 - drudec(freqc(i)) ; enddo
            endif
          endif
          call writeout(dielec,dielecc,nbas,job1%out,'dielec','Dielectric matrix')
          if(ikpt==1) then
            dielec (1,:)   = 1
            dielecc(1,1,:) = 1
          endif
          Rend
          Nfence(matrix)
          Nfence(matrixc)
        endif

c
c       Update RPA correlation energy (Drude term not taken into account)
        if(job1%type==J_RPA) then
          if(.not.allocated(rpa_int)) then ; allocate ( rpa_int(nfreq) ) ; rpa_int = 0 ; endif
          Nfence(matrix)
          Oif(ikpt==1) dielec(1,:) = [( trace(head(:,i))/3 , i=1,nfreq )]
          Nfence(matrix)
          allocate ( eig(nbas) )
          do i = 1,nfreq ; Pcycle(i,nfreq,.true.)
            beginSplit(i)
            call Mdiagonalize(eig,dielec(:,i))
            ifR rpa_int(i) = rpa_int(i) + sum( log(eig) + 1-eig ) * count(kptp==ikpt) / nkpt
            endSplit
          enddo
          deallocate ( eig )
          Nfence(matrix)
          Oif(ikpt==1) dielec(1,:) = 1
          Nfence(matrix)
          nullify(dielec)
          nullify(dielecc)
          goto 4
        endif
        
c
c       Invert dielectric matrices
        Rwrite(6,'(A'NoA) 'Invert dielectric matrices... '

        call cpu_time(cputime)
        
        do i = 1,nfreq ; Pcycle(i,nfreq,.true.)
          beginSplitNodes(i)
          call Minverse(dielec(:,i) MpiC(win_matrix) )
          endSplitNodes
        enddo

# ifdef MPI
        if(nfreq>0) then
          do i = 1,nfreq ; Pcycle_not(i,nfreq,.true.)
            m = Nsize
            beginSplitNodes(i)            ! if Nsize=m, then this complete instance of dielec(:,i) has not entered
            Oif(Nsize==m) dielec(:,i) = 0 ! Minverse, and we must set dielec(:,i) to zero here.
            endSplitNodes
          enddo
          Nfence(matrix)
          Mpi( Ocall Msum(matrix,comm=Ocomm) )
          Nfence(matrix)          
        endif
# endif

        do i = 1,nfreqc ; Pcycle(i,nfreqc,.true.)
          beginSplitNodes(i)
          call Minverse(dielecc(:,:,i) MpiC(win_matrixc) )
          endSplitNodes
        enddo

# ifdef MPI
        if(nfreqc>0) then
          do i = 1,nfreqc ; Pcycle_not(i,nfreqc,.true.)
            m = Nsize
            beginSplitNodes(i)               ! if m==Nsize, then this complete instance of dielecc(:,:,i) has not entered
            Oif(Nsize==m) dielecc(:,:,i) = 0 ! Minverse, and we must set dielecc(:,:,i) to zero here.
            endSplitNodes
          enddo
          Nfence(matrixc)
          Mpi( Ocall Msum(matrixc,comm=Ocomm) )
          Nfence(matrixc)
        endif
# endif

        Rcall cpu_done(cputime)

c       Head and wings of the inverse
        if(ikpt==1) then
          call cpu_time(cputime)
          Rwrite(6,'(A'NoA) 'Calculate head and wings of the inverse... '
          allocate ( wing2(3,nbas) )
          do i = 1,nfreq
            wing2       = MCONJG(wing(:,:,i))
            wing(:,:,i) = - matmat(wing(:,:,i),MCONJG(dielec(:,i)))
            head(:,i)   = head(:,i) + packmat ( matmat ( wing2,transpose(wing(:,:,i)) ) )
          enddo
          deallocate ( wing2 )
          do i = 1,nfreqc
            wingc(:,:,1,i) = - matmat(wingc(:,:,1,i),transpose(dielecc(:,:,i)))
            headc(:,:,i)   = headc(:,:,i) + matmat ( wingc(:,:,2,i),transpose(wingc(:,:,1,i)) )
            wingc(:,:,2,i) = - matmat(wingc(:,:,2,i),dielecc(:,:,i))
          enddo
          Rcall cpu_done(cputime)
          ldum = .false.
          if     (any(freq ==0)) then ; i = minloc(abs(freq), 1) ; rdum = trace(head   (:,i))/3 ; ldum = .true.
          else if(any(freqc==0)) then ; i = minloc(abs(freqc),1) ; rdum = trace(headc(:,:,i))/3 ; ldum = .true.
          endif
          Rif(ldum) then
            if(metal) then ; write(6,'(A,F11.3)') 'Dielectricity constant (w/o Drude term):',rdum
            else           ; write(6,'(A,F11.3)') 'Dielectricity constant:',rdum
            endif
          endif
        endif

        if(iand(restart,W_cor)/=0.and.typ<=0) call writecor(3)

 3      if(wrtdielec.or.job1%type==J_DIEL) then
          if(ikpt==1) then
            if(lkptadd) then
              Rcall writeout([(dotprod(vec,matvec(head   (:,i),vec))-drude (freq (i)),i=1,nfreq )],
     &                       [(dotprod(vec,matmul(headc(:,:,i),vec))-drudec(freqc(i)),i=1,nfreqc)],1,1,'dielecR',
     &                       'Dielectric constant (LF)')
            else
              Rcall writeout([(trace(head   (:,i))/3-drude (freq (i)),i=1,nfreq )],
     &                       [(trace(headc(:,:,i))/3-drudec(freqc(i)),i=1,nfreqc)],1,1,'dielecR','Dielectric constant (LF)')
            endif
          else
            Rcall writeout(1/dielec(1,:),1/dielecc(1,1,:),1,1,'dielecR','Dielectric matrix (LF)')
          endif
          if(job1%type==J_DIEL) goto 4
        endif

c
c       Calculate renormalized susceptibility if JOB = SUSCEPR
        if(job1%type==J_SUSR) then
          Nfence(matrix)
          Nfence(matrixc)
          if(ikpt==1) then ! avoid division by zero
            coul(1)   = 1
            sqcoul(1) = 1
          endif
          suscep  => matrix
          suscepc => matrixc
          do i = 1,nfreq ; Ncycle(i-1)
            call p_plusdiag   ( -1d0,     dielec(:,i) )              !          eps - 1
            call p_diagmatdiag( 1/sqcoul, dielec(:,i) )              ! v^-1/2 ( eps - 1 ) v^-1/2
          enddo
          do i = 1,nfreqc ; Ncycle(i-1)
            call p_plusdiag   ( -1d0,     dielecc(:,:,i) )           !          eps - 1
            call p_diagmatdiag( 1/sqcoul, dielecc(:,:,i), 1/sqcoul ) ! v^-1/2 ( eps - 1 ) v^-1/2            
          enddo
          nullify(dielec)
          nullify(dielecc)
          if(ikpt==1) then
            coul(1)   = 0
            sqcoul(1) = 0
          endif
          Nfence(matrix)
          Nfence(matrixc)
          Rcall writeout(suscep,suscepc,nbas,job1%out,'suscepR','Renormalized susceptibility matrix')
          Nfence(matrix)
          Nfence(matrixc)
          nullify(suscep)
          nullify(suscepc)
          goto 4
        endif

c
c       Replace plasma frequency by PLASMA argument
        Rif(ikpt==1.and.metal) then
          call getkey(inp,'PLASMA', line, section='SUSCEP',  default=' ')
          if(line=='METAL')  then ! use metallic screening (scaling of head of W), flag: plasma = -1
            plasma = -1d0
            write(6,'(A)') 'Plasma frequency set to zero. Metallic screening enforced.'
          else if(line/=' ') then ! replace plasma frequency by argument
            plasma = ch2r(line,'*PLASMA')
            if(plasma<0) Error('Argument to PLASMA must not be negative.')
            write(6,'(A,F10.5,A)') 'Plasma frequency set to',plasma*hartree,' eV.'
          endif
        endif
        Mpi( call Mcast(plasma) )

c
c       Add Drude term
        if(ikpt==1.and.metal.and.plasma>0) then
          Rif(nfreq>=2.and.job1%type==J_GW) call test_plasma
          do i = 1,nfreq
            k = 0
            do j = 1,3
              k         = k + j
              head(k,i) = head(k,i) - drude(freq(i))
            enddo
          enddo
          do i = 1,nfreqc
            do j = 1,3
              headc(j,j,i) = headc(j,j,i) - drudec(freqc(i))
            enddo
          enddo
        endif

        Nfence(matrix)
        Nfence(matrixc)

c
c       Calculate screened interaction
        screen => matrix
        do i = 1,nfreq ; Ncycle(i-1)
          call p_diagmatdiag( sqcoul, dielec(:,i) )            ! v^1/2 eps v^1/2
        enddo
        nullify ( dielec )

        screenc => matrixc
        do i = 1,nfreqc ; Ncycle(i-1)
          call p_diagmatdiag( sqcoul, dielecc(:,:,i), sqcoul ) ! v^1/2 eps v^1/2
        enddo

        nullify ( dielecc )
        if(ikpt==1) then
          do i = 1,nbas
            wing(:,i,:) = sqrt(4*pi) * sqcoul(i) * wing(:,i,:)
          enddo
          if(nfreqc>0) then
            do i = 1,nbas
              wingc(:,i,:,:) = sqrt(4*pi) * sqcoul(i) * wingc(:,i,:,:)
            enddo
          endif
        endif

        if(wrtscreen.or.job1%type==J_SCR) then
          Nfence(matrix)
          Nfence(matrixc)
          Rif(ikpt==1) then
            screen (1,:)   = [( (trace(head   (:,i))/3)**(-1) ,i=1,nfreq)]
            screenc(1,1,:) = [( (trace(headc(:,:,i))/3)**(-1) ,i=1,nfreqc)]
          endif
          Rcall writeout(screen,screenc,nbas,job1%out,'screen','Screened potential')
          Rif(ikpt==1) then
            screen (1,:)   = 0
            screenc(1,1,:) = 0
          endif
          Nfence(matrix)
          Nfence(matrixc)
          if(job1%type==J_SCR) goto 4
        endif

        if(job1%type==J_SCRW) then
          Nfence(matrix)
          Rwrite(6,'(A'NoA) 'Add contribution to Wannier-projected bare and screened interaction... '
          call cpu_time(cputime)
          call coulomb_wannier(screenw,barew,screenc,headc,ctrafo,nbas,ikpt,nfreqc,rsite,nsite,spinw)
          Rif(iand(restart,W_wcou)/=0) call write_dump(3,[ikpt],1,'spex.wcou',job1 MpiC(Mcolor) MpiC(Maxcolor) )
          Rcall cpu_done(cputime)
          Nfence(matrix)
          goto 4
        endif

c        Rif(job1%out==-1) Error('matrices not in eigenspace of Coulomb.') ! disabled because it is obsolete (?)

# ifdef gammatest
        if(optm/=2) Error('use OPTIMIZE PW ... for gammatest.')
        Mpi( if(Mcolor/=0) Error('do not use MPIKPT for gammatest.') )
        do k = 1,nfreq
          screen(:,k) = 0
          j           = 0
          do i = 1,nbas
            vec         = matmul ( rlat , kpt(:,ikpt) + gptm(:,pgptm(gpw(i),ikpt)) )
            rdum        = dot_product(vec,matvec(head(:,k),vec))
            j           = j + i
            screen(j,k) = 4*pi/rdum
          enddo
          if(ikpt==1) screen(1,k) = 0
        enddo
        wing = 0
        do k = 1,nfreqc
          screenc(:,:,k) = 0
          do i = 1,nbas
            vec            = matmul ( rlat , kpt(:,ikpt) + gptm(:,pgptm(gpw(i),ikpt)) )
            cdum           = dot_product(vec,matvec(headc(:,:,k),vec))
            screenc(i,i,k) = 4*pi/cdum
          enddo
          if(ikpt==1) screenc(1,1,k) = 0
        enddo
        if(nfreqc>0) wingc = 0
#   ifndef noWARN
#     warning gammatest
#   endif
# endif

c
c       "Metallic" screening (PLASMA METAL)
        if(metal.and.plasma==-1d0) then
          plasma = 0
          head0  = head(:,1)
          head1  = unpackmat(head0)
          do i = 1,nfreq  ; if(freq(i)==0)  then ; head(:,i) = packmat(identity(3) * infinity()) ; cycle ; endif
            mat          = matmat ( head(:,i)    , matmat ( invert(head0-head(:,i))    , head1-identity(3) ) )
            head(:,i)    = packmat(mat+transpose(mat))/2
          enddo
          do i = 1,nfreqc ; if(freqc(i)==0) then ; headc(:,:,i) =      identity(3) * infinity()  ; cycle ; endif
            mat          = matmat ( headc(:,:,i) , matmat ( invert(head1-headc(:,:,i)) , head1-identity(3) ) )
            headc(:,:,i) =        (mat+transpose(mat))/2
          enddo
          Rwrite(6,'(A,F10.5)') 'Metallic screening enforced by scaling of W(head):',dble(1/(1-3/trace(head0)))
        endif

c
c       Preparations for Gamma point
        if(ikpt==1) then
          ! wrtinfo: Write spex.diel
          Rif(wrtinfo) then
            iunit = fopen('spex.diel',status='unknown')
            do i = 1,nfreq
              write(iunit,'(A,F12.8)') 'Dielectric tensor for w = i *',freq(i)
              write(iunit,ifInv('(3F13.7)','(6F13.7)')) transpose(unpackmat(head(:,i)))
              write(iunit,*)
            enddo
            do i = 1,nfreqc
              write(iunit,'(A,F12.8)') 'Dielectric tensor for w =',freqc(i)
              write(iunit,'(6F13.7)') transpose(headc(:,:,i))
              write(iunit,*)
            enddo
            call fclose(iunit)
            write(6,'(A)') 'Dielectric tensors written to spex.diel.'
          endif
          ! Invert dielectric tensor. Add k averaged contribution to body of screen.
          call cpu_time(cputime)
          Rwrite(6,'(A'NoA) 'Invert dielectric tensor... '
          Nfence(matrix)
          do i = Mrange(1,nfreq)
            h0   = unpackmat(head(:,i)) ; call invert_angular(h0)
            cdum = trace(h0)
            m    = 0
            do j = 1,nbas
              do k = 1,j
                m           = m + 1
                cdum1       = dot_product ( wing(:,j,i) , wing(:,k,i) )                ! trace ( wingwing )
                cdum2       = dot_product ( wing(:,j,i) , matmul( h0 , wing(:,k,i) ) ) ! trace ( h0 * wingwing )
                screen(m,i) = screen(m,i) + ( cdum*cdum1 + 2*cdum2 ) / 15
              enddo
            enddo
          enddo
# ifdef MPI
          Nfence(matrix)
          do j = 0,Msize-1
            Ocall Mcast(screen(:,j*nfreq/Msize+1:(j+1)*nfreq/Msize),Opnt(j),comm=Ocomm)
          enddo
# endif
          if(nfreqc>0) then
            Nfence(matrixc)
            do i = Mrange(1,nfreqc)
              h0   = headc(:,:,i) ; call invert_angular(h0)
              cdum = trace(h0)
              do j = 1,nbas
                do k = 1,nbas
                  cdum1          = sum ( wingc(:,j,2,i) * wingc(:,k,1,i) )
                  cdum2          = sum ( wingc(:,j,2,i) * matmul( h0 , wingc(:,k,1,i) ) )
                  screenc(k,j,i) = screenc(k,j,i) + ( cdum*cdum1 + 2*cdum2 ) / 15
                enddo
              enddo
            enddo
# ifdef MPI
            Nfence(matrixc)
            do j = 0,Msize-1
              Ocall Mcast(screenc(:,:,j*nfreqc/Msize+1:(j+1)*nfreqc/Msize),Opnt(j),comm=Ocomm)
            enddo
# endif
          endif
          Rcall cpu_done(cputime)
        endif

c
c       Update contracted screened interaction (for CORES, COSX, IBC)
        Rif(eval(3)) then
          call screen_contract
          if(iand(restart,W_ccou)/=0) call write_dump(4,[ikpt],1,'spex.ccou',job1 MpiC(Mcolor) MpiC(Maxcolor) )
        endif

c
c       Update self-energy
        Nfence(matrix)
        Nfence(matrixc)
        if(any(eval(:2))) call selfenergy(job1,ikpt,eval,nbas,coul,Win(ctrafo),Win(matrix), head, wing,
     &                                                                         Win(matrixc),headc,wingc,plasma)

# ifdef CHECK_SENERGY
#   ifndef noWARN
#     warning CHECK_SENERGY defined
#   endif
        if(eval(1)) then ; Rwrite(*,*) selfx ; Rread(*,*) ; Mpi(call mpi_barrier(Mcomm,i);) selfx = 0 ; endif
        if(eval(2)) then ; Rwrite(*,*) selfc ; Rread(*,*) ; Mpi(call mpi_barrier(Mcomm,i);) Nfence(selfc);ifO selfc=0;Nfence(selfc)
        endif
# endif

        Rif(iand(restart,W_sigx)/=0.and.eval(1)) call write_dump(1,[ikpt],1,'spex.sigx',job1 MpiC(Mcolor) MpiC(Maxcolor) )
        Rif(iand(restart,W_sigc)/=0.and.eval(2)) call write_dump(2,[ikpt],1,'spex.sigc',job1 MpiC(Mcolor) MpiC(Maxcolor) )

        nullify ( screen  )
        nullify ( screenc )

 4      Ndeallocate( matrix  )
        Ndeallocate( matrixc )
# ifndef gammatest
        deallocate ( head,wing   )
        deallocate ( headc,wingc )
# endif
        deallocate ( coul )
        Ndeallocate( ctrafo )
        Ndeallocate( coulomb0 )
        if(allocated(sqcoul)) deallocate (sqcoul)
        if(allocated(gpw))    deallocate (gpw)
 10     if(allocated(cblock)) deallocate (cblock)
        Mpi( call cpu_time(time_idle) ; Rif(Mcolor/=0) call Mwrite_n(1,IKPT,Mnkpt) )

# ifdef BREAKPNT
#   warning BREAKPNT set
        if(option==ikpt) Error('breakpoint reached.')
# endif

      enddo ! ikpt1 loop

      Mpi( Rif(Mcolor/=0) call Mwrite_n(2, MPI_order ,Mnkpt) )      

c
c     MPIKPT: Finish node-parallel run
# ifdef MPI
      if(Mcolor/=0) then
        Oif(Maxcolor<size(kpt_list)) call mpi_file_close(Mfile,Merr)
        call end_split_nodes     ; time_maxidle = time_idle
        call cpu_time(time_idle) ; time_idle    = time_idle - time_maxidle
        call mpi_reduce(time_idle,time_maxidle,1,mpi_real,mpi_max,0,Mcomm,Merr)
        Rwrite(6,'(/A,//A)') '=========','MPI idle times: '//chr(time_idle,'F15.5')//' '//chr(time_maxidle,'F15.5')        
      endif
# endif

c
c     Contribution from contracted W
      if(allocated(screenk_mt)) then
        Mpi (                           call Msum(screenk_mt) )
        Mpi ( if(allocated(screenk_pw)) call Msum(screenk_pw) )
        eval = .not.lkpt(nkpti+1,:)
        call selfenergy_contract(eval(:2),job1)
        Rif(iand(restart,W_sigx)/=0.and.eval(1)) call write_dump(1,[nkpti+1],1,'spex.sigx',job1 MpiC(Mcolor) MpiC(Maxcolor) )
        Rif(iand(restart,W_sigc)/=0.and.eval(2)) call write_dump(2,[nkpti+1],1,'spex.sigc',job1 MpiC(Mcolor) MpiC(Maxcolor) )
      endif

c
c     Contribution from core exchange
      if(allocated(selfx)) then
        if(.not.lkpt(nkpti+2,1)) then
          call exchange_core(job1)
          Rif(iand(restart,W_sigx)/=0) call write_dump(1,[nkpti+2],1,'spex.sigx',job1 MpiC(Mcolor) MpiC(Maxcolor) )
        endif
      endif

c
c     MPI: Sum all contributions
# ifdef MPI
      if(allocated(selfx))          call Msum(selfx)
      if(job1%type==J_GW)   then ; Ocall Msum(selfc,comm=Ocomm) ; Nfence(selfc) ; endif
      if(job1%type==J_RPA)  then ;  call Msum(rpa_int)                          ; endif
      if(job1%type==J_SCRW) then ;  call Msum(barew) ; call Msum(screenw)       ; endif
# endif

c
c     Write Coulomb matrix
      Rif(job1%type==J_SCRW) then
        call writecou(nkpti)
        call writecou_eff
      endif

      if(allocated(screenw))   tDeallocate ( screenw )
      if(allocated(barew))     tDeallocate ( barew )
      if(allocated(rsite))      deallocate ( rsite )
      if(allocated(freqc))      deallocate ( freqc )
      if(allocated(screenk_mt)) deallocate ( screenk_mt )
      if(allocated(screenk_pw)) deallocate ( screenk_pw )
      deallocate(lkpt)

      contains

c     ------------

      ! Drude function
      ! (1) for imaginary frequencies
      function drude(frq)
      implicit none
      real_dp, intent(in) :: frq
      real_dp             :: drude
      drude = 0
      if(plasma>0) then
        if(frq==0) then ; drude = infinity()
        else            ; drude = -plasma**2 / frq**2
        endif
      endif
      end function drude
      ! (2) for complex frequencies
      function drudec(frq)
      implicit none
      complex_dp, intent(in) :: frq
      complex_dp             :: drudec
      drudec = 0
      if(plasma>0) then
        if(frq==0) then ; drudec = infinity()
        else            ; drudec = plasma**2 / frq**2
        endif
      endif
      end function drudec

c     ------------

c
c     Returns infinity.
      function infinity()
# ifdef noARITHM
      real_dp :: infinity
      infinity = huge(0d0)
# else
      use, intrinsic :: ieee_arithmetic
      implicit none
      real_dp :: infinity
      if(ieee_support_inf(infinity)) then
        infinity = ieee_value(infinity,ieee_positive_inf)
      else
        infinity = huge(0d0)
      endif
# endif
      end function infinity

c     ------------

c
c     Test if the frequency meshes are sufficiently dense. If not, print a warning.
      subroutine test_plasma
      implicit none
      real_dp :: rdum
      if(plasma==-1d0) return
      rdum = plasma/sqrt(dble(head(1,1)+head(3,1)+head(6,1))/3)
      if(rdum<freq(2)) then
        Warn('Small plasma frequency.')
        write(0,'(A)') '             Drude term will produce large variations close to w=0.'
        write(0,'(A)') '             It is recommended to use additional frequency points there.'
        write(0,'(A)') '             Or use PLASMA 0 (or METAL) as an alternative.'
        write(6,'(A)') '--- Copy the following into a file and use MESH "<file>" ---'
        write(6,'(A,F12.8,A)') '# Frequency mesh modified for small plasma frequency',plasma*hartree,' eV'
        write(6,'(F22.15)') 0d0
        rdum = rdum / 2
        if(   rdum<freq(2)) write(6,'(F22.15,A)')    rdum,' # additional point'
        if( 2*rdum<freq(2)) write(6,'(F22.15,A)')  2*rdum,' # additional point'
        if( 4*rdum<freq(2)) write(6,'(F22.15,A)')  4*rdum,' # additional point'
        if( 6*rdum<freq(2)) write(6,'(F22.15,A)')  6*rdum,' # additional point'
        if( 8*rdum<freq(2)) write(6,'(F22.15,A)')  8*rdum,' # additional point'
        if(12*rdum<freq(2)) write(6,'(F22.15,A)') 12*rdum,' # additional point'
        if(16*rdum<freq(2)) write(6,'(F22.15,A)') 16*rdum,' # additional point'
        if(20*rdum<freq(2)) write(6,'(F22.15,A)') 20*rdum,' # additional point'
        if(30*rdum<freq(2)) write(6,'(F22.15,A)') 30*rdum,' # additional point'
        if(40*rdum<freq(2)) write(6,'(F22.15,A)') 40*rdum,' # additional point'
        write(6,'(F22.15)') freq(2:)
        write(6,'(A)') '----------------------------------------------------------'
      endif
      if(nfreqc>1) then
        if(rdum<4*contour_par(4)) then
          Warn('W mesh possibly too coarse for Drude pole.')
          write(0,'(A)')         '             Make W mesh finer or switch off Drude term (with PLASMA 0/METAL) '//
     &                                        'or use Pade CONTOUR or CONTINUE.'
        endif
      endif
      end subroutine test_plasma

c     ------------

      subroutine read_rsite
      implicit none
      character(:),  allocatable :: charr(:)
      character(10)              :: parse_out(3)
      character(4)               :: parse1_out(2)
      real_dp                    :: vec(3),rdum
      logical                    :: ldum
      integer                    :: i,j,m,n
      integer                    :: nn1(3),nn2(3),n1,n2,n3
      integer                    :: ch2i
      Rbegin
      call getkey(inp,'RSITE', charr, section='WANNIER', default=['(0,0,0)'] )
      n    = size(charr) ; if(n==0) Error('Arguments missing after RSITE.')
      ldum = .false.
      do
        nsite = 0
        do i = 1,n
          m = index(charr(i),')-(')
          if(m/=0) then
            call parser(parse_out,charr(i)(:m),  '(1,2,3)','RSITE') ; do j = 1,3 ; nn1(j) = ch2i(parse_out(j),'RSITE') ; enddo
            call parser(parse_out,charr(i)(m+2:),'(1,2,3)','RSITE') ; do j = 1,3 ; nn2(j) = ch2i(parse_out(j),'RSITE') ; enddo
          else
            call parser(parse_out,charr(i),'(1,2,3)','RSITE')
            do j = 1,3
              if(index(parse_out(j),':')/=0) then
                call parser(parse1_out,parse_out(j),'1:2','RSITE')
                nn1(j) = ch2i(parse1_out(1),'RSITE')
                nn2(j) = ch2i(parse1_out(2),'RSITE')
              else
                nn1(j) = ch2i(parse_out(j),'RSITE')
                nn2(j) = nn1(j)
              endif
            enddo
          endif
          vec = nn2-nn1 ; rdum = sqrt(1d0*sum((nn2-nn1)**2)) ; if(rdum>1d-10) vec = vec / rdum
          do n3 = nn1(3),nn2(3),sign(1,nn2(3)-nn1(3))
            do n2 = nn1(2),nn2(2),sign(1,nn2(2)-nn1(2))
              do n1 = nn1(1),nn2(1),sign(1,nn2(1)-nn1(1))
                if(m/=0.and.maxval(abs([n1,n2,n3]-nn1-vec*sqrt(1d0*sum(([n1,n2,n3]-nn1)**2))))>1d-10) cycle
                nsite = nsite + 1 ; if(ldum) rsite(:,nsite) = [n1,n2,n3]
              enddo
            enddo
          enddo
        enddo
        if(ldum) exit
        ldum = .true.
        allocate ( rsite(3,nsite) )
      enddo
      deallocate ( charr )
      if(any(rsite(:,1)/=0)) then
        nsite       = nsite + 1 ; call reallocate(rsite,3,nsite)
        rsite(:,2:) = rsite(:,:nsite-1)
        rsite(:,1)  = 0
      endif
      Rend
      Mpi( call Mcast(nsite); call Mcastl(rsite) )
      end subroutine read_rsite

c     ------------

c     Write Coulomb matrices (barew and screenw) to spex.cou
      subroutine writecou(ikpt)
      implicit none
      integer,     intent(in) :: ikpt
      integer                 :: isite
      integer                 :: iunit
      integer                 :: i,j,m,mm,n
      complex_dp              :: cdum
      character(3)            :: cjob
      character(5), parameter :: spin_ch(0:4) = [ 'uu+dd','uu   ','dd   ','ud   ','du   ' ]
      cjob  = ' ' ; if(maxval(job%indx)>1) cjob = chr(job1%indx)
      iunit = fopen('spex'//trim(cjob)//'.cou',form='formatted',status='unknown')
      write(iunit,'(A,I3,2X,A)') 'Spin index:',spinw,spin_ch(spinw)
      write(iunit,'(A,I3)')      'Number of Wannier functions:',nwan
      write(iunit,'(A,I15)')     'Number of sites:',            nsite
      write(iunit,'(A,I9)')      'Number of frequencies:',      nfreqc
      write(iunit,'(A,I13)')     'Completed k point:',          ikpt
      write(iunit,'(/A)') 'Wannier-projected bare interaction (in eV):'
      do isite = 1,nsite
        write(iunit,'(7A)') 'Site (',trim(chr(rsite(1,isite))),',',trim(chr(rsite(2,isite))),',',trim(chr(rsite(3,isite))),')'
        do i = 1,nwan ; do j = 1,nwan ; do m = 1,nwan ; do n = 1,nwan
          cdum = barew(i,j,m,n,isite) * hartree
          if(abs(cdum)>1d-10) write(iunit,'(4I4,2F17.10)') i,j,m,n,cdum
        enddo ; enddo ; enddo ; enddo
      enddo
      write(iunit,'(/A)') 'Wannier-projected screened interaction (in eV):'
      do isite = 1,nsite
        write(iunit,'(7A)') 'Site (',trim(chr(rsite(1,isite))),',',trim(chr(rsite(2,isite))),',',trim(chr(rsite(3,isite))),')'
        do mm = 1,nfreqc ; do i = 1,nwan ; do j = 1,nwan ; do m = 1,nwan ; do n = 1,nwan
          cdum = screenw(i,j,m,n,mm,isite) * hartree
          if(abs(cdum)>1d-10) write(iunit,'(2F15.10,4I4,2F20.10)') freqc(mm)*hartree,i,j,m,n,cdum
        enddo ; enddo ; enddo ; enddo ; enddo
      enddo
      write(iunit,*)
      call fclose(iunit)
      write(6,'(/A)') '---------'
      write(6,'(/A)') 'Full matrix written to spex'//trim(cjob)//'.cou.'
      end subroutine writecou

c     ------------

c     Effective Coulomb parameters written to output
      subroutine writecou_eff
      implicit none
# include "interface/orbitals.inc"
      character(:), allocatable :: charr(:)
      real_dp,      allocatable :: eulerwan(:,:)
      integer,      allocatable :: lwan(:),lmwan(:),centwan(:)
      integer                   :: i,j,m,mm,n
      complex_dp                :: cdum,cdum1
      barew   = barew   * hartree
      screenw = screenw * hartree
      freqc   = freqc   * hartree
      ! (a) get orbital characters from ORBITALS
      call getkey(inp,'ORBITALS',charr,section='WANNIER')
      call orbitals(lwan,lmwan,centwan,eulerwan,i,charr(3:),size(charr)-2)
      deallocate ( charr )
      centwan = abs(centwan)
      if(l_soc) then
        call reallocate(lwan,      2*nwan) ; lwan      (nwan+1:2*nwan) = lwan      (:nwan)
        call reallocate(lmwan,     2*nwan) ; lmwan     (nwan+1:2*nwan) = lmwan     (:nwan)
        call reallocate(centwan,   2*nwan) ; centwan   (nwan+1:2*nwan) = centwan   (:nwan)        
        call reallocate(eulerwan,3,2*nwan) ; eulerwan(:,nwan+1:2*nwan) = eulerwan(:,:nwan)
        nwan = nwan * 2
      endif
      call getkey(inp,'SUBSET',charr,section='WANNIER',status=j)
      if(j/=0) then
        j              = len_trim(charr(1)) - 1 ; call str2iarr(iarr,charr(1)(2:j),'SUBSET')
        i              = size(iarr)
        lwan(:i)       = lwan(iarr)             ; call reallocate(lwan,i)
        lmwan(:i)      = lmwan(iarr)            ; call reallocate(lmwan,i)
        centwan(:i)    = centwan(iarr)          ; call reallocate(centwan,i)
        eulerwan(:,:i) = eulerwan(:,iarr)       ; call reallocate(eulerwan,3,i)
        deallocate ( charr,iarr )
      endif
      if(i/=nwan) Bug('Count error.')
      ! (b) Effective parameters
      if(all(lwan>=100.and.lwan/=202.and.lwan/=302)) then
        write(6,'(/A)') 'No effective parameters given due to absence of full atomic shells in ORBITALS.'
      else
        write(6,'(/A'NoA) 'Effective Coulomb and exchange parameters for'
        if(any(lwan<100))                                        write(6,'(A'NoA) ' full shells'
        if(any(lwan<100).and.(any(lwan==202).or.any(lwan==302))) write(6,'(A'NoA) ' and'
        if(                   any(lwan==202).or.any(lwan==302))  write(6,'(A'NoA) ' t2g/eg orbitals'
        write(6,'(/A)')   '(Orbital character according to ORBITALS def, actual character could differ)'
        m = 1
        do while(m<=nwan)
          if(lwan(m)<100) then ; n = 2*lwan(m) + 1
          else                 ; n = lwan(m) / 100
          endif
          mm = m + n - 1
          if(lwan(m)<100) then
            write(6,'(/A,I3'NoA) 'Atom',centwan(m)
            if(lwan(m)==0) then ; write(6,'(A,I2)') '  --  Hubbard parameter for l =',lwan(m)
            else                ; write(6,'(A,I2)') '  --  Hubbard-Hund parameters for l =',lwan(m)
            endif
            cdum = sum( [ ((barew(i,j,i,j,1), i=m,mm),j=m,mm) ] ) / n**2
            write(6,'(A,F10.5'NoA) '   V  =',real(cdum)
            if(n>1) then
              cdum = cdum - sum( [ ((barew(i,j,i,j,1)-barew(i,j,j,i,1), i=m,mm),j=m,mm) ] ) / (n*(n-1))
              if(nwan>1) write(6,'(A,F10.5'NoA) ',  Jb =',real(cdum)
            endif
            write(6,*)
            cdum = sum( [ ((screenw(i,j,i,j,1,1),i=m,mm),j=m,mm) ] ) / n**2
            write(6,'(A,F10.5'NoA) '   U  =',real(cdum)
            if(n>1) then
              cdum = cdum - sum( [ ((screenw(i,j,i,j,1,1)-screenw(i,j,j,i,1,1),i=m,mm),j=m,mm) ] ) / (n*(n-1))
              if(nwan>1) write(6,'(A,F10.5'NoA) ',  J  =',real(cdum)
            endif
            write(6,*)
          else if(any(lwan(m)==[202,302])) then
            write(6,'(/A,I3'NoA) 'Atom',centwan(m)
            if(lwan(m)==302) then ; write(6,'(A)') '  --  Kanamori parameters for t2g orbitals'
            else                  ; write(6,'(A)') '  --  Kanamori parameters for eg orbitals'
            endif
            cdum = sum( [ (barew(i,i,i,i,1), i=m,mm) ] ) / n ; cdum1 = cdum
            write(6,'(A,F10.5'NoA) '   V  =',real(cdum)
            cdum = sum( [ ((barew(i,j,i,j,1), i=m,mm),j=m,mm) ] ) / (n*(n-1)) - cdum1 / (n-1)
            write(6,'(A,F10.5'NoA) ',  V'' =',real(cdum)
            if(n>1) then
              cdum = sum( [ ((barew(i,j,j,i,1), i=m,mm),j=m,mm) ] ) / (n*(n-1)) - cdum1 / (n-1)
              write(6,'(A,F10.5'NoA) ',  Jb =',real(cdum)
            endif
            write(6,*)
            cdum = sum( [ (screenw(i,i,i,i,1,1), i=m,mm) ] ) / n ; cdum1 = cdum
            write(6,'(A,F10.5'NoA) '   U  =',real(cdum)
            cdum = sum( [ ((screenw(i,j,i,j,1,1), i=m,mm),j=m,mm) ] ) / (n*(n-1)) - cdum1 / (n-1)
            write(6,'(A,F10.5'NoA) ',  U'' =',real(cdum)
            if(n>1) then
              cdum = sum( [ ((screenw(i,j,j,i,1,1), i=m,mm),j=m,mm) ] ) / (n*(n-1)) - cdum1 / (n-1)
              write(6,'(A,F10.5'NoA) ',  J  =',real(cdum)
            endif
            write(6,*)
          endif
          m = mm + 1
        enddo
      endif
      deallocate ( lwan,lmwan,centwan,eulerwan )
      barew   = barew   / hartree
      screenw = screenw / hartree
      freqc   = freqc   / hartree
      end subroutine writecou_eff

c     ------------

c
c     Determines basis set - defines n, ctrafo, coul, trafo, symon, (gpw)

      subroutine get_basis
      implicit none
      MCOMPLEX_dp, pointer_cnt :: olap(:,:)
      integer                  :: leng,i
      Mpi( integer             :: win_olap )

      trafo = optm/=0

      if(optm<=1) then

        ! Use eigenvectors of Coulomb matrix as basis

        nbas  = nbasm(ikpt)
        symon = .true.

        !   allocation
        allocate  ( coul(nbas) )
        Nallocate0( ctrafo, (S_ nbas,nbas S_) )
        Nallocate ( olap,   (S_ ngptm(ikpt),ngptm(ikpt) S_) )

        !   diagonalize Coulomb matrix
        if(optm==1) nbas = max(0,nint(optcutm))
        call diagonalize_coulomb(Win(ctrafo),coul,Win(olap),coulomb0,ikpt,nbas,-optcutm)
        Ndeallocate ( olap )
        if(optm==0.and.nbas/=nbasm(ikpt)) then
          if(noapw) Error('Numbers of Coulomb eigenvectors and MPB functions disagree.')
          trafo = .true.
        endif
        leng = 5+max(2,len_trim(chr(ceiling(coul(nbas)))))
        Rwrite(6,'(A,F'//chr(leng)//'.4,''   (MB'',F6.2,'')'')') 'Minimum Coulomb eigenvalue:',coul(nbas),sqrt(4*pi/coul(nbas))
        call reallocate ( coul,nbas )

      else if(optm==2) then

        ! Use plane waves as basis

        symon = .false. ; RWarn('Symmetry switched off for OPTIMIZE PW.')

        !   allocation
        allocate  ( coul(ngptm(ikpt)) )
        Nallocate0( ctrafo,(S_ nbasm(ikpt),ngptm(ikpt) S_) )

        !   determine G point set
        nbas = max(0,nint(optcutm))
        allocate ( gpw(ngptm(ikpt)) )
        call getgpw(gpw,nbas,-optcutm)

        !   determine transformation coefficients
        call getctrafo(Win(ctrafo),nbas,gpw,ikpt)

        !   calculate diagonal elements of Coulomb matrix (in basis of plane waves)
        do i = 1,nbas
          coul(i) = 4*pi / sum(matmul(rlat,kpt(:,ikpt)+gptm(:,pgptm(gpw(i),ikpt)))**2)
        enddo
        if(ikpt==1) coul(1) = 0

      else
        Error(' optm not known.')
      endif

      end subroutine get_basis

c     ------------

      subroutine getgpw(gpw,nbas,optcutm_r)
      implicit none
      integer, intent(inout) :: nbas
      real_dp, intent(in)    :: optcutm_r
      integer                :: gpw(:)
      real_dp                :: rarr(ngptm(ikpt)),kp(3)
      integer                :: i
      kp = kpt(:,ikpt)
      if(nbas>ngptm(ikpt)) then
        Warn('Number of plane waves exceeds number of G points. Is reduced')
        nbas = ngptm(ikpt)
      endif
      do i = 1,ngptm(ikpt)
        rarr(i) = sqrt(sum(matmul(rlat,kp+gptm(:,pgptm(i,ikpt)))**2))
      enddo
      call rorderp(gpw,rarr,ngptm(ikpt))
      rarr = rarr(gpw)
      if(ngptm(ikpt)==1) then ; nbas = 1 ; return ; endif
      if(nbas==0) then
        nbas = 1
        do while(4*pi/rarr(nbas+1)**2>optcutm_r)
          nbas = nbas + 1 ; if(nbas==ngptm(ikpt)) exit
        enddo
      else
        do while(nbas<ngptm(ikpt))
          if(abs(rarr(nbas+1)-rarr(nbas))>1d-8) exit
          nbas = nbas + 1
        enddo
      endif
      end subroutine getgpw

c     ------------

      subroutine writeout(mat,matc,dim,out,filename,description)
      implicit none
      integer,     intent(in) :: dim,out
      MCOMPLEX_dp, intent(in) :: mat(dim*(dim+1)/2,nfreq)
      complex_dp,  intent(in) :: matc(dim,dim,nfreqc)
      character(4)            :: eunit
      character*(*)           :: filename,description
      real_dp                 :: kp(3),rdum1,rdum2
      integer                 :: iunit,out1,g(3)
      integer                 :: ifreq,i,j,ij,k,ikindx

      if(escale==hartree) then ; eunit = '(eV)'
      else if(escale==1)  then ; eunit = '(Ha)'
      else                     ; eunit = '    '
      endif

      if(nfreq+nfreqc<=0) return
      kp = kpt(:,ikpt)

      if(out==-1) then
        iunit = fopen(filename//trim(cjob)//'.bin','unknown',form='unformatted')
        write(iunit) kp,description
        do ifreq = 1,nfreq
          write(iunit) freq(ifreq)
          write(iunit) mat(:,ifreq)
        enddo
        do ifreq = 1,nfreqc
          write(iunit) freqc(ifreq)
          write(iunit) matc(:,:,ifreq)
        enddo
        call fclose(iunit)
        write(6,'(A)') description//' written to binary file '//filename//trim(cjob)//'.bin .'
        return
      endif

      out1 = out
      if(out>dim) then
        out1 = dim
        Warn('Submatrix for output cannot be larger than total matrix.')
        write(0,'(2(A,I3))') '          Submatrix is reduced to',out1,' x',out1
      endif

      if(ikpt<=nkpt) then ; ikindx = ikpt
      else                ; ikindx = 0 ; do k = 1,nkpt ; if(sum(abs(kpt(:,ikpt)-kpt(:,k)))<1d-8) ikindx = -k ; enddo
      endif
      iunit = fopen(filename//trim(cjob),'unknown')
      write(iunit,'(A)')        '# Projected '//description
      write(iunit,'(A,9F10.5)') '# lattvec:',lat
      write(iunit,'(A,3F10.5)') '# k point:',kpt(:,ikpt)
      write(iunit,'(A,I6'NoA)   '# k index:',abs(ikindx) ; if(ikindx<=0) write(iunit,'(A'NoA) '  (additional)'
      write(iunit,'(/A,I2)')    '# spin:',spin
      do j = 1,out1
        do i = 1,out1
          ij = j*(j-1)/2 + i
          if(nfreq>0.and.i<=j.or.nfreqc>0) then
            if(out1>1) write(iunit,'(/A,2I5/)') '#',i,j
            if(allocated(gpw)) then ;                                                     g=gptm(:,pgptm(gpw(i),ikpt))
              write(iunit,'(A,3I3,F10.5)') '# k+G1 =',g,sqrt(sum(matmul(rlat,kp+g)**2)) ; g=gptm(:,pgptm(gpw(j),ikpt))
              write(iunit,'(A,3I3,F10.5)') '# k+G2 =',g,sqrt(sum(matmul(rlat,kp+g)**2))
            else
              if(coul(i)==0) then ; rdum1 = 0 ; else ; rdum1 = sqrt(4*pi)/coul(i) ; endif
              if(coul(j)==0) then ; rdum2 = 0 ; else ; rdum2 = sqrt(4*pi)/coul(j) ; endif
              write(iunit,'(A,F10.5)')     '# |k+G1| =',rdum1
              write(iunit,'(A,F10.5)')     '# |k+G2| =',rdum2
            endif
          endif
          if(nfreq>0.and.i<=j) then
            write(iunit,'(/A)') '#   im. Frequency '//eunit//'      Real Part           ' NoInv(//'Imaginary Part')
            do ifreq = 1,nfreq
              write(iunit,'(3F20.10)') freq(ifreq)*escale,dble(mat(ij,ifreq)) NoInvC(imag(mat(ij,ifreq)))
            enddo
          endif
          if(nfreqc>0) then
            if(all(imag(freqc)==0)) then
              write(iunit,'(/A)') '#       Frequency '//eunit//'      Real Part           Imaginary Part'
              do ifreq = 1,nfreqc
                write(iunit,'(3F20.10)') real(freqc(ifreq))*escale,matc(i,j,ifreq)
              enddo
            else
              write(iunit,'(/A)') '#       Frequency '//eunit//'      Real Part           Imaginary Part'
              do ifreq = 1,nfreqc
                write(iunit,'(4F20.10)') freqc(ifreq)*escale,matc(i,j,ifreq)
              enddo
            endif
          endif
        enddo
      enddo
      call fclose(iunit)
      write(6,'(A)') description//' written to '//filename//trim(cjob)//' .'
      end subroutine writeout

c     ------------

c     readcor:  Reads  spex.cor
c     writecor: Writes spex.cor

# define ERRSTOP if(Herr/=0) Error('Fatal HDF5 error.')

      subroutine readcor(typ)
      use Hwrapper
      use hdf5
      implicit none
      integer, intent(inout)  :: typ
      integer(HID_T)          :: Hfile,Hkpt
      integer                 :: Herr
      character(38)           :: Hname,fname
      character(2)            :: cjob
      character(6)            :: label1
      integer                 :: iarr(2),i,j
      integer                 :: iunit,version
      real_dp                 :: rarr(4)
      real_dp,    allocatable :: freq1(:)
      complex_dp, allocatable :: freqc1(:)
      logical                 :: ldum
      logical, save           :: first = .true.
      write(Hname,'(F12.10,'','',F12.10,'','',F12.10)') kpt(:,ikpt)
      i    = 0
      cjob = ' ' !; if(maxval(job%indx)>1) write(cjob,'(I2.2)') job1%indx
      Rbegin
      inquire(file='spex'//trim(cjob)//'.cor.map/present',exist=ldum)
      if(ldum) then
        inquire(file='spex'//trim(cjob)//'.cor.map/'//Hname,exist=ldum)
        if(ldum) then
          iunit = fopen('spex'//trim(cjob)//'.cor.map/'//Hname,status='old')
          read(iunit,*,iostat=j) i
          if(j/=0) Error('Could not read corfile handle.')
          if(i<=0) Bug('corfile handle non-positive.')
          call fclose(iunit)
        endif
      endif
      if(i==0) then ; fname = 'spex'//trim(cjob)//'.cor'
      else          ; fname = 'spex'//trim(cjob)//'.cor.'//chr(i)
      endif
      inquire(file=trim(fname),exist=ldum)
      Rend
      Mpi( call Mcast(ldum) ; call Mcast(fname) )
      if(.not.ldum) then
        typ = -1
        return
      endif
      ! Open file
      call hdf_fopen(Hfile,trim(fname),0)
      ! Read version number
      version = 0 ! default
      call h5aexists_f(Hfile,'version',ldum,Herr) ; ERRSTOP
      if(ldum) call hdf_rdwr_a(Hfile,'version',0,version)
      if(all(version/=[0,1])) Error('COR version number unknown: '//Chr(version))
      ! Open k-point group (return if it does not exist)
      call h5lexists_f(Hfile,Hname,ldum,Herr) ; ERRSTOP
      if(.not.ldum) then
        typ = -1
        call hdf_fclose(Hfile)
        return
      endif
      call h5gopen_f(Hfile,Hname,Hkpt,Herr) ; ERRSTOP
      ! Read coulomb0
      call hdf_rdwr_a(Hkpt,'nbasm',0,i) ; Rif(i/=nbasm(ikpt)) Error('Current basis size differs from the one in '//trim(fname))
      Nallocate0 ( coulomb0,(S_ nbasm(ikpt)*(nbasm(ikpt)+1)/2 S_) )
      call hdf_rdwr(Hkpt,'coulomb',0,coulomb0) ; Nfence(coulomb0)
      Nfence(coulomb0)
      Rwrite(6,'(A)') 'Coulomb matrix read from '//trim(fname)
      if(typ>0) then
        call h5aexists_f(Hkpt,'type',ldum,Herr) ; ERRSTOP
      endif
      if(typ==0.or..not.ldum) then
        typ = 0
        call h5gclose_f(Hkpt,Herr) ; ERRSTOP
        call hdf_fclose(Hfile)
        return
      endif
      ! Read file attributes
      if(first) then
        first = .false.
        call hdf_rdwr_a(Hfile,'nfreq,nfreqc',0,iarr)
        if(nfreq/=0) then
          if(iarr(1)==0) Error('No imaginary frequency mesh in '//trim(fname))
          allocate(freq1(iarr(1)))
          call hdf_rdwr_a(Hfile,'freq',0,freq1)
          ldum = .false.
          if(iarr(1)/=nfreq)        then ; ldum = .true. ; nfreq = iarr(1) ; deallocate(freq) ; allocate(freq(nfreq))
          else if(any(freq1/=freq)) then ; ldum = .true.
          endif
          if(ldum) then
            freq = freq1
            RWarn('Imaginary frequency mesh replaced by the one from '//trim(fname)//': '//)
     &        Chr(freq(1))//' : '//Chr(freq(nfreq))
            Rif(any(job1%type==[J_SCRW,J_RPA,J_SX,J_COSX]).or.job1%type==J_GW.and.oselfc==1)
     &                        Error('Cannot modify mesh for SCREENW, RPAENE, SX/COSX, or GW(CONTINUE).')
            Rif(nfreq==1) then
              if(job1%type==J_GW) then ; Warn('Imaginary frequency mesh contains only single point.')
              else                     ; Warn('Imaginary frequency mesh contains a single point.')
              endif
            endif
          endif
        endif
        if(nfreqc/=0) then
          if(iarr(2)==0) Error('No complex frequency mesh in '//trim(fname))
          allocate(freqc1(iarr(2)))
          call hdf_rdwr_a(Hfile,'freqc',0,freqc1)
          ldum = .false.
          if(iarr(2)/=nfreqc)         then ; ldum = .true. ; nfreqc = iarr(2) ; deallocate(freqc) ; allocate(freqc(nfreqc))
          else if(any(freqc1/=freqc)) then ; ldum = .true.
          endif
          if(ldum) then
            freqc = freqc1
            RWarn('Complex frequency mesh replaced by the one from '//trim(fname)//': '//)
     &        Chr(real(freqc(1)))     //'+img*'//Chr(imag(freqc(1)))//' : '//
     &        Chr(real(freqc(nfreqc)))//'+img*'//Chr(imag(freqc(nfreqc)))
            Rif(any(job1%type==[J_SCRW,J_RPA])) Error('Must not modify mesh for SCREENW or RPAENE.')
            Rif(nfreqc==1) Warn('Complex frequency mesh contains a single point.')
          endif
        endif
      endif
      ! Read k-point attributes
      call hdf_rdwr_a(Hkpt,'type',   0,i) ; Rif(i>typ) Error('Matrix type inconsistent.')
      typ = i
      call hdf_rdwr_a(Hkpt,'label',  0,label1)
      call hdf_rdwr_a(Hkpt,'dim',    0,nbas)
      call hdf_rdwr_a(Hkpt,'contour',0,rarr)
      call hdf_rdwr_a(Hkpt,'oselfc', 0,i)
      Rif(i/=oselfc) then
        if(i==1.and.oselfc>1.and.nfreqc/=0) Bug('oselfc(file)==1 but oselfc>1 and nfreqc/=0.')
        Warn('Data taken from different self-energy type.')
      endif
      Rif(any(rarr/=contour_par).and.oselfc>1) Warn('Contour parameters different.')
      if(ikpt==1) call hdf_rdwr_a(Hkpt,'plasma',0,plasma)
      call hdf_rdwr_a(Hkpt,'ncblock',0,i)
      ! Allocate arrays
      if(i/=0) allocate ( cblock(i) )
      Nallocate0 ( ctrafo,  (S_ nbasm(ikpt),nbas       S_) )
      Nallocate0 ( matrix,  (S_ nbas*(nbas+1)/2,NFREQ  S_) )
      Nallocate0 ( matrixc, (S_ nbas,nbas,NFREQC       S_) )
      allocate   ( coul(nbas) )
      allocate   ( head(6,nfreq)     ) ; allocate ( wing(3,nbas,nfreq) )
      allocate   ( headc(3,3,nfreqc) ) ; allocate ( wingc(3,nbas,2,nfreqc) )
      ! Read data
      call hdf_rdwr(Hkpt,'coul',   0,coul)
      call hdf_rdwr(Hkpt,'ctrafo', 0,ctrafo) ; Nfence(ctrafo)
      if(allocated(cblock)) call hdf_rdwr(Hkpt,'cblock',0,cblock)
      if(version==0) then
        if(nfreq>0) then
          call hdf_rdwr(Hkpt,'matrix',0,matrix)
          Nfence(matrix)
        endif
        if(nfreqc>0) then
          call hdf_rdwr(Hkpt,'matrixc',0,matrixc)
          Nfence(matrixc)
        endif
      else        
        do i = 1,nfreq
          call hdf_rdwr(Hkpt,'matrix_'//Chr(i),0,matrix(:,i))
        enddo        
        do i = 1,nfreqc
          call hdf_rdwr(Hkpt,'matrixc_'//Chr(i),0,matrixc(:,:,i))
        enddo
        Nfence(matrix)
        Nfence(matrixc)
      endif
      if(nfreq>0) then
        if(ikpt==1) then
          call hdf_rdwr(Hkpt,'head',0,head)
          call hdf_rdwr(Hkpt,'wing',0,wing)
        endif
      endif
      if(nfreqc>0) then
        if(ikpt==1) then
          call hdf_rdwr(Hkpt,'headc',0,headc)
          call hdf_rdwr(Hkpt,'wingc',0,wingc)
        endif
      endif
      if(typ==1) then
        suscep => matrix ; suscepc => matrixc
        Rwrite(6,'(A'NoA) 'Susceptibility '
      else if(typ==2.or.typ==3) then
        dielec => matrix ; dielecc => matrixc ; allocate(sqcoul(nbas)) ; sqcoul = sqrt(coul)
        Rif(typ==2) write(6,'(A'NoA) 'Dielectric '
        Rif(typ==3) write(6,'(A'NoA) 'Inverse dielectric '
      else
        RBug('Unknown type index: '//trim(chr(typ))//'.')
      endif
      Rwrite(6,'(A'NoA) 'matrices read from '//trim(fname)
      Rif(label1/=job1%label) write(6,'(A'NoA) ' (job '//trim(label1)//')'
      Rwrite(6,*)
      ! Close all
      call h5gclose_f(Hkpt,Herr) ; ERRSTOP
      call hdf_fclose(Hfile)
      if(ikpt==1.and.oibc/=0) call ibc1_contraction(freq*img,nfreq,.false.)
      end subroutine readcor

      subroutine writecor(typ)
      use Hwrapper
      use hdf5
      implicit none
      integer, intent(in) :: typ
      integer(HID_T)      :: Hfile,Hkpt
      integer             :: Herr
      character(38)       :: Hname,fname
      character(80)       :: err
      character(2)        :: cjob
      logical             :: ldum
      integer             :: version
      integer             :: ncblock MpiC(ios)
      if(typ/=0) then
        Nfence(matrix)
        Nfence(matrixc)
      endif
# ifdef HDF5ser
      MnoR( goto 20 ) ! return if not root (only root writes in case of HDF5ser)
# endif
      cjob  = ' ' !; if(maxval(job%indx)>1) write(cjob,'(I2.2)') job1%indx
      fname = 'spex'//trim(cjob)//'.cor' ; Mpi( if(Mcolor/=0) fname(len_trim(fname)+1:) = '.'//trim(chr(Mcolor)) )
      ! Open file
      call hdf_fopen(Hfile,trim(fname),1)
      ! Write version number
      version = 1
      call hdf_rdwr_a(Hfile,'version',2,version)
      ! Open k-point group
      write(Hname,'(F12.10,'','',F12.10,'','',F12.10)') kpt(:,ikpt)
      call h5lexists_f(Hfile,Hname,ldum,Herr)   ; ERRSTOP
      if(ldum) then
        call h5gopen_f(Hfile,Hname,Hkpt,Herr)   ; ERRSTOP
      else
        call h5gcreate_f(Hfile,Hname,Hkpt,Herr) ; ERRSTOP
      endif
      select case(typ)
        case (0)     ; Rwrite(6,'(A'NoA) 'Coulomb matrix '
        case (1)     ; Rwrite(6,'(A'NoA) 'Susceptibility matrices '
        case (2)     ; Rwrite(6,'(A'NoA) 'Dielectric matrices '
        case (3)     ; Rwrite(6,'(A'NoA) 'Inverse dielectric matrices '
        case default ; RBug('Unknown type index: '//trim(chr(typ))//'.')
      end select
      Rwrite(6,'(A)') 'written to '//trim(fname)
      ! Write coulomb0 (and return if typ==0)
      call hdf_rdwr_a(Hkpt,'nbasm',  2,nbasm(ikpt))
      call hdf_rdwr  (Hkpt,'coulomb',2,coulomb0)
      if(typ==0) goto 10
      ! Write file attributes
      call hdf_rdwr_a(Hfile,'nfreq,nfreqc',2,[nfreq,nfreqc])
      if(nfreq >0) call hdf_rdwr_a(Hfile,'freq', 2,freq)
      if(nfreqc>0) call hdf_rdwr_a(Hfile,'freqc',2,freqc)
      ! Write k-point attributes
      call h5aexists_f(Hkpt,'type',ldum,Herr) ; ERRSTOP
      if(ldum) then
        RInfo('Dataset '//Hname//' already exists in '//trim(fname)//'. Skipping ...')
        call h5gclose_f(Hkpt,Herr) ; ERRSTOP
        call hdf_fclose(Hfile)
        goto 20 ! return
      endif
      ncblock = 0 ; if(allocated(cblock)) ncblock = size(cblock)
      call hdf_rdwr_a(Hkpt,'type',   1,typ)
      call hdf_rdwr_a(Hkpt,'label',  1,job1%label)
      call hdf_rdwr_a(Hkpt,'dim',    1,nbas)
      call hdf_rdwr_a(Hkpt,'ncblock',1,ncblock)
      call hdf_rdwr_a(Hkpt,'oselfc', 1,oselfc)
      call hdf_rdwr_a(Hkpt,'contour',1,contour_par)
      if(ikpt==1) call hdf_rdwr_a(Hkpt,'plasma',1,plasma)
      ! Write data
      call hdf_rdwr(Hkpt,'coul',  1,coul(:nbas))
      call hdf_rdwr(Hkpt,'ctrafo',1,ctrafo(:,:nbas))
      if(ncblock>0) call hdf_rdwr(Hkpt,'cblock',1,cblock)
      do i = 1,nfreq
        call hdf_rdwr(Hkpt,'matrix_'//Chr(i),1,matrix(:,i))
      enddo
      do i = 1,nfreqc
        call hdf_rdwr(Hkpt,'matrixc_'//Chr(i),1,matrixc(:,:,i))
      enddo
      if(nfreq>0) then
        !call hdf_rdwr(Hkpt,'matrix',1,matrix) ! obsolete COR version 0
        if(ikpt==1) then
          call hdf_rdwr(Hkpt,'head',1,head)
          call hdf_rdwr(Hkpt,'wing',1,wing)
        endif
      endif
      if(nfreqc>0) then
        !call hdf_rdwr(Hkpt,'matrixc',1,matrixc) ! obsolete COR version 0
        if(ikpt==1) then
          call hdf_rdwr(Hkpt,'headc',1,headc)
          call hdf_rdwr(Hkpt,'wingc',1,wingc)
        endif
      endif
      ! Close all
 10   call h5gclose_f(Hkpt,Herr) ; ERRSTOP
      call hdf_fclose(Hfile)
      Rbegin
      inquire(file='spex'//trim(cjob)//'.cor.map/present',exist=ldum)
# ifdef MPI
      if(.not.ldum) then
        if(Mcolor>0) then
          call shell_command('mkdir -p spex'//trim(cjob)//'.cor.map',ios)
          if(ios/=0) Error('Could not create directory spex'//trim(cjob)//'.cor.map')
          call shell_command('touch spex'//trim(cjob)//'.cor.map/present',ios)
          if(ios/=0) Error('Could not create empty file in directory spex'//trim(cjob)//'.cor.map')
          do ! make sure that the file just created is seen from Fortran
            inquire(file='spex'//trim(cjob)//'.cor.map/present',exist=ldum)
            if(ldum) exit
            Info('File creation not ready. Rank '//Chr(Mrank0)//' will sleep for 0.1 seconds.')
            call asleep(0.1)
          enddo
        endif
      endif
# endif
      if(ldum) then
        Mpi( if(Mcolor==0) ) inquire(file='spex'//trim(cjob)//'.cor.map/'//Hname,exist=ldum)
        if(ldum) then
          iunit = fopen('spex'//trim(cjob)//'.cor.map/'//Hname,status='unknown')
          Mpi( if(Mcolor>0)  then ; write(iunit,*) Mcolor ; call fclose(iunit) ; endif )
          Mpi( if(Mcolor==0) )                              call fclose(iunit,status='delete')
        endif
      endif
      Rend
 20   if(typ/=0) then
        Nfence(matrix)
        Nfence(matrixc)
      endif
      end subroutine writecor

c     ------------

# ifdef MPI
#   define CB call Mcast
#   define CA call Mcastl
      subroutine Mbroadcast0
      implicit none
      CB(nfreq); CB(nfreqr); CB(nfreqc); CB(oselfc)
      end subroutine

      subroutine Mbroadcast1
      implicit none
      CB(lkpt)
      CA(freq) ; CA(freqr) ; CA(freqc) ; CB(contour_par)
      end subroutine Mbroadcast1

      subroutine Mbroadcast3
      implicit none
      CB(head);  CB(wing)
      CB(headc); CB(wingc)
      end subroutine Mbroadcast3

      subroutine Mbroadcast4
      implicit none
      CB(plasma); CB(nbas); CB(coul); CB(sqcoul); CB(ctrafo)
      CA(cblock); CB(head); CB(wing); CB(headc) ; CB(wingc)
      Obegin
      CB(matrix,comm=Ocomm)
      CB(ctrafo,comm=Ocomm)
      Oend
      Nfence(matrix)
      Nfence(ctrafo)
      end subroutine Mbroadcast4
#   undef CA
#   undef CB

# ifdef old_mpikpt

c
c     Distributes k points onto subcommunicators (and opens unit 6 for output)
c     -> Mcolor  : Color for current node (-> subcommunicators)
c     -> Msub(k) : Index (color) of subcommunicator that will calculate k
c
c     Distribution type (see below)
# define DISTTYPE 2
      subroutine Mdistribute_kpt
      implicit none
      integer              :: kpt1(nkpt),nkpt1,nkpts(nkpt),sym1(nsym),nsym1,pnt(nkpti),nsub,isub,i,n,k
      real_dp              :: time0(nkpti),t,t0
      real_dp, allocatable :: time(:)
# if DISTTYPE <= 1
      real_dp, allocatable :: time1(:)
# endif
      Rbegin
      n = min ( count([(any(eval.and..not.lkpt(i,:)),i=1,nkpti)]) , Osize )
      allocate ( time(n) )
# if DISTTYPE <=1
      allocate ( time1(2:n) )
# endif      
      do ikpt = 1,nkpti
        call getkpt1(kpt1,nkpt1,nkpts,sym1,nsym1,ikpt,0,.false.)
        time0(ikpt) = nkpt1 !; write(*,*) time0(ikpt)
      enddo
      call rorderpf(pnt,time0,nkpti)
      pnt = pnt(nkpti:1:-1)
      do nsub = 2,size(time)
        time(:nsub) = [ (nsub-i,i=1,nsub) ] * 1d-10
        do k = 1,nkpti
          ikpt = pnt(k) ; if(all(lkpt(ikpt,:))) cycle
          t0   = huge(0d0)
          do i = 1,nsub
            n = Osize / nsub ; if(i<=Osize-n*nsub) n = n + 1
            t = time(i) + time0(ikpt) / n
            if(t<t0) then ; t0 = t ; isub = i ; endif
          enddo
          time(isub) = t0
        enddo
# if   DISTTYPE == 0
        time1(nsub) = maxval(time(:nsub)) - minval(time(:nsub))   ! DISTTYPE=0: minimize maximal idle time
# elif DISTTYPE == 1
        time1(nsub) = maxval(time(:nsub)) - sum(time(:nsub))/nsub ! DISTTYPE=1: minimize average idle time
# endif
      enddo
# if DISTTYPE == 2
      nsub        = size(time)                                    ! DISTTYPE=2: always take maximum number of subcommunicators
# else
      nsub        = minloc(time1,1) + 1
# endif
      time(:nsub) = [ (nsub-i,i=1,nsub) ] * 1d-10
      do k = 1,nkpti
        ikpt = pnt(k) ; if(all(lkpt(ikpt,:))) cycle
        t0   = huge(0d0)
        do i = 1,nsub
          n = Osize / nsub ; if(i<=Osize-n*nsub) n = n + 1
          t = time(i) + time0(ikpt) / n
          if(t<t0) then ; t0 = t ; isub = i ; endif
        enddo
        time(isub) = t0 !; write(*,*) ikpt,isub,n,time(isub)
        Msub(ikpt) = isub
      enddo
      deallocate ( time )
# if DISTTYPE <=1
      deallocate ( time1 )
# endif      
      
      Rend
      call Mcast(Msub)
      call Mcast(nsub)
      i = 0
      do isub = 1,nsub
        n = Osize / nsub ; if(isub<=Osize-n*nsub) n = n + 1
        i = i + n        ; if(Orank<i) exit
      enddo
      Mcolor = isub
      end subroutine Mdistribute_kpt

# else

c
c     Simulates MPIKPT run and determines "best" division in subcommunicators
c     Returns  color    : Color for current node (-> subcommunicators)
c              maxcolor : maximal color index (=number of subcommunicators)
c     Updates  kpt_list : list of kpoints to be calculated (will be ordered from high to low computional cost)
c     Parameter       x : 0<x<1, x is the "node-parallelized portion" of the code,
c                         n nodes calculate one kpoint [x/n+1-x]^(-1) times faster than one node
c
      subroutine Mdistribute_kpt(color,maxcolor,kpt_list)
      implicit none
      integer, intent(out)   :: color,maxcolor
      integer, intent(inout) :: kpt_list(:)
      real_dp, parameter     :: x = 0.5d0
      real_dp                :: subtime(min(Osize,size(kpt_list))),tottime(min(Osize,size(kpt_list)))
      real_dp                :: time(size(kpt_list))
      integer                :: pnt(size(kpt_list))
      integer                :: kpt1(nkpt),nkpt1,nkpts(nkpt),sym1(nsym),nsym1
      integer                :: ikpt,nsub,isub,k,n,ix
      if(x<0.or.x>1) Bug('scale_ratio not between 0 and 1.')
      ! order kpt_list according to computational time
      maxcolor = min(Osize,size(kpt_list))
      do k = 1,size(kpt_list)
        ikpt = kpt_list(k)
        call getkpt1(kpt1,nkpt1,nkpts,sym1,nsym1,ikpt,0,.false.)
        time(k) = nkpt1 - k*1d-10
      enddo
      call rorderpf(pnt,time,size(kpt_list))
      kpt_list = kpt_list(pnt(size(kpt_list):1:-1))
      time     = time(pnt(size(kpt_list):1:-1))
      do nsub = 1,maxcolor ! number of subcommunicators
        subtime(:nsub) = [ (k,k=1,nsub) ] * 1d-10
        do k = 1,size(kpt_list)
          isub          = minloc(subtime(:nsub),1)
          n             = Osize / nsub ; if(isub<=Osize-n*nsub) n = n + 1 ! n = number of nodes in subcomm. isub
          subtime(isub) = subtime(isub) + time(k) * ( x/n + (1-x) )
        enddo
        tottime(nsub) = maxval(subtime(:nsub))
      enddo
      n = maxcolor
      do i = n-1,1,-1
        if(tottime(i)<tottime(maxcolor)-1d-10) maxcolor = i ! search backwards to prefer higher maxcolor
      enddo
      if(maxcolor>1) then
        color = mod(Orank,maxcolor) + 1
      else
        color    = 0
        maxcolor = 0
      endif
      end subroutine Mdistribute_kpt

# endif

# endif

c     ------------

c
c     Updates contracted screened interaction for
c     CORES/IBC : add screen matrices for contraction  term (->screenk_mt)
c     COSX      : add screen matrices for coulomb-hole term (->screenk_mt,screenk_pw)
c
c     MPI: Called only by root process(es).
      subroutine screen_contract
      implicit none
      MCOMPLEX_dp, allocatable :: mat(:,:,:)
      complex_dp               :: matc(maxlmindxm,maxlmindxm,ncent) InvC(ctrafo1(nbasp,nbas))
      real                     :: cputime
      integer                  :: g(3)
      integer                  :: ifreq,ic,m,itype,ieq,nlmp
      integer                  :: i,j,k,l
      write(6,'(A'NoA) 'Calculate W contraction... '
      call cpu_time(cputime)
# ifdef INV
      ctrafo1 = ctrafo(:nbasp,:nbas)
      call desymmetrize(ctrafo1,nbasp,nbas,1)
# else
#   define CTRAFO1 ctrafo
# endif
      ! define screenk_mt
      do ifreq = 1,nfreq
        ic = 0
        m  = 0
        do itype = 1,ntype
          nlmp = sum( [ ((2*l+1)*nindxm(l,itype),l=0,lcutm(itype)) ] )
          do ieq = 1,neq(itype)
            ic                   = ic + 1
            matc(:nlmp,:nlmp,ic) = matmac ( matmat ( CTRAFO1(m+1:m+nlmp,:nbas) , screen(:,ifreq) ) -
     &                                      matdiag( CTRAFO1(m+1:m+nlmp,:nbas) , coul            ) ,
     &                                               CTRAFO1(m+1:m+nlmp,:nbas) )
            m                    = m + nlmp
          enddo
        enddo
        if(m/=nbasp) Bug('MT count error.')
        do k = 1,nkpt
          if(kptp(k)==ikpt) call mtrafo_mt(screenk_mt(:,:,:,ifreq),matc,maxlmindxm,maxlmindxm,ikpt,symkpt(k),0,3,.true.)
        enddo
      enddo
      if(job1%type==J_COSX) then
        ! define screenk_pw
        allocate ( mat(ngptm(ikpt),ngptm(ikpt),2) )
        call unitarytrafo( mat(:,:,1) , screen(:,1) - diagonalp(coul) , ctrafo(nbasp+1:,:nbas) , 2 )
        do k = 1,nkpt
          if(kptp(k)==ikpt) then
            call mtrafo_pw Inv(_r) (mat(:,:,2),mat(:,:,1),ngptm(ikpt),ngptm(ikpt),ikpt,symkpt(k),3,.false.)
            do j = 1,ngptm(k)
              do i = 1,ngptm(k)
                g                          = gptm(:,pgptm(i,k)) - gptm(:,pgptm(j,k))
                screenk_pw(g(1),g(2),g(3)) = screenk_pw(g(1),g(2),g(3)) + mat(i,j,2)
              enddo
            enddo
          endif
        enddo
        deallocate ( mat )
      endif
      call cpu_done(cputime)
      end subroutine screen_contract

c     ------------

      end

c     ------------

c     Read Coulomb matrices (barew and screenw) from spex.cou
c     (out) lkpt       - barew/screenw contain contributions from
c     (in)  defb(defs) - read bare (screened) Coulomb matrix
      subroutine readcou(ikpt,barew,screenw,defb,defs,spin,rsite,nsite,freqc,nfreqc)
      use global, only: nwan,img,hartree
      use file
      use, intrinsic :: iso_fortran_env
      implicit none
      integer,    intent(in)  :: nsite,nfreqc,spin
      integer,    intent(out) :: ikpt
      integer,    intent(in)  :: rsite(3,nsite)
      logical,    intent(in)  :: defb,defs
      complex_dp, intent(in)  :: freqc(nfreqc)
      complex_dp, intent(out) :: barew(nwan,nwan,nwan,nwan,nsite)
      complex_dp, intent(out) :: screenw(nwan,nwan,nwan,nwan,nfreqc,nsite)
      logical                 :: foundb(nsite),founds(nfreqc,nsite)
      integer                 :: spn,nsit,nfrq,rsit(3)
      integer                 :: isite,isit,ifrq
      integer                 :: iunit
      integer                 :: i,j,k,l,ios
      logical                 :: ldum
      real_dp                 :: re,im,fr,fi
      character(120)          :: line
      ikpt   = 0
      foundb = .false.
      founds = .false.
      if(defb) barew   = 0
      if(defs) screenw = 0
      inquire(file='spex.cou',exist=ldum)
      if(.not.ldum) return
      iunit = fopen('spex.cou',form='formatted',status='old')
      call nextline ; read(line(12:),*,iostat=ios) spn  ; if(ios/=0)                call fatal('spin')
      call nextline ; read(line(29:),*,iostat=ios) i    ; if(ios/=0.or.i/=nwan)     call fatal('Different Wannier set')
      call nextline ; read(line(17:),*,iostat=ios) nsit ; if(ios/=0.or.nsit<nsite)  call fatal('Two few sites')
      call nextline ; read(line(23:),*,iostat=ios) nfrq ; if(ios/=0.or.nfrq<nfreqc) call fatal('Too few frequencies')
      call nextline ; read(line(19:),*,iostat=ios) ikpt ; if(ios/=0)                call fatal('k point')
      call nextline ; if(line/=' ') call fatal('Expected empty line.')
      call nextline ; if(line==' ') call fatal('Did not expect empty line.')
      do isit = 1,nsit
        call nextline ; i = index(line,')') ; read(line(7:i-1),*,iostat=ios) rsit ; if(ios/=0) call fatal('site (bare)')
        do
          call nextline ; if(line==' ') exit
          read(line,*,iostat=ios) i,j,k,l,re,im ; if(ios/=0) call fatal('bare Coulomb data')
          do isite = 1,nsite
            if(defb.and.all(rsit==rsite(:,isite))) then
              barew(i,j,k,l,isite) = ( re + img*im ) / hartree
              foundb(isite)        = .true.
            endif
          enddo
        enddo
      enddo
      call nextline ; if(line==' ') call fatal('Did not expect empty line.')
      do isit = 1,nsit
        call nextline ; read(line(7:index(line,')')-1),*,iostat=ios) rsit ; if(ios/=0) call fatal('site (screened)')
        do
          call nextline ; if(line==' ') exit
          read(line,*,iostat=ios) fr,fi,i,j,k,l,re,im ; if(ios/=0) call fatal('screened Coulomb data')
          do isite = 1,nsite
            do ifrq = 1,nfreqc
              if(defs.and.all(rsit==rsite(:,isite)).and.abs(freqc(ifrq)*hartree-(fr+img*fi))<1d-10) then
                screenw(i,j,k,l,ifrq,isite) = ( re + img*im ) / hartree
                founds(ifrq,isite)          = .true.
              endif
            enddo
          enddo
        enddo
      enddo
      if(defb.and.any(.not.foundb)) call fatal('Missing sites')
      if(defs.and.any(.not.founds)) call fatal('Missing sites or frequencies')
      call fclose(iunit)
      if(spn/=spin) then
        if(spin==3.and.spn==4.or.spin==4.and.spn==3) then
          write(6,'(/A,I2,A,I2,''.'')') 'Read W with spin index',spn,'. Will be transformed to',spin
          do isite = 1,nsite
            if(defb) call transpose2(barew(:,:,:,:,isite),2)
            do ifrq = 1,nfreqc
              if(defs) call transpose2(screenw(:,:,:,:,ifrq,isite),2)
            enddo
          enddo
        else
          Error('Cannot convert Coulomb matrix read from spex.cou to requested spin index.')
        endif
      endif
      contains
      subroutine nextline
      implicit none
      read(iunit,'(A)',iostat=ios) line ; if(ios/=0) call fatal('next line')
      end subroutine nextline
      subroutine fatal(string)
      character(*), intent(in) :: string
      if(ios/=0) then ; Error('Error while reading '//trim(string)//'. File: spex.cou')
      else            ; Error(trim(string)//'. File: spex.cou')
      endif
      end subroutine fatal
      end

c     ------------

c     Define Frequencies
c     freq (:nfreq)  - (global) Imaginary mesh for self-energy (i*freq)
c     freqr(:nfreqr) - (global) Real      mesh for self-energy
c     freqc(:nfreqc) - (global) Real      mesh for W (or T)
c                               Complex   mesh for spectra
c     contour_par    - Contour parameters
c
c     ltrs = true  : return freqc = 0, ..., max(abs(freqmin),abs(freqmax))
c     ltrs = false : return freqc = freqmin, ..., freqmax
c     freqmin, freqmax obtained from freqc_range
c
c     for testing
c# define freqc_pad
c# define freqc_regular

c begin interface
      subroutine define_freq(job1,contour_par,ospin,ltrs,lwrite)
      use global !inc
      use arrays !inc
      use file
      use key
      use, intrinsic :: iso_fortran_env
      implicit none
      type(jobtype), intent(in)   :: job1
      real_dp,       intent(out)  :: contour_par(4)
      logical,       intent(in)   :: lwrite,ltrs
      integer,       intent(in)   :: ospin
c end interface
      logical                     :: define
      integer                     :: i,j,m
      integer                     :: iunit,istat,freqadd
      real_dp                     :: rdum,freqmax,freqmin
      character(256)              :: line
      character(:),   allocatable :: arg(:)
      character(80)               :: parse_out(3)
      real_dp                     :: ch2r
      integer                     :: ch2i
# ifdef maxfreqc_test
#  warning maxfreqc test included
      real_dp                     :: maxfreqc
# endif

      oselfc = 0
      nfreq  = 0
      nfreqr = 0
      nfreqc = 0
      contour_par = 0

c     (a) Catch trivial cases and return
      ! Nothing to be done
      if(any(job1%type==[J_HF,J_HFE,J_PBE0])) return
      ! Only need w=0
      if(any(job1%type==[J_SX,J_COSX])) then
        allocate ( freq(1),freqc(0) )
        freq  = 0d0
        nfreq = 1
        return
      endif
      ! Get complex mesh for spectra
      if(any(job1%type==[J_SPEC,J_SCRW])) then
        allocate(freqc(0))
        call getfreq(freqc,nfreqc,job1) ; deallocate(freqc) ; allocate(freqc(nfreqc))
        call getfreq(freqc,nfreqc,job1) ; allocate(freq(0))
        return
      endif

c     (b) Read spex.inp (->contour_par,oselfc)
      if(any(job1%type==[J_GW,J_GT])) then
        allocate ( character(80) :: arg(1) )
        call getkey(inp,'CONTINUE', arg, section='SENERGY', status=istat)
        if(istat>=1) then ! CONTINUE
          if(size(arg)>=2) Error('More than one argument after CONTINUE.')
          deallocate ( arg )
          call getkey(inp,'CONTOUR', arg, section='SENERGY', status=istat)
          if(istat/=0) Error('Keywords CONTINUE and CONTOUR given at the same time.')
          oselfc = 1
        else              ! CONTOUR
          deallocate ( arg )
          call getkey(inp,'CONTOUR', arg, section='SENERGY', status=istat)
          if(istat==0) then
            oselfc = 1    ! CONTINUE is default
          else
            if(istat==1)    Error('Missing arguments after CONTOUR.')
            if(size(arg)>2) Error('Expected one or two arguments after CONTOUR.')
            if(size(arg)==2) contour_par(4) = ch2r(arg(2),'*CONTOUR')
            if(contour_par(4)<0) Error('Second argument to CONTOUR must be positive.')
            if(index(arg(1),'{')/=0) then
              if(arg(1)(:2)=='+-') then
                oselfc = 3 ; if(lwrite) write(6,'(/A)') 'Parameters for contour integration:'
                if(index(arg(1),'..')/=0) then
                  call parser(parse_out,arg(1), '+-{1..2,3}' ,'CONTOUR')
                  Info('".." obsolete; use ":" instead.')
                else
                  call parser(parse_out,arg(1), '+-{1:2,3}' ,'CONTOUR')
                endif
                if(any(parse_out(:2)=='*')) Error('Wildcards can only be used for CONTOUR {...}.')
              else
                oselfc = 4 ; if(lwrite) write(6,'(/A)') 'Parameters for contour integration:'
                if(index(arg(1),'..')/=0) then
                  call parser(parse_out,arg(1),'{1..2,3}','CONTOUR')
                  Info('".." obsolete; use ":" instead.')
                else
                  call parser(parse_out,arg(1),'{1:2,3}','CONTOUR')
                endif
              endif
              if(parse_out(1)/='*') then
                contour_par(1) = ch2r(parse_out(1),'*CONTOUR')
                if(lwrite) write(6,'(A,F10.5)') 'User-defined lower bound of contour:',contour_par(1)
              else
                contour_par(1) = minval( [ ( ene(job1%band(i),job1%kpt(i),job1%spin(i)), i=1,size(job1%band) ) ] ) - efermi - 0.07
                if(lwrite) write(6,'(A,F10.5)') 'Automatic lower bound of contour:',contour_par(1)
              endif
              if(parse_out(2)/='*') then
                contour_par(2) = ch2r(parse_out(2),'*CONTOUR')
                if(lwrite) write(6,'(A,F10.5)') 'User-defined upper bound of contour:',contour_par(2)
              else
                contour_par(2) = maxval( [ ( ene(job1%band(i),job1%kpt(i),job1%spin(i)), i=1,size(job1%band) ) ] ) - efermi + 0.07
                if(lwrite) write(6,'(A,F10.5)') 'Automatic upper bound of contour:',contour_par(2)
              endif
              contour_par(3) = ch2r(parse_out(3),'*CONTOUR')
              if(contour_par(3)<=0)             Error('increment must be positive.')
              if(contour_par(1)>contour_par(2)) Error('lower > upper bound.')
            else
              oselfc           = 2
              contour_par(1)   = ch2r(arg(1),'*CONTOUR')
              contour_par(2:3) = 0
              if(contour_par(1)<=0) Error('first argument to CONTOUR must be positive.')
            endif
          endif
        endif
        if(allocated(arg)) deallocate ( arg )

c       (c) Real frequencies (contour integration)
        if(oselfc==1) then
          nfreqr   = 1 ; allocate ( freqr(nfreqr) )
          freqr(1) = 0d0
        else if(oselfc==2) then
          nfreqr   = 2 ; allocate ( freqr(nfreqr) )
          freqr(1) = - contour_par(1)
          freqr(2) =   contour_par(1)
        else
          rdum = ( contour_par(2) - contour_par(1) ) / contour_par(3) ; nfreqr = nint(rdum) + 1 ; allocate ( freqr(nfreqr) )
          rdum = ( contour_par(2) - contour_par(1) ) / (nfreqr-1)
          if(nfreqr==1) Error('Frequency range for self-energy only contains one point.')
          do i = 1,nfreqr
            freqr(i) = contour_par(1) + rdum * (i-1)
          enddo
          if(lwrite) write(6,'(A)') 'Number of frequencies for interpolation:'//trim(chr(nfreqr,'I5'))//
     &                              '  ('//trim(chr(rdum,'F8.5'))//' )'
        endif
      endif

c     (c) Imaginary frequencies
      call getkey(inp,'MESH',arg,section='SENERGY',default=['10','10'])
      if(size(arg)==0) then
        Error('Missing arguments after MESH.')
      else if(size(arg)==1) then
        if(arg(1)(:1)//arg(1)(len_trim(arg(1)):)/='""') Error('Filenames have to be given in quotes "..." (MESH)')
        arg(1) = arg(1)(2:len_trim(arg(1))-1)
        if(lwrite) write(6,'(/A)') 'Frequency mesh on imaginary axis is taken from file '//trim(arg(1))//'.'
        iunit  = fopen(arg(1),status='old')
        define = .false.
 1      i      = 0
        do
          read(iunit,'(A)',iostat=j) line
          if(j/=0) exit
          if(index(line,'#')/=0) line(index(line,'#'):) = ' '
          if(line/=' ') then
            i = i + 1
            if(define) read(line,*) freq(i)
          endif
        enddo
        nfreq = i
        if(.not.define) then
          define = .true.
          allocate ( freq(nfreq) )
          rewind(iunit)
          goto 1
        endif
        call rorder(freq,nfreq)
        if(freq(1)/=0) Error('imaginary frequency mesh must contain zero frequency.')
        call fclose(iunit)
      else if(size(arg)==2) then
        freqadd = 0
        i       = index(arg(1),'+')
        if(i/=0) then
          freqadd    = ch2i(arg(1)(i+1:),'MESH')
          arg(1)(i:) = ' '
        endif
        nfreq   = ch2i(arg(1), 'MESH') ; if(nfreq<=0) Error('First argument to MESH must be positive.')
        freqmax = ch2r(arg(2),'*MESH')
        if(define) then
        endif
        allocate ( freq(nfreq+freqadd*(freqadd+1)/2) )
        do i = 1,nfreq
          if(freqmax>0) then ; rdum = 0.9d0 * (i-1)/(nfreq-1)
          else               ; rdum = 1.0d0 * (i-1)/ nfreq
          endif
          freq(i) = rdum / (1-rdum)
        enddo
        if(nfreq>1) then
          rdum         = abs(freqmax) / freq(nfreq)
          freq(:nfreq) = rdum * freq(:nfreq)
        endif
        j = nfreq
        do m = 1,freqadd
          do i = 1,m
            j       = j + 1
            freq(j) = freq(freqadd-m+1) + (freq(freqadd-m+2)-freq(freqadd-m+1)) / (m+1) * i
          enddo
        enddo
        if(j/=size(freq)) Bug('Count error in freq.')
        nfreq = j
        call rorder(freq,nfreq)
      else
        Error('Expected maximally two arguments after MESH.')
      endif
      deallocate ( arg )

c     (d) Frequencies for screened interaction (->freqc)
      if(any(job1%type==[J_GW,J_GT]).and.oselfc>1.and.contour_par(4)/=0) then
        call freqc_range(freqmin,freqmax,job1,ospin) 
# ifdef freqc_pad
#   warning freqmin freqmax padded
        Warn('freqmin/freqmax range padded (0.3)')
        freqmin = freqmin - .3
        freqmax = freqmax + .3
# endif
# ifdef freqc_regular
#   warning regularized freqmin/freqmax
        freqmin = floor(freqmin/contour_par(4))*contour_par(4)
        freqmax = ceiling(freqmax/contour_par(4))*contour_par(4)
# endif
        if(ltrs) then
          freqmax = max(abs(freqmin),abs(freqmax))
          freqmin = 0
        endif
# ifdef maxfreqc_test
        write(6,'(A'NoA) 'maxfreqc test... '
        rdum   = maxfreqc(job1)
        if(abs(freqmax-rdum)>1d-10) Bug('maxfreqc test failed: '//Chr(rdum)//' vs. '//Chr(freqmax))
        write(6,'(A)') 'passed'
# endif
        nfreqc = 1 + max ( 1 , nint ( (freqmax-freqmin)/contour_par(4) ) ) ; allocate ( freqc(nfreqc) )
        freqc  = [ ( ( (nfreqc-1-i)*freqmin + i*freqmax ) / (nfreqc-1) , i=0,nfreqc-1 ) ]
        if(lwrite) write(6,'(A)') 'Number of frequencies for interpolation:'//trim(chr(nfreqc,'I5'))//
     &                            '  ('//trim(chr((freqmax-freqmin)/(nfreqc-1),'F8.5'))//' )'
        if(lwrite) write(6,'(A)') 'Frequency range for interpolation from '//Chf(freqmin,'F9.5')//' to '//Chf(freqmax,'F9.5')
      else
        allocate ( freqc(0) )
      endif

      end

c     ------------

      subroutine getfreq(frq,nfrq,job1)
      use global
      use file
      use, intrinsic :: iso_fortran_env
      implicit none
      integer,       intent(inout) :: nfrq
      complex_dp,    intent(out)   :: frq(nfrq)
      type(jobtype), intent(in)    :: job1
      character(80)                :: line
      integer                      :: iunit,i,ios
      real_dp                      :: rdum
      if(job1%freqfile/=' ') then
        iunit = fopen(job1%freqfile)
        i = 0
        do
          read(iunit,'(A)',iostat=ios) line
          if(ios/=0) exit
          if(index(line,'#')/=0) line = line(:index(line,'#')-1)
          if(line/=' ') then
            i = i + 1
            if(nfrq==0) then
              read(line,*,iostat=ios)
            else
              if(index(line,'(')/=0) then ; read(line,*,iostat=ios) frq(i)
              else                        ; read(line,*,iostat=ios) rdum ; frq(i) = rdum
              endif
            endif
            if(ios/=0) Error('Error in frequency file.')
          endif
        enddo
        nfrq = i
        call fclose(iunit)
      else
        if(nfrq==0) then
          nfrq   = nint((job1%freq(2)-job1%freq(1))/job1%freq(3)) + 1
        else
          frq(1) = job1%freq(1)
          do i = 2,nfrq
            frq(i) = frq(i-1) + job1%freq(3)
            if(abs(frq(i))<1d-12) frq(i) = 0
          enddo
        endif
      endif
      end

c     ------------

      subroutine getctrafo(Win(ctrafo),nn,gpw,ikpt)
      use global
      use, intrinsic :: iso_fortran_env
      implicit none
      integer,    intent(in)  :: nn,gpw(nn),ikpt MpiC(win_ctrafo)
      MCOMPLEX_dp,intent(out) :: ctrafo(nbasm(ikpt),nn)
      complex_dp              :: ctrafo1(nbasm(ikpt))
      integer                 :: igpt,igptp,i,j,itype,ieq,ic,ng,l,m,n,lm MpiC(Merr)
      real_dp                 :: kp(3),q(3),qnorm,sphbes(maxgrid,0:maxlcutm),csphbes(maxindxm,0:maxlcutm)
      complex_dp              :: y((maxlcutm+1)**2),cexp
      real_dp                 :: intgrf
      Mpi( if(win_ctrafo/=MPI_UNDEFINED) then ; Nfence(ctrafo) ; endif )
      kp = kpt(:,ikpt)
      do igpt = 1,nn ; Mpi( if(win_ctrafo/=MPI_UNDEFINED) then ; Ncycle(igpt-1) ; endif )
        igptp = pgptm(gpw(igpt),ikpt)
        q     = matmul(rlat,kp+gptm(:,igptp))
        qnorm = sqrt(sum(q**2))
        j     = 0
        ic    = 0
        call harmonicsr(y,q,maxlcutm)
        y = conjg(y)
        do itype = 1,ntype
          ng = grid(itype)%number
          do i = 1,ng
            call sphbessel(sphbes(i,:lcutm(itype)),rgrid(i,itype)*qnorm,lcutm(itype))
          enddo
          do l = 0,lcutm(itype)
            do n = 1,nindxm(l,itype)
              csphbes(n,l) = intgrf( rgrid(:ng,itype)*basm(:ng,n,l,itype)*sphbes(:ng,l), itype )
            enddo
          enddo
          do ieq = 1,neq(itype)
            ic   = ic + 1
            cexp = 4*pi * exp(img* 2*pi *dot_product(gptm(:,igptp),cent(:,ic)))
            lm   = 0
            do l = 0,lcutm(itype)
              do m = -l,l
                lm = lm + 1
                do n = 1,nindxm(l,itype)
                  j          = j + 1
                  ctrafo1(j) = cexp * img**l * y(lm) * csphbes(n,l) / sqrt(vol)
                enddo
              enddo
            enddo
          enddo
        enddo
        do i = 1,ngptm(ikpt)
          if(i==gpw(igpt)) then ; ctrafo1(nbasp+i) = 1
          else                  ; ctrafo1(nbasp+i) = 0
          endif
        enddo
        Inv( call symmetrize(ctrafo1,nbasm(ikpt),1,1) )
        Inv( if(any(abs(imag(ctrafo1))>1d-12)) Error('Nonzero imaginary part.') )
        ctrafo(:,igpt) = ctrafo1
      enddo
      Mpi( if(win_ctrafo/=MPI_UNDEFINED) then ; Nfence(ctrafo) ; endif )

      end
