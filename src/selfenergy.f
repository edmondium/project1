c--------------------------------------------------------------------------------
c Copyright (c) 2021 Peter Grünberg Institut, Forschungszentrum Jülich, Germany
c This file is part of SPEX and available as free software under the conditions
c of the MIT license as expressed in the LICENSE file in more detail.
c--------------------------------------------------------------------------------

c Adds the contribution of the current k-point (ikpt) (and its symmetry-equivalent points) to the matrix elements
c of the correlation self-energy (-> selfc).
c
c job1%full = .true.  : full self-energy matrix including off-diagonal elements
c                       in block-matrix form according to block(:,:)
c job1%full = .false. : only diagonal elements
c
c (1) If the self-energy is later analytically continued from the imaginary to the real frequency axis, we evaluate
c     the matrix elements <n|SIGMA(iw)|m> by frequency convolution (along the imaginary frequency axis)
c
c          s         s          s        -OMEGA        3   u/o  +inf             s                  c               s          s  ~         ~       s       s
c     < phi   | SIGMA (iw) | phi   >  =  ------ INT   d k  SUM  INT dw' (iw+iw'-E      )^(-1) SUM  W  (k,iw')  < phi      | phi   M    >  < M    phi   | phi      >
c          nq        c          mq       16pi^4    BZ       n'  -inf             n'q+k         IJ   IJ              n'q+k      mq  k,I       k,J    nq      n'q+k
c
c (2) If the frequency integration is performed by contour integration, we calculate
c
c          s         s          s
c     < phi   | SIGMA (w0) | phi   >  by summing
c          nq        c          mq
c
c     (a) the integral along the imaginary frequency axis using the formula (1) above with the substitution
c
c         iw -> w0     (replacing a purely imaginary by a purely real argument)
c
c     (b) the sum over the residues, i.e., the integrated values of W at the poles of the Green function
c
c             OMEGA        3  u/o             s                 c     s               s         s  ~        ~       s       s
c         +/- ----- INT   d k SUM  theta[+/-(E     +w0)]  SUM  W  (k,E     -w0)  < phi     | phi   M    > < M    phi   | phi      >
c             8pi^3    BZ      n'             n'q+k        IJ   IJ    n'q+k           n'q+k     nq  k,I      k,J    nq      n'q+k
c
c         for w0 < 0 / w0 > 0. Both contributions are calculated at the frequencies freqr1(:nfreqr) which are calculated
c         from freqr(:nfreqr) and the band energies (see below).
c         theta is the Heaviside function.
c         Contributions (a) and (b) are summed.
c
c Input:
c   screen = screened interaction W
c   head   = head of W
c   ctrafo = transformation coefficients of eigenvectors of v wrt the mixed basis (screen is defined in the basis of eigenvectors)
c   n      = rank of W
c   ikpt   = current summation k-point (from the irreducible zone)
c
c Imported from global:
c   block  = definition of the block matrices (taking into account irreps, see irrep.f)
c   nfreq  = number of mesh points along imaginary frequency axis
c   freq   = mesh points along imaginary frequency axis
c   nfreqr = number of (real) frequencies w0 at which the self-energy is evaluated
c   freqr  = list of those frequencies (relative to efermi or the KS energy, see correlation.f)
c   nfreqc = number of (real) frequencies which are used to interpolate W (i.e., screenc and headc)
c            use Pade extrapolation for nfreqc = 0
c   freqc  = list of those frequencies
c
c The divergence at the Gamma point is integrated as in exchange.f (but anisotropy is taken into account in divergence_h/c).
c
c
c
c USE OF SYMMETRY:
c
c
c We explain the calculation of <i|GW|j>(ikpt), which is the contribution of the current k point ikpt to <i|GW|j>.
c i(r) and j(r) are states at q. (Prefactors such as "img" are omitted for simplicity.)
c
c <i|GW|j>(ikpt) = SUM(n) SUM(k0) INT dw <i (nk) | W | (nk) j >
c                = SUM(n) SUM(k0) INT dw INT j(r) nk(r)* W(r,r') nk(r') i(r')* dr dr',
c
c where we have used the symmetry W(r,r')=W(r',r). INT dw is the frequency convolution. The k0 summation runs over all k0 that
c are symmetry equivalent to ikpt. The definition is such that k=q+k0 rather than k=q-k0.
c
c Using symmetry, we can restrict the summation to the subset {k1} of {k0}. The set {k1} [kpt1(:nkpt1)] includes only those
c points k1 that are NOT pairwise related by a symop of the little group, which is defined as the set of symops P that leave
c q invariant [sym1(:nsym1)]. In other words, we have removed all symmetry equivalents (wrt sym1) from the set {k0} and have
c retained only one representive of each group [size nkpts(k1)] of symmetry-equivalent k0 points. The contribution of those
c that have been removed is taken into account by symmetrization. (The k1 are elements of EIBZ and equivalent to ikpt.)
c The little group also contains the inverses P^(-1). We use that W is invariant under the spatial part R of P, i.e.,
c W=RWR. (If P does not include time reversal (tr), R=P.)
c
c <i|GW|j>(ikpt) = SUM(n) SUM(k0) INT dw                                < i (nk) | W | (nk) j >
c                = SUM(n) SUM(k1) INT dw nkpts(k1)/nsym1   SUM(P)  <i P^(-1)(nk) | W | P^(-1)(nk) j>
c                = SUM(n) SUM(k1) INT dw nkpts(k1)/nsym1 { SUM(P,no tr) <Pi (nk) | W | (nk) Pj>
c                                                        + SUM(P,tr)    <Pj (nk) | W | (nk) Pi> } ,
c
c where k stands for q+k0 and q+k1, respectively.
c
c The tr case follows from
c <i P^(-1)(nk) | W | P^(-1)(nk) j> = <Ri (nk)* | W | (nk)* Rj> = <Pi* (nk)* | W | (nk)* Pj* > = <Pj (nk) | W | (nk) Pi> .
c Note that, if W is Hermitian (e.g., W=W(iw) or W=v), one can also write
c <i P^(-1)(nk) | W | P^(-1)(nk) j> = <Pi (nk) | W | (nk) Pj>* .
c
c Introducing the irreps P(i'i) = <i'|Pi> (defined in the subspace including state i), we have
c <i|GW|j>(ikpt) = SUM(n) SUM(k1) INT dw nkpts(k)/nsym1 SUM(i'j') { [ SUM(P,no tr) P(i'i)* P(j'j) ]  <i' (nk) | W | (nk) j'>
c                                                                 + [ SUM(P,tr)    P(i'i)* P(j'j) ]* <j' (nk) | W | (nk) i'> }
c                = SUM(n) SUM(k1) INT dw nkpts(k)/nsym1 SUM(i'j') {       irrep_contr(ij,i'j')       <i' (nk) | W | (nk) j'>
c                                                                 +       irrep_contt(ij,j'i')       <j' (nk) | W | (nk) i'> } ,
c which defines the contracted irreps irrep_contr and irrep_contt.
c
c In the following, we discuss the evaluation of <i(nk)|W|(nk)j>. The operator W is given in an auxiliary Bloch basis,
c <a|W|b>. The auxiliary basis can be the mixed basis (W=v) or the basis of Coulomb eigenvectors (W=W(w)).
c The matrices of W are only known in the irreducible wedge. Therefore, we use a suitable transformation [symkpt(k1)]
c P: IBZ->EIBZ, ikpt->k1, PI=I', Pa=a', etc.. In the following, summing over I, J, a, b is implicit.
c
c First case: P does not include tr:
c   <i|(nk)W(nk)|j> = <nk|jI'> <I'|a'> <a'|W|b'> <b'|J'> <J'i|nk>
c                   = <nk|jPI> <PI|Pa> <Pa|W|Pb> <Pb|PJ> <PJi|nk>
c                   = <nk|jPI>  <I|a>   <a|W|b>   <b|J>  <PJi|nk> .
c
c Second case: P includes tr:
c   <i|(nk)W(nk)|j> = <nk|jI'> <I'|a'> <a'|W|b'> <b'|J'> <J'i|nk>
c                   = <nk|jPI> <PI|Pa> <Pa|W|Pb> <Pb|PJ> <PJi|nk>
c                   = <nk|jPI>  <a|I>   <b|W|a>   <J|b>  <PJi|nk>
c                   = <nk|iPJ>* <J|b>   <b|W|a>   <a|I>  <PIj|nk>* ,
c
c which differs from the first case only by a complex conjugation of the vectors and a transpose operator. If W is Hermitian
c (i.e, for W(iw) and v), we may also use
c
c   <i|(nk)W(nk)|j> =   <nk|jPI> <PI|Pa> <Pa|W|Pb> <Pb|PJ> <PJi|nk>
c                   =   <nk|jPI>  <I|a>*  <a|W|b>*  <b|J>* <PJi|nk>
c                   = [ <nk|jPI>* <I|a>   <a|W|b>   <b|J>  <PJi|nk>* ]* ,
c
c i.e., we may replace the transpose operator by a complex conjugate.
c
c When testing symmetry, note that Pade interpolation is not a linear operation. This means that CONTOUR with Pade and FREQINT PADE
c will produce slightly different results if symmetry is switched off or not (switch_off_symmetry macro in susceptibility.f).
c
c For the same reason, one can get different results in the serial and parallelized version: In the serial version, the Green-function
c band summation is (degenerate-)subspace-wise and equivalent-k-point-wise. In the parallelized version, the summation is over
c degenerate bands and equivalent kpoints only if MPISYM is set. (By default, this summation is not made because of additional MPI
c communication). The sum is identical, but the non-linear Pade operation on summands can lead to deviations.
c
c
c MPI flowchart of loop over blocks:
c      
c do ! block loop
c   MPIBLK:  cycle if another subcommunicator calculates this block
c   COPYDEG: cycle if another degenerate state has been calculated already (disabled for FULL calculations)
c   Mblocks determines subdivision of (1:size(band),1:nkpt1,1:nbnd) onto processes of current subcommunicator:
c     (b1:b2,k1:k2,u1:nbnd:us), so third dimension is not a block but "u-strides" over full range
c     Mblocks may decide to take few processes idle, which are then split from subcommunicator
c   MEM: if memory not sufficient, divide work into packets that are carried out consecutively
c     division done collectively to make sure that all processes run the same number of packets      
c
c   do ! MEM packet loop
c     Exchange self-energy: If (band) is divided among processes, array cprod needs to be bcasted so that
c       full self-energy matrix can be calculated -> matx(1:size(band),1:size(band))
c
c     GW SELF-ENERGY
c       Determine sets of {ikpt2}. Each set contains only equivalent ikpt2 (so matrices can be summed before w integration)
c           
c     do ! ikpt2 set loop
c       GW self-energy: If (band) is divided among processes, array cprod needs to be bcasted so that
c         each process can calculate full matrix -> mat(1:size(band),1:size(band),:)
c         From now on, there is no (band) division anymore and the u-stride gets increased accordingly -> us1
c
c       do ! loop over u-bands
c         MPISYM: Average over degenerate states in same process and other processes
c         iw integral (along imaginary axis) -> selfc_hlp
c       enddo
c
c       --- same for residues (-> selfc_hlp) ---
c
c     enddo ! ikpt2 set loop
c   enddo ! MEM packet loop
c
c   MPI reduce selfc_hlp and update selfc(:,:)
c
c enddo ! block loop
c
c Note:
c - Several dummy beginSplit/endSplit pairs are necessary to avoid deadlock for ikpt2 sets (nk_dead) and MPISPLIT (nu_dead).
c
c
c Note:
c - wings of W are used only to give the Gamma point contribution to the body of W
c - on exit, the arrays screen and screenc are destroyed if cblock is allocated (see macro WBLOCK in susceptibility).
c - on exit, the array coulomb0 is destroyed if enable_multipole is true (see keyword MULTIPOLE).
c

c Only calculate one of the identical diagonal self-energy matrix elements of
c degenerate states and copy the others (only serial version and .not.full).
c This can be disabled here by commenting out the following.
# define COPYDEG

# ifndef COPYDEG
#   warning COPYDEG unset for testing
# endif

# ifdef MPI
#   define NFENCE Nfence(screen) ; if(nfreqc>0) then ; Nfence(screenc) ; endif
# else
#   define NFENCE continue
# endif

# ifndef M_PART
#   define M_PART 1
#   include __FILE__
#   undef M_PART
#   define M_PART 2
#   include __FILE__
#   undef M_PART
#   define M_PART 3
# endif

c -----------------------------
# if M_PART == 1
c -----------------------------

# include "cppmacro.h"
# include "jobtype.h"

c begin interface
      subroutine selfenergy(job1,ikpt,eval,nbas,coul,Win(ctrafo),Win(screen), head, wing,
     &                                                           Win(screenc),headc,wingc,plasma)

      use global !inc
      use arrays !inc
      use util
      use timer_util
      use wrapper
      use freq_integral
      use file
      use, intrinsic :: iso_fortran_env
      Mpi ( use Mwrapper )
      Mpi ( use key )
      Mpi2( use, intrinsic :: iso_c_binding )
      Load( use readwrite ) LoadC( only: read_wavef0 ) LoadC( read_wavef2 )

      implicit none
      type(jobtype), intent(in)           :: job1
      integer,       intent(in)           :: nbas,ikpt
      logical,       intent(in)           :: eval(2)
      real_dp,       intent(in), optional :: coul(nbas)
      complex_dp,    intent(in), optional :: headc(3,3,nfreqc),wingc(3,nbas,2,nfreqc)
      MCOMPLEX_dp,   intent(in), optional :: head(6,nfreq),wing(3,nbas,nfreq)
      MCOMPLEX_dp,               optional :: ctrafo(nbasm(ikpt),nbas)      ! can be modified but leaves subroutine unchanged
      MCOMPLEX_dp,               optional :: screen(nbas*(nbas+1)/2,nfreq) ! modified by cblock
      complex_dp,                optional :: screenc(nbas,nbas,nfreqc)     ! modified by cblock
      real_dp,                   optional :: plasma                        ! changed if plasma=-1d0 (PLASMA METAL)
# ifdef MPI
      integer,       intent(in), optional :: win_screen,win_screenc,win_ctrafo
# endif
c end interface
      complex_dp,    allocatable         :: matc(:,:,:),matc_(:,:,:),matz(:,:),cvec(:)
      complex_dp,    allocatable         :: cprod_ibc(:,:,:,:)
      complex_dp,    allocatable, target :: cprod2c(:,:)
      MCOMPLEX_dp,   allocatable         :: mat(:,:,:),matx(:,:),mat_(:,:,:,:),vec(:)
      MCOMPLEX_dp,   allocatable, target :: cprod(:,:,:),cprod2(:,:)
      MCOMPLEX_dp,   pointer_cnt         :: cprodp(:,:,:)
      MCOMPLEX_dp,   allocatable         :: irrep_contr(:,:,:,:),irrep_contt(:,:,:,:)
      MCOMPLEX_dp,   allocatable         :: moment(:,:,:),moment_diag(:,:)
      MCOMPLEX_dp,   pointer_cnt         :: coulmat(:,:),coulloc(:),screen_(:,:),screen2_(:,:,:)
      MCOMPLEX_dp,   pointer_cnt         :: eigv(:,:),olap(:,:)
      real_dp,       allocatable         :: eig(:)
      real_dp,       allocatable         :: wintgrc(:,:,:,:)
      integer                            :: nwblock
      integer,       allocatable         :: band(:),pnt(:),wblock(:),bfreqc(:,:,:),bu(:,:,:)
      complex_dp,    allocatable         :: wfreqintegral(:,:,:),wfreqintegral1(:,:,:),wfreqintegral2(:,:,:)
      complex_dp,    allocatable         :: pm(:,:,:,:),pv(:,:,:),zm(:,:,:,:),cpade_head(:,:,:)
      complex_dp,    allocatable         :: selfc_ibc(:,:,:),selfc_hlp(:,:,:)
      complex_dp                         :: h0(3,3),h1(3,3),h2(3,3)
      complex_dp                         :: cpade(nfreq+1),pole(nfreq*smooth(2)),resid(nfreq*smooth(2)),wfreq(0:3,nfreq)
      complex_dp                         :: cdum,cdum1,cdum2
      complex_dp                         :: cmat(3,3),hw(3,3,max(nfreq,nfreqr),0:2),divergence_c(nfreqc),iheadc(3,3,nfreqc)
      MCOMPLEX_dp                        :: aspline(0:3,nfreq)
      MCOMPLEX_dp                        :: ihead(6,nfreq),mom(3)
      real_dp                            :: dene(3),dene2(3,3),ddene(3,3),divfac_x,divergence_x,divergence_h(nfreq)
      real_dp                            :: eneq,enek,enediff,freqr1(nfreqr),weight,wght,rdum,rdum1,coeff(nfreqc)
      real_dp                            :: frq(27),pfrq(0:3,26),emin,emax,f1,f2,logdiv1,logdiv2
      real_dp                            :: memo,memd
      integer                            :: maxpack
      integer                            :: minb(nfreqr,nspin1),maxb(nfreqr,nspin1),npole1
      integer                            :: nkpt1,kpt1(nkpt),sym1(nsym),nsym1,nkpts(nkpt)
      integer                            :: nfrq,nfreq1,ifreq,ifreqr,ifreqc,bfreqc1,bfreqc2
      integer                            :: iself,iselfc,iselfx,iblock,ib,ibandq,ikptq,ispin,s,iband,sizeband,id
      integer                            :: ndeg,deg0,deg1,band1,band2,nbnd,nbnd1,nbnd0,ibnd
      integer                            :: isub,jsub,isym
      integer                            :: ikpt1,ikpt2,ikpt2_
      integer                            :: i,j,k,l,m,p,ka,kb,iuf,sgn,itype
      logical                            :: drude_separate,newkptq,evalx,evalc,ldum,enable_multipole
      real                               :: time_mt,time_pw,time_trafo,time_exch,time_gamma,time_mat,time_equiv,time_freq
      real                               :: time_ibc,time_tot,time_tetra MpiC(time_idle) MpiC(time_maxidle)
      integer                            :: b1,b2,nb,k1,k2,nk,u1,us,nu,nu0,nu1,nu2,u2,u,iu,iu2,bu1,bu2
      integer                            :: kptsum
      complex_dp                         :: pade_func
      real_dp                            :: logdiv,error_coulomb_multipole
      real_dp, external                  :: scaling_selfenergy_cprod,scaling_linear
# ifdef LOAD
      integer                            :: kpt2(nkpt)
# endif
# ifdef MPI      
      integer                            :: win_eigv,win_olap,win_coulmat,win_coulloc,win_screen_,win_screen2_
      logical, save                      :: first = .true.,mpisym = .false.
      real_dp, save                      :: Mover
      integer                            :: Mcolor,Merr,Msub(nblock),Mstat(mpi_status_size),rank,Mfile
      integer                            :: Mblock(nblock)
      complex_dp,  allocatable, target   :: cprod3c(:,:)
      complex_dp,  pointer_cnt           :: cprod_c(:,:)
      MCOMPLEX_dp, allocatable, target   :: cprod3(:,:)
      MCOMPLEX_dp, pointer_cnt           :: cprod_(:,:),Ninit_coulmat(:)
      type(c_ptr)                        :: ptr
      integer                            :: bounds(4,3)
      integer                            :: b3,b4,nb1,irank,urank,Mrank1
      integer                            :: u1_,us_,nu_,nu_dead,nk_dead,nsplit
# endif

      interface equivalent_kpoints
      procedure equivalent_kpoints_r,equivalent_kpoints_c
      end interface

      call timer_start('Routine selfenergy SLF')

      evalx = eval(1)
      evalc = eval(2)

      nullify(coulmat)
      nullify(coulloc)
      nullify(screen_)
      nullify(screen2_)

      enable_multipole = l_multipole .and. evalx .and. all(job1%type/=[J_SX,J_COSX])

# ifdef MPI
      if(first) then
        Rbegin
        Mover = 0
        if(job1%type==J_GW.and.job1%full) then ; Mover = 0  ! only parallelize over n'     (default for GW FULL)
        else                                   ; Mover = 10 ! also parallelize over blocks (default otherwise)
        endif
        call getkey(inp,'MPIBLK',Mover,section='SENERGY',status=i,mini=0d0)
        if(i==1) Mover = 10 ! default
        if(nfreqc==0) call getkey(inp,'MPISYM',mpisym,section='SENERGY',default=.false.)
        if(mpisym.and.nfreqc/=0.and.freqint==2) then
          Info('MPISYM disabled for present case (nfreqc/=0 and FREQINT PADE)')
          mpisym = .false.
        endif
        Rend
        first = .false.
        call Mcast(Mover)
        call Mcast(mpisym)
      endif
# endif

# ifdef CHECK_SENERGY
#   warning CHECK_SENERGY defined
                    Rwrite(6,'(''screen: '',F30.14)') sum(abs(screen)**2)
      Rif(ikpt==1)   write(6,'(''head:   '',F30.14)') sum(abs(head)**2)
      Rif(nfreqc>0) then
                     write(6,'(''screen:'',F30.14)') sum(abs(screenc)**2)
        Rif(ikpt==1) write(6,'(''head:  '',F30.14)') sum(abs(headc)**2)
      endif
      Rwrite(6,*)
# endif

      if(all(job1%type/=[J_HF,J_GW,J_RPA,J_HFE,J_SX,J_COSX,J_PBE0]))
     &                                  Bug('Wrong type of calculation (job1%type).')
      if(.not.(evalx.or.evalc))         Bug('evalx and evalc both false.')
      if(evalc.and.nfreqr==0)           Bug('nfreqr is zero.')
      if(nfreqc>0) then
        if(any(abs(imag(freqc))>1d-10)) Bug('Nonzero imaginary part in freqc.')
      endif

      Rif(ikpt==1.and.evalc.and.freqint>=2) then
        if(ozero)     Info('FREQINT SPLINE assumed for zero-order corrections.')
        if(job1%full) Warn('FREQINT NONLIN never tested for FULL calculations. Use with care!')
      endif

      if(oibc/=0.and.allocated(cblock)) then
        Info('cblock deallocated for keyword IBC')
        deallocate ( cblock )
      endif

      maxpack    = 1
      time_mt    = 0
      time_pw    = 0
      time_trafo = 0
      time_exch  = 0
      time_gamma = 0
      time_mat   = 0
      time_equiv = 0
      time_freq  = 0
      time_ibc   = 0
      time_tetra = 0

c
c     Define inverse head (-> ihead/c)
c
      if(ikpt==1) then
        if(present(screen)) then
          Mpi( ihead = 0 )
          do i = Mrange(1,nfreq)
            h0         = unpackmat(head(:,i)) ; call invert_angular(h0)
            ihead(:,i) = 4*pi * packmat(h0)
          enddo
          Mpi( call Msum(ihead) )
          if(nfreqc>0) then
            Mpi( iheadc = 0 )
            do i = Mrange(1,nfreqc)
              h0            = headc(:,:,i) ; call invert_angular(h0)
              iheadc(:,:,i) = 4*pi * h0
            enddo
            Mpi( call Msum(iheadc) )
          endif
        endif
      endif

c
c     Define Gamma divergence (->divergence/_h/_c)
      if(ikpt==1) then
        call timer_start('SLF divergence')
        if(divergence==0) then
          RWarn('Contribution of Coulomb divergence has not been calculated yet. This should not happen!')
          call gamma_divergence(.false.)
        endif
        if(present(screen)) then
          do i = 1,nfreq  ; rdum = trace(ihead(:,i))/3    ; call gamma_divergence_h(divergence_h(i),head(:,i),   rdum) ; enddo
          do i = 1,nfreqc ; cdum = trace(iheadc(:,:,i))/3 ; call gamma_divergence_c(divergence_c(i),headc(:,:,i),cdum) ; enddo
        endif
        call timer_stop('SLF divergence',times=time_gamma)
      endif

c
c     Prepare Coulomb matrix for exchange calculation (->coulomb0)
c     - case SX and COSX: calculate static screened interaction
c     - if needed, multiply with the inverse of the overlap matrix for current k point
# define NC Ncol(1,ngptm(ikpt))
      if(evalx) then
        call timer_start('SLF exchange prep')
        Nfence(coulomb0)
        Obegin
        if(any(job1%type==[J_SX,J_COSX])) then
          allocate(coulmat(nbasm(ikpt),nbasm(ikpt)))
          call unitarytrafo ( coulmat , screen(:,1) , ctrafo , 2 )
          coulomb0 = packmat(coulmat)
          deallocate(coulmat)
        endif
        if(ikpt==1) then
          if(any(job1%type==[J_SX,J_COSX])) then ; divergence_x =   divergence_h(1) / vol ; divfac_x = trace(ihead(:,1))/3 / vol
          else                                   ; divergence_x = 4*pi * divergence / vol ; divfac_x = 4*pi / vol
          endif
        endif
        Oend
        Nfence(coulomb0)
        Mpi( call Mcast(divergence_x) ; call Mcast(divfac_x) )
        if(    fullpw .and. all(job1%type/=[J_SX,J_COSX]) .or.
     &    .not.fullpw .and. any(job1%type==[J_SX,J_COSX]) ) then
          allocate  ( eig(ngptm(ikpt)) )
          Nallocate ( eigv,(S_ ngptm(ikpt),ngptm(ikpt) S_) )
          Nallocate ( olap,(S_ ngptm(ikpt),ngptm(ikpt) S_) )
          Ocall olap_gptm(olap,ngptm(ikpt),ikpt,ikpt)        ; Nfence(olap)
          ! Calculate pseudoinverse of olap -> olap**(-1)
          call Mdiagonalize(Win(eigv),eig,olap)
          Nfence(eigv)
          Mpi( Ocall Msum(eigv,comm=Ocomm) )
          Nfence(eigv)
          if(any(eig< 0)) Error('Overlap matrix not positive definite.')
          if(any(eig==0)) Error('Overlap matrix singular.')
          where(eig<cutzero) eig = 0
          where(eig/=0)      eig = 1/eig
          do i = Nrange(1,ngptm(ikpt))
            eigv(:,i) = sqrt(eig(i)) * eigv(:,i)
          enddo
          Nfence(eigv)
          Nfence(olap) ; if(size(eigv(NC,:))>0) olap(:,NC) = matmac(eigv,eigv(NC,:))
          Nfence(olap)
          Ndeallocate ( eigv )
          ! Multiply coulomb0 with olap**(-1) from both sides
          Nallocate ( coulmat,(S_ nbasm(ikpt),nbasm(ikpt) S_) ) ! unpacked coulomb0
          ifO call p_unpackmat(coulmat,coulomb0)
# ifndef INV
          Rbegin
          rdum = sum( [ (abs(imag(coulmat(i,i))),i=1,nbasm(ikpt)) ] ) / nbasm(ikpt)
          if(rdum>1d-6) Error('Very large imaginary part on Coulomb diagonal: '//Chf(rdum,'F20.12'))
          if(rdum>1d-12) Warn('Large imaginary part on Coulomb diagonal: '//Chf(rdum,'F20.12'))
          Rend
# endif
# undef NC
# define NC Ncol(1,nbasm(ikpt))
          Nfence(coulmat) ; coulmat(nbasp+1:,NC) = matmat( olap , coulmat(nbasp+1:,NC) )
          Nfence(coulmat) ; coulmat(NC,nbasp+1:) = matmat( coulmat(NC,nbasp+1:) , olap )
          Nfence(coulmat)
          Rbegin
          rdum = 0
          do j = 1,nbasm(ikpt)
            do i = 1,j
              rdum = rdum + abs(coulmat(i,j)-MCONJG(coulmat(j,i)))
            enddo
          enddo
          rdum = rdum / (nbasm(ikpt)*(nbasm(ikpt)+1)/2)
          write(6,'(A,ES8.1)') 'Hermiticity of transformed Coulomb matrix:',rdum
          if(rdum>1d-6) then
            Warn('Large average deviation from Hermiticity in Coulomb: '//Chf(rdum,'F11.8')//'. CUTZERO might help.')
          endif
          Rend
          deallocate (eig)
          Ndeallocate(olap)
          if(enable_multipole) then ! copy coulmat to coulomb0, coulmat is redefined below
            Nfence(coulomb0)
            Ocall p_packmat(coulomb0,coulmat,.false.)
            Nfence(coulomb0)
            Ndeallocate(coulmat)
          endif
        else if(.not.enable_multipole) then ! unpack coulomb0 -> coulmat
          Nallocate0 ( coulmat,(S_ nbasm(ikpt),nbasm(ikpt) S_) )
          ifO call p_unpackmat(coulmat,coulomb0)
          Nfence(coulmat)
        endif
        if(enable_multipole) then ! coulomb0 -> coulloc, coulmat
          if(ikpt==1) then
            if(fullpw) then ; Ocall coulomb_sphaverage(coulomb0,5)
            else            ; Ocall coulomb_sphaverage(coulomb0,1)
            endif
            Nfence(coulomb0)
          endif
          rdum = error_coulomb_multipole(coulomb0,nbasm(ikpt),0)
          Rwrite(6,'(A,ES8.1)') 'Residual interaction after transformation:',rdum
          if     (rdum>1d-6) then ; RError('Residual multipole interaction too large: '//Chf(rdum,'ES8.1'))
          else if(rdum>1d-8) then ; RWarn ('Residual multipole interaction large: '    //Chf(rdum,'ES8.1'))
          endif
          l = nbasm(ikpt) - nbasp + sum( [ (neq(itype) * (lcutm(itype)+1)**2, itype=1,ntype) ] )
          m = sum( [ (( neq(itype)*(2*l+1) * nindxm(l,itype)**2, l=0,lcutm(itype)),itype=1,ntype) ] )
          Nallocate( coulmat,(S_ l,l S_) )
          Nallocate( coulloc,(S_ m S_) )
          Ocall coulomb_multipole(coulloc,coulmat,coulomb0,nbasm(ikpt))
          Nfence(coulmat)
          Nfence(coulloc)
        endif
        call timer_stop('SLF exchange prep',times=time_exch)
      endif
# undef NC

c
c     Define  W^c = W - v  (->screen,ihead/w,divergence_h/c) in order to remove the divergence
      if(evalc.or.job1%type==J_COSX) then
        NFENCE
        if(nfreq>0) then
          do i = 1,nbas
            ifO screen(i*(i+1)/2,:) = screen(i*(i+1)/2,:) - coul(i)
          enddo
          if(ikpt==1) then
            ihead(1,:) = ihead(1,:) - 4*pi
            ihead(3,:) = ihead(3,:) - 4*pi
            ihead(6,:) = ihead(6,:) - 4*pi
          endif
        endif
        if(oselfc/=1.and.nfreqc>0) then
          do i = 1,nbas
            ifO screenc(i,i,:) = screenc(i,i,:) - coul(i)
          enddo
          if(ikpt==1) then
            iheadc(1,1,:) = iheadc(1,1,:) - 4*pi
            iheadc(2,2,:) = iheadc(2,2,:) - 4*pi
            iheadc(3,3,:) = iheadc(3,3,:) - 4*pi
          endif
        endif
        if(ikpt==1) then
          if(nfreq >0) divergence_h = divergence_h - 4*pi * divergence
          if(nfreqc>0) divergence_c = divergence_c - 4*pi * divergence
        endif
        NFENCE
      endif

      if(any(job1%type==[J_SX,J_COSX])) then
        if(nfreq/=1) Error('nfreq not set to 1. Frequency mesh read from spex.cor?')
        nfreq = 0 ! we only need the static screened interaction, which is stored in coulomb0; avoid loops over ifreq=1,..,nfreq
      endif

c
c     Test if screenc is symmetric (case inversion symmetry)
# ifdef INV
      if(evalc) then
        do ifreqc = Mrange1(nfreqc)
          if(sum(abs(screenc(:,:,ifreqc)-transpose(screenc(:,:,ifreqc))))>1d-10)
     &      Bug('Screened interaction not symmetric (screenc).')
        enddo
      endif
# endif

c
c     WBLOCK (see susceptibility): Decompose W into blocks (redefine screen,screenc,ctrafo)
      if(evalc) then
        if(allocated(cblock)) then
          call timer_start('SLF cblock')
          NFENCE ; Nfence(ctrafo)
          nwblock = maxval(cblock) ; allocate ( wblock(nwblock),pnt(nbas) )
          k       = 0
          do i = 1,nwblock
            wblock(i) = count(cblock==i)
            do j = 1,nbas
              if(cblock(j)==i) then
                k      = k + 1
                pnt(k) = j
              endif
            enddo
          enddo
          do i = Nrange(1,nbasm(ikpt))
            ctrafo(i,:) = ctrafo(i,pnt)
          enddo
          Nallocate0(screen_,(S_ sum(wblock**2),nfreq S_))        
          do i = Nrange(1,nfreq)
            call block_diag_packed(screen_(:,i),screen(:,i))
          enddo
          Nfence(screen_)
          do i = 1,nfreqc
            do j = Nrange(1,nbas) ; screenc(:,j,i) = screenc(pnt,j,i) ; enddo
          enddo
          NFENCE
          do i = 1,nfreqc
            do j = Nrange(1,nbas) ; screenc(j,:,i) = screenc(j,pnt,i) ; enddo
          enddo
          deallocate ( pnt )
          NFENCE
          do i = Nrange(1,nfreqc)
            call collapse_screen_cblock_c(screenc(:,:,i))
          enddo
          NFENCE ; Nfence(ctrafo)
          call timer_stop('SLF cblock')
        else
          Nallocate0(screen2_,(S_ nbas,nbas,nfreq S_) )
          do i = Nrange(1,nfreq)
            call p_unpackmat(screen2_(:,:,i),screen(:,i))
          enddo
          Nfence(screen2_)
        endif
      endif

      if(.not.evalc)     then ; nfreq1 = 0
      else if(oselfc/=1) then ; nfreq1 = 1     ! contour integration:   self-energy arguments  w + img*freq(:1)      (freq(1)=0)
      else                    ; nfreq1 = nfreq ! analytic continuation: self-energy arguments      img*freq(:nfreq)
      endif

      if(evalc.and.(freqint<=1.or.ozero)) then
        if(oselfc==1) then ; allocate ( wfreqintegral(0:3,nfreq,nfreq)  )
        else               ; allocate ( wfreqintegral(0:3,nfreq,nfreqr) )
        endif        
      endif

c
c     Remove Drude term, which is treated analytically
      if(ikpt==1.and.evalc) then
        ! Does Drude term exist?
        drude_separate = metal.and.plasma/=0
        ! Write head of screened interaction to spex.head
        Rif(wrtinfo) then
          i = fopen('spex.head',status='unknown',numbered=.true.)
          write(i,'(A)')         '# Head element of screened interaction on imaginary axis'
          if(drude_separate)
     &      write(i,'(A,F10.5)') '# Compare with -4pi * plasma**2 / (plasma**2 + freq**2),  plasma =',plasma
          write(i,'(F10.5,F15.10)') (freq(ifreq),dble((trace(ihead(:,ifreq)))/3),ifreq=1,nfreq)
          call fclose(i)
        endif
        ! If Padé is used, we cannot separate off the Drude term.
        drude_separate = drude_separate.and.freqint<2.and.(oselfc==1.or.nfreqc/=0)
        ! Check if separating off makes the head smoother
        if(drude_separate) then
          rdum  = 0
          rdum1 = 0
          do ifreq = 1,nfreq-1
            enediff = 4*pi * plasma**2 * ( ( freq(ifreq+1)**2 + plasma**2 )**(-1) -
     &                                     ( freq(ifreq)  **2 + plasma**2 )**(-1) )
            rdum    = rdum  + abs( trace(ihead(:,ifreq+1)-ihead(:,ifreq))/3           )
            rdum1   = rdum1 + abs( trace(ihead(:,ifreq+1)-ihead(:,ifreq))/3 + enediff )
          enddo
          Rwrite(6,'(A,F5.1,A,F5.1)') 'Variation in W(head):',rdum,' ->',rdum1
          if(rdum1>rdum) drude_separate = .false.
        endif
        ! Remove Drude term from head and Gamma divergence
        if(drude_separate) then
          Rwrite(6,'(A)') 'Drude term is separated off and treated analytically.'
          do ifreq = 1,nfreq
            rdum                = -plasma**2 / ( plasma**2 + freq(ifreq)**2 )
            ihead(:,ifreq)      = ihead(:,ifreq) - packmat ( 4*pi * rdum * identity(3) )
            divergence_h(ifreq) = divergence_h(ifreq) - divergence * 4*pi * rdum
          enddo
        endif
      endif

      Rbegin
      write(6,'(/A'NoA) 'Add contribution to self-energy'
      Rend

c
c     Pade coefficients from head (-> cpade_head)
      if(ikpt==1.and.ozero.and.evalc.and.oselfc/=1.and.nfreqc==0) then
        allocate ( cpade_head(nfreq+1,3,3) )
        cpade_head = 0
        k          = 0
        do j = 1,3
          do i = 1,j
            k = k + 1
            if(all(abs(ihead(k,:))>1d-10)) then
                       call pade_init(cpade_head(:,i,j),img*freq,       ihead(k,:) *(1d0,0d0),nfreq,-1)
              if(i/=j) call pade_init(cpade_head(:,j,i),img*freq,MCONJG(ihead(k,:))*(1d0,0d0),nfreq,-1)
            endif
          enddo
        enddo
      endif

c
c     Distribute blocks over processes (Mover/=0)
# ifdef MPI
      Mrank1 = Mrank
      Mcolor = 0
      Msub   = 0
      if(Mover>0) then
        call Mdistribute
#   ifndef old_mpiblk
        if(Mcolor/=0) then          
          Mrank1 = -1             ! avoid printing progress bar (and incrementation of iblock)
          Rwrite(6,'(A'NoA) '...' ! just three dots instead of progress bar
          call Mtask_init(nblock,maxval(Msub),Mfile,Mrank,Mcomm)
        endif
#   endif
        if(Mcolor/=0) call begin_split(Mcolor)
      endif
# endif

c
c
c     LOOP OVER BLOCKS   ( job1%full = .true.  : block matrices of full self-energy matrix;
c                          job1%full = .false. : dummy definition: each state is a "block" of its own )

      ikptq  = 0
      ispin  = 0
      iblock = 0

      do ! start loop blocks

c       Next block
        Mpi( if(Mrank1>=0) then ) ! for MPIBLK (then Mrank1=-1), iblock is assigned by Mtask below
        iblock = iblock + 1       ! take next block
        if(iblock>nblock) exit    ! exit if all done
        Mpi( endif )

c       Progress bar
        Mpi( if(Mrank1==0) then ) ! no progress bar for MPIBLK
        do i = 1,(49*iblock)/nblock-(49*(iblock-1))/nblock
          write(6,'(''.'''NoA)
        enddo
        Mpi( endif )

c       MPIBLK
# ifdef MPI
#   ifdef old_mpiblk
        if(Msub(iblock)/=Mcolor) cycle ! cycle here if current block is calculated by different subcommunicator
#   else
        if(Mcolor>0) then
          Rcall Mtask_indx(i,0,0,Mfile)
          call Mcast(i)
          if(i>nblock) exit
          iblock = Mblock(i)
        endif
#   endif
# endif
        
        ib = block(1,iblock) ! first band of block        

c       COPYDEG: If a state degenerate with the present one has been calculated already, we cycle here.
c       The self-energy is simply copied (after block loop).
# ifdef COPYDEG
        if(.not.job1%full) then
          if(any( [ (job1%kpt(i)==job1%kpt(ib) .and. job1%spin(i)==job1%spin(ib) .and.
     &              same_eigenspace(job1%band(ib),job1%band(i),job1%kpt(i),job1%spin(i)), i=1,ib-1) ] )) cycle
        endif
# endif       

        newkptq  = job1%kpt(ib)/=ikptq.or.job1%spin(ib)/=ispin
        ibandq   = job1%band(ib)
        ikptq    = job1%kpt(ib)
        ispin    = job1%spin(ib)

        iselfc   = sum( sizeblock(:iblock-1)**2 )
        iselfx   = sum( sizeblock(:iblock-1)*(sizeblock(:iblock-1)+1)/2 )
        sizeband = sizeblock(iblock) ! equals 1 for diag. and equals size(band) for FULL

c
c       Get irreducible k points kpt1(:nkpt1) wrt current q point ikptq (EIBZ)
        call getkpt1(kpt1,nkpt1,nkpts,sym1,nsym1,ikptq,ikpt,.false.)

c
c       Determine bands that are contained in current block (->band)
        if(job1%full) then
          ! according to block (full self-energy calculations)
          allocate ( band(sizeband) )
          band = job1%band(block(:sizeband,iblock))
        else
          ! Determine first and last state of degenerate subspace (only diagonal self-energy elements)
          call getdeg(deg0,deg1,ibandq,ikptq,ispin)
          ndeg = deg1 - deg0 + 1
          ! double spin degeneracies (SOC+inv.sym.) are not due to spatial symmetry and don't have to be averaged over: leave out.
          if(l_soc.and.invsym/=0.and.ndeg==2) then
            ndeg = 1
            deg1 = deg0
          endif
          ! define array band accordingly
          if(ndeg*nkpt1>=count(kptp==ikpt) Load(.and..false.) ) then ! in this case summing over all k points is faster (disabled for LOAD)
            allocate ( band(1) )
            band = ibandq
            deg0 = ibandq
            deg1 = ibandq
            call getkpt1_fullBZ(kpt1,nkpt1,nkpts,sym1,nsym1,ikpt)
          else
            allocate ( band(ndeg) )
            band = [ (i,i=deg0,deg1) ]
          endif
        endif

c
c       Define frequency arguments w of self-energy SIGMA(w)  (-> freqr1)
        if(evalc) then
          if(oselfc==1.or.oselfc==4) then
            freqr1 = freqr
          else
            eneq = ene(ibandq,ikptq,ispin)
            if(job1%full)
     &        Error('Full self-energy calculations only implemented for CONTINUE and CONTOUR [{...}] <...>')
            if(eneq<=efermi) then ; freqr1 = eneq - efermi - freqr
            else                  ; freqr1 = eneq - efermi + freqr
            endif
          endif
          allocate(selfc_hlp(sizeband,sizeband,size(selfc,2)))
          selfc_hlp = 0
        endif

c
c       Initialize addition of equivalent k points (define contracted irreps -> irrep_contr/t)
        if(use_sym.and.job1%full) then!.and.any(nkpts(:nkpt1)/=1)) then
          call timer_start('SLF equivalent prep')
          ndeg = deg(ibandq,ikptq,ispin) - ibandq + 1 ! dimension of irrep
          if(.not.trsoff) then
            allocate ( irrep_contt(ndeg**2,ndeg**2,size(band)/ndeg,size(band)/ndeg) )
            irrep_contt = 0
          endif
          allocate ( irrep_contr(ndeg**2,ndeg**2,size(band)/ndeg,size(band)/ndeg) )
          irrep_contr = 0
          do i = 1,size(band)/ndeg
            do j = 1,size(band)/ndeg
              isub = psub(block((i-1)*ndeg+1,iblock))
              jsub = psub(block((j-1)*ndeg+1,iblock))
              do k = 1,ndeg ; do l = 1,ndeg ; do m = 1,ndeg ; do p = 1,ndeg
                deg0 = k+(p-1)*ndeg
                deg1 = l+(m-1)*ndeg
                do isym = 1,nsym1
                  cdum = conjg ( irrep_sub(l,k,isub,sym1(isym)) ) * irrep_sub(m,p,jsub,sym1(isym))
                  if(sym1(isym)>nsymt) then
                    irrep_contt(m+(l-1)*ndeg,k+(p-1)*ndeg,i,j) = irrep_contt(m+(l-1)*ndeg,k+(p-1)*ndeg,i,j) + conjg(cdum)
                  else
                    irrep_contr(l+(m-1)*ndeg,k+(p-1)*ndeg,i,j) = irrep_contr(l+(m-1)*ndeg,k+(p-1)*ndeg,i,j) + cdum
                  endif
                enddo
              enddo ; enddo ; enddo ; enddo
            enddo
          enddo
          call timer_stop('SLF equivalent prep',times=time_equiv)
        endif

c
c       Gamma point
        if(ikpt==1) then
          call timer_start('SLF gamma')
          if(ozero) then
c           Calculate momentum matrix <nq|-i\nabla|n'q>
            Load( allocate ( cmtq(maxlmindx,ncent,nband(ikptq,ispin),nspin3) ) )
            Load( allocate ( cpwq(maxgpt,         nband(ikptq,ispin),nspin3) ) )
            Load( call read_wavef2([(i,i=1,nband(ikptq,ispin))],nband(ikptq,ispin),ikptq,ispin,cmtq,cpwq)  )
            allocate ( moment(nband(ikptq,ispin),size(band),3) )
            do i = Mrange1(size(band))
              call momentum_matrix(moment(:,i,:),[ikptq],1,ispin,ispin,band(i),band(i),1,nband(ikptq,ispin),.false.
     &                             MpiC(.false.) LoadC(.true.) )
            enddo
            MrangeDistr( moment (:, McolD1(size(band),i) ,:) ,i )
            if(evalc.and.newkptq) then ! all diagonal elements needed for correlation
              if(allocated(moment_diag)) deallocate ( moment_diag )
              allocate ( moment_diag(nband(ikptq,ispin),3) )
              do iband = Mrange1(nband(ikptq,ispin))
                call momentum_matrix(moment_diag(iband,:),[ikptq],1,ispin,ispin,iband,iband,iband,iband,.false.
     &                               MpiC(.false.) LoadC(.true.) )
              enddo
              MrangeDistr( moment_diag ( McolD1(nband(ikptq,ispin),i) ,:) ,i)
            endif
          endif
          Rbegin
c         (A) Pole contribution 1/k**2 at k->0 (head element)
          ! Correlation
          if(evalc) then
            if(freqint>=2.or.nfreqc==0) call pade_init(cpade,img*freq,divergence_h*(1d0,0d0),nfreq,-1)
            if(freqint>=2) call pade_poles(pole,resid,npole1,img*freq,divergence_h*(1d0,0d0),cpade,nfreq,smooth(2),-1,.false.)
          endif
          do ib = 1,sizeband
            eneq = ene(band(ib),ikptq,ispin) - efermi
            if(ologdiv) then
              if(allocated(moment)) then
                mom = moment(band(ib),band(ib),:)
              else
                call momentum_matrix(mom,[ikptq],1,ispin,ispin,band(ib),band(ib),band(ib),band(ib),.false. MpiC(.false.) )
              endif
              logdiv1 = logdiv(eneq,dble(mom))
            endif
            ! Correlation
            if(evalc) then
              do ifreqr = 1,nfreqr
                enediff = eneq - freqr1(ifreqr)
                if(ologdiv) logdiv2 = logdiv(enediff,dble(mom))
                ! Frequency convolution
                do ifreq = 1,nfreq1
                  if     (freqint<=1) then ; cdum = freqintegral(divergence_h,freq,nfreq,1,freq(ifreq)+img*enediff,0)
                  else if(freqint>=2) then ; cdum = freqintegral_poles(pole,resid,npole1,1,freq(ifreq)+img*enediff,0)
                  endif
                  if(drude_separate) cdum = cdum + img*plasma/2 * drudefunc(freq(ifreq)+img*enediff,1) * 4*pi * divergence
                  if(ologdiv)        cdum = real(cdum) * abs(2*logdiv2-1) + img * imag(cdum) ! 2*logdiv-1 switches from 1 (enediff<<0) to -1 (enediff>>0); abs() because, if 2*logdiv-1<0, the frequency integral has switched sign of the real part
                  if(oselfc==1) then ; selfc_hlp(ib,ib,ifreq)  = selfc_hlp(ib,ib,ifreq)  - cdum/vol ! The minus sign (-cdum/vol) comes from the prefactor i 
                  else               ; selfc_hlp(ib,ib,ifreqr) = selfc_hlp(ib,ib,ifreqr) - cdum/vol ! and the integration differential d(iw) = i*dw.
                  endif
                enddo
                ! Residues / eneq*enediff<=0: (un)occupied states require a positive (negative) argument to W, i.e., eneq and enediff must not have identical signs.
                if(oselfc/=1 .and. ( eneq*enediff<=0 .or. (ologdiv.and.abs(logdiv1-logdiv2)>1d-8) ) ) then
                  if(nfreqc==0) then
                    cdum = pade_func((1d0,0d0)*abs(enediff),img*freq,cpade,nfreq)
                  else
                    call getcoeff(coeff,abs(enediff),real(freqc),nfreqc)
                    cdum = sum(coeff*divergence_c)
                  endif
                  if(ologdiv) then
                    cdum = cdum * (logdiv1-logdiv2) ! also contains factor -1 for clockwise integration                    
                  else
                    if(enediff==0)       cdum =  cdum / 2 ! the other half comes from the integral over iw                    
                    if(freqr1(ifreqr)>0) cdum = -cdum     ! integration runs the perimeter of the rectangle clock-wise (factor -1)
                  endif
                  selfc_hlp(ib,ib,ifreqr) = selfc_hlp(ib,ib,ifreqr) - cdum/vol ! The minus sign comes from the prefactor i/(2*pi) times the factor 2*pi*i from the residues
                endif
              enddo
            endif
            ! Exchange
            if(evalx) then
              if(lomit(1)) then ; if(any(omit==band(ib))) cycle ; endif
              iself = iselfx + ib*(ib+1)/2
              if(ologdiv)     then ; selfx(iself) = selfx(iself) - divergence_x * logdiv1
              else if(eneq<0) then ; selfx(iself) = selfx(iself) - divergence_x
              endif                
              if(job1%type==J_COSX)  selfx(iself) = selfx(iself) + (divergence_x-4*pi/vol*divergence)/2
            endif
          enddo
          Rend
c         (B) Zero-order contribution 1/k**2 * k**2 at k->0 (OZERO)
          if(ozero) then
            ! Loop over bands
            if(oselfc==1) then ; allocate ( wfreqintegral1(0:3,nfreq,nfreq), wfreqintegral2(0:3,nfreq,nfreq) )
            else               ; allocate ( wfreqintegral1(0:3,nfreq,nfreqr),wfreqintegral2(0:3,nfreq,nfreqr) )
            endif
            allocate ( pm(3,3,size(band),size(band)),pv(3,size(band),size(band)),zm(3,3,size(band),size(band)) )
            ! Calculate second order term of <E phi_i | phi_n> = O(1) + O(k) + k*zm*k (->zm)
            zm = 0
            do j = 1,size(band)   ! phi_n (zm needed only if n is one of band)
              do i = 1,size(band) ! phi_i
                ldum = same_eigenspace( band(i) , band(j) , ikptq,ispin ,edeg)
                do iband = 1,nband(ikptq,ispin)
                  cmat = reshape ( [ ((moment(iband,i,k)*MCONJG(moment(iband,j,l)),k=1,3),l=1,3) ] , [3,3] )
                  if(ldum) then ! band(i), band(j) degenerate
                    if(.not.same_eigenspace( iband , band(j) , ikptq,ispin ,edeg)) ! contributes only if iband not degenerate to band(i),band(j)
     &                zm(:,:,i,j) = zm(:,:,i,j) - cmat /   (ene(band(j),ikptq,ispin) - ene(iband, ikptq,ispin))**2 / 2
                  else          ! band(i), band(j) not degenerate                     
                    if(same_eigenspace( iband , band(j) , ikptq,ispin ,edeg)) then ! contribution if iband, band(j) degenerate
                      zm(:,:,i,j) = zm(:,:,i,j) - cmat /   (ene(band(j),ikptq,ispin) - ene(band(i),ikptq,ispin))**2
                    else                                                           ! contribution if iband, band(j) not degenerate
                      zm(:,:,i,j) = zm(:,:,i,j) + cmat / ( (ene(band(j),ikptq,ispin) - ene(band(i),ikptq,ispin)) *
     &                                                     (ene(band(j),ikptq,ispin) - ene(iband, ikptq,ispin)) )
                    endif
                  endif
                enddo
              enddo
            enddo
            if(evalc) then ; nbnd = nband(ikptq,ispin)
            else           ; nbnd = bando
            endif
            m    = 0
            deg0 = 1
            do while(deg0<=nbnd)
              enek = ene(deg0,ikptq,ispin) - efermi
              deg1 = deg0
              do while(deg1<nbnd)
                if(abs(ene(deg1+1,ikptq,ispin)-ene(deg0,ikptq,ispin))>edeg) exit
                deg1 = deg1 + 1
              enddo
              if(deg1<deg0) Bug('Count error.')
              if(deg1>nbnd) Bug('deg1 exceeded upper bound.')
              ifMODP(m)
              if(evalc) then
                ! Initialize frequency convolutions
                if(oselfc==1) then ; allocate ( cvec(nfreq)  ) ; cvec = freq + img*enek
                else               ; allocate ( cvec(nfreqr) ) ; cvec = img*(enek-freqr1)
                endif
                call freqintegral_init(wfreqintegral1,freq,nfreq,cvec-1d-4,size(cvec))
                call freqintegral_init(wfreqintegral2,freq,nfreq,cvec+1d-4,size(cvec))
                call freqintegral_init(wfreqintegral, freq,nfreq,cvec,     size(cvec))
                deallocate ( cvec )
                wfreqintegral2 =   wfreqintegral2 + wfreqintegral1
                wfreqintegral1 = ( wfreqintegral2 - wfreqintegral1*2 ) / 2d-4
                wfreqintegral2 = ( wfreqintegral2 - wfreqintegral *2 ) / 1d-8
                ! Frequency convolutions (->hw)
                do j = 1,3
                  do i = 1,3
                    k = min(i,j) + max(i,j)*(max(i,j)-1)/2
                    call spline_init(aspline, freq, ihead(k,:), nfreq)
                    if(i>j) aspline = MCONJG(aspline)
                    do ifreq = 1,size(selfc,2)
                      hw(i,j,ifreq,0) = sum(wfreqintegral (:,:,ifreq)*aspline) / ( 2*pi*img)
                      hw(i,j,ifreq,1) = sum(wfreqintegral1(:,:,ifreq)*aspline) / ( 2*pi)
                      hw(i,j,ifreq,2) = sum(wfreqintegral2(:,:,ifreq)*aspline) / (-2*pi*img)
                    enddo
                  enddo
                enddo
                if(drude_separate) then
                  if(oselfc==1) then ; allocate ( cvec(nfreq)  ) ; cvec = freq + img* enek
                  else               ; allocate ( cvec(nfreqr) ) ; cvec =        img*(enek-freqr1)
                  endif
                  do ifreq = 1,size(selfc,2)
                    hw(:,:,ifreq,0) = hw(:,:,ifreq,0) +     identity(3) * (2*pi*plasma) * drudefunc(cvec(ifreq),1) * img
                    hw(:,:,ifreq,1) = hw(:,:,ifreq,1) +     identity(3) * (2*pi*plasma) * drudefunc(cvec(ifreq),2)
                    hw(:,:,ifreq,2) = hw(:,:,ifreq,2) - 2 * identity(3) * (2*pi*plasma) * drudefunc(cvec(ifreq),3) * img
                  enddo
                  deallocate ( cvec )
                endif
                ! Calculate derivatives of KS energies (->dene, dene2=dene*dene^T, ddene)
# if defined(__GFORTRAN__) && __GNUC__ < 5
#   warning rewritten to avoid gfortran 4.9.1 bug
                dene  = [ ( sum ( [ (moment_diag(id,i),id=deg0,deg1) ] ) , i=1,3 ) ]
                dene  = dene / (deg1-deg0+1)
# else
                dene  = [ ( sum ( [ (moment_diag(id,i),id=deg0,deg1) ] ) , i=1,3 ) ] / (deg1-deg0+1) ! average over degenerate states
# endif
                dene2 = reshape ( [ ((dene(i)*dene(j),i=1,3),j=1,3) ] , [3,3] )
                do k = 1,size(band) ; band1 = band(k) ! ddene only needed if one of band(:) falls inside deg0..deg1
                  if(band1>=deg0.and.band1<=deg1) then
                    ddene = identity(3) + 2 * reshape ( (/ ((
     &                sum   ( moment(:nband(ikptq,ispin),k,i) *
     &                MCONJG (moment(:nband(ikptq,ispin),k,j)) / (ene(band1,ikptq,ispin)-ene(:nband(ikptq,ispin),ikptq,ispin)),
     &                abs(ene(band1,ikptq,ispin)-ene(:nband(ikptq,ispin),ikptq,ispin))>edeg ), i=1,3),j=1,3) /) , [3,3] )
                    exit
                  endif
                enddo
              endif
              ! Calculate k expansion of <E phi_i | phi_n> <phi_n | phi_j E> = O(1) + k*pv + k*pm*k (->pm,pv)
              do iband = deg0,deg1 ! phi_n
                allocate ( cvec(3) )
                do j = 1,size(band)
                  do i = 1,size(band)
                    cvec = 0
                    cmat = 0
                    if( (band(i)<deg0.or.band(i)>deg1) .and. (band(j)<deg0.or.band(j)>deg1) ) then
                      cmat = cmat + reshape ( [((moment(iband,i,k)*MCONJG(moment(iband,j,l)),k=1,3),l=1,3)], [3,3] )
     &                              / ( (ene(iband,ikptq,ispin) - ene(band(i),ikptq,ispin)) *
     &                                  (ene(iband,ikptq,ispin) - ene(band(j),ikptq,ispin)) )
                    endif
                    if(band(i)==iband) then
                      if(band(j)<deg0.or.band(j)>deg1)
     &                cvec = cvec + MCONJG(moment(iband,j,:)) / (ene(iband,ikptq,ispin) - ene(band(j),ikptq,ispin))
                      cmat = cmat + conjg(zm(:,:,j,i))
                    endif
                    if(band(j)==iband) then
                      if(band(i)<deg0.or.band(i)>deg1)
     &                cvec = cvec + moment(iband,i,:) / (ene(iband,ikptq,ispin) - ene(band(i),ikptq,ispin))
                      cmat = cmat + zm(:,:,i,j)
                    endif
                    pm(:,:,i,j) = cmat
                    pv(:,  i,j) = cvec
                  enddo
                enddo
                deallocate ( cvec )
                ! Exchange self-energy
                if(evalx.and.enek<0) then
                  iself = iselfx
                  do j = 1,sizeband
                    do i = 1,j
                      iself        = iself + 1
                      selfx(iself) = selfx(iself) - divfac_x/(3*nkpt) * ( pm(1,1,i,j)+pm(2,2,i,j)+pm(3,3,i,j) )
                    enddo
                  enddo
                endif
                ! Correlation self-energy
                if(evalc) then
                  do ifreq = 1,size(selfc,2)
# ifdef CMP_0405adj
                    Error('ddene undefined due to change in selfenergy!')
                    h0     =        hw(:,:,ifreq,0)
                    h2     = matmul(hw(:,:,ifreq,1),ddene) / 2 + matmul(hw(:,:,ifreq,2),dene2) / 2
# else
                    h0     = hw(:,:,ifreq,0)
                    h1     = hw(:,:,ifreq,1)
                    h2     = hw(:,:,ifreq,2)
# endif
                    if(oselfc/=1) then
                      enediff = enek - freqr1(ifreq)
                      if(enek*enediff<=0) then
                        p = 1
                        if(enediff==0) p = 2
                        if(enek>0)     p = -p
                        do i = 1,3 ; do j = 1,3
                          cdum = 0 ; cdum1 = 0 ; cdum2 = 0
                          if(nfreqc==0) then
                            if(any(cpade_head(:,i,j)/=0)) then
                              cdum  = pade_func((1d0,0d0)* abs(enediff),      img*freq,cpade_head(:,i,j),nfreq)
                              cdum1 = pade_func((1d0,0d0)*(abs(enediff)-1d-4),img*freq,cpade_head(:,i,j),nfreq)
                              cdum2 = pade_func((1d0,0d0)*(abs(enediff)+1d-4),img*freq,cpade_head(:,i,j),nfreq)
                            endif
                          else
                            call getcoeff(coeff,abs(enediff),     real(freqc),nfreqc) ; cdum  = sum(coeff*iheadc(i,j,:))
                            call getcoeff(coeff,abs(enediff-1d-4),real(freqc),nfreqc) ; cdum1 = sum(coeff*iheadc(i,j,:))
                            call getcoeff(coeff,abs(enediff+1d-4),real(freqc),nfreqc) ; cdum2 = sum(coeff*iheadc(i,j,:))
                          endif
                          cdum2   = cdum1 + cdum2
                          cdum1   = ( cdum2 - cdum1*2 ) / 2d-4
                          cdum2   = ( cdum2 - cdum *2 ) / 1d-8
                          h0(i,j) = h0(i,j) + cdum  / p
# ifndef CMP_0405adj
c                          h1(i,j) = h1(i,j) + cdum1 / p ! currently disabled because a Taylor expansion of W along the real axis is probably unreliable
c                          h2(i,j) = h2(i,j) + cdum2 / p
# endif
                        enddo ; enddo
                      endif
                    endif
                    do j = 1,sizeband
                      do i = 1,sizeband
                        cdum  = 0
# ifdef CMP_0405adj
                        if(i==j.and.iband==band(i))
     &                  cdum = cdum - ( h2(1,1)+h2(2,2)+h2(3,3) ) / (3*vol*nkpt)
                        cdum = cdum - ( sum(h0*pm(:,:,i,j))     ) / (3*vol*nkpt)
# else
                        cmat = reshape ( [((dene(k)*pv(l,i,j),k=1,3),l=1,3)], [3,3] )
                        if(i==j.and.iband==band(i))
     &                  cdum = cdum - ( avg(h1,(1d0,0d0)*ddene) + avg(h2,(1d0,0d0)*dene2) ) / (vol*nkpt*2)
                        cdum = cdum - ( avg(h1,cmat)            + avg(h0,pm(:,:,i,j))     ) / (vol*nkpt)
# endif
                        selfc_hlp(i,j,ifreq) = selfc_hlp(i,j,ifreq) + cdum
                      enddo
                    enddo
                  enddo
                endif
              enddo
              endMOD
              deg0 = deg1 + 1
            enddo
            deallocate ( pm,pv,zm,wfreqintegral1,wfreqintegral2 )
          endif
          call timer_stop('SLF gamma',times=time_gamma)
        endif

c
c       Define tetrahedron weights
        if(evalc.and.freqint/=3.and..not.allocated(wintgrc)) then
          if(oselfc==1) then
            ! copy from wintgr
            minb(1,:) = 1
            maxb(1,:) = bando
            allocate ( wintgrc(nkpti+nkpti2,bando,1,nspin1) )
            wintgrc(:nkpti,:,1,:) = wintgr(:nkpti,:,:)
            if(lkptadd) wintgrc(nkpti+1:,:,1,:) = wintgr(nkpt+1:nkpt+nkpti2,:,:)
          else
            call timer_start('SLF tetrahedron')
            ! determine minimum and maximum bands
            minb = 0
            maxb = 0
            do ifreqr = 1,nfreqr
              do s = 1,nspin1 ; if(oselfc/=4.and.s/=ispin) cycle
                ! (1) minimum
                i = 1
                do while(all(ene(i,:,s)-freqr1(ifreqr)<=efermi).and.i/=maxeband) ! all states below efermi(+freqr1) will get weight 1/nkpt
                  i = i + 1
                enddo
                j = i
                do k = 1,size(ene,2) ! be sure to catch all degenerate states
                  if(j<=nband(k,s)) i = min(i,deg(j,k,s))
                enddo
                minb(ifreqr,s) = i
                ! (2) maximum
                do while(any(ene(i,:,s)-freqr1(ifreqr)<efermi).and.i/=maxeband) ! all states above efermi(+freqr1) will get weight 0
                  i = i + 1
                enddo
                i = i - 1
                if(i>0) then
                  j = i
                  do k = 1,size(ene,2) ! be sure to catch all degenerate states
                    if(j<=nband(k,s)) then
                      l = deg(j,k,s) ; if(l<j) l = deg(l,k,s)
                      i = max(i,l)
                    endif
                  enddo
                endif
                maxb(ifreqr,s) = i
              enddo
            enddo
            ! calculate weights
            i = maxval(maxb-minb)+1
            if(oselfc/=4) then ; allocate ( wintgrc(nkpti+nkpti2,i,nfreqr,ispin:ispin) )
            else               ; allocate ( wintgrc(nkpti+nkpti2,i,nfreqr,nspin1) )
            endif
            if(i>0) then              
              wintgrc = 0
              do s = 1,nspin1 ; if(oselfc/=4.and.s/=ispin) cycle
                do ifreqr = 1,nfreqr
                  call tetrahedron_init(wintgrc(:nkpti,:,ifreqr,s),nkpti,minb(ifreqr,s),maxb(ifreqr,s),
     &                                  efermi+freqr1(ifreqr),s,.false.)
                  if(lkptadd)
     &              call tetrahedron_init(wintgrc(nkpti+1:,:,ifreqr,s),nkpti2,minb(ifreqr,s),maxb(ifreqr,s),
     &                                    efermi+freqr1(ifreqr),s,.true.)
                enddo
              enddo
            endif
            call timer_stop('SLF tetrahedron',times=time_tetra)
          endif
        endif

c       IBC: allocate auxiliary selfc_ibc array
        if(oibc/=0.and.evalc) then
          Mpi( Error('IBC&MPI not implemented') )
          allocate ( selfc_ibc(size(band),size(band),nfreq) )
          selfc_ibc = 0
        endif

c       MPI: read cmtq/cpwq
# ifdef LOAD
        Error('LOAD disabled. Please compile without -DLOAD.')
        if(allocated(cmtq)) then
          cmtq(:,:,:size(band),:) = cmtq(:,:,band,:)
          cpwq(:,  :size(band),:) = cpwq(:,  band,:)
          call reallocate ( cmtq,maxlmindx,ncent,size(band),nspin3 )
          call reallocate ( cpwq,maxgpt,         size(band),nspin3 )
        else
          allocate ( cmtq(maxlmindx,ncent,size(band),nspin3) )
          allocate ( cpwq(maxgpt,         size(band),nspin3) )
          call read_wavef2(band,size(band),ikptq,ispin,cmtq,cpwq)
        endif
        band = [ (i,i=1,size(band)) ]
# endif

        if(evalc) then ; nbnd = maxval ( [ (nband(kptsum(ikptq,kpt1(k)),ispin),k=1,nkpt1) ] )
        else           ; nbnd = bando
        endif
        nbnd0 = nbnd ! nbnd can later be smaller than nbnd0 if MEM separation is required

c
c       LOAD: allocate cmt/cpw
# ifdef LOAD
        if(pack/=bandup.or.newkptq) then                    ! MEM separation or newkptq
          newkptq = .true.                                  ! force execution of read_wavef below
          nbnd    = ceiling( (1d0*pack-1d-12) Mpi(/Msize) ) ! maximal number of bands per process
          kindx   = 0                                       ! prepare kindx and kpt2
          i       = 0
          do k = 1,nkpt1
            ikpt1 = kpt1(k)
            ikpt2 = kptsum(ikpt1,ikptq)
            if(storeibz) ikpt2 = kptp(ikpt2)
            if(kindx(ikpt2)==0) then ; i = i + 1 ; kindx(ikpt2) = i ; kpt2(i) = ikpt2 ; endif
            if(storeibz) kindx(kptsum(ikpt1,ikptq)) = kindx(ikpt2)
          enddo
          nkpt_2 = i
          if(associated(cmt)) deallocate(cmt,cpw)
          allocate ( cmt(maxlmindx,ncent,nbnd,nkpt_2,nspin3) )
          allocate ( cpw(maxgpt,         nbnd,nkpt_2,nspin3) )
        endif

        if(newkptq) call read_wavef0([(i,i=band1,band2)],kpt2(:nkpt_2),ispin,cmt,cpw,.true.)
# endif

c
c       Wave-function products at kpt1

        if(evalc.and.oibc/=0) then
          Error('IBC currently disabled')
          allocate ( cprod_ibc(nbasp,maxlmindx,size(band),nkpt1) )
          cprod_ibc = 0
        endif

c
c       Determine for which real frequencies {w} we must calculate matrices matc(w) = <nn''|W(w)|n''n'>, n''=iband  (->bfreqc))
        nfrq = 0
        if(nfreqc/=0) then
          call timer_start('SLF bfreqc')
          allocate(bfreqc(2,nbnd,nkpt1)) ! bfreqc(1:2,iu,k) is first and last freqc-index for band iu (at kpoint k)
          bfreqc = 0
          Mpi( m = -1 )
          do k = 1,nkpt1
            ikpt1  = kpt1(k)
            ikpt2  = kptsum(ikpt1,ikptq)
            ikpt2_ = kptp(ikpt2) ; if(ikpt2_>nkpt) ikpt2_ = ikpt2_ - nkpt + nkpti ! for wintgrc
            do iband = 1,min(nbnd,nband(ikpt2,ispin)) ; McycleP(m)
              enek    = ene(iband,ikpt2,ispin) - efermi
              bfreqc1 = nfreqc+1
              bfreqc2 = 0
              if(freqint==3) then ! FREQINT NONLIN
                call tetrahedron_nonlin(nfrq,frq,pfrq,enek,iband,ikpt2,ispin,.false.)
                f1 = frq(1)    + enek + 1d-12
                f2 = frq(nfrq) + enek - 1d-12
                do ifreqr = 1,nfreqr
                  emin = min(0d0,freqr1(ifreqr))
                  emax = max(0d0,freqr1(ifreqr))
                  sgn  = nint( -sign(1d0,freqr1(ifreqr)) )
                  if(emax-emin>1d-8.and.f1<emax.and.f2>emin) then
                    do ifreqc = 1,nfreqc-1
                      rdum  = sgn*real(freqc(ifreqc))   + freqr1(ifreqr)
                      rdum1 = sgn*real(freqc(ifreqc+1)) + freqr1(ifreqr)
                      if(min(rdum,rdum1)<min(emax,f2).and.max(rdum,rdum1)>max(emin,f1)) then
                        bfreqc1 = min(bfreqc1,ifreqc)
                        bfreqc2 = max(bfreqc2,ifreqc+1)
                      endif
                    enddo
                  endif
                enddo
              else ! FREQINT PADE/SPLINE
                do ifreqr = 1,nfreqr
                  enediff = enek - freqr1(ifreqr)
                  weight  = 0
                  if     (iband<=bando)                    weight = -wintgr(ikpt2,iband,ispin)
                  if     (iband< minb(ifreqr,ispin)) then; weight = weight + 1d0/nkpt
                  else if(iband<=maxb(ifreqr,ispin)) then; weight = weight + wintgrc(ikpt2_,iband-minb(ifreqr,ispin)+1,ifreqr,ispin)
                  endif
                  if(abs(weight)>1d-10) then
                    call getcoeff(coeff,abs(enediff),real(freqc),nfreqc)
                    do ifreqc = 1,nfreqc
                      if(coeff(ifreqc)/=0) then
                        bfreqc1 = min(bfreqc1,ifreqc)
                        bfreqc2 = max(bfreqc2,ifreqc)
                      endif
                    enddo
                  endif
                enddo
              endif
              if(bfreqc1<=bfreqc2) bfreqc(:,iband,k) = [bfreqc1,bfreqc2]
            enddo            
          enddo
          Mpi( call Msum(bfreqc) )
          where( bfreqc(1,:,:) == 0 ) bfreqc(1,:,:) = 1 ! empty ranges are defined as [0,1]
          allocate(bu(2,nfreqc,nkpt1))   ! bu(1:2,ifreqc,k) is first and last u-index, where frequency ifreqc is needed
          bu(1,:,:) = maxband + 1
          bu(2,:,:) = 0
          do k = 1,nkpt1
            ikpt1 = kpt1(k)
            ikpt2 = kptsum(ikpt1,ikptq)
            do iband = 1,min(nbnd,nband(ikpt2,ispin))
              bfreqc1                 = bfreqc(1,iband,k)
              bfreqc2                 = bfreqc(2,iband,k)
              bu(1,bfreqc1:bfreqc2,k) = min( iband , bu(1,bfreqc1:bfreqc2,k) )
              bu(2,bfreqc1:bfreqc2,k) = max( iband , bu(2,bfreqc1:bfreqc2,k) )
            enddo
          enddo
          where( bu(1,:,:) == maxband+1 ) bu(1,:,:) = 1 ! empty ranges are defined as [0,1]
          call timer_stop('SLF bfreqc')
        endif

c
c MPI work distribution by Mblocks: band(b1:b2), kpt1(k1:k2), [u1,u1+us,u1+2*us,etc.]
c 
c Macro COLOR_EQ(i,j) returns same "color" for all processes that work on the same ith and jth block (i,j=1,2,3)
c        
# ifdef MPI
#   define COLOR_EQ(i,j) bounds(3,i)+bounds(3,j)*bounds(4,i)
        call Mblocks(bounds,[size(band),nkpt1,nbnd],3,scaling_linear)       
        b1 = bounds(1,1)
        b2 = bounds(2,1)
        nb = b2 - b1 + 1
        k1 = bounds(1,2)
        k2 = bounds(2,2)
        nk = k2 - k1 + 1
        u1 = bounds(3,3) + 1
        us = bounds(4,3)
        nu = ( nbnd - u1 ) / us + 1
        if(any([nb,nk,nu]<=0)) then ! cycle here if current rank is taken idle for this round
          beginSplit(mpi_undefined)
          endSplit
          cycle
        else
          beginSplit(1)
        endif
# else
        b1 = 1 ; b2 = size(band) ; nb = b2 - b1 + 1
        k1 = 1 ; k2 = nkpt1      ; nk = k2 - k1 + 1
        u1 = 1 ; us = 1          ; nu = nbnd
# endif        

c NUPACK: maximal packet size
# define NUFACTOR 0.2d0
# define NUPACK   int(maxgpt*NUFACTOR)
        nu1 = max(1, min(nu,int(NUPACK/size(band))) ) ! u-dimension for help array cprod2

c
c       MEM work packages
        memo = (maxmem-mem) Mpi(/Nsize0) ! memo: physical memory available for current process
        do
          call memory_demand(memd)       ! memd: memory demand of current process
          ldum = memd<memo
          Mpi( call mpi_allreduce(memd<memo,ldum,1,mpi_logical,mpi_land,Mcomm,Merr) ) ! all processes should run the same number of MEM loops
          if(ldum) exit
          if(nu<=1) then
            call memory_demand(memd,output=.true.)
            Error('Memory demand exceeds available memory (see list in stdout).')
          endif
          nbnd = (nbnd+1) / 2 ; Mpi( nbnd = max(nbnd,us) ) ! halve u bands if memory not sufficient
          nu   = (nbnd-u1)/us + 1
          nu1  = (nu1+1) / 2
        enddo
        maxpack = max( maxpack , (nbnd0-1) / nbnd + 1 )
        nbnd1   = nbnd
        Mpi( nsplit = 0 )


        do ! MEM loop starts here

c
c       Calculate wavefunction products (->cprod)
        allocate ( cprod(nbasm(ikpt),nb*nu,k1:k2) ) ! main cprod array

        beginSplit( COLOR_EQ(1,2) )
        call timer_start('SLF wavefproducts1_mt')
        call wavefproducts1_mt(cprod,nbasm(ikpt),band(b1:b2),nb,ikptq,ispin,kpt1(k1:k2),nk,u1,nbnd,us)
        call timer_stop('SLF wavefproducts1_mt',times=time_mt)

        call timer_start('SLF wavefproducts1_pw')
        call wavefproducts1_pw(cprod,nbasm(ikpt),band(b1:b2),nb,ikptq,ispin,kpt1(k1:k2),nk,u1,nbnd,us)
        call timer_stop('SLF wavefproducts1_pw',times=time_pw)
        endSplit

c
c       Rotate to kptp(ikpt1), and transform to eigenbasis
        do k = k1,k2
          call timer_start('SLF trafo')
          if(lomit(3)) cprod(:,(omit-1)*nb+1,k) = 0
          ikpt1 = kpt1(k)
          ikpt2 = kptsum(ikpt1,ikptq)
          call mtrafo1 Inv(_r) (cprod(:,:,k),nbasm(ikpt),nb*nu,kptp(ikpt1),-symkpt(ikpt1),1)
          call timer_stop('SLF trafo',times=time_trafo)
        enddo

c
c --------------------------------
c       Exchange self-energy
c --------------------------------

        if(.not.evalx) goto 1 ! skip if HF not requested

        beginSplit( COLOR_EQ(2,3) )
        if(u1<=bando) then
          call timer_start('SLF exchange')
          nu0 = ( min(bando,nbnd) - u1 ) / us + 1 ! nu0: special nu for exchange (nu0/=nu for JOB=GW)          
          allocate( cprod2(nbasm(ikpt),nb*min(nu0,nu1)) ) ! help cprod2 array for dgemm
          allocate( matx(size(band),size(band)) )
          matx = 0          
          do k = k1,k2
            ikpt1 = kpt1(k)
            ikpt2 = kptsum(ikpt1,ikptq)
            do u2 = 0,nu0-1,nu1 ! use cprod2 help array in packets
              nu2 = min(nu0-u2,nu1) ! actual packet u2+1:u2+nu2
              call timer_start('SLF exchange-mat')
              if(enable_multipole) then
                call matmat_coulomb_multipole(cprod2,coulloc,coulmat,cprod(:,u2*nb+1:(u2+nu2)*nb,k),ikpt,nb*nu2)
              else
                cprod2(:,:nb*nu2) = MATMAT( coulmat, cprod(:,u2*nb+1:(u2+nu2)*nb,k) )
              endif
              call timer_stop('SLF exchange-mat')
              call timer_start('SLF exchange-dist')
# ifdef MPI
              do irank = 0,Msize-1
                if(irank==Mrank) then ; b3 = b1 ; b4 = b2 ; endif
                call Mcast(b3,rank=irank)
                call Mcast(b4,rank=irank) ; nb1 = b4 - b3 + 1
                if(irank==Mrank) then ;                                         cprod_ => cprod2(:,:nb1*nu2)
                else                  ; allocate(cprod3(nbasm(ikpt),nb1*nu2)) ; cprod_ => cprod3
                endif
                call Mcast(cprod_,rank=irank)
# else
#   define b3 b1
#   define b4 b2
#   define nb1 nb                
#   define cprod_ cprod2                
# endif
                do u = 0,nu2-1
                  iu    = u*nb1          ! offset for cprod_
                  iu2   = (u2+u)*nb      ! offset for cprod
                  iband = (u2+u)*us + u1 ; if(iband>bando) Bug('>') ; if(lomit(1)) then ; if(any(omit==iband)) cycle ; endif
                  wght  = nkpts(k) * wintgr(kptp(ikpt2),iband,ispin)
                  if(symkpt(ikpt1)<=nsymt) then
                    matx(b3:b4,b1:b2) = matx(b3:b4,b1:b2) + wght * MCONJG( macmat( cprod_(:,iu+1:iu+nb1), cprod(:,iu2+1:iu2+nb,k) ))
                  else
                    matx(b3:b4,b1:b2) = matx(b3:b4,b1:b2) + wght *         macmat( cprod_(:,iu+1:iu+nb1), cprod(:,iu2+1:iu2+nb,k) )
                  endif
                enddo
# ifdef MPI
                nullify(cprod_)
                if(irank/=Mrank) deallocate(cprod3)
              enddo
# else
#   undef b3
#   undef b4
#   undef nb1
#   undef cprod_
# endif
              call timer_stop('SLF exchange-dist')
            enddo
          enddo
          deallocate(cprod2)
          call timer_stop('SLF exchange',times=time_exch)

          call timer_start('SLF equiv')
          ! Add contribution of equivalent k points
          call equivalent_kpoints(matx)
          ! Add to selfx
          iself = iselfx
          do j = 1,sizeband
            do i = 1,j
              iself        = iself + 1
              selfx(iself) = selfx(iself) - matx(i,j)
            enddo
          enddo
          
          deallocate(matx)
          call timer_stop('SLF equiv')
        endif
        endSplit

 1      if(.not.evalc) goto 2 ! skip the rest if GW not requested

c
c --------------------------
c       GW self-energy 
c --------------------------

c
c       Transform to eigenfunctions of Coulomb
        call timer_start('SLF trafo')
        do k = k1,k2
          do i = 1,size(cprod,2),NUPACK
            j                  = min(i+NUPACK-1,size(cprod,2))
            cprod(:nbas,i:j,k) = macmat ( ctrafo , cprod(:,i:j,k) )
          enddo
        enddo

c
c       Change leading order nbasm(ikpt)->nbas and define cprodp pointer accordingly
        if(nbas<nbasm(ikpt)) call collapse(cprod,nbas,nbasm(ikpt),nb*nu*nk)
        cprodp(1:nbas,1:nb*nu,k1:k2) => cprod
        call timer_stop('SLF trafo',times=time_trafo)        

c
c       MPI: Determine maximal number of equivalent k points among processes (to avoid deadlock) (->nk_dead)
# ifdef MPI
        k  = 0
        ka = k1
        do while(ka<=k2)
          ikpt2 = kptp(kptsum(kpt1(ka),ikptq))
          kb    = ka
          do while(kb<k2)
            if(kptp(kptsum(kpt1(kb+1),ikptq))/=ikpt2) exit
            kb = kb + 1
          enddo
          k  = k  + 1
          ka = kb + 1
        enddo
        call mpi_allreduce(k,nk_dead,1,mpi_integer,mpi_max,Mcomm,Merr)
# endif

c
c       Loop over the equivalent k points of EIBZ (Green function)
c       (We will sum up the matrices mat(c) for ikpt1=kpt1(ka:kb) because they all lead to symmetry-equiv. ikpt2, which makes the frequ. convolution identical.)
        ka = k1
        do while(ka<=k2)
          ikpt2 = kptp(kptsum(kpt1(ka),ikptq))
          kb    = ka          
          do while(kb<k2)
            if(kptp(kptsum(kpt1(kb+1),ikptq))/=ikpt2) exit
            kb = kb + 1
          enddo
          ikpt2_ = ikpt2 ; if(ikpt2>nkpt) ikpt2_ = ikpt2 - nkpt + nkpti ! for wintgrc

          Mpi( nk_dead = nk_dead - 1 )

c
c -------------------------------------------------
c         (A) iw integral (and Padé residues)
c -------------------------------------------------

          beginSplit( COLOR_EQ(2,3) )

          ! If b-range MPI-divided, define new u1_, us_, nu_ for after matrix-vector product
# ifdef MPI
          u1_ = u1 + us * Mrank
          us_ = us * Msize
          nu_ = ( nbnd - u1_ ) / us_ + 1 
# else
#   define u1_ u1
#   define us_ us
#   define nu_ nu
# endif

          allocate( mat_(size(band),size(band),0:nu_-1,nfreq), matx(size(band),size(band)) )
          mat_ = 0

c
c         Vector-matrix-vector products: mat_=cprod(i)*screen*cprod(j)
          call timer_start('SLF mat')
          allocate ( cprod2(nbas,nb*nu1) ) ! help cprod2 array for dgemm
          do ifreq = 1,nfreq
            do k = ka,kb ! we can sum mat over k because all ikpt2 are symmetry-equivalent
              ikpt1 = kpt1(k)
              do u2 = 0,nu-1,nu1 ! use cprod2 help array in packets
                nu2 = min(nu-u2,nu1) ! actual packet u2+1:u2+nu2
                if(allocated(cblock)) then
                  call matmat_cblock_ ifInv(r,c) (cprod2,screen_(:,ifreq),cprodp(:,u2*nb+1:(u2+nu2)*nb,k),nb*nu2)
                else                  
                  cprod2(:,:nb*nu2) = MATMAT( screen2_(:,:,ifreq) , cprodp(:,u2*nb+1:(u2+nu2)*nb,k) )
                endif
# ifdef MPI
                do irank = 0,Msize-1
                  if(irank==Mrank) then ; b3 = b1 ; b4 = b2 ; endif
                  call Mcast(b3,rank=irank)
                  call Mcast(b4,rank=irank) ; nb1 = b4 - b3 + 1
                  if(irank==Mrank) then ;                                  cprod_ => cprod2(:,:nb1*nu2)
                  else                  ; allocate(cprod3(nbas,nb1*nu2)) ; cprod_ => cprod3
                  endif
                  call Mcast(cprod_,rank=irank)
# else
#   define b3 b1
#   define b4 b2
#   define nb1 nb
#   define cprod_ cprod2
# endif
                  do u = 0,nu2-1
                    matx  = 0
                    iu    = u*nb1          ! offset for cprod_
                    iu2   = (u2+u)*nb      ! offset for cprodp
                    iband = (u2+u)*us + u1 ! band index
                    if(iband>nband(ikpt2,ispin)) exit ; if(lomit(2)) then ; if(any(omit==iband)) cycle ; endif
                    if(symkpt(ikpt1)<=nsymt) then
                      matx(b3:b4,b1:b2) = nkpts(k) * MCONJG( macmat( cprod_(:,iu+1:iu+nb1), cprodp(:,iu2+1:iu2+nb,k) ) )
                    else
                      matx(b3:b4,b1:b2) = nkpts(k) *         macmat( cprod_(:,iu+1:iu+nb1), cprodp(:,iu2+1:iu2+nb,k) )
                    endif
# ifdef MPI
                    urank = mod(u2+u,Msize) ! rank that will calculate current band
                    call Msum(matx,rank=urank)
                    if(Mrank==urank) mat_(:,:,(u2+u)/Msize,ifreq) = mat_(:,:,(u2+u)/Msize,ifreq) + matx
                  enddo ! u
                  nullify(cprod_)
                  if(irank/=Mrank) deallocate(cprod3)                  
                enddo ! irank
# else
                    mat_(:,:,u2+u,ifreq) = mat_(:,:,u2+u,ifreq) + matx
                  enddo ! u
#   undef b3
#   undef b4
#   undef nb1
#   undef cprod_
# endif
              enddo ! u2           
            enddo ! k
          enddo ! ifreq
          deallocate(cprod2,matx)
          call timer_stop('SLF mat',times=time_mat)

          endSplit

c
c         Loop over bands (Green function)

          Mpi( if(mpisym) call mpi_allreduce(nu_,nu_dead,1,mpi_integer,mpi_max,Mcomm,Merr) )

          allocate(mat(size(band),size(band),nfreq))
          mat = 0
          
          do u = 0,nu_-1            

            iband = u*us_ + u1_
            if(iband>nband(ikpt2,ispin)) then
              Mpi( if(mpisym) then ; beginSplit(mpi_undefined) ; endSplit ; endif ) ! avoid deadlock
              cycle
            endif
            
            mat = mat + mat_(:,:,u,:)
            call getdeg(deg0,deg1,iband,ikpt2,ispin)
            if(u+1<nu_.and.(u+1)*us_+u1_<=deg1) then
              Mpi( if(mpisym) then ; beginSplit(mpi_undefined) ; endSplit ; endif ) ! avoid deadlock
              cycle ! simply sum mat(:,:,:) if next state is degenerate with current
            endif

#ifdef MPI
            if(mpisym) then
              beginSplit((deg0-1)*(nkpti+nkpti2)+ikpt2_) ! processes working on same degenerate subspace are grouped in same subcommunicator
              call Msum(mat,rank=0)
              if(Mrank==0) then ; endSplit
              else              ; endSplit ; mat = 0 ; cycle
              endif
            endif
# endif

            enek = ene(iband,ikpt2,ispin) - efermi

c
c           Add contribution of equivalent kpoints
            do ifreq = 1,nfreq
              call equivalent_kpoints(mat(:,:,ifreq))
            enddo            

c
c           Multiply Green function and perform frequency integration (along imaginary axis)

            call timer_start('SLF freq prep')

c           Initialize frequency convolutions
            if(freqint<=1) then
              if(oselfc==1) then ; allocate ( cvec(nfreq)  ) ; cvec = freq + img*enek
              else               ; allocate ( cvec(nfreqr) ) ; cvec = img*(enek-freqr1)
              endif
              call freqintegral_init(wfreqintegral,freq,nfreq,cvec,size(cvec))
              deallocate ( cvec )
            endif

c           FREQINT NONLIN: Get tetrahedron weight function
            if(freqint==3) call tetrahedron_nonlin(nfrq,frq,pfrq,enek,iband,ikpt2,ispin,.false.)

            call timer_stop('SLF freq prep')

c           Loop over self-energy matrix elements i,j
            call timer_start('SLF freq')
            do j = 1,sizeband
              do i = 1,j             
                if(maxval(abs(mat(i,j,:)))>1d-10) then
                  if(freqint>=2.or.nfreqc==0) call pade_init(cpade,img*freq,(1d0,0d0)*mat(i,j,1:),nfreq,-1)
                  if(freqint>=2) then
                    call pade_poles(pole,resid,npole1,img*freq,(1d0,0d0)*mat(i,j,1:),cpade,nfreq,smooth(2),-1,.false.)
                  else
                    call spline_init(aspline,freq,mat(i,j,1:),nfreq)
                  endif

                  call iw_integral                              ! Imaginary-frequency integration
                  if(oselfc/=1.and.nfreqc==0) call add_residues ! Residues are calculated from Padé approximants of mat (can be called here)

                endif
              enddo
            enddo
            call timer_stop('SLF freq',times=time_freq)

            mat = 0

          enddo ! u-loop

# ifdef MPI
          if(mpisym) then
            do u = nu_,nu_dead-1 ! avoid deadlock
              beginSplit(mpi_undefined)
              endSplit
            enddo
          endif
# endif

          deallocate(mat,mat_)

c
c -----------------------------------------
c         (B) Real-frequency residues
c -----------------------------------------

          if(oselfc/=1.and.nfreqc>0) then ! Residues are calculated from W on the real axis (new uloop necessary)

            beginSplit( COLOR_EQ(2,3) )

            ! If b-range MPI-divided, define new u1_, us_, nu_ for after matrix-vector product
# ifdef MPI
            u1_ = u1 + us * Mrank
            us_ = us * Msize
            nu_ = ( nbnd - u1_ ) / us_ + 1 
# endif

            nfrq = sum(bfreqc(2,u1_:nbnd:us_,ka)-bfreqc(1,u1_:nbnd:us_,ka)+1)
            allocate( matc_(size(band),size(band),nfrq), matz(size(band),size(band)) )
            matc_ = 0

c
c           Vector-matrix-vector products: matc_=cprod(i)*screenc*cprod(j)
            call timer_start('SLF mat')
            allocate ( cprod2c(nbas,nb*nu1) ) ! help cprod2 array for dgemm
            do ifreqc = 1,nfreqc
              do k = ka,kb ! we can sum mat over k because all ikpt2 are symmetry-equivalent
                ikpt1 = kpt1(k)
                bu1   = max(    0 , ( bu(1,ifreqc,k) - u1 + us-1 ) / us )
                bu2   = min( nu-1 , ( bu(2,ifreqc,k) - u1        ) / us )
                do u2 = bu1,bu2,nu1 ! use cprod2 help array in packets
                  nu2 = min(bu2-u2+1,nu1) ! actual packet u2+1:u2+nu2
                  if(allocated(cblock)) then
                    call matmat_cblock_c(cprod2c,screenc(:,:,ifreqc),cprodp(:,u2*nb+1:(u2+nu2)*nb,k),nb*nu2)
                  else
                    cprod2c(:,:nb*nu2) = MATMAT( screenc(:,:,ifreqc) , cprodp(:,u2*nb+1:(u2+nu2)*nb,k) )
                  endif
# ifdef MPI
                  do irank = 0,Msize-1
                    if(irank==Mrank) then ; b3 = b1 ; b4 = b2 ; endif
                    call Mcast(b3,rank=irank)
                    call Mcast(b4,rank=irank) ; nb1 = b4 - b3 + 1
                    if(irank==Mrank) then ;                                   cprod_c => cprod2c(:,:nb1*nu2)
                    else                  ; allocate(cprod3c(nbas,nb1*nu2)) ; cprod_c => cprod3c
                    endif
                    call Mcast(cprod_c,rank=irank)
# else
#   define b3 b1
#   define b4 b2
#   define nb1 nb
#   define cprod_c cprod2c                    
# endif
                    do u = 0,nu2-1
                      iu    = u*nb1          ! offset for cprod_c
                      iu2   = (u2+u)*nb      ! offset for cprodp
                      iband = (u2+u)*us + u1 ! band index
                      if(iband>nband(ikpt2,ispin)) exit ; if(lomit(2)) then ; if(any(omit==iband)) cycle ; endif
                      if(ifreqc<bfreqc(1,iband,k).or.ifreqc>bfreqc(2,iband,k)) cycle
                      matz  = 0
                      if(symkpt(ikpt1)<=nsymt) then
                        matz(b3:b4,b1:b2) = nkpts(k) * transpose( macmat( cprodp(:,iu2+1:iu2+nb,k), cprod_c(:,iu+1:iu+nb1) ) )
                      else
                        matz(b1:b2,b3:b4) = nkpts(k) *            macmat( cprodp(:,iu2+1:iu2+nb,k), cprod_c(:,iu+1:iu+nb1) )
                      endif
                      Mpi( urank = mod(u2+u,Msize)    ) ! rank that will calculate current band
                      Mpi( call Msum(matz,rank=urank) )
                      Mpi( if(Mrank/=urank) cycle     )
                      iuf = sum(bfreqc(2,u1_:iband-1:us_,k)-bfreqc(1,u1_:iband-1:us_,k)+1) + ifreqc - bfreqc(1,iband,k) + 1 ! composite u-freqc-index for matc_
                      if(iuf<1.or.iuf>nfrq) Bug('iuf out of range 1-'//Chr(nfrq)//': '//Chr(iuf))
                      matc_(:,:,iuf) = matc_(:,:,iuf) + matz
                    enddo ! u
# ifdef MPI
                    nullify(cprod_c)
                    if(irank/=Mrank) deallocate(cprod3c)
                  enddo ! irank
# else
#   undef b3
#   undef b4
#   undef nb1
#   undef cprod_c            
# endif
                  
                enddo ! u2
              enddo ! k                
            enddo ! ifreqc
            deallocate(cprod2c,matz)
            call timer_stop('SLF mat',times=time_mat)

            endSplit

c
c           Loop over bands (Green function)

            Mpi( if(mpisym) call mpi_allreduce(nu_,nu_dead,1,mpi_integer,mpi_max,Mcomm,Merr) )

            iuf = 1
            do u = 0,nu_-1
              iband   = u*us_ + u1_
              bfreqc1 = bfreqc(1,iband,ka)
              bfreqc2 = bfreqc(2,iband,ka)
              if(iband>nband(ikpt2,ispin).or.bfreqc1>bfreqc2) then
                Mpi( if(mpisym) then ; beginSplit(mpi_undefined) ; endSplit ; endif ) ! avoid deadlock
                cycle
              endif
              
              if(.not.allocated(matc)) then ! if we cycle in degenerate subspace, matc is still allocated
                allocate(matc(size(band),size(band),bfreqc1:bfreqc2))
                matc = 0
              endif
              matc = matc + matc_(:,:,iuf:iuf+bfreqc2-bfreqc1)
              call getdeg(deg0,deg1,iband,ikpt2,ispin)
              if(u+1<nu_.and.(u+1)*us_+u1_<=deg1) then ! simply sum matc(:,:,:) if next state is degenerate with current
                iuf = iuf + bfreqc2 - bfreqc1 + 1
                Mpi( if(mpisym) then ; beginSplit(mpi_undefined) ; endSplit ; endif ) ! avoid deadlock
                cycle
              endif

# ifdef MPI
              if(mpisym) then
                beginSplit((deg0-1)*(nkpti+nkpti2)+ikpt2_) ! processes working on same degenerate subspace are grouped in same subcommunicator
                call Msum(matc,rank=0)
                if(Mrank==0) then ; endSplit
                else              ; endSplit ; deallocate(matc) ; cycle
                endif
              endif
# endif
              
              enek = ene(iband,ikpt2,ispin) - efermi

c
c             Add contribution of equivalent kpoints
              do ifreqc = bfreqc1,bfreqc2
                call equivalent_kpoints(matc(:,:,ifreqc))
              enddo

c
c             Multiply Green function and sum residues

              call timer_start('SLF freq residues')

c             FREQINT NONLIN: Get tetrahedron weight function
              if(freqint==3) call tetrahedron_nonlin(nfrq,frq,pfrq,enek,iband,ikpt2,ispin,.false.)

c             Loop over self-energy matrix elements i,j
              do j = 1,sizeband
                do i = 1,j             
                  if(maxval(abs(matc(i,j,:)))>1d-10) call add_residues
                enddo
              enddo
              call timer_stop('SLF freq residues',times=time_freq)
              
              iuf = iuf + bfreqc2 - bfreqc1 + 1
              deallocate(matc)

            enddo ! u-loop

            deallocate(matc_)

# ifdef MPI
            if(mpisym) then
              do u = nu_,nu_dead-1 ! avoid deadlock
                beginSplit(mpi_undefined)
                endSplit
              enddo
            endif
# endif

          endif
            
          ka = kb + 1
        enddo ! ka..kb loop

# ifdef MPI
        if(oselfc/=1.and.nfreqc>0) nk_dead = nk_dead * 2 ! there are two loops if CONTOUR (and nfreqc>0) is specified
        do k = 1,nk_dead ! avoid deadlock
          beginSplit(mpi_undefined) ; endSplit ! finish subcommunicator of matrix multiplications for iw integral (and CONTOUR)
          if(mpisym) then
            call mpi_allreduce(0,nu_dead,1,mpi_integer,mpi_max,Mcomm,Merr)
            do u = 0,nu_dead-1
              beginSplit(mpi_undefined) ; endSplit ! finish all subcommunicators opened above for deg. subspaces (only MPISYM)
            enddo
          endif
        enddo
# endif        

        nullify(cprodp)
 2      deallocate ( cprod )

c
c       MEM loop cycles here
        if(nbnd==nbnd0) exit ! exit if all MEM packages done, otherwise cycle
        u1 = nbnd + 1 Mpi(+bounds(3,3)) ! MPI: it can happen that u1>nbnd, then the process has to be taken idle
# ifdef MPI
        call mpi_allreduce(u1>nbnd0,ldum,1,mpi_logical,mpi_lor,Mcomm,Merr) ! check if any process has to be taken idle
        if(ldum) nsplit = nsplit + 1
        if(u1>nbnd0)  then ; beginSplit(mpi_undefined) ; exit ! take idle
        else if(ldum) then ; beginSplit(1)                    ! the others are grouped into a subcommunicator
        endif
# endif
        nbnd = min( nbnd + nbnd1 , nbnd0 )
        nu   = ( nbnd - u1 ) / us + 1 ; if(nu==0) Bug('Empty MEM work package: nu==0.')
        nu1  = min (nu , nu1 )

        enddo ! MEM loop
        Mpi2( do i = 1,nsplit ; endSplit ; enddo ) ! finish all subcommunicators opened above in MEM loop cycle

        endSplit
      
# ifndef MPI
#   undef u1_
#   undef us_
#   undef nu_
# endif

c
c       IBC: Add equivalent k points and add to ibc_selfc
        if(oibc/=0.and.evalc) then
          write(*,*) 'Add equivalent k points for IBC.'
          if(job1%full.and.any(nkpts(:nkpt1)/=1)) then
            Error('Not implemented job1%full and IBC.')
          endif
          if(.not.job1%full) then
            ! use great orthogonality theorem for diagonal elements
            do ifreq = 1,nfreq
              ibc_selfc(iselfc+1,ifreq) = sum( [ (selfc_ibc(i,i,ifreq),i=1,size(band)) ] ) / size(band)
            enddo
            write(*,'(3F20.10)') (freq(ifreq),-ibc_selfc(iselfc+1,ifreq),ifreq=1,nfreq)
            write(800,'(3F20.10)') (freq(ifreq),-ibc_selfc(iselfc+1,ifreq),ifreq=1,nfreq)
          endif
          deallocate ( selfc_ibc )
        endif

        if(allocated(wintgrc)) then
          if(oselfc==2.or.oselfc==3) deallocate ( wintgrc )
        endif
        if(allocated(irrep_contr)) deallocate ( irrep_contr )
        if(allocated(irrep_contt)) deallocate ( irrep_contt )

        Load(deallocate(cmtq,cpwq))

        deallocate ( band )
        if(allocated(moment)) deallocate ( moment )
        if(allocated(bfreqc)) deallocate ( bfreqc )
        if(allocated(bu))     deallocate ( bu     )

c       Update selfc with current self-energy
        if(evalc) then
          Mpi( call Msum(selfc_hlp,rank=0) )
          Rbegin
          do ifreq = 1,size(selfc,2)
            selfc(iselfc+1:iselfc+sizeband**2,ifreq) = selfc(iselfc+1:iselfc+sizeband**2,ifreq) +
     &                                                 reshape( selfc_hlp(:,:,ifreq), [sizeband**2] )
          enddo
          Rend
          deallocate(selfc_hlp)    
        endif

      enddo ! loop blocks

      call timer_start('SLF final')

# ifdef MPI
      call timer_start('SLF MPI idle')
      if(Mcolor/=0) then
        call end_split ! finish subcommunicator of MPIBLK
#   ifndef old_mpiblk
        call mpi_file_close(Mfile,Merr)
#   endif
      endif
      if(evalc) then
        Nfence(selfc)
        Ocall Msum(selfc,0,comm=Ocomm) ; MnoR(ifO selfc = 0 )
        Nfence(selfc)
      endif
      if(evalx) then
        call Msum(selfx,0) ; MnoR( selfx = 0 )
      endif
      call timer_stop('SLF MPI idle',time=time_idle)
      call mpi_reduce(time_idle,time_maxidle,1,mpi_real,mpi_max,0,Mcomm,Merr)
# endif

c     COPYDEG: Copy self-energy to degenerate states that come later in the JOB list.
# ifdef COPYDEG
      if(.not.job1%full) then
        Rbegin
        do i = 1,size(job1%band)          
          ibandq = job1%band(i)
          ikptq  = job1%kpt (i)
          ispin  = job1%spin(i)
          do j = i+1,size(job1%band)
            if(job1%kpt(j)==ikptq .and. job1%spin(j)==ispin .and. same_eigenspace(ibandq,job1%band(j),ikptq,ispin) ) then
              if(evalx) selfx(j)   = selfx(i)
              if(evalc) selfc(j,:) = selfc(i,:)
            endif
          enddo
        enddo
        Rend
        Mpi( if(evalc) then ; Nfence(selfc) ; endif )
      endif
# endif        

      if(allocated(moment_diag))   deallocate ( moment_diag )
      if(allocated(wfreqintegral)) deallocate ( wfreqintegral )
      if(allocated(cpade_head))    deallocate ( cpade_head )
      if(allocated(wintgrc))       deallocate ( wintgrc )
      Load( if(associated(cmt))    deallocate ( cmt,cpw ) )

      if(associated(coulmat))    tNdeallocate ( coulmat )
      if(associated(coulloc))    tNdeallocate ( coulloc )
      if(associated(screen_))    tNdeallocate ( screen_ )
      if(associated(screen2_))   tNdeallocate ( screen2_ )

      if(any(job1%type==[J_SX,J_COSX])) nfreq = size(freq) ! nfreq was temporarily set to 0

      if(allocated(wblock)) deallocate ( wblock )

      Mpi( m = maxpack ; call mpi_allreduce(m,maxpack,1,mpi_integer,mpi_max,Mcomm,Merr) )

      call timer_stop('SLF final')
      call timer_stop('Routine selfenergy SLF',time=time_tot)

      Rbegin
      write(6,'(/31X,A)')  'Timings'
      write(6,'(A,F13.5)') '  MT wave-func. products:',time_mt
      write(6,'(A,F13.5)') '  PW wave-func. products:',time_pw
      if(evalc)
     &write(6,'(A,F13.5)') '  Matrix-vector products:',time_mat
      write(6,'(A,F13.5)') '  Transformations:       ',time_trafo
      if(evalx)
     &write(6,'(A,F13.5)') '  Exchange self-energy:  ',time_exch
      if(job1%full)
     &write(6,'(A,F13.5)') '  Equivalent k points:   ',time_equiv
      if(evalc)
     &write(6,'(A,F13.5)') '  Frequency integration: ',time_freq
      if(oselfc/=0)
     &write(6,'(A,F13.5)') '  Tetrahedron weights:   ',time_tetra
      if(ikpt==1)
     &write(6,'(A,F13.5)') '  Gamma-point correction:',time_gamma
      if(evalc.and.oibc/=0)
     &write(6,'(A,F13.5)') '  IBC double-counting:   ',time_ibc
# ifdef MPI
      write(6,'(A,2F13.5)')'  MPI idle time:         ',time_idle,time_maxidle
# endif
      write(6,'(A,F13.5)') '  Other:                 ',time_tot-(time_mt+time_pw+time_mat+time_trafo+time_exch+time_equiv+time_freq+
     &                                                           time_tetra+time_gamma+time_ibc Mpi(+time_idle))
      write(6,'(A)')       '  ------------------------------------'
      write(6,'(A,F13.5)') '  Total:                 ',time_tot
      if(maxpack>1 Mpi(.or.Mover>0) ) write(6,*)
      if(maxpack>1) write(6,'(A)') '  Maximal number of MEM packets: '//Chr(maxpack)
# ifdef MPI
      if(Mover>0)   write(6,'(A)') '  Number of MPIBLK workers / processes / blocks: '//
     &                             Chr(maxval(Msub))//'/'//Chr(Msize)//'/'//Chr(nblock)
# endif
      if(job1%type==J_RPA) write(6,*)

      if(otimer==3) then
        call timer_print(['SLF '],select=.true.,title='TIMING SELFENERGY',      prnt=.true. andR )
        call timer_print(['WFP '],select=.true.,title='TIMING WAVEFPRODUCTS_MT',prnt=.true. andR )
      else if(otimer==4) then
        call timer_print(['SLF '],select=.true.,title='TIMING SELFENERGY',      prnt=.true., unit= ifMpi(Mrank,0) )
        call timer_print(['WFP '],select=.true.,title='TIMING WAVEFPRODUCTS_MT',prnt=.true., unit= ifMpi(Mrank,0) )
      endif

      Rend

      contains

c -----------------------------

      function drudefunc(c,n)
      implicit none
      complex_dp             :: drudefunc
      complex_dp, intent(in) :: c
      integer,    intent(in) :: n
      real_dp                :: c2
      c2 = imag(c)
      if (plasma==0) then ; drudefunc =   0
      else if(c2>0)  then ; drudefunc =   1/( img*plasma+c)**n
      else if(c2<0)  then ; drudefunc =   1/(-img*plasma+c)**n
      else                ; drudefunc = ( 1/( img*plasma+c)**n + 1/(-img*plasma+c)**n ) / 2
      endif
      end function drudefunc

c -----------------------------

      subroutine block_diag_packed(matout,matin)
      implicit none
      MCOMPLEX_dp, intent(in)  :: matin(nbas*(nbas+1)/2)
      MCOMPLEX_dp, intent(out) :: matout(*)
      MCOMPLEX_dp              :: mat1(nbas,nbas)
      integer                  :: i,j,jj,n
      call p_unpackmat(mat1,matin)
      do i = 1,nbas
        mat1(:,i) = mat1(pnt,i)
      enddo
      do i = 1,nbas
        mat1(i,:) = mat1(i,pnt)
      enddo
      j  = 0
      jj = 0
      do i = 1,nwblock
        n                    = wblock(i)
        matout(jj+1:jj+n**2) = reshape( mat1(j+1:j+n,j+1:j+n) , [n**2] )
        j                    = j  + n
        jj                   = jj + n**2
      enddo
      end subroutine block_diag_packed

      subroutine matvec_sym(res,mat,vec)
      implicit none
      MCOMPLEX_dp, intent(out) :: res(nbas)
      MCOMPLEX_dp, intent(in)  :: mat(nbas*(nbas+1)/2),vec(nbas)
      integer                  :: i,n0,n1,m0,m1
      res = 0
      n0  = 1
      m0  = 1
      do i = 1,nwblock
        n1         = n0 + wblock(i)                 - 1
        m1         = m0 + wblock(i)*(wblock(i)+1)/2 - 1
        res(n0:n1) = res(n0:n1) + matvec ( mat(m0:m1) , vec(n0:n1) )
        n0         = n1 + 1
        m0         = m1 + 1
      enddo
      end subroutine matvec_sym

      subroutine matvec_gen(res,mat,vec)
      implicit none
      complex_dp,  intent(out) :: res(nbas)
      complex_dp,  intent(in)  :: mat(nbas,nbas)
      MCOMPLEX_dp, intent(in)  :: vec(nbas)
      integer                  :: i,n0,n1
      res = 0
      n0  = 1
      do i = 1,nwblock
        n1         = n0 + wblock(i) - 1
        res(n0:n1) = res(n0:n1) + matvec ( mat(n0:n1,n0:n1) , vec(n0:n1) )
        n0         = n1 + 1
      enddo
      end subroutine matvec_gen

c ------------------------------

# ifdef MPI
c     Distribute job over processes. The current communicator (Mcomm) is split into subgroups s of P(s) processors.
c     ( sum(s) P(s) is the total number of processors. )
c
c     Assume computation time for each block b : T(b) = N(b) * [ M + K(b)*N'/P(s) ]
c     N(b) - number of bands in block b
c     K(b) - number of Green-function k-points for block b (nkpt1)
c     N'   - number of Green-function bands
c     M    - relative cost of N' independent term (Mover, keyword MPIBLK)
c
c     We determine
c     - the number of subgroups
c     - the number of processors of each subgroup P(s)
c     - the subgroup index S(b) for the block b
c     so as to minimize
c     R(S,P) = max(s) [ SUM[b,S(b)=s] T(b) ]
c
c     Output
c     - Mcolor    : Subgroup index of current rank
c     - Msub(i)   : Subgroup index of ith block                                                    ! for old_mpiblk
c     - Mblock(i) : ith block go be calculated, ordered according to decreasing computational time ! for new implementation
      subroutine Mdistribute
      implicit none
      integer :: s,proc(Msize),sub(nblock),i,j,iblock,nn,nk(nblock),nb(nblock)
      integer :: pnt(nblock)
      real_dp :: time(nblock),stime(Msize),mintime,maxtime
      if(Mover==0) then ! do not distribute jobs
        Mcolor = 0
        Msub   = 0
        return
      endif
      mintime = huge(0d0)
      ! information from blocks
      do iblock = 1,nblock
        ikptq      = job1%kpt(block(1,iblock)) ; call getkpt1(kpt1,nk(iblock),nkpts,sym1,nsym1,ikptq,ikpt,.false.) ! -> nk(iblock)
        nb(iblock) = sizeblock(iblock)
      enddo
      ! number of Green-function bands
      if(job1%type==J_GW) then ; nn = maxband
      else                     ; nn = bando
      endif
      ! loop over subgroups
      do s = 1,Msize
        ! distribute Msize processors over subgroups -> proc
        proc = Msize / s
        if(sum(proc(:s))/=Msize) then
          i        = Msize - Msize/s*s
          proc(:i) = proc(:i) + 1
        endif
        if(sum(proc(:s))/=Msize) Bug('Count error.')
        ! computation time for each block -> time(:)
        time = nb * ( Mover + 1d0 * nk * nn / Msize * s )
        call rorderp(pnt,time,nblock)
        ! distribute blocks over subgroups
        stime = -proc*1d-15
        do j = nblock,1,-1
          iblock      = pnt(j)
          i           = minloc(stime(:s),1)
          stime(i)    = stime(i) + nb(iblock) * ( Mover + 1d0 * nk(iblock) * nn / proc(i) )
          sub(iblock) = i
        enddo
        maxtime = maxval(stime(:s))
        if(maxtime<mintime) then
          mintime = maxtime
          Msub    = sub              ! for old_mpiblk
          Mblock  = pnt(nblock:1:-1) ! for new implementation
          Mcolor  = 0 ; do while(sum(proc(:Mcolor))<=Mrank) ; Mcolor = Mcolor + 1 ; enddo
        endif
      enddo
      end subroutine Mdistribute
# endif

c ------------------------------

c
c Returns memory demand for self-energy calculation in memd
c If output=true, a corresponding list is written in MB.
      subroutine memory_demand(memd,output)      
      implicit none
      real_dp, intent(out)          :: memd
      logical, intent(in), optional :: output
      logical                       :: output1
      real_dp                       :: memc
      integer                       :: nu0,u1_,us_,nu_,nfrq,k
      output1 = .false. ; if(present(output)) output1 = output
      if(output1) write(6,'(A)') 'MEMORY DEMAND:'
      memd = 0
      if(evalx.and.u1<=bando) then
        nu0  = (bando-u1)/us+1
        memd = MBYTES * nbasm(ikpt) * nb*min(nu0,nu1) Mpi(*2) ! for cprod2 (and cprod3)
        if(enable_multipole) then
          k    = sum( [ (neq(itype)*(lcut(itype)+1)**2,itype=1,ntype) ] )
          memd = memd + MBYTES * (ngptm(ikpt)+k) * nb*min(nu0,nu1) ! help array in matmat_coulomb_multipole
        endif
        if(output1) write(6,'(A)') 'for exchange: '//Chr(memd/1024d0**2)//' MB'
      endif
      if(evalc) then
        beginSplit( COLOR_EQ(2,3) )
        u1_  = u1 Mpi( + us * Mrank )
        us_  = us Mpi( * Msize      )
        nu_  = ( nbnd - u1_ ) / us_ + 1
        endSplit
        memc = MBYTES * size(band)**2*nu_*nfreq ! for mat_
     &       + MBYTES * nbas*nb*nu1 Mpi(*2)     ! for cprod2 (and cprod3)
        if(allocated(cblock)) then
          memc =memc + MBYTES * maxval(wblock(:nwblock)) * nb*nu1 ! for array slice in matmat of matmat_cblock
        endif
        memd = max( memd , memc )
        if(output1) write(6,'(A)') 'for correlation (iw integral): '//Chr(memc/1024d0**2)//' MB'
        if(oselfc/=1.and.nfreqc>0) then
          nfrq = maxval( [ (sum(bfreqc(2,u1_:nbnd:us_,k)-bfreqc(1,u1_:nbnd:us_,k)+1),k=k1,k2) ] )
          memc = 16d0 * size(band)**2*nfrq      ! for matc_
     &         + 16d0 * nbas*nb*nu1 Mpi(*2)     ! for cprod2c (and cprod3c)
          if(allocated(cblock)) then
            memc =memc + 16d0 * maxval(wblock(:nwblock)) * nb*nu1 ! for array slice in matmat of matmat_cblock
          endif
          memd = max( memd , memc )
          if(output1) write(6,'(A)') 'for correlation (residues): '//Chr(memc/1024d0**2)//' MB'
        endif
      endif
      end subroutine memory_demand
      
c ------------------------------

c Returns  1/(4*pi) * INT (kAk) (kBk) dÂ²k = [ tr(a)*tr(b) + 2*tr(ab) ] / 15
c - Simplifies to tr(ab) if a (or b) is a multiple of the identity matrix.
c - At least one of the matrices must be symmetric.
      function avg(a,b)
      implicit none
      complex_dp             :: avg
      complex_dp, intent(in) :: a(3,3),b(3,3)
      if(maxval(abs(a-transpose(a)))>1d-8.and.maxval(abs(b-transpose(b)))>1d-8) Error('both matrices asymmetric')
      avg = ( (a(1,1)+a(2,2)+a(3,3)) * (b(1,1)+b(2,2)+b(3,3)) + 2*sum(a*b) ) / 15
      end function avg

c ------------------------------

c
c Updates selfc with contribution from iw integration.
      subroutine iw_integral
      implicit none      

      if(freqint<=2) then ! FREQINT SPLINE/PADE

        do ifreq = 1,size(selfc,2)
          if(oselfc==1) then ; ifreqr = 1     ; rdum = freq(ifreq) ! the latter only for Pade integration
          else               ; ifreqr = ifreq ; rdum = 0           !              - " -
          endif
          if      (iband<minb(ifreqr,ispin)) then ; weight = 1d0/nkpt
          else if (iband>maxb(ifreqr,ispin)) then ; weight = 0
          else                                    ; weight = wintgrc(ikpt2_,iband-minb(ifreqr,ispin)+1,ifreqr,ispin)
          endif
          enediff = enek - freqr1(ifreqr)
          ! occ. states
          if(weight>1d-10) then
            if(freqint<=1) then
              if(enediff>0) then ; wfreq = conjg(wfreqintegral(:,:,ifreq)) ! protrudes above Fermi energy
              else               ; wfreq =       wfreqintegral(:,:,ifreq)  ! normal case
              endif
              cdum = sum(wfreq*aspline) / (2*pi*img)
            else
              cdum = freqintegral_poles(pole,resid,npole1,1,rdum-img*abs(enediff),0)
            endif
            if(rdum==0.and.enediff==0) cdum = cdum + mat(i,j,1)/2 ! add delta function (due to negative infinitesimal imaginary part)
            selfc_hlp(i,j,ifreq) = selfc_hlp(i,j,ifreq) - cdum * weight
            if(i/=j) then
# ifndef INV
              if(freqint<=1)then; cdum = sum(wfreq*conjg(aspline)) / (2*pi*img)
              else              ; cdum = freqintegral_poles(-conjg(pole),-conjg(resid),npole1,1,rdum-img*abs(enediff),0)
              endif
# endif
              selfc_hlp(j,i,ifreq) = selfc_hlp(j,i,ifreq) - cdum * weight
            endif
          endif
          ! unocc. states
          weight = 1d0/nkpt - weight
          if(weight>1d-10) then
            if(freqint<=1) then
              if(enediff<0) then ; wfreq = conjg(wfreqintegral(:,:,ifreq)) ! protrudes below Fermi energy
              else               ; wfreq =       wfreqintegral(:,:,ifreq)  ! normal case
              endif
              cdum = sum(wfreq*aspline) / (2*pi*img)
            else
              cdum = freqintegral_poles(pole,resid,npole1,1,rdum+img*abs(enediff),0)
            endif
            if(rdum==0.and.enediff==0) cdum = cdum - mat(i,j,1)/2 ! add delta function (due to positive infinitesimal imaginary part)
            selfc_hlp(i,j,ifreq) = selfc_hlp(i,j,ifreq) - cdum * weight
            if(i/=j) then
# ifndef INV
              if(freqint<=1)then; cdum = sum(wfreq*conjg(aspline)) / (2*pi*img)
              else              ; cdum = freqintegral_poles(-conjg(pole),-conjg(resid),npole1,1,rdum+img*abs(enediff),0)
              endif
# endif
              selfc_hlp(j,i,ifreq) = selfc_hlp(j,i,ifreq) - cdum * weight
            endif
          endif
        enddo
                    
      else ! FREQINT NONLIN

        if(npole1>0) then
          do ifreq = 1,size(selfc,2)
            if(oselfc==1) then ; cdum1 = freq(ifreq) * img
            else               ; cdum1 = freqr1(ifreq)
            endif
            cdum  = freqintegral_poles_nonlin(pole,resid,npole1,1,cdum1,0,enek,nfrq,frq,pfrq)
            selfc_hlp(i,j,ifreq) = selfc_hlp(i,j,ifreq) - cdum
            if(i/=j) then
              cdum  = freqintegral_poles_nonlin(-conjg(pole),-conjg(resid),npole1,1,cdum1,0,enek,nfrq,frq,pfrq)
              selfc_hlp(j,i,ifreq) = selfc_hlp(j,i,ifreq) - cdum
            endif
          enddo
        endif
                    
      endif

      end subroutine iw_integral

c ------------------------------

c
c Updates selfc with contribution from residues.
      subroutine add_residues
      implicit none
      complex_dp :: residues,residues_pade

      if(oselfc==1) Bug('oselfc==1.')
      if(nfreqc>0.and..not.allocated(matc)) Bug('matc not allocated.')
                  
      if(freqint<=2) then ! FREQINT SPLINE/PADE

        do ifreqr = 1,nfreqr
          weight = 0
          if     (iband<=bando)                     weight = -wintgr(ikpt2,iband,ispin)
          if     (iband< minb(ifreqr,ispin)) then ; weight = weight + 1d0/nkpt
          else if(iband<=maxb(ifreqr,ispin)) then ; weight = weight + wintgrc(ikpt2_,iband-minb(ifreqr,ispin)+1,ifreqr,ispin)
          endif
          if(abs(weight)>1d-10) then
            enediff = enek - freqr1(ifreqr)
            if(nfreqc==0) then
              cdum = pade_func((1d0,0d0)*abs(enediff),img*freq,cpade,nfreq)
            else
              call getcoeff(coeff,abs(enediff),real(freqc),nfreqc)
              cdum = sum(coeff(bfreqc1:bfreqc2)*matc(i,j,:))
            endif
            selfc_hlp(i,j,ifreqr) = selfc_hlp(i,j,ifreqr) + weight * cdum
            if(i/=j) then
# ifndef INV
              if(nfreqc==0) then
                cdum = pade_func(-(1d0,0d0)*abs(enediff),-img*freq,conjg(cpade),nfreq)
              else
                call getcoeff(coeff,abs(enediff),real(freqc),nfreqc)
                cdum = sum(coeff(bfreqc1:bfreqc2)*matc(j,i,:))
              endif
# endif
              selfc_hlp(j,i,ifreqr) = selfc_hlp(j,i,ifreqr) + weight * cdum
            endif
          endif
        enddo

      else ! FREQINT NONLIN

        allocate ( cvec(nfreqc) ) ; cvec = 0
        f1 = frq(1)    + enek + 1d-10
        f2 = frq(nfrq) + enek - 1d-10
        do ifreqr = 1,nfreqr
          emin = min(0d0,freqr1(ifreqr))
          emax = max(0d0,freqr1(ifreqr))
          sgn  = nint( -sign(1d0,freqr1(ifreqr)) )
          if(emax-emin>1d-8.and.f1<emax.and.f2>emin) then
            if(nfreqc>0) then
              rdum  = sgn*real(freqc(bfreqc1)) + freqr1(ifreqr)
              rdum1 = sgn*real(freqc(bfreqc2)) + freqr1(ifreqr)
              if(max(emin,f1)<min(rdum,rdum1)) then
                write(*,*) 'maxof',emin,f1
                write(*,*) min(rdum,rdum1),sgn,'lower',bfreqc1,bfreqc2
                write(*,*) 'test',sgn*real(freqc(bfreqc2)) + freqr1(ifreqr)
                read(*,*)
                Warn('Lower bound error (freqc).')
              endif
              if(min(emax,f2)>max(rdum,rdum1)) then
                write(*,*) min(emax,f2),max(rdum,rdum1),sgn,'upper'
                read(*,*)
                Warn('Upper bound error (freqc).')
              endif
              cvec(bfreqc1:bfreqc2) = matc(i,j,:)
              cdum = residues(emin,emax,enek,nfrq,frq,pfrq,nfreqc,real(freqc)+sgn*freqr1(ifreqr),cvec,sgn==-1)
            else
              cdum = residues_pade(emin,emax,enek,nfrq,frq,pfrq,npole1,sgn*pole+freqr1(ifreqr),sgn*resid)
            endif
            if(emax==0d0) cdum = -cdum
            selfc_hlp(i,j,ifreqr) = selfc_hlp(i,j,ifreqr) + cdum
            if(i/=j) then
              if(nfreqc>0) then
                cvec(bfreqc1:bfreqc2) = matc(j,i,:)
                cdum = residues(emin,emax,enek,nfrq,frq,pfrq,nfreqc,real(freqc)+sgn*freqr1(ifreqr),cvec,sgn==-1)
              else
                cdum = residues_pade(emin,emax,enek,nfrq,frq,pfrq,npole1,sgn*(-conjg(pole))+freqr1(ifreqr),sgn*(-conjg(resid)))
              endif
              if(emax==0d0) cdum = -cdum
              selfc_hlp(j,i,ifreqr) = selfc_hlp(j,i,ifreqr) + cdum
            endif
          endif
        enddo
        deallocate(cvec)

      endif
                  
      end subroutine add_residues

c ------------------------------

# if 0
c     Converts data from double precision to single precision (to regularize Pade extrapolation). Currently not used.
      subroutine regularize(val,n)
      implicit none
      integer,     intent(in)    :: n
      MCOMPLEX_dp, intent(inout) :: val(n)
      integer                    :: i
      do i = 1,n
        val(i) = ifInv( real , cmplx ) (val(i))
      enddo
      end subroutine regularize
# endif

c ------------------------------

c     Removes empty space from array a(:n,:nn) -> a(:m,:nn)  
      subroutine collapse(a,m,n,nn)
      implicit none
      MCOMPLEX_dp, intent(inout) :: a(*)
      integer,     intent(in)    :: m,n,nn
      integer                    :: i1,i2,j
      i1 = 0
      i2 = 0
      do j = 1,nn
        a(i1+1:i1+m) = a(i2+1:i2+m)
        i1           = i1 + m
        i2           = i2 + n
      enddo
      end subroutine collapse
      
c ------------------------------

c ------------------------------
# elif M_PART == 2
c ------------------------------

# ifndef Datatyp
#   define Datatyp real_dp
#   define Sufx    _r
#   include __FILE__
#   define Datatyp complex_dp
#   define Sufx    _c
# endif

      ! Collapse for cblock (collapsing diagonal blocks to contiguous array)
      subroutine collapse_screen_cblock Sufx (array)
      implicit none
      Datatyp, intent(inout) :: array(*)
      integer                :: iwblock,i,j,k,l
      k = 0
      l = 0
      do iwblock = 1,nwblock
        do j = 1,wblock(iwblock)
          do i = 1,wblock(iwblock)
            k        = k + 1
            l        = l + 1 ; if(l>nbas**2.or.k>nbas**2) Bug('counters k or l out of bounds.')
            array(l) = array(k)
          enddo
          k = k + sum(wblock) - wblock(iwblock)
        enddo
        k = k + wblock(iwblock)
      enddo
      end subroutine collapse_screen_cblock Sufx

c     ------------------------------

      ! Performs matrix-matrix product cprodout(:,:dim) = screen * cprodin(:,:dim)
      ! screen array must be cblock-collapsed, cprodin must be cblock-ordered, cprodout is cblock-ordered
      subroutine matmat_cblock Sufx (cprodout,screen,cprodin,dim)
      implicit none
      integer,     intent(in)         :: dim
      Datatyp,     intent(out)        :: cprodout(nbas,dim)
      Datatyp,     intent(in), target :: screen(*)
      MCOMPLEX_dp, intent(in)         :: cprodin(nbas,dim)
      Datatyp,     pointer_cnt        :: screenp(:,:)
      integer                         :: iwblock
      integer                         :: n,n1,nn
      nn = 0 ! index for screenp1(:)
      n1 = 0 ! index for cprodp
      do iwblock = 1,nwblock
        n                     =  wblock(iwblock)
        screenp(1:n,1:n)      => screen(nn+1:nn+n**2)
        cprodout(n1+1:n1+n,:) =  matmat(screenp,cprodin(n1+1:n1+n,:))
        n1                    =  n1 + n
        nn                    =  nn + n**2
        nullify(screenp)
      enddo
      end subroutine matmat_cblock Sufx

c     ------------------------------
      
      ! Add contribution of equivalent k points
      subroutine equivalent_kpoints Sufx (mat)
      implicit none
      Datatyp, intent(inout) :: mat(size(band),size(band))
      Datatyp                :: vec1(ndeg**2),vec2(ndeg**2)
      integer                :: i,j,l,m
      ! use contracted irreps for full matrix
      if(job1%full.and.any(nkpts(:nkpt1)/=1)) then
        call timer_start('SLF equivalent')
        mat = mat / nsym1
        do j = 1,size(band)/ndeg
          do i = 1,size(band)/ndeg
            l = (i-1)*ndeg
            m = (j-1)*ndeg
            if(trsoff) then ! no time-reversal symmetry
              vec1                       = reshape ( mat(l+1:l+ndeg,m+1:m+ndeg)        , [ ndeg**2     ] )
              mat(l+1:l+ndeg,m+1:m+ndeg) = reshape ( matmul(vec1,irrep_contr(:,:,i,j)) , [ ndeg , ndeg ] )
            else            ! with time-reversal symmetry
              if(i>j) cycle
              vec1                       = reshape ( mat(l+1:l+ndeg,m+1:m+ndeg)        , [ ndeg**2     ] )
              vec2                       = reshape ( mat(m+1:m+ndeg,l+1:l+ndeg)        , [ ndeg**2     ] )
              mat(l+1:l+ndeg,m+1:m+ndeg) = reshape ( matmul(vec1,irrep_contr(:,:,i,j)) , [ ndeg , ndeg ] )
     &                                   + reshape ( matmul(vec2,irrep_contt(:,:,i,j)) , [ ndeg , ndeg ] )
              if(i/=j)
     &        mat(m+1:m+ndeg,l+1:l+ndeg) = reshape ( matmul(vec2,irrep_contr(:,:,j,i)) , [ ndeg , ndeg ] )
     &                                   + reshape ( matmul(vec1,irrep_contt(:,:,j,i)) , [ ndeg , ndeg ] )
            endif
          enddo
        enddo
        call timer_stop('SLF equivalent',times=time_equiv)
      endif
      ! use great orthogonality theorem for diagonal elements
      if(.not.job1%full) mat(1,1) = sum( [ (mat(i,i),i=1,size(band)) ] ) / size(band)      
      end subroutine equivalent_kpoints Sufx

# undef Datatyp
# undef Sufx

c ------------------------------
# elif M_PART == 3
c ------------------------------

      end

c -----------------------------

      function scaling_selfenergy_cprod(dim,ndim)
      use global, only: ntype,neq,nindx,lcut,ncent
      use, intrinsic :: iso_fortran_env
      implicit none
      real_dp             :: scaling_selfenergy_cprod
      integer, intent(in) :: ndim,dim(ndim)
      integer             :: itype,l
      real_dp             :: overhead
      overhead = 1d0/6 * sum( [ ((neq(itype)*(2*l+1)*nindx(l,itype),l=0,lcut(itype)),itype=1,ntype) ] ) / ncent ! taken from scaling_suscep (to be optimized...)
      if(ndim/=3) Bug('Array dimension is not three.')
      scaling_selfenergy_cprod = dim(1) * ( overhead + dim(2) * dim(3) )
      end
      
c -----------------------------
# endif
c -----------------------------
