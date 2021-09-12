c--------------------------------------------------------------------------------
c Copyright (c) 2021 Peter Grünberg Institut, Forschungszentrum Jülich, Germany
c This file is part of SPEX and available as free software under the conditions
c of the MIT license as expressed in the LICENSE file in more detail.
c--------------------------------------------------------------------------------

c Calculates the susceptibility matrices (suscep) at k point ikpt (including head, wings for ikpt=1) and
c frequencies frqh(:nfrqh) (imaginary) and frqc(:nfrqc) (general):
c
c The spectral function
c
c  ss'           occ unocc      s'        s  _         _       s       s'                      s     s'                     s     s'
c S  (k,w) = SUM SUM SUM   < phi     | phi   M    >  < M    phi   | phi      >  [ delta ( w + e   - e      ) - delta ( w - e   + e      ) ]
c  IJ         q   n   n'        q+kn'     qn  k,J       k,I    qn      q+kn'                   qn    q+kn'                  qn    q+kn'
c
c is calculated and integrated to obtain the susceptibility
c
c       ss'            infinity       ss'                                                                                 ss'
c suscep   (k,w) = 2  P  INT     w'  S   (k,w) / ( w - w' +i sgn(w') eta )  d w' = SUM [ wghtSr(w ,w)+img*wghtSi(w ,w) ] S  (k,w )
c       IJ            -infinity       IJ                                            i            i                i       IJ    i
c
c where the integral is evaluated with a weighted sum over frequencies w = frqs as indicated. (s=s'=1 for spin=1 etc.)
c                                                                       i
c If trafo=.true., the matrices suscep(:,:,i) are returned in the basis {B} defined by
c                  _
c | B    > = SUM | M    > ctrafo(J,I) ,  I = 1,...,dim
c    k,I      J     k,J
c        _
c where {M} is the biorthogonal set to the mixed basis.
c
c
c Symmetry usage: <i|CHI|j> = SUM(P) <i|P^-1*CHI'*P|j>
c                           = SUM(P,no tr) <Pi|CHI'|Pj> + SUM(P,tr) <Pi*|CHI'|Pj*>
c                           = SUM(P,no tr) <Pi|CHI'|Pj> + SUM(P,tr) <Pj|CHI'|Pi>
c                           = { SUM(P,no tr) <Pi|k><l|Pj> } <k|CHI'|l>             | irrep_contr(ij,kl) = SUM(P,no tr) <Pi|k><l|Pj>
c                           + { SUM(P,tr)    <Pj|l><k|Pi> } <l|CHI'|k>             | irrep_contt(ij,kl) = SUM(P,tr)    <Pj|l><k|Pi>
c CHI' = summation only over EIBZ.
c
c
c Case multdiff_>0 multiplies <M i|j> by <M i|j>*(ej-ei). This is particularly reasonable for the wings and head at the Gamma point (but also at small k vectors)
c because there <M i|j> \propto 1/(ej-ei) varies a lot for small energy differences. The factor 1/(ej-ei) is later considered as a factor 1/w (with w being a
c Hilbert frequency) in the weighted sum of the Hilbert trafo (multdiff_=1) or in the integrations for the weighting factor (multdiff_=2).
c
c Usage of symmetry can be switched off here.
c # define switch_off_symmetry
c
c PREORDER: If defined, the arrays cprod and wtetra are reordered such that nonzero terms are (as) continuous (as possible) so as to make usage of ddot possible.
c Otherwise the old algorithm (using a pointer to nonzero elements) is used.
c PREORDER 1: Order according to energy difference
c PREORDER 2: Order according to weights
# define PREORDER 2
c
c WBLOCK: If defined, information about the block structure of suscep (and dielec, screen etc.) is stored in the array cblock. The information will later be used
c in the calculation of the self-energy (if cblock is allocated).
# define WBLOCK
c
c ACCESS_CONTIG: If defined, the arrays susceph and suscepc are defined with favorable contiguous memory access.
c However, this requires a subsequent reordering of the array, which is computationally very expensive.
c# define ACCESS_CONTIG
c
c SPECF_MSUM: If defined, the final summation for the susceptibility (step "suscep") is done in parallel using
c collective MPI reduce operations. The default is that each process updates its own contribution into the
c shared array(s) suscep(c) with load/store operations without using MPI calls.
c# define SPECF_MSUM      

# ifdef SPECF_MSUM
#   warning SPECF_MSUM defined
# endif

# include "cppmacro.h"

# ifdef switch_off_symmetry
#   warning switch_off_symmetry set
# endif

# ifdef MPI
#   define OUTP line(len_trim(line)+1:)
#   define NOA //')'
#   define MpiR_WRITE if(Mrank0==0)then;write(6,'(A'NoA)trim(line);line='';endif
# else
#   define OUTP 6
#   define NOA NoA
#   define MpiR_WRITE
# endif

# ifdef LOAD
#   define BANDDEF 1,no,no+1,no+nu
# else
#   define BANDDEF bo1,bo2,bu1,bu2
# endif

# ifdef MPI
#   define NFENCE if(nfrqc>0) then ; Nfence(suscepc) ; endif ; if(nfrqh>0) then ; Nfence(susceph) ; endif
# else
#   define NFENCE
# endif

c# define Test_cprod

      subroutine susceptibility ( ikpt,spin,dim,ctrafo,coul,trafo,symon,plasma,
     &                            nfrqc,frqc,Win(suscepc),headc,wingc , ! non-Hermitian matrices
     &                            nfrqh,frqh,Win(susceph),headh,wingh ) ! Hermitian matrices (imaginary frequency axis)

      use global
# ifdef ACCESS_CONTIG      
#   warning ACCESS_CONTIG defined
      use util, only: chr,reshape_permute_r,reshape_permute_c,array_resort_c
# else
      use util, only: chr
# endif      
      use wrapper
      use file
      use key
      use timer_util
      use, intrinsic :: iso_fortran_env
      Mpi( use Mwrapper )
      Mpi( use iso_c_binding )
      Load( use readwrite ) LoadC( only: read_wavef1 ) LoadC( read_wavef3 )

      implicit none

      integer,     intent(in)  :: dim,ikpt,spin
      integer,     intent(in)  :: nfrqc,nfrqh MpiC(win_suscepc) MpiC(win_susceph)
      complex_dp,  intent(in)  :: frqc(nfrqc),frqh(nfrqh)
      logical,     intent(in)  :: trafo,symon
      real_dp,     intent(in)  :: coul(dim)

      MCOMPLEX_dp, intent(in)  :: ctrafo(nbasm(ikpt),dim)
# ifdef ACCESS_CONTIG
      MCOMPLEX_dp, intent(out), target :: susceph(dim*(dim+1)/2,nfrqh)
      complex_dp,  intent(out), target :: suscepc(dim,dim,nfrqc)
      MCOMPLEX_dp, intent(out) :: headh(6,nfrqh),wingh(3,dim,nfrqh)
      complex_dp,  intent(out) :: headc(3,3,nfrqc),wingc(3,dim,2,nfrqc)
      complex_dp,  pointer_cnt :: suscepc_(:,:)
      MCOMPLEX_dp, pointer_cnt :: susceph_(:,:)
      external                 :: get_index_suscepc  
# else
      MCOMPLEX_dp, intent(out) :: susceph(dim*(dim+1)/2,nfrqh),headh(6,nfrqh),wingh(3,dim,nfrqh)
      complex_dp,  intent(out) :: suscepc(dim,dim,nfrqc),headc(3,3,nfrqc),wingc(3,dim,2,nfrqc)
# endif
      real_dp,     intent(out) :: plasma

      MCOMPLEX_dp, pointer_cnt :: irrep_contr(:),irrep_contt(:)
      MCOMPLEX_dp, allocatable :: cprod(:,:),dcprod(:,:),momentum(:,:,:,:),ctrafo1(:,:)
# ifdef FIBC
      complex_dp,  allocatable :: ibc_suscep1(:,:,:)
# else
      complex_dp,  allocatable :: ibc_suscep1(:,:,:,:)
# endif
      complex_dp,  allocatable :: irrep(:,:,:,:)
      complex_dp,  allocatable :: wing0(:,:,:),wing1(:,:,:),wing2(:,:,:)
      complex_dp,  allocatable :: suscep0(:,:),suscep1(:,:),suscep0_(:,:),suscep1_(:,:)
      complex_dp               :: head0(3,3),head2(3,3)
      complex_dp               :: carr(nbasm(ikpt)),cdum NoInvC(cdum1)
      MCOMPLEX_dp              :: mdum1,mdum2

      real_dp,     allocatable :: wghtSr(:,:),wghtSi(:,:),wtetra(:,:),re(:),im(:),sf(:)
      real_dp,     allocatable :: sus1(:),susc1(:) NoInvC(sus2(:)) NoInvC(susc2(:))
      real_dp                  :: rot(3,3),wdrude(bando,nkpt),integrand(bando,nkpt)
      real_dp                  :: rdum,err,memo,memc,memp,mnc,mxc
      real                     :: time1,time2

      integer,     allocatable :: nindex(:),index(:,:),pnt(:)
      integer                  :: kpt1(nkpt),nkpt1,kpt1p(nkpt),nkpts(nkpt),sym1(nsym),nsym1,kpt1p_(nkpt)
      integer                  :: ifrqS,ifrq,ikpt1,ikpt2,isym,iunit
      integer                  :: ispin1,ispin2
      integer                  :: numax,nu,nomax,no,nkmax,bu1,bu2,bo1,bo2,no1,no2,nu1,nu2,nk,nk1,nk2,bk1,bk2
      integer                  :: nirrep,dirrep(dim),maxdirrep,d,d1,d2,d3
      integer                  :: i,i1,i2,j,k,k1,k2,l,ll,n,m,ij,iterm
      integer                  :: ii,ii1,ii2
      integer                  :: multdiff_
      real,        pointer_cnt :: sirrep(:,:)
      real,        parameter   :: sirrep_threshold = 1e-9
      real                     :: sdum
      logical                  :: first_term
      logical                  :: ldum
      character(80)            :: line
      character(22)            :: wrt
      integer                  :: nfrqs
# ifndef SPECF_MSUM
      integer                  :: nfrqs1,frqs1
# endif
      real_dp                  :: fspec(2)
      real_dp,     allocatable :: frqs(:)
      real_dp                  :: ddot    
      integer                  :: kptsum

# ifdef MPI
      integer                  :: ggt,rank,msiz,bounds(4,3)
      integer                  :: Merr,win_irrep_contr,win_irrep_contt,win_sirrep,nskip
      integer                  :: mm,m1,m2
      type(c_ptr)              :: ptr
      real_dp, external        :: scaling_suscep
# endif

c     Parameters generated from spex.inp
      logical, save            :: firstcall=.true.,lhub
      integer, save            :: multdiff
      real_dp, save            :: disorder,wghtthr
      real_dp, save, allocatable :: hub(:,:,:)

# include "interface/getfrqs.inc"

      call timer_start('Routine susceptibility SUS')
      Rcall timer_start('SUS preparations')

      CHKMEM0(entered sus)

      if(spin/=0.and.l_soc) Error('spin/=0 and SOC not implemented.')

      NFENCE

      Oif(nfrqh>0) susceph = 0
      Oif(nfrqc>0) suscepc = 0
      if(ikpt==1) then
        plasma = 0
        if(nfrqh>0) headh = 0
        if(nfrqh>0) wingh = 0
        if(nfrqc>0) headc = 0
        if(nfrqc>0) wingc = 0
      endif

      NFENCE

# ifdef ACCESS_CONTIG
      if(nfrqh>0) susceph_(1:nfrqh,1:dim*(dim+1)/2) => susceph ! Help pointers for faster memory access
      if(nfrqc>0) suscepc_(1:nfrqc,1:dim**2)        => suscepc ! (Matrices will be resorted by reshape_permute)
# endif

c
c     IBC : Calculate contraction term (->ibc_suscep0)
      if(oibc/=0) then
        if(nfrqc>0) Error('IBC not implemented for real frequencies.')
        if(ikpt==1) call ibc1_contraction(frqh,nfrqh,.true.)
# ifdef FIBC
        allocate ( ibc_suscep1(dim,nbasp,nfrqh) )
# else
        allocate ( ibc_suscep1(maxlmindxm,maxlmindxm,ncent,nfrqh) )
# endif
        ibc_suscep1 = 0
      endif

c
c     Determine Hilbert mesh
      Rcall getfspec(fspec,nfrqh,nfrqc,frqc,.true.)
      Mpi( call Mcast(fspec) )

c
c     Input from spex.inp
      if(firstcall) then
        Rbegin
        call getkey(inp,'WGHTTHR', wghtthr, section='SUSCEP', default=1d-10, mini=0d0)
        call getkey(inp,'MULTDIFF',   line, section='SUSCEP', default='   ')
        if     (line==' '  )   then ; multdiff = 0 ! Default: MULTDIFF ON for Gamma point, else OFF
        else if(line=='ON ')   then ; multdiff = 1 ! always ON
        else if(line=='OFF')   then ; multdiff = 2 ! always OFF
        else if(line=='INT')   then ; multdiff = 3 ! Default   with integrated 1/w
        else if(line=='INTON') then ; multdiff = 4 ! always on with integrated 1/w
        else                        ; Error('Choose MULTDIFF ON, OFF, INT or INTON.')
        endif
        call getkey(inp,'DISORDER', disorder, section='SUSCEP',  default=0d0, mine=0d0 )
        if(disorder/=0) disorder = 1/(2*disorder) ! eta = 1/(2*tau)
        call getkey(inp,'HUBBARD', line, section='SUSCEP', status=j)
        if(j==0) then
          lhub = .false.
        else
          lhub = .true.
          if(nwan==0) Error('No Wannier orbitals defined. (Needed for HUBBARD.)')
          if(l_soc)   Error('SOC not implemented for HUBBARD.')
          allocate ( hub(maxband,nkpt,nspin) ) ; hub = 0
          ldum = .false.
          if(j==2) then
            inquire(file=line,exist=ldum)
          endif
          if(ldum) then
            iunit = fopen(line,form='formatted',status='old')
            read(iunit,*) i1,i2
            do ispin1 = 1,nspin ; do k = 1,nkpt ; do i = i1,i2
              read(iunit,*) hub(i,k,ispin1)
            enddo ; enddo ; enddo
            call fclose(iunit)
            write(6,'(A)') 'Hubbard subspace read from file '//trim(line)
          else
            do ispin1 = 1,nspin ; do k = 1,nkpt ; do i = wanbandi,wanbandf
              hub(i,k,ispin1) = sum(abs(uwan(i,:,k,ispin1))**2)
            enddo ; enddo ; enddo
            if(j==2) then
              iunit = fopen(line,form='formatted',status='new')
              write(iunit,*) wanbandi,wanbandf
              do ispin1 = 1,nspin ; do k = 1,nkpt ; do i = wanbandi,wanbandf
                write(iunit,*) hub(i,k,ispin1)
              enddo ; enddo ; enddo
              call fclose(iunit)
              write(6,'(A)') 'File '//trim(line)//' written.'
            else
              ldum = .true.
            endif
          endif
        endif
        Rend
        Mpi( call Mbroadcast1 )
        if(lhub.and..not.ldum) call finish
      endif

      if(any(real(frqh)/=0)) Error('Found nonzero real part in frequencies (frqh).')

      Rbegin

      write(6,'(/A)') '---------'
      write(6,'(/A)')            'Calculation of susceptibility'
      write(6,'(/A,I7)')         'Order of array suscep:',dim
      write(6,'(A,F10.1," MB")') 'Size of array suscep:',((MBYTES*nfrqh+16d0*nfrqc)*dim**2)/megabyte
      if(dim<10) Warn('Array dimension lower than 10. Is this correct?')

c
c     Get irreducible k points kpt1(:nkpt1) (wrt current k point ikpt)
      if(symon) then ; call getkpt1       (kpt1,nkpt1,nkpts,sym1,nsym1,ikpt,0,.true.)
      else           ; call getkpt1_fullBZ(kpt1,nkpt1,nkpts,sym1,nsym1,0)
      endif
      call getkpt1p(kpt1p,kpt1,nkpt1,sym1,nsym1,nkpts,.false.)

      Rend

      first_term = .true.
      Mpi(Load( rank=Mrank;msiz=Msize ))
      Mpi(Load( Rif(firstcall) Info('Parallelization over k disabled for LOAD.') ))
      Mpi( call Mbroadcast2 )

c
c     Prepare irreducible representations inside eigenspaces of coulomb (->cblock,irrep,dirrep,nirrep)
      Nallocate(sirrep,(S_ dim,dim S_))
      ifO sirrep = 1. + sirrep_threshold ! sirrep > sirrep_threshold by default
      Nfence(sirrep)
      if(nsym1>1.and.trafo) then
        Rcall timer_stop('SUS preparations')
        call timer_start('SUS irreps')
        ! (a) Look for "degenerate" eigenvalues of coulomb
c        if(ikpt==1) then ; call deg_limit(mnc,mxc,coul(2:),dim-1,.true.)
c        else             ; call deg_limit(mnc,mxc,coul,    dim,  .true.)
c        endif
c        write(*,*) mnc,mxc
        mnc    =  huge(0d0)
        mxc    = -huge(0d0)
        err    = 0
        n      = 1
        nirrep = 0
        do while(n<=dim)
          nu = n
          do
            nu = nu + 1
            if(nu>dim)                     exit
            if(abs(coul(nu)-coul(n))>cdeg) exit
          enddo
          if(n>1)    mnc = min(mnc,abs(coul(n)-coul(n-1)))
          if(nu-1>n) mxc = max(mxc,abs(coul(nu-1)-coul(n)))
          nirrep         = nirrep + 1
          dirrep(nirrep) = nu - n
          n              = nu
        enddo
        maxdirrep = maxval(dirrep(:nirrep))
        if(mxc>=0) then ; Rwrite(6,'(A,2ES8.1)') 'Min/max eigenvalue diffs:  ',mnc,mxc
        else            ; Rwrite(6,'(A,2ES8.1)') 'Min eigenvalue diff:       ',mnc
        endif
        Allocate_ ( irrep,(maxdirrep,maxdirrep,nirrep,nsym1) ) ; irrep = 0
        ! (b) Calculate irreps
        allocate ( ctrafo1(nbasm(ikpt),1) )
        irrep = 0
        n     = 0
        CHKMEM0(before loop)
        do k = 1,nirrep
          d = dirrep(k)
          ifMOD(k)
          do i1 = 1,d
            ctrafo1(:nbasp,1) = ctrafo(:nbasp,n+i1)
            call olapvec(ctrafo1(nbasp+1:,1),ctrafo(nbasp+1:,n+i1),gptm(:,pgptm(:ngptm(ikpt),ikpt)),ngptm(ikpt))
            do isym = 1,nsym1
              call mtrafo(carr,(1d0,0d0)*ctrafo1(:,1),nbasm(ikpt),1,ikpt,sym1(isym),1,.false.)
              do i2 = 1,d
                irrep(i1,i2,k,isym) = dotprod(ctrafo(:,n+i2),carr)
              enddo
            enddo
          enddo
          do isym = 1,nsym1
            rdum = maxval ( abs ( identity(d) - matmat(conjg(transpose(irrep(:d,:d,k,isym))),irrep(:d,:d,k,isym)) ) )
            err  = max(rdum,err)
          enddo
          endMOD
          n = n + d
        enddo
        CHKMEM0(after loop)
        Mpi( call Msum(irrep) )
        deallocate ( ctrafo1 )
c        Rcall cpu_done(time1)
        ! (c) Construct contracted irreps
        n = sum ( [ (((dirrep(k1)*dirrep(k2))**2,k1=1,k2),k2=1,nirrep) ] )
        CHKMEM0(before irrep_contr alloc)
        Nallocate( irrep_contr, (S_ n S_) )
        ifO irrep_contr = 0 ; Nfence(irrep_contr)
        if(.not.trsoff) then
          Nallocate( irrep_contt, (S_ n S_) )
          ifO irrep_contt = 0 ; Nfence(irrep_contt)
        endif
        CHKMEM0(done irrep_contr alloc)
        l = 0
        do k2 = 1,nirrep ; d2 = dirrep(k2) !
          do k1 = 1,k2   ; d1 = dirrep(k1) ! Actually irrep_contr(t)=0 if d1/=d2, but near degeneracies can lead to problems
            n = (d1*d2)**2
            ifMOD(l)
            if(trsoff) then
              irrep_contr(l+1:l+n) = irrep_contract(0)
            else
              irrep_contr(l+1:l+n) = irrep_contract(1)
              irrep_contt(l+1:l+n) = irrep_contract(2)
            endif
            endMOD
            l = l + n
          enddo
        enddo
        CHKMEM0(after irrep contr/t 0)
        Nfence(irrep_contr) ; Mpi( Ocall Msum(irrep_contr,comm=Ocomm) ) ; Nfence(irrep_contr)
        if(.not.trsoff) then
          Nfence(irrep_contt) ; Mpi( Ocall Msum(irrep_contt,comm=Ocomm) ) ; Nfence(irrep_contt)
        endif
        CHKMEM0(after irrep contr/t 1)
c        Rcall cpu_done(time1)
        ! (d) Determine non-zero elements (->sirrep, for susceptibility)
        Nfence(sirrep)
        l  = 0
        i2 = 0
        do k2 = 1,nirrep ; d2 = dirrep(k2) ; i1 = 0
          do k1 = 1,k2   ; d1 = dirrep(k1)
            ifNMOD(k2*(k2-1)/2+k1)
            sdum = sum(abs(irrep_contr(l+1:l+(d1*d2)**2))) / (d1*d2)**2
            if(.not.trsoff) sdum = sdum + sum(abs(irrep_contt(l+1:l+(d1*d2)**2))) / (d1*d2)**2
            sirrep(i1+1:i1+d1,i2+1:i2+d2) = sdum
            sirrep(i2+1:i2+d2,i1+1:i1+d1) = sdum
            endMOD          
            l  = l + (d1*d2)**2
            i1 = i1 + d1
          enddo
          i2 = i2 + d2
        enddo
        Nfence(sirrep)
        ! (e) Find irrep indices (nirrep)        
        allocate ( nindex(nirrep) ) ! nindex = index for irrep
        nindex = 0
        i      = 0
        i2     = 1
        do k2 = 1,nirrep ; i1 = 1
          if(nindex(k2)==0) then
            i          = i + 1
            nindex(k2) = i
          endif
          do k1 = 1,nirrep
            sdum = sirrep(i1,i2)
            if(sdum>=sirrep_threshold) then
              if(nindex(k1)==0) then
                nindex(k1) = nindex(k2)
              else if(nindex(k1)/=nindex(k2)) then
                l = 0
                m = 0
                n = nindex(k1)                
                do k = 1,nirrep
                  if     (nindex(k)==nindex(k2)) then ; l = l + 1
                  else if(nindex(k)==n)          then ; m = m + 1 ; nindex(k) = nindex(k2)
                  endif
                enddo
                RInfo('Must merge blocks: '//Chf(sdum,'ES8.1')//', sizes: '//Chr(l)//'+'//Chr(m))
              endif
            endif
            i1 = i1 + dirrep(k1)
          enddo
          i2 = i2 + dirrep(k2)
        enddo
        n = i
        j = 0
        do i = 1,n
          if(count(nindex==i)>0) then
            j = j + 1
            where(nindex==i) nindex = j
          endif
        enddo
        n = j
        do i = 1,n
          if(count(nindex==i)==0) Error('Block definition failed.')
        enddo
        CHKMEM0(after block def)
        Rwrite(6,'(A,I12'NoA) 'Number of blocks:',n
# ifdef WBLOCK
        ! (f) Define blocks (->cblock, for selfenergy)
        if(n>1) then
          allocate ( cblock(dim) )
          j = 0
          do k = 1,nirrep
            do i = 1,dirrep(k)
              j         = j + 1
              cblock(j) = nindex(k)
            enddo
          enddo
          Rbegin
          write(6,'(A'NoA) ' : '
          do i = 1,n
            j = count(cblock==i) ; if(i/=1) write(6,'(A'NoA) '/'
            write(6,'(I'//chr(int(log10(j+.1d0))+1)NoA) j
          enddo
          Rend
        endif
        Rwrite(6,*)
# elif !defined(noWARN)
#   warning WBLOCK disabled
# endif
        deallocate ( nindex )
        call timer_stop('SUS irreps',time=time1)
        Rbegin
        write(6,'(A,F13.1,A)') 'Non-zero elements:',count(sirrep>=sirrep_threshold)*100d0/dim**2,' %   (Timing:'//Chr(time1)//')'
        if(err>1d-2) Error('Too large deviation of irrep from unitarity: '//Chr(err)//'. Try to increase LEXP or use CUTZERO.')
        if(err>1d-8) Warn('Large deviation of irrep from unitarity: '//Chf(err,'F11.8'))
        Rend
        Rcall timer_start('SUS preparations')
      endif

      if     (spin==0) then ; ispin1 = 1 ; ispin2 = 1
      else if(spin==1) then ; ispin1 = 1 ; ispin2 = 1
      else if(spin==2) then ; ispin1 = 2 ; ispin2 = 2
      else if(spin==3) then ; ispin1 = 2 ; ispin2 = 1
      else if(spin==4) then ; ispin1 = 1 ; ispin2 = 2
      else                  ; Error('unknown spin index.')
      endif

      Rcall timer_stop('SUS preparations')

c
c     Spin loop
      iterm = 0
 2    iterm = iterm + 1

      CHKMEM0(bef Hilbert)

c
c     Define Hilbert mesh
      Rcall timer_start('SUS preparations')
      nfrqs = 0
      if(fspec(1)/=0) then        
        allocate ( re(2) )
        Rbegin
        ! Get frequency range -> re(1)..re(2)        
        if(gauss(1)/=0.or.fspec(1)<0) then
          re(1) = 0 ; re(2) = maxval(ene(:maxeband,:nkpt,:)) - minval(ene(:maxeband,:nkpt,:))
        else
          call tetrahedron5_init(re,[0d0],0,ikpt,kpt1,kpt1p,nkpt1,1,ispin1,ispin2,1,bando,bandu,maxeband,1,2,1,1)
        endif
        if(re(1)>=re(2)) then
          if(all(ene(1,kpt1(:nkpt1),ispin1)>=efermi)) then
            Warn('No virtual transitions because all states unoccupied.')
          else if(all( [ (ene(maxeband,kptsum(kpt1(ikpt1),ikpt),ispin2),ikpt1=1,nkpt1) ] <= efermi)) then
            Warn('No virtual transitions because all states occupied.')
          else
            Bug('Tetrahedron initialization returned empty range.')
          endif
        else
          ! Define
          call getfrqs(frqs,nfrqs,fspec,re(1),re(2),.true.)
          if(first_term) then
            write(6,'(A'NoA) 'Number of frequencies: '
            if(nfrqh> 0) write(6,'(I6'NoA) nfrqh
            if(nfrqc> 0) write(6,'(I6'NoA) nfrqc
            if(nfrqs/=0) write(6,'(I6'NoA) nfrqs
            write(6,*)
            ! Check
            if(nfrqc>0) then
              if(all(imag(frqc)==0)) then
                i = 0
                do ifrqs = 1,nfrqs-1
                  i = max(i,count(real(frqc)>frqs(ifrqs).and.real(frqc)<frqs(ifrqs+1)))
                enddo
                if(i>10) Warn('Coarse Hilbert mesh. Line segments in imaginary part of length: '//chr(i))
              endif
            endif
          endif
        endif
        Rend
        Mpi( call Mcast(re) )
        if(re(1)>=re(2)) then ; deallocate ( re ) ; goto 20 ; endif ! skip
        deallocate( re )
        Mpi( call Mcast(nfrqs) ; call Mcastl(frqs) )
      endif

      multdiff_ = 0 ! off
      if(spin<=2.and.nfrqs/=0) then
        if(multdiff==0.and.ikpt==1 .or. multdiff==1) multdiff_ = 1 ! old behavior : multiply 1/w in weighted sum
        if(multdiff==3.and.ikpt==1 .or. multdiff==4) multdiff_ = 2 ! new behavior : 1/w included in weighting factors
      endif
      Rif(all(multdiff/=[0,2])) then
        if(.not.trafo.and.ikpt/=1) Info('MULTDIFF has no effect for k/=0 with NOAPW.')
        if(nfrqs==0)               Info('MULTDIFF has no effect with HILBERT NONE.')
        if(spin>2)                 Info('MULTDIFF has no effect for magnetic susceptibility (spins different).')
      endif

c
c     Allocate arrays sus1, susc1, sus2, and susc2
      Allocate_( sus1,(nfrqh+nfrqc) ) ; NoInv( Allocate_ ( sus2,(nfrqh+nfrqc) ) )      
      Allocate_( susc1,(nfrqc)      ) ; NoInv( Allocate_ ( susc2,(nfrqc)      ) )      
      
c
c     Determine weights for Hilbert transformation (Hilbert mesh -> requested mesh)
      if(nfrqs/=0) then
        allocate ( wghtSr(nfrqs,nfrqh+nfrqc) )
        allocate ( wghtSi(nfrqs,nfrqh+nfrqc) )
        if(multdiff_/=2) call get_wghtS(0)
      endif

      Rcall timer_stop('SUS preparations')

c
c     Calculate packet sizes
c
c         MPI package of current rank: nk1->nk2|nk, no1->no2|no, nu1->nu2|nu     (lower->upper|number)
c     max MEM package of current rank:       nkmax,       nomax,       numax
c         MEM package of current loop: bk1->bk2|nk, bo1->bo2|no, bu1->bu2|nu     (nk, no, nu redefined, kpt1 and nkpts adjusted to kpoint packet)

      ! MPI work package (no1:no2,nu1:nu2,nk1:nk2)
      if(first_term) then
        nk1 = 1
        nk2 = nkpt1
        no1 = 1
        no2 = bando
        nu1 = 1
        if(nfrqs==0) then ; nu2 = maxband
        else              ; nu2 = maxeband
        endif
        if(oibc==0) nu2 = nu2 - bandu + 1
# ifdef MPI
#   ifdef LOAD
        call Mblocks(bounds,[bando,nu2-nu1+1],2,scaling_suscep)
        bounds(:,2:) = bounds(:,:2)
        bounds(:,1)  = [1,nkpt1,0,1]
#   else
        call Mblocks(bounds,[nkpt1,bando,nu2-nu1+1],3,scaling_suscep)          
#   endif
        nk1 = bounds(1,1) ; no1 = bounds(1,2) ; nu1 = bounds(1,3)
        nk2 = bounds(2,1) ; no2 = bounds(2,2) ; nu2 = bounds(2,3)
# endif
        if(oibc==0) then        
          nu1 = nu1 + bandu - 1
          nu2 = nu2 + bandu - 1
        endif
      endif
      no = no2 - no1 + 1
      nu = nu2 - nu1 + 1
      nk = nk2 - nk1 + 1

# ifdef MPI
      Rif(first_term) then
        write(6,'(A,10X,A)') 'MPI work packages:',Chr(bounds(4,1))//'*'//Chr(bounds(4,2))//'*'//Chr(bounds(4,3))//
     &                                      ' ('//Chr(nk)//'*'//Chr(no)//'*'//Chr(nu)//')'
      endif
# endif

      ! MEM work packages
      Nfence(mem)
      d     = 0 ; Mpi( if(no*nu*nk==0) goto 4 ) ! idle ranks skip
      nomax = no ! full packet if MEM sufficient
      numax = nu
      nkmax = nk
# ifdef PREORDER
      memp = (MBYTES*maxbasm+8*nfrqs        ) ! memory demand for wavefproducts, wtetra, and index (unless PREORDER) is memp*No*Nu*Nk
# else
      memp = (MBYTES*maxbasm+8*nfrqs+4*nfrqs)
# endif
# ifdef LOAD
      if(ikpt==1) then ; memc =    (16d0*maxlmindx+MBYTES*maxgpt)*nspin3 ! LOAD: memory demand for wavefunctions is memc*(No+Nu)*Nk
      else             ; memc = max(16d0*maxlmindx,MBYTES*maxgpt)*nspin3
      endif
      memo = ( maxmem - mem ) Mpi(/Nsize0) ! memory available for current rank
      if( memo < (memc*(no+nu)+memp*no*nu)*nk ) then ! we need MEM work packages
        rdum  = 2*memc/(3*memp)
        nomax = min ( max ( int ( sqrt(rdum**2+memo/(3*memp)) - rdum ), 1 ) , no ) ! from  no*nu*memp + (no+nu)*memc = memo  with nu = 3*no
        d1    = (no-1)/nomax + 1
        nomax = (no+d1-1)/d1
        numax = (memo-nomax*memc) / (memc+nomax*memp)
      endif
# else
      memo = ( maxmem - mem ) Mpi(/Nsize0) ! memory available for current rank
      if( memo < memp*no*nu*nk ) then ! we need MEM work packages
        if( memo >= memp*no*nu ) then ! kpoint division sufficient
          nkmax = int( memo / (memp*no*nu) )
        else                          ! ... insufficient
          nkmax = 1
          if( memo >= memp*nu ) then ! occ-states devision sufficient
            nomax = int( memo / (memp*nu) )
          else                       ! ... insufficient
            if( memo < memp )
     &        Error('Maximal allowed storage exceeded. Increase MEM to at least '//chr(mem+memp)//' or reduce parameters.')
            nomax = 1
            numax = int( memo / memp )
          endif
        endif
      endif
      d1    = (no-1)/nomax + 1
      nomax = (no+d1-1)/d1      
# endif
      d2    = (nu-1)/numax + 1
      numax = (nu+d2-1)/d2
      d3    = (nk-1)/nkmax + 1
      nkmax = (nk+d3-1)/d3
      Rif(first_term) then
        if(nomax/=no.or.numax/=nu.or.nkmax/=nk) then
          write(6,'(A,10X,A)') 'MEM work packages:',Chr(d3)//'*'//Chr(d1)//'*'//Chr(d2)//
     &                                     ' ('//Chr(nkmax)//'*'//Chr(nomax)//'*'//Chr(numax)//')'
        endif
        write(6,'(A,F13.1,A)') 'Auxiliary storage:',memo/megabyte,' MB'
      endif
      if(nomax*numax*nkmax==0) Error('Maximal allowed storage exceeded. Increase MEM or reduce number of frequencies!')

      ! Number of MEM work packages (d=d1*d2*d3) can be different in different ranks -> some ranks must skip (->nskip)
# ifdef MPI
      d = d1 * d2 * d3
 4    Nfence(mem)
      call mpi_allreduce(d,d1,1,mpi_integer,mpi_max,Mcomm,Merr)
      nskip = d1 - d
      if(d==0) goto 5
# endif

      if(trafo.and.wrtinfo.and.nfrqs/=0) allocate ( sf(nfrqs) )

      Rbegin
      if(nspin/=1) then
        if(ispin1==1.and.ispin2==1) write(6,'(/A'NoA) 'up/up    '
        if(ispin1==2.and.ispin2==2) write(6,'(/A'NoA) 'down/down'
        if(ispin1==1.and.ispin2==2) write(6,'(/A'NoA) 'up/down  '
        if(ispin1==2.and.ispin2==1) write(6,'(/A'NoA) 'down/up  '
      else
        write(6,'(/A'NoA) '         '
      endif
      write(6,'(A)')            '        Bands | Timings'
      write(6,'(A'NoA) '                       |      MT      PW'
# ifdef FIBC
      if(trafo)   write(6,'(A'NoA) '   trafo'
      if(oibc/=0) write(6,'(A'NoA) '     IBC'
# else
      if(oibc/=0) write(6,'(A'NoA) '     IBC'
      if(trafo)   write(6,'(A'NoA) '   trafo'
# endif
      if(ikpt==1.and.spin<=2) write(6,'(A'NoA) '      kp'
      if(nfrqs/=0) write(6,'(A'NoA) '   wghts'
      write(6,'(A)') '  suscep'
      Rend
      if(trafo.and.wrtinfo.and.nfrqs/=0) sf = 0

c
c     Loop over states in packets of nk kpoints and no and nu bands
      line = ' '
      bo1  = no1
      bu1  = nu1
      bk1  = nk1
      do while(bk1<=nk2) ; bk2 = min(bk1+nkmax-1,nk2) ; nk = bk2 - bk1 + 1 ; bo1 = no1
      do while(bo1<=no2) ; bo2 = min(bo1+nomax-1,no2) ; no = bo2 - bo1 + 1 ; bu1 = nu1
      do while(bu1<=nu2) ; bu2 = min(bu1+numax-1,nu2) ; nu = bu2 - bu1 + 1
        wrt = Chr(bk1)//'/'//Chr(bo1)//'/'//Chr(bu1)//' '//Chr(nk)//'*'//Chr(no)//'*'//Chr(nu)
        write(OUTP,'(A'NoA) adjustr(wrt)

        MpiR_WRITE

c                           s       s'
c       Calculate < M    phi   | phi     > (Here, M is biorthogonal if .not.fullpw, n=occ, n'=unocc, s=ispin1, s'=ispin2)
c                    k,I    qn      q+kn'

        if(dble(nu)*no*nk>huge(0)) Bug('Integer overflow for cprod. Please report.')
        Allocate_ ( cprod,(nu*no*nk,nbasm(ikpt)) )

        Load( Rcall timer_start('SUS read_wavef1') )
        Load( if(ikpt==1) call read_wavef1(bo1,bo2,bu1,bu2,kpt1(bk1),nk,ikpt,ispin1,ispin2,.true.,.true.)  ) ! k=0  : we need cmt,cpw for later
        Load( if(ikpt/=1) call read_wavef1(bo1,bo2,bu1,bu2,kpt1(bk1),nk,ikpt,ispin1,ispin2,.true.,.false.) ) ! k/=0 : only read cmt
        Load( Rcall timer_stop('SUS read_wavef1') )

        beginSplit(bo1*nkpt+kpt1(bk1))
        call timer_start('SUS wavefproducts_mt')
        call wavefproducts2_mt(cprod,            kpt1(bk1),nk,ikpt,ispin1,ispin2, BANDDEF )
        call timer_stop('SUS wavefproducts_mt',time=time1)
        write(OUTP,'(F10.2'NOA) time1 ; MpiR_WRITE
        endSplit

        Load( Rcall timer_start('SUS read_wavef1') )
        Load( if(ikpt/=1) deallocate(cmt) )
        Load( if(ikpt/=1) call read_wavef1(bo1,bo2,bu1,bu2,kpt1(bk1),nk,ikpt,ispin1,ispin2,.false.,.true.) ) ! k/=0 : only read cpw
        Load( Rcall timer_stop('SUS read_wavef1') )
        
        beginSplit(bo1*nkpt+kpt1(1))
        call timer_start('SUS wavefproducts_pw')
        call wavefproducts2_pw(cprod(:,nbasp+1:),kpt1(bk1),nk,ikpt,ispin1,ispin2, BANDDEF )
        call timer_stop('SUS wavefproducts_pw',time=time1)
        write(OUTP,'(F8.2'NOA) time1 ; MpiR_WRITE
        endSplit

        Load( if(ikpt/=1) deallocate(cpw) )

# ifndef FIBC
        if(oibc/=0) then
          call double_counting ; if(.not.allocated(cprod)) cycle
        endif
# endif

c
c       Transform to eigenvectors of Coulomb matrix
# define STEP (dim/4)
        if(trafo) then
          call timer_start('SUS trafo')
          NoInv( cprod = conjg(cprod) )
          if(fullpw) then ! in this case wavefproducts returns <M  phi|phi>
            do i = 1,nu*no*nk,STEP
              j               = min(i+STEP-1,nu*no*nk)
              cprod(i:j,:dim) = matmat ( cprod(i:j,:),ctrafo  )
            enddo
          else            ! in this case wavefproducts returns <M~ phi|phi>
            allocate ( ctrafo1(nbasm(ikpt),dim) )
            ctrafo1(:nbasp,:) = ctrafo(:nbasp,:dim)
            do i = 1,dim
              call olapvec(ctrafo1(nbasp+1:,i),ctrafo(nbasp+1:,i),gptm(:,pgptm(:ngptm(ikpt),ikpt)),ngptm(ikpt))
            enddo
            do i = 1,nu*no*nk,STEP
              j               = min(i+STEP-1,nu*no*nk)
              cprod(i:j,:dim) = matmul ( cprod(i:j,:),ctrafo1 )
            enddo
            deallocate ( ctrafo1 )
          endif
          NoInv( cprod = conjg(cprod) )
          call timer_stop('SUS trafo',time=time1)
          write(OUTP,'(F8.2'NOA) time1 ; MpiR_WRITE
# ifdef Test_cprod
          call test_cprod(142,cprod,sqrt(sum(matmul(rlat,kpt(:,ikpt))**2)))
# endif
        endif
        call checkmem('trafo done',0d0)

# ifdef FIBC
        if(oibc/=0) then
          call double_counting ; if(.not.allocated(cprod)) cycle
        endif
# endif

c
c       Calculate derivative of < exp(ikr) phi_n(q) | phi_n'(q+k) > for head and wings
c       which is given by < phi_n | -i nabla | phi_n' > / ( e_n' - e_n ) .
c       If multdiff_>0, we calculate only the numerator and later divide the weights
c       by the denominator. This is necessary to converge the spectral function for
c       transition metals (e.g., Ni) at a small (or zero) k vector.
        if(ikpt==1.and.ispin1==ispin2) then
          allocate ( dcprod (nu*no*nk,3) )
          call timer_start('SUS momentum')
          call momentum_matrix( dcprod,kpt1(bk1),nk,ispin1,ispin2 , BANDDEF , .false. MpiC(.false.) LoadC(.false.) )
          call timer_stop('SUS momentum',time=time1)
          write(OUTP,'(F8.2'NOA) time1 ; MpiR_WRITE
# ifdef Test_cprod
          allocate(ctrafo1(nu*no*nk,1))
          rot(:,1) = matmul(rlat,kptadd) ; rdum = sqrt(sum(rot(:,1)**2))
          rot(:,1) = rot(:,1) /rdum
# endif
          dcprod = dcprod / sqrt(vol)
          if(multdiff_==0) then
            i = 0
            do k = bk1,bk2
              ikpt1 = kpt1(k)
              do i1 = bo1,bo2
                do i2 = bu1,bu2
                  i    = i + 1
                  rdum = ene(i2,ikpt1,ispin2) - ene(i1,ikpt1,ispin1)
                  if(abs(rdum)>1d-8) then ; dcprod(i,:) = dcprod(i,:) / rdum
                  else                    ; dcprod(i,:) = 0
                  endif
# ifdef Test_cprod
                  ctrafo1(i,1) = dot_product(dcprod(i,:),rot(:,1))
# endif
                enddo
              enddo
            enddo
          endif
          Load( deallocate(cmt,cpw) )
# ifdef Test_cprod
          call test_cprod(140,ctrafo1,1d0)
          deallocate(ctrafo1)
# endif
        endif

        call checkmem('momentum done',0d0)

c       Leave out certain transitions (keywords HUBBARD and OMIT)
        if(lhub.or.lomit(2)) then
          i = 0
          do k = bk1,bk2
            ikpt1 = kpt1(k)
            ikpt2 = kptsum(ikpt,ikpt1)
            do i1 = bo1,bo2
              do i2 = bu1,bu2
                i    = i + 1
                rdum = 1
                if(lhub) then
                  rdum = 1 - hub(i1,kptp(ikpt1),ispin1) * hub(i2,kptp(ikpt2),ispin2)
                endif
                if(lomit(2)) then
                  if(any(omit==i1).or.any(omit==i2)) rdum = 0
                endif
                if(rdum/=1) then
                  cprod(i,:) = cprod(i,:) * rdum
                  if(ikpt==1.and.ispin1==ispin2) dcprod(i,:) = dcprod(i,:) * rdum
                endif
              enddo
            enddo
          enddo
        endif

c
c       multdiff_: Multiply cprod(:,1) by enediff (only relevant for case ikpt/=1)
        if(multdiff_/=0.and.trafo) then
          i1 = 0
          do k = bk1,bk2
            ikpt1 = kpt1(k)
            ikpt2 = kptsum(ikpt,ikpt1)
            do i = bo1,bo2
              do j = bu1,bu2
                i1          = i1 + 1
                rdum        = ene(j,ikpt2,ispin2)-ene(i,ikpt1,ispin1)
                cprod(i1,1) = cprod(i1,1) * rdum
              enddo
            enddo
          enddo
        endif

        n = nu*no*nk
        allocate ( re(n) )
# ifndef INV
        allocate ( im(n) )
# endif

c
c       Loop over spectral-function frequencies and calculate BZ integration weights
        if(nfrqs/=0) then
          call timer_start('SUS weights')
          Allocate_ ( wtetra,(n,nfrqs) )
          call getkpt1p(kpt1p_,kpt1(bk1),nk,sym1,nsym1,nkpts(bk1),.false.)
          wtetra = 0
          if(gauss(1)/=0) then
            call gauss4_init      (wtetra,frqs,nfrqs,ikpt,kpt1(bk1),kpt1p_,nk,1,ispin1,ispin2,bo1,bo2,bu1,bu2)
          else
            call tetrahedron5_init(wtetra,frqs,nfrqs,ikpt,kpt1(bk1),kpt1p_,nk,1,ispin1,ispin2,bo1,bo2,bu1,bu2,bo1,bo2,bu1,bu2)
          endif
          wtetra = wtetra / nsym1
          call timer_stop('SUS weights',time=time1)
          write(OUTP,'(F8.2'NOA) time1 ; MpiR_WRITE
        else
          ! HILBERT NONE: simply sum over virtual transitions
          if(metal)          Error('System is a metal. Hilbert mesh needed (HILBERT).')
          if(any(frqc/=0d0)) Error('Real frequencies present. Hilbert mesh needed (HILBERT).')
          allocate ( wghtSr(n,nfrqh+nfrqc) )
          wghtSr = 0
          i      = 0
          do k = bk1,bk2
            ikpt1 = kpt1(k)
            ikpt2 = kptsum(ikpt1,ikpt)
            do i1 = bo1,bo2
              do i2 = bu1,bu2
                i = i + 1
                if(i2<=nband(ikpt2,ispin2)) then
                  rdum = ene(i2,ikpt2,ispin2) - ene(i1,ikpt1,ispin1)
                  do ifrq = 1,nfrqh
                    wghtSr(i,ifrq) = -2*rdum / ((frqh(ifrq)/img)**2 + rdum**2) / nkpt / nsym1 * nkpts(k)
                  enddo
                  do ifrq = 1,nfrqc
                    wghtSr(i,nfrqh+ifrq) = -2*rdum / ((frqc(ifrq)/img)**2 + rdum**2) / nkpt / nsym1 * nkpts(k)
                  enddo
                endif
              enddo
            enddo
          enddo
          if(i/=n) Bug('count error.')
        endif

        call timer_start('SUS suscep')

c
c       Spectral function
        if(trafo.and.wrtinfo.and.nfrqs/=0) then
          rdum = nsym1
          do ifrqS = 1,nfrqs
            if(multdiff_/=0.and.frqS(ifrqS)/=0) rdum = nsym1 / frqs(ifrqS)**2
            if(ikpt==1.and.ispin1==ispin2) then
              sf(ifrqS) = sf(ifrqS) + rdum*sum(wtetra(:,ifrqS)*(abs(dcprod(:,1))**2+abs(dcprod(:,2))**2+abs(dcprod(:,3))**2))/3
            else
              sf(ifrqS) = sf(ifrqS) + rdum*sum(wtetra(:,ifrqS)*abs(cprod(:n,1))**2)
            endif
          enddo
        endif

c
c          PREORDER: Preorder arrays wtetra and (d)cprod to enable use of ddot in spectral_function
c       NO PREORDER: Define array index of nonzero weights wtetra.
        if(nfrqs/=0) then
          Timing( call timer_start('SUS suscep-prep') )
          allocate ( nindex(nfrqs) )
# ifdef PREORDER
          allocate ( index(1,nfrqs), pnt(n) )
#   if PREORDER == 1
          i = 0
          do k = bk1,bk2
            ikpt1 = kpt1(k)
            ikpt2 = kptsum(ikpt1,ikpt)
            do i1 = bo1,bo2
              do i2 = bu1,bu2
                i     = i + 1
                re(i) = ene(i2,ikpt2,ispin2) - ene(i1,ikpt1,ispin1)
              enddo
            enddo
          enddo
          call rorderpf(pnt,re,n)
#   elif PREORDER == 2
          pnt = 0
          l   = 0
          do ifrqS = 1,nfrqs
            do i = 1,n
              if(wtetra(i,ifrqS)>wghtthr.and.pnt(i)==0) then
                l      = l + 1
                pnt(i) = l
              endif
            enddo
          enddo
          do i = 1,n
            if(pnt(i)==0) then
              l      = l + 1
              pnt(i) = l
            endif
          enddo
          if(l/=n) Bug('count error.')
          pnt(pnt) = [ (i,i=1,n) ] ! take the inverse
#   else
#     error 'unknown PREORDER type'
#   endif
          do ifrqS = 1,nfrqs
            wtetra(:,ifrqS) = wtetra(pnt,ifrqS)
            do i1 = 1,n
              if(wtetra(i1,ifrqS)>wghtthr) exit
            enddo
            do i2 = n,1,-1
              if(wtetra(i2,ifrqS)>wghtthr) exit
            enddo
            index(1,ifrqS) = i1
            nindex(ifrqS)  = i2 - i1 + 1
          enddo
#   ifndef SPECF_MSUM          
          if(any(nindex>0)) then
            frqs1  = minval( [ (i,i=1,nfrqs) ] , nindex>0 )
            nfrqs1 = maxval( [ (i,i=1,nfrqs) ] , nindex>0 ) - frqs1 + 1
          else
            frqs1  = 1
            nfrqs1 = 0
          endif
#   endif          
          do i = 1,dim
            cprod(:n,i) = cprod(pnt,i)
          enddo
          if(ikpt==1.and.ispin1==ispin2) then
            do i = 1,3
              dcprod(:n,i) = dcprod(pnt,i)
            enddo
          endif
          deallocate ( pnt )
# else
          allocate ( index(n,nfrqs) )
          do ifrqS = 1,nfrqs
            k = 0
            do i = 1,n
              if(wtetra(i,ifrqS)>wghtthr) then
                k = k + 1 ; index(k,ifrqS) = i ; wtetra(k,ifrqS) = wtetra(i,ifrqS)
              endif
            enddo
            nindex(ifrqS) = k
          enddo
# endif
          Timing( call timer_stop('SUS suscep-prep') )
        endif

c
c       ----------------------------------------
c       Calculate contribution of current packet
c       ----------------------------------------

        beginSplitNodes(1) ! split off active processes (to avoid deadlock)

        l = 0
        do ll = 0,2 ! loops 0->2 only if multdiff_==2

        if(multdiff_==2) call get_wghtS(ll)

        if(ikpt==1.and.ispin1==ispin2) then
          call timer_start('SUS headwings')
          call checkmem('head...',0d0)
c         (A0) Head (only Gamma point)
          if(multdiff_/=2.or.ll==2) then
            if(multdiff_==1) l = 2
            ij = 0
            do k2 = 1,3
              do k1 = 1,k2
                ij = ij + 1
                do k = 1,n
# ifdef INV
                  re(k) = dcprod(k,k1) * dcprod(k,k2)
# else
                  cdum  = dcprod(k,k1) * conjg(dcprod(k,k2)) ; re(k) = real(cdum) ; im(k) = imag(cdum)
# endif
                enddo
                call spectral_function(sus1,susc1 NoInvC(sus2) NoInvC(susc2) ,l)
# ifdef SPECF_MSUM
                MnoR( cycle )
# endif
                ij = (k2-1)*k2/2 + k1
                do ifrq = 1,nfrqh
                  headh(ij,ifrq) = headh(ij,ifrq) + sus1(ifrq) NoInv( + img * sus2(ifrq) )
                enddo
                do ifrq = 1,nfrqc
                  headc(k1,k2,ifrq) = headc(k1,k2,ifrq) + sus1(nfrqh+ifrq) NoInv( + img * sus2(nfrqh+ifrq) )
     &                                            + img * susc1(ifrq)      NoInv( -       susc2(ifrq)      )
                  if(k1/=k2)
     &            headc(k2,k1,ifrq) = headc(k2,k1,ifrq) + sus1(nfrqh+ifrq) NoInv( - img * sus2(nfrqh+ifrq) )
     &                                            + img * susc1(ifrq)      NoInv( +       susc2(ifrq)      )
                enddo
              enddo
            enddo
          endif
c         (A1) Wings (only Gamma point)
          if(multdiff_/=2.or.ll==1) then
            call checkmem('wings...',0d0)
            if(multdiff_==1) l = 1
            do k2 = 1,3
              do i = 1,dim
                do k = 1,n
# ifdef INV
                  re(k) = cprod(k,i) * dcprod(k,k2)
# else
                  cdum  = cprod(k,i) * conjg(dcprod(k,k2)) ; re(k) = real(cdum) ; im(k) = imag(cdum)
# endif
                enddo
                call spectral_function(sus1,susc1 NoInvC(sus2) NoInvC(susc2) ,l)
# ifdef SPECF_MSUM
                MnoR( cycle )
# endif
                do ifrq = 1,nfrqh
                  wingh(k2,i,ifrq) = wingh(k2,i,ifrq) + sus1(ifrq) NoInv( + img * sus2(ifrq) )
                enddo
                do ifrq = 1,nfrqc
                  wingc(k2,i,1,ifrq) = wingc(k2,i,1,ifrq) + sus1(nfrqh+ifrq) NoInv( + img * sus2(nfrqh+ifrq) )
     &                                              + img * susc1(ifrq)      NoInv( -       susc2(ifrq)      )
                  wingc(k2,i,2,ifrq) = wingc(k2,i,2,ifrq) + sus1(nfrqh+ifrq) NoInv( - img * sus2(nfrqh+ifrq) )
     &                                              + img * susc1(ifrq)      NoInv( +       susc2(ifrq)      )
                enddo
              enddo
            enddo
          endif
          call timer_stop('SUS headwings')
        endif

c       (B) Body

        call checkmem('body...',0d0)
        if(.not.trafo.and.ll/=0) cycle ! .not.trafo : for pure MB functions (NOAPW) only body (multdiff_ not used)

# if defined(MPI) && !defined(SPECF_MSUM)
        mm = 0
        do j = 1,dim
          do i = 1,j
            if(trafo.and.multdiff_==2) then      ! multdiff_ = 2
              if(ll==0.and.(i==1.or.j==1)) cycle !   only body,          loop i = 2,j
              if(ll> 0.and.j> 1)           cycle !   only head or wings, loop i = 1,1
            endif
            if(sirrep(i,j)<sirrep_threshold) cycle
            mm = mm + 1
          enddo
        enddo
        do rank = 0,Nsize-1
          m  = 0
          m1 =   1 + mod(rank+Nrank,Nsize)   * mm / Nsize
          m2 = ( 1 + mod(rank+Nrank,Nsize) ) * mm / Nsize
# endif

        l = 0
        do j = 1,dim
          do i = 1,j
            if(trafo.and.multdiff_==2) then      ! multdiff_ = 2
              if(ll==0.and.(i==1.or.j==1)) cycle !   only body,          loop i = 2,j
              if(ll> 0.and.j> 1)           cycle !   only head or wings, loop i = 1,1
            endif
            if(sirrep(i,j)<sirrep_threshold) cycle
# ifndef SPECF_MSUM            
            Mpi( m = m + 1 ; if(m<m1.or.m>m2) cycle )
# endif
            if(trafo.and.multdiff_==1) l = count([i,j]==1)
            Timing( call timer_start('SUS suscep-cprods:'//Chr(n)) )
            do k = 1,n
# ifdef INV
              re(k) = cprod(k,i) * cprod(k,j)
# else
              cdum  = cprod(k,i) * conjg(cprod(k,j)) ; re(k) = real(cdum) ; im(k) = imag(cdum)
# endif
            enddo
            Timing( call timer_stop('SUS suscep-cprods:'//Chr(n)) )
            Timing( call timer_start('SUS specf') )
            call spectral_function(sus1,susc1 NoInvC(sus2) NoInvC(susc2) ,l)
            Timing( call timer_stop('SUS specf') )
# ifdef SPECF_MSUM
            MnoR( cycle )
# endif
            Timing( call timer_start('SUS suscep-distrib') )
            ij = j*(j-1)/2 + i
# ifdef ACCESS_CONTIG
            do ifrq = 1,nfrqh
              susceph_(ifrq,ij) = susceph_(ifrq,ij) + sus1(ifrq) NoInv( + img * sus2(ifrq) )
            enddo
            ij = (j-1)**2 + 2*i - 1
            do ifrq = 1,nfrqc
              suscepc_(ifrq,ij) = suscepc_(ifrq,ij) + sus1(nfrqh+ifrq) NoInv( + img * sus2(nfrqh+ifrq) )
     &                                        + img * susc1(ifrq)      NoInv( -       susc2(ifrq)      )          
            enddo
            if(i/=j) then
              ij = ij + 1
              do ifrq = 1,nfrqc
                suscepc_(ifrq,ij) = suscepc_(ifrq,ij) + sus1(nfrqh+ifrq) NoInv( - img * sus2(nfrqh+ifrq) )
     &                                          + img * susc1(ifrq)      NoInv( +       susc2(ifrq)      )
              enddo
            endif
# else
            do ifrq = 1,nfrqh
              susceph(ij,ifrq)  = susceph(ij,ifrq) + sus1(ifrq) NoInv( + img * sus2(ifrq) )
            enddo
            do ifrq = 1,nfrqc
              cdum              = cmplx(  sus1(nfrqh+ifrq) , susc1(ifrq)      )
              NoInv( cdum1      = cmplx( -susc2(ifrq)      , sus2(nfrqh+ifrq) ) )
              suscepc(i,j,ifrq) = suscepc(i,j,ifrq) + cdum NoInv( +cdum1 )
              NoInv( if(i/=j) suscepc(j,i,ifrq) = suscepc(j,i,ifrq) + cdum-cdum1 )
            enddo
# endif
            Timing( call timer_stop('SUS suscep-distrib') )
          enddo
        enddo

        Timing( call timer_start('SUS suscep-Nfence') )
        NFENCE
        Timing( call timer_stop('SUS suscep-Nfence') )

# if defined(MPI) && !defined(SPECF_MSUM)
        enddo
# endif
        
        if(multdiff_/=2) exit
        enddo
      
        Deallocate_ ( cprod )
        deallocate  ( re )
        if(nfrqs==0)           deallocate ( wghtSr )
        if(allocated(im))      deallocate ( im     )
        if(allocated(index))   deallocate ( index  )
        if(allocated(nindex))  deallocate ( nindex )
        if(allocated(wtetra)) tDeallocate ( wtetra )
        if(allocated(dcprod))  deallocate ( dcprod )

        call timer_stop('SUS suscep',time=time1)
        write(OUTP,'(F8.2)') time1 ; MpiR_WRITE

        endSplitNodes ! end split active processes

        Mpi ( call Mwrite(line) ; line = ' ')

        NFENCE

      bu1 = bu2 + 1
      enddo
      bo1 = bo2 + 1
      enddo
      bk1 = bk2 + 1
      enddo ! End loop over packets

# ifdef MPI
c     Inactive processes "simulate" a packet to avoid deadlocks      
 5    do d = 1,nskip  ! inactive processes (avoid deadlock)
        beginSplit(0) ! around wavefproducts_mt
        endSplit
        beginSplit(0) ! around wavefproducts_pw
        endSplit
        m = Nsize
        beginSplitNodes(0) ! around suscep summation
        do ll = 0,2
          do rank = 1,m-Nsize
            NFENCE
          enddo
          if(multdiff_/=2) exit
        enddo
        endSplitNodes
        call Mwrite('skipline')
        NFENCE
      enddo
# endif

      first_term = .false.

c
c     Write spectral function (wrtinfo)
      if(trafo.and.wrtinfo.and.nfrqs/=0) then
        Mpi( call Msum(sf,0) )
        Rbegin
        i = fopen('spex.sf',status='new',numbered=.true.)
        write(i,'(A)')               '# Spectral function on Hilbert frequency mesh (head element of susceptibility).'
        write(i,'(A,I4,A,3F10.5,A)') '# k point:',ikpt,'(',kpt(:,ikpt),')'
        write(i,'(A,I4)')            '# spin:',spin
        do ifrqS = 1,nfrqs
          write(i,'(2F20.10)') frqs(ifrqS),sf(ifrqS)
        enddo
        call fclose(i)
        Rend
        deallocate(sf)
      endif

      Rif(wrtinfo.and..not.trafo) Info('Files ''spex.sf'' not written because trafo unset.')
      Rif(wrtinfo.and.nfrqs==0)   Info('Files ''spex.sf'' not written because of HILBERT NONE.')

      Deallocate_(sus1)
      NoInv( Deallocate_(sus2) )
      Deallocate_(susc1)
      NoInv( Deallocate_(susc2) )
      if(allocated(frqs))   deallocate ( frqs   )
      if(allocated(wghtSr)) deallocate ( wghtSr )
      if(allocated(wghtSi)) deallocate ( wghtSi )

c
c     Loop over spins if necessary
 20   if(.not.l_soc.and.spin==0) then
        if(nspin==2) then
          if(ispin1==1) then ; ispin1 = 2 ; ispin2 = 2 ; goto 2 ; endif
        else
          NFENCE
          Oif(nfrqh>0) susceph = susceph * 2
          Oif(nfrqc>0) suscepc = suscepc * 2
          if(ikpt==1) then
            if(nfrqh>0) then ; headh = headh * 2 ; wingh = wingh * 2 ; endif
            if(nfrqc>0) then ; headc = headc * 2 ; wingc = wingc * 2 ; endif
          endif
        endif
      endif
      if(spin==3) then
        if(ispin1==2) then ; ispin1 = 1 ; ispin2 = 2 ; goto 2 ; endif
      else if(spin==4) then
        if(ispin1==1) then ; ispin1 = 2 ; ispin2 = 1 ; goto 2 ; endif
      endif

      NFENCE  

# ifdef ACCESS_CONTIG
c
c     Reorder susceph(c) from frequency-major into element-major order
      call timer_start('SUS reorder')
      if(nfrqh>0) then
        Ocall reshape_permute_ ifInv(r,c) (susceph,[nfrqh,dim*(dim+1)/2],[2,1])
        Nfence(susceph)
        nullify(susceph_)
      endif        
      if(nfrqc>0) then
        Ocall reshape_permute_c(suscepc,[nfrqc,dim**2],[2,1])
        Nfence(suscepc)
        nullify(suscepc_)
        do ifrq = Nrange(1,nfrqc)
          call array_resort_c(suscepc(:,:,ifrq),dim**2,get_index_suscepc,[dim])
        enddo
      endif
      call timer_stop('SUS reorder')              
# elif defined(INV)
c
c     Add lower matrix triangle (suscepc is symmetric for INV)      
      if(nfrqc>0) then
        call timer_start('SUS transpose')
        do ifrq = Nrange(1,nfrqc)
          do j = 1,dim
            do i = j+1,dim
              suscepc(i,j,ifrq) = suscepc(j,i,ifrq)
            enddo
          enddo
        enddo
        call timer_stop('SUS transpose')      
      endif
# endif      
      

      NFENCE

c
c     Calculate Drude term (plasma frequency)
      if(ikpt==1.and.metal.and.spin<=2) then
        Rcall timer_start('SUS plasma')
        if(symon) then ; call getkpt1       (kpt1,nkpt1,nkpts,sym1,nsym1,ikpt,0,.false.)
        else           ; call getkpt1_fullBZ(kpt1,nkpt1,nkpts,sym1,nsym1,0)
        endif
        ispin1 = 1 ; if(spin==2) ispin1 = 2
        do
          call tetrahedron4_init(wdrude,ispin1)
          k1 = 1     ; Mpi( k1 = Mrank*nkpt1/Msize+1   )
          k2 = nkpt1 ; Mpi( k2 = (Mrank+1)*nkpt1/Msize )
          Load( call read_wavef3(1,bando,kpt1,nkpt1,ispin1) )
          allocate ( momentum(bando,bando,k1:k2,3) )
          call momentum_matrix(momentum,kpt1(k1:k2),k2-k1+1,ispin1,ispin1,1,bando,1,bando,.false. MpiC(.false.) LoadC(.false.) )
          Load( deallocate ( cmt,cpw ) )
          integrand = 0
          ! loop over "little group" of k points
          allocate ( re(3) )
          do ikpt1 = k1,k2
            ikpt2 = kpt1(ikpt1)
            i1    = 1
            rdum  = 1
            do while(i1<=bando)
              if(lhub) rdum = 1 - hub(i1,kptp(ikpt2),ispin1)**2
              i2    = deg(i1,ikpt2,ispin1) ; if(i2>bando) exit
              re(1) = sum( [ (momentum(i,i,ikpt1,1),i=i1,i2) ] ) / (i2-i1+1) ! take the trace of the momentum matrix which is invariant under rotation in the degenerate subspace (to make the plasma fr
              re(2) = sum( [ (momentum(i,i,ikpt1,2),i=i1,i2) ] ) / (i2-i1+1) !
              re(3) = sum( [ (momentum(i,i,ikpt1,3),i=i1,i2) ] ) / (i2-i1+1) !
              integrand(i1:i2,ikpt2) = sum(abs(re)**2) * rdum
              i1                     = i2 + 1
            enddo
            ! add symmetry-equivalent k points
            do isym = 2,nsym1
              ikpt2              = kptsym(kpt1(ikpt1),sym1(isym))
              integrand(:,ikpt2) = integrand(:,kpt1(ikpt1))
            enddo
          enddo
          deallocate ( momentum,re )
          Mpi( call Msum(integrand) )
          if(sum(abs(integrand))<1d-8) then
            RWarn('All integrand points are nearly zero.')
            Rwrite(0,'(A)') '                Possibly only high-symmetry points in the BZ.'
            Rwrite(0,'(A)') '                No contribution to plasma frequency calculated.'
          else
            plasma = plasma + 4*pi/3/vol * sum(wdrude*integrand)
          endif
          if(spin/=0.or.l_soc.or.nspin==1.or.ispin1==2) exit
          ispin1 = 2
        enddo
        if(spin==0.and..not.l_soc.and.nspin==1) plasma = plasma * 2
        if(plasma<-1d-10) Error('Plasma frequency undefined.')
        plasma = max(0d0,plasma)
        Rcall timer_stop('SUS plasma',time=time1)
        Rwrite(6,'(A,F8.2,A)') 'Plasma frequency timing:',time1
      endif

c
c     Reduce from processes
      NFENCE
      Mpi( call Mreduce )
      NFENCE

c
c     IBC : Hermitianize and add to suscep
      Rif(oibc/=0) call add_to_suscep

c
c     Add contribution of symmetry-equivalent q points
      if(nsym1>1) then
        Rwrite(6,'(A'NoA) 'Add contribution of symmetry-equivalent k points... '
        Rcall timer_start('SUS symmetry-equiv')
c       (A) Body
        if(trafo) then
          ! Case: transformation to Coulomb eigenvectors (use contracted irreps calculated above)
c Define OLD_ALGO for old algorithm. (Kept for testing.)
c# define OLD_ALGO
# ifdef OLD_ALGO
#   warning OLD_ALGO defined
          Allocate_(suscep0,(dim,dim))
          Allocate_(suscep1,(dim,dim))
# endif
          do ifrq = Nrange(1,nfrqh+nfrqc)
# ifdef OLD_ALGO
            if(ifrq<=nfrqh) then ; suscep0 = unpackmat(susceph(:,ifrq))
            else                 ; suscep0 = suscepc(:,:,ifrq-nfrqh)
            endif
            suscep1 = 0
# endif
            i2      = 0
            l       = 0
            do k2 = 1,nirrep ; d2 = dirrep(k2) ; i1 = 0               
# ifdef OLD_ALGO
              do k1 = 1,k2 ; d1 = dirrep(k1) ; n  = (d1*d2)**2
                suscep1(i1+1:i1+d1,i2+1:i2+d2)   = multiply_contract( irrep_contr(l+1:l+n) ,  suscep0(i1+1:i1+d1,i2+1:i2+d2)  )
                if(.not.trsoff) then
                  suscep1(i1+1:i1+d1,i2+1:i2+d2) = suscep1(i1+1:i1+d1,i2+1:i2+d2) +
     &                                    multiply_contract( irrep_contt(l+1:l+n) , transpose(suscep0(i2+1:i2+d2,i1+1:i1+d1)) )
                endif
                if(k1/=k2) then
                  suscep1(i2+1:i2+d2,i1+1:i1+d1)   = transpose(
     &                            multiply_contract( irrep_contr(l+1:l+n) ,  transpose(suscep0(i2+1:i2+d2,i1+1:i1+d1))  ) )
                  if(.not.trsoff) then
                    suscep1(i2+1:i2+d2,i1+1:i1+d1) = suscep1(i2+1:i2+d2,i1+1:i1+d1) + transpose(
     &                                      multiply_contract( irrep_contt(l+1:l+n) , suscep0(i1+1:i1+d1,i2+1:i2+d2) ) )
                  endif
                endif
# else
              ! get matrix slice (->suscep0,suscep0_)
              allocate( suscep0(i2+d2,d2) )
              allocate( suscep1(i2+d2,d2) )
              if(ifrq<=nfrqh) then
                do i = 1,d2
                  ii               = (i2+i) * (i2+i-1) / 2
                  suscep0(:i2+i,i) = susceph(ii+1:ii+i2+i,ifrq)
                enddo
                do i = 1,d2     ; ii2 = i2  + i
                  do j = 1,d2-i ; ii1 = ii2 + j
                    ii                = ii1*(ii1-1)/2 + ii2
                    suscep0(i2+i+j,i) = MCONJG(susceph(ii,ifrq))
                  enddo
                enddo                
              else
                allocate( suscep0_(i2+d2,d2) )
                allocate( suscep1_(i2,   d2) )                
                suscep0  =           suscepc(    :i2+d2 , i2+1:i2+d2 , ifrq-nfrqh)
                suscep0_ = transpose(suscepc(i2+1:i2+d2 ,     :i2+d2 , ifrq-nfrqh))
              endif
              ! symmetrization for current frequency ifrq
              do k1 = 1,k2 ; d1 = dirrep(k1) ; n = (d1*d2)**2
                if(sirrep(i1+1,i2+1)<sirrep_threshold) then
                  suscep1(i1+1:i1+d1,:) = 0
                  if(ifrq>nfrqh) suscep1_(i1+1:i1+d1,:) = 0
                else
                  suscep1(i1+1:i1+d1,:) =  multiply_contract( irrep_contr(l+1:l+n) , suscep0(i1+1:i1+d1,:) )
                  if(.not.trsoff) then
                    if(ifrq<=nfrqh) then ; suscep1(i1+1:i1+d1,:) = suscep1(i1+1:i1+d1,:) +
     &                                     multiply_contract( irrep_contt(l+1:l+n) , MCONJG(suscep0 (i1+1:i1+d1,:)) )
                    else                 ; suscep1(i1+1:i1+d1,:) = suscep1(i1+1:i1+d1,:) +
     &                                     multiply_contract( irrep_contt(l+1:l+n) ,        suscep0_(i1+1:i1+d1,:)  )
                    endif
                  endif
                  if(ifrq>nfrqh.and.k1/=k2) then
                    suscep1_(i1+1:i1+d1,:) = multiply_contract( MCONJG(irrep_contr(l+1:l+n)) , suscep0_(i1+1:i1+d1,:) )
                    if(.not.trsoff)
     &              suscep1_(i1+1:i1+d1,:) = suscep1_(i1+1:i1+d1,:) +
     &                                       multiply_contract( MCONJG(irrep_contt(l+1:l+n)) , suscep0 (i1+1:i1+d1,:) )
                  endif
                endif
# endif                                
                i1 = i1 + d1
                l  = l + n
              enddo
# ifndef OLD_ALGO
              if(ifrq<=nfrqh) then
                do i = 1,d2
                  m                        = (i2+i) * (i2+i-1) / 2
                  susceph(m+1:m+i2+i,ifrq) = suscep1(:i2+i,i)
                enddo
              else
                suscepc(:i2+d2,i2+1:i2+d2,ifrq-nfrqh) =           suscep1
                suscepc(i2+1:i2+d2,:i2,   ifrq-nfrqh) = transpose(suscep1_)
              endif              
              deallocate ( suscep0 )
              deallocate ( suscep1 )              
              if(ifrq>nfrqh) then
                deallocate ( suscep0_ )
                deallocate ( suscep1_ )
              endif              
# endif
              i2 = i2 + d2
            enddo
            if(i1/=dim.or.i2/=dim) Bug('Count error.')
# ifdef OLD_ALGO
            if(ifrq<=nfrqh) then ; susceph(:,ifrq)         = packmat(suscep1))
            else                 ; suscepc(:,:,ifrq-nfrqh) = suscep1
            endif
# endif
          enddo
          Ndeallocate ( irrep_contr )
          if(.not.trsoff) tNdeallocate ( irrep_contt )
# ifdef OLD_ALGO
          Deallocate_(suscep0)
          Deallocate_(suscep1)
# endif
        else
          ! Case: No transformation to Coulomb eigenvectors (use mtrafo for body)
          Allocate_(suscep0,(dim,dim))
          Allocate_(suscep1,(dim,dim))
          do ifrq = Nrange(1,nfrqh+nfrqc)
            if(ifrq<=nfrqh) then ; suscep0 = unpackmat(susceph(:,ifrq))
            else                 ; suscep0 = suscepc(:,:,ifrq-nfrqh)
            endif
            suscep1 = suscep0
            do isym = 2,nsym1
              call mtrafo(suscep1,suscep0,dim,dim,ikpt,-sym1(isym),3,.true.)
            enddo
            if(ifrq<=nfrqh) then ; susceph(:,ifrq)         = packmat(suscep1)
            else                 ; suscepc(:,:,ifrq-nfrqh) = suscep1
            endif
          enddo
          Deallocate_(suscep0)
          Deallocate_(suscep1)
        endif
c       (B) Head and wings
        if(ikpt==1) then
          allocate ( wing0(3,dim,2),wing1(3,dim,2),wing2(3,dim,2) )
          do ifrq = 1,nfrqh+nfrqc
            if(ifrq<=nfrqh) then
              head0        = unpackmat(headh(:,ifrq))
              wing0(:,:,1) = wingh(:,:,ifrq)
              wing0(:,:,2) = MCONJG(wingh(:,:,ifrq))
            else
              head0        = headc(:,:,  ifrq-nfrqh)
              wing0        = wingc(:,:,:,ifrq-nfrqh)
            endif
            head2 = head0
            wing2 = wing0
            do isym = 2,nsym1
              if(trafo) then
                do i = 1,3
                  n = 0
                  do k = 1,nirrep
                    j                  = dirrep(k)
                    wing1(i,n+1:n+j,1) = matmul(conjg(irrep(:j,:j,k,isym)),wing0(i,n+1:n+j,1))
                    wing1(i,n+1:n+j,2) = matmul(wing0(i,n+1:n+j,2),transpose(irrep(:j,:j,k,isym)))
                    n                  = n + j
                  enddo
                enddo
              else
                do i = 1,3
                  call mtrafo(wing1(i,:,1),wing0(i,:,1),nbasm(ikpt),1,ikpt,-sym1(isym),1,.false.)
                  call mtrafo(wing1(i,:,2),wing0(i,:,2),1,nbasm(ikpt),ikpt,-sym1(isym),2,.false.)
                enddo
              endif
              i   = sym(sym1(isym))%inv
              rot = matmul(rlat,matmul(sym(i)%rrot,transpose(lat)))/(2*pi)
              if(sym1(isym)>nsymt) then ! time-reversal symmetry
                head2        = head2        + transpose ( matmat(rot,matmat(head0,transpose(rot))) )
                wing2(:,:,1) = wing2(:,:,1) + matmat(rot,wing1(:,:,2))
                wing2(:,:,2) = wing2(:,:,2) + matmat(rot,wing1(:,:,1))
              else
                head2        = head2        +             matmat(rot,matmat(head0,transpose(rot)))
                wing2(:,:,1) = wing2(:,:,1) + matmat(rot,wing1(:,:,1))
                wing2(:,:,2) = wing2(:,:,2) + matmat(rot,wing1(:,:,2))
              endif
            enddo
            if(ifrq<=nfrqh) then
              headh(:,  ifrq) = packmat(head2)
              wingh(:,:,ifrq) = wing2(:,:,1)
            else
              headc(:,:,  ifrq-nfrqh) = head2
              wingc(:,:,:,ifrq-nfrqh) = wing2
            endif
          enddo
          deallocate ( wing0,wing1,wing2 )
        endif
        Rcall timer_stop('SUS symmetry-equiv',time=time1)
        Rwrite(6,'(A,F7.2,A)') 'done   ( Timing: ',time1,' )'
        if(trafo) tDeallocate ( irrep )
        NFENCE
      endif

c
c     Add core contribution
      if(any(cores)) then
        Rcall timer_start('SUS add_core')
        call add_core
        Rcall timer_stop('SUS add_core')
      endif

# if 0
      do ifrq = 1,nfrqh
        allocate ( re(dim) )
        call diagonalize(re,susceph(:,ifrq))
        write(*,*) ifrq
        write(*,*) re
        write(*,*)
c        read(*,*)
        deallocate ( re )
      enddo
      do ifrq = 1,nfrqh
        allocate ( re(dim) )
        call diagonalize(re,susceph(:,ifrq))
        write(*,*) ifrq,sum(re)
        deallocate ( re )
      enddo
c      read(*,*)
      Error(' ')
# endif

c# ifdef MPI
c      if(ikpt==1) then
c        if(nfrqh/=0) then
c          call Mcast(headh)
c          call Mcast(wingh)
c        endif
c        if(nfrqc/=0) then
c          call Mcast(headc)
c          call Mcast(wingc)
c        endif
c      endif
c# endif

      if(ikpt==1.and.metal.and.plasma/=0) then
        plasma = sqrt(plasma)
        Rwrite(6,'(A,F9.6,A'NoA) 'Plasma frequency: ',27.2113831243217d0*plasma,' eV'
        if(plasma*hartree<1d-6) then
          plasma = 0
          Rwrite(6,'(A'NoA) ' < 10^-6 eV, set to zero'
        endif
        Rwrite(6,*)
      endif

      Rwrite(6,'(/A/)') '---------'

      if(allocated(ibc_suscep1)) deallocate ( ibc_suscep1 )
      Ndeallocate(sirrep)

# if 0
      Rbegin
      if(nfrqh/=0) then
        do ifrq = 1,nfrqh
          write(*,*) 's', sum(abs(susceph(:,ifrq))**2)
          if(ikpt==1) then
            write(*,*) 'h',sum(abs(headh(:,ifrq))**2)
            write(*,*) 'w',sum(abs(wingh(:,:,ifrq))**2)
          endif
        enddo
      endif
      do ifrq = 1,nfrqc
        write(*,*) ifrq,'s',sum(abs(suscepc(:,:,ifrq))**2)
        if(ikpt==1) then
          write(*,*) 'h',sum(abs(headc(:,:,ifrq))**2)
          write(*,*) 'w',sum(abs(wingc(:,:,:,ifrq))**2)
        endif
      enddo
      Rread(*,*)
      Rend
      Mpi( call mpi_barrier(Mcomm,i) )
# endif

      firstcall = .false.

      call timer_stop('Routine susceptibility SUS')

      if(otimer==3) then
        call timer_print(['SUS '],select=.true.,title='TIMING SUSCEPTIBILITY',  prnt=.true. andR )
        call timer_print(['WFP '],select=.true.,title='TIMING WAVEFPRODUCTS_MT',prnt=.true. andR )
      else if(otimer==4) then
        call timer_print(['SUS '],select=.true.,title='TIMING SUSCEPTIBILITY',  prnt=.true.,unit= ifMpi(Mrank,0) )
        call timer_print(['WFP '],select=.true.,title='TIMING WAVEFPRODUCTS_MT',prnt=.true.,unit= ifMpi(Mrank,0) )
      endif

      contains

c     ---------------

      subroutine get_wghtS(l)
      implicit none
      integer, intent(in) :: l
      integer             :: term
      term = 0 ; if(spin==3.or.spin==4) term = iterm
      call getwghtS(wghtSr,           wghtSi,           frqh,nfrqh,frqs,nfrqs,disorder,l,term)
      if(nfrqc>0)
     &call getwghtS(wghtSr(:,nfrqh+1),wghtSi(:,nfrqh+1),frqc,nfrqc,frqs,nfrqs,disorder,l,term)
      if(spin==3.or.spin==4) then
        wghtSr = -wghtSr
        wghtSi = -wghtSi
      endif
      end subroutine get_wghtS

c     ---------------

c     Integrate spectral function over frqs and return contributions to susceptibility (sus1,sus2,susc1,susc2)
c     SPECF_MSUM unset (default): Each process calculates its contribution and returns it
c     SPECF_MSUM set:             Contributions are collectively summed (Msum) and the sum is returned in rank 0
      subroutine spectral_function(sus1, susc1 NoInvC(sus2) NoInvC(susc2) ,mm)
      implicit none
      integer, intent(in)  :: mm
      real_dp, intent(out) :: sus1(nfrqh+nfrqc) NoInvC(sus2(nfrqh+nfrqc))
      real_dp, intent(out) :: susc1(nfrqc)      NoInvC(susc2(nfrqc))
      real_dp              :: reS(nfrqs) NoInvC(imS(nfrqs))
      integer              :: ifrqs,ifrq,i1
      real_dp              :: rdum1 NoInvC(rdum2)
      real_dp              :: ddot
# ifndef PREORDER
      integer              :: k
# endif

      ! HILBERT NONE (simple summation)
      if(nfrqs==0) then
        do ifrq = 1,nfrqh+nfrqc
                 sus1(ifrq) = ddot(n,re,1,wghtSr(1,ifrq),1)
          NoInv( sus2(ifrq) = ddot(n,im,1,wghtSr(1,ifrq),1) )
        enddo
# ifdef SPECF_MSUM
        call Msum(sus1,rank=0) ; NoInv( call Msum(sus2,rank=0) )
# endif
        return
      endif

      ! Calculate spectral function for HILBERT
# ifdef PREORDER
      do ifrqs = 1,nfrqs
        if(nindex(ifrqs)>0) then
                 i1         = index(1,ifrqs)
                 reS(ifrqs) = ddot(nindex(ifrqs),re(i1),1,wtetra(i1,ifrqs),1)
          NoInv( imS(ifrqs) = ddot(nindex(ifrqs),im(i1),1,wtetra(i1,ifrqs),1) )
        else
                 reS(ifrqs) = 0
          NoInv( imS(ifrqs) = 0 )
        endif
      enddo
# else
      do ifrqs = 1,nfrqs
               rdum1 = 0
        NoInv( rdum2 = 0 )
        do k = 1,nindex(ifrqs)
                 i1    = index(k,ifrqs)
                 rdum1 = rdum1 + re(i1) * wtetra(k,ifrqs)
          NoInv( rdum2 = rdum2 + im(i1) * wtetra(k,ifrqs) )
        enddo
               reS(ifrqs) = rdum1
        NoInv( imS(ifrqs) = rdum2 )
      enddo
# endif

      if(mm/=0) then
        do ifrqs = 1,nfrqs
          if(frqs(ifrqs)>0) then
                   rdum1      = frqs(ifrqs)**mm
                   reS(ifrqs) = reS(ifrqs) / rdum1
            NoInv( imS(ifrqs) = imS(ifrqs) / rdum1 )
          endif
        enddo
      endif

# ifdef SPECF_MSUM
      call Msum(reS) ; NoInv( call Msum(imS) )
      sus1  = 0      ; NoInv( sus2  = 0 )
      susc1 = 0      ; NoInv( susc2 = 0 )
      do ifrq = Mrange1(nfrqh+nfrqc)
               sus1(ifrq)  = ddot(nfrqs,reS,1,wghtSr(1,ifrq),1)
        NoInv( sus2(ifrq)  = ddot(nfrqs,imS,1,wghtSr(1,ifrq),1) )
      enddo
      do ifrq = Mrange1(nfrqc)
               susc1(ifrq) = ddot(nfrqs,reS,1,wghtSi(1,nfrqh+ifrq),1)
        NoInv( susc2(ifrq) = ddot(nfrqs,imS,1,wghtSi(1,nfrqh+ifrq),1) )
      enddo
                         call Msum(sus1, rank=0) ; NoInv( call Msum(sus2, rank=0)  )
      if(nfrqc>0) then ; call Msum(susc1,rank=0) ; NoInv( call Msum(susc2,rank=0); ) endif
# else 
      do ifrq = 1,nfrqh+nfrqc
               sus1(ifrq)  = ddot(nfrqs1,reS(frqs1),1,wghtSr(frqs1,ifrq),1)
        NoInv( sus2(ifrq)  = ddot(nfrqs1,imS(frqs1),1,wghtSr(frqs1,ifrq),1) )
      enddo
      do ifrq = 1,nfrqc
               susc1(ifrq) = ddot(nfrqs1,reS(frqs1),1,wghtSi(frqs1,nfrqh+ifrq),1)
        NoInv( susc2(ifrq) = ddot(nfrqs1,imS(frqs1),1,wghtSi(frqs1,nfrqh+ifrq),1) )
      enddo
# endif

      end subroutine spectral_function

c     ---------------

      function irrep_contract(mode)
      implicit none
      MCOMPLEX_dp            :: irrep_contract((d1*d2)**2)
      complex_dp             :: cdum
      integer,    intent(in) :: mode
      integer                :: i,j,k,l,m
      do i = 1,d1
        do j = 1,d2
          do k = 1,d1
            do l = 1,d2
              if     (mode==0) then ; cdum = sum ( conjg(irrep(i,k,k1,:)) * irrep(j,l,k2,:)                          )
              else if(mode==1) then ; cdum = sum ( conjg(irrep(i,k,k1,:)) * irrep(j,l,k2,:) , sym1(:nsym1)<=nsym/2 )
              else                  ; cdum = sum ( conjg(irrep(j,l,k2,:)) * irrep(i,k,k1,:) , sym1(:nsym1)>nsym/2 )
              endif
# ifdef INV
              if(abs(imag(cdum))>1d-10) Error('Nonzero imaginary part in contracted irrep.')
# endif
              m                 = i+d1*(j-1)+d1*d2*(k-1)+d1*d1*d2*(l-1) ; if(m>(d1*d2)**2) Error('bug?')
              irrep_contract(m) = cdum
            enddo
          enddo
        enddo
      enddo
      end function irrep_contract

c     ---------------

      function multiply_contract(irrep_contr,suscep0)
      implicit none
      complex_dp              :: multiply_contract(d1,d2)
      MCOMPLEX_dp, intent(in) :: irrep_contr((d1*d2)**2)
      MCOMPLEX_dp             :: contract(d1*d2,d1*d2)
      complex_dp,  intent(in) :: suscep0(d1,d2)
      complex_dp              :: suscep(d1*d2)
      suscep            = reshape ( suscep0                      , [ d1*d2 ] )
      contract          = reshape ( irrep_contr                  , [ d1*d2,d1*d2 ] )
      multiply_contract = reshape ( matmul ( contract , suscep ) , [ d1,d2 ] )
      end function multiply_contract

c     ---------------

c
c     IBC : Calculate double-counting term (->ibc_suscep1)
      subroutine double_counting
      use, intrinsic :: iso_fortran_env
      implicit none
      integer :: i,j,k,i1,i2
      call timer_start('SUS IBC double_counting')
# ifdef FIBC
      call ibc1_doublecount(ibc_suscep1,cprod,dim,ikpt,kpt1,nkpt1,nkpts,nfrqh,bo1,bo2,bu1,bu2)
# else
      call ibc1_doublecount(ibc_suscep1,cprod,    ikpt,kpt1,nkpt1,nkpts,nfrqh,bo1,bo2,bu1,bu2)
# endif
      call timer_stop('SUS IBC double_counting',time=time1)
      write(6,'(F8.2'NoA) time1
      if(bu2<bandu) then
        deallocate ( cprod )
        return
      endif
      if(bu1<bandu) then
        i = 0
        j = 0
        do k = 1,nkpt1 ; do i1 = bo1,bo2 ; do i2 = bu1,bu2
          i = i + 1
          if(i2>=bandu) then
            j          = j + 1
            cprod(j,:) = cprod(i,:)
          endif
        enddo ; enddo ; enddo
        bu1 = bandu
        nu  = bu2-bu1+1
      endif
      end subroutine double_counting

c     ---------------

c
c     IBC : Hermitianize and add to suscep
      subroutine add_to_suscep
      use, intrinsic :: iso_fortran_env
      implicit none
      complex_dp,  allocatable :: ctrafo2(:,:)
      integer                  :: ifrq,itype,ieq,ic,i,n
# ifdef INV
#   define CTRAFO ctrafo2
# endif
# ifdef FIBC
      complex_dp,  allocatable :: suscep2(:,:)
      write(6,'(A'NoA) 'Add IBC to susceptibility... '
      call timer_start('SUS IBC add_to_suscep')
#   ifdef INV
      if(trafo) then
        allocate ( ctrafo2(nbasp,dim) )
        ctrafo2 = ctrafo(:nbasp,:)
        call desymmetrize(ctrafo2,nbasp,dim,1)
      endif
#   endif
      allocate ( suscep2(dim,nbasp) )
      do ifrq = 1,nfrqh
        ! Transform to Coulomb eigenvectors from left (case trafo)
        suscep2 = 0
        i       = 0
        ic      = 0
        do itype = 1,ntype
          n = sum ( [ ((2*l+1)*nindxm(l,itype),l=0,lcutm(itype)) ] )
          do ieq = 1,neq(itype)
            ic = ic + 1
            if(trafo) then
              suscep2(:,i+1:i+n)       = matmul ( transpose(conjg(CTRAFO(i+1:i+n,:))) , ibc_suscep0(:n,:n,ic,ifrq) )
            else
              suscep2(i+1:i+n,i+1:i+n) = ibc_suscep0(:n,:n,ic,ifrq)
            endif
            i = i + n
          enddo
        enddo
#   ifdef INV
        if(.not.trafo) then
          call symmetrize(suscep2,dim,nbasp,1)
        endif
#   endif
        ! Add contraction and double-counting term
        suscep2 = ( suscep2 + ibc_suscep1(:,:,ifrq) ) / nsym1
        ! Transform to Coulomb eigenvectors from right (case trafo)
        if(trafo) then
          suscep1           = matmul ( suscep2 , CTRAFO(:nbasp,:) )
        else
          suscep1           = 0
          suscep1(:,:nbasp) = suscep2
#   ifdef INV
          call symmetrize(suscep1,dim,nbasp,2)
#   endif
        endif
#   ifdef INV
        if(sum(abs(imag(suscep1)))>1d-10) Bug('Found nonzero imaginary part in ibc_suscep.')
#   endif
        ! Hermitianize
        suscep1         = ( suscep1 + conjg(transpose(suscep1)) ) / 2
        ! Add to suscep
        susceph(:,ifrq) = susceph(:,ifrq) + packmat(suscep1)
      enddo
      deallocate ( suscep2,ibc_suscep1 )
      if(allocated(ctrafo2)) deallocate ( ctrafo2 )        
      call timer_stop('IBC add_to_suscep',time=time1)
      write(6,'(A,F7.2,A)') 'done   ( Timing: ',time1,' )'
# else
      write(6,'(A'NoA) 'Add IBC to susceptibility... (no FIBC)'
      call timer_start('SUS IBC add_to_suscep')
      ! Add contraction and double-counting term
      ibc_suscep1 = ( ibc_suscep0 + ibc_suscep1 ) / nsym1
      ! Set matrix elements to zero that involve the constant function
      ibc_suscep1(1,:,:,:) = 0
      ibc_suscep1(:,1,:,:) = 0
      ! Hermitianize
      do ifrq = 1,nfrqh
        do ic = 1,ncent
          ibc_suscep1(:,:,ic,ifrq) = ( ibc_suscep1(:,:,ic,ifrq) + conjg(transpose(ibc_suscep1(:,:,ic,ifrq))) ) / 2
        enddo
      enddo
      ! Add to suscep
#   ifdef INV
      if(trafo) then
        allocate ( ctrafo2(nbasp,dim) )
        ctrafo2 = ctrafo
        call desymmetrize(ctrafo2,nbasp,dim,1)
      endif
#   endif
      if(allocated(suscep1)) deallocate(suscep1)
      allocate(suscep1(dim,dim))
      do ifrq = 1,nfrqh
        suscep1 = 0
        i       = 0
        ic      = 0
        do itype = 1,ntype
          n = sum ( [ ((2*l+1)*nindxm(l,itype),l=0,lcutm(itype)) ] )
          do ieq = 1,neq(itype)
            ic = ic + 1
            if(trafo) then
              suscep1 = suscep1 + matmul ( transpose(conjg(CTRAFO(i+1:i+n,:))) ,
     &                            matmul(ibc_suscep1(:n,:n,ic,ifrq),CTRAFO(i+1:i+n,:)) )
            else
              suscep1(i+1:i+n,i+1:i+n) = suscep1(i+1:i+n,i+1:i+n) + ibc_suscep1(:n,:n,ic,ifrq)
            endif
            i = i + n
          enddo
        enddo
        write(*,*) 'imaginary part = ', sum(abs(imag(suscep1)))
#   ifdef INV
        write(*,*) 'trafo=', trafo
        if(.not.trafo) then
          call symmetrize(suscep1,nbasm(ikpt),nbasp,3)
        endif
        if(sum(abs(imag(suscep1)))>1d-10) then
          write(*,*) 'imaginary part = ', sum(abs(imag(suscep1)))
          Bug('Found nonzero imaginary part in ibc_suscep.')
        endif
#   endif
        susceph(:,ifrq) = susceph(:,ifrq) + packmat(suscep1)
      enddo
      deallocate ( ibc_suscep1, suscep1 )
      if(allocated(ctrafo2)) deallocate ( ctrafo2 )
      call timer_stop('IBC add_to_suscep',time=time1)
      write(6,'(A,F7.2,A)') 'done   ( Timing: ',time1,' )'
# endif
# ifdef INV
#   undef CTRAFO
# endif
      end subroutine add_to_suscep

c     ---------------

      subroutine add_core
      implicit none
      complex_dp, allocatable :: suscep(:,:),wing1(:,:),wing2(:,:),ctrafo1(:,:)
      complex_dp              :: sus(dim,dim),wng(3,dim)
      real                    :: cputime
      integer                 :: iunit,ios
      integer                 :: itype,ieq,ifrq,nlm,lm,lm0,i,j,k MpiC(f1) MpiC(f2) MpiC(ff1(0:Msize-1)) MpiC(ff2(0:Msize-1))
      Rwrite(6,'(A'NoA) 'Add contribution of core susceptibility... '
c      suscepc = 0
c      if(ikpt==1) then
c        headc = 0
c        wingc = 0
c      endif
c# warning suscepc set to zero
      Mpi( MrangeDef(f1,f2,1,nfrqh+nfrqc) ; ff1=0;ff2=0 ; ff1(Mrank)=f1;ff2(Mrank)=f2 ; call Msum(ff1) ; call Msum(ff2) )
      call cpu_time(cputime)
      iunit = fopen('spex.core',form='unformatted',status='old')
      read(iunit) i,j,k ; if(any([i,j,k]/=[spin,ntype,nfrqh+nfrqc])) Error('Wrong parameters in spex.core.')
      lm0 = 0
      do itype = 1,ntype
        nlm = sum( [ ((2*l+1)*nindxm(l,itype),l=0,lcutm(itype)) ] )
        read(iunit) i   ; if(i==0) then ; lm0 = lm0 + nlm * neq(itype) ; cycle ; endif
        read(iunit) i,j ; if(any([i,j]/=[neq(itype),nlm])) Error('Wrong atom type parameters in spex.core.')
        allocate( suscep(nlm,nlm),wing1(nlm,3),wing2(3,nlm) )
        allocate( ctrafo1(nlm*neq(itype),dim) )
        ctrafo1 = ctrafo(lm0+1:lm0+nlm*neq(itype),:)
        Inv( call desymmetrize(ctrafo1,nlm*neq(itype),dim,3*itype+1) )
        do ifrq = 1,nfrqh+nfrqc ; Mpi2( if(ifrq<f1.or.ifrq>f2) then ; do ieq = 1,neq(itype) ; read(iunit) ; enddo ; cycle ; endif )
          lm = 0
          do ieq = 1,neq(itype)
            if(ikpt==1) then ; read(iunit,iostat=ios) suscep,wing1,wing2
            else             ; read(iunit,iostat=ios) suscep
            endif
            if(ios/=0) Error('Read error (suscep in spex.core).')
            sus = matmat( transpose(conjg(ctrafo1(lm+1:lm+nlm,:))) , matmat( suscep , ctrafo1(lm+1:lm+nlm,:) ) )
            if(ifrq<=nfrqh) then ; susceph(:,ifrq)         = susceph(:,ifrq)         + packmat(sus)
            else                 ; suscepc(:,:,ifrq-nfrqh) = suscepc(:,:,ifrq-nfrqh) + sus
            endif
            if(ikpt==1) then
              wng = transpose( matmul( transpose(conjg(ctrafo1(lm+1:lm+nlm,:))) , wing1 ) )
              if(ifrq<=nfrqh) then ; wingh(:,:,ifrq)         = wingh(:,:,ifrq)         + wng
              else                 ; wingc(:,:,1,ifrq-nfrqh) = wingc(:,:,1,ifrq-nfrqh) + wng
                                     wingc(:,:,2,ifrq-nfrqh) = wingc(:,:,2,ifrq-nfrqh) + matmul( wing2 , ctrafo1(lm+1:lm+nlm,:) )
              endif
            endif
            lm = lm + nlm
          enddo
        enddo
        lm0 = lm0 + nlm * neq(itype)
        deallocate( suscep,wing1,wing2,ctrafo1 )
      enddo
      if(ikpt==1) then
        do ifrq = 1,nfrqh+nfrqc ; Mpi( if(ifrq<f1.or.ifrq>f2) then ; read(iunit) ; cycle ; endif )
          read(iunit) head0
          if(ifrq<=nfrqh) then ; headh(:,ifrq)         = headh(:,ifrq)         + packmat(head0)
          else                 ; headc(:,:,ifrq-nfrqh) = headc(:,:,ifrq-nfrqh) + head0
          endif
        enddo
      endif
      call fclose(iunit)
      Rcall cpu_done(cputime)
      MpiR( write(6,'(A'NoA) 'MPI distribute... ' )
      NFENCE
      if(nfrqh>0) then
        OrangeDistr(susceph(:,  ff1(i):min(nfrqh,ff2(i))),i)
        if(ikpt==1) then
          MrangeDistr(wingh(:,:,ff1(i):min(nfrqh,ff2(i))),i)
          MrangeDistr(headh(:,  ff1(i):min(nfrqh,ff2(i))),i)
        endif
      endif
      if(nfrqc>0) then
        OrangeDistr(suscepc(:,:,  max(ff1(i)-nfrqh,1):ff2(i)-nfrqh),i)
        if(ikpt==1) then
          MrangeDistr(wingc(:,:,:,max(ff1(i)-nfrqh,1):ff2(i)-nfrqh),i)
          MrangeDistr(headc(:,:,  max(ff1(i)-nfrqh,1):ff2(i)-nfrqh),i)
        endif
      endif
      NFENCE
      MpiR( call cpu_done(cputime) )
      end subroutine add_core

c     ---------------

# ifdef Test_cprod
        subroutine test_cprod(iunit,cprod,rkp)
        use, intrinsic :: iso_fortran_env
        implicit none
        integer,     intent(in) :: iunit
        real_dp,     intent(in) :: rkp
        MCOMPLEX_dp, intent(in) :: cprod(bu1:bu2,bo1:bo2,nkpt1)
        rewind(iunit)
        write(iunit,'(5I5)') bo1,bo2,bu1,bu2,nkpt1
        write(iunit,*) cprod(29,1,1)/rkp,rkp
        i = 0
        do k = 1,nkpt1
          ikpt1 = kpt1(k)
          ikpt2 = kptsum(kpt1(k),ikpt)
          do i1 = bo1,bo2
            do i2 = bu1,bu2
              i  = i + 1
              d1 = deg(i1,ikpt1,ispin1) ; if(d1<i1) cycle
              d2 = deg(i2,ikpt2,ispin2) ; if(d2<i2) cycle
              write(iunit,'(5I5,F14.9)') i,i2,i1,d2-i2+1,d1-i1+1,sum(abs(cprod(i2:d2,i1:d1,k))**2)/rkp**2
            enddo
          enddo
        enddo
        end subroutine test_cprod
#   warning Test_cprod set
# endif

# ifdef MPI
#   define CB call Mcast
#   define CA call Mcastl
      subroutine Mbroadcast1
      use Mwrapper
      use, intrinsic :: iso_fortran_env
      implicit none
      CB(wghtthr); CB(multdiff); CB(disorder); CB(lhub)
      CA(hub) ; CB(ldum)
      end subroutine Mbroadcast1

      subroutine Mbroadcast2
      use Mwrapper
      use, intrinsic :: iso_fortran_env
      implicit none
      CB(nkpt1); CB(nsym1); CB(kpt1); CB(nkpts); CB(sym1); CB(kpt1p)
      end subroutine Mbroadcast2

      ! returns work packet for current rank (rank, msiz) -> no1,no2,nu1,nu2
      subroutine Mdistribute(n1,n2)
      use util, only: chr
      use, intrinsic :: iso_fortran_env
      implicit none
      integer, intent(in)  :: n1,n2
      integer, allocatable :: prim(:)
      integer              :: i,j,n,m1,m2,k1,k2,min1,d1,d2
      real_dp              :: q1,q2,dv,sp
      real_dp, parameter   :: aspect = 10d0 ! aspect ratio of optimal block
c# warning aspect changed
      allocate ( prim(ceiling(log(1d0*msiz)/log(2d0))) )
      call primfac(prim,n,msiz)
      if(product(prim(:n))/=msiz) Bug('Prime factorization failed.')
      min1 = huge(0)
      do i = 0,2**n-1
        k1 = product(prim(:n),[(iand(i,2**j)/=0,j=0,n-1)]) ! # processes for occupied states
        k2 = msiz / k1                                     ! # processes for unoccupied states
        if(k1*k2/=msiz) Error('Bug.')
        if(k1>n1) cycle
        if(k2>n2) cycle
        m1 = ceiling ( 1d0 * n1 / k1 - 1d-12 )
        m2 = ceiling ( 1d0 * n2 / k2 - 1d-12 )
        q1 = sqrt(m1*m2/aspect)
        q2 = aspect*q1
        dv = (abs(m1-q1)*q2+abs(m2-q2)*q1/10)*msiz ! deviation from "optimal" block (q1,q2) with correct aspect ratio
        sp = abs(m1*m2*msiz-n1*n2)                 ! spill of block area
c        if(Mrank==0) then
c          write(*,*) k1,k2,m1,m2
c          write(*,*) q1,q2,dv/5,sp
c        endif
        if(dv/5+sp<min1) then ! distribution criterion (minimize dv/5+sp)
          min1 = dv/5+sp
          d1   = k1
          d2   = k2
        endif
      enddo
      if(d1==0) Error('Could not determine MPI work packets. (Decrease processes!)')

      m1 = ceiling ( 1d0 * n1 / d1 - 1d-12 )
      m2 = ceiling ( 1d0 * n2 / d2 - 1d-12 )
      Rif(first_term) then
        write(6,'(A,10X,A)') 'MPI work packages:',trim(chr(d1))//'*'//trim(chr(d2))//
     &                                      ' ('//trim(chr(m1))//'*'//trim(chr(m2))//')'
      endif

      q1 = 0
      m1 = 1
      i  = 0
      do k1 = 1,d1
        q1 = q1 + 1d0 * n1 / d1
        q2 = 0
        m2 = 1
        do k2 = 1,d2
          q2 = q2 + 1d0 * n2 / d2
          if(i==rank) then
            no1 = m1
            no2 = nint(q1)
            nu1 = m2
            nu2 = nint(q2)
            return
          endif
          i  = i + 1
          m2 = nint(q2)+1
        enddo
        m1 = nint(q1)+1
      enddo
      Bug('Work distribution to processes failed.')
      end subroutine Mdistribute

      subroutine Mreduce
      use Mwrapper
      use, intrinsic :: iso_fortran_env
      implicit none
      integer                  :: i
      real                     :: cputime
      call mpi_barrier(Mcomm,i)
      Rwrite(6,'(/A'NoA) 'MPI reduce... '
      if(Nrank==0) then
        call cpu_time(cputime)
        do i = 1,nfrqh
          call Msum(susceph(:,i),comm=Ocomm)
        enddo
        do i = 1,nfrqc
          call Msum(suscepc(:,:,i),comm=Ocomm)
        enddo
      endif
      if(ikpt==1) then
        do i = 1,nfrqh
          call Msum(headh(:,i))
          call Msum(wingh(:,:,i))
        enddo
        do i = 1,nfrqc
          call Msum(headc(:,:,i))
          call Msum(wingc(:,:,1,i))
          call Msum(wingc(:,:,2,i))
        enddo
      endif
      Rcall cpu_done(cputime)
      end subroutine Mreduce
#   undef CA
#   undef CB
# endif

c     ---------------

      end

c     ---------------

c
c     Returns symmetry operations in sym1(:nsym1) that leave kvec(:3) (in internal coordinates) invariant.
c     (Dimension of sym1 must be large enough, nsym always sufficient.)      
      subroutine getsym1(sym1,nsym1,kvec)
      use global, only: sym,nsym,modulo1r
      use, intrinsic :: iso_fortran_env
      implicit none
      integer, intent(out) :: nsym1,sym1(*)
      real_dp, intent(in)  :: kvec(3)
      real_dp              :: kvec1(3)
      integer              :: isym,i
      i = 0
      do isym = 1,nsym
        kvec1 = matmul(sym(isym)%rrot,kvec)        
        if(all( abs(modulo1r(kvec1-kvec)) < 1d-10 )) then
          i       = i + 1
          sym1(i) = isym
        endif
      enddo
      nsym1 = i
      end

c
c     Look for symmetry operations that leave kpt(:,ikptq) invariant
c     Define corresponding reduced k-point set EIBZ (if ikpt/=0, only k points that are equivalent to ikpt)
c
c     sym1(:nsym1)  : symmetry operations
c     kpt1(:nkpt1)  : indices of irreducible k-point set wrt sym1 (EIBZ)
c     nkpts(:nkpt1) : multiplicity of k points
      subroutine getkpt1(kpt1,nkpt1,nkpts,sym1,nsym1,ikptq,ikpt,writeout)
      use global
      logical, intent(in)  :: writeout
      integer, intent(in)  :: ikptq,ikpt
      integer, intent(out) :: kpt1(nkpt),nkpt1,nkpts(nkpt),sym1(nsym),nsym1
      integer, allocatable :: pnt(:),kpt2(:)
      integer              :: i,i1,k,isym
      logical              :: check(nkpt)
      integer              :: kptsum      
# ifdef switch_off_symmetry
      call getkpt1_fullBZ(kpt1,nkpt1,nkpts,sym1,nsym1,ikpt)
      Warn('Symmetry is switched off.')
      return
# endif
      if(.not.use_sym) then
        call getkpt1_fullBZ(kpt1,nkpt1,nkpts,sym1,nsym1,ikpt)
        return
      endif
      i = 0
      do isym = 1,nsym
        if(kptsym(ikptq,isym)==ikptq) then
          i       = i + 1
          sym1(i) = isym
        endif
      enddo
      nsym1 = i ; if(writeout) write(6,'(A,I4)') 'Number of sym. op.:      ',nsym1
      check = .false.
      nkpts = 0
      i1    = 0
      do i = 1,nkpt
        if(ikpt/=0) then
          if(kptp(i)/=kptp(ikpt)) cycle
        endif
        if(check(i)) cycle
        i1 = i1 + 1 ; kpt1(i1) = i
        do isym = 1,nsym1
          k = kptsym(i,sym1(isym))
          if(.not.check(k)) then
            nkpts(i1) = nkpts(i1) + 1
            check(k)  = .true.
          endif
        enddo
      enddo
      nkpt1 = i1 ; if(writeout) write(6,'(A,I5)') 'Number of irr. k-points:',nkpt1
      allocate ( pnt(nkpt1),kpt2(nkpt1) ) ! order kpt1 such that kptp(kpt1+ikptq) is size-ordered (i.e., k points with the same kptp(ikpt2) are consecutive).
      do i = 1,nkpt1
        kpt2(i) = kptp(kptsum(kpt1(i),ikptq))
      enddo
      call rorderp(pnt,kpt2*1d0,nkpt1)
      kpt1(:nkpt1)  = kpt1(pnt)
      nkpts(:nkpt1) = nkpts(pnt)
      deallocate ( pnt,kpt2 )
      end

c     Replace EIBZ by full BZ.
      subroutine getkpt1_fullBZ(kpt1,nkpt1,nkpts,sym1,nsym1,ikpt)
      use global
      integer, intent(in)  :: ikpt
      integer, intent(out) :: kpt1(nkpt),nkpt1,nkpts(nkpt),sym1(nsym),nsym1
      integer              :: i,j
      if(ikpt/=0) then
        j = 0
        do i = 1,nkpt
          if(kptp(i)==kptp(ikpt)) then
            j       = j + 1
            kpt1(j) = i
          endif
        enddo
        nkpt1 = j
      else
        nkpt1 = nkpt
        kpt1  = [(i,i=1,nkpt)]
      endif
      nkpts   = 1
      nsym1   = 1
      sym1(1) = 1
      end

      ! Construct kpt1p (parent k points of EIBZ); kpt1p(k) is zero if k point does not have a parent in kpt1
      ! lsymkpt1 = true : also write corresponding symops to symkpt1
      subroutine getkpt1p(kpt1p,kpt1,nkpt1,sym1,nsym1,symkpt1,lsymkpt1)
      use global, only: kptsym,nkpt
      implicit none
      integer, intent(in)  :: nkpt1,nsym1
      integer, intent(out) :: kpt1p(nkpt),symkpt1(*)
      integer, intent(in)  :: kpt1(nkpt1),sym1(nsym1)
      logical, intent(in)  :: lsymkpt1
      integer              :: i,isym,ikpt
      kpt1p = 0
      if(lsymkpt1) symkpt1(:nkpt) = 0
      do i = 1,nkpt1        
        do isym = 1,nsym1
          ikpt        = kptsym(kpt1(i),sym1(isym))
          kpt1p(ikpt) = i
          if(lsymkpt1) then
            if(symkpt1(ikpt)==0) symkpt1(ikpt) = sym1(isym)
          endif
        enddo
      enddo
      end

c     ---------------

      subroutine dipole
      use global
      use key
      implicit none
      integer, allocatable :: bnd(:)
      logical              :: occup(maxband)
      MCOMPLEX_dp          :: momentum(3)
      real_dp              :: rdum,ediff
      integer              :: ikpt,ikpt0,ispin
      integer              :: i,ibando,ibandu
      Load( Error('DIPOLE not implemented for -DLOAD.') )
      allocate(bnd(2))
      call getkey(inp,'DIPOLE', bnd, status=i)
      if     (i==0) then ; Bug('DIPOLE not found.')
      else if(i==1) then ; bnd = [ 1,huge(0) ]
      else               ; if(any(bnd<1)) Error('Arguments to DIPOLE must be positive integers.')
      endif
      write(6,'(//A)') 'Dipole matrix elements:'
      do ispin = 1,nspin1
        if     (nspin1==1) then ; write(6,*)
        else if(ispin ==1) then ; write(6,'(/A)') 'Spin up'
        else if(ispin ==2) then ; write(6,'(/A)') 'Spin down'
        endif
        write(6,'(A)') 'K point  Transition  Matrix element (squared)  (transition energy)'
        i = 0
        do ikpt = 1,nkpti+1
          ikpt0 = ikpt
          if(ikpt==nkpti+1) then
            if(.not.lkptadd) exit
            ikpt0 = nkpt+1
          endif
          occup = ene(:,ikpt0,ispin)<=efermi
          do ibando = bnd(1),min(nband(ikpt0,ispin),bnd(2))
            if(occup(ibando)) then
              do ibandu = bnd(1),min(nband(ikpt0,ispin),bnd(2))
                if(.not.occup(ibandu)) then
                  call momentum_matrix(momentum,[ikpt0],1,ispin,ispin,ibando,ibando,ibandu,ibandu,.false. MpiC(.false.) )
                  ediff = ene(ibandu,ikpt,ispin) - ene(ibando,ikpt,ispin)
                  rdum  = 0
                  if(ediff>1d-8) then
                    rdum = sum(abs(momentum)**2) / ediff**2
                  endif
                  if(rdum>=.5d-6) then
                    write(6,'(I7,I6,'' ->'',I3,F16.6,12X,''('',F8.3,'' eV )'')')
     &                ikpt0,ibando,ibandu,rdum,hartree*(ene(ibandu,ikpt,ispin)-ene(ibando,ikpt,ispin))
                    i = i + 1
                  endif
                endif
              enddo
            endif
          enddo
        enddo
        if(i==0) write(6,'(A/,A)') '... no values ...','Either all occupied or all unoccupied or all absolute values below 5e-6.'
      enddo
      end

c     ---------------

# ifdef MPI

      function scaling_suscep(dim,ndim)
      use global, only: neq,nindx,lcut,ntype,ncent
      use, intrinsic :: iso_fortran_env
      implicit none
      real_dp             :: scaling_suscep
      integer, intent(in) :: ndim,dim(ndim)
      integer             :: itype,l
      real_dp             :: overhead
      ! prefactor 2/3 experimentally determined for wavefproducts2_mt; assuming that this routine makes up 25% of the calculation time gives a prefactor 1/6
      overhead = 1d0/6 * sum( [ ((neq(itype)*(2*l+1)*nindx(l,itype),l=0,lcut(itype)),itype=1,ntype) ] ) / ncent            
      if(ndim==2) then
        scaling_suscep = dim(2) * ( overhead + dim(3) )
      else if(ndim==3) then
        scaling_suscep = dim(1) * dim(2) * ( overhead + dim(3) )
      else
        Bug('Array dimension neither two nor three.')
      endif
      end
# endif

c     ---------------

# ifdef ACCESS_CONTIG

c
c     Transforms matrix indices ( 11, 12, 21, 22, 13, 31, 23, 32, 33, ... )
c                          into ( 11, 21, 31, 41, ... )
      subroutine get_index_suscepc(index,args)
      implicit none
      integer, intent(inout) :: index
      integer, intent(in)    :: args(1)
      integer                :: dim,i,j
      dim = args(1)
      j   = 0
      do while(j**2<index)
        j = j + 1
      enddo
      index = index - (j-1)**2 + 1 ! = 2i (2i+1) for ij (ji)
      i     = index / 2 ; if(mod(index,2)==1) then ; index = i ; i = j ; j = index ; endif ! transpose ij->ji     
      index = (j-1)*dim + i
      end

# endif

c     ---------------

# undef OUTP
# undef MpiR_WRITE
