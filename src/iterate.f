c--------------------------------------------------------------------------------
c Copyright (c) 2021 Peter Gr�nberg Institut, Forschungszentrum J�lich, Germany
c This file is part of SPEX and available as free software under the conditions
c of the MIT license as expressed in the LICENSE file in more detail.
c--------------------------------------------------------------------------------

# ifndef INV
#   define INV
#   define INV_SOC
#   include "iterate.f"
#   undef INV_SOC
#   undef INV
# endif

# include "cppmacro.h"
# include "restype.h"

c 
c Solve KS equations (one-shot calculation).
c
c it_mode = 1 : (NR) non-relativistic
c it_mode = 2 : (SR) scalar-relativistic
c it_mode = 3 : (FR) full-relativistic
c it_stop = 0 : continue after calculation
c it_stop = 1 : stop after calculation
c it_stop = 2 : band-structure calculation (and stop)
c

c MPI: Define kind of ordering of output (MPI_order 0, 1, 2 for ii==0, ii>0, ii<0, respectively, see Mwrite_n)
# define MPI_order 2

# if MPI_order < 0 || MPI_order > 2
#   error unknown MPI_order
# endif

# ifdef INV_SOC
      subroutine iterate_inv_soc(it_mode,it_stop) ! special version for INV and l_soc: constructs real symmetric Hamiltonian in 1st variation
# else
      subroutine iterate(it_mode,it_stop)
# endif      

# ifndef TDDFT

      use global
      use wrapper
      use key
      use util
      use readwrite
      use Hwrapper
      use hdf5
      use, intrinsic :: iso_fortran_env
      Mpi ( use Mwrapper )
      Mpi2( use, intrinsic :: iso_c_binding )
      implicit none
      integer,       intent(in)  :: it_mode,it_stop
      complex_dp,    allocatable :: hpwlo(:,:,:),opwlo(:,:,:),hlolo(:,:,:,:),carr(:,:)
      complex_dp,    allocatable :: cmt1(:,:),cmt2(:,:),cmt3(:,:),cmt4(:,:),hlp(:),hlp2(:,:)
      complex_dp,    allocatable :: cmt0(:,:,:)
      MCOMPLEX_dp,   allocatable :: cpw0(:,:),mcmt(:,:)
      MCOMPLEX_dp,   pointer_cnt :: hamiltonian(:),overlap(:),evec(:,:)
      complex_dp,    pointer_cnt :: hamiltonian_soc(:,:),evec_soc(:,:)
      real_dp,       allocatable :: eig(:)
      complex_dp                 :: ylm((maxlcut+1)**2),ylm2((maxlcut+1)**2),ylm1(0:maxlcut)
      complex_dp,    allocatable :: harm(:,:),phas(:,:),phas0(:),phas1(:),phas2(:)
# ifdef CHECK_ifort18bug
#   warning hmtpw,hmtpw0 made allocatable
      complex_dp, allocatable, target :: hmtpw(:,:,:),hmtpw0(:,:)
# else
      complex_dp, target              :: hmtpw(2,(maxlcut+1)**2,maxval(neq)),hmtpw0(2,(maxlcut+1)**2)
# endif
      complex_dp,    pointer_cnt :: hmtpw_p(:,:),hmt(:,:,:,:,:)
      complex_dp                 :: dwgn1(-maxlcut:maxlcut,-maxlcut:maxlcut,0:maxlcut)
      complex_dp                 :: cdum,cexp,cexp1,cexp2(maxval(neq)),kphas(ncent),imgl
      MCOMPLEX_dp                :: mcdum,mexp
      real_dp                    :: pi4vol
      real_dp                    :: omt(maxindx,maxindx,0:maxlcut,ntype),omtpw(2,0:maxlcut)
      real_dp                    :: omt0(maxindx,maxindx,0:maxlcut,ntype)
      real_dp                    :: apwmat(2,2,0:maxlcut,ntype)
      real_dp                    :: integ(maxindx,maxindx,maxlh,0:maxlcut,0:maxlcut),integr(maxindx,maxindx),trafo(maxindx,maxindx)
      real_dp                    :: dbas1(maxgrid,maxindx),dbas2(maxgrid,maxindx),dvmt(maxgrid),dvmt_so(maxgrid)
      real_dp                    :: ulo,dulo,a(2),rdum,rdum1,x,rad
      real_dp                    :: ubas1,ubas2,dubas1,dubas2
      real_dp                    :: kvec(3),kg1(3),kg2(3),rot(3,3),kg1n,kg2n,jl(0:maxlcut+1),dj,gaunt1,leg(0:maxlcut)
      real_dp,       allocatable :: kgn(:),apw(:,:,:,:),gtrafo(:,:),avg(:),qpath(:)
      real_dp,       allocatable :: proj(:,:,:),project(:,:,:,:)
      integer,       allocatable :: pnt(:,:)
      integer,       parameter   :: packet = 10
      integer                    :: neq1
      integer                    :: g(3),nband00,nband1,nband_soc
      integer                    :: nbas1,nlotot,maxlmnlo,ngrid
      integer                    :: isym,ilow,ifac
      integer                    :: ispin,ikpt
      integer                    :: itype,ieq,ic,ilh,jlh
      integer                    :: itype1,ieq1,ic1,itype2,ieq2,ic2
      integer                    :: i,i1,j,k,k1,l,l1,l2,m,m1,m2,n,n1,n2,lm,lm1,lm2,lmn,lmn1,lmn2
      integer                    :: ix,iy,iz
      integer                    :: iunit
      logical                    :: found,ldum,linv
      real                       :: cputime,cputime0
      real_dp                    :: intgrf,gaunt
      character(:),  allocatable :: lines(:)
# ifdef INV_SOC
      complex_dp                 :: stepfunction
# else
      MCOMPLEX_dp                :: stepfunction
# endif
# ifdef INV
      logical,       allocatable :: invp(:)
# endif      
# ifdef MPI
#   if MPI_order == 2
      integer                    :: ikpt0
#   endif
      type(c_ptr)                :: ptr
      integer                    :: win_hamiltonian,win_overlap,win_evec,win_hamiltonian_soc,win_evec_soc,win_hmt
      integer                    :: Merr,Mcolor,Munit,Mstat(mpi_status_size),ios,Mwrt_n,Mcoms(1:Msize-1)
      real_dp                    :: xmem
      logical                    :: pcycle_skip
      character(256)             :: line
      complex_dp,    pointer_cnt :: Ninit_cmt(:),Ninit_hamiltonian_soc(:),Ninit_hmt(:)
#   ifdef INV_SOC
      complex_dp,    pointer_cnt :: Ninit_cpw(:),Ninit_hamiltonian(:),Ninit_overlap(:)      
#   else
      MCOMPLEX_dp,   pointer_cnt :: Ninit_cpw(:),Ninit_hamiltonian(:),Ninit_overlap(:)
#   endif
# endif
# include "interface/getqlist.inc"
# include "interface/band_info.inc"

c#   define Mdo(arg) Mpi( if(Mrank==0) then ; ) arg Mpi( ; endif )
#   define Mdo(arg)

# ifdef LOAD
      Error('ITERATE not implemented for -DLOAD.')
# else

      Rwrite(6,'(//A)') '### subroutine iterate ###'

      if     (it_mode==2.and.     l_soc) then ; Bug('ITERATE SR but l_soc=.true.')
      else if(it_mode==3.and..not.l_soc) then ; Bug('ITERATE FR but l_soc=.false.')
      endif

      call cpu_time(cputime0)

      pi4vol = 4*pi / sqrt(vol)

c
c     Prepare band structure calculation
      if(it_stop==2) then
        deallocate(kpt)
        Rbegin
        call getkey(inp,'ITERATE', lines, default=[' ']) ; if(lines(size(lines))=='BANDS') lines(size(lines)) = ' '
        call getqlist(kpt,nkpt,20,lines(size(lines)))
        deallocate(lines)
        write(6,'(/A)') 'Band structure calculation with '//Chr(nkpt)//' k points.'
        if(bandinfo) then
          call getkey(inp,'PROJECT',lines,status=i)
          if(i==0)      then ; Bug('PROJECT not found.')
          else if(i==1) then ; allocate(character(80) :: lines(ntype)) ; lines = '[s,p,d,f,g]'
          endif
        endif
        Rend
# ifdef MPI        
        call Mcast(nkpt)
        call Mcastl(kpt)
        if(bandinfo) then
          ifR n = size(lines)
          call Mcast(n) ; MpiR( l = len(lines(1)) )
          call Mcast(l) ; MnoR( allocate(character(l) :: lines(n)) )
          do i = 1,n
            call Mcast(lines(i))
          enddo
        endif
# endif
        nkpti = nkpt
        nkpt1 = nkpt
        nkpt2 = nkpt
      endif

c
c     Determine number of local orbitals
      i = 0
      j = 0
      do itype = 1,ntype
        if(any(nindx(:lcut(itype),itype)<2)) Error('missing u and/or udot.')
        i = i + sum( [ (neq(itype) * (nindx(l,itype)-2) * (2*l+1), l=0,lcut(itype)) ] )
        j = max(j,           sum( [ ((nindx(l,itype)-2) * (2*l+1), l=0,lcut(itype)) ] ) )
      enddo
      nlotot   = i
      maxlmnlo = j

      symm_maxindx = max(symm_maxindx,3) ! the maximal index used in this routine is indx=2(=3-1) (see trafo.f)

      allocate ( nband(nkpt2,nspin2),ngpt(nkpt2) ) ; nband = 0

c
c     Construct list of G vectors for APWs
      Mdo( call cpu_time(cputime) )
      ngpt = 0
 1    do ikpt = Mrange1(nkpt1)
        if(ikpt<=nkpti) then ; kvec = kpt(:,ikpt)
        else                 ; kvec = kpt(:,ikpt+nkpt-nkpti)
        endif
        i     = 0
        n     = 0
        found = .true.
        do while(found)
          found = .false.
          do ix = -n,n
            do iy = -(n-abs(ix)),n-abs(ix)
              iz   = n - abs(ix) - abs(iy)
 2            g    = [ ix,iy,iz ] - nint(kvec)
              kg1  = matmul(rlat,kvec+g)
              kg1n = sum(kg1**2)
              if(kg1n<=gcut**2) then
                i     = i + 1
                found = .true.
                if(allocated(gpt1)) gpt1(:,i,ikpt) = g
              endif
              if(iz>0) then
                iz = -iz
                goto 2
              endif
            enddo
          enddo
          n = n + 1
        enddo
        ngpt(ikpt) = i
      enddo
      if(.not.allocated(gpt1)) then
        Mpi( call Msum(ngpt) )
        maxgpt = maxval(ngpt) ; allocate ( gpt1(3,maxgpt,nkpt2) )        
        gpt1   = 0
        goto 1
      endif
      Mpi( call Msum(gpt1) )
      Mdo( Rwrite(6,'(A'NoA) 'prep0 ' ; Rcall cpu_done(cputime) )

      g = [ ( maxval(gpt1(i,:,:)) - minval(gpt1(i,:,:)) + 1 , i = 1,3 ) ]
      Allocate_ ( cstep, ( -g(1):g(1),-g(2):g(2),-g(3):g(3) ) )
      cstep = 0
      l     = -1
      do k = -g(3),g(3)
      do j = -g(2),g(2)
      do i = -g(1),g(1) ; McycleP(l)
        cstep(i,j,k) = stepfunction( [i,j,k] )
      enddo
      enddo
      enddo
      Mpi( call Msum(cstep) )

      Mdo( Rwrite(6,'(A'NoA) 'cstep ' ; Rcall cpu_done(cputime) )

c
c     Allocate arrays
      Rif(maxene/=-huge(0d0)) Warn('Max. energy defined after NBAND. Cannot determine number of bands. Large memory demand.')
      if(nband0==-1) then ; nband00 = maxgpt + nlotot ;  if(l_soc) nband00 = 2 * nband00
      else                ; nband00 = nband0          ; Rif(l_soc.and.mod(nband00,2)/=0) Error('Use even number of bands for SOC.')
      endif
      allocate ( ene(nband00,nkpt2,nspin) ) ; maxeband = huge(0)
      Mpi( call Mcast(nband00) )
      ene   = 0
      nband = 0
      if(it_stop==0) then
        if(storeibz) then
          Nallocate0 ( cpw,(S_ maxgpt,         nband00,nkpt1,nspin2 S_) )
          Nallocate0 ( cmt,(S_ maxlmindx,ncent,nband00,nkpt1,nspin2 S_) )
        else
          Nallocate0 ( cpw,(S_ maxgpt,         nband00,nkpt2,nspin2 S_) )
          Nallocate0 ( cmt,(S_ maxlmindx,ncent,nband00,nkpt2,nspin2 S_) )
        endif
      endif

      Mdo( Rwrite(6,'(A'NoA) 'alloc ' ; Rcall cpu_done(cputime) )

c
c     Loop over spins
c
      do ispin = 1,nspin1

      Rif(nspin1==2) then
        if(ispin==1) write(6,'(/A)') 'SPIN UP'
        if(ispin==2) write(6,'(/A)') 'SPIN DOWN'
      endif

c
c     Precalculation (k independent)

      Mdo ( call cpu_time(cputime) )
      
c     (A) Calculate matrix for construction of APWs and define lo radial functions with the properties ulo(R)=ulo'(R)=0 (bas1/2 is redefined!)
      do itype = 1,ntype
        do l = 0,lcut(itype)
          apwmat(1,1,l,itype) =  ubas(1,l,itype,ispin)
          apwmat(1,2,l,itype) =  ubas(2,l,itype,ispin)
          apwmat(2,1,l,itype) = dubas(1,l,itype,ispin)
          apwmat(2,2,l,itype) = dubas(2,l,itype,ispin)
          call inverse(apwmat(:,:,l,itype))
          do n1 = 1,nindx(l,itype)
            do n2 = 1,nindx(l,itype)
              omt0(n1,n2,l,itype) = intgrf( bas1(:,n1,l,itype,ispin)*bas1(:,n2,l,itype,ispin)
     &                                    + bas2(:,n1,l,itype,ispin)*bas2(:,n2,l,itype,ispin) , itype )
            enddo
          enddo
          do n = 3,nindx(l,itype)
            ulo                     =  ubas(n,l,itype,ispin)
            dulo                    = dubas(n,l,itype,ispin)
            a                       = matmul(apwmat(:,:,l,itype),-[ulo,dulo])
            bas1(:,n,l,itype,ispin) = a(1) * bas1(:,1,l,itype,ispin) +
     &                                a(2) * bas1(:,2,l,itype,ispin) +
     &                                       bas1(:,n,l,itype,ispin)
            bas2(:,n,l,itype,ispin) = a(1) * bas2(:,1,l,itype,ispin) +
     &                                a(2) * bas2(:,2,l,itype,ispin) +
     &                                       bas2(:,n,l,itype,ispin)
          enddo
          do n1 = 1,nindx(l,itype)
            do n2 = 1,nindx(l,itype)
              omt(n1,n2,l,itype) = intgrf( bas1(:,n1,l,itype,ispin)*bas1(:,n2,l,itype,ispin) +
     &                                     bas2(:,n1,l,itype,ispin)*bas2(:,n2,l,itype,ispin) , itype )
            enddo
          enddo
        enddo
      enddo
      Mdo( write(6,'(A'NoA) 'prep ' ; call cpu_done(cputime) )

c     (B) MT part of kinetic energy (->hmt)
      Nallocate0( hmt , (S_ maxindx,maxindx,(maxlcut+1)**2,(maxlcut+1)**2,ncent S_) )
      Obegin
      if(it_mode>=2) then
        ! We calculate the scalar-relativistic kinetic energy including spherical part of effective potential energy
        do itype = 1,ntype
          lm = 0
          do l = 0,lcut(itype)
            do n1 = 1,nindx(l,itype)
              do n2 = 1,nindx(l,itype)
                rdum1 = 0
                if(n1==2.and.n2/=2) rdum1 = omt0(1,n2,l,itype)
                if(n2==2.and.n1/=2) rdum1 = omt0(n1,1,l,itype)
                rdum          = omt0(n1,n2,l,itype)
                integr(n1,n2) = ( ( ebas(n1,l,itype,ispin) + ebas(n2,l,itype,ispin) ) * rdum + rdum1 ) / 2
              enddo
            enddo
            trafo = identity(maxindx)
            do n = 3,nindx(l,itype)
              ulo         =  ubas(n,l,itype,ispin)
              dulo        = dubas(n,l,itype,ispin)
              trafo(:2,n) = matmul(apwmat(:,:,l,itype),-[ulo,dulo])
            enddo
            n             = nindx(l,itype)
            integr(:n,:n) = matmul ( transpose(trafo(:n,:n)) , matmul(integr(:n,:n),trafo(:n,:n)) )
            do m = -l,l
              lm                     = lm + 1
              hmt(:n,:n,lm,lm,itype) = integr(:n,:n)
            enddo
          enddo
        enddo
      else
        ! we calculate the non-relativistic kinetic energy
        ! - 1/4 INT (         u_l(r) Y^*_lm(r) laplace u_l'(r) Y_l'm'(r)
        !           + laplace u_l(r) Y^*_lm(r)         u_l'(r) Y_l'm'(r) ) d^3 r
        do itype = 1,ntype
          ngrid = grid(itype)%number
          rad   = grid(itype)%radius
          lm    = 0
          do l = 0,lcut(itype)
            do n = 1,nindx(l,itype)
              ! construct first derivative u'(r)
              call derivative(dbas1(:,n),bas1(:,n,l,itype,ispin),itype)
              call derivative(dbas2(:,n),bas2(:,n,l,itype,ispin),itype)
              dbas1(:,n) = dbas1(:,n) - bas1(:,n,l,itype,ispin)/rgrid(:,itype)
              dbas2(:,n) = dbas2(:,n) - bas2(:,n,l,itype,ispin)/rgrid(:,itype)
              if(n<=2) then ; dbas1(ngrid,n) = rad * dubas(n,l,itype,ispin) ! replace derivative at last point by the more accurate dubas value
              else          ; dbas1(ngrid,n) = 0                            !
              endif
            enddo
            do n1 = 1,nindx(l,itype)
              if(n1<=2) then ; ubas1 = ubas(n1,l,itype,ispin) ; dubas1 = dubas(n1,l,itype,ispin)
              else           ; ubas1 = 0                      ; dubas1 = 0
              endif
              do n2 = 1,n1
                if(n2<=2) then ; ubas2 = ubas(n2,l,itype,ispin) ; dubas2 = dubas(n2,l,itype,ispin)
                else           ; ubas2 = 0                      ; dubas2 = 0
                endif
                integr(n1,n2) = rad**2 * ( ubas1*dubas2 + dubas1*ubas2 ) / 2 -
     &                          intgrf( dbas1(:,n1)*dbas1(:,n2) + dbas2(:,n1)*dbas2(:,n2) ,itype) -
     &                l*(l+1) * intgrf( ( bas1(:,n1,l,itype,ispin)*bas1(:,n2,l,itype,ispin) +
     &                                    bas2(:,n1,l,itype,ispin)*bas2(:,n2,l,itype,ispin) ) / rgrid(:,itype)**2 ,itype)
                integr(n2,n1) = integr(n1,n2)
              enddo
            enddo
            n = nindx(l,itype)
            do m = -l,l
              lm                     =   lm + 1
              hmt(:n,:n,lm,lm,itype) = - integr(:n,:n) / 2
            enddo
          enddo
        enddo
      endif      
      Mdo( Rwrite(6,'(A'NoA) 'prep2 ' ; Rcall cpu_done(cputime) )
      Oend

c     (C) MT part of total potential (->hmt)
      do itype = 1,ntype
        Mpi( integ = 0 ; i = 0 )
        ! precalculate radial integrals
        do l2 = 0,lcut(itype)
          do l1 = 0,l2
            do ilh = 1,nlh(itype)
              l = llh(ilh,itype)
              if(it_mode>=2.and.ilh==1) cycle ! leave out spherical part
              if(mod(l1+l2+l,2)/=0)     cycle ; McycleP(i)
              do n2 = 1,nindx(l2,itype)
                do n1 = 1,nindx(l1,itype)
                  integ(n1,n2,ilh,l1,l2) = intgrf( vmt(:,ilh,itype,ispin) * (
     &                                             bas1(:,n1,l1,itype,ispin)*bas1(:,n2,l2,itype,ispin) +
     &                                             bas2(:,n1,l1,itype,ispin)*bas2(:,n2,l2,itype,ispin) ) , itype )                  
                  if(ilh==1) integ(n1,n2,1,l1,l2)   = integ(n1,n2,1,l1,l2) * sqrt(4*pi) ! multiply with 1/Y_00 such that V_00=v_00*Y_00
                  if(l1/=l2) integ(n2,n1,ilh,l2,l1) = integ(n1,n2,ilh,l1,l2)            ! offdiagonal element
                enddo                
              enddo
            enddo
          enddo
        enddo
        Mpi( call Msum(integ) )
        Obegin
        lm2 = 0
        do l2 = 0,lcut(itype)
          n2 = nindx(l2,itype)
          do m2 = -l2,l2
            lm2 = lm2 + 1
            lm1 = 0
            do l1 = 0,l2
              n1 = nindx(l1,itype)
              do m1 = -l1,l1
                lm1 = lm1 + 1
                do ilh = 1,nlh(itype)
                  if(it_mode>=2.and.ilh==1)   cycle ! leave out spherical part
                  l = llh(ilh,itype)
                  if(mod(l1+l2+l,2)/=0)       cycle ! Gaunt condition                  
                  if(l<abs(l1-l2).or.l>l1+l2) cycle ! Gaunt condition
                  do jlh = 1,nmlh(ilh,itype)
                    m = mlh(jlh,ilh,itype)
                    if(m/=m1-m2)              cycle ! Gaunt condition
                    gaunt1                     = gaunt(l,l1,l2,m,m1,m2)
                    hmt(:n1,:n2,lm1,lm2,itype) = hmt(:n1,:n2,lm1,lm2,itype) + integ(:n1,:n2,ilh,l1,l2) * gaunt1 * clh(jlh,ilh,itype)
                  enddo
                enddo
                if(lm1<lm2) then ! offdiagonal elements
                  n1                         = nindx(l1,itype)
                  n2                         = nindx(l2,itype)
                  hmt(:n2,:n1,lm2,lm1,itype) = conjg(transpose( hmt(:n1,:n2,lm1,lm2,itype) ) )
                endif
              enddo
            enddo
          enddo
        enddo
        Oend
      enddo
      Mdo( Rwrite(6,'(A'NoA) 'prep3 ' ; Rcall cpu_done(cputime) )

c     (C) make hmt atom-dependent and rotate symmetry-equivalent parts to the global coordinate system
      Obegin
      ic = ncent
      do itype = ntype,1,-1
        l = lcut(itype)
        do ieq = neq(itype),1,-1
          isym = symcent(ic) ; if(isym<0) Bug('symcent < 0')
          if(isym/=1) then
            rot = matmul(lat,matmul(sym(abs(isym))%rot,transpose(rlat)))/(2*pi) ! this is according to Fleur's "symcent"
            if(isym<0) rot = -rot                                               ! in Spex's symcent, isym is always positive
            call dwigner(dwgn1(-l:l,-l:l,0:l),rot,l)
            do l2 = 0,lcut(itype)
              lm2 = l2**2 + 1
              do l1 = 0,lcut(itype)
                lm1 = l1**2 + 1
                do n2 = 1,nindx(l2,itype)
                  do n1 = 1,nindx(l1,itype)
                    hmt(n1,n2,lm1:lm1+2*l1,lm2:lm2+2*l2,ic) =
     &                matmul( conjg(transpose(dwgn1(-l1:l1,-l1:l1,l1))), matmul (
     &                hmt(n1,n2,lm1:lm1+2*l1,lm2:lm2+2*l2,itype) , dwgn1(-l2:l2,-l2:l2,l2) ) )
                  enddo
                enddo
              enddo
            enddo
          else
            hmt(:,:,:,:,ic) = hmt(:,:,:,:,itype)
          endif
          ic = ic - 1
        enddo
      enddo
      Oend
      Nfence(hmt)

      Mdo( Rwrite(6,'(A'NoA) 'prep4 ' ; Rcall cpu_done(cputime) )

c
c     MPI: Estimate memory demand per subgroup (for pcycle)
# ifdef MPI
      Nfence(mem)
#   if MPI_order < 2
      Mwrt_n = nkpt1 ! for Mwrite_n
#   endif
      nbas1 = maxgpt + nlotot ! maximal nbas1 for xmem
      if(maxene/=-huge(0d0)) then
        nband1 = nbas1
      else
        if(l_soc) then ; nband1 = nband00/2 + 3 ! number of bands required for first variation (plus a few extra)
        else           ; nband1 = nband00   + 1
        endif
        nband1 = min( nband1 , nbas1 )
      endif
      xmem = (nbas1*(nbas1+1)) * MBYTES +                               ! hamiltonian & overlap
     &       maxgpt * (ncent+(maxlcut+1)**2+ntype*(maxlcut+1)) * 16d0 + ! harm & phas & apw
     &       maxgpt * ntype * 4d0 + nbas1 * nband1 * MBYTES +           ! pnt & evec
     &       nbas1**2 * MBYTES * 5.5d0                                  ! ScaLAPACK memory cost: 3 full matrices for mat1, olap1, and evec1, work arrays estimated to amount to 2 full matrices + 0.5 for safety
      xmem = xmem / ( maxmem - mem )
      call Mcast(xmem)
      Rwrite(6,'(/A/)') 'Auxiliary storage required per k point: '//Chr(nint(xmem*100))//'%'
      Rif(xmem>1) Warn('Memory demand likely too high.')
#   ifdef ALL_ON_ONE      
      xmem = 1 ! enforce single group
#     warning ALL_ON_ONE defined
#   endif
# endif

c
c     Loop over k points

      do ikpt = 1,nkpt1 ; PcycleM(ikpt,nkpt1,xmem,.false.)

# ifdef MPI
#   if   MPI_order == 0
        call begin_split_nodes(ikpt)
        Rcall Mwrite_n(0,0,Mwrt_n)
#   elif MPI_order == 1
        call begin_split_nodes(ikpt)
        Rcall Mwrite_n(0,ikpt,Mwrt_n)
#   else
        Mwrt_n = Msize
        ikpt0  = ikpt
        call Mcast(ikpt0)
        call begin_split_nodes(ikpt)
        Mwrt_n = Mwrt_n / Msize
        Rcall Mwrite_n(0,-(ikpt-ikpt0+1),Mwrt_n)
#   endif
# endif

        Rwrite(6,'(/A,A)') 'K point '//Chr(ikpt) MpiC('  (MPI:'//Chr(Msize)//')')
        if(ikpt<=nkpti) then
          kvec = kpt(:,ikpt)
        else
          kvec = kpt(:,ikpt+nkpt-nkpti)
        endif

        nbas1 = ngpt(ikpt) + nlotot
        Rwrite(6,'(A,I5,''  ('',I3,'')'')') '  Number of APWs (LOs): ',ngpt(ikpt),nlotot

        Nallocate0 (     overlap, (S_ nbas1*(nbas1+1)/2 S_) )
        Nallocate0 ( hamiltonian, (S_ nbas1*(nbas1+1)/2 S_) )

        call cpu_time(cputime)

c
c       Precalculation (k dependent)
        allocate  ( kgn(ngpt(ikpt)) )
        Allocate_ ( harm,((maxlcut+1)**2,ngpt(ikpt)) )   ; Wrtmem( harm = 0 )
        Allocate_ ( phas,(ncent,ngpt(ikpt)) )            ; Wrtmem( phas = 0 )
        Allocate_ ( apw,(2,0:maxlcut,ngpt(ikpt),ntype) ) ; Wrtmem( apw  = 0 )
        Allocate_ ( pnt,(ngpt(ikpt),ntype) )             ; Wrtmem( pnt  = 0 )
        do j = 1,ngpt(ikpt)
          ! harmonics * (-img)**l (->harm)
          kg1    = matmul(rlat,kvec+gpt1(:,j,ikpt))
          kgn(j) = sqrt(sum(kg1**2))
          call harmonicsr(harm(:,j),kg1,maxlcut)
          lm   = 0
          imgl = 1
          do l = 0,maxlcut
            harm(lm+1:lm+2*l+1,j) = harm(lm+1:lm+2*l+1,j) * imgl
            imgl                  = imgl * (-img)
            lm                    = lm + 2*l+1
          enddo
          ! APW matching coefficients (->apw)
          do itype = 1,ntype
            x = kgn(j) * grid(itype)%radius
            call sphbessel(jl,x,lcut(itype)+1)
            do l = 0,lcut(itype)
              dj = 0
              if(x==0) then ; if(l==1) dj = 1d0/3
              else          ;          dj = jl(l)*l/x - jl(l+1)
              endif
              apw(:,l,j,itype) = matmul(apwmat(:,:,l,itype),[jl(l),kgn(j)*dj])
            enddo
          enddo
          ! phases exp(iGRa) (->phas)
          do ic = 1,ncent
            phas(ic,j) = exp( img * 2*pi * dot_product(gpt1(:,j,ikpt),cent(:,ic)) )
          enddo
        enddo
        do ic = 1,ncent
          kphas(ic) = exp( img * 2*pi * dot_product(kvec,cent(:,ic)) )
        enddo
        ! Determine unique phases (exp(iGRa))
        Mdo( write(6,'(A'NoA) '1a' ; call cpu_done(cputime) )
        allocate(avg(ntype))
        avg = 0
        ic  = 0
        do itype = 1,ntype          
c          call rnorderpf(pnt(:,itype),[ transpose(real(phas(ic+1:ic+neq(itype),:))) ,
c     &                                  transpose(imag(phas(ic+1:ic+neq(itype),:))) ],ngpt(ikpt),2*neq(itype))
          call rnorderpf(pnt(:,itype), phas(ic+1,1), 2*ncent,ngpt(ikpt), ngpt(ikpt),2*neq(itype), .true. )
# ifdef test_phasegroup
#   warning trivial pointer for testing
          pnt(:,itype) = [(i,i=1,ngpt(ikpt))]
# endif
          k = 0
          i = 1
          allocate(phas0(neq(itype)))
          do while(i<=ngpt(ikpt))
            j     = i
            phas0 = phas(ic+1:ic+neq(itype),pnt(i,itype))
            do while(j<ngpt(ikpt))
              if( any(abs(phas(ic+1:ic+neq(itype),pnt(j+1,itype))-phas0)>1d-12) ) exit
              j = j + 1
            enddo
            avg(itype) = avg(itype) + j-i+1
            k          = k + 1
            i          = j + 1
          enddo
          deallocate(phas0)
          avg(itype) = avg(itype) / k
          ic         = ic + neq(itype)
        enddo        
        Mdo( write(6,'(A'NoA) '1' ; call cpu_done(cputime) )

        if(ntype==1) then ; Rwrite(6,'(A'NoA) '  Average size of G groups: '
        else              ; Rwrite(6,'(A'NoA) '  Average sizes of G groups: '
        endif
        do itype = 1,ntype ; Rwrite(6,'(A'NoA) Chr(nint(avg(itype)))//' ' ; enddo
        Rwrite(6,'(/A'NoA) '  Build FLAPW Hamiltonian... ' ; Mdo( write(6,*) )

c
c       IPW part of kinetic energy and overlap
        k = 0
        do j = Mstride(1,ngpt(ikpt))
          kg2n = kgn(j)
          do i = 1,j
            kg1n           = kgn(i)
            k              = j*(j-1)/2 + i
            mcdum          = cstep( gpt1(1,i,ikpt)-gpt1(1,j,ikpt), gpt1(2,i,ikpt)-gpt1(2,j,ikpt), gpt1(3,i,ikpt)-gpt1(3,j,ikpt) )
            overlap(k)     = mcdum
            hamiltonian(k) = mcdum * 0.25d0 * ( kg1n**2 + kg2n**2 )
          enddo
        enddo
        Mdo( write(6,'(A'NoA) '2' ; call cpu_done(cputime) )

c
c       IPW part of total potential
        k = 0
        do j = Mstride(1,ngpt(ikpt))
          do i = 1,j
            k              = j*(j-1)/2 + i
            g              = gpt1(:,i,ikpt) - gpt1(:,j,ikpt)
            hamiltonian(k) = hamiltonian(k) + vpw(g(1),g(2),g(3),ispin)
          enddo
        enddo
        Mdo( write(6,'(A'NoA) '4' ; call cpu_done(cputime) )

        Nfence( hamiltonian )
        Nfence( overlap     )

c
c       Add hmt to Hamiltonian

# ifdef CHECK_ifort18bug
        allocate(hmtpw(2,(maxlcut+1)**2,maxval(neq)),hmtpw0(2,(maxlcut+1)**2))
# endif

# ifdef INV
        ic = 0
        do itype = 1,ntype
          ic   = ic + 1
          linv = pcent(ic,invsym)==ic
          do ieq = 2,neq(itype)
            ic = ic + 1
            if( linv .neqv. (pcent(ic,invsym)==ic) )
     &        Bug('Found atom type with both inversion-invariant and inversion noninvariant sites. This is unexpected.')
          enddo
        enddo
# endif

        ! APW-APW        
        nullify(hmtpw_p)
        Mdo( lmn1 = 0 )
        Mdo( lmn2 = 0 )
        do j = Mstride(1,ngpt(ikpt))
c          Mdo( write(6,'(A'NoA) '8run '//trim(chr(j)) ; call cpu_done(cputime) )
          kg2  = matmul(rlat,kvec+gpt1(:,j,ikpt))
          kg2n = sum(kg2**2)
          if(kg2n>1d-25) then ; kg2 = kg2 / sqrt(kg2n)
          else                ; leg(0) = 1 ; leg(1:) = 0
          endif
          do itype = 1,ntype
            ylm2  = conjg( harm(:,j) )
            ic    = sum(neq(:itype-1))
            lm2   = (lcut(itype)+1)**2
            lm    = 0
            hmtpw = 0
            omtpw = 0
# ifdef INV
            allocate(invp(neq(itype)))
            invp = (/ ( pcent(ic+ieq,invsym)>=ic+ieq , ieq=1,neq(itype) ) /) ! INV: invp = inversion parents
            linv = pcent(ic+1,invsym)==ic+1                                  !      linv = atoms inversion-invariant (valid for the whole atom type, is checked above)
            neq1 = count(invp)
# else
            neq1 = neq(itype)
# endif
            if(neq1==1) then ; hmtpw_p => hmtpw (:,:,1)
            else             ; hmtpw_p => hmtpw0
            endif
            allocate(phas0(neq1),phas1(neq1),phas2(neq1))
            phas2 = ifInv( pack( phas(ic+1:ic+neq(itype),j) , invp ) , phas(ic+1:ic+neq(itype),j) )
            do l = 0,lcut(itype)
              a = apw(:,l,j,itype) * pi4vol**2
              do m = -l,l
                lm   = lm + 1
                ieq1 = 0
                do ieq = 1,neq(itype) ; ifInv( if(.not.invp(ieq)) cycle ; ieq1 = ieq1 + 1 , ieq1 = ieq )
                  hmtpw(:,:lm2,ieq1) = hmtpw(:,:lm2,ieq1) + (
     &                                 a(1) * hmt(:2,1,:lm2,lm,ic+ieq) +
     &                                 a(2) * hmt(:2,2,:lm2,lm,ic+ieq) ) * ( ylm2(lm) * phas2(ieq1) )
                enddo
              enddo
              omtpw(:,l) = (2*l+1) * ( a(1) * omt(:2,1,l,itype) +
     &                                 a(2) * omt(:2,2,l,itype) )
            enddo
            phas0 = 0
            do i1 = 1,ngpt(ikpt)
              i     = pnt(i1,itype) ; if(i>j) cycle
              phas1 = ifInv( pack( phas(ic+1:ic+neq(itype),i) , invp ) , phas(ic+1:ic+neq(itype),i) )
              Mdo( lmn1 = lmn1 + 1 )
              if( any(abs(phas0-phas1)>1d-12) ) then
                Mdo( lmn2 = lmn2 + 1 )
                phas0 = phas1
                mexp  = conjg(phas1(1)) * phas2(1)
                if(neq1>1) then
                  hmtpw0(:,:lm2) = hmtpw(:,:lm2,1)
                  do ieq = 2,neq1
                    cexp =        conjg(phas1(ieq)) * phas1(1) ! phas2(ieq) included in hmtpw, phas1(1) cancelled below                    
                    mexp = mexp + conjg(phas1(ieq)) * phas2(ieq)
                    if(abs(imag(cexp))<1d-12) then
                      if(abs(1+real(cexp))<1d-12) then ; hmtpw0(:,:lm2) = hmtpw0(:,:lm2) - hmtpw(:,:lm2,ieq) ! cexp = -1
                      else                             ; hmtpw0(:,:lm2) = hmtpw0(:,:lm2) + hmtpw(:,:lm2,ieq) ! cexp =  1
                      endif
                    else                               ; hmtpw0(:,:lm2) = hmtpw0(:,:lm2) + hmtpw(:,:lm2,ieq) * cexp
                    endif
                  enddo
                endif
              endif
              if(kg2n>1d-25) then
                kg1  = matmul(rlat,kvec+gpt1(:,i,ikpt))
                kg1n = sum(kg1**2)                
                if(kg1n>1d-25) then
                  kg1 = kg1 / sqrt(kg1n)
                  call legendre( leg , dot_product(kg1,kg2) , lcut(itype) )
                else
                  leg(0) = 1 ; leg(1:lcut(itype)) = 0
                endif
              endif
              lm    = 0
              cdum  = 0
              mcdum = 0
              do l = 0,lcut(itype)
                a     = apw(:,l,i,itype)
                mcdum = mcdum + leg(l) * ( a(1) * omtpw(1,l) + a(2) * omtpw(2,l) )
# ifdef INV
                if(linv) then
                  do m = -l,-1
                    lm   = lm + 1
                    cdum = cdum + harm(lm,i) * ( a(1) * hmtpw_p(1,lm) + a(2) * hmtpw_p(2,lm) )
                  enddo
                  lm   = lm + 1
                  cdum = cdum + harm(lm,i) * ( a(1) * hmtpw_p(1,lm) + a(2) * hmtpw_p(2,lm) ) * 0.5d0
                  lm   = lm + l
                else
# endif
                do m = -l,l
                  lm   = lm + 1
                  cdum = cdum + harm(lm,i) * ( a(1) * hmtpw_p(1,lm) + a(2) * hmtpw_p(2,lm) )
                enddo
                Inv( endif )
              enddo
              Inv( if(.not.linv) mcdum = mcdum * 2 )
              k              = (j-1)*j/2 + i
              hamiltonian(k) = hamiltonian(k) + cdum * conjg( phas1(1) )    Inv( *2 )
              overlap(k)     = overlap(k)     + mcdum * mexp * (4*pi)**(-1)
            enddo
            deallocate(phas0,phas1,phas2 InvC(invp) )            
            nullify(hmtpw_p)
          enddo
        enddo
        Mdo( write(6,'(A'NoA) '8' ; call cpu_done(cputime) )
        Mdo( write(6,'(A,2I6)') 'phase'//Chr(lmn1)//' '//Chr(lmn2) ; call cpu_done(cputime) )        

        ! APW-lo
        allocate ( hpwlo(Mcol1(ngpt(ikpt)),maxlmnlo,ncent) )
        allocate ( opwlo(Mcol1(ngpt(ikpt)),maxlmnlo,ncent) )
        hpwlo = 0
        opwlo = 0
        ic    = 0
        j     = 0
        do itype = 1,ntype
          do ieq = 1,neq(itype)
            ic   = ic + 1
            lm2  = 0
            lmn2 = 0
            do l2 = 0,lcut(itype)
              do m2 = -l2,l2
                lm2 = lm2 + 1
                do n2 = 3,nindx(l2,itype)
                  hmtpw0(:,:) = hmt(:2,n2,:,lm2,ic)
                  lmn2        = lmn2 + 1
                  j           = j + 1
                  k1          = (ngpt(ikpt)+j-1)*(ngpt(ikpt)+j)/2
                  do i = Mrange1(ngpt(ikpt))
                    k    = k1 + i
                    cexp = conjg( kphas(ic) * phas(ic,i) )
                    ylm  = harm(:,i)
                    lm   = 0
                    do l = 0,lcut(itype)
                      a    = apw(:,l,i,itype) * pi4vol
                      do m = -l,l
                        lm               = lm + 1
                        hpwlo(i,lmn2,ic) = hpwlo(i,lmn2,ic) + cexp * ylm(lm) *
     &                                     ( a(1) * hmtpw0(1,lm) + a(2) * hmtpw0(2,lm) )
                        if(l==l2.and.m==m2) then
                          opwlo(i,lmn2,ic) = opwlo(i,lmn2,ic) + cexp * ylm(lm) *
     &                                       ( a(1) * omt(1,n2,l2,itype) + a(2) * omt(2,n2,l2,itype) )
                        endif
                      enddo
                    enddo
                  enddo
                enddo
              enddo
            enddo
          enddo
        enddo
# ifdef CHECK_ifort18bug
        deallocate(hmtpw,hmtpw0)
# endif
        Mdo( write(6,'(A'NoA) '9' ; call cpu_done(cputime) )
        
# ifdef INV
        do ic = 1,ncent
          hpwlo(:,:,ic) = hpwlo(:,:,ic) * kphas(ic)
          opwlo(:,:,ic) = opwlo(:,:,ic) * kphas(ic)
        enddo
        ! symmetrize LOs
# if 1
        call symmetrize(hpwlo,size(hpwlo,1),maxlmnlo*ncent,2 + 3*(ntype+1)*2 + 3*(ntype+1)*symm_maxindx*maxlmnlo )
        call symmetrize(opwlo,size(opwlo,1),maxlmnlo*ncent,2 + 3*(ntype+1)*2 + 3*(ntype+1)*symm_maxindx*maxlmnlo )        
# else
        ! old version kept for testing...
        rdum = sqrt(0.5d0)!1/sqrt(2d0)
        ic   = 0
        do itype = 1,ntype
          do ieq = 1,neq(itype)
            ic  = ic + 1
            ic1 = pcent(ic,invsym)
            if(ic1>=ic) then
              lmn1 = 0
              do l = 0,lcut(itype)
                if(mod(l,2)==0) then ; imgl = 1
                else                 ; imgl = img
                endif
c                imgl = img**l
                do m = -l,l
                  ifac = (-1)**(l+m)
                  do n = 3,nindx(l,itype)
                    lmn1 = lmn1 + 1
                    lmn2 = lmn1 - 2*m*(nindx(l,itype)-2)
                    if(ic==ic1.and.m==0) then
                      do i = Mrange1(ngpt(ikpt))
                        hpwlo(i,lmn1,ic)  = hpwlo(i,lmn1,ic) * imgl
                        opwlo(i,lmn1,ic)  = opwlo(i,lmn1,ic) * imgl
                      enddo
                    else if(ic<ic1.or.m<0) then
                      do i = Mrange1(ngpt(ikpt))
                        cdum              = hpwlo(i,lmn1,ic)
                        hpwlo(i,lmn1,ic)  = ( cdum + ifac * hpwlo(i,lmn2,ic1) ) * rdum
                        hpwlo(i,lmn2,ic1) = ( cdum - ifac * hpwlo(i,lmn2,ic1) ) * (rdum * img)
                        cdum              = opwlo(i,lmn1,ic)
                        opwlo(i,lmn1,ic)  = ( cdum + ifac * opwlo(i,lmn2,ic1) ) * rdum
                        opwlo(i,lmn2,ic1) = ( cdum - ifac * opwlo(i,lmn2,ic1) ) * (rdum * img)
                      enddo
                    endif
                  enddo
                enddo
              enddo
            endif
          enddo
        enddo
# endif        
        Mdo( write(6,'(A'NoA) '10' ; call cpu_done(cputime) )
# endif
        ! copy hpwlo -> hamiltonian
        j  = ngpt(ikpt)
        ic = 0
        do itype = 1,ntype
          do ieq = 1,neq(itype)
            ic   = ic + 1 
            lmn2 = 0
            do l = 0,lcut(itype)
              do m = -l,l
                do n = 3,nindx(l,itype)
                  lmn2 = lmn2 + 1
                  j    = j + 1
                  k1   = j*(j-1)/2
                  do i = Mrange1(ngpt(ikpt))
                    k              = k1 + i
                    hamiltonian(k) = hpwlo(i,lmn2,ic)
                    overlap(k)     = opwlo(i,lmn2,ic)
                  enddo
                enddo
              enddo
            enddo
          enddo
        enddo
        deallocate ( hpwlo )
        deallocate ( opwlo )
        Mdo( write(6,'(A'NoA) '11' ; call cpu_done(cputime) )
        ! lo-lo
        allocate ( hlolo(maxlmnlo,ncent,maxlmnlo,ncent) )
        hlolo = 0
        j     = ngpt(ikpt)
        ic    = 0
        do itype = 1,ntype
          do ieq = 1,neq(itype)
            ic   = ic + 1 
            lm2  = 0
            lmn2 = 0
            do l2 = 0,lcut(itype)
              do m2 = -l2,l2
                lm2 = lm2 + 1
                do n2 = 3,nindx(l2,itype)
                  lmn2 = lmn2 + 1
                  j    = j + 1 ; Mcycle(j)
                  i    = ngpt(ikpt)
                  ! go to proper atom
                  do l = 1,itype-1
                    i = i + sum( [ (neq(l) * (nindx(l1,l)-2) * (2*l1+1), l1=0,lcut(l)) ] )
                  enddo
                  i = i + (ieq-1) * sum( [ ( (nindx(l1,itype)-2) * (2*l1+1), l1=0,lcut(itype)) ] )
                  !
                  lm1  = 0
                  lmn1 = 0
                  do l1 = 0,lcut(itype)
                    do m1 = -l1,l1
                      lm1 = lm1 + 1
                      do n1 = 3,nindx(l1,itype)
                        lmn1                   = lmn1 + 1
                        i                      = i + 1
                        hlolo(lmn1,ic,lmn2,ic) = hmt(n1,n2,lm1,lm2,ic)
                        if(l1==l2.and.m1==m2.and.i<=j) then
                          k          = j*(j-1)/2 + i
                          overlap(k) = overlap(k) + omt(n1,n2,l1,itype)
                        endif
                      enddo
                    enddo
                  enddo
                enddo
              enddo
            enddo
          enddo
        enddo
        Mpi( call Msum(hlolo) )
        Mdo( write(6,'(A'NoA) '12' ; call cpu_done(cputime) )
# ifdef INV
        ! symmetrize LOs
# if 1
        call symmetrize ( hlolo , maxlmnlo*ncent,maxlmnlo*ncent , 3 + 3*(ntype+1)*2 + 3*(ntype+1)*symm_maxindx*maxlmnlo )
# else
        ! old version kept for testing ...
        allocate ( carr(maxlmnlo,ncent) )
        rdum = sqrt(0.5d0)
        ic   = 0
        do itype = 1,ntype
          do ieq = 1,neq(itype)
            ic  = ic + 1
            ic1 = pcent(ic,invsym)
            if(ic1>=ic) then
              lmn1 = 0
              do l = 0,lcut(itype)
                if(mod(l,2)==0) then ; imgl = 1
                else                 ; imgl = -img
                endif
c                imgl = img**(-l)
                do m = -l,l
                  ifac = (-1)**(l+m)
                  do n = 3,nindx(l,itype)
                    lmn1 = lmn1 + 1
                    lmn2 = lmn1 - 2*m*(nindx(l,itype)-2)
                    if(ic==ic1.and.m==0) then
                      hlolo(lmn1,ic ,:,:) = hlolo(lmn1,ic,:,:) * imgl
                    else if(ic<ic1.or.m<0) then
                      carr                = hlolo(lmn1,ic,:,:)
                      hlolo(lmn1,ic ,:,:) =   ( carr + ifac * hlolo(lmn2,ic1,:,:) ) *   rdum
                      hlolo(lmn2,ic1,:,:) = - ( carr - ifac * hlolo(lmn2,ic1,:,:) ) * ( rdum * img )
                    endif
                  enddo
                enddo
              enddo
            endif
          enddo
        enddo
        ic = 0
        do itype = 1,ntype
          do ieq = 1,neq(itype)
            ic  = ic + 1
            ic1 = pcent(ic,invsym)
            if(ic1>=ic) then
              lmn1 = 0
              do l = 0,lcut(itype)
                if(mod(l,2)==0) then ; imgl = 1
                else                 ; imgl = img
                endif
c                imgl = img**l
                do m = -l,l
                  ifac = (-1)**(l+m)
                  do n = 3,nindx(l,itype)
                    lmn1 = lmn1 + 1
                    lmn2 = lmn1 - 2*m*(nindx(l,itype)-2)
                    if(ic==ic1.and.m==0) then
                      hlolo(:,:,lmn1,ic ) = hlolo(:,:,lmn1,ic) * imgl
                    else if(ic<ic1.or.m<0) then
                      carr                = hlolo(:,:,lmn1,ic)
                      hlolo(:,:,lmn1,ic ) = ( carr + ifac * hlolo(:,:,lmn2,ic1) ) *   rdum
                      hlolo(:,:,lmn2,ic1) = ( carr - ifac * hlolo(:,:,lmn2,ic1) ) * ( rdum * img )
                    endif
                  enddo
                enddo
              enddo
            endif
          enddo
        enddo
        deallocate ( carr )
# endif
        Mdo( write(6,'(A'NoA) '13' ; call cpu_done(cputime) )
# endif
        ! copy hlolo -> hamiltonian
        j   = ngpt(ikpt)
        ic1 = 0
        do itype1 = 1,ntype
          do ieq1 = 1,neq(itype1)
            ic1  = ic1 + 1 
            lmn1 = 0
            do l1 = 0,lcut(itype1)
              do m1 = -l1,l1
                do n1 = 3,nindx(l1,itype1)
                  lmn1 = lmn1 + 1
                  j    = j + 1 ; Mcycle(j)
                  i    = ngpt(ikpt)
                  ic2  = 0
                  do itype2 = 1,ntype
                    do ieq2 = 1,neq(itype2)
                      ic2  = ic2 + 1
                      lmn2 = 0
                      do l2 = 0,lcut(itype2)
                        do m2 = -l2,l2
                          do n2 = 3,nindx(l2,itype2)
                            lmn2           = lmn2 + 1
                            i              = i + 1 ; if(i>j) cycle
                            k              = j*(j-1)/2 + i
                            hamiltonian(k) = hlolo(lmn1,ic1,lmn2,ic2)
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
        deallocate ( hlolo )
        Mdo( write(6,'(A'NoA) '14' ; call cpu_done(cputime) )
        Rcall cpu_done(cputime)

        Nfence( hamiltonian ) ; MpiO( call Msum(hamiltonian,comm=Ocomm) ) ; Nfence( hamiltonian )
        Nfence( overlap     ) ; MpiO( call Msum(    overlap,comm=Ocomm) ) ; Nfence( overlap     )

c       Diagonalize Hamiltonian
        Rwrite(6,'(A'NoA) '  Diagonalize Hamiltonian... '
        if(maxene/=-huge(0d0)) then
          nband1 = nbas1
          ldum   = .false.
        else
          if(l_soc) then ; nband1 = nband00/2 + 3 ! number of bands required for first variation (plus a few extra)
          else           ; nband1 = nband00   + 1
          endif
          nband1 = min( nband1 , nbas1 )
          ldum   = nband1<nbas1
        endif
        allocate ( eig(nbas1) )
        Nallocate ( evec, (S_ nbas1,nband1 S_) )
        call Mdiagonalize(Win(evec),eig,hamiltonian,overlap,1,nband1,-huge(0d0),maxene)
        Nfence(evec)
        Mpi( Ocall Msum(evec,comm=Ocomm) )
        Nfence(evec)        
        Rif(minval(eig(:nband1))<minebas(ispin)-0.5d0) then
          Warn('Lowest eigenenergy more than 0.5 Ha below lowest energy parameter. You might have picked up a ghost state.')
        endif
        Ndeallocate ( hamiltonian )
        Ndeallocate ( overlap     )
        if(.not.l_soc.and.nband1<nband0) Error('Did not find enough bands. Decrease NBAND.')
        ilow = count(eig(:nband1)<it_elow)
        Nfence(evec)
        if(ilow>0) then
          nband1              = nband1 - ilow
          eig(:nband1)        = eig(ilow+1:ilow+nband1)
          ifO evec(:,:nband1) = evec(:,ilow+1:ilow+nband1)
        endif
        Nfence(evec)
        if(ldum) nband1 = nband1 - 1
        if(it_mode<=2) then
          maxeband = min(maxeband,nband1)
          ifR ene(:nband1,ikpt,ispin) = eig(:nband1)
        endif
        if(ldum) then
          call cut_deg(nband1,eig)
          if(l_soc.and.nband1>nband00/2) then
            do while(nband1>=nband00/2)
              n      = nband1
              nband1 = nband1 - 1
              call cut_deg(nband1,eig)
            enddo
            nband1 = n
          endif
        endif
        Rcall cpu_done(cputime)
        Rif(ilow>0) write(6,'(I4,A,F10.5)') ilow,' eigenvalues below',it_elow
        Rwrite(6,'(A)') '  Kohn-Sham energies ('//Chr(nband1)//'):'
        Rwrite(6,'(5F15.10)') eig(:nband1)
        nband(ikpt,ispin) = nband1
        if(it_stop/=0.and.it_mode<=2.and..not.bandinfo) goto 4

c       Calculate cmt array
        Rwrite(6,'(A'NoA) '  Construct cmt array... '
        Mdo( write(6,'(A'NoA) '15_' ; call cpu_done(cputime) )
        ! APWs
        allocate(cmt0(Mcol1(nband1),maxlmindx,ncent))
        allocate(cpw0(Mcol1(nband1),ngpt(ikpt)))
        allocate(mcmt(Mcol1(nband1),-2*maxlcut:2*maxlcut+1))

c        if(size(cmt0,1)==0) goto 5

# ifdef test_phasegroup
        do n1 = 1,50
# endif
        cmt0 = 0        
        ic   = 0
        do itype = 1,ntype

          if(avg(itype)>2) then

          !
          ! algorithm with phase groups
          !
          
          allocate(phas0(neq(itype)))
# ifdef __INTEL_COMPILER
          do i = 1,ngpt(ikpt) ! workaround to avoid a segfault caused by ifort in the diagonalization routines below.
            cpw0(Mcol1(nband1),i) = evec(pnt(i,itype),Mcol1(nband1))
          enddo
# else
          cpw0(Mcol1(nband1),:) = transpose(evec(pnt(:,itype),Mcol1(nband1)))
# endif          
# ifdef test_phasegroup
          do i = 1,ngpt(ikpt),n1
            j     = min( i-1+n1 , ngpt(ikpt) )
            phas0 = phas(ic+1:ic+neq(itype),pnt(i,itype))
# else
          i = 1
          do while(i<=ngpt(ikpt))
            ! G groups with same phases (->i:j)
            phas0 = phas(ic+1:ic+neq(itype),pnt(i,itype))
            j     = i
            do while(j<ngpt(ikpt))
              if( any(abs(phas(ic+1:ic+neq(itype),pnt(j+1,itype))-phas0)>1d-12) ) exit
              j = j + 1
            enddo
# endif
            ! Determine transformation matrix G -> lmn (->gtrafo)
            imgl = 1
            lmn  = 1
            allocate(gtrafo(i:j,-2*maxlcut:2*maxlcut+1))                         
            do l = 0,lcut(itype) 
              lm1 = (l+1)*l + 1 ! m=0
              lm2 = (l+1)**2    ! m=l
              do k = i,j
                k1       = pnt(k,itype)
                a        = apw(:,l,k1,itype) * pi4vol
                ylm1(:l) = harm(lm1:lm2,k1) * imgl ! note that harm is defined as (-img)**l * Ylm; multiplying with imgl=img**l eliminates (-img)**l
                do m = 0,l
                  gtrafo(k,2*m:2*m+1)   = a * real( ylm1(m) ) ! m>=0 | we use real harmonics, which are the real/imag. parts of the 
                  if(m>0)                                            ! complex harmonics Ylm with m>=0; factor sqrt(2)*(-1)**m omitted
     &            gtrafo(k,-2*m:-2*m+1) = a * imag( ylm1(m) ) ! m<0
                enddo
              enddo              
              ! perform trafo (->mcmt)
              if(size(cpw0,1)>0) mcmt(:,-2*l:2*l+1) = matmat ( cpw0(:,i:j) , gtrafo(:,-2*l:2*l+1) )
              ! distribute phases (->cmt0)
              do ieq = 1,neq(itype)
                cexp = kphas(ic+ieq) * phas0(ieq) * imgl ! includes factor img**l
                lmn1 = lmn
                lmn2 = -2*l
                do m = -l,l
                  cmt0(:,lmn1:lmn1+1,ic+ieq) = cmt0(:,lmn1:lmn1+1,ic+ieq) + mcmt(:,lmn2:lmn2+1) * cexp
                  lmn1                       = lmn1 + nindx(l,itype)
                  lmn2                       = lmn2 + 2
                enddo
              enddo
              lmn  = lmn + (2*l+1) * nindx(l,itype)
              imgl = imgl * img
            enddo
            deallocate(gtrafo)
# ifndef test_phasegroup            
            i = j + 1
# endif
          enddo
          deallocate(phas0)
          ! transform real harmonics to complex harmonics (->cmt0)
          allocate(cmt1(Mcol1(nband1),2))
          do ieq = 1,neq(itype)
            lmn1 = 0
            do l = 0,lcut(itype)
              do m = l,1,-1
                lmn2                         = lmn1 + 2*m * nindx(l,itype)  ! lmn1 = offset at -m ; lmn2 = offset at +m
                cmt1(:,:2)                   = cmt0(:,lmn1+1:lmn1+2,ic+ieq) ! store coefficients at -m
                cmt0(:,lmn1+1:lmn1+2,ic+ieq) = (-1)**m * ( cmt0(:,lmn2+1:lmn2+2,ic+ieq) + img * cmt0(:,lmn1+1:lmn1+2,ic+ieq) ) ! transform to get conjg(Y(l,-m)) = Y(l,m)  * (-1)**m
                cmt0(:,lmn2+1:lmn2+2,ic+ieq) =             cmt0(:,lmn2+1:lmn2+2,ic+ieq) - img * cmt1(:,:2)                     ! transform to get conjg(Y(l,m))  = Y(l,-m) * (-1)**m
                lmn1                         = lmn1 + nindx(l,itype)
              enddo
              lmn2 = lmn1 + (l+1) * nindx(l,itype)
              lmn1 = lmn2 ! skip m=0..l
            enddo
          enddo
          deallocate(cmt1)
c          write(*,*) sum(abs(cmt0)),sum(cmt0)

          else
c          cmt0 = 0

          !
          ! simple algorithm
          !

          cpw0(Mcol1(nband1),:) = transpose(evec(:ngpt(ikpt),Mcol1(nband1)))
          do j = 1,ngpt(ikpt)
            do ieq = 1,neq(itype)
              cexp2(ieq) = kphas(ic+ieq) * phas(ic+ieq,j)
            enddo
            ylm  = conjg( harm(:,j) )
            lm   = 0
            lmn1 = 0
            do l = 0,lcut(itype)
              a = apw(:,l,j,itype) * pi4vol
              do m = -l,l
                lm = lm + 1
                do ieq = 1,neq(itype)
                  cexp                  = ylm(lm) * cexp2(ieq)
                  cmt0(:,lmn1+1,ic+ieq) = cmt0(:,lmn1+1,ic+ieq) + cpw0(:,j) * ( a(1) * cexp )
                  cmt0(:,lmn1+2,ic+ieq) = cmt0(:,lmn1+2,ic+ieq) + cpw0(:,j) * ( a(2) * cexp )
                enddo
                lmn1 = lmn1 + nindx(l,itype)                
              enddo
            enddo
          enddo

          endif ! end algorithm

c          write(*,*) sum(abs(cmt0)),sum(cmt0)
c          call finish

          ic = ic + neq(itype)          
        enddo ! end loop over atom types
        
        Mdo( write(6,'(A'NoA) '15' ; call cpu_done(cputime) )

        deallocate(mcmt)
        cpw0(Mcol1(nband1),:) = transpose(evec(:ngpt(ikpt),Mcol1(nband1)))

# ifdef test_phasegroup
        enddo
# endif
        
        ! LOs
        ! (Desymmetrization explicitly done here instead of by "desymmetrize" because of the array structure of cmt0)
        j  = ngpt(ikpt)
        ic = 0
        do itype = 1,ntype
          do ieq = 1,neq(itype)
            ic   = ic + 1                  ; Inv( ic1   = pcent(ic,invsym) )
            cexp = kphas(ic) * sqrt(0.5d0) ; Inv( cexp1 = kphas(ic1) * sqrt(0.5d0) )
            lmn1 = 0
            do l = 0,lcut(itype)
              if(mod(l,2)==0) then ; cdum = kphas(ic)
              else                 ; cdum = kphas(ic) * img
              endif
              do m = -l,l
                ifac = (-1)**(l+m)
                do n = 1,nindx(l,itype)
                  lmn1 = lmn1 + 1
                  lmn2 = lmn1 - 2*m*nindx(l,itype)
                  if(n>2) then
                    j = j + 1
                    do k = Mrange1(nband1)
# ifdef INV
                      if(ic==ic1.and.m==0) then
                        cmt0(k,lmn1,ic ) = cmt0(k,lmn1,ic ) + cdum  *   evec(j,k)
                      else if(ic<ic1.or.ic==ic1.and.m<0) then
                        cmt0(k,lmn1,ic ) = cmt0(k,lmn1,ic ) + cexp  *   evec(j,k)
                        cmt0(k,lmn2,ic1) = cmt0(k,lmn2,ic1) + cexp1 * ( evec(j,k) * ifac )
                      else if(ic>ic1.or.ic==ic1.and.m>0) then
                        cmt0(k,lmn1,ic ) = cmt0(k,lmn1,ic ) - cexp  * ( evec(j,k) * ifac ) * img
                        cmt0(k,lmn2,ic1) = cmt0(k,lmn2,ic1) + cexp1 *   evec(j,k)          * img
                      endif
# else
                      cmt0(k,lmn1,ic) = evec(j,k)
# endif
                    enddo
                  endif
                enddo
              enddo
            enddo
          enddo
        enddo

 5      Rcall cpu_done(cputime)

c       Second variation for spin-orbit coupling
        if(it_mode==3) then
          Rif(ispin==2) Error('SOC not implemented for spin-polarized systems.')
          Rwrite(6,'(A'NoA) '  Build SOC Hamiltonian  ... '          
          Nallocate0 ( hamiltonian_soc, (S_ 2*nband1,2*nband1 S_) )
          allocate ( hlp(nband1) )
          
          do i = Mrange1(nband1)
            hamiltonian_soc(i,i)               = eig(i)
            hamiltonian_soc(nband1+i,nband1+i) = eig(i)
          enddo

          ! Define dwgn1 for rotation if sqaxis != zaxis
          if(any(sqa/=0)) then
            call rotation_sph(rot,sqa(1),sqa(2))
            call dwigner(dwgn1,rot,maxlcut)
          endif

          do itype = 1,ntype
            call derivative(dvmt,vmt(:,1,itype,1),itype)
            do i = 1,grid(itype)%number
              dvmt(i) = dvmt(i) / rgrid(i,itype)
            enddo

            do l = 1,lcut(itype)
              n = nindx(l,itype)              
              allocate ( cmt1(n,nband1),cmt2(n,nband1),cmt3(n,nband1),cmt4(n,nband1) )
              ! radial integral INT r² u1(r) u2(r) 1/(2Mc)² r² dV/dr dr
              integr = 0
              rdum   = 2 * clight**2
              do i = 1,nindx(l,itype)
                do j = 1,nindx(l,itype)
                  dvmt_so     = dvmt * ( 0.5d0 / (1+(ebas(i,l,itype,1)-vmt(:,1,itype,1))/rdum)**2
     &                                 + 0.5d0 / (1+(ebas(j,l,itype,1)-vmt(:,1,itype,1))/rdum)**2 ) / (2*rdum)
                  integr(i,j) = intgrf( dvmt_so * ( bas1(:,i,l,itype,1) * bas1(:,j,l,itype,1) ) , itype )
                enddo
              enddo

              ! perform L*S and multply with cmt
              ic = sum(neq(:itype-1))
              do ieq = 1,neq(itype)
                ic = ic + 1
                lm = sum( [ ((2*i+1)*nindx(i,itype),i=0,l-1) ] )
                do m = -l,l
                  do n1 = Mrange1(nband1)
                    if(m/=-l) then ; lm1 = lm - n  ; cmt1(:n,n1) = conjg( cmt_trafo(lm1,m-1) ) ; endif
                                     lm1 = lm      ; cmt2(:n,n1) = conjg( cmt_trafo(lm1,m)   )
                    if(m/=l)  then ; lm1 = lm + n  ; cmt3(:n,n1) = conjg( cmt_trafo(lm1,m+1) ) ; endif
                                                     cmt4(:n,n1) = matmul ( integr(:n,:n) , conjg(cmt2(:n,n1)) )
                  enddo
                  MrangeDistr( cmt1(:, McolD1(nband1,i) ) , i )
                  MrangeDistr( cmt2(:, McolD1(nband1,i) ) , i )
                  MrangeDistr( cmt3(:, McolD1(nband1,i) ) , i )
                  MrangeDistr( cmt4(:, McolD1(nband1,i) ) , i )

                  do n2 = Mrange1(nband1)
                    ! up/up + down/down
                    hlp                                  = matmul ( cmt4(:,n2) , cmt2 )                   
                    hamiltonian_soc(:nband1,n2)          = hamiltonian_soc(:nband1,n2)          + m * hlp
                    hamiltonian_soc(nband1+1:,nband1+n2) = hamiltonian_soc(nband1+1:,nband1+n2) - m * hlp
                    ! up/down
                    if(m/=-l) then
                      hamiltonian_soc(:nband1,nband1+n2) = hamiltonian_soc(:nband1,nband1+n2) + sqrt( l*(l+1d0) - m*(m-1) ) *
     &                                                     matmul ( cmt4(:,n2) , cmt1 )
                    endif
                    ! down/up
                    if(m/=l) then
                      hamiltonian_soc(nband1+1:,n2)      = hamiltonian_soc(nband1+1:,n2)      + sqrt( l*(l+1d0) - m*(m+1) ) *
     &                                                     matmul ( cmt4(:,n2) , cmt3 )
                    endif
                  enddo
                    
                  lm = lm + nindx(l,itype)

                enddo
              enddo

              deallocate ( cmt1,cmt2,cmt3,cmt4 )
            enddo
          enddo          
          deallocate ( eig,hlp )
          Nfence(hamiltonian_soc) ; MpiO( call Msum(hamiltonian_soc,comm=Ocomm) ) ; Nfence(hamiltonian_soc)          
          Rcall cpu_done(cputime)
          Rwrite(6,'(A'NoA) '  Diagonalize Hamiltonian... '
          nband_soc = 2 * nband1
          Nallocate ( evec_soc, (S_ nband_soc,nband_soc S_) )
          allocate ( eig(nband_soc) )
          Mdo( write(6,'(A'NoA) 'before diag' ; call cpu_done(cputime) )
          call Mdiagonalize ( Win(evec_soc),eig,hamiltonian_soc,1,nband_soc,-huge(0d0),maxene)
          Nfence(evec_soc)
          Mpi( Ocall Msum(evec_soc,comm=Ocomm) )
          Nfence(evec_soc)          
          Ndeallocate( hamiltonian_soc )
          Rcall cpu_done(cputime)
          if(nband_soc>nband00) then
            nband_soc = nband00 + 1
            call cut_deg(nband_soc,eig)
          endif
          nband(ikpt,ispin) = nband_soc
          Rwrite(6,'(A)') '  Kohn-Sham energies with SOC ('//Chr(nband_soc)//'):'
          Rwrite(6,'(5F15.10)') eig(:nband_soc)
          maxeband = min(maxeband,nband_soc)
          ifR ene(:nband_soc,ikpt,ispin) = eig(:nband_soc)
        endif

        ! Undo local-orbital construction in coefficients
        if(.not.it_stop==2) then ! if it_stop=2, we need the coefficients before (!) the lo construction of basis is undone
          do itype = 1,ntype
            do l = 0,lcut(itype)
              do n = 3,nindx(l,itype)
                ulo  =  ubas(n,l,itype,ispin)
                dulo = dubas(n,l,itype,ispin)
                a    = matmul(apwmat(:,:,l,itype),-[ulo,dulo])
                ic   = sum(neq(:itype-1))
                do ieq = 1,neq(itype)
                  ic = ic + 1
                  do m = 0,2*l
                    lmn1 = sum( [ ( (2*l1+1)*nindx(l1,itype),l1=0,l-1 ) ] ) + m*nindx(l,itype)
                    do k = Mrange1(nband1)
                      cmt0(k,lmn1+1,ic) = cmt0(k,lmn1+1,ic) + a(1) * cmt0(k,lmn1+n,ic)
                      cmt0(k,lmn1+2,ic) = cmt0(k,lmn1+2,ic) + a(2) * cmt0(k,lmn1+n,ic)
                    enddo
                  enddo
                enddo
              enddo
            enddo
          enddo
        endif

        if(it_stop==0.or.it_stop==2.and.bandinfo) then
          Rif(it_stop==2) allocate(cmt(maxlmindx,ncent,nband(ikpt,ispin),1,nspin3))
          if(it_mode==3) then
            Rwrite(6,'(A'NoA) 'Construct cmt/cpw arrays... '
            if(it_stop/=2) then
              allocate(hlp2(ngpt(ikpt),packet),carr(size(cpw0,2),size(cpw0,1)))
              carr = transpose(cpw0)
              do n = 1,nband_soc,packet
                i = min(nband_soc-n+1,packet)
                if(size(carr)==0) then ; hlp2(:,:i)  = 0
                else                   ; hlp2(:,:i) = matmat ( carr , evec_soc(        Mcol1(nband1)        ,n:n+i-1) )              
                endif
                Mpi( call Msum(hlp2(:,:i)) )
                ifR cpw(:ngpt(ikpt),n:n+i-1,ikpt,1) = hlp2(:,:i)
              enddo
              do n = 1,nband_soc,packet
                i = min(nband_soc-n+1,packet)
                if(size(carr)==0) then ; hlp2(:,:i) = 0
                else                   ; hlp2(:,:i) = matmat ( carr , evec_soc( nband1+Mcol1(nband1)+nband1 ,n:n+i-1) )
                endif
                Mpi( call Msum(hlp2(:,:i)) )
                ifR cpw(:ngpt(ikpt),n:n+i-1,ikpt,2) = hlp2(:,:i)
              enddo
              deallocate(hlp2,carr)
              k = ikpt
            else
              k = 1
            endif
            allocate(hlp2(maxlmindx,10),carr(size(cmt0,2),size(cmt0,1)))
            do ic = 1,ncent
              carr = transpose(cmt0(:,:,ic))
              do n = 1,nband_soc,packet
                i = min(nband_soc-n+1,packet)
                if(size(carr)==0) then ; hlp2(:,:i) = 0
                else                   ; hlp2(:,:i) = matmat ( carr , evec_soc(        Mcol1(nband1)        ,n:n+i-1) )
                endif
                Mpi( call Msum(hlp2(:,:i)) )
                ifR cmt(:,ic,n:n+i-1,k,1) = hlp2(:,:i)
              enddo
              do n = 1,nband_soc,packet
                i = min(nband_soc-n+1,packet)                                
                if(size(carr)==0) then ; hlp2(:,:i) = 0
                else                   ; hlp2(:,:i) = matmat ( carr , evec_soc( nband1+Mcol1(nband1)+nband1 ,n:n+i-1) )
                endif
                Mpi( call Msum(hlp2(:,:i)) )
                ifR cmt(:,ic,n:n+i-1,k,2) = hlp2(:,:i)
              enddo
            enddo
            deallocate(hlp2,carr)
            Rcall cpu_done(cputime)
          else
            if(it_stop==2) then ! only cmt is needed, which is allocated only in root
# ifdef MPI
              m = maxlmindx * ncent ! all other procs have to send their share (cmt0) to root
              do i = 0,Msize-1
                do n = 1+i*nband1/Msize,(1+i)*nband1/Msize
                  if(Mrank/=0) then
                    if(i==Mrank)    call mpi_send(cmt0(n,:,:),   m,mpi_double_complex,0,1,Mcomm,      Merr)
                  else
                    if(i/=0) then ; call mpi_recv(cmt(:,:,n,1,1),m,mpi_double_complex,i,1,Mcomm,Mstat,Merr)
                    else          ; cmt(:,:,n,1,1) = cmt0(n,:,:)
                    endif
                  endif
                enddo
              enddo
# else
              do n = 1,nband1
                cmt(:,:,n,1,1) = cmt0(n,:,:)
              enddo
# endif
            else
              do n = Mrange1(nband1)
                cpw(:ngpt(ikpt),n,ikpt,ispin) = cpw0(n,:)
                cmt(:,:,        n,ikpt,ispin) = cmt0(n,:,:)
              enddo
            endif
          endif
          Mdo( write(6,'(A'NoA) '16' ; call cpu_done(cputime) )
          if(it_stop==2) then
            Rbegin
            call band_info(lines,cmt(:,:,:,1,:),spin=ispin,result=proj)
            if(.not.allocated(project)) then
              allocate(project(size(proj,1),nspin1,size(ene,1),nkpt))
              project = 0
            endif
            project(:,ispin,:size(proj,3),ikpt) = proj(:,1,:)
            deallocate(proj,cmt)
            Rend
          endif
        endif

        MnoR( nband(ikpt,ispin) = 0 ) ! non-root processes must "forget" nband(ikpt,ispin) because of Msum below

        if(it_mode==3) tNdeallocate ( evec_soc )
        deallocate ( cmt0,cpw0 )
 4      deallocate ( eig,kgn,avg )
        Deallocate_( harm )
        Deallocate_( phas )
        Deallocate_( apw  )
        Deallocate_( pnt  )
        Ndeallocate( evec )

# ifdef MPI
#   if   MPI_order == 0
        Rcall Mwrite_n(1,0,Mwrt_n)
#   elif MPI_order == 1
        Rcall Mwrite_n(1,ikpt,Mwrt_n)
#   else
        Rcall Mwrite_n(1,-(ikpt-ikpt0+1),Mwrt_n)
        Rcall Mwrite_n(2,-(ikpt-ikpt0+1),Mwrt_n)
#   endif
        call end_split_nodes        
# endif
        
      enddo ! end loop over k points

# if defined(MPI) && MPI_order < 2
      call Mwrite_n(2, MPI_order ,Mwrt_n)
# endif

      Ndeallocate(hmt)

c     Undo local-orbital construction in basis functions
      do itype = 1,ntype
        do l = 0,lcut(itype)
          do n = 3,nindx(l,itype)
            ulo                     =  ubas(n,l,itype,ispin)
            dulo                    = dubas(n,l,itype,ispin)
            a                       = matmul(apwmat(:,:,l,itype),-[ulo,dulo])
            bas1(:,n,l,itype,ispin) = - a(1) * bas1(:,1,l,itype,ispin)
     &                                - a(2) * bas1(:,2,l,itype,ispin)
     &                                +        bas1(:,n,l,itype,ispin)
            bas2(:,n,l,itype,ispin) = - a(1) * bas2(:,1,l,itype,ispin)
     &                                - a(2) * bas2(:,2,l,itype,ispin)
     &                                +        bas2(:,n,l,itype,ispin)
          enddo
        enddo
      enddo
      
      enddo ! end loop over spins

      Rcall cpu_time(cputime)
      Rwrite(6,'(/A,F12.5)') 'Timing (iterate):',cputime-cputime0

# ifdef MPI
      call Msum(nband)
      call Msum(ene)
      call mpi_allreduce(mpi_in_place,maxeband,1,mpi_integer,mpi_min,Mcomm,Merr)

      if(it_stop==0) then
        Nfence(cpw)
        Nfence(cmt)
        Ocall Msum(cpw,comm=Ocomm)
        Ocall Msum(cmt,comm=Ocomm)
        Nfence(cpw)
        Nfence(cmt)
      endif

      if(it_stop==2.and.bandinfo) then ! Collect contribution of each "former" root
        do i = 1,Msize-1
          Mcolor = mpi_undefined ; if(any(Mrank==[0,i])) Mcolor = 0
          call mpi_comm_split(Mcomm,Mcolor,Mrank,Mcoms(i),Merr)
        enddo
        if(Mrank/=0) then
          ldum = allocated(project)
          call mpi_send(ldum,1,mpi_logical,0,0,Mcomm,Merr)
          if(ldum) then
            call Msum(project,rank=0,comm=Mcoms(Mrank))
            deallocate(project)
          endif
        else
          do i = 1,Msize-1
            call mpi_recv(ldum,1,mpi_logical,i,0,Mcomm,Mstat,Merr)
            if(ldum) call Msum(project,rank=0,comm=Mcoms(i))
          enddo
        endif
        do i = 1,Msize-1
          if(Mcoms(i)/=mpi_comm_null) call mpi_comm_free(Mcoms(i),Merr)
        enddo
      endif
# endif

      maxband  = maxval(nband)
      maxeband = min(maxeband,maxband)

c
c     Write out band structure      
      Rif(it_stop==2) then
        iunit = fopen('bands',status='unknown')
        allocate(qpath(0:nkpt))
        call bands_header(iunit,'# Band structure (ITERATE)',qpath,kpt,nkpt)        
        if(bandinfo) then
          write(iunit,'(A,9X,A,8X,A'NoA) '#','     ','       ' ; if(nspin1==2) write(iunit,'(8X,A'NoA) '       '
          call write_orbitals_header(iunit,lines,size(lines),nspin1,10,1) ; write(iunit,*)
          write(iunit,'(A,9X,A,8X,A'NoA) '#','kpath',' energy' ; if(nspin1==2) write(iunit,'(8X,A'NoA) '       '
          call write_orbitals_header(iunit,lines,size(lines),nspin1,10,2) ; write(iunit,*)
        else
          write(iunit,'(A,9X,A,8X,A)')   '#','kpath',' energy'
        endif
        if(bandinfo) then ; m = minval(nband(:,:nspin1)) ! maximal common band index of calculated wavefunctions
        else              ; m = maxeband                 ! maximal band index of calculated energies
        endif
        if(m==0) Bug('Minimal number of bands is zero: minval(nband)==0')
        do i = 1,m
          do ikpt = 1,nkpt
            write(iunit,'('//Chr(1+nspin1)//'F15.9'NoA) qpath(ikpt),ene(i,ikpt,:)*escale
            if(bandinfo) then
              write(iunit,'('//Chr(nspin1*size(project,1))//'F10.4'NoA) transpose(project(:,:,i,ikpt))
            endif
            write(iunit,*)
          enddo
          write(iunit,*)
          write(iunit,*)
        enddo
        call fclose(iunit)
        if(bandinfo) deallocate(project)
        write(6,'(/A)') 'Band structure written to bands'
      endif

      if(it_stop==2.and.bandinfo) deallocate(lines)

      Rwrite(6,*)
      if(it_stop/=0) then
        Rwrite(6,'(A)') 'ITERATE finished.'
        call finish
      endif

c
c     RESTART: Prepare output file
c
      if(iand(restart,W_KS)/=0) then
        beginSingle
        call kpt_reorder1(0)
        call write_param(l_soc)
        do ispin = 1,nspin1
          do ikpt = 1,nkpt1
            n = nband(ikpt,ispin)
            if(l_soc) then
              call write_wavef(ikpt,ispin,ngpt(ikpt),gpt1(:,:ngpt(ikpt),ikpt),
     &                         ene(:n,ikpt,ispin),cmt(:,:,:n,ikpt,:),    cpw(:,:n,ikpt,:))
            else
              call write_wavef(ikpt,ispin,ngpt(ikpt),gpt1(:,:ngpt(ikpt),ikpt),
     &                         ene(:n,ikpt,ispin),cmt(:,:,:n,ikpt,ispin),cpw(:,:n,ikpt,ispin))
            endif
          enddo
        enddo
        call kpt_reorder2(0)
        endSingle
      endif

      MnoR( deallocate(ene) )

# endif

c --------------------------------------

      contains

      function cmt_trafo(off,m)
      implicit none
      integer, intent(in) :: off,m
      complex_dp          :: cmt_trafo(n)
      integer             :: i
      if(all(sqa==0)) then
        cmt_trafo = cmt0(n1,off+1:off+n,ic)
      else
        do i = 1,n
          cmt_trafo(i) = dot_product ( dwgn1(-l:l,m,l) , cmt0(n1,off-(l+m)*n+i:off+(l-m)*n+i:n,ic) )
        enddo
      endif
      end function cmt_trafo

c --------------------------------------

      subroutine cut_deg(nband1,eig)
      implicit none
      real_dp, intent(in)    :: eig(*)
      integer, intent(inout) :: nband1
      do
        if(nband1==0) Error('Reached zero bands while cutting degenerate subspace.')
        if(abs(eig(nband1+1)-eig(nband1))>edeg) exit
        nband1 = nband1 - 1
      enddo
      end subroutine cut_deg

c --------------------------------------

      end

c --------------------------------------

# else
      Error('ITERATE not implemented with -DTDDFT.')
      end
# endif

# ifndef INV_SOC

c --------------------------------------

c     Sums vxc and qsgw matrices and returns the result in vxcmat(:nb,:nb).
c
      subroutine calc_vxc(vxcmat,band,nb,ikpt,ispin,l_qsgw)
      use readwrite, only: read_qsgw
      use wrapper,   only: unpackmat
      use, intrinsic :: iso_fortran_env
      implicit none
      integer,     intent(in)  :: nb,band(nb),ikpt,ispin
      logical,     intent(in)  :: l_qsgw
      MCOMPLEX_dp, intent(out) :: vxcmat(nb,nb)
      MCOMPLEX_dp              :: vxchlp(nb*(nb+1)/2)
      call calc_vxc_p(vxchlp,band,nb,ikpt,ispin)
      vxcmat = unpackmat(vxchlp)
      if(l_qsgw) call read_qsgw(vxcmat,band,nb,ikpt,ispin)      
      end

c --------------------------------------

c     Calculates vxc matrix and returns it in packed storage: vxcmat(:nb*(nb+1)/2). (Does not add qsgw matrix.)
c
      subroutine calc_vxc_p(vxcmat,band,nb,ikpt,ispin)
      use readwrite
      use wrapper, only: matvec
      use, intrinsic :: iso_fortran_env
      implicit none
      integer,     intent(in)  :: nb,band(nb),ikpt,ispin
      MCOMPLEX_dp, intent(out) :: vxcmat(nb*(nb+1)/2)
      MCOMPLEX_dp              :: vpwxc(ngpt(ikpt)*(ngpt(ikpt)+1)/2),vpwhlp(ngpt(ikpt))
      complex_dp               :: vrmat(maxindx,maxindx),vrhlp(maxindx),cmt1(maxlmindx,maxval(neq),nb),cmt2(-maxlcut:maxlcut)
      real_dp                  :: integ(maxindx,maxindx,maxlh,0:maxlcut,0:maxlcut)
      real_dp                  :: intgrf,gaunt,gaunt1
      integer                  :: g(3),isym,ispin1,iband,jband
      integer                  :: i,j,ij,itype,ieq,ic,ic0,l,l1,l2,m,m1,m2,n,n1,n2,lm,lm1,lm2,ilh,jlh,ikindx
      if(storeibz.and.kptp(ikpt)/=ikpt) Error('K-point index out of range (storeibz).')
      vxcmat = 0
# ifdef LOAD
      allocate ( cmt(maxlmindx,ncent,nb,1,nspin3) )
      allocate ( cpw(maxgpt,         nb,1,nspin3) )
      call read_wavef2(band,size(band),ikpt,ispin,cmt,cpw)
      ikindx = 1
      ispin1 = 1
# else
      ikindx = kindx(ikpt)
      ispin1 = ispin
# endif
      ! MT
 1    integ = 0
      do itype = 1,ntype
        ! precalculate radial integrals
        do l2 = 0,lcut(itype)
          do l1 = 0,lcut(itype)
            do ilh = 1,nlh(itype)
              do n2 = 1,nindx(l2,itype)
                do n1 = 1,nindx(l1,itype)
                  integ(n1,n2,ilh,l1,l2) = intgrf( vmt_xc(:,ilh,itype,ispin) * (
     &                                             bas1(:,n1,l1,itype,ispin)*bas1(:,n2,l2,itype,ispin) +
     &                                             bas2(:,n1,l1,itype,ispin)*bas2(:,n2,l2,itype,ispin) ) , itype )
                enddo
              enddo
            enddo
          enddo
        enddo

        ! rotate cmt (->cmt1) (using a symmetry operation that relates the current atom to its parent; since the potential is locally symmetric, we can use any suitable symmetry operation isym.)
        ic0 = sum(neq(:itype-1)) + 1
        ic  = sum(neq(:itype-1))
        do ieq = 1,neq(itype)
          ic = ic + 1
          do isym = 1,nsym
            if(pcent(ic,isym)==ic0) exit
          enddo
          if(pcent(ic,isym)/=ic0)
     &      Bug('Did not find a suitable symmetry operation relating the current atom to its parent.')
          lm = 0
          do l = 0,lcut(itype)
            n1 = nindx(l,itype)
            do n = 1,n1
              do i = 1,nb ; iband = ifLoad(i,band(i))
                lm1                    = lm + n
                lm2                    = lm + n + n1*2*l
                if(l_soc) then
                  cmt2(-l:l)           = cmt(lm1:lm2:n1,ic,iband,ikindx,1) * sym(isym)%esoc(1,ispin1) +
     &                                   cmt(lm1:lm2:n1,ic,iband,ikindx,2) * sym(isym)%esoc(2,ispin1)
                else
                  cmt2(-l:l)           = cmt(lm1:lm2:n1,ic,iband,ikindx,ispin1)
                endif
                cmt1(lm1:lm2:n1,ieq,i) = matmul(dwgn(-l:l,-l:l,l,isym),cmt2(-l:l))!cmt(lm1:lm2:n1,ic,band(i),ikpt,ispin1))
              enddo
            enddo
            lm = lm + n1*(2*l+1)
          enddo
        enddo

        ! multiply cmt
        do l2 = 0,lcut(itype)
          do m2 = -l2,l2
            lm2 = sum( [ (nindx(l,itype)*(2*l+1),l=0,l2-1) ] ) + nindx(l2,itype) * (l2+m2)
            n2  = nindx(l2,itype)

            do l1 = 0,lcut(itype)
              do m1 = -l1,l1
                lm1 = sum( [ (nindx(l,itype)*(2*l+1),l=0,l1-1) ] ) + nindx(l1,itype) * (l1+m1)
                n1  = nindx(l1,itype)

                vrmat = 0
                do ilh = 1,nlh(itype)
                  l = llh(ilh,itype)
                  if(l<abs(l1-l2).or.l>l1+l2) cycle ! Gaunt condition
                  do jlh = 1,nmlh(ilh,itype)
                    m = mlh(jlh,ilh,itype)
                    if(m/=m1-m2)              cycle ! Gaunt condition
                    gaunt1         = gaunt(l,l1,l2,m,m1,m2) ! = INT Y_l1m1^* Y_lm Y_l2m2
                    vrmat(:n1,:n2) = vrmat(:n1,:n2) + gaunt1 * clh(jlh,ilh,itype) * integ(:n1,:n2,ilh,l1,l2)
                  enddo
                enddo

                ic = sum(neq(:itype-1))
                do ieq = 1,neq(itype)
                  ic = ic + 1
                  ij = 0
                  do j = 1,nb
                    vrhlp(:n1) = matmul( vrmat(:n1,:n2), cmt1(lm2+1:lm2+n2,ieq,j) )
                    do i = 1,j
                      ij         = ij + 1
                      vxcmat(ij) = vxcmat(ij) + dot_product ( cmt1(lm1+1:lm1+n1,ieq,i) , vrhlp(:n1) )
                    enddo
                  enddo
                enddo

              enddo
            enddo
          enddo
        enddo
      enddo

      ! PW
      ij = 0
      do j = 1,ngpt(ikpt)
        do i = 1,j
          g         = gpt(:,pgpt(i,ikpt)) - gpt(:,pgpt(j,ikpt))
          ij        = ij + 1
          vpwxc(ij) = vpw_xc(g(1),g(2),g(3),ispin)
        enddo
      enddo
      ij = 0
      do j = 1,nb ; jband = ifLoad(j,band(j))
        vpwhlp = matvec( vpwxc, cpw(:ngpt(ikpt),jband,ikindx,ispin1) )
        do i = 1,j ; iband = ifLoad(i,band(i))
          ij         = ij + 1
          vxcmat(ij) = vxcmat(ij) + dot_product ( cpw(:ngpt(ikpt),iband,ikindx,ispin1) , vpwhlp )
        enddo
      enddo
      if(l_soc.and.ispin1==1) then
        ispin1 = 2
        goto 1
      endif
      Load( deallocate(cmt,cpw) )
      end

c --------------------------------------

      subroutine plussoc
      use global
      use file
      use wrapper
      use Hwrapper
      use readwrite
      use hdf5
      use, intrinsic :: iso_fortran_env
      implicit none
      real_dp                 :: kpt_new(3,nkpt)
      integer                 :: kptp_new(nkpt),symkpt_new(nkpt),nkpti_new
      real_dp                 :: rot(3,3),kvec(3)
      complex_dp, allocatable :: cmt1(:,:),cmt2(:,:),cmt3(:,:),cmt4(:,:)
      complex_dp, allocatable :: hamiltonian(:,:),eigv(:,:),cpw_1(:,:,:),cmt_1(:,:,:,:)
      complex_dp, allocatable :: dwgn1(:,:,:)
      real_dp,    allocatable :: eig(:)
      complex_dp              :: cdum
      integer                 :: isym,n,sym1(nsym),iunit
      integer                 :: i,j,ikpt,ivec(3),ikptp
      integer                 :: ispin,itype,ieq,ic,l,lm,lm1,m
      integer                 :: s1,s2,ms1,ms2,nb1,n1,n2
      real_dp                 :: dvmt(maxgrid),dvmt_so(maxgrid),integr(maxindx,maxindx),rdum,sqaxis(3)
      real_dp                 :: determinant,intgrf

      beginSingle

      Load(       Error('PLUSSOC not implemented for -DLOAD.') )
      if(l_soc)   Error('SOC already included. Do not use PLUSSOC.')
      if(lkptadd) Error('not implemented yet for additional k point')

      sqaxis = [ sin(sqa(1)) * cos(sqa(2)) , sin(sqa(1)) * sin(sqa(2)) , cos(sqa(1)) ]
      nspin3 = 2 ! write out spinors

      write(6,'(/A)')      'Second variation for spin-orbit coupling.'
      write(6,'(A,3F8.4)') '  Spin quantization axis:',sin(sqa(1))*cos(sqa(2)),sin(sqa(1))*sin(sqa(2)),cos(sqa(1))

      if(nspin==2) then
        write(6,'(A)') '  Non-collinear case'
        iunit = fopen('spex.sym',status='unknown')
        n     = 0
        do isym = 1,nsym
          if(isym>nsymt) cycle
          rot = matmul(lat,matmul(sym(isym)%rot,transpose(rlat)))/(2*pi)
          rot = rot / determinant(rot)
          if(any(abs(matmul(rot,sqaxis)-sqaxis)>1d-12)) cycle
          n       = n + 1
          sym1(n) = isym
        enddo
        do isym = 1,n
          write(iunit,'(''!'',I3)') sym1(isym)
          write(iunit,'(3I5,F14.5)') (sym(sym1(isym))%rot(i,:),sym(sym1(isym))%transl(i),i=1,3)
        enddo
        call fclose(iunit)
        write(6,'(A)') '  New symmetry file for non-collinear case was written to spex.sym.'
        call ibz(kpt_new,kptp_new,symkpt_new,nkpti_new,[0d0,0d0,0d0],nspin==2)
        write(6,'(A,I3,A)') '  Corresponding IBZ has',nkpti_new,' k points:'
        do ikpt = 1,nkpti_new
          write(6,'(''  ('',F7.5,'','',F7.5,'','',F7.5,'')'')') kpt_new(:,ikpt)
        enddo
        if(storeibz.and.nkpti_new>nkpti) Error('PLUSSOC + spinpol. implemented only with STOREBZ.')
      else
        nkpti_new  = nkpti
        kpt_new    = kpt(:,:nkpt)
        kptp_new   = kptp(:nkpt)
        symkpt_new = symkpt(:nkpt)
      endif

      ! Rewrite parameter file
      call write_param(.true.)

      ! Define dwgn1 for rotation if sqaxis != zaxis
      if(any(sqa/=0)) then
        allocate ( dwgn1(-maxlcut:maxlcut,-maxlcut:maxlcut,0:maxlcut) )
        call rotation_sph(rot,sqa(1),sqa(2))
        call dwigner(dwgn1,rot,maxlcut)
      endif

      ! Second variation and write data to file

      do ikpt = 1,nkpti_new
        kvec  = kpt_new(:,ikpt)
        ivec  = nint(kvec*nkpt3)
        ikptp = pkpt(ivec(1),ivec(2),ivec(3))
        if(storeibz.and.ikptp>nkpti) Bug('ikptp > nkpti.')
        n = sum(nband(ikptp,:)) * 2/nspin
        allocate ( hamiltonian(n,n),eigv(n,n),eig(n) )
        hamiltonian = 0
        ! SOC term
        nb1 = nband(ikptp,1)
        n   = maxval(nband(ikptp,:))
        allocate ( cmt1(maxindx,n),cmt2(maxindx,n),cmt3(maxindx,n),cmt4(maxindx,n) )
        do itype = 1,ntype
          call derivative(dvmt,vmt(:,1,itype,1),itype)
          do i = 1,grid(itype)%number
            dvmt(i) = dvmt(i) / rgrid(i,itype)
          enddo
          do l = 1,lcut(itype)
            ! radial integral INT r² u1(r) u2(r) 1/(2Mc)² r² dV/dr dr
            integr = 0
            rdum   = 2 * clight**2
            do i = 1,nindx(l,itype)
              do j = 1,nindx(l,itype)
              dvmt_so     = dvmt * ( 0.5d0 / (1+(ebas(i,l,itype,1)-vmt(:,1,itype,1))/rdum)**2
     &                             + 0.5d0 / (1+(ebas(j,l,itype,1)-vmt(:,1,itype,1))/rdum)**2 ) / (2*rdum)
              integr(i,j) = intgrf( dvmt_so * ( bas1(:,i,l,itype,1) * bas1(:,j,l,itype,1) ) , itype )
              enddo
            enddo
            n = nindx(l,itype)

            ! perform L*S and multply with cmt
            ic = sum(neq(:itype-1))
            do ieq = 1,neq(itype)
              ic = ic + 1

              do s1 = 1,2 ; ms1 = min(s1,nspin)
              do s2 = 1,2 ; ms2 = min(s2,nspin)

              lm = sum( [ ((2*i+1)*nindx(i,itype),i=0,l-1) ] )

              do m = -l,l
                do n1 = 1,maxval(nband(ikptp,:))
                  if(m/=-l) then ; lm1 = lm - n  ; cmt1(:n,n1) = cmt_trafo(lm1,m-1,ms1) ; endif
                                   lm1 = lm      ; cmt2(:n,n1) = cmt_trafo(lm1,m,  ms1)
                  if(m/=l)  then ; lm1 = lm + n  ; cmt3(:n,n1) = cmt_trafo(lm1,m+1,ms1) ; endif
                                   lm1 = lm      ; cmt4(:n,n1) = cmt_trafo(lm1,m,  ms2)
                                                   cmt4(:n,n1) = matmul ( integr(:n,:n) , cmt4(:n,n1) )
                enddo
                do n1 = 1,nband(ikptp,ms1)
                  do n2 = 1,nband(ikptp,ms2)
                    if     (s1==1.and.s2==1) then ! up/up
                      cdum                       = dot_product ( cmt2(:n,n1),cmt4(:n,n2) )
                      hamiltonian(n1,n2)         = hamiltonian(n1,n2)         + m * cdum
                    else if(s1==2.and.s2==2) then ! down/down
                      cdum                       = dot_product ( cmt2(:n,n1),cmt4(:n,n2) )
                      hamiltonian(nb1+n1,nb1+n2) = hamiltonian(nb1+n1,nb1+n2) - m * cdum
                    else if(s1==1.and.s2==2) then ! up/down
                      if(m/=-l) then
                        hamiltonian(n1,nb1+n2)   = hamiltonian(n1,nb1+n2)     + sqrt( l*(l+1d0) - m*(m-1) ) *
     &                                             dot_product ( cmt1(:n,n1),cmt4(:n,n2) )
                      endif
                    else if(s1==2.and.s2==1) then ! down/up
                      if(m/=l) then
                        hamiltonian(nb1+n1,n2)   = hamiltonian(nb1+n1,n2)     + sqrt( l*(l+1d0) - m*(m+1) ) *
     &                                             dot_product ( cmt3(:n,n1),cmt4(:n,n2) )
                      endif
                    endif
                  enddo
                enddo
                lm = lm + nindx(l,itype)

              enddo

              enddo ! s1
              enddo ! s2

            enddo

          enddo

        enddo

        deallocate ( cmt1,cmt2,cmt3,cmt4 )

        ! diagonal : eigenvalues
        do ispin = 1,2
          n = sum(nband(ikptp,:ispin-1))
          do i = 1,nband(ikptp,min(ispin,nspin))
            hamiltonian(n+i,n+i) = hamiltonian(n+i,n+i) + ene(i,ikptp,min(ispin,nspin))
          enddo
        enddo

        call diagonalize(eigv,eig,hamiltonian)

        write(6,'(/A,I5,A,I5)') 'K point',ikpt,' ->',ikptp
        do i = 1,nkpt
          if(kptp_new(i)==ikpt) write(6,'(I4,3F10.5)') i,kpt(:,i)
        enddo
        write(6,'(5F12.6)') eig

        n = size(eig)

        allocate ( cmt_1(size(cmt,1),size(cmt,2),n,2) )
        allocate ( cpw_1(size(cpw,1),            n,2) )

        cpw_1(:ngpt(ikptp),:,1) = matmul(cpw(:ngpt(ikptp),:nband(ikptp,1),    ikptp,1),    eigv(:nband(ikptp,1)   ,:))
        cpw_1(:ngpt(ikptp),:,2) = matmul(cpw(:ngpt(ikptp),:nband(ikptp,nspin),ikptp,nspin),eigv( nband(ikptp,1)+1:,:))

        do ic = 1,ncent
          cmt_1(:,ic,:,1) = matmul(cmt(:,ic,:nband(ikptp,1),    ikptp,1),    eigv(:nband(ikptp,1)   ,:))
          cmt_1(:,ic,:,2) = matmul(cmt(:,ic,:nband(ikptp,nspin),ikptp,nspin),eigv( nband(ikptp,1)+1:,:))
        enddo

        call write_wavef(ikpt,1,ngpt(ikptp),gpt(:,pgpt(:ngpt(ikptp),ikptp)),eig,cmt_1,cpwin_c=cpw_1) ! spin index is always 1

        deallocate ( cmt_1,cpw_1 )
        deallocate ( hamiltonian,eigv,eig )

      enddo

      if(any(sqa/=0)) deallocate(dwgn1)

      write(6,'(/A)') 'New input data written.'
      if(nspin==2) write(6,'(A)') 'New symmetry file written to spex.sym.'
# ifdef noHDF
      Info('Input data can only be read by Spex when compiled without -DnoHDF.')
# endif

      endSingle

      call finish

      contains

      ! returns cmt coefficients in the coordinate system with sqaxis as the z axis
      ! (dwgn1 is defined for rotation z->sqaxis; therefore, its inverse is used)
      function cmt_trafo(off,m,ispin)
      implicit none
      integer, intent(in) :: off,m,ispin
      complex_dp          :: cmt_trafo(n)
      integer             :: i
      if(all(sqa==0)) then
        cmt_trafo = cmt(off+1:off+n,ic,n1,ikptp,ispin)
      else
        do i = 1,n
          cmt_trafo(i) = dot_product ( dwgn1(-l:l,m,l) , cmt(off-(l+m)*n+i:off+(l-m)*n+i:n,ic,n1,ikptp,ispin) )
        enddo
      endif
      end function cmt_trafo

      end

c --------------------------------------

c     Returns list of q vectors in qvec(:,:nq) as obtained from KPTPATH or a file (qfile includes quotes "...")
c
c begin interface
      subroutine getqlist(qvec,nq,mkpt_def,qfile)
      use global
      use key
      use file
      use util
      use, intrinsic :: iso_fortran_env !inc
      implicit none
      integer, intent(in)               :: mkpt_def
      integer, intent(out)              :: nq
      real_dp, intent(out), allocatable :: qvec(:,:)
      character(*), intent(in)          :: qfile
c end interface
      real_dp                           :: qfac,qmin,qlen
      integer                           :: ikpt,i,n,iq
      integer                           :: iunit,istat
      character(80)                     :: line
      if(qfile==' ') then
        if(nkpt_path==0) Error('KPTPATH definition missing.')
        if(mkpt_path==0) then ; qmin = rvol**(1d0/3) / mkpt_def  ! shortest step size (default)
        else                  ; qmin = rvol**(1d0/3) / mkpt_path ! shortest step size (user defined)
        endif
        do
          nq = 0
          do ikpt = 1,nkpt_path-1
            qlen = sqrt(sum(matmul(rlat,kpt_path(:,ikpt+1)-kpt_path(:,ikpt))**2))
            n    = max(1,nint(qlen/qmin))
            do i = 0,n ; if(i==0.and.ikpt/=1) cycle
              nq = nq + 1              
              if(allocated(qvec)) qvec(:,nq) = ( kpt_path(:,ikpt)*(n-i) + kpt_path(:,ikpt+1)*i ) / n
            enddo
          enddo
          if(.not.allocated(qvec)) then ; allocate(qvec(3,nq))
          else                          ; exit
          endif
        enddo
      else
        if(qfile(:1)//qfile(len_trim(qfile):)/='""') Error('Filenames must be given in quotes "..." (qpoint file)')
        iunit = fopen(qfile(2:len_trim(qfile)-1),status='old')
        read(iunit,*,iostat=istat) nq,qfac ; if(istat/=0) Error('Read error: '//trim(qfile)//', line: 1')
        allocate(qvec(3,nq))
        do iq = 1,nq
          read(iunit,'(A)',iostat=istat) line  ; if(istat/=0) Error('Read error: '//trim(qfile))
          read(line,*,iostat=istat) qvec(:,iq) ; if(istat/=0) Error('Read error: '//trim(qfile)//', line: '//chr(iq+1))
          qvec(:,iq) = qvec(:,iq) / qfac
        enddo
        call fclose(iunit)
      endif
      end

c --------------------------------------

c     Writes header of "bands*" file and also returns qpath(0:nq) (k axis of band structure).
      subroutine bands_header(iunit,title,qpath,qvec,nq)
      use global, only: nkpt_path,kpt_path,nkpt_label,kpt_label,ckpt_label,rlat
      use, intrinsic :: iso_fortran_env
      implicit none
      integer,      intent(in)  :: iunit,nq
      character(*), intent(in)  :: title
      real_dp,      intent(out) :: qpath(0:nq)
      real_dp,      intent(in)  :: qvec(3,nq)
      real_dp                   :: qvec0(3)
      integer                   :: ikpt,iq,i
      write(iunit,'(A)') title
      ikpt     = 1
      qpath(0) = 0
      qvec0    = qvec(:,1)
      do iq = 1,nq
        qpath(iq) = qpath(iq-1) + sqrt(sum(matmul(rlat,qvec(:,iq)-qvec0)**2))
        if(nkpt_path>0) then
          if(all(abs(qvec(:,iq)-kpt_path(:,ikpt))<1d-12)) then
            write(iunit,'(''#'',F15.10,3F10.5,''  '''NoA) qpath(iq),qvec(:,iq)
            do i = 1,nkpt_label
              if(all(abs(kpt_label(:,i)-kpt_path(:,ikpt))<1d-12)) write(iunit,'(A'NoA) ' ('//trim(ckpt_label(i))//')'
            enddo
            write(iunit,*)
            ikpt = ikpt + 1
          endif
        endif
        qvec0 = qvec(:,iq)
      enddo
      end

c --------------------------------------

# endif

