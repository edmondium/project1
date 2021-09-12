c--------------------------------------------------------------------------------
c Copyright (c) 2021 Peter Grünberg Institut, Forschungszentrum Jülich, Germany
c This file is part of SPEX and available as free software under the conditions
c of the MIT license as expressed in the LICENSE file in more detail.
c--------------------------------------------------------------------------------

c     Defines the mixed basis:
c
c     1) products of interstitial plane waves up to a cutoff gcutm
c     2) MT functions in the form U_il(r)*Y_lm(r) where the radial functions U_il(r)
c        are linear combinations of basis-function products. The corresponding
c        coefficients are the components of eigenvectors obtained from
c        the diagonalization of the overlap matrix of basis-function products.
c        Eigenvectors with eigenvalues below a given parameter (tolerance) are
c        omitted, because they do not contribute (much) to the flexibility of the
c        basis.
c
c     Explicit input:
c       Mode 1 = | val  * val  >
c       Mode 2 = | core * val  > + | val  * val  >
c       Mode 3 = | core * core > + | val  * val  >
c       Mode 4 = | core * val  > + | core * core > + | val  * val >
c
c     Selected output (global variables and arrays):
c       ngptmall         = number of unique G points
c       gptm(:,:)        = unique G vectors
c       ngptm(:)         = number of G vector at k points
c       maxgptm          = maximal ngptm
c       pgptm(:,:)       = points to G vectors for given k point (list of G vectors)
c       nbasp            = number of mixed-basis MT functions
c       nbasm(:)         = number of mixed-basis functions at k points (=nbasp+ngptm)
c       maxbasm          = maximal nbasm
c       nindxm(:,:)      = number of radial functions for given lm
c       basm(:,:,:,:)    = radial part of MT functions ( r*U_il(r) )

# include "cppmacro.h"
# include "jobtype.h"
# include "restype.h"
      
      subroutine mixedbasis(mode)

      use global
      use util
      use wrapper
      use readwrite
      use key
      Mpi( use Mwrapper )

      use, intrinsic :: iso_fortran_env
      implicit none
      integer,       intent(in)  :: mode
      integer                    :: ikpt,g(3)
      logical                    :: selecmat(max(maxindx,maxindxc),0:(max(maxlcut,maxlcutc)),
     &                                       max(maxindx,maxindxc),0:(max(maxlcut,maxlcutc)))
      integer                    :: spin,ispin,ispin1,ispin2,itype,ieq,icent,iband,
     &                              l1,l2,l,n,n1,n2,nn,i,j,m,m1,m2,lm1_0,lm2_0,lm,lm1,lm2,ll,ng,no,nu,x,y,z
      integer,       allocatable :: ihelp(:),nindxp(:,:,:)
      real_dp,       allocatable :: olap(:,:),olap2(:,:),work(:),eig(:),eigv(:,:),error(:,:),norm(:,:),
     &                              basprod(:,:),basp(:,:),dbasp(:),basm1(:,:,:,:),dbasm1(:,:,:),prodm(:,:)
c      MCOMPLEX_dp,   allocatable :: olappw(:,:),cprod(:,:,:)
      complex_dp,    allocatable :: chelp1(:,:,:),chelp2(:,:,:)
      complex_dp                 :: cdum
      real_dp                    :: rdum,rdum1,kvec(3),mp1,mp2
      logical                    :: offdiag,ldum,ldum1,first
      logical,       allocatable :: seleco(:,:,:),selecu(:,:,:)
      character(8),  parameter   :: chmode(4) = [ 'vv      ','cv+vv   ','cc+vv   ','cc+cv+vv' ]
      character(80)              :: parse_out(6),line
      character(:),  allocatable :: charr(:)
      logical                    :: chkprod,addbas
      real_dp                    :: intgrf,gaunt
      real                       :: cputime
      MCOMPLEX_dp                :: stepfunction
      integer                    :: ch2i
      type prodtype
        integer                  :: l1,l2,n1,n2
      end type prodtype
      type(prodtype), allocatable :: prod(:)

c      kpt(:,nkpti) = [1d-4,0d-4,0d-4]/sqrt(1d0) ; write(*,*) 'SMALL KPT APPENDED!!!'

      Rbegin

      call cpu_time(cputime)

      write(6,'(//A,A,A)') '### subroutine: mixedbasis ( Mode: ',trim(chmode(mode)),' ) ###'

      if(all(mode/=[1,2,3,4]))                     Error('Mode is unknown.')
      if((mode==2.or.mode==4).and..not.any(cores)) Error('No cores specified.')

c     Read spex.inp
      call getkey(inp,'TOL',     tolerance, section='MBASIS', default=0.0001d0, mini=0d0, writeout=.false.)
      call getkey(inp,'ADDBAS',     addbas, section='MBASIS', default=.false.,            writeout=.false.)
      call getkey(inp,'CHKPROD',   chkprod, section='MBASIS', default=.false.,            writeout=.false.)
      Mpi ( if(chkprod) then ; write(6,'(A)') 'MPI: CHKPROD disabled.' ; chkprod = .false. ; endif )
      if(any(job(:)%type==[J_COSX]).and..not.addbas) Info('It is recommendable to use ADDBAS for COSX calculations.')
      if(maxindx>0) then
        allocate( seleco(maxindx,0:maxlcut,ntype),selecu(maxindx,0:maxlcut,ntype) )
        allocate( character(80) ::  charr(ntype) )
        call getkey(inp,'SELECT',    charr, section='MBASIS', status=i,                   writeout=.false.)
        if(i==1) Error('Missing arguments after SELECT.')
        if(i==0) then
          do itype = 1,ntype
            write(charr(itype),'(I1,A,I1)') lcutm(itype)/2,';',lcutm(itype)/2+1 ! define default
          enddo
        endif
        allocate ( ihelp(0:max(4,maxlcut)) )
        seleco = .false.
        selecu = .false.
        do itype = 1,ntype
          call parser(parse_out,charr(itype),'1<,2><,5>;3<,4><,6>','SELECT')
          do l = 0,lcut(itype)
            do i = 3,nindx(l,itype)
              seleco(i,l,itype) = any(ebas(i,l,itype,:)<ebas(1,l,itype,:)) ! The default is to include semicore LOs
              selecu(i,l,itype) = any(ebas(i,l,itype,:)<ebas(1,l,itype,:)) ! but not LOs at higher energies.
            enddo
          enddo
          if(parse_out(2)(:1)==',') then
            if(parse_out(2)==','.and.maxindx>=3) seleco(3:,:,itype) = .false. ! empty LO definition -> exclude LOs
            parse_out(5) = parse_out(2)(2:)
            parse_out(2) = ' '
          endif
          if(parse_out(4)(:1)==',') then
            if(parse_out(4)==','.and.maxindx>=3) selecu(3:,:,itype) = .false.
            parse_out(6) = parse_out(4)(2:)
            parse_out(4) = ' '
          endif
          ihelp    = 0
          ihelp(1) = ch2i(parse_out(1),'SELECT') ; if(parse_out(2)/=' ') ihelp(2) = ch2i(parse_out(2),'SELECT')
          ihelp(3) = ch2i(parse_out(3),'SELECT') ; if(parse_out(4)/=' ') ihelp(4) = ch2i(parse_out(4),'SELECT')
          if(any(ihelp<0))           Error('values after keyword SELECT must not be negative.')
          if(any(ihelp>lcut(itype))) Error('values after keyword SELECT must not be larger than l cutoff.')
          seleco(1,:ihelp(1),itype) = .true. ; if(parse_out(2)/=' ') seleco(2,:ihelp(2),itype) = .true.
          selecu(1,:ihelp(3),itype) = .true. ; if(parse_out(4)/=' ') selecu(2,:ihelp(4),itype) = .true.
          if(parse_out(5)/=' ') then
            if(maxindx<=2) Error('LOs SELECTed but no LOs defined in DFT input.')
            call lodef(line,parse_out(5))
            if(len_trim(line)/=nlo(itype)) Error('Incorrect number of LOs after SELECT.')
            do i = 1,len_trim(line)
              j = count(llo(:i,itype)==llo(i,itype))
              if     (line(i:i)=='0') then ; seleco(2+j,llo(i,itype),itype) = .false.
              else if(line(i:i)/='1') then ; Error('expected series of 0 and 1 after SELECT (LO def)')
              endif
            enddo
          endif
        enddo
        deallocate ( charr,ihelp )
      endif

      call getkey(inp,'WFADJUST', ldum, section='MBASIS',  default=.false.)
      if(ldum) then
# if defined(LOAD) || defined(MPI)
        Error('WFADJUST disabled for macros LOAD or MPI.')
# endif
        Warn('WFADJUST has been rewritten and is untested.')
        maxlmindx = 0
        icent     = 0
        do itype = 1,ntype
          do ieq = 1,neq(itype)
            icent = icent + 1
            lm  = 0
            lm1 = 0
            do l = 0,lcut(itype)
              do m = -l,l
                do n = 1,nindx(l,itype)
                  lm = lm + 1
                  if(seleco(n,l,itype).or.selecu(n,l,itype)) then
                    lm1                  = lm1 + 1
                    cmt(lm1,icent,:,:,:) = cmt(lm,icent,:,:,:)
                  endif
                enddo
              enddo
            enddo
          enddo
          maxlmindx = max(maxlmindx,lm1)
        enddo
        do itype = 1,ntype
          ll = -1
          do l = 0,lcut(itype)
            nn = 0
            do n = 1,nindx(l,itype)
              if(seleco(n,l,itype).or.selecu(n,l,itype)) then
                nn                   = nn + 1
                ll                   = max(ll,l)
                bas1(:,nn,l,itype,:) = bas1(:,n,l,itype,:)
                bas2(:,nn,l,itype,:) = bas2(:,n,l,itype,:)
              endif
            enddo
            nindx(l,itype) = nn
          enddo
          lcut(itype) = ll
        enddo
        maxindx = maxval(nindx)
      endif

c
c     Define G-point set gptm for mixed basis with cutoff gcutm and corresponding pointer pgptm

      allocate ( ngptm(nkpt2),nbasm(nkpt2) )
c      i = 4*pi/3*(1.1*gcutm)**3 * vol/(8*pi**3) * nkpt2 ! 2
c      write(*,*) i,nkpt2,i*nkpt2
      first = .true.
 1    ngptm = 0
      i     = 0
      n     = -1
      kvec  = [ 0d0,0d0,0d0 ] ! workaround for the xlf compiler
      rdum1 = maxval ( [ (sqrt(sum(matmul(rlat,kpt(:,ikpt)+kvec)**2)),ikpt=1,nkpt2) ] )
      do
        n    = n + 1
        ldum = .false.
        do x = -n,n
          n1 = n-abs(x)
          do y = -n1,n1
            n2 = n1-abs(y)
            do z = -n2,n2,max(2*n2,1)
              g     = [x,y,z]
              rdum  = sqrt(sum(matmul(rlat,g)**2))-rdum1 ; if(rdum>gcutm) cycle
              ldum1 = .false.
              do ikpt = 1,nkpt2
                kvec = kpt(:,ikpt)
                rdum = sum(matmul(rlat,kvec+g)**2)
                if(rdum<=gcutm**2) then
                  if(.not.ldum1) then
                    i     = i + 1 !; if(i>size(gptm,2)) call reallocate ( gptm,3,i+100 )
                    ldum1 = .true.
                    if(.not.first) gptm(:,i) = g
                  endif
                  !if(ngptm(ikpt)+1>size(pgptm,1)) call reallocate ( pgptm,ngptm(ikpt)+100,nkpt2 )
                  ngptm(ikpt) = ngptm(ikpt) + 1
                  ldum        = .true.
                  if(.not.first) pgptm(ngptm(ikpt),ikpt) = i
                endif
              enddo
            enddo
          enddo
        enddo
        if(.not.ldum) exit
      enddo
      if(first) then
        first    = .false.
        ngptmall = i
        maxgptm  = maxval(ngptm)
        allocate ( gptm(3,ngptmall) )
        allocate ( pgptm(maxgptm,nkpt2) )
        goto 1
      endif
      allocate ( work(maxgptm),ihelp(maxgptm) )
      do ikpt = 1,nkpt2
        do i = 1,ngptm(ikpt)
          work(i) = sqrt(sum(matmul(rlat,gptm(:,pgptm(i,ikpt)))**2)) + i*1d-12 ! "+i*1d-12" prevents rounding errors from changing the order
        enddo
        call rorderp(ihelp,work,ngptm(ikpt))
        pgptm(:ngptm(ikpt),ikpt) = pgptm(ihelp(:ngptm(ikpt)),ikpt)
      enddo
      deallocate ( work,ihelp )

      write(6,'(/A)')   'Mixed basis'
      write(6,'(A)')    '  Plane-wave basis:'
      write(6,'(A,I6)') '    Number of unique G vectors:',ngptmall
      write(6,'(A,I5)') '    Maximal number of G vectors:',maxgptm

      g = [ (maxval(abs(gptm(i,:))),i=1,3) ]
      allocate ( pntgptm(-g(1):g(1),-g(2):g(2),-g(3):g(3),nkpt2) )
      pntgptm = 0
      do ikpt = 1,nkpt2
        do i = 1,ngptm(ikpt)
          j                            = pgptm(i,ikpt)
          g                            = gptm(:,j)
          pntgptm(g(1),g(2),g(3),ikpt) = i
        enddo
      enddo

      ! Precalculate Fourier series of step function (cstep), will be redefined in mixedbasis
      g = 2 * [ (maxval(abs(gpt(i,:))),i=1,3) ] + [ (maxval(abs(gptm(i,:))),i=1,3) ] + 1
      if(allocated(cstep)) tDeallocate(cstep)
      Allocate_ ( cstep,(-g(1):g(1),-g(2):g(2),-g(3):g(3)) )
      do x = -g(1),g(1)
      do y = -g(2),g(2)
      do z = -g(3),g(3)
        cstep(x,y,z) = stepfunction([x,y,z])
      enddo
      enddo
      enddo

      if(iand(restart,R_mbas)/=0) then
        call read_mixedbasis(ldum)
        if(ldum) goto 3
      endif

c
c     Set up MT product basis

      write(6,'(A)') '  MT product basis:'
      write(6,'(A)') '    Reduction due to overlap (quality of orthonormality, should be < 1.0E-06)'
      allocate (nindxm(0:maxlcutm,ntype),basm1(maxgrid,1,0:maxlcutm,ntype),dbasm1(1,0:maxlcutm,ntype))
      nindxm = 0
      basm1  = 0
      do itype = 1,ntype
        ng = grid(itype)%number
        do l = 0,lcutm(itype)
          first = .true.
 2        n     = 0
          m     = 0
          do ispin = 1,nspin ! we use products | up up > and | down down >
            !  products  | val * val >
            selecmat = .true.
            do l1 = 0,lcut(itype)
              do n1 = 1,nindx(l1,itype)
                rdum = intgrf(bas1(:ng,n1,l1,itype,ispin)**2+bas2(:ng,n1,l1,itype,ispin)**2,itype)
                do l2 = 0,lcut(itype)
                  do n2 = 1,nindx(l2,itype)
                    if(l<abs(l1-l2).or.l>l1+l2)    cycle ! Fulfill triangular inequality
                    if(mod(l+l1+l2,2)/=0)          cycle ! Use that the Gaunt coeff. vanish for l+l1+l2 = odd
                    if(.not.selecmat(n1,l1,n2,l2)) cycle ! Prevent double counting of products a*b = b*a
                    m = m + 1                               ! m not used at the moment
                    if(.not.seleco(n1,l1,itype))   cycle ! Basis-function selection in spex.inp (occ.)
                    if(.not.selecu(n2,l2,itype))   cycle ! Basis-function selection in spex.inp (unocc.)
                    n = n + 1 ; selecmat(n2,l2,n1,l1) = .false.
                    !if(n>size(basm,2)) call reallocate(basm,maxgrid,n,maxlcutm,ntype)
                    if(.not.first) then
                      basp(:ng,n) = ( bas1(:ng,n1,l1,itype,ispin)*bas1(:ng,n2,l2,itype,ispin)   ! Calculate those basis-function products
     &                              + bas2(:ng,n1,l1,itype,ispin)*bas2(:ng,n2,l2,itype,ispin) ) ! from which the product basis will be constructed.
     &                              / rgrid(:ng,itype)                                          ! Note that we do not need the more general two-component product functions!
                      dbasp(n)    =  ubas(n1,l1,itype,ispin) * dubas(n2,l2,itype,ispin) +
     &                              dubas(n1,l1,itype,ispin) *  ubas(n2,l2,itype,ispin)
                      rdum1       = sqrt(rdum*intgrf(bas1(:ng,n2,l2,itype,ispin)**2+bas2(:ng,n2,l2,itype,ispin)**2,itype))
                      basp(:ng,n) = basp(:ng,n) / rdum1
                      dbasp(n)    = dbasp(n)    / rdum1
                    endif
                  enddo
                enddo
              enddo
            enddo
            if(mode==2.or.mode==4) then
              !  products  | core * val >
              do l1 = 0,lcutc(itype)
                do n1 = 1,nindxc(l1,itype)
                  if(.not.cores(n1,l1,itype)) cycle
                  do l2 = 0,lcut(itype)
                    do n2 = 1,nindx(l2,itype)
                      if(.not.selecu(n2,l2,itype)) cycle ! Basis-function selection in spex.inp (unocc.)
                      if(l<abs(l1-l2).or.l>l1+l2)  cycle ! Fulfill triangular inequality
                      if(mod(l+l1+l2,2)/=0)        cycle ! Use that the Gaunt coeff. vanish for l+l1+l2 = odd
                      n = n + 1
                      !if(n>size(basm,2)) call reallocate(basm,maxgrid,n,maxlcutm,ntype)
                      if(.not.first) then
                        basp(:ng,n) = ( core1(:ng,n1,l1,itype,ispin)*bas1(:ng,n2,l2,itype,ispin)   ! Calculate those basis-function products
     &                                + core2(:ng,n1,l1,itype,ispin)*bas2(:ng,n2,l2,itype,ispin) ) ! from which the product basis will be constructed.
     &                                / rgrid(:ng,itype)                                           ! Note that we do not need the more general two-component product functions!
                        dbasp(n)    = 0
                        rdum        = sqrt(intgrf(bas1(:ng,n2,l2,itype,ispin)**2+bas2(:ng,n2,l2,itype,ispin)**2,itype))
                        basp(:ng,n) = basp(:ng,n) / rdum
                      endif
                    enddo
                  enddo
                enddo
              enddo
            endif
            if(mode==3.or.mode==4) then
              !  products  | core * core >
              selecmat = .true.
              do l1 = 0,lcutc(itype)
                do n1 = 1,nindxc(l1,itype)
                  do l2 = 0,lcutc(itype)
                    do n2 = 1,nindxc(l2,itype)
                      if(l<abs(l1-l2).or.l>l1+l2)    cycle ! Fulfill triangular inequality
                      if(mod(l+l1+l2,2)/=0)          cycle ! Use that the Gaunt coeff. vanish for l+l1+l2 = odd
                      if(.not.selecmat(n1,l1,n2,l2)) cycle ! Prevent double counting of products a*b = b*a
                      selecmat(n2,l2,n1,l1) = .false.
                      n = n + 1
                      !if(n>size(basm,2)) call reallocate(basm,maxgrid,n,maxlcutm,ntype)
                      if(.not.first) then
                        basp(:ng,n) = ( core1(:ng,n1,l1,itype,ispin)*core1(:ng,n2,l2,itype,ispin)   ! Calculate those basis-function products
     &                                + core2(:ng,n1,l1,itype,ispin)*core2(:ng,n2,l2,itype,ispin) ) ! from which the product basis will be constructed.
     &                                / rgrid(:ng,itype)                                            ! Note that we do not need the more general two-component product functions!
                        dbasp(n)    = 0
                      endif
                    enddo
                  enddo
                enddo
              enddo
            endif
            ! add bas1 if required
            if(addbas) then
              do i = 1,nindx(l,itype)
                if(seleco(i,l,itype)) then
                  n = n + 1
                  if(.not.first) then
                    basp(:ng,n) = bas1(:ng,i,l,itype,ispin)
                    dbasp(n)    = dubas(i,l,itype,ispin)
                  endif
                endif
              enddo
            endif
          enddo
          ! stop if there aren't any product-basis functions
          if(n==0)
     &      Error('No MT products found for l='//Chr(l)//' and atom type='//Chr(itype)//'. Check LCUT and SELECT (MBASIS).')
          ! cycle for second loop
          if(first) then
            first = .false.
            allocate ( basp(maxgrid,n),dbasp(n) )
            goto 2
          endif
          ! normalize product functions
          do i = 1,n
            rdum        = sqrt ( intgrf(basp(:ng,i)**2,itype) )
            basp(:ng,i) = basp(:ng,i) / rdum
            dbasp(i)    = dbasp(i) / rdum
          enddo
          ! calculate overlap
          allocate ( olap(n,n),eigv(n,n),eig(n),ihelp(n) )
          olap = 0
          do n2 = 1,n
            do n1 = 1,n2
              olap(n1,n2) = intgrf(basp(:ng,n1)*basp(:ng,n2),itype)
              olap(n2,n1) = olap(n1,n2)
            enddo
          enddo
          ! diagonalize
          call diagonalize(eigv,eig,packmat(olap))
          ! get rid of linear dependencies (eigenvalue <= tolerance)
          nn = 0
          do i = n,1,-1
            if(eig(i)>tolerance) then
              nn        = nn + 1
              ihelp(nn) = i
            endif
          enddo
          nindxm(l,itype) = nn
          eig(:nn)        = eig(ihelp(:nn))
          eigv(:,:nn)     = eigv(:,ihelp(:nn))
          ! Loewdin orthonormalization of product basis (start with constant function for l=0)
          m = 0
          if(l==0) then
            basm1(:ng,1,0,itype) = rgrid(:ng,itype) / sqrt(grid(itype)%radius**3/3)
            dbasm1(   1,0,itype) = 0
            m = 1
          endif
          if(size(basm1,2)<nn+m) then
            call reallocate ( basm1,maxgrid,nn+m,maxlcutm,ntype )
            call reallocate ( dbasm1,       nn+m,maxlcutm,ntype )
          endif
          do j = 1,nn
            m = m + 1
            do i = 1,ng
              basm1(i,m,l,itype) = dot_product(basp(i,:n),eigv(:,j)) / sqrt(eig(j))
            enddo
            dbasm1(m,l,itype) = dot_product(dbasp(:n),eigv(:,j)) / sqrt(eig(j))
          enddo
          if(l==0) then
            ! Do a Gram-Schmidt orthonormalization to the constant function
            do i = 2,nn+1
              do j = 1,i-1
                rdum                 = intgrf(basm1(:ng,i,0,itype)*basm1(:ng,j,0,itype),itype)
                basm1(:ng,i,0,itype) = basm1(:ng,i,0,itype) - rdum * basm1(:ng,j,0,itype)
                dbasm1(   i,0,itype) = dbasm1(   i,0,itype) - rdum * dbasm1(   j,0,itype)
              enddo
              rdum                 = sqrt(intgrf(basm1(:ng,i,0,itype)**2,itype))
              basm1(:ng,i,0,itype) = basm1(:ng,i,0,itype) / rdum
              dbasm1(   i,0,itype) = dbasm1(   i,0,itype) / rdum
            enddo
            nn              = nn + 1 ; deallocate ( olap ) ; allocate ( olap(nn,nn) )
            nindxm(l,itype) = nn
          else if(l_multipole)  then
            ! Make MT functions multipole-free except the first (Gram-Schmidt for l=0 has the same effect)
            do i = 2,nn
              mp1  = intgrf(basm1(:ng,1,l,itype)*rgrid(:ng,itype)**(l+1),itype)
              mp2  = intgrf(basm1(:ng,i,l,itype)*rgrid(:ng,itype)**(l+1),itype)
              rdum = sqrt(mp1**2+mp2**2)
              basm1(:ng,[1,i],l,itype) = matmul( basm1(:ng,[1,i],l,itype) , reshape( [mp1,mp2,mp2,-mp1] , [2,2] ) ) / rdum
              dbasm1([1,i],l,itype)    = matmul( dbasm1([1,i],l,itype),     reshape( [mp1,mp2,mp2,-mp1] , [2,2] ) ) / rdum
            enddo
          endif
c          do i = 1,nn
c            write(214,*) i,l,intgrf(basm1(:ng,i,l,itype)*rgrid(:ng,itype)**(l+1),itype)
c          enddo
          ! Check orthonormality of product basis
          do n2 = 1,nn
            do n1 = 1,nn
              olap(n1,n2) = intgrf(basm1(:ng,n1,l,itype)*basm1(:ng,n2,l,itype),itype)
            enddo
          enddo
          rdum = sum ( (identity(nn)-olap(:nn,:nn))**2 )
          rdum = sqrt(rdum)/nn
          if(rdum>1d-6) Error('Bad orthonormality of '//lchar(l)//'-product basis ('//Chr(rdum)//'). Increase tolerance.')
          write(6,'(6X,A,I4,'' ->'',I4,''   ('',ES8.1,'' )'')') lchar(l)//':',n,nn,rdum
          deallocate(olap,eigv,eig,ihelp,basp,dbasp)
        enddo
      enddo

      maxindxm = maxval(nindxm)
      allocate ( basm(maxgrid,maxindxm,0:maxlcutm,ntype) )
      allocate ( dbasm(       maxindxm,0:maxlcutm,ntype) )
      basm  = basm1
      dbasm = dbasm1
      deallocate ( basm1,dbasm1 )

      if(iand(restart,W_mbas)/=0) call write_mixedbasis
      deallocate ( seleco,selecu )

 3    spin = nspin
      if(nspin>1) then
        do i = 1,njob
          if(allocated(job(i)%spin)) then
            if(size(job(i)%spin)>0) then
              if(job(i)%spin(1)>=3) spin = 4
            endif
          endif
        enddo
      endif

      ! determine number of basis-function products (nindxp,nindxp2)
      allocate ( nindxp(0:maxlcutm,ntype,0:1) )
      nindxp = 0
      i      = 0
      do itype = 1,ntype
        do l = 0,lcutm(itype)
          nn = 0
          do l1 = 0,lcut(itype)
            do l2 = abs(l-l1),min(l+l1,lcut(itype))
              do n1 = 1,nindx(l1,itype)
                do n2 = 1,nindx(l2,itype)
                  nindxp(l,itype,1) = nindxp(l,itype,1) + 1   !__  all products (for up/down, down/up)
                  if(l2<=l1.and..not.(l2==l1.and.n2>n1)) then !
                    nindxp(l,itype,0) = nindxp(l,itype,0) + 1 !    products without a*b / b*a (for up/up, down/down)
                  endif
                enddo
              enddo
            enddo
          enddo
        enddo
      enddo

      maxindxp = maxval(nindxp)
      allocate ( prodm(maxindxm,maxindxp),prod(maxindxp) )

c One day I should take some time to rewrite the following ...

      if(chkprod) write(6,'(2X,A)') 'Accuracy of MT product basis:'
      allocate ( basprod(maxgrid,maxindxp),
     &           olap(maxindxp,maxindxp),olap2(maxindxp,maxindxp),work(maxindxp),
     &           error(maxband-bandu+1,bando),norm(maxband-bandu+1,bando) )
      prod%l1 = 0
      prod%l2 = 0
      prod%n1 = 0
      prod%n2 = 0

      do ispin=1,spin
        ldum = chkprod.and.spin>1
        select case(ispin)
          case(1) ; ispin1 = 1 ; ispin2 = 1 ; if(ldum) write(6,'(2X,A)') 'Spin up/up'
          case(2) ; ispin1 = 2 ; ispin2 = 2 ; if(ldum) write(6,'(2X,A)') 'Spin down/down'
          case(3) ; ispin1 = 1 ; ispin2 = 2 ; if(ldum) write(6,'(2X,A)') 'Spin up/down'
          case(4) ; ispin1 = 2 ; ispin2 = 1 ; if(ldum) write(6,'(2X,A)') 'Spin down/up'
          case default ; Error('unknown spin index.')
        end select
        no = bando
        nu = min(nband(1,ispin1),nband(1,ispin2))-bandu+1
        allocate ( chelp1(maxindxp,nu,no),chelp2(maxindxp,nu,no) )
        do itype=1,ntype
          if(chkprod) then
            if(ntype>1) write(6,'(4X,A,I3)') 'Atom type',itype
            write(6,'(6X,A)') 'Average (maximum) relative error wrt basis-function products'
          endif
          ng    = grid(itype)%number
          error = 0
          norm  = 0
          do l = 0,lcutm(itype)
            prodm = 0
            i     = 0
            do l1 = 0,lcut(itype)
              ll = lcut(itype)
              if(ispin<=2) ll = l1
              do l2 = 0,ll
                if(l>=abs(l1-l2).and.l<=l1+l2) then

                  do n1 = 1,nindx(l1,itype)
                    nn = nindx(l2,itype)
                    if(l1==l2.and.ispin<=2) nn = n1
                    do n2=1,nn
                      i = i + 1 ; if(i>maxindxp) Error('maximal product index exceeded.')
                      prod(i)%l1 = l1
                      prod(i)%l2 = l2
                      prod(i)%n1 = n1
                      prod(i)%n2 = n2

                      basprod(:ng,i) = ( bas1(:ng,n1,l1,itype,ispin1)*bas1(:ng,n2,l2,itype,ispin2)   ! Calculate all basis-function products to obtain
     &                                 + bas2(:ng,n1,l1,itype,ispin1)*bas2(:ng,n2,l2,itype,ispin2) ) ! the overlaps with the product-basis functions (prodm)
     &                                 / rgrid(:ng,itype)
                      work(i) = intgrf(basprod(:ng,i)**2,itype)
                      do j=1,nindxm(l,itype)
                        prodm(j,i) = intgrf(basprod(:ng,i)*basm(:ng,j,l,itype),itype)
                      enddo
                    enddo
                  enddo

                endif
              enddo
            enddo
            if(.not.chkprod) cycle
            n  = i
            nn = nindxm(l,itype)
            ! calculate olap2 = overlap of basis-function products
            do i=1,n
              do j=1,i
                olap2(i,j) = intgrf(basprod(:ng,i)*basprod(:ng,j),itype)
                olap2(j,i) = olap2(i,j)
              enddo
            enddo
            ! calculate difference of basis-function products and their approximations with the product basis
            do i=1,n
              do j=1,ng
                basprod(j,i) = dot_product(prodm(1:nn,i),basm(j,1:nn,l,itype)) - basprod(j,i)
              enddo
            enddo
            ! calculate olap = overlap of differences
            do i=1,n
              do j=1,i
                olap(i,j) = intgrf(basprod(:ng,i)*basprod(:ng,j),itype)
                olap(j,i) = olap(i,j)
              enddo
              work(i) = olap(i,i) / work(i)
            enddo
            write(6,'(6X,A,F9.6,''  ('',F9.6,'' )'')') lchar(l)//':',sqrt(sum(work(1:n)**2)/n),maxval(work(1:n))

c           Accuracy wrt wave functions

            icent = sum(neq(1:itype-1))+1 ! take first of equivalent centers

            do m=-l,l

              n      = nindxp(l,itype,spin/3)
              chelp1 = 0

              do i=1,n

                l1 = prod(i)%l1
                l2 = prod(i)%l2
                n1 = prod(i)%n1
                n2 = prod(i)%n2

                offdiag = l1/=l2.or.n1/=n2 ! offdiag=true means that a*b is different from b*a.

                lm1_0 = sum ( [ (nindx(ll,itype)*(2*ll+1),ll=0,l1-1) ] )
                lm2_0 = sum ( [ (nindx(ll,itype)*(2*ll+1),ll=0,l2-1) ] )

                do m1=-l1,l1
                  m2 = m1 - m ! condition: m1 - m2 = m
                  if(abs(m2)<=l2) then
                    lm1  = lm1_0 + n1
                    lm2  = lm2_0 + n2 + (m2+l2)*nindx(l2,itype)
                    rdum = gaunt(l2,l1,l,m2,m1,m)
                    do iband = 1,no
                      cdum              = rdum * conjg(cmt(lm1,icent,iband,1,ispin1))
                      chelp1(i,:,iband) = chelp1(i,:,iband) + cmt(lm2,icent,no+1:,1,ispin2) * cdum
                    enddo
                  endif
                  m2 = m1 + m ! condition: - m1 + m2 = m
                  if(abs(m2)<=l2.and.offdiag) then
                    lm1  = lm1_0 + n1
                    lm2  = lm2_0 + n2 + (m2+l2)*nindx(l2,itype)
                    rdum = gaunt(l1,l2,l,m1,m2,m)
                    do iband = 1,no
                      cdum              = rdum * conjg(cmt(lm2,icent,iband,1,ispin1))
                      chelp1(i,:,iband) = chelp1(i,:,iband) + cmt(lm1,icent,no+1:,1,ispin2) * cdum
                    enddo
                  endif
                  lm1_0 = lm1_0 + nindx(l1,itype)
                enddo

              enddo

              do j=1,no
                do i=1,nu
                  chelp2(:n,i,j) = matmul  (olap(:n,:n),chelp1(:n,i,j))
                  error(i,j)     = error(i,j) + dot_product(chelp1(:n,i,j),chelp2(:n,i,j))
                  chelp2(:n,i,j) = matmul (olap2(:n,:n),chelp1(:n,i,j))
                  norm(i,j)      = norm(i,j)  + dot_product(chelp1(:n,i,j),chelp2(:n,i,j))
                enddo
              enddo

            enddo


          enddo

          if(chkprod) then
            error = error / norm
            write(6,'(6X,A)') 'Average (maximum) relative error wrt wave-function products','(occ-unocc comb. at Gamma-point)'
            write(6,'(5X,F12.9,''  ('',F12.9,'' )'')') sqrt(sum(error(:nu,:no)**2)/no/nu),maxval(error(:nu,:no))
          endif

        enddo

        deallocate (chelp1,chelp2)

      enddo

# if 0
 123  read(*,*) i
      if(i==0) goto 124
      do itype = 1,ntype
        ng = grid(itype)%number
        do l = 0,lcutm(itype)
          if(nindxm(l,itype)==0) cycle
          rdum = 0
          do i = 1,nindxm(l,itype)
            rdum1 = intgrf(rgrid(:ng,itype)**(l+1)*basm(:ng,i,l,itype),itype)
            if(abs(rdum1)>rdum) then ; n = i ; rdum = rdum1 ; endif
          enddo
          do i = 1,nindxm(l,itype)
            if(i==n) cycle
            rdum1             = intgrf(rgrid(:ng,itype)**(l+1)*basm(:ng,i,l,itype),itype) ; write(*,*) 'rdum1',rdum1
            rdum              = intgrf(rgrid(:ng,itype)**(l+1)*basm(:ng,n,l,itype),itype)
            basprod(:,1)      = basm(:,n,l,itype)
            basm(:,n,l,itype) = rdum /sqrt(rdum**2+rdum1**2) * basprod(:,1) + rdum1/sqrt(rdum**2+rdum1**2) * basm(:ng,i,l,itype)
            basm(:,i,l,itype) = rdum1/sqrt(rdum**2+rdum1**2) * basprod(:,1) - rdum /sqrt(rdum**2+rdum1**2) * basm(:ng,i,l,itype)
          enddo
          do i = 1,nindxm(l,itype)
            write(*,'(F15.10$)') intgrf(rgrid(:ng,itype)**(l+1)*basm(:ng,i,l,itype),itype)
            do j = 1,nindxm(l,itype)
              write(6,'(F10.5$)') intgrf(basm(:ng,i,l,itype)*basm(:ng,j,l,itype),itype)
            enddo
            write(*,*)
          enddo
        enddo
      enddo
  124 continue
# endif

      nbasp   = 0
      maxbasp = 0
      do itype=1,ntype
        do i = 1,neq(itype)
          n = 0
          do l = 0,lcutm(itype)
            do m = -l,l
              do j = 1,nindxm(l,itype)
                nbasp = nbasp + 1
                n     = n     + 1
              enddo
            enddo
          enddo
          if(n>maxbasp) maxbasp = n
        enddo
      enddo
      maxbasm    = nbasp  + maxgptm
      nbasm      = nbasp  + ngptm
      maxlmindxm = 0
      do itype = 1,ntype
        maxlmindxm = max ( maxlmindxm , sum ( [ ((2*l+1)*nindxm(l,itype),l=0,lcutm(itype)) ] ) )
      enddo
      write(6,'(4X,A,I5)') 'Total number of product-basis functions:',nbasp

      deallocate (basprod,olap,olap2)

# if 0
      if(bandinfo==0.and.(any(job%type==J_HF).or.any(job%type==J_GW))) then
        call teststop('bandinfo with new routines')
        write(6,'(/A)') 'Band info'
        write(6,'(A)')  '  Square integral of wave-function products  INT |conjg(phi_{n,k}(r))*phi_{m,k} (r)|^2 dr'
        write(6,'(A)')  '  n,k (rows)    = all bands'
        write(6,'(A)')  '  m,k (columns) = bands as defined by JOB QP/HF'
        n = 0
        do m = 1,njob
          if(job(m)%type/=J_HF.and.job(m)%type/=J_GW) cycle
          n = n + size(job(m)%band)
        enddo
        allocate ( cprod(nbasm(1),maxband,n),olappw(ngptm(1),ngptm(1)) )
        call olap_pw(olappw,gptm(:,pgptm(:ngptm(1),1)),ngptm(1))
        i = 0
        do m = 1,njob
          if(job(m)%type/=J_HF.and.job(m)%type/=J_GW) cycle
          do j = 1,size(job(m)%band)
            i = i + 1
            call wavefproducts1_mt(cprod(:,:,i),nbasm(1),job(m)%band(j),1,job(m)%kpt(j),job(m)%spin(j),[1],1,1,maxband)
            call wavefproducts1_pw(cprod(:,:,i),nbasm(1),job(m)%band(j),1,job(m)%kpt(j),job(m)%spin(j),[1],1,1,maxband)
          enddo
        enddo
        do j = 1,maxband
          write(6,'(I3'NoA) j
          do i = 1,n
            rdum = dotprod( cprod(:nbasp,j,i) , cprod(:nbasp,j,i) ) +
     &             dotprod( cprod(nbasp+1:nbasm(1),j,i) , matmul ( olappw , cprod(nbasp+1:nbasm(1),j,i) ) )
            write(6,'(F7.4'NoA) rdum
          enddo
          write(6,*)
        enddo
        deallocate ( cprod,olappw )
      endif
# endif
       
      call cpu_time(rdum) ; cputime = rdum - cputime
      write(6,'(/A,F12.5)') 'Timing:',cputime

      Rend

# ifdef MPI
#   define CB call Mcast
#   define CA call Mcastl
      MnoR( if(allocated(cstep)) tDeallocate(cstep) )
      CB(maxbasm);  CA(nbasm)
      CB(maxindxm); CB(maxlmindxm); CB(nbasp); CB(maxbasp); CB(maxindxp)
      CA(nindxm);   CA(basm);       CA(dbasm)
      CB(maxgptm);  CB(ngptmall);   CA(ngptm); CA(pntgptm); CA(gptm)
      CA(pgptm);    CA(cstep);      CB(tolerance)
      CB(l_multipole) ! might have been redefined from "spex.mb"
      MnoR( CHKMEM(cstep) )
#   undef CA
#   undef CB
# endif

      contains

c     ----------------

      subroutine lodef(out,in)
      implicit none
      character(80), intent(out) :: out
      character(80), intent(in)  :: in
      integer                    :: i,j,k,n
      out = ' '
      j   = 0
      i   = 1
      do while(i<=len(in).and.in(i:i)/=' ')
        if(in(i:i)=='0'.or.in(i:i)=='1') then
          j        = j + 1
          out(j:j) = in(i:i) ; if(j>len(out)) Error('Out of dimension.')
        else
          n = ch2i(in(i:i),'SELECT (lodef)')
          i = i + 1
          if(in(i:i)/='0'.and.in(i:i)/='1') Error('Wrong LO definition (expected 0 or 1).')
          do k = 1,n
            j        = j + 1 ; if(j>len(out)) Error('Out of dimension.')
            out(j:j) = in(i:i)
          enddo
        endif
        i = i + 1
      enddo
      end subroutine lodef

c     ----------------

      subroutine read_mixedbasis(lsuccess)
      use file
      implicit none
      logical, intent(out) :: lsuccess
      integer, allocatable :: selo(:,:,:),selu(:,:,:)
      integer              :: iunit,n0,n1,iarr(ntype)
      real_dp              :: rdum,rmat(3,3)
      integer              :: l1,l2,version,mode1
      inquire(file='spex.mb',exist=lsuccess)
      if(lsuccess) then
        iunit = fopen('spex.mb',form='unformatted',status='old')
        mode1 = 0
        ! Preamble
        read(iunit) version   ; if(version<1) then ; read(iunit) mode1 ; else ; backspace(iunit) ; endif
        ! MB parameters
        read(iunit) n0        ; if(n0/=ntype) goto 1
        read(iunit) rmat      ; if(any(rmat/=lat)) goto 1
        read(iunit) rdum,iarr ; if(rdum/=tolerance.or.any(iarr/=lcutm)) goto 1
        read(iunit) n0,n1     ; allocate ( selo(n0,0:n1,ntype),selu(n0,0:n1,ntype) )
        read(iunit) selo,selu ; if(any((selo/=0).neqv.seleco).or.any((selu/=0).neqv.selecu)) goto 1 ! see (*) in readwrite
        read(iunit) n0,rdum   ; if(n0/=optm.or.rdum/=optcutm) goto 1
        read(iunit) l1,l2     ; if((l1/=0).or.(l2/=0).neqv.noapw) goto 1
        read(iunit) n0        ; if(n0/=maxgrid) goto 1
        deallocate ( selo,selu )
        ! MB MT functions
        allocate ( nindxm(0:maxlcutm,ntype) )
        read(iunit) nindxm ; maxindxm = maxval(nindxm)
        allocate ( basm(maxgrid,maxindxm,0:maxlcutm,ntype),dbasm(maxindxm,0:maxlcutm,ntype) )
        read(iunit) basm
        read(iunit) dbasm
        if(version==-2) then
          read(iunit) l1
          if(l1==1.and..not.l_multipole) then
            Info('MT basis from "spex.mb" MULTIPOLE-ordered but not used because MULTIPOLE OFF specified in input.')
          else if(l1==0.and.l_multipole) then
            Info('MT basis from "spex.mb" not MULTIPOLE-ordered; MULTIPOLE must be switched OFF.')
            l_multipole = .false.
          endif
        else
          l_multipole = .false.
        endif
        call fclose(iunit)
        if(mode1==0) then ; write(6,'(/A)')     '  MT basis read in from spex.mb (old version).'
        else              ; write(6,'(/A,A,A)') '  MT basis read in from spex.mb (mode: ',trim(chmode(mode1)),').'
        endif
      endif
      return
 1    call fclose(iunit)
      write(6,'(A)') '  MT parameters changed from spex.mb.'
      lsuccess = .false.
      end subroutine read_mixedbasis

c     ----------------

      subroutine write_mixedbasis
      use file
      implicit none
      integer :: selo(maxindx,0:maxlcut,ntype),selu(maxindx,0:maxlcut,ntype)
      integer :: iunit,i
      where(seleco) ; selo = 1 ; elsewhere ; selo = 0 ; endwhere
      where(selecu) ; selu = 1 ; elsewhere ; selu = 0 ; endwhere
      iunit = fopen('spex.mb',form='unformatted',status='unknown')
      ! Preamble
      write(iunit) -2 ! version number (non-positive)
      write(iunit) mode
      ! MB parameters
      write(iunit) ntype
      write(iunit) lat
      write(iunit) tolerance,lcutm
      write(iunit) maxindx,maxlcut
      write(iunit) selo,selu
      write(iunit) optm,optcutm ; i = 0 ; if(noapw) i = 1
      write(iunit) 0,i
      write(iunit) maxgrid
      ! MB MT functions
      write(iunit) nindxm
      write(iunit) basm
      write(iunit) dbasm
      if(l_multipole) then ; write(iunit) 1
      else                 ; write(iunit) 0
      endif
      call fclose(iunit)
      write(6,'(A)') '    MT basis written to spex.mb.'
      end subroutine write_mixedbasis

c     ----------------

      end

c     ----------------

      subroutine wavefunction_test(psi,iband,ikpt,ispin,itype,ieq,theta,phi) ! for testing
      use global
      use, intrinsic :: iso_fortran_env
      implicit none
      integer    :: iband,ikpt,itype,ieq,ispin
      real_dp    :: theta,phi
      complex_dp :: psi(maxgrid,2),y((lcut(itype)+1)**2),cdum
      integer    :: l,m,i,lm,lmi,icent
      call harmonics(y,theta,phi,lcut(itype))
      icent = sum(neq(1:itype-1))+ieq
      psi = 0
      lm  = 0
      lmi = 0
      do l=0,lcut(itype)
        do m=-l,l
          lm = lm + 1
          do i=1,nindx(l,itype)
            lmi = lmi + 1
            cdum = cmt(lmi,icent,iband,ikpt,ispin) * y(lm)
            psi(:,1) = psi(:,1) + cdum*bas1(:,i,l,itype,ispin)
            psi(:,2) = psi(:,2) + cdum*bas2(:,i,l,itype,ispin)
          enddo
        enddo
      enddo
      end

c     ------------------

      function wavefunc(iband,ikpt,ispin,itype,ieq,r)
      use global
      use, intrinsic :: iso_fortran_env
      implicit none
      complex_dp :: wavefunc,y((maxlcut+1)**2),bas(maxlmindx)
      integer    :: iband,ikpt,ispin,itype,ieq,ic,i
      real_dp    :: r(3),rr,interpolated
      integer    :: l,m,lm,n,lm1
      ic = 0
      do i = 1,itype-1
        ic = ic + neq(i)
      enddo
      ic = ic + ieq
      rr = sqrt(sum(r**2))
      if(rr>grid(itype)%radius) Error('MT radius exceeded.')
      if(rr==0) Error('point falls into MT center.')
      call harmonicsr(y,r,lcut(itype))
      bas = 0
      lm  = 0
      lm1 = 0
      do l = 0,lcut(itype)
        do m = -l,l
          lm  = lm + 1
          do n = 1,nindx(l,itype)
            lm1 = lm1 + 1
            bas(lm1) = interpolated(rr,rgrid(:,itype),bas1(:,n,l,itype,ispin),grid(itype)%number)/rr * y(lm)
          enddo
        enddo
      enddo
      wavefunc = sum(cmt(:lm1,ic,iband,ikpt,ispin)*bas(:lm1))
      end

      function interpolated(rr,r,f,n)
      use, intrinsic :: iso_fortran_env
      implicit none
      integer, intent(in) :: n
      real_dp, intent(in) :: r(n),f(n),rr
      real_dp             :: interpolated
      real_dp             :: r1(0:n),f1(0:n)
      real_dp             :: interpol
      integer             :: i
      r1(0) = 0 ; r1(1:) = r
      f1(0) = 0 ; f1(1:) = f
      do i = 1,n
        if(rr<=r1(i)) goto 1
      enddo
      Error('MT radius exceeded')
 1    i = max(2,i)
      i = min(i,n-1)
      interpolated = interpol(rr,r1(i-2),r1(i-1),r1(i),r1(i+1),f1(i-2),f1(i-1),f1(i),f1(i+1))
      end

      function interpol(rr,r1,r2,r3,r4,f1,f2,f3,f4)
      use, intrinsic :: iso_fortran_env
      implicit none
      real_dp rr,r1,r2,r3,r4,f1,f2,f3,f4,interpol
      interpol = ((rr-r2)*(rr-r3)*(rr-r4)) / ((r1-r2)*(r1-r3)*(r1-r4)) * f1 +
     &           ((rr-r1)*(rr-r3)*(rr-r4)) / ((r2-r1)*(r2-r3)*(r2-r4)) * f2 +
     &           ((rr-r1)*(rr-r2)*(rr-r4)) / ((r3-r1)*(r3-r2)*(r3-r4)) * f3 +
     &           ((rr-r1)*(rr-r2)*(rr-r3)) / ((r4-r1)*(r4-r2)*(r4-r3)) * f4
      end

