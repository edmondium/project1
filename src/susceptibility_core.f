c--------------------------------------------------------------------------------
c Copyright (c) 2021 Peter Grünberg Institut, Forschungszentrum Jülich, Germany
c This file is part of SPEX and available as free software under the conditions
c of the MIT license as expressed in the LICENSE file in more detail.
c--------------------------------------------------------------------------------

c Calculates the core contribution to the susceptibility in a pure MT basis set and stores it on hard disk ("spex.core").
c
c It is later read and transformed to the usual mixed basis in the routines susceptibility1/2.
c
c Here, the basis functions are identical to the corresponding (unsymmetrized) MT functions M_aLMP(r) but lack the phase
c factor exp(ikR_a). In this way the core contribution is k-independent and must only be calculated once.
c
c The input parameters are the list of (complex) frequencies frq(1:nfrq) and
c spin = 0 :                "uu+dd"
c spin = 1 : s = 1, s' = 1  "uu"
c spin = 2 : s = 2, s' = 2  "dd"
c spin = 3 : s = 1, s' = 2  "ud"
c spin = 4 : s = 2, s' = 1  "du" .
c
c The projections cprod that are calculated below are
c
c   _          c,s1       s2                 L+l            s2                           R_MT  c,s1     s2
c < M       phi      | phi    > = delta      SUM    SUM  cmt               C              INT u    (r) u     (r) M   (r) dr / sqrt(N)
c    k,aLMP    n,q-k      n'q          aa' l'=|L-l|  p'     n'q,al'(m+M)p'  lm,l'(m+M),LM  0   alp      al'p'     aLP
c
c with the core state n=(a'lmp) .
c
c If set, the core states are orthogonalized to the basis functions.
# define Orthog_core
c
c Use of symmetry can be switched off here.
c # define switch_off_symmetry_core
c
c Write detailed timing information.
c # define timing

# ifndef Orthog_core
#   warning Orthog_core undefined
# endif      

# ifdef switch_off_symmetry_core
#   warning switch_off_symmetry_core defined
#   define NKPTI nkpt
#   define KPTP [(i,i=1,nkpt)]
# endif
# ifdef timing
#   define Time(arg) arg
# else
#   define Time(arg)
# endif

# include "cppmacro.h"
# include "jobtype.h"
# include "restype.h"

      subroutine susceptibility_core(spin,frq,nfrq)

      use global
      use file
      use key
      use wrapper
      use, intrinsic :: iso_fortran_env

      implicit none
      integer,    intent(in)  :: spin,nfrq
      complex_dp, intent(in)  :: frq(nfrq)
      complex_dp, allocatable :: suscep(:,:,:,:),wing1(:,:,:,:),wing2(:,:,:,:),head(:,:,:)
      complex_dp, allocatable :: s0(:,:,:),w10(:,:,:),w1(:,:,:),w20(:,:,:),w2(:,:,:)
      complex_dp              :: tocart(3,3)
      complex_dp              :: head0(3,3)
      complex_dp, allocatable :: c(:),p(:,:),q(:,:,:),qq(:,:,:),qq0(:,:,:),ww(:,:)
      real_dp,    allocatable :: weight(:,:,:),wghtSr(:,:),wghtSi(:,:)
      real_dp                 :: rot(3,3),vv(3)
      real_dp                 :: integral(maxindx,maxindxm,0:maxlcut,0:maxlcutm),hlp(maxindxm,maxindx)
      real_dp                 :: integralc(maxindx,2)
      real_dp                 :: rdum,rdum1,rdum2,gnt,svol,clebsch
      real_dp                 :: eval(maxindxc,2*maxlcutc+2,0:maxlcutc),evec(2,maxindxc,2*maxlcutc+2,maxlcutc),ham(2,2)
      real_dp,    parameter   :: sq2 = 1/sqrt(2d0)
      integer,    allocatable :: nindex(:,:),index(:,:,:)
      integer                 :: spin1(2),spin2(2),ispin1,ispin2,term(2),nterm,iterm
      integer                 :: iunit,ios
      integer                 :: isym
      integer                 :: ifrqs,itype,ieq,ic,lcut1
      integer                 :: lc,nc,mc,ll,nn,mm,l,n,m,llm,lm0,lm,l_1,l_2,ifrq,nlm,nl
      integer                 :: i,j,k,ij,lk
      logical                 :: define
      real                    :: cputime
      real_dp                 :: intgrf,gaunt
      integer                 :: nfrqs
      real_dp,    allocatable :: frqs(:)

c     Parameters generated from spex.inp
      logical, save           :: firstcall=.true.
      real_dp, save           :: disorder,wghtthr,fspec(2)

# ifdef Orthog_core
      real_dp                 :: core1_tmp(size(core1,1),size(core1,2),size(core1,3),size(core1,4),size(core1,5))
      real_dp                 :: core2_tmp(size(core1,1),size(core1,2),size(core1,3),size(core1,4),size(core1,5))
      real_dp                 :: bas1o(size(bas1,1),size(bas1,2)),bas2o(size(bas1,1),size(bas1,2)),olap(maxindx,maxindx)
# endif

# include "interface/getfrqs.inc"

      svol = 1/sqrt(vol)

      write(6,*)

c     RESTART: Read spex.core if possible.
      inquire(file='spex.core',exist=define)
      if(iand(restart,R_core)/=0.and.define) then
        iunit = fopen('spex.core',form='unformatted',status='old')
        read(iunit,iostat=ios) i,j,k
        if(any([i,j,k,ios]/=[spin,ntype,nfrq,0])) goto 1
        do itype = 1,ntype
          n   = count(cores(:,:,itype))
          nlm = sum( [ ((2*l+1)*nindxm(l,itype),l=0,lcutm(itype)) ] )
          read(iunit) i
          if(i/=n) write(6,'(A,I3,A,I3)') 'spex.core: Number of core shells of atom type',itype,' differs:',n
          if(i==0) cycle
          read(iunit,iostat=ios) i,j
          if(any([i,j,ios]/=[neq(itype),nlm,0])) goto 1
          do i = 1,neq(itype)*nfrq
            read(iunit,iostat=ios) ; if(ios/=0) goto 1
          enddo
        enddo
        call fclose(iunit)
        write(6,'(A)') 'File "spex.core" exists and will be read.'
        return
 1      call fclose(iunit)
        if(ios==0) then ; write(6,'(A/)') 'File "spex.core" exists but parameters have changed.'
        else            ; write(6,'(A/)') 'File "spex.core" exists but cannot be read.'
        endif
      endif

      Load(         Error('Not implemented for LOAD yet.') )
      if(l_soc)     Error('Not implemented for SOC yet.')
# ifdef switch_off_symmetry_core
      if(storeibz)  Error('"switch_off_symmetry_core" requires STOREBZ.')
# else
      if(lcore_soc.and.nspin==2)
     &              Error('CORESOC & nspin==2; please compile with -Dswitch_off_symmetry_core.')
# endif

# ifdef Orthog_core
      core1_tmp = core1
      core2_tmp = core2
# else
      write(0,'(A)') 'susceptibility_core: Note! Orthog_core unset; residual head and wings might be large.'
# endif

# ifdef Orthog_core
      write(6,'(A)')   'Calculation of core susceptibility '
      write(6,'(A)')   '  Core states are orthonormalized wrt basis functions:'
# else
      write(6,'(A'NoA) 'Calculation of core susceptibility... '
# endif

c     Input from spex.inp
      if(firstcall) then
        firstcall = .false.
        call getkey(inp,'WGHTTHR',   wghtthr, section='SUSCEP', default=1d-10, mini=0d0)
        call getkey(inp,'DISORDER', disorder, section='SUSCEP', default=0d0,   mine=0d0)
        if(disorder/=0) disorder = 1/(2*disorder) ! eta = 1/(2*tau)
        call getfspec(fspec,0,nfrq,frq,.false.)
        if(fspec(1)==0) Error('CORES not implemented for HILBERT NONE.')
      endif
      if(gauss(1)/=0) then ; rdum1 = 0
      else                 ; rdum1 = eip + egap - efermi ! dispersionless core state set at Fermi energy (see below)
      endif
      rdum2 = maxval(ene(:maxeband,:nkpt,:)) - efermi
      call getfrqs(frqs,nfrqs,fspec,rdum1,rdum2,.true.)

      allocate(wghtSr(nfrqs,nfrq))
      allocate(wghtSi(nfrqs,nfrq))

      call cpu_time(cputime)
      iunit = fopen('spex.core',form='unformatted',status='unknown')
      write(iunit) spin,ntype,nfrq

      ! Spin algebra
      select case(spin)
        case(0) ; if(nspin1==2) then ; spin1 = [1,2] ; spin2 = [1,2] ; term = [0,0] ; nterm = 2
                  else               ; spin1 =  1    ; spin2 =  1    ; term =  0    ; nterm = 1
                  endif
        case(1) ; spin1 =  1    ; spin2 =  1    ; term =  0    ; nterm = 1
        case(2) ; spin1 =  2    ; spin2 =  2    ; term =  0    ; nterm = 1
        case(3) ; spin1 = [1,2] ; spin2 = [2,1] ; term = [1,2] ; nterm = 2
        case(4) ; spin1 = [2,1] ; spin2 = [1,2] ; term = [1,2] ; nterm = 2
      end select

      ! Determine BZ integration weights (e.g., tetrahedron weights)
      allocate (      c((maxeband-bandu+1)*NKPTI) )
      allocate ( weight((maxeband-bandu+1)*NKPTI,nfrqs,nterm), nindex(nfrqs,nterm) ) ; weight = 0
      do iterm = 1,nterm
        ispin1 = spin1(iterm)
        ispin2 = spin2(iterm)
        if(gauss(1)/=0) then
          call gauss4_init      (weight(:,:,iterm),frqs,nfrqs,1,[(i,i=1,NKPTI)],KPTP,NKPTI,1,ispin1,ispin2,0,0,bandu,maxeband)
        else
          call tetrahedron5_init(weight(:,:,iterm),frqs,nfrqs,1,[(i,i=1,NKPTI)],KPTP,NKPTI,1,ispin1,ispin2,0,0,bandu,maxeband, ! "0" is flag for core states (see tetrahedron.f)
     &                                                                                                     0,0,bandu,maxeband) ! (only ispin2 is relevant)
        endif
        nindex(:,iterm) = [ (count(weight(:,ifrqs,iterm)>wghtthr),ifrqs=1,nfrqs) ]
      enddo
      allocate ( index(maxval(nindex),nfrqs,nterm) )
      do iterm = 1,nterm
        do ifrqs = 1,nfrqs
          n = 0
          do i = 1,size(weight,1)
            if(weight(i,ifrqs,iterm)>wghtthr) then
              n                     = n + 1
              index(n,ifrqs,iterm)  = i
              weight(n,ifrqs,iterm) = weight(i,ifrqs,iterm)
            endif
          enddo
        enddo
      enddo

      allocate(head(3,3,nfrq))
      head = 0

      ! Transformation matrix "tocart": Y1m->cartesian basis (for head/wings)
      tocart      = reshape ( [ sq2,-sq2,0d0 , 0d0,0d0,1d0 , -sq2,-sq2,0d0 ] , [3,3] ) * 2 * sqrt(pi/3)
      tocart(2,:) = tocart(2,:) * (-img)

      ! Loop over atom types
      ic = 0
      do itype = 1,ntype

        n = count(cores(:,:,itype))
        write(iunit) n
        if(n==0) then
          ic = ic + neq(itype)
          cycle
        endif

        ! Calculate Zeeman-splitted core states (lcore_soc & nspin==2)
        if(lcore_soc.and.nspin==2) then
          eval = 0
          do lc = 0,lcutc(itype)
            do nc = 1,nindxc(lc,itype),min(lc+1,2)
              vv(1) = vmat(nc,nc)
              if(lc>0) then
                vv(2) = vmat(nc+1,nc+1)
                vv(3) = vmat(nc+1,nc)
              endif
              mm = 0
              do m = -(2*lc+1),2*lc+1,2 ; mm = mm + 1
                if(abs(m)==2*lc+1) then
                  eval(  nc,mm,lc) = ecore(nc,lc,itype,1) + m*vv(1) / (2*lc+1)
                  evec(:,nc,mm,lc) = [ 1d0,0d0 ]
                  if(l>0) evec(:,nc+1,mm,lc) = 0
                else
                  ham(1,1) = ecore(nc,  lc,itype,1) + m*vv(1) / (2*lc+1)
                  ham(2,2) = ecore(nc+1,lc,itype,1) - m*vv(2) / (2*lc+1)
                  ham(1,2) = -sqrt(1d0*(2*lc+1-m)*(2*lc+1+m))*vv(3) / (2*lc+1)
                  ham(2,1) = ham(1,2)
                  call diagonalize(evec(:,nc:nc+1,mm,lc),eval(nc:nc+1,mm,lc),ham)
                  rdum = eval(  nc,mm,lc) ; eval(  nc,mm,lc) = eval(  nc+1,mm,lc) ; eval(  nc+1,mm,lc) = rdum
                  rdum = evec(1,nc,mm,lc) ; evec(1,nc,mm,lc) = evec(1,nc+1,mm,lc) ; evec(1,nc+1,mm,lc) = rdum
                  rdum = evec(2,nc,mm,lc) ; evec(2,nc,mm,lc) = evec(2,nc+1,mm,lc) ; evec(2,nc+1,mm,lc) = rdum
                endif
              enddo
            enddo
          enddo
        endif

#       ifdef Orthog_core
        ! Orthogonalize core to basis
        if(spin<3) then
          write(6,'(A,I3,A'NoA) '  Atom type',itype,':'
          do iterm = 1,nterm
            ispin1 = spin1(iterm)
            do lc = 0,lcutc(itype)
              if(any(cores(:,lc,itype))) then
                do n = 1,nindx(lc,itype) ; do nn = 1,nindx(lc,itype)
                  olap(n,nn) = intgrf( bas1(:,n,lc,itype,ispin1) * bas1(:,nn,lc,itype,ispin1) +
     &                                 bas2(:,n,lc,itype,ispin1) * bas2(:,nn,lc,itype,ispin1) , itype )
                enddo ; enddo
                n           = nindx(lc,itype) ; call sqrtmat(olap(:n,:n))
                olap        = transpose(invert(olap(:n,:n)))
                bas1o(:,:n) = matmul ( bas1(:,:n,lc,itype,ispin1) , olap(:n,:n) )
                bas2o(:,:n) = matmul ( bas2(:,:n,lc,itype,ispin1) , olap(:n,:n) )
                do n = 1,nindx(lc,itype) ; do nn = 1,nindx(lc,itype)
                  if(abs(intgrf(bas1o(:,n)*bas1o(:,nn)+bas2o(:,n)*bas2o(:,nn),itype)-min(n,nn)/max(n,nn))>1d-10)
     &              Bug('(Orthog_core) Orthogonalization of radial basis failed.')
                enddo ; enddo
                do nc = 1,nindxc(lc,itype)
                  do n = 1,nindx(lc,itype)                    
                    rdum = intgrf ( core1(:,nc,lc,itype,ispin1) * bas1o(:,n) +
     &                              core2(:,nc,lc,itype,ispin1) * bas2o(:,n) , itype )
                    core1(:,nc,lc,itype,ispin1) = core1(:,nc,lc,itype,ispin1) - rdum * bas1o(:,n)
                    core2(:,nc,lc,itype,ispin1) = core2(:,nc,lc,itype,ispin1) - rdum * bas2o(:,n)
                  enddo
                  rdum = sqrt ( intgrf ( core1(:,nc,lc,itype,ispin1)**2 + core2(:,nc,lc,itype,ispin1)**2 , itype ) )
                  if(rdum<.7d0) Warn('Orthogonalization to valence basis reduced core norm to '//Chr(rdum))
                  if(cores(nc,lc,itype)) write(6,'(F6.1,''%'''NoA) rdum*100
                  core1(:,nc,lc,itype,ispin1) = core1(:,nc,lc,itype,ispin1) / rdum
                  core2(:,nc,lc,itype,ispin1) = core2(:,nc,lc,itype,ispin1) / rdum
                  do n = 1,nindx(lc,itype)
                    if(abs( intgrf ( core1(:,nc,lc,itype,ispin1) * bas1(:,n,lc,itype,ispin1) +
     &                               core2(:,nc,lc,itype,ispin1) * bas2(:,n,lc,itype,ispin1) , itype ) ) >1d-10)
     &                Bug('(Orthog_core) Core orthogonalization failed.')
                  enddo
                enddo
              endif
            enddo
          enddo
          write(6,*)
        endif
#       endif

        nlm = sum( [ ((2*l+1)*nindxm(l,itype),l=0,lcutm(itype)) ] )
        allocate( suscep(nlm,nlm,neq(itype),nfrq) )
        allocate( wing1 (nlm,3,  neq(itype),nfrq) )
        allocate( wing2 (3,nlm,  neq(itype),nfrq) )
        suscep = 0
        wing1  = 0
        wing2  = 0
        do ieq = 1,neq(itype)
          ic = ic + 1

          ! Determine p(ij,ws) = sum(n=unocc) c(qn,i) c(qn,j) delta(ws-eqn), ws = Hilbert frequencies
          lcut1 = min(lcut(itype),lcutc(itype)+max(1,lcutm(itype)))
          nl    = sum( [ ((2*l+1)*nindx(l,itype),l=0,lcut1) ] )
          allocate( p(nl*(nl+1)/2,nfrqs) )
          allocate( q(nl,nl,nfrq),qq(nlm,nl,nfrq),qq0(3,nl,nfrq) )

          ! Loop over spins
          do iterm = 1,nterm
          ispin1 = spin1(iterm)
          ispin2 = spin2(iterm)

          Time( call cpu_time(cputime) )
          ij = 0
          do j = 1,nl
          do i = 1,j
            ij = ij + 1
            lk = 0
            c = 0
            do k = 1,NKPTI
              m = min(maxeband,nband(k,ispin2))
              call getdeg(n,nn,m,k,ispin2)
              if(nn>m) m = n - 1
              do l = bandu,m
                lk    = lk + 1
                c(lk) = cmt(i,ic,l,k,ispin2) * conjg(cmt(j,ic,l,k,ispin2))
              enddo
              m            = maxeband - m
              c(lk+1:lk+m) = 0
              lk           = lk + m
            enddo
            do ifrqs = 1,nfrqs
              nn          = nindex(ifrqs,iterm)
              p(ij,ifrqs) = dot_product( weight(:nn,ifrqs,iterm),c(index(:nn,ifrqs,iterm)) )
            enddo
          enddo
          enddo
          Time( write(*,'(''p:'''NoA) ; call cpu_done(cputime) )

          ! Loop over core states at current atom (ic)
          do lc = 0,lcutc(itype)
          do nc = 1,nindxc(lc,itype)
            if(cores(nc,lc,itype)) then

              if(.not.lcore_soc.or.nspin==1) then ! otherwise done later
                call perform_hilbert(ecore(nc,lc,itype,ispin1)) ! Hilbert transformation frqs->frq (-> q)
                call calc_radial(nc,[1d0,0d0])                  ! Calculate radial integrals (-> integral(c))
                Time( write(*,'(''q:'''NoA) ; call cpu_done(cputime) )
              endif

              ! Loop over core mc and multiply from left and right with MPBs and nabla
              do mc = -lc,lc

                ! lcore_soc: LS Clebsch-Gordon coefficients
                clebsch = 1
                if(lcore_soc) then
                  j = max(mod(nc,2),1-lc) ! j=0 for lc-1/2   and j=1  for lc+1/2
                  if(nspin==1) then
                    if(spin==0.and.nterm==1) then ; clebsch = 1d0 * (lc                          +j) / (2*lc+1) ! spin-summed
                    else                          ; clebsch = 1d0 * (lc+ mc*(3-ispin1*2)*(2*j-1) +j) / (2*lc+1)
                    endif
                    if(clebsch==0) cycle
                  else
                    m  = 2*mc + (3-2*ispin1) ! 2*mj
                    mm = (m+1)/2 + lc + 1    ! 2*mj -> Counting index
                    if(all(evec(:,nc,mm,lc)==0)) cycle
                    if(eval(nc,mm,lc)==0) Bug('Core energy undefined.')
                    call perform_hilbert(eval(nc,mm,lc))      ! Hilbert transformation (-> q) for lcore_soc and nspin==2
                    call calc_radial(nc-1+j,evec(:,nc,mm,lc)) ! Calculate radial integrals (-> integral(c))
                    Time( write(*,'(''q:'''NoA) ; call cpu_done(cputime) )
                  endif
                endif

                ! Left-multiplication (-> qq) for wing/suscep
                qq  = 0
                llm = 0
                do ll = 0,lcutm(itype)
                  do mm = -ll,ll
                    m   = mc + mm
                    l_1 = max(abs(ll-lc),abs(m)) ; l_1 = l_1 + mod(lc+ll+l_1,2)
                    l_2 = min(ll+lc,lcut1)
                    do l = l_1,l_2,2
                      lm  = sum( [ ((2*i+1)*nindx(i,itype),i=0,l-1) ] ) + (m+l) * nindx(l,itype)
                      n   = nindx(l,itype)
                      nn  = nindxm(ll,itype)
                      hlp = transpose(integral(:,:,l,ll))
                      gnt = gaunt(ll,l,lc,mm,m,mc)
                      do ifrq = 1,nfrq
                        qq(llm+1:llm+nn,:,ifrq) = qq(llm+1:llm+nn,:,ifrq) + matmul( hlp(:nn,:n) , q(lm+1:lm+n,:,ifrq) ) * gnt
                      enddo
                    enddo
                    llm = llm + nindxm(ll,itype)
                  enddo
                enddo
                ! Left-multiplication (-> qq0) for head/wing
                qq0 = 0
                if(lc<lcut(itype)) then
                  lm0 = sum( [((2*l+1)*nindx(l,itype),l=0,lc)] )
                  l   = lc + 1
                  do mm = -1,1
                    lm = lm0 + (l+mc+mm)*nindx(l,itype)
                    do n = 1,nindx(l,itype)
                      lm            = lm + 1
                      gnt           = gaunt(1,lc+1,lc,mm,mc+mm,mc)
                      qq0(mm+2,:,:) = qq0(mm+2,:,:) + q(lm,:,:) * integralc(n,1) * gnt * svol
                    enddo
                  enddo
                endif
                if(lc>0) then
                  lm0 = sum( [((2*l+1)*nindx(l,itype),l=0,lc-2)] )
                  l   = lc - 1
                  do mm = -1,1 ; if(abs(mc+mm)>l) cycle
                    lm = lm0 + (l+mc+mm)*nindx(l,itype)
                    do n = 1,nindx(l,itype)
                      lm            = lm + 1
                      gnt           = gaunt(1,lc-1,lc,mm,mc+mm,mc)
                      qq0(mm+2,:,:) = qq0(mm+2,:,:) + q(lm,:,:) * integralc(n,2) * gnt * svol
                    enddo
                  enddo
                endif
                ! Right-multiplication (-> suscep,wing2)
                llm = 0
                do ll = 0,lcutm(itype)
                  do mm = -ll,ll
                    m   = mc + mm
                    l_1 = max(abs(ll-lc),abs(m)) ; l_1 = l_1 + mod(lc+ll+l_1,2)
                    l_2 = min(ll+lc,lcut1)
                    do l = l_1,l_2,2
                      lm  = sum( [ ((2*i+1)*nindx(i,itype),i=0,l-1) ] ) + (m+l) * nindx(l,itype)
                      n   = nindx(l,itype)
                      nn  = nindxm(ll,itype)
                      gnt = gaunt(ll,l,lc,mm,m,mc) * clebsch
                      do ifrq = 1,nfrq
                        suscep(:,llm+1:llm+nn,ieq,ifrq) = suscep(:,llm+1:llm+nn,ieq,ifrq) +
     &                                                    matmul ( qq(:,lm+1:lm+n,ifrq)  , integral(:n,:nn,l,ll) ) * gnt
                      enddo
                      do ifrq = 1,nfrq
                        wing2(:,llm+1:llm+nn,ieq,ifrq)  = wing2(:,llm+1:llm+nn,ieq,ifrq)  +
     &                                                    matmul ( qq0(:,lm+1:lm+n,ifrq) , integral(:n,:nn,l,ll) ) * gnt
                      enddo
                    enddo
                    llm = llm + nindxm(ll,itype)
                  enddo
                enddo
                ! Right-multiplication (-> head/wing1)
                if(lc<lcut(itype)) then
                  lm0 = sum( [((2*l+1)*nindx(l,itype),l=0,lc)] )
                  l   = lc + 1
                  do mm = -1,1
                    lm = lm0 + (l+mc+mm)*nindx(l,itype)
                    do n = 1,nindx(l,itype)
                      lm                  = lm + 1
                      gnt                 = gaunt(1,lc+1,lc,mm,mc+mm,mc) * clebsch
                      head (:,mm+2,:)     = head (:,mm+2,:)     + qq0(:,lm,:) * integralc(n,1) * gnt * svol
                      wing1(:,mm+2,ieq,:) = wing1(:,mm+2,ieq,:) +  qq(:,lm,:) * integralc(n,1) * gnt * svol
                    enddo
                  enddo
                endif
                if(lc>0) then
                  lm0 = sum( [((2*l+1)*nindx(l,itype),l=0,lc-2)] )
                  l   = lc - 1
                  do mm = -1,1 ; if(abs(mc+mm)>l) cycle
                    lm = lm0 + (l+mc+mm)*nindx(l,itype)
                    do n = 1,nindx(l,itype)
                      lm                  = lm + 1
                      gnt                 = gaunt(1,lc-1,lc,mm,mc+mm,mc) * clebsch
                      head (:,mm+2,:)     = head (:,mm+2,:)     + qq0(:,lm,:) * integralc(n,2) * gnt * svol
                      wing1(:,mm+2,ieq,:) = wing1(:,mm+2,ieq,:) +  qq(:,lm,:) * integralc(n,2) * gnt * svol
                    enddo
                  enddo
                endif
                Time( if(lcore_soc.and.nspin==2) then ; write(*,'(''s:'''NoA) ; call cpu_done(cputime) ; endif )

              enddo

              Time( if(.not.lcore_soc.or.nspin==1) then ; write(*,'(''s:'''NoA) ; call cpu_done(cputime) ; endif )

            endif
          enddo
          enddo

          enddo ! spin loop

          deallocate( p,q,qq,qq0 )

        enddo ! atom equiv

        ! Spin factor 2
        if(spin==0.and.nterm==1) then
          suscep = suscep * 2
          wing1  = wing1  * 2
          wing2  = wing2  * 2
        endif

        ! Transform wings to cartesian coordinates
        do ifrq = 1,nfrq
          do lm = 1,nlm
            wing2(:,lm,:,ifrq) = matmul ( conjg(tocart) , wing2(:,lm,:,ifrq) ) *   img  ! factor img from momentum operator: -img*nabla
            wing1(lm,:,:,ifrq) = matmul (       tocart  , wing1(lm,:,:,ifrq) ) * (-img) ! (and -1 because we differentiated the core state)
          enddo
        enddo

        ! Symmetrize and write to file
        Time( call cpu_time(cputime) )
        write(iunit) neq(itype),nlm
        allocate(s0(nlm,nlm,neq(itype)))
        allocate(w1(nlm,3,neq(itype)),w10(nlm,3,neq(itype)))
        allocate(w2(3,nlm,neq(itype)),w20(3,nlm,neq(itype)),ww(3,nlm))
        do ifrq = 1,nfrq
# ifndef switch_off_symmetry_core
          s0  = suscep(:,:,:,ifrq)
          w10 = wing1 (:,:,:,ifrq)
          w20 = wing2 (:,:,:,ifrq)
          do isym = 2,nsym
            rot = matmul(rlat,matmul(sym(isym)%rrot,transpose(lat)))/(2*pi)
            do ieq = 1,neq(itype)
              ww          = matmul( rot , w20(:,:,ieq)                  )
              w1(:,:,ieq) = matmul(       w10(:,:,ieq) , transpose(rot) )
              if(isym>nsymt) then ; w2(:,:,ieq) = conjg(transpose(w1(:,:,ieq))) ; w1(:,:,ieq) = conjg(transpose(ww))
              else                ; w2(:,:,ieq) = ww
              endif
            enddo
            call mtrafo_mt( wing1(:,:,:,ifrq),w1,  nlm,3,1,isym,itype,1,.true.)
            call mtrafo_mt( wing2(:,:,:,ifrq),w2,  3,nlm,1,isym,itype,2,.true.)
            call mtrafo_mt(suscep(:,:,:,ifrq),s0,nlm,nlm,1,isym,itype,3,.true.)
          enddo
          suscep(:,:,:,ifrq) = suscep(:,:,:,ifrq) / nsym
          wing1 (:,:,:,ifrq) = wing1 (:,:,:,ifrq) / nsym
          wing2 (:,:,:,ifrq) = wing2 (:,:,:,ifrq) / nsym
# endif
          do ieq = 1,neq(itype)
            if(spin<=2) then ; write(iunit) suscep(:,:,ieq,ifrq),wing1(:,:,ieq,ifrq),wing2(:,:,ieq,ifrq)
            else             ; write(iunit) suscep(:,:,ieq,ifrq)
            endif
          enddo
        enddo
        Time( write(*,'(''t:'''NoA) ; call cpu_done(cputime) )

        deallocate( s0,w1,w2,w10,w20,ww )
        deallocate( suscep,wing1,wing2 )

      enddo ! atom type

      if(spin==0.and.nterm==1) head = head * 2 ! spin

      ! Transform head to cartesian coordinates, symmetrize, and write to file
      if(spin<=2) then
        do ifrq = 1,nfrq
c          write(*,*) ifrq
c          write(*,'(6F15.10)') head(:,:,ifrq)
          head(:,:,ifrq) = matmul ( conjg(tocart) , matmul ( head(:,:,ifrq) , transpose(tocart) ) )
          head0          = head(:,:,ifrq)
# ifndef switch_off_symmetry_core
          do isym = 2,nsym
            rot = matmul(rlat,matmul(sym(isym)%rrot,transpose(lat)))/(2*pi)
            if(isym>nsymt) then ; head(:,:,ifrq) = head(:,:,ifrq) + transpose ( matmul(rot,matmul(head0,transpose(rot))) )
            else                ; head(:,:,ifrq) = head(:,:,ifrq) +             matmul(rot,matmul(head0,transpose(rot)))
            endif
          enddo
          head(:,:,ifrq) = head(:,:,ifrq) / nsym
# endif
c          write(*,'(6F15.10)') head(:,:,ifrq)
c          write(*,*) sum(head(:,:,ifrq))
          write(iunit) head(:,:,ifrq)
        enddo
      endif
c      stop

      deallocate(head)
      deallocate(wghtSr,wghtSi)
      deallocate(weight,c,index,nindex)

# ifdef Orthog_core
      rdum = cputime
      call cpu_time(cputime) ; write(6,'(A,F8.2)') '  Timing:',cputime-rdum
# else
      call cpu_done(cputime)
# endif
      call fclose(iunit)
      write(6,'(A)') 'Core susceptibility written to "spex.core".'

# ifdef Orthog_core
      core1 = core1_tmp
      core2 = core2_tmp
# endif

      contains

c     ---------------

c     Returns <nc1|dV|dc2>.
      function vmat(nc1,nc2)
      implicit none
      integer, intent(in) :: nc1,nc2
      real_dp             :: vmat
      real_dp             :: intgrf
      vmat = intgrf( ( core1(:,nc1,lc,itype,1) * core1(:,nc2,lc,itype,1) +
     &                 core2(:,nc1,lc,itype,1) * core2(:,nc2,lc,itype,1) ) * (vmt(:,1,itype,1)-vmt(:,1,itype,2))/2, itype )
      end function vmat

c     ---------------

c     Perform Hilbert transformation for core energy ecor.
      subroutine perform_hilbert(ecor)
      implicit none
      real_dp, intent(in) :: ecor
      integer             :: ifrq
      call getwghtS(wghtSr,wghtSi,frq,nfrq,frqs-ecor+efermi,nfrqs,disorder,0,term(iterm))
      if(spin>=3) then ; wghtSr = -wghtSr ; wghtSi = -wghtSi ; endif
      do ifrq = 1,nfrq
        if(abs(real(frq(ifrq)))<1d-10) then
          q(:,:,ifrq) = unpackmat( matvec(p,wghtSr(:,ifrq)) )
        else
          q(:,:,ifrq) = unpackmat( matvec(p,wghtSr(:,ifrq)) ) + img * unpackmat( matvec(p,wghtSi(:,ifrq)) )
        endif
      enddo
      end subroutine perform_hilbert

c     ---------------

c     Calculate radial integrals
      subroutine calc_radial(nc,evec)
      implicit none
      integer, intent(in) :: nc
      real_dp, intent(in) :: evec(2)
      real_dp             :: cor1(maxgrid),cor2(maxgrid),clebsch(2)
      integer             :: n,nn,l,ll,s
      real_dp             :: intgrf
      if(all(evec==[1d0,0d0])) then
        cor1 = core1(:,nc,lc,itype,ispin1)
        cor2 = core2(:,nc,lc,itype,ispin1)
      else
        s          = 3 - 2*ispin1
        clebsch(1) = sqrt( 0.5d0 * (2*lc+1+s*m) / (2*lc+1) )
        clebsch(2) = sqrt( 0.5d0 * (2*lc+1-s*m) / (2*lc+1) ) * (-s)
        cor1       = matmul( core1(:,nc:nc+1,lc,itype,ispin1) , evec*clebsch )
        cor2       = matmul( core2(:,nc:nc+1,lc,itype,ispin1) , evec*clebsch )
      endif
      ! Calculate radial integrals (-> integralc) for head/wing
      if(lc<lcut(itype)) then
        do n = 1,nindx(lc+1,itype)
          integralc(n,1) = intgrf( ( bas1(:,n,lc+1,itype,ispin2) * cor1 +
     &                               bas2(:,n,lc+1,itype,ispin2) * cor2 ) * rgrid(:,itype) , itype)
        enddo
      endif
      if(lc>0) then
        do n = 1,nindx(lc-1,itype)
          integralc(n,2) = intgrf( ( bas1(:,n,lc-1,itype,ispin2) * cor1 +
     &                               bas2(:,n,lc-1,itype,ispin2) * cor2 ) * rgrid(:,itype) , itype)
        enddo
      endif
      ! Calculate radial integrals (-> integral) for wing/suscep
      do ll = 0,lcutm(itype)
        do l = abs(ll-lc),min(ll+lc,lcut(itype)),2
          do nn = 1,nindxm(ll,itype)
            do n = 1,nindx(l,itype)
              integral(n,nn,l,ll) = intgrf ( basm(:,nn,ll,itype) / rgrid(:,itype) *
     &                                     ( bas1(:,n,l,itype,ispin2) * cor1 +
     &                                       bas2(:,n,l,itype,ispin2) * cor2 ) , itype )
            enddo
          enddo
        enddo
      enddo
      end subroutine calc_radial

c     ---------------

      end

c     ---------------

# if 0

c     Calculates k-independent of susceptibility contribution for extrapolar approximation and stores the result in spex.xsu.
c     Also calculates contractions for self-energy and stores them in spex.xse.

# define OMIT_CORE

# ifdef switch_off_symmetry
#   define NKPTI nkpt
# else
#   define NKPTI nkpti
# endif
      subroutine susceptibility_extra(spin,frq,nfrq,job1)
      use global
      use file
      use wrapper
      use key
      use, intrinsic :: iso_fortran_env
      implicit none
      integer,     intent(in)  :: spin,nfrq
      complex_dp,  intent(in)  :: frq(nfrq)
      type(jobtype)            :: job1
      complex_dp,  allocatable :: cprod(:,:,:,:,:),suscep_mt(:,:,:)
      complex_dp,  allocatable :: contract_mt(:,:,:),suscep0(:,:,:),suscep1(:,:,:)
      MCOMPLEX_dp, allocatable :: suscep_pw(:,:),suscep2(:),suscep3(:)
      MCOMPLEX_dp, pointer     :: contract_pw(:,:)
      MCOMPLEX_dp              :: mcarr(nfrq)
      complex_dp               :: carr(nfrq),cdum,hlp(nbasp,nbasp)
      real_dp                  :: integral(maxindx,0:maxlcut,maxindxm,0:maxlcutm)
      real_dp                  :: bracket(nfrq,NKPTI,bando)
      real_dp                  :: gaunt1,rdum,wght,enediff
      integer,     allocatable :: gptm2(:,:)
      integer                  :: pnt( -2*maxval(abs(gptm(1,:))) : 2*maxval(abs(gptm(1,:))) ,
     &                                 -2*maxval(abs(gptm(2,:))) : 2*maxval(abs(gptm(2,:))) ,
     &                                 -2*maxval(abs(gptm(3,:))) : 2*maxval(abs(gptm(3,:))) )
      integer                  :: iunit,g(3),isym
      integer                  :: nkpts(NKPTI),ispin
      integer                  :: ic,itype,ieq,lc,nc,ll,nn,l,n,mc,mm,m,ibas,ibasc,ibas1,i,j,i1,i2,nbasis
      integer                  :: iband,k,ifrq,ncore,icore,ngptm2,ikpt
      integer                  :: ibandq,ikptq,ispinq
      logical                  :: ldum,same_eigenspace
      real_dp                  :: disorder
      real                     :: cputime
      real_dp                  :: intgrf,gaunt
      call getkey(inp,'DISORDER', disorder, section='SUSCEP',  default=0d0, mine=0d0 )
      if(disorder/=0) disorder = 1/(2*disorder) ! eta = 1/(2*tau)
      if(NKPTI/=nkpti.and.NKPTI/=nkpt) Error('wrong definition of NKPTI.')
      if(NKPTI==nkpt.and.storeibz) Error('STOREIBZ not implemented.')
      ! set up G-point set and pointer pnt
      pnt = 0
      do ikpt = 1,nkpt
        do i = 1,ngptm(ikpt)
          do j = 1,ngptm(ikpt)
            g = gptm(:,pgptm(i,ikpt)) - gptm(:,pgptm(j,ikpt))
            pnt(g(1),g(2),g(3)) = 1
            if(any(g/=0)) pnt(-g(1),-g(2),-g(3)) = 0 ! hermiticity
          enddo
        enddo
      enddo
      ngptm2 = sum(pnt)
      allocate ( gptm2(3,ngptm2) )
      pnt = 0
      k   = 0
      do ikpt = 1,nkpt
        do i = 1,ngptm(ikpt)
          do j = 1,ngptm(ikpt)
            g = gptm(:,pgptm(i,ikpt)) - gptm(:,pgptm(j,ikpt))
            if(pnt(g(1),g(2),g(3))==0.and.pnt(-g(1),-g(2),-g(3))==0) then ! only store one of the pair (g,-g) (hermiticity)
              k                   = k + 1
              gptm2(:,k)          = g
              pnt(g(1),g(2),g(3)) = k
            endif
          enddo
        enddo
      enddo
      ! restart?
      write(6,'(/A)') '---------'
      if(restart>0) then
        inquire(file='spex.xsu',exist=ldum)
        if(ldum) then
          write(6,'(/A)') 'File spex.xsu already exists and will be used.'
          goto 3
        endif
      endif
      write(6,'(/A)') 'Calculation of k-independent extrapolar contribution to susceptibility ...'
      if(ene_extra<=0) Error('ene_extra has wrong value.')
      if(disorder/=0)  Warn('Disorder not implemented for extrapolar approximation.')
      if(spin==0.or.spin==1.or.spin==3) then ; ispin = 1
      else                                         ; ispin = 2
      endif
      do ifrq = 1,nfrq
        if(real(frq(ifrq))/=0.and.imag(frq(ifrq))/=0) then
          Error('extrapolar approximation only implemented for real or purely imaginary frequencies.')
        endif
      enddo
      call cpu_time(cputime)
      write(6,'(A'NoA) '  MT contribution ... '
      ! precalculate frequency-dependent expression (in brackets)
      do k = 1,NKPTI
        nkpts(k) = count(kptp(:nkpt)==k)
      enddo
# ifdef switch_off_symmetry
      nkpts = nsym ! case for testing
# endif
      ncore = 0
      do itype = 1,ntype
        i = sum ( [ ( (2*l+1)*nindxc(l,itype),l=0,lcutc(itype)) ] )
        if(i>ncore) ncore = i
      enddo
      allocate ( cprod (ncore,NKPTI,bando,maxbasp,ncent) )
      allocate ( contract_mt(maxbasp*(maxbasp+1)/2,NKPTI,ncent) )
      allocate ( suscep_mt(nfrq,maxbasp*(maxbasp+1)/2,ncent) )
 1    do k = 1,NKPTI
        do iband = 1,bando
          wght    = wintgr(k,iband,ispin) * nkpts(k) / nsym
          enediff = ene_extra - ene(iband,k,ispin)
          do ifrq = 1,nfrq
            bracket(ifrq,k,iband) = ( 1/(frq(ifrq)-enediff) - 1/(frq(ifrq)+enediff) ) * wght
          enddo
        enddo
      enddo
      ! calculate valence-core wavefunction products (for core contribution)
      cprod = 0
      ic    = 0
      do itype = 1,ntype
        do ieq = 1,neq(itype)
          ic    = ic + 1
          icore = 0
          do lc = 0,lcutc(itype)
            do nc = 1,nindxc(lc,itype)
              ! Calculate radial integrals (-> integral)
              do ll = 0,lcutm(itype)
                do nn = 1,nindxm(ll,itype)
                  do l = abs(ll-lc),min(ll+lc,lcut(itype))
                    do n = 1,nindx(l,itype)
                      integral(n,l,nn,ll) = intgrf ( basm(:,nn,ll,itype) / rgrid(:,itype) *
     &                                      ( core1(:,nc,lc,itype,ispin) * bas1(:,n,l,itype,ispin) +
     &                                        core2(:,nc,lc,itype,ispin) * bas2(:,n,l,itype,ispin) ) , itype )
                    enddo
                  enddo
                enddo
              enddo
              do mc = -lc,lc
                ! Calculate projection < M core | cond > (-> cprod)
                icore = icore + 1
                ibasc = 0
                do ll = 0,lcutm(itype)
                  do mm = -ll,ll
                    do nn = 1,nindxm(ll,itype)
                      ibasc = ibasc + 1
                      m     = mc - mm
                      do l = max(abs(ll-lc),abs(m)),min(ll+lc,lcut(itype))
                        gaunt1 = gaunt(l,lc,ll,m,mc,mm)
                        if(gaunt1/=0) then
                          ibas1 = sum ( [ ((2*j+1)*nindx(j,itype),j=0,l-1) ] ) + (l+m) * nindx(l,itype)
                          do n = 1,nindx(l,itype)
                            ibas1 = ibas1 + 1
                            rdum  = integral(n,l,nn,ll) * gaunt1
                            do iband = 1,bando
                              do k = 1,NKPTI
                                cprod(icore,k,iband,ibasc,ic) = cprod(icore,k,iband,ibasc,ic) +
     &                                                          rdum * conjg(cmt(ibas1,ic,iband,k,ispin))
                              enddo
                            enddo
                          enddo
                        endif
                      enddo
                    enddo
                  enddo
                enddo
              enddo
            enddo
          enddo
        enddo
      enddo
      ! calculate MT correction (contractions + core)
      suscep_mt = 0
      do iband = 1,bando
c        call contraction(contract_mt,contract_pw,gptm2,ngptm2,[(i,i=1,NKPTI)],NKPTI,ispin,iband,1,[.false.,.false.])
        ic = 0
        do itype = 1,ntype
          nbasis = sum ( [ ( (2*l+1)*nindxm(l,itype),l=0,lcutm(itype)) ] )
          ncore  = sum ( [ ( (2*l+1)*nindxc(l,itype),l=0,lcutc(itype)) ] )
          do ieq = 1,neq(itype)
            ic = ic + 1
            l  = 0
            do j = 1,nbasis
              do i = 1,j
                l    = l + 1
                carr = 0
                do k = 1,NKPTI
                  cdum = contract_mt(l,k,ic)
# ifndef OMIT_CORE
                  do icore = 1,ncore
                    cdum = cdum - cprod(icore,k,iband,i,ic) * conjg(cprod(icore,k,iband,j,ic))
                  enddo
# endif
                  do ifrq = 1,nfrq
                    carr(ifrq) = carr(ifrq) + cdum * bracket(ifrq,k,iband)
                  enddo
                enddo
                do ifrq = 1,nfrq
                  suscep_mt(ifrq,l,ic) = suscep_mt(ifrq,l,ic) + carr(ifrq)
                enddo
              enddo
            enddo
          enddo
        enddo
      enddo
      ! loop over spins
      if(spin==0) then
        if(ispin==1) then
          if(nspin==1) then ; suscep_mt = suscep_mt * 2
          else                ; ispin = 2 ; goto 1
          endif
        else
          ispin = 1
        endif
      endif
      call cpu_done(cputime)
      ! symmetrize MT contribution and write to spex.xsu
      write(6,'(A'NoA) '  symmetrize      ... '
      iunit = fopen('spex.xsu',form='unformatted',status='unknown')
      write(iunit) ene_extra
      do ifrq = 1,nfrq
        ic = 0
        do itype = 1,ntype
          nbasis = sum ( [ ( (2*l+1)*nindxm(l,itype),l=0,lcutm(itype)) ] )
# ifndef switch_off_symmetry
          allocate ( suscep0(nbasis,nbasis,neq(itype)) )
          allocate ( suscep1(nbasis,nbasis,neq(itype)) )
          do ieq = 1,neq(itype)
            suscep0(:,:,ieq) = unpackmat(suscep_mt(ifrq,:nbasis*(nbasis+1)/2,ic+ieq))
          enddo
          suscep1 = suscep0
          do isym = 2,nsym
            call trafo_mt(suscep1,suscep0,nbasis,nbasis,isym,itype)
          enddo
          do ieq = 1,neq(itype)
            suscep_mt(ifrq,:nbasis*(nbasis+1)/2,ic+ieq) = packmat(suscep1(:,:,ieq))
          enddo
          deallocate ( suscep0,suscep1 )
# endif
          do ieq = 1,neq(itype)
            ic = ic + 1
            write(iunit) nbasis
            write(iunit) suscep_mt(ifrq,:nbasis*(nbasis+1)/2,ic)
          enddo
        enddo
      enddo
      deallocate ( suscep_mt,cprod )
      deallocate ( contract_mt )
      call cpu_done(cputime)
      ! calculate PW correction
      write(6,'(A'NoA) '  PW contribution ... '
      allocate ( suscep_pw(nfrq,ngptm2) )
      allocate ( contract_pw(ngptm2,NKPTI) )
      suscep_pw = 0
 2    do iband = 1,bando
c        call contraction(contract_mt,contract_pw,gptm2,ngptm2,[(i,i=1,NKPTI)],NKPTI,ispin,iband,2,[.false.,.false.])
        do i = 1,ngptm2
          mcarr = 0
          do k = 1,NKPTI
            do ifrq = 1,nfrq
              mcarr(ifrq) = mcarr(ifrq) + contract_pw(i,k) * bracket(ifrq,k,iband)
            enddo
          enddo
          do ifrq = 1,nfrq
            suscep_pw(ifrq,i) = suscep_pw(ifrq,i) + mcarr(ifrq)
          enddo
        enddo
      enddo
      ! loop over spins
      if(spin==0) then
        if(ispin==1) then
          if(nspin==1) then ; suscep_pw = suscep_pw * 2
          else                ; ispin = 2 ; goto 2
          endif
        else
          ispin = 1
        endif
      endif
      ! symmetrize PW contribution and write to spex.xsu
      write(iunit) ngptm2,pnt
      allocate ( suscep2(ngptm2) )
      allocate ( suscep3(ngptm2) )
      do ifrq = 1,nfrq
# ifndef switch_off_symmetry
        suscep2 = suscep_pw(ifrq,:)
        suscep3 = suscep2
        do isym = 2,nsym
          call trafo_pw(suscep3,suscep2,gptm2,ngptm2,pnt,isym)
        enddo
        suscep_pw(ifrq,:) = suscep3
# endif
        write(iunit) suscep_pw(ifrq,:)
      enddo
      deallocate ( contract_pw )
      deallocate ( suscep_pw )
      call fclose(iunit)
      call cpu_done(cputime)
 3    if(job1%type==J_GW) then
        call cpu_time(cputime)
        write(6,'(/A'NoA) 'Calculation of k-independent contractions for self-energy ...'
        iunit = fopen('spex.xse',form='unformatted',status='unknown')
        allocate ( contract_mt(maxbasp*(maxbasp+1)/2,ncent,1) )
        allocate ( contract_pw(ngptm2,1) )
        do iband = 1,size(job1%band)
          if(iband/=1) then
            if(ikptq==job1%kpt(iband).and.ispinq==job1%spin(iband)) then
              i = job1%band(iband)
              if(same_eigenspace(ibandq,i,ikptq,ispinq)) cycle ! skip degenerate band
            endif
          endif
          ibandq = job1%band(iband)
          ikptq  = job1%kpt(iband)
          ispinq = job1%spin(iband)
c          call contraction(contract_mt,contract_pw,gptm2,ngptm2,[ikptq],1,ispinq,ibandq,3,[.false.,.false.])
          write(iunit) ibandq,ikptq,ispinq
          ic = 0
          do itype = 1,ntype
            nbasis = sum ( [ ( (2*l+1)*nindxm(l,itype),l=0,lcutm(itype)) ] )
            do ieq = 1,neq(itype)
              ic = ic + 1
              write(iunit) nbasis
              write(iunit) contract_mt(:nbasis*(nbasis+1)/2,ic,1)!((contract_mt(i,j,ic,1),i=1,j),j=1,nbasis)
            enddo
          enddo
          write(iunit) ngptm2,pnt
          write(iunit) gptm2
          write(iunit) contract_pw(:ngptm2,1)
        enddo
        deallocate ( contract_mt,contract_pw )
        call fclose(iunit)
        call cpu_done(cputime)
      endif
      deallocate ( gptm2 )
      end

      ! transform to <P(isym) M|...> = <M | P(isym)^-1...>
      subroutine trafo_pw(vec1,vec0,gptm2,ngptm2,pnt,isym)
      use global
      use, intrinsic :: iso_fortran_env
      implicit none
      integer,     intent(in)    :: isym,ngptm2,gptm2(3,ngptm2)
      integer,     intent(in)    :: pnt( -2*maxval(abs(gptm(1,:))) : 2*maxval(abs(gptm(1,:))) ,
     &                                   -2*maxval(abs(gptm(2,:))) : 2*maxval(abs(gptm(2,:))) ,
     &                                   -2*maxval(abs(gptm(3,:))) : 2*maxval(abs(gptm(3,:))) )
      MCOMPLEX_dp, intent(in)    :: vec0(ngptm2)
      MCOMPLEX_dp, intent(inout) :: vec1(ngptm2)
      MCOMPLEX_dp                :: mcdum
      integer                    :: i,j,g(3)
      do i = 1,ngptm2
        g = matmul(sym(isym)%rrot,gptm2(:,i))
        j = pnt(g(1),g(2),g(3))
        if(j==0) then
          j     = pnt(-g(1),-g(2),-g(3)) ; if(j==0) Error('missing pointer.')
          mcdum = MCONJG(vec0(j))
        else
          mcdum = vec0(j)
        endif
        if(isym>nsymt) then
          vec1(i) = vec1(i) + conjg ( mcdum * exp(img*2*pi*dot_product(g,sym(isym)%transl)) )
        else
          vec1(i) = vec1(i) +         mcdum * exp(img*2*pi*dot_product(g,sym(isym)%transl))
        endif
      enddo
      end

      ! transform to <P(isym)M |...| P(isym)M> = < M | P(isym)^-1 ... | M >
      subroutine trafo_mt(matrix1,matrix0,nbasis1,nbasis0,isym,itype)
      use global
      use wrapper
      use, intrinsic :: iso_fortran_env
      implicit none
      integer,    intent(in)    :: isym,nbasis1,nbasis0,itype
      complex_dp, intent(inout) :: matrix1(nbasis1,nbasis1,neq(itype))
      complex_dp, intent(in)    :: matrix0(nbasis0,nbasis0,neq(itype))
      complex_dp                :: matrix(nbasis1,nbasis1)
      integer                   :: ic0,ic,ieq,ieq1,ibas,l,n,nn
      if(nbasis1>nbasis0) Error('wrong dimensions: nbasis1 > nbasis0.')
      ic0 = sum(neq(:itype-1))
      ic  = ic0
      do ieq = 1,neq(itype)
        ic   = ic + 1
        ieq1 = pcent(ic,isym) - ic0
        ibas = 1
        do l = 0,lcutm(itype)
          nn = nindxm(l,itype)
          do n = 1,nn
            matrix(:,ibas:ibas+nn*2*l:nn) = matmat(matrix0(:nbasis1,ibas:ibas+nn*2*l:nn,ieq1),dwgn(-l:l,-l:l,l,isym))
            ibas                          = ibas + 1
          enddo
          ibas = ibas + 2*l*nn
        enddo
        ibas = 1
        do l = 0,lcutm(itype)
          nn = nindxm(l,itype)
          do n = 1,nn
            if(isym>nsymt) then
              matrix1(ibas:ibas+nn*2*l:nn,:,ieq) = matrix1(ibas:ibas+nn*2*l:nn,:,ieq) + conjg( ! note that the "operator" is real
     &                                             matmat(conjg(transpose(dwgn(-l:l,-l:l,l,isym))),matrix(ibas:ibas+nn*2*l:nn,:)) )
            else
              matrix1(ibas:ibas+nn*2*l:nn,:,ieq) = matrix1(ibas:ibas+nn*2*l:nn,:,ieq) +
     &                                             matmat(conjg(transpose(dwgn(-l:l,-l:l,l,isym))),matrix(ibas:ibas+nn*2*l:nn,:))
            endif
            ibas                                 = ibas + 1
          enddo
          ibas = ibas + 2*l*nn
        enddo
      enddo
      end

# endif
