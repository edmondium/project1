c--------------------------------------------------------------------------------
c Copyright (c) 2021 Peter Grünberg Institut, Forschungszentrum Jülich, Germany
c This file is part of SPEX and available as free software under the conditions
c of the MIT license as expressed in the LICENSE file in more detail.
c--------------------------------------------------------------------------------

c Version for FLEUR v26b (2019.03)
c
c ----------------------------------
c
      
# include "cppmacro.h"

      module readwrite

      use global
      use file
      use Hwrapper
      use hdf5
      use, intrinsic :: iso_fortran_env

      character(*), parameter :: dftlong    = 'FLEUR v0.26b (2019.03)'
      character(*), parameter :: dftshort   = 'fleur'
      character(*), parameter :: file_wavef = 'KS.hdf'
c Parameters that must also be present in routine getinput
      logical                 :: invs
      integer, allocatable    :: gpt1(:,:,:)
      integer, allocatable    :: nlo(:),llo(:,:)
      integer                 :: nkpt1,nband0,kpt_order,nldau=0
      real_dp                 :: fac_kptadd,maxene
c Parameters that are needed for file input and output
      integer                 :: neigd,nlod,lcutd
      integer, private        :: ntypd
      logical                 :: invs_inp

      contains

c ----------------------------------

      subroutine read_param
      use util, only: chr
      implicit none
      real_dp              :: rdum
      real_dp, allocatable :: rarr(:)
      integer, allocatable :: iarr(:)      
      integer              :: iunit
      integer              :: version,i,l,idum,l1,l2,l3
      ! Read gwa : LAPW parameters
      iunit = fopen('gwa',form='unformatted')
      read(iunit) version
      if(version>-1) Error('Fleur-GW version is '//Chr(-version)//', but I need version 1 or 2')
      if(version==-1) Warn('Fleur-GW version is 1. Same compiler should be used for Fleur and Spex!')
      read(iunit) nspin,ncent,ntype,l,nlod
      ntypd = ntype
      allocate ( lcut(ntype),ztype(ntype),cent(3,ncent),nlo(ntype),
     &           llo(nlod,ntype),neq(ntype),rarr(ntype) )
      llo = -1
      allocate(iarr(ncent))
      read(iunit) iarr,lcut,(idum,i=1,ntype*(l+1)),rarr,
     &            cent,latpar,lat,rdum,neigd,lcutd,nlo,
     &            (llo(1:nlo(i),i),i=1,ntype)
      ztype = nint(rarr)
      if(version==-1) then
        read(iunit,iostat=i) invs,invs_inp,l_soc
      else
        read(iunit,iostat=i) l1,l2,l3 ; invs=l1/=0;invs_inp=l2/=0;l_soc=l3/=0 ! (*) Different compilers store logicals differently, which can lead to undefined behavior if a logical is read from a file
      endif
      sqa = 0
      if(i/=0) then
        l_soc = .false. ! older version
      else if(version==-2) then
        backspace(iunit)
        read(iunit,iostat=i) l1,l2,l3,sqa ; invs=l1/=0;invs_inp=l2/=0;l_soc=l3/=0
        if(i/=0) sqa = 0
      endif
      neq = [ (count(iarr==i), i=1,ntype) ] ! Obtain number of equivalent atoms per atom type
      deallocate(iarr)                      ! => neq.
      call fclose(iunit)
      ! Read LATTC : gcut (which is for no apparent reason in this file.)
      iunit = fopen('LATTC')
      read(iunit,*) ; read(iunit,*) ; read(iunit,*) ; read(iunit,*) ; read(iunit,*) gcut
      call fclose(iunit)
      end subroutine read_param

c ----------------------------------

      subroutine read_symmetry
      implicit none
      integer :: i,j
      integer :: iunit
      logical :: ldum
      inquire ( file='spex.sym',exist=ldum )
      if(ldum) then ; iunit = fopen('spex.sym',form='formatted')
      else          ; iunit = fopen('sym.out', form='formatted')
      endif
      read(iunit,*) nsym
      allocate ( sym(nsym) )
      do i = 1,nsym
        read(iunit,*)
        read(iunit,*) ( sym(i)%rot(j,:),sym(i)%transl(j), j=1,3 )
        if(any(abs(2*sym(i)%transl-nint(2*sym(i)%transl))>1d-15))
     &    Error('INV and non-sym. vec. is not half a reciprocal lattice vector.')
      enddo
      call fclose(iunit)
      end subroutine read_symmetry

c ----------------------------------

      subroutine read_radial
      implicit none
      real_dp, allocatable :: rarr(:)
      integer              :: i,l,ispin,itype,nrec
      integer              :: iunit,indx(0:maxlcut)
      iunit    = fopen('radfun',form='unformatted')
      read(iunit) grid(:)%number
      maxgrid = maxval(grid(:)%number)
      allocate ( bas1(maxgrid,maxindx,0:maxlcut,ntype,nspin2),
     &           bas2(maxgrid,maxindx,0:maxlcut,ntype,nspin2),
     &           ubas (maxindx,0:maxlcut,ntype,nspin2),
     &           dubas(maxindx,0:maxlcut,ntype,nspin2),
     &           ebas( maxindx,0:maxlcut,ntype,nspin2) )
      allocate ( rarr(maxgrid) )
      do ispin = 1,nspin
        do itype = 1,ntype
          read(iunit) rarr(1:grid(itype)%number),nrec
          if(nrec/=lcut(itype)) Warn('lnonsph <> lmax in Fleur. Wavefunctions are not orthonormal.')
          indx                  = 0
          grid(itype)%first     = rarr(1)                  !
          grid(itype)%increment = log(rarr(2)/rarr(1))     ! Define grid
          grid(itype)%radius    = rarr(grid(itype)%number) !
          do l = 0,nrec
            if(l<=lcut(itype)) then
              read(iunit)
     &          bas1(1:grid(itype)%number,indx(l)+1,l,itype,ispin),
     &          bas1(1:grid(itype)%number,indx(l)+2,l,itype,ispin),
     &          bas2(1:grid(itype)%number,indx(l)+1,l,itype,ispin),
     &          bas2(1:grid(itype)%number,indx(l)+2,l,itype,ispin),
     &          ubas(indx(l)+1,l,itype,ispin),dubas(indx(l)+1,l,itype,ispin),
     &          ubas(indx(l)+2,l,itype,ispin),dubas(indx(l)+2,l,itype,ispin),
     &          ebas(indx(l)+1,l,itype,ispin)
              ebas(indx(l)+2,l,itype,ispin) = ebas(indx(l)+1,l,itype,ispin)
              indx(l) = 2
            else
              read(iunit)
            endif
          enddo
          do i = 1,nlo(itype)
            l       = llo(i,itype)
            indx(l) = indx(l) + 1
            read(iunit)
     &        bas1(1:grid(itype)%number,indx(l),l,itype,ispin),
     &        bas2(1:grid(itype)%number,indx(l),l,itype,ispin),
     &        ubas(indx(l),l,itype,ispin),dubas(indx(l),l,itype,ispin),
     &        ebas(indx(l),l,itype,ispin)
          enddo
        enddo
      enddo
      call fclose(iunit)
      deallocate ( rarr )
# ifdef IPOLTEST
      allocate ( rgrid(maxgrid,ntype) )
      do itype = 1,ntype
        rgrid(1,itype) = grid(itype)%first
        do i = 2,grid(itype)%number
          rgrid(i,itype) = rgrid(i-1,itype) * exp(grid(itype)%increment)
        enddo
      enddo
# endif      
      end subroutine read_radial

c ----------------------------------

      subroutine read_core
      implicit none
      real_dp, allocatable :: rarr(:),rarr2(:)
      real_dp              :: wcore1,wcore2,ecore1,ecore2
      integer              :: iarr(0:maxlcut)
      integer              :: ncore
      integer              :: i,l,ll,ispin,itype,idum
      integer              :: iunit
      ! first determine maximal l-number of cores (lcutc,nindxc), then read core states
      lcore_soc =  .true. ! only for testing: set to .false. to test reading of non-relativistic core states
      iunit     =  fopen('ecore',form='unformatted')
      maxindxc  =  0
      lcutc     = -1
      do ispin = 1,nspin
        do itype = 1,ntype
          read(iunit) ncore
          iarr = 0
          i    = 0
          do while(i<ncore)
            i = i + 1
            read(iunit) l
            if(l>maxlcut) Error('l(core) > l(basis).')
            if(l>lcutc(itype)) lcutc(itype) = l
            iarr(l) = iarr(l) + 1
            if(l>=1) then
              i = i + 1
              read(iunit) ll
              if(ll/=l) Error('SO splitted state not found.')
              if(lcore_soc) iarr(l) = iarr(l) + 1
            endif
            if(iarr(l)>maxindxc) maxindxc = iarr(l)
          enddo
        enddo
      enddo
      maxlcutc = maxval(lcutc)
      allocate ( nindxc(0:maxlcutc,ntype),rarr(maxgrid),rarr2(maxgrid),
     &           ecore(maxindxc,0:maxlcutc,ntype,nspin),
     &           core1(maxgrid,maxindxc,0:maxlcutc,ntype,nspin),
     &           core2(maxgrid,maxindxc,0:maxlcutc,ntype,nspin))
      core1 = 0
      core2 = 0
      rewind(iunit)
      do ispin = 1,nspin
        nindxc = 0 ! It is assumed that the numbers of spin-up and spin-down core states are equal.
        do itype = 1,ntype
          read(iunit) ncore
          i = 0
          do while (i<ncore)
            i = i + 1
            read(iunit) l
            nindxc(l,itype) = nindxc(l,itype) + 1
            idum            = nindxc(l,itype)
            if(nindxc(l,itype)>maxindxc) Bug('maximal core function index exceeded.')
            backspace(iunit)
            read(iunit) ll,
     &        wcore1,ecore1,
     &        core1(1:grid(itype)%number,idum,l,itype,ispin),
     &        core2(1:grid(itype)%number,idum,l,itype,ispin)
            if(ll==0) then
              ecore(idum,l,itype,ispin) = ecore1
            else
              i = i + 1
              if(lcore_soc) then
                ecore(idum,l,itype,ispin) = ecore1
                nindxc(l,itype)           = nindxc(l,itype) + 1
                idum                      = nindxc(l,itype)
                read(iunit) ll,
     &            wcore2,ecore(idum,l,itype,ispin),
     &            core1(1:grid(itype)%number,idum,l,itype,ispin),
     &            core2(1:grid(itype)%number,idum,l,itype,ispin)
              else
                read(iunit) ll,
     &            wcore2,ecore2,
     &            rarr (1:grid(itype)%number),
     &            rarr2(1:grid(itype)%number)
                ecore(idum,l,itype,ispin) =
     &            ( wcore1 * ecore1 + wcore2 * ecore2 ) /
     &            ( wcore1 + wcore2 )
                core1(1:grid(itype)%number,idum,l,itype,ispin) =
     &            ( wcore1 * core1(1:grid(itype)%number,idum,l,itype,ispin)
     &            + wcore2 *  rarr(1:grid(itype)%number) ) /
     &            ( wcore1 + wcore2 )
                core2(1:grid(itype)%number,idum,l,itype,ispin) =
     &            ( wcore1 * core2(1:grid(itype)%number,idum,l,itype,ispin)
     &            + wcore2 * rarr2(1:grid(itype)%number) ) /
     &            ( wcore1 + wcore2 )
              endif
            endif
          enddo
        enddo
      enddo
      call fclose(iunit)
      if(lcore_soc) then
        do ispin = 1,nspin ! Order: 2p3/2 2p1/2 3p3/2 3p1/2 ...
          do itype = 1,ntype
            do l = 1,lcutc(itype)
              do i = 1,nindxc(l,itype),2
                ecore1                     = ecore(  i,  l,itype,ispin)
                ecore(  i,  l,itype,ispin) = ecore(  i+1,l,itype,ispin)
                ecore(  i+1,l,itype,ispin) = ecore1
                rarr                       = core1(:,i,  l,itype,ispin)
                rarr2                      = core2(:,i,  l,itype,ispin)
                core1(:,i,  l,itype,ispin) = core1(:,i+1,l,itype,ispin)
                core2(:,i,  l,itype,ispin) = core2(:,i+1,l,itype,ispin)
                core1(:,i+1,l,itype,ispin) = rarr
                core2(:,i+1,l,itype,ispin) = rarr2
              enddo
            enddo
          enddo
        enddo
      endif
      deallocate (rarr,rarr2)
      end subroutine read_core

c ----------------------------------

      subroutine read_gpt
      use hdf5
      use Hwrapper
      implicit none
      integer        :: dim(2),Herr
      character(8)   :: ckpt
      real_dp        :: rarr(3),rdum
      integer        :: i,ispin,ikpt
      logical        :: first,ldum
      integer(HID_T) :: Hfile,Hkpt
# define ERRSTOP if(Herr/=0) Error('Fatal HDF5 error.')
      call kpt_reorder1(0)
      call hdf_fopen(Hfile,'KS.hdf',0)
      call hdf_rdwr_a(Hfile,'version',0,i) ; if(i/=1) Error('Incompatible KS.hdf file.')
      first = .true.
 1    ispin = 1
      do ikpt = 1,nkpt1
        write(ckpt,'(I1,''-'',I6.6)') ispin,ikpt
        call h5lexists_f(Hfile,ckpt,ldum,Herr) ; ERRSTOP
        if(.not.ldum) Error('K data set '//ckpt//' not found in KS.hdf.')
        call h5gopen_f(Hfile,ckpt,Hkpt,Herr) ; ERRSTOP
        if(first) then
          call hdf_rdwr_a(Hkpt,'kpt',       0,rarr)
          call hdf_rdwr_a(Hkpt,'ngpt,nband',0,dim) ; ngpt(ikpt) = dim(1)
          rdum = maxval(abs(modulo1(rarr-kpt(:,ikpt))))
          if(rdum>1d-12) then
            write(0,'(A,I4,A)') 'read_gpt: k-point',ikpt,' inconsistent'
            write(0,'(7X,3F10.5)') kpt(:,ikpt)
            write(0,'(7X,3F10.5)') rarr
            Error('k-point sets inconsistent (check BZ and KS.hdf).')
          endif
        else
          call hdf_rdwr(Hkpt,'gpt',0,gpt1(:,:ngpt(ikpt),ikpt),0) ; gpt1(:,ngpt(ikpt)+1:,ikpt) = 0
        endif
        call h5gclose_f(Hkpt,Herr) ; ERRSTOP
      enddo
      if(first) then
        allocate(gpt1(3,maxval(ngpt(:nkpt1)),nkpt2))
        gpt1  = 0
        first = .false.
        goto 1
      endif
      call hdf_fclose(Hfile)
      call kpt_reorder2(0)
      end subroutine read_gpt

c ----------------------------------

      subroutine read_ene
      use hdf5
      use Hwrapper
      implicit none
      integer        :: k,ispin
      logical        :: first
      integer(HID_T) :: Hfile,Hkpt
      integer        :: dim(2),Herr
      character(8)   :: ckpt
      call kpt_reorder1(0)
      call hdf_fopen(Hfile,'KS.hdf',0)
      first    = .true.
 1    do ispin = 1,nspin1
        do k = 1,nkpt1
          write(ckpt,'(I1,''-'',I6.6)') ispin,k
          call h5gopen_f(Hfile,ckpt,Hkpt,Herr) ; ERRSTOP
          call hdf_rdwr_a(Hkpt,'ngpt,nband',0,dim) ; nband(k,ispin) = dim(2)
          if(.not.first)
     &      call hdf_rdwr(Hkpt,'energies',0,ene(:nband(k,ispin),k,ispin),0)
          call h5gclose_f(Hkpt,Herr) ; ERRSTOP
        enddo
      enddo
      if(first) then
        first   = .false.
        maxband = maxval(nband)
        allocate(ene(maxband,nkpt1,nspin1))
        goto 1
      endif
      call hdf_fclose(Hfile)
      call kpt_reorder2(0)
      end subroutine read_ene

c ----------------------------------

      subroutine write_kpt
      implicit none
      real_dp :: scale
      integer :: iunit
      integer :: i,mult
      integer :: kgv
      call kpt_reorder1(0)
      iunit = fopen('kpts',status='unknown')
      scale = kgv(nkpt3,3) ; scale = scale * fac_kptadd
      write(iunit,'(I5,F20.10)') nkpt1,scale
      do i = 1,nkpt1
        if(i<=nkpti) then ; mult = count(kptp(:nkpt)==i)
        else              ; mult = 0
        endif
        write(iunit,'(4F10.5)') kpt(:,i)*scale,mult*1d0
      enddo
      call fclose(iunit)
      call kpt_reorder2(0)
      write(0,'(A)') 'File kpts written.'
      end subroutine write_kpt

      subroutine write_qpt(qpt,nqpt,label)
      implicit none
      integer,      intent(in) :: nqpt
      real_dp,      intent(in) :: qpt(3,nqpt)
      character(*), intent(in) :: label(nqpt)
      end subroutine write_qpt

c ----------------------------------

      subroutine read_pot(mode)
      implicit none
      integer,     intent(in)  :: mode
      complex_dp,  allocatable :: vstar_c(:)
      real_dp,     allocatable :: vstar_r(:)
      real_dp,     allocatable :: phasestar(:,:,:)
      real_dp                  :: rdum,gmax
      integer,     allocatable :: pointstar(:,:,:)
      integer                  :: nstar,nstar2,k1d,k2d,k3d,maxmlh
      integer                  :: iunit,ispin,itype,ieq,ic,idum,ilh,istar
      integer                  :: i,j,k,isym
      logical                  :: read_tot,doloop,lfirst
      character(2)             :: chdum
      if(mode==0) Error('mode=0 disabled.')
      if(all(mode/=[1,2])) Bug('Mode for read_pot unknown.')
      lfirst = .not.allocated(llh)
      if(mode==1.or.mode==2) then
        ! read stars
        iunit = fopen('stars',form='unformatted',status='old')
        read(iunit) gmax,(idum,i=1,11),nstar,nstar2,k1d,k2d,k3d
        if(invs_inp) then ; allocate ( vstar_r(nstar) )
        else              ; allocate ( vstar_c(nstar) )
        endif
        allocate ( phasestar(-k1d:k1d,-k2d:k2d,-k3d:k3d) )
        allocate ( pointstar(-k1d:k1d,-k2d:k2d,-k3d:k3d) )
        if(lfirst) then
          allocate ( vpw   (-k1d:k1d,-k2d:k2d,-k3d:k3d,nspin) )
          allocate ( vpw_xc(-k1d:k1d,-k2d:k2d,-k3d:k3d,nspin) )
        endif
        read(iunit) (idum,i=1,nstar+nstar2),phasestar,
     &    (rdum,i=1,  nstar+  nstar2),
     &    (idum,i=1,3*nstar+2*nstar2),
     &    pointstar
        call fclose(iunit)
      else if(mode/=0) then
        Error('unknown mode!')
      endif
      ! read lattice harmonics
      if(lfirst) then
        iunit = fopen('latharm',form='unformatted',status='old')
        read(iunit) maxlh,maxmlh
        maxlh = maxlh + 1 ! Fleur starts at 0 for nlh
        allocate ( nlh(ntype),llh(maxlh,ntype),nmlh(maxlh,ntype),
     &             mlh(maxmlh,maxlh,ntype),clh(maxmlh,maxlh,ntype) )
        ic = 0
        do itype = 1,ntype
          read(iunit) nlh(itype),llh(:nlh(itype)+1,itype),nmlh(:nlh(itype)+1,itype),
     &                mlh(:maxmlh,:nlh(itype)+1,itype),clh(:maxmlh,:nlh(itype)+1,itype)
          do ieq = 1,neq(itype)
            ic = ic + 1
            read(iunit) isym ; if(isym<0) isym = symtab(invsym,-isym)
            if(pcent(sum(neq(:itype-1))+1,sym(isym)%inv)/=ic) Error('symcent error.')            
          enddo
        enddo
        nlh = nlh + 1     ! Fleur starts at 0 for nlh
        call fclose(iunit)
        if(mode==0) then
          allocate ( vmt(maxgrid,1,ntype,nspin) )
        else
          allocate ( vmt   (maxgrid,maxlh,ntype,nspin) )
          allocate ( vmt_xc(maxgrid,maxlh,ntype,nspin) )
        endif
      endif
      ! read potential
      if(mode==0) then
        iunit    = fopen('pottot', form='unformatted',status='old')
        read_tot = .true.
        doloop   = .false.
      else if(mode==1) then
        iunit    = fopen('potcoul', form='unformatted',status='old')
        read_tot = .false.
        doloop   = .true.
      else if(mode==2) then
        iunit    = fopen('potx', form='unformatted',status='old')
        read_tot = .true.
        doloop   = .true.
      endif
 1    read(iunit)
      read(iunit)
      do ispin = 1,nspin
        read(iunit)
        read(iunit) itype
        if(itype/=ntype) Error('wrong number of atom types in pot file.')
        do itype = 1,ntype
          read(iunit) chdum,idum,idum
          if(idum/=grid(itype)%number) Error('wrong number of grid points in pot file.')
          read(iunit) idum,idum
          if     (idum<nlh(itype)-1) then
            Warn('Less lattice harmonics needed (pot file) than defined in file latharm.')
            nlh(itype) = idum + 1
          else if(idum>nlh(itype)-1) then ;
            Error('More lattice harmonics needed (pot file) than defined in file latharm.')
          endif
          do ilh = 1,nlh(itype)
            read(iunit)
            if(mode==0.and.ilh/=1) then ; read(iunit)
            else                        ; read(iunit) vmt(:grid(itype)%number,ilh,itype,ispin)
            endif
          enddo
          ! pottot  : spherical term includes the factor Y_00=1/sqrt(4pi) and is multiplied with r 
          ! potx    : spherical term does NOT include the factor Y_00 but is multiplied with r
          ! potcoul : spherical term does NOT include the factor Y_00 and is NOT multiplied with r
          if(read_tot) then
            vmt(:,1,itype,ispin) = vmt(:,1,itype,ispin) * sqrt(4*pi)
          else
# ifdef IPOLTEST            
            vmt(:,1,itype,ispin) = vmt(:,1,itype,ispin) * rgrid0(:,itype)
# else
            vmt(:,1,itype,ispin) = vmt(:,1,itype,ispin) * rgrid(:,itype)
# endif            
          endif
        enddo
        if(mode/=0) then
          read(iunit) idum
          if(idum/=nstar) Error('Number of stars in files pottot, potcoul, and stars does not agree.')
          if(invs_inp) then ; read(iunit) vstar_r
          else              ; read(iunit) vstar_c
          endif
          vpw(:,:,:,ispin) = 0
          do k = -k3d,k3d
            do j = -k2d,k2d
              do i = -k1d,k1d
                istar = pointstar(i,j,k)
                if(istar/=0.and.mode/=0) then
                  if(sum(matmul(rlat,[i,j,k])**2)>gmax**2+1d-6)
     &              Bug('G point outside gmax sphere but pointer is nonzero.')
                  if(invs_inp) then ; vpw(i,j,k,ispin) = phasestar(i,j,k) * vstar_r(istar)
                  else              ; vpw(i,j,k,ispin) = phasestar(i,j,k) * vstar_c(istar)
                  endif
                endif
              enddo
            enddo
          enddo
        else
          read(iunit)
          read(iunit)
        endif
      enddo
      call fclose(iunit)
      if(doloop) then
        vmt_xc   = vmt
        vpw_xc   = vpw
        iunit    = fopen('pottot',form='unformatted',status='old')
        read_tot = .true.
        doloop   = .false.
        goto 1
      endif
      if(mode==1) then
        vmt_xc = vmt - vmt_xc
        vpw_xc = vpw - vpw_xc
      else if(mode==2) then
        vmt_xc(:,1,:,:) = vmt_xc(:,1,:,:) / sqrt(4*pi) ! undo multiplication (see above)
      endif
      if(mode/=0) then
        deallocate ( phasestar,pointstar )
        if(invs_inp) then ; deallocate ( vstar_r )
        else              ; deallocate ( vstar_c )
        endif
      endif
      end subroutine read_pot

c ----------------------------------

      subroutine read_wannier
      implicit none
      integer                 :: ikpt,ispin
      integer                 :: iunit,ix,iy,iz,i,j,iwan,ikpt1,iband
      real_dp                 :: r1,r2
      logical                 :: found
      if(any(modulo1r(2*sym(symkpt(ikpt))%transl)/=0))
     &  Error('Translation vector coefficient neither 0.0 nor 0.5. Fleur cannot do that!(?)')
      do ispin = 1,nspin2
        if(ispin==1) then
          iunit = fopen('WF1.umn',status='old')
        else
          iunit = fopen('WF2.umn',status='old')
        endif        
        do ikpt = 1,nkpt
          rewind(iunit)
          read(iunit,*,end=2,err=3)
          read(iunit,*,end=2,err=3) iband,ikpt1,iwan
          if(iband/=nwanband) Error('number of bands in WF?.umn_k are insonsistent with proj')
          if(ikpt1/=nkpt    ) Error('number of k points in WF?.umn_k are insonsistent with '//trim(inpfile))
          if(iwan /=nwan    ) Error('number of Wannier functions in WF?.umn_k are insonsistent with proj')
          if(.not.allocated(uwan)) then
            allocate(uwan(nwanband,nwan,nkpt,nspin2))
            Error('Check for wanbandi and wanbandf in read_wannier!')
          endif
          found = .false.
          do ix = 0,nkpt3(1)-1
            do iy = 0,nkpt3(2)-1
              do iz = 0,nkpt3(3)-1
                if(pkpt(ix,iy,iz)/=ikpt) then
                  do i = 1,nwan*nwanband ; read(iunit,*,end=2,err=3) ; enddo
                else
                  found = .true.
                  do i = 1,nwan
                    do j = 1,nwanband
                      read(iunit,*,end=2,err=3) iband,iwan,ikpt1,r1,r2
                      if(iband/=j.or.iwan/=i) Error('wrong index.')
                      uwan(j,i,ikpt,ispin) = r1 + img*r2
                    enddo
                  enddo
                  goto 1
                endif
              enddo
            enddo
          enddo
 1        if(.not.found) Error('Did not find k-point index in file WF?.umn_k')
        enddo
        call fclose(iunit)
      enddo
      return
 2    Error('Reached end of file WF?.umn_k')
 3    Error('Read error from file WF?.umn_k')
      end subroutine read_wannier

c ----------------------------------

      subroutine read_vxc(vxcmat,band,nb,ikpt_in,ispin0,l_qsgw)
      use, intrinsic :: iso_fortran_env
      implicit none
      integer,     intent(in)  :: nb,band(nb),ikpt_in,ispin0
      logical,     intent(in)  :: l_qsgw
      MCOMPLEX_dp, intent(out) :: vxcmat(nb,nb)
      integer                  :: iunit
      integer                  :: ispin,ikpt,ikpt0
      integer                  :: i,j,k,nband1(2),n,i1,i2,ib1,ib2
      MCOMPLEX_dp, allocatable :: vxcfull(:,:)
      real_dp,     allocatable :: vxchelp(:,:)
      complex_dp,  allocatable :: socvec(:,:)
      call kpt_reorder1(3)
      if     (ikpt_in<=nkpti)  then ; ikpt0 = ikpt_in
      else if(ikpt_in==nkpt+1) then ; ikpt0 = nkpti + 1
      else                          ; Error('k-point index out of range.')
      endif
      if(l_soc) then
        iunit = fopen('SOCVEC', form='unformatted')
        n = 0
        do ikpt = 1,nkpt1
          if(ikpt/=ikpt0) then
            read(iunit)
          else
            allocate ( socvec(2*neigd,2*neigd) )
            read(iunit) n,n,j,nband1,k,socvec
            if(n/=sum(nband1)) Error('transformation matrix in SOCVEC not quadratic.')
            if(j/=nspin)       Error('number of spins in SOCVEC inconsistent.')
            if(k/=2*neigd)     Error('array dimension in SOCVEC inconsistent.')
            call fclose(iunit)
            exit
          endif
        enddo
        if(n==0)           Error('Not a parent k point.')
        if(maxval(band)>n) Error('band indices out of range.')
      else
        nband1(:nspin1) = nband(ikpt0,:)
        if(maxval(band)>nband1(ispin0)) Error('band indices out of range.')
      endif
      vxcmat = 0
      iunit  = fopen('vxcfull' ,form='unformatted')
      do ispin = 1,nspin
        allocate ( vxcfull(nband1(ispin),nband1(ispin)) )
# ifndef INV
        if(invs) then
          allocate ( vxchelp(nband1(ispin),nband1(ispin)) )
          vxchelp = 0
        endif
# endif
        vxcfull = 0
        do ikpt = 1,nkpt1
          if(ikpt/=ikpt0.or.ispin/=ispin0) then
            read(iunit)
          else
            if(allocated(vxchelp)) then
              read(iunit) ((vxchelp(i,j),i=1,j),j=1,nband1(ispin))
              vxcfull = vxchelp
            else
              read(iunit) ((vxcfull(i,j),i=1,j),j=1,nband1(ispin))
            endif
            do j = 1,nband1(ispin)
              do i = 1,j-1
                vxcfull(j,i) = MCONJG(vxcfull(i,j))
              enddo
            enddo
            if(l_soc) then
              i = sum(nband1(:ispin-1)) + 1
              j = sum(nband1(:ispin))
              do i2 = 1,nb
                do i1 = 1,nb
                  ib1           = band(i1)
                  ib2           = band(i2)
                  vxcmat(i1,i2) = vxcmat(i1,i2) + dot_product(socvec(i:j,ib1),matmul(vxcfull,socvec(i:j,ib2)))
                enddo
              enddo
              if(nspin==1) then
                i = sum(nband1(:ispin)) + 1
                j = sum(nband1(:ispin)) * 2
                do i2 = 1,nb
                  do i1 = 1,nb
                    ib1           = band(i1)
                    ib2           = band(i2)
                    vxcmat(i1,i2) = vxcmat(i1,i2) + dot_product(socvec(i:j,ib1),matmul(vxcfull,socvec(i:j,ib2)))
                  enddo
                enddo
              endif
            else
              vxcmat = vxcfull(band,band)
            endif
            if(.not.l_soc.or.ispin==nspin) then
              deallocate ( vxcfull )
              if(l_soc) deallocate ( socvec )
              goto 1
            endif
          endif
        enddo
        deallocate ( vxcfull )
        if(allocated(vxchelp)) deallocate ( vxchelp )
      enddo
 1    call fclose(iunit)
      call kpt_reorder2(3)
      if(l_qsgw) call read_qsgw(vxcmat,band,nb,ikpt_in,ispin0)
      end subroutine read_vxc

c ----------------------------------

      subroutine read_qsgw(vxcmat,band,nb,ikpt_in,ispin0)
      use, intrinsic :: iso_fortran_env
      implicit none
      integer,     intent(in)    :: nb,band(nb),ikpt_in,ispin0
      MCOMPLEX_dp, intent(inout) :: vxcmat(nb,nb)
      integer                    :: iunit
      integer                    :: ispin,ikpt,ikpt0
      integer                    :: i,j,n
      MCOMPLEX_dp, allocatable   :: sigfull(:,:)
      call kpt_reorder1(3)
      ikpt0 = ikpt_in
      if(ikpt0==nkpt+1) ikpt0 = nkpti + 1
      iunit = fopen('qsgw',form='unformatted')
      do ispin = 1,nspin
        do ikpt = 1,nkpt1
          if(ikpt/=ikpt0.or.ispin/=ispin0) then
            read(iunit)
            read(iunit)
          else
            read(iunit) n
            if(n<=maxval(band)) Error('Did not find enough bands in qsgw.')
            allocate ( sigfull(n,n) )
            read(iunit) ((sigfull(i,j),i=1,j),j=1,n)
            do j = 1,n
              do i = 1,j-1
                sigfull(j,i) = MCONJG(sigfull(i,j))
              enddo
            enddo
            vxcmat = vxcmat + sigfull(band,band)
            deallocate ( sigfull )
            goto 1
          endif
        enddo
      enddo
 1    call fclose(iunit)
      call kpt_reorder2(3)
      end subroutine read_qsgw

c ----------------------------------

      subroutine write_param(lsoc)
      implicit none
      logical, intent(in) :: lsoc
      integer             :: iunit,i,itype,l,l1,l2,l3
      iunit = fopen('gwa',form='unformatted',status='unknown')
      write(iunit) -2
      write(iunit) nspin,ncent,ntype,maxlcut,nlod
      write(iunit) ((itype,i=1,neq(itype)),itype=1,ntype),lcut,((l+1,l=0,maxlcut),itype=1,ntype),1d0*ztype,cent,latpar,lat,vol,
     &             neigd,lcutd,nlo,(llo(1:nlo(i),i),i=1,ntype)
      l1 = 0 ; if(invs)     l1 = 1
      l2 = 0 ; if(invs_inp) l2 = 1
      l3 = 0 ; if(lsoc)     l3 = 1
      write(iunit) l1,l2,l3,sqa
      call fclose(iunit)
      end subroutine write_param

c ----------------------------------

      subroutine write_wavef(ikpt,ispin,ngpt1,gpt1,eig,cmtin,cpwin,cpwin_c)
      use global
      use Hwrapper
      use hdf5
      implicit none
      integer,     intent(in)           :: ikpt,ispin,ngpt1
      integer,     intent(in)           :: gpt1(3,ngpt1)
      real_dp,     intent(in)           :: eig(:)
      complex_dp,  intent(in)           :: cmtin(maxlmindx,ncent,size(eig),nspin3)
      MCOMPLEX_dp, intent(in), optional :: cpwin(maxgpt,size(eig),nspin3)
      complex_dp,  intent(in), optional :: cpwin_c(maxgpt,size(eig),nspin3)
      integer(HID_T)                    :: Hfile,Hkpt,Horb
      integer                           :: nbnd
      integer                           :: Herr
      character(8)                      :: ckpt
      if(present(cpwin).eqv.present(cpwin_c)) Bug('cpwin and cpwin_c both present or missing.')
      if(ngpt1>maxgpt)                        Bug('ngpt1 > maxgpt.')
      nbnd  = size(eig)
      if(ikpt==1.and.ispin==1) then ; call hdf_fopen(Hfile,'KS.hdf',2) ; call hdf_rdwr_a(Hfile,'version',1,1)
      else                          ; call hdf_fopen(Hfile,'KS.hdf',1)
      endif
      write(ckpt,'(I1,''-'',I6.6)') ispin,ikpt
      call h5gcreate_f(Hfile,ckpt,Hkpt,Herr) ; ERRSTOP
      call hdf_rdwr_a(Hkpt,'kpt',       1,kpt(:,ikpt))
      call hdf_rdwr_a(Hkpt,'ngpt,nband',1,[ngpt1,nbnd])
      call hdf_rdwr  (Hkpt,'gpt',       1,gpt1(:,:ngpt1))
      call hdf_rdwr  (Hkpt,'energies',  1,eig)
      call h5gcreate_f(Hkpt,'orbitals',Horb,Herr) ; ERRSTOP
      if(present(cpwin)) then ; call hdf_rdwr(Horb,'pw',1,cpwin)
      else                    ; call hdf_rdwr(Horb,'pw',1,cpwin_c)
      endif
      call hdf_rdwr(Horb,'mt',1,cmtin)
      call h5gclose_f(Horb,Herr) ; ERRSTOP
      call h5gclose_f(Hkpt,Herr) ; ERRSTOP
      call hdf_fclose(Hfile)
      end subroutine write_wavef

c ----------------------------------

# ifndef LOAD

c ----------------------------------

      subroutine read_wavef
      use util, only: chr
      Mpi( use Mwrapper )
      use, intrinsic :: iso_fortran_env
      implicit none
      integer                  :: ispin,k,ng MpiC(Merr)
      integer(HID_T)           :: Hfile
      character(8)             :: ckpt
      complex_dp,  allocatable :: cmthlp(:,:,:,:)
      MCOMPLEX_dp, allocatable :: cpwhlp(:,:,:)
      call kpt_reorder1(1)
      call hdf_fopen(Hfile,'KS.hdf',0)
      Hpos = 0 ; if(l_soc) Hpos = 1
      ! Read states
      do ispin = 1,nspin1
        do k = 1,nkpti+nkpti2
          write(ckpt,'(I1,''-'',I6.6)') ispin,k
          ng = ngpt(k)
          ! read file
          Allocate_(cmthlp,(maxlmindx,ncent,nband(k,ispin),nspin3))
          call hdf_rdwr(Hfile,'/'//ckpt//'/orbitals/mt',0, cmthlp,0)
          if(l_soc) then
            ifO cmt(:,:,:nband(k,ispin),k,:)     = cmthlp
          else
            ifO cmt(:,:,:nband(k,ispin),k,ispin) = cmthlp(:,:,:,1)
          endif
          Deallocate_(cmthlp)
          Allocate_(cpwhlp,(ng,nband(k,ispin),nspin3))
          call hdf_rdwr(Hfile,'/'//ckpt//'/orbitals/pw',0, cpwhlp,0)
          if(l_soc) then
            ifO cpw(:ng,:nband(k,ispin),k,:)     = cpwhlp
          else
            ifO cpw(:ng,:nband(k,ispin),k,ispin) = cpwhlp(:,:,1)
          endif
          Deallocate_(cpwhlp)
          call checkmem('read:'//Chr(k)//','//Chr(ispin),0d0)
        enddo
      enddo
      Nfence(cmt)
      Nfence(cpw)
      call hdf_fclose(Hfile)
      call kpt_reorder2(1)
      end subroutine read_wavef

c ----------------------------------

# else

c ----------------------------------

      subroutine read_wavef(band,kpt1,ispin,cmtout,cpwout,str)
      use, intrinsic :: iso_fortran_env
      implicit none
      integer,     intent(in)            :: band(:),kpt1(:),ispin,str
      complex_dp,  intent(out), optional :: cmtout(:,:,:,:,:)
      MCOMPLEX_dp, intent(out), optional :: cpwout(:,:,:,:)
      logical                            :: lmt,lpw,lbcnt
      integer                            :: i,k,ib,b1,b2,nb,ikpt1
      integer(HID_T)                     :: Hfile
      character(8)                       :: ckpt
      if(str/=1 Mpi(.and.str/=Msize) ) Bug('Wrong stride (str) argument.')
      call kpt_reorder1(3)
      lmt   = present(cmtout) ; if(lmt) cmtout = 0
      lpw   = present(cpwout) ; if(lpw) cpwout = 0
      lbcnt = all ( [( band(i+1)-band(i)==1 , i=1,size(band)-1 )] )
      if(lbcnt) then
        b1 = band(1)
        b2 = band(size(band))
      else if(str/=1) then
        Error('not implemented: str/=1 and band(:).')
      endif
      Hpos  = 0 ; if(l_soc) Hpos = 1
      call hdf_fopen(Hfile,'KS.hdf',0)
      do k = 1,size(kpt1)
        ikpt1 = kpt1(k)
        write(ckpt,'(I1,''-'',I6.6)') ispin,ikpt1
        if(lbcnt) then
          nb = min( nband(ikpt1,ispin) , b2 ) - b1 + 1
          if(str/=1) nb = (nb-1) / str + 1
          if(l_soc) then
            if(lmt) call hdf_rdwr(Hfile,'/'//ckpt//'/orbitals/mt',0,cmtout(:,:,         :nb,k,:),b1-1,str)
            if(lpw) call hdf_rdwr(Hfile,'/'//ckpt//'/orbitals/pw',0,cpwout(:ngpt(ikpt1),:nb,k,:),b1-1,str)
          else
            if(lmt) call hdf_rdwr(Hfile,'/'//ckpt//'/orbitals/mt',0,cmtout(:,:,         :nb,k,1),b1-1,str)
            if(lpw) call hdf_rdwr(Hfile,'/'//ckpt//'/orbitals/pw',0,cpwout(:ngpt(ikpt1),:nb,k,1),b1-1,str)
          endif
        else
          do ib = 1,size(band)
            if(band(ib)>nband(ikpt1,ispin)) cycle
            if(l_soc) then
              if(lmt) call hdf_rdwr(Hfile,'/'//ckpt//'/orbitals/mt',0,cmtout(:,:,         ib:ib,k,:),band(ib)-1,str)
              if(lpw) call hdf_rdwr(Hfile,'/'//ckpt//'/orbitals/pw',0,cpwout(:ngpt(ikpt1),ib:ib,k,:),band(ib)-1,str)
            else
              if(lmt) call hdf_rdwr(Hfile,'/'//ckpt//'/orbitals/mt',0,cmtout(:,:,         ib:ib,k,1),band(ib)-1,str)
              if(lpw) call hdf_rdwr(Hfile,'/'//ckpt//'/orbitals/pw',0,cpwout(:ngpt(ikpt1),ib:ib,k,1),band(ib)-1,str)
            endif
          enddo
        endif
      enddo
      call hdf_fclose(Hfile)
      call kpt_reorder2(3)
      end subroutine read_wavef

c ----------------------------------

# endif

c ----------------------------------

# include "readwrite.inc"      

c ----------------------------------

      end module readwrite
