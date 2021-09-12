c--------------------------------------------------------------------------------
c Copyright (c) 2021 Peter Grünberg Institut, Forschungszentrum Jülich, Germany
c This file is part of SPEX and available as free software under the conditions
c of the MIT license as expressed in the LICENSE file in more detail.
c--------------------------------------------------------------------------------

c Version for FLEUR (FleurMaX Release 5 version 33)
c
c -------------------------------------------------
c

# include "cppmacro.h"

      module readwrite

      use global
      use file
      use Hwrapper
      use hdf5
      use, intrinsic :: iso_fortran_env

      character(*), parameter :: dftlong    = 'FleurMaX Release 5 version 33'
      character(*), parameter :: dftshort   = 'fleurR5'
      character(*), parameter :: file_wavef = 'eig_gw.hdf'
c Parameters that must also be present in routine getinput
      logical                 :: invs
      integer, allocatable    :: gpt1(:,:,:)
      integer, allocatable    :: nlo(:),llo(:,:)
      integer                 :: nkpt1,nband0,kpt_order,nldau=0
      real_dp                 :: fac_kptadd,maxene
c Parameters that are needed for file input and output
      logical                 :: invs_inp
      integer                 :: neigd,nlod,lcutd
      integer, private        :: ntypd

      contains

c ----------------------------------

      subroutine read_param
      implicit none
      integer(HID_T) :: Hfile
      integer        :: i,l,itype
      call hdf_fopen (Hfile,'basis.hdf',0)
      call hdf_rdwr_a(Hfile,'meta/version',  0,i) ; if(i/=1) Error('Wrong version number.')
      call hdf_rdwr_a(Hfile,'general/invs',  0,i) ; invs_inp = i==1
      call hdf_rdwr_a(Hfile,'general/l_real',0,i) ; invs     = i==1
      call hdf_rdwr_a(Hfile,'general/jspins',0,nspin)
      call hdf_rdwr_a(Hfile,'general/l_soc', 0,i) ; l_soc    = i==1
      call hdf_rdwr_a(Hfile,'general/rkmax', 0,gcut)
      call hdf_rdwr_a(Hfile,'atoms/nAtoms',  0,ncent)
      call hdf_rdwr_a(Hfile,'atoms/nTypes',  0,ntype)
      allocate ( lcut(ntype),ztype(ntype),cent(3,ncent),nlo(ntype),neq(ntype) )
      call hdf_rdwr  (Hfile,'atoms/neq',     0,neq)
      call hdf_rdwr  (Hfile,'atoms/ztype',   0,cent(1,:ntype)) ; ztype = nint(cent(1,:ntype))
      call hdf_rdwr  (Hfile,'atoms/ratoms',  0,cent)
      call hdf_rdwr  (Hfile,'atoms/lmax',    0,lcut)
      call hdf_rdwr_a(Hfile,'atoms/lmaxd',   0,i) ; if(i/=maxval(lcut)) Bug('lmaxd/=maxval(lmax).')
      latpar = -1 ! flag for undefined
      call hdf_rdwr  (Hfile,'cell/amat',     0,lat)
      call hdf_rdwr  (Hfile,'atoms/nlo',     0,nlo)
      call hdf_rdwr_a(Hfile,'atoms/nlod',    0,nlod) ; if(nlod<maxval(nlo)) Bug('nlod<maxval(nlo)')
      allocate ( llo(nlod,ntype) )
      if(nlod>0) call hdf_rdwr(Hfile,'atoms/llo',0,llo)
      call hdf_fclose(Hfile)
      end subroutine read_param

c ----------------------------------

      subroutine read_symmetry
      implicit none
      integer(HID_T)       :: Hfile
      integer, allocatable :: mrot(:,:,:)
      real_dp, allocatable :: tau(:,:)
      integer              :: i
      call hdf_fopen (Hfile,'pot.hdf',0)
      call hdf_rdwr_a(Hfile,'structure-1/nop', 0,nsym) ; allocate ( sym(nsym),mrot(3,3,nsym),tau(3,nsym) )
      call hdf_rdwr  (Hfile,'structure-1/mrot',0,mrot)
      call hdf_rdwr  (Hfile,'structure-1/tau', 0,tau)
      do i = 1,nsym
        sym(i)%rot    = mrot(:,:,i)
        sym(i)%transl = tau(:,i)
        if(any(abs(2*sym(i)%transl-nint(2*sym(i)%transl))>1d-15))
     &    Error('INV and non-sym. vec. is not half a reciprocal lattice vector.')
      enddo
      deallocate(mrot,tau)
      call hdf_fclose(Hfile)
      end subroutine read_symmetry

c ----------------------------------

      subroutine read_radial
      use util, only: chr
      implicit none
      real_dp, allocatable :: bas(:,:,:),rarr(:,:),u(:),du(:),e(:)
      integer(HID_T)       :: Hfile
      integer              :: indx(0:maxlcut)
      integer              :: i,l,n,itype,ispin
      Hpos = 0
      call hdf_fopen(Hfile,'basis.hdf',0)
      call hdf_rdwr (Hfile,'atoms/jri',0,grid(:)%number)
      maxgrid = maxval(grid(:)%number)
      allocate ( bas1(maxgrid,maxindx,0:maxlcut,ntype,nspin2),
     &           bas2(maxgrid,maxindx,0:maxlcut,ntype,nspin2),
     &           ubas (maxindx,0:maxlcut,ntype,nspin2),
     &           dubas(maxindx,0:maxlcut,ntype,nspin2),
     &           ebas( maxindx,0:maxlcut,ntype,nspin2) )
      allocate ( rarr(maxgrid,ntype) )
      call hdf_rdwr(Hfile,'atoms/rmsh',0,rarr)
      do itype = 1,ntype
        grid(itype)%first     = rarr(1,itype)                    !
        grid(itype)%increment = log(rarr(2,itype)/rarr(1,itype)) ! Define grid
        grid(itype)%radius    = rarr(grid(itype)%number,itype)   !
      enddo
      deallocate(rarr)
      call hdf_rdwr_a(Hfile,'atoms/jmtd',0,n)
      if(n/=maxgrid) Error('Grid dimension (maxgrid,jmtd) does not agree.')
      do ispin = 1,nspin
        do itype = 1,ntype
          indx = 0
          allocate(bas(maxgrid,2,0:maxlcut))
          call hdf_rdwr(Hfile,'jsp_'//trim(chr(ispin))//'/itype_'//trim(chr(itype))//'/f',0,bas)
          bas1(:,1,:,itype,ispin) = bas(:,1,:)
          bas2(:,1,:,itype,ispin) = bas(:,2,:)
          call hdf_rdwr(Hfile,'jsp_'//trim(chr(ispin))//'/itype_'//trim(chr(itype))//'/g',0,bas)
          bas1(:,2,:,itype,ispin) = bas(:,1,:)
          bas2(:,2,:,itype,ispin) = bas(:,2,:)
          call hdf_rdwr(Hfile,'jsp_'//trim(chr(ispin))//'/itype_'//trim(chr(itype))//'/us',0,ubas(1,:,itype,ispin))
          call hdf_rdwr(Hfile,'jsp_'//trim(chr(ispin))//'/itype_'//trim(chr(itype))//'/dus',0,dubas(1,:,itype,ispin))
          call hdf_rdwr(Hfile,'jsp_'//trim(chr(ispin))//'/itype_'//trim(chr(itype))//'/uds',0,ubas(2,:,itype,ispin))
          call hdf_rdwr(Hfile,'jsp_'//trim(chr(ispin))//'/itype_'//trim(chr(itype))//'/duds',0,dubas(2,:,itype,ispin))
          call hdf_rdwr(Hfile,'jsp_'//trim(chr(ispin))//'/itype_'//trim(chr(itype))//'/el0',0,ebas(1,:,itype,ispin))
          ebas(2,:,itype,ispin) = ebas(1,:,itype,ispin)
          deallocate(bas)
          indx = 2
          if(nlo(itype)>0) then
            allocate(bas(maxgrid,2,nlo(itype)),u(nlo(itype)),du(nlo(itype)),e(nlo(itype)))
            call hdf_rdwr(Hfile,'jsp_'//trim(chr(ispin))//'/itype_'//trim(chr(itype))//'/flo',0,bas,0)
            call hdf_rdwr(Hfile,'jsp_'//trim(chr(ispin))//'/itype_'//trim(chr(itype))//'/ulos',0,u,0)
            call hdf_rdwr(Hfile,'jsp_'//trim(chr(ispin))//'/itype_'//trim(chr(itype))//'/dulos',0,du,0)
            call hdf_rdwr(Hfile,'jsp_'//trim(chr(ispin))//'/itype_'//trim(chr(itype))//'/ello0',0,e,0)
            do i = 1,nlo(itype)
              l                       = llo(i,itype)
              indx(l)                 = indx(l) + 1
              n                       = indx(l)
              bas1(:,n,l,itype,ispin) = bas(:,1,i)
              bas2(:,n,l,itype,ispin) = bas(:,2,i)
              ubas  (n,l,itype,ispin) = u(i)
              dubas (n,l,itype,ispin) = du(i)
              ebas  (n,l,itype,ispin) = e(i)
            enddo
            deallocate(bas,u,du,e)
          endif
        enddo
      enddo
      call hdf_fclose(Hfile)
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
      lcore_soc = .true.
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
              iarr(l) = iarr(l) + 1
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
              ecore(idum,l,itype,ispin) = ecore1
              nindxc(l,itype)           = nindxc(l,itype) + 1
              idum                      = nindxc(l,itype)
              read(iunit) ll,
     &          wcore2,ecore(idum,l,itype,ispin),
     &          core1(1:grid(itype)%number,idum,l,itype,ispin),
     &          core2(1:grid(itype)%number,idum,l,itype,ispin)
            endif
          enddo
        enddo
      enddo
      call fclose(iunit)
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
      deallocate (rarr,rarr2)
      end subroutine read_core

c ----------------------------------

      subroutine read_gpt
      implicit none
      integer(HID_T)       :: Hfile
      integer              :: n,ispin,ikpt,Herr
      integer, allocatable :: gpt2(:,:)
      logical              :: lexist
      character(38)        :: ckpt
      call kpt_reorder1(0)
      call hdf_fopen(Hfile,'basis.hdf',0)
      do ikpt = 1,nkpt1
        write(ckpt,'(F12.10,'','',F12.10,'','',F12.10)') kpt(:,ikpt)
        call h5lexists_f(Hfile,'jsp_1/kpt_'//ckpt,lexist,Herr)
        if(Herr/=0)     Error('Existence check of group failed.')
        if(.not.lexist) Error('Group jsp_1/kpt_'//ckpt//' does not exist in basis.hdf')
        call hdf_rdwr_a(Hfile,'jsp_1/kpt_'//ckpt//'/nv',0,ngpt(ikpt))
        if(nspin1==2) then
          call hdf_rdwr_a(Hfile,'jsp_2/kpt_'//ckpt//'/nv',0,n)
          if(n/=ngpt(ikpt)) Bug('Spin-dependent ngpt.')
        endif
      enddo
      n = maxval(ngpt(:nkpt1))
      allocate(gpt1(3,n,nkpt2),gpt2(3,n))
      gpt1 = 0
      do ikpt = 1,nkpt1
        write(ckpt,'(F12.10,'','',F12.10,'','',F12.10)') kpt(:,ikpt)
        call hdf_rdwr(Hfile,'jsp_1/kpt_'//ckpt//'/gvec',0,gpt1(:,:ngpt(ikpt),ikpt))
        if(nspin1==2) then
          call hdf_rdwr(Hfile,'jsp_2/kpt_'//ckpt//'/gvec',0,gpt2(:,:ngpt(ikpt)))
          if(any(gpt2(:,:ngpt(ikpt))/=gpt1(:,:ngpt(ikpt),ikpt))) Bug('Spin-dependent gpt.')
        endif
      enddo
      deallocate(gpt2)
      call hdf_fclose(Hfile)
      call kpt_reorder2(0)
      end subroutine read_gpt

c ----------------------------------

      subroutine read_ene
      use util, only: chr
      implicit none
      integer(HID_T)       :: Hfile
      integer              :: ispin,ikpt,n,Herr
      real_dp, allocatable :: enehlp(:)
      logical              :: lexist
      character(38)        :: ckpt
      call kpt_reorder1(0)
      call hdf_fopen(Hfile,'eig_gw.hdf',0)
      Hpos = 0
      ! Read number of eigenvalues in file
      do ispin = 1,nspin1
        do ikpt = 1,nkpt1
          write(ckpt,'(F12.10,'','',F12.10,'','',F12.10)') kpt(:,ikpt)
          call h5lexists_f(Hfile,'jsp_1/kpt_'//ckpt,lexist,Herr)
          if(Herr/=0)     Error('Existence check of group failed.')
          if(.not.lexist) Error('Group jsp_1/kpt_'//ckpt//' does not exist in eig_gw.hdf')
          call hdf_rdwr_a(Hfile,'jsp_'//trim(chr(ispin))//'/kpt_'//ckpt//'/numbands',0,nband(ikpt,ispin))
        enddo
      enddo
      maxband = maxval(nband(:nkpt1,:nspin1))
      allocate(enehlp(maxband))
      ! Read eigenvalues
      allocate(ene(maxband,nkpt1,nspin)) ; ene = 0
      do ispin = 1,nspin1
        do ikpt = 1,nkpt1
          write(ckpt,'(F12.10,'','',F12.10,'','',F12.10)') kpt(:,ikpt)
          call hdf_rdwr_a(Hfile,'jsp_'//trim(chr(ispin))//'/kpt_'//ckpt//'/numbands',0,n)
          call hdf_rdwr  (Hfile,'jsp_'//trim(chr(ispin))//'/kpt_'//ckpt//'/eig',0,enehlp(:n),0)
          n                  = min(n,maxband)
          ene(:n,ikpt,ispin) = enehlp(:n)
        enddo
      enddo
      call hdf_fclose(Hfile)
      call kpt_reorder2(0)
      end subroutine read_ene

c ----------------------------------

      subroutine write_kpt_xml(filename,listname,typename,kpts,nkpts,wghts,scale,label)
      use util, only: chr
      implicit none
      character(*), intent(in)           :: filename,listname,typename
      integer,      intent(in)           :: nkpts
      real_dp,      intent(in)           :: kpts(3,nkpts),wghts(nkpts),scale
      character(*), intent(in), optional :: label(:)
      integer                            :: iunit,iunit_s
      integer                            :: i,ios
      integer                            :: indent
      logical                            :: lwrite
      character(256)                     :: line
      iunit   = fopen('kpts.xml',status='old',action='readwrite')
      iunit_s = fopen(status='scratch',action='readwrite')
      lwrite  = .true.
      indent  = -1
      do
        read(iunit,'(A)',end=1,iostat=ios) line
        if(ios/=0) Error('Read error while reading from "kpts.xml".')
        if(line/=' '.and.indent==-1) indent = index(trim(line),' ',back=.true.)
        if(adjustl(line)=='</kPointLists>') exit
        i = index(line,'<kPointList ')
        if(lwrite) then
          if(i/=0.and.line(:i-1)==' '.and.index(line,'name="'//listname//'"')/=0) then
            lwrite = .false.
          else
            write(iunit_s,'(A)') trim(line)
          endif
        else
          if(adjustl(line)=='</kPointList>') lwrite = .true.
        endif
      enddo      
 1    line = ' '
      write(iunit_s,'(A,A)') line(:indent),'   <kPointList name="'//listname//'" count="'//Chr(nkpts)//'" type="'//typename//'">'
      do i = 1,nkpts
        write(iunit_s,'(A,A'NoA) line(:indent),'      <kPoint weight="'//trim(chr(dble(wghts(i)),'F10.5'))//'" '
        if(present(label)) then
          if(label(i)/=' ') write(iunit_s,'(A'NoA) 'label="'//trim(label(i))//'"' 
        endif
        write(iunit_s,'(A)') '>   '//
     &    Chr(kpts(1,i)*scale)//'/'//Chr(scale)//'   '//
     &    Chr(kpts(2,i)*scale)//'/'//Chr(scale)//'   '//
     &    Chr(kpts(3,i)*scale)//'/'//Chr(scale)//'   </kPoint>'
      enddo
      write(iunit_s,'(A,A)') line(:indent),'   </kPointList>'
      write(iunit_s,'(A,A)') line(:indent),'</kPointLists>'
      rewind(iunit_s)
      rewind(iunit)
      do
        read(iunit_s,'(A)',end=2,iostat=ios) line
        if(ios/=0) Error('Read error while reading from "kpts.xml" scratch file.')
        write(iunit,'(A)',iostat=ios) trim(line)
        if(ios/=0) Error('Write error while writing to "kpts.xml".')        
      enddo
 2    call fclose(iunit_s)
      call fclose(iunit)  
      end subroutine write_kpt_xml

      subroutine write_kpt
      implicit none
      real_dp, allocatable :: wghts(:)
      real_dp              :: scale
      logical              :: lwrite
      integer              :: i
      integer              :: kgv
      call kpt_reorder1(0)
      allocate(wghts(nkpt1))
      scale = kgv(nkpt3,3)
      scale = scale * fac_kptadd
      do i = 1,nkpt1
        if(i<=nkpti) then ; wghts(i) = count(kptp(:nkpt)==i)
        else              ; wghts(i) = 0
        endif
      enddo
      call write_kpt_xml('kpts.xml','spex','SPEX-mesh',kpt,nkpt1,wghts,scale)
      deallocate(wghts)
      write(6,'(/A)') 'File kpts.xml complemented with Spex k-point set.'
      write(6,'(A)')  '  To use it, set listName="spex" in inp.xml.'
      write(6,'(A/)') '  Also, do not forget to set gw="2" and numbands="(max number of bands)" '//
     &                'if Fleur should prepare input data for Spex.'
      call kpt_reorder2(0)
      end subroutine write_kpt

      subroutine write_qpt(qpt,nqpt,label)
      implicit none
      integer,      intent(in) :: nqpt
      real_dp,      intent(in) :: qpt(3,nqpt)
      character(*), intent(in) :: label(nqpt)
      real_dp                  :: scale
      integer                  :: i,kgv
      scale = kgv(nkpt3,3)
      call write_kpt_xml('kpts.xml','spex_path','path',qpt,nqpt,[(0d0,i=1,nqpt)],scale,label)
      write(6,'(A)')  'File kpts.xml complemented with Spex q-point path.'
      write(6,'(A/)') '  To use it, set listName="spex_path" in inp.xml.'
      end subroutine write_qpt
      
c ----------------------------------

      subroutine read_pot(mode)
      implicit none
      integer,     intent(in) :: mode
      integer(HID_T)          :: Hfile
      complex_dp, allocatable :: phasestar(:,:,:),vstar(:,:)
      real_dp                 :: rdum,gmax
      integer,    allocatable :: pointstar(:,:,:),ntypsy(:),nlh_tmp(:),llh_tmp(:,:),nmlh_tmp(:,:),mlh_tmp(:,:,:)
      complex_dp, allocatable :: clh_tmp(:,:,:)
      integer                 :: mx1,mx2,mx3,nstar,maxmlh,ispin,itype,j,ntypsd,iatom,nsym,nat,n
      logical                 :: lfirst
      if(all(mode/=[1,2])) Bug('Mode for read_pot unknown.')
      lfirst = .not.allocated(llh)
      call hdf_fopen(Hfile,'pot.hdf',0)
      ! vmt
      if(lfirst) then
        call hdf_rdwr_a(Hfile,'latharms-1/nlhd',0,maxlh) ; maxlh = maxlh + 1
        call hdf_rdwr_a(Hfile,'latharms-1/memd',0,maxmlh)
        call hdf_rdwr_a(Hfile,'latharms-1/ntypsd',0,ntypsd)
        call hdf_rdwr_a(Hfile,'structure-1/nat',0,nat)
        allocate (ntypsy(nat),nlh_tmp(ntypsd),llh_tmp(maxlh,ntypsd),nmlh_tmp(maxlh,ntypsd),
     &            mlh_tmp(maxmlh,maxlh,ntypsd),clh_tmp(maxmlh,maxlh,ntypsd))
        call hdf_rdwr(Hfile,'structure-1/ntypsy',0,ntypsy)
        call hdf_rdwr(Hfile,'latharms-1/nlh',0,nlh_tmp)
        call hdf_rdwr(Hfile,'latharms-1/llh',0,llh_tmp)
        call hdf_rdwr(Hfile,'latharms-1/nmem',0,nmlh_tmp)
        call hdf_rdwr(Hfile,'latharms-1/mlh',0,mlh_tmp)
        call hdf_rdwr(Hfile,'latharms-1/clnu',0,clh_tmp)
        allocate ( nlh(ntype),llh(maxlh,ntype),nmlh(maxlh,ntype),
     &             mlh(maxmlh,maxlh,ntype),clh(maxmlh,maxlh,ntype) )
        nlh  = 0 ; nmlh = 0 ; clh = 0
        llh  = 0 ; mlh  = 0
        do itype = 1,ntype
          iatom           = sum(neq(:itype-1))+1
          nsym            = ntypsy(iatom)
          nlh(itype)      = nlh_tmp(nsym)+1
          n               = nlh(itype)
          llh(:n,itype)   = llh_tmp(:n,nsym)
          nmlh(:n,itype)  = nmlh_tmp(:n,nsym)
          mlh(:,:n,itype) = mlh_tmp(:,:n,nsym)
          clh(:,:n,itype) = clh_tmp(:,:n,nsym)
        enddo
        deallocate(ntypsy,nlh_tmp,llh_tmp,nmlh_tmp,mlh_tmp,clh_tmp)
        if(maxval(nlh)/=maxlh)   Error('nlhd/=maxlh.')
        if(maxval(nmlh)/=maxmlh) Error('memd/=maxmlh.')
        allocate ( vmt   (maxgrid,maxlh,ntype,nspin) ,
     &             vmt_xc(maxgrid,maxlh,ntype,nspin) )
      endif
      vmt    = 0
      vmt_xc = 0
      call hdf_rdwr(Hfile,'pottot/in/fr',0,vmt,0)
      do ispin = 1,nspin
        do itype = 1,ntype
          vmt(:,1,itype,ispin) = vmt(:,1,itype,ispin) * sqrt(4*pi)
        enddo
      enddo
      if(mode==1) then
        call hdf_rdwr(Hfile,'potcoul/in/fr',0,vmt_xc,0)
        do ispin = 1,nspin
          do itype = 1,ntype
            vmt_xc(:,1,itype,ispin) = vmt_xc(:,1,itype,ispin) * rgrid(:,itype)
c            vmt_xc(:,:,itype,ispin) = vmt(:,:,itype,ispin) - vmt_xc(:,:,itype,ispin)
          enddo
        enddo
        vmt_xc = vmt - vmt_xc
      else if(mode==2) then
        call hdf_rdwr(Hfile,'potx/in/fr',0,vmt_xc,0)
      endif
      ! vpw
      call hdf_rdwr_a(Hfile,'stars-1/gmax',0,gmax)
      call hdf_rdwr_a(Hfile,'stars-1/mx1',0,mx1)
      call hdf_rdwr_a(Hfile,'stars-1/mx2',0,mx2)
      call hdf_rdwr_a(Hfile,'stars-1/mx3',0,mx3)
      allocate ( phasestar(-mx1:mx1,-mx2:mx2,-mx3:mx3) )
      allocate ( pointstar(-mx1:mx1,-mx2:mx2,-mx3:mx3) ) ; pointstar = -1
      if(lfirst) then
        allocate ( vpw   (-mx1:mx1,-mx2:mx2,-mx3:mx3,nspin) )
        allocate ( vpw_xc(-mx1:mx1,-mx2:mx2,-mx3:mx3,nspin) )
      endif
      call hdf_rdwr(Hfile,'stars-1/rgphs',0,phasestar)
      call hdf_rdwr(Hfile,'stars-1/ig',0,pointstar)
      if(any(pointstar<0)) Bug('Check stars-1/ig dimensions.')
      if(any(abs(imag(phasestar))>1d-12)) then
        Inv( Error('Star phases have nonzero imaginary part with inversion symmetry.') )
        Info('Star phases are complex-valued.')
      endif
      call hdf_rdwr_a(Hfile,'stars-1/ng3',  0,nstar) ; allocate(vstar(nstar,nspin))
      call hdf_rdwr  (Hfile,'pottot/in/fpw',0,vstar)
      call v_distrib(vpw,vstar)
      if(mode==1) then
        call hdf_rdwr(Hfile,'potcoul/in/fpw',0,vstar) ; if(nspin==2) vstar(:,2) = vstar(:,1)
        call v_distrib(vpw_xc,vstar) ; vpw_xc = vpw - vpw_xc
c        do ispin = 1,nspin
c          do j = -mx3,mx3
c            vpw_xc(:,:,j,ispin) = vpw(:,:,j,ispin) - vpw_xc(:,:,j,ispin) 
c          enddo
c        enddo
      else
        call hdf_rdwr(Hfile,'potx/in/fpw',   0,vstar)
        call v_distrib(vpw_xc,vstar)
      endif
      call hdf_fclose(Hfile)
      contains
      subroutine v_distrib(vpw,vstar)
      implicit none
      MCOMPLEX_dp, intent(out) :: vpw(-mx1:mx1,-mx2:mx2,-mx3:mx3,nspin)
      complex_dp,  intent(in)  :: vstar(:,:)
      complex_dp               :: cdum
      integer                  :: ispin,istar,i,j,k
      do ispin = 1,nspin
        do k = -mx3,mx3
          do j = -mx2,mx2
            do i = -mx1,mx1
              istar = pointstar(i,j,k)
              if(istar/=0) then
                if(sum(matmul(rlat,[i,j,k])**2)>gmax**2+1d-6)
     &            Bug('G point outside gmax sphere but pointer is nonzero.')
                cdum = phasestar(i,j,k) * vstar(istar,ispin)
                Inv( if(abs(imag(cdum))>1d-12) Error('Potential coefficient has imaginary part.') )
                vpw(i,j,k,ispin) = cdum
              else
                vpw(i,j,k,ispin) = 0d0
              endif
            enddo
          enddo
        enddo
      enddo
      end subroutine v_distrib
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

      subroutine read_vxc(vxcmat,band,nb,ikpt_in,ispin,l_qsgw)
      use, intrinsic :: iso_fortran_env
      implicit none
      integer,     intent(in)  :: nb,band(nb),ikpt_in,ispin
      logical,     intent(in)  :: l_qsgw
      MCOMPLEX_dp, intent(out) :: vxcmat(nb,nb)
      integer                  :: ikpt
      ! ...
      call kpt_reorder1(3)
      if     (ikpt_in<=nkpti)  then ; ikpt = ikpt_in
      else if(ikpt_in==nkpt+1) then ; ikpt = nkpti + 1
      else                          ; Error('k-point index out of range.')
      endif
      ! ...
      Error('Routine is not available: read_vxc.')
      ! ...
      call kpt_reorder2(3)
      if(l_qsgw) call read_qsgw(vxcmat,band,nb,ikpt_in,ispin)
      end subroutine read_vxc

c ----------------------------------

      subroutine read_qsgw(vxcmat,band,nb,ikpt_in,ispin)
      use, intrinsic :: iso_fortran_env
      implicit none
      integer,     intent(in)    :: nb,band(nb),ikpt_in,ispin
      MCOMPLEX_dp, intent(inout) :: vxcmat(nb,nb)
      integer                    :: ikpt
      ! ...
      call kpt_reorder1(3)
      if     (ikpt_in<=nkpti)  then ; ikpt = ikpt_in
      else if(ikpt_in==nkpt+1) then ; ikpt = nkpti + 1
      else                          ; Error('k-point index out of range.')
      endif
      ! ...
      Error('Routine is not available: read_qsgw.')
      ! ...
      call kpt_reorder2(3)
      end subroutine read_qsgw

c ----------------------------------

      subroutine write_param(lsoc)
      use util, only: chr
      use hdf5
      use Hwrapper
      implicit none
      logical, intent(in)  :: lsoc
      real_dp, allocatable :: bas(:,:,:),rarr(:,:),u(:),du(:),e(:)
      integer(HID_T)       :: Hfile,Hgrp,Hspin,Htype,Hkpt
      integer              :: Herr
      integer              :: indx(0:maxlcut)
      integer              :: i,l,n,itype,ispin,ikpt
      character(38)        :: ckpt
      Hpos = 0
      call hdf_fopen(Hfile,'basis.hdf',2)
      call h5gcreate_f(Hfile,'meta',   Hgrp,Herr) ; call h5gclose_f(Hgrp,Herr)
      call h5gcreate_f(Hfile,'general',Hgrp,Herr) ; call h5gclose_f(Hgrp,Herr)
      call h5gcreate_f(Hfile,'atoms',  Hgrp,Herr) ; call h5gclose_f(Hgrp,Herr)
      call h5gcreate_f(Hfile,'cell',   Hgrp,Herr) ; call h5gclose_f(Hgrp,Herr)
      i = 1
      call hdf_rdwr_a(Hfile,'meta/version',  1,i)
      if(invs_inp) then ; i = 1
      else              ; i = 0
      endif
      call hdf_rdwr_a(Hfile,'general/invs',  1,i)
      if(invs) then ; i = 1
      else          ; i = 0
      endif
      call hdf_rdwr_a(Hfile,'general/l_real',1,i)
      call hdf_rdwr_a(Hfile,'general/jspins',1,nspin)
      if(l_soc) then ; i = 1
      else           ; i = 0
      endif
      call hdf_rdwr_a(Hfile,'general/l_soc', 1,i)
      call hdf_rdwr_a(Hfile,'general/rkmax', 1,gcut)
      call hdf_rdwr_a(Hfile,'atoms/nAtoms',  1,ncent)
      call hdf_rdwr_a(Hfile,'atoms/nTypes',  1,ntype)
      call hdf_rdwr  (Hfile,'atoms/neq',     1,neq)
      call hdf_rdwr  (Hfile,'atoms/ztype',   1,ztype)
      call hdf_rdwr  (Hfile,'atoms/ratoms',  1,cent)
      call hdf_rdwr  (Hfile,'atoms/lmax',    1,lcut)
      call hdf_rdwr_a(Hfile,'atoms/lmaxd',   1,maxval(lcut))
      call hdf_rdwr_a(Hfile,'cell/scaleCell',1,latpar)
      call hdf_rdwr  (Hfile,'cell/amat',     1,lat)
      call hdf_rdwr  (Hfile,'atoms/nlo',     1,nlo)
      call hdf_rdwr_a(Hfile,'atoms/nlod',    1,nlod)
      call hdf_rdwr  (Hfile,'atoms/llo',     1,llo)
      call hdf_rdwr  (Hfile,'atoms/jri',     1,grid(:)%number)
      maxgrid = maxval(grid(:)%number)
      call hdf_rdwr_a(Hfile,'atoms/jmtd',    1,maxgrid)
      allocate ( rarr(maxgrid,ntype) )
      rarr = 0
      do itype = 1,ntype
        rarr(1,itype)                  = grid(itype)%first
        rarr(2,itype)                  = rarr(1,itype)*exp(grid(itype)%increment)
        rarr(grid(itype)%number,itype) = grid(itype)%radius
      enddo
      call hdf_rdwr  (Hfile,'atoms/rmsh',    1,rarr)
      deallocate(rarr)
      do ispin = 1,nspin
        call h5gcreate_f(Hfile,'jsp_'//trim(chr(ispin)),Hspin,Herr)
        do itype = 1,ntype
          call h5gcreate_f(Hspin,'itype_'//trim(chr(itype)),Htype,Herr)
          allocate(bas(maxgrid,2,0:maxlcut))
          bas(:,1,:) = bas1(:,1,:,itype,ispin)
          bas(:,2,:) = bas2(:,1,:,itype,ispin)
          call hdf_rdwr(Htype,'f',1,bas)
          bas(:,1,:) = bas1(:,2,:,itype,ispin)
          bas(:,2,:) = bas2(:,2,:,itype,ispin)
          call hdf_rdwr(Htype,'g',1,bas)
          deallocate(bas)
          call hdf_rdwr(Htype,'us',1,ubas(1,:,itype,ispin))
          call hdf_rdwr(Htype,'dus',1,dubas(1,:,itype,ispin))
          call hdf_rdwr(Htype,'uds',1,ubas(2,:,itype,ispin))
          call hdf_rdwr(Htype,'duds',1,dubas(2,:,itype,ispin))
          call hdf_rdwr(Htype,'el0',1,ebas(1,:,itype,ispin))
          indx = 2
          if(nlo(itype)>0) then
            allocate(bas(maxgrid,2,nlo(itype)),u(nlo(itype)),du(nlo(itype)),e(nlo(itype)))
            do i = 1,nlo(itype)
              l          = llo(i,itype)
              indx(l)    = indx(l) + 1
              n          = indx(l)
              bas(:,1,i) = bas1(:,n,l,itype,ispin)
              bas(:,2,i) = bas2(:,n,l,itype,ispin)
              u(i)       = ubas(n,l,itype,ispin)
              du(i)      = dubas(n,l,itype,ispin)
              e(i)       = ebas(n,l,itype,ispin)
            enddo
            call hdf_rdwr(Htype,'flo',1,bas)
            call hdf_rdwr(Htype,'ulos',1,u)
            call hdf_rdwr(Htype,'dulos',1,du)
            call hdf_rdwr(Htype,'ello0',1,e)
            deallocate(bas,u,du,e)
          endif
          call h5gclose_f(Htype,Herr)
        enddo
        do ikpt = 1,nkpt1
          write(ckpt,'(F12.10,'','',F12.10,'','',F12.10)') kpt(:,ikpt)
          call h5gcreate_f(Hspin,'kpt_'//ckpt,Hkpt,Herr)
          call hdf_rdwr_a(Hkpt,'nv',1,ngpt(ikpt))
          call hdf_rdwr(Hkpt,'gvec',1,gpt1(:,:ngpt(ikpt),ikpt))
          call h5gclose_f(Hkpt,Herr)
        enddo
        call h5gclose_f(Hspin,Herr)
      enddo
      call hdf_fclose(Hfile)
      end subroutine write_param

c ----------------------------------

      subroutine write_wavef(ikpt,ispin,ngpt1,gpt1,eig,cmtin,cpwin,cpwin_c)
      use, intrinsic :: iso_fortran_env
      use util, only: chr
      use global
      use Hwrapper
      use hdf5
      Mpi( use Mwrapper )
      implicit none
      integer,     intent(in)           :: ikpt,ispin,ngpt1
      integer,     intent(in)           :: gpt1(3,ngpt1)
      real_dp,     intent(in)           :: eig(:)
      complex_dp,  intent(in)           :: cmtin(maxlmindx,ncent,size(eig),nspin3)
      MCOMPLEX_dp, intent(in), optional :: cpwin(maxgpt,size(eig),nspin3)
      complex_dp,  intent(in), optional :: cpwin_c(maxgpt,size(eig),nspin3)
      integer        :: Herr
      integer(HID_T) :: Hfile, Hkpt, Hspin
      character(80)  :: filename, spath, kpath
      character(38)  :: ckpt
      logical        :: l_exist, link_exists

      if(present(cpwin).eqv.present(cpwin_c)) Bug('cpwin and cpwin_c both present or missing.')
      if(ngpt1>maxgpt)                        Bug('ngpt1 > maxgpt.')

      write(ckpt,'(F12.10,'','',F12.10,'','',F12.10)') kpt(:,ikpt)
      spath = 'jsp_'//trim(chr(ispin))
      kpath = trim(spath)//'/kpt_'//ckpt

      filename = 'eig_gw.hdf'
      if((ispin==1).and.(ikpt==1)) then
        call hdf_fopen(Hfile,trim(filename),2) ; call hdf_rdwr_a(Hfile,'version',1,1)
      else
        call hdf_fopen(Hfile,trim(filename),1)
      endif
      ! Group of spin
      if(ikpt==1) call h5gcreate_f(Hfile,trim(spath),Hspin,Herr)
      call h5gcreate_f(Hfile,trim(kpath),Hkpt,Herr)
      call hdf_rdwr_a(Hkpt,'numbands',1,size(eig))
      call hdf_rdwr(Hkpt,'eig',1,eig)
      if(present(cpwin)) then ; call hdf_rdwr(Hkpt,'pw',1,cpwin)
      else                    ; call hdf_rdwr(Hkpt,'pw',1,cpwin_c)
      endif
      call hdf_rdwr(Hkpt,'mt',1,cmtin)
      call h5gclose_f(Hkpt,Herr)
      if(ikpt==1) call h5gclose_f(Hspin,Herr)
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
      integer        :: ispin,k,ng,nb MpiC(Merr)
      integer(HID_T) :: Hfile
      integer        :: k1
      character(38)  :: ckpt
      complex_dp,  allocatable :: cmthlp(:,:,:,:)!maxlmindx,ncent,maxband,nspin3)
      MCOMPLEX_dp, allocatable :: cpwhlp(:,:,:)!maxgpt,maxband,nspin3)
      call kpt_reorder1(1)
      call hdf_fopen(Hfile,'eig_gw.hdf',0)
      Hpos = 0 ; if(l_soc) Hpos = 1
      ! Read states
      do ispin = 1,nspin1
        do k = 1,nkpti+nkpti2
          write(ckpt,'(F12.10,'','',F12.10,'','',F12.10)') kpt(:,k)
          ng = ngpt(k)
          nb = nband(k,ispin)
          ! read file
          Allocate_(cmthlp,(maxlmindx,ncent,nb,nspin3))
          call hdf_rdwr(Hfile,'jsp_'//trim(chr(ispin))//'/kpt_'//ckpt//'/mt',0, cmthlp,0)
          if(l_soc) then
            ifO cmt(:,:,:nb,k,:)     = cmthlp
          else
            ifO cmt(:,:,:nb,k,ispin) = cmthlp(:,:,:,1)
          endif
          Deallocate_(cmthlp)
          Allocate_(cpwhlp,(ng,nb,nspin3))
          call hdf_rdwr(Hfile,'jsp_'//trim(chr(ispin))//'/kpt_'//ckpt//'/pw',0, cpwhlp,0)
          if(l_soc) then
            ifO cpw(:ng,:nb,k,:)     = cpwhlp
          else
            ifO cpw(:ng,:nb,k,ispin) = cpwhlp(:,:,1)
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
      integer                            :: b1,b2
      b1    = 0
      b2    = 0
      lmt   = present(cmtout) ; if(lmt) then ; cmtout = 0 ; b1 =     size(cmtout,3)     ; b2 =     size(cmtout,4)     ;  endif
      lpw   = present(cpwout) ; if(lpw) then ; cpwout = 0 ; b1 = min(size(cpwout,2),b1) ; b2 = min(size(cpwout,3),b2) ; endif
      lbcnt = all ( [( band(i+1)-band(i)==1 , i=1,size(band)-1 )] )
      if(str/=1 Mpi(.and.str/=Msize) ) Bug('Wrong stride (str) argument.')
      if(b1<size(band))                Bug('Band dimension too small.')
      if(b2<size(kpt1))                Bug('kpoint dimension too small.')
      call kpt_reorder1(3)
      if(lbcnt) then
        b1 = band(1)
        b2 = band(size(band))
      else if(str/=1) then
        Error('not implemented: str/=1 and band(:).')
      endif
      ! Hpos = ?
      ! ...
      Error('Routine is not available: read_wavef (LOAD version).')
      ! ...
      call kpt_reorder2(3)
      end subroutine read_wavef

c ----------------------------------

# endif

c ----------------------------------

# include "readwrite.inc"

c ----------------------------------

      end module readwrite
