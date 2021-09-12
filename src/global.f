c--------------------------------------------------------------------------------
c Copyright (c) 2021 Peter Grünberg Institut, Forschungszentrum Jülich, Germany
c This file is part of SPEX and available as free software under the conditions
c of the MIT license as expressed in the LICENSE file in more detail.
c--------------------------------------------------------------------------------

c Global module for parameters and arrays.
c
c Special datatypes (see "cppmacro.h")
c   real_dp     = double precision (8 bytes)
c   complex_dp  = double complex (16 bytes)      
c   integer_dp  = large-range integer (8 bytes)
c   MCOMPLEX_dp = real_dp if INV is defined, otherwise complex_dp

# include "cppmacro.h"

      module global

      use, intrinsic :: iso_fortran_env
      use, intrinsic :: iso_c_binding

      implicit none
      public

      Mpi( include 'mpif.h' )

      integer                    :: inp,uin
      character(128)             :: inpfile

c Input parameters
      logical                    :: metal,lkptadd,storeibz,alignbd,obloechl,ozero,trsoff,
     &                              tetraf,coulstore,fullpw,ologdiv,noapw,l_soc,l_qsgw,
     &                              lcore_soc,use_sym,bandinfo,l_multipole
      integer                    :: restart,ovxc,option,otimer
      integer                    :: njob,oselfc,freqint,ogrid
      integer,       allocatable :: smooth(:),omit(:)      
      real_dp                    :: tolerance
      real_dp                    :: maxmem,ene_extra,gcutf,cutzero,deltaex,it_elow
      logical                    :: lomit(3),wrtsuscep,wrtdielec,wrtscreen,wrtinfo,wrtext MpiC(mpiprod(2))
      logical,       allocatable :: cores(:,:,:)
      integer                    :: nkpt_label
      real_dp,       allocatable :: kpt_label(:,:)
      character(10), allocatable :: ckpt_label(:)
      type jobtype
        integer                  :: indx,type,out
        integer,     allocatable :: band(:),kpt(:),spin(:)
        real_dp                  :: freq(3)
        character                :: label*6,freqfile*256,kern*4
        logical                  :: full
      end type jobtype
      type(jobtype), allocatable :: job(:)

c System
      integer                  :: ntype,ncent,nelec
      integer                  :: nspin,nspin1,nspin2,nspin3 ! nspin1/2/3 only for SOC
      integer,     allocatable :: neq(:),ztype(:),atype(:)
      real_dp                  :: latpar,lat(3,3),rlat(3,3),vol,rvol,sqa(2)
      real_dp,     allocatable :: cent(:,:)
      integer,     allocatable :: pcent(:,:),tcent(:,:,:)

c Wave functions
      real_dp                  :: checksum
      ! core
      integer,     allocatable :: lcutc(:),nindxc(:,:)
      integer                  :: maxlcutc,maxindxc
      real_dp,     allocatable :: ecore(:,:,:,:),core1(:,:,:,:,:),core2(:,:,:,:,:)
      ! valence
      integer                  :: maxband,maxeband,bando,bandu
      integer,     allocatable :: nband(:,:)
      MCOMPLEX_dp, pointer_cnt :: cpw(:,:,:,:)
      complex_dp,  pointer_cnt :: cmt(:,:,:,:,:)
      complex_dp,  allocatable :: phase(:)
      real_dp,     allocatable :: ene(:,:,:)
      integer,     allocatable :: deg(:,:,:)
      real_dp                  :: efermi,egap,eip
      integer                  :: maxlmindx
# ifdef LOAD
      MCOMPLEX_dp, allocatable, target :: cpwq(:,:,:)
      complex_dp,  allocatable, target :: cmtq(:,:,:,:)
# elif defined(MPI)
      integer                  :: win_cpw,win_cmt
# endif
      ! basis
      integer                  :: maxlcut,maxindx,maxgpt,ngptall
      integer,     allocatable :: lcut(:),nindx(:,:),lcutp(:)
      real_dp                  :: gcut
      real_dp,     allocatable :: bas1(:,:,:,:,:),bas2(:,:,:,:,:),ubas(:,:,:,:),dubas(:,:,:,:),ebas(:,:,:,:),minebas(:)
      integer,     allocatable :: gpt(:,:),pgpt(:,:),ngpt(:),pntgpt(:,:,:)

      ! subspaces / irreps
      complex_dp,  allocatable :: irrep_sub(:,:,:,:)
      integer,     allocatable :: psub(:),block(:,:),sizeblock(:)
      integer                  :: nsub,nblock

c Mixed basis
      integer                  :: maxbasm
      integer,     allocatable :: nbasm(:)
      integer                  :: optm
      real_dp                  :: optcutm
      ! MT
      integer                  :: maxlcutm,maxindxm,maxlmindxm,nbasp,maxbasp,maxindxp
      integer,     allocatable :: lcutm(:)
      integer,     allocatable :: nindxm(:,:)
      real_dp,     allocatable :: basm(:,:,:,:)
      real_dp,     allocatable :: dbasm(:,:,:)
      ! IPW
      real_dp                  :: gcutm
      integer                  :: maxgptm,ngptmall
      integer,     allocatable :: ngptm(:),pntgptm(:,:,:,:)
      integer,     allocatable :: gptm(:,:),pgptm(:,:)

c K-point set
      integer                  :: nkpt,nkpt2,nkpti,nkpti2
      real_dp,     allocatable :: kpt(:,:)
      integer,     allocatable :: kptp(:),kindx(:)
      integer,     allocatable :: nkpt3(:)
      integer,     allocatable :: symkpt(:),kptsym(:,:),gkptsym(:,:,:)
      real_dp                  :: kptadd(3)
      real_dp,     allocatable :: wintgr(:,:,:),gauss(:)
      integer                  :: nkpt_path,mkpt_path
      real_dp,     allocatable :: kpt_path(:,:)
      integer, allocatable, target :: pkpt(:,:,:),pkpt1(:,:,:)

c Symmetry
      integer                  :: invsym,nsym,nsymt,ngensym
      integer,     allocatable :: symtab(:,:),gensym(:),symgen(:,:)
      type symtype
        integer                :: rot(3,3),rrot(3,3)
        real_dp                :: transl(3)
        complex_dp             :: esoc(2,2)
        integer                :: inv
        logical                :: trs
      end type symtype
      type(symtype), allocatable :: sym(:)

c Potentials
# ifndef TDDFT
      integer                  :: maxlh      
      integer,     allocatable :: nlh(:),llh(:,:),nmlh(:,:),mlh(:,:,:),symcent(:)
      complex_dp,  allocatable :: clh(:,:,:)
      real_dp                  :: vavg

      real_dp,     pointer_cnt :: vmt (:,:,:,:),vmt_xc (:,:,:,:)
      MCOMPLEX_dp, pointer_cnt :: vpw (:,:,:,:),vpw_xc (:,:,:,:)
      MCOMPLEX_dp, pointer_cnt :: vpw_(:,:,:,:),vpw_xc_(:,:,:,:)
#   ifdef MPI
      integer                  :: win_vmt,win_vmt_xc,win_vpw,win_vpw_xc,win_vpw_,win_vpw_xc_
#   endif
# endif

c Coulomb matrix
      MCOMPLEX_dp, pointer_cnt :: coulomb(:,:)      
      MCOMPLEX_dp, pointer_cnt :: coulomb0(:)
      MCOMPLEX_dp, allocatable :: constfunc(:)
      integer,     allocatable :: cblock(:)
      integer                  :: lexp MpiC(win_coulomb0) MpiC(win_coulomb)
      complex_dp,  allocatable :: screenk_mt(:,:,:,:)
      MCOMPLEX_dp, allocatable :: screenk_pw(:,:,:)

c IBC
      integer                  :: oibc
      real_dp,     allocatable :: ibc_proj(:,:,:,:,:,:)
      complex_dp,  allocatable :: ibc_suscep0(:,:,:,:)
      complex_dp,  allocatable :: ibc_selfc(:,:)

c Memory bookkeeping
      real_dp,     pointer     :: mem  ! mem = total memory used in virtual node
      real_dp,     target      :: mem_ !
      integer                  :: nalloc_ MpiC(win_mem)
# if defined(MPI) && defined(WRTMEM)
      real_dp,     pointer     :: mem0 ,mem_prv0  ! mem0 (mem_prv0) = total (private) memory used in physical node (only MPI & WRTMEM)
      real_dp,     target      :: mem0_,mem_prv0_ !
      integer                  :: win_mem0,win_mem_prv0
# endif

c Miscellaneous
      real_dp                  :: radshmin,divergence
      real_dp,     allocatable :: rpa_int(:)
      real_dp,     allocatable :: fac(:),sfac(:),facfac(:)
      complex_dp,  allocatable :: dwgn(:,:,:,:)
      MCOMPLEX_dp, allocatable :: olap_prod(:,:,:),cstep(:,:,:)
      integer,     allocatable :: olap_dim(:),pnt_prod(:,:,:,:)
      integer                  :: maxolap_dim
      integer                  :: symm_maxindx = 8 ! must be increased if necessary (see trafo.f)
      character, parameter     :: lchar(0:38) =
     &                            (/'s','p','d','f','g','h','i','j','k','l','m','n','o',
     &                              'x','x','x','x','x','x','x','x','x','x','x','x','x',
     &                              'x','x','x','x','x','x','x','x','x','x','x','x','x' /)
      type statetype
        integer                :: b,k,s
      end type statetype

c Radial grids
      type gridtype
        integer                :: number
        real_dp                :: first,increment,radius
      end type gridtype
      type(gridtype), allocatable :: grid(:)
      integer                     :: maxgrid
      real_dp, allocatable        :: gridf(:,:),rgrid(:,:)
      real_dp, allocatable        :: rgrid0(:,:) ! only for interpolations

c Accuracy of degenerate levels (used as a criterion for defining degenerate subspaces)
      real_dp,     parameter   :: edeg0 = 1d-8, cdeg = 1d-8
      real_dp                  :: edeg  = edeg0 ! will be refined in getinput

c Other accuracies
      real_dp,     parameter   :: acc_checksum = 1d-8

c Constants
      real_dp,     parameter   :: pi = 3.1415926535897932384626433d0
      real_dp,     parameter   :: hartree = 27.21138386d0! (old value: 27.211608d0)
      real_dp,     parameter   :: escale = hartree ! unit for data files (in eV if escale=hartree)
      real_dp,     parameter   :: clight = 137.0359895d0
      real_dp,     parameter   :: angstrom = 0.5291772108d0      
      integer,     parameter   :: kilobyte = 1024 , megabyte = 1048576
      complex_dp,  parameter   :: img = (0d0,1d0)

c Wannier
      integer                  :: nwan,nwanband,wanbandi,wanbandf
      complex_dp, pointer_cnt  :: cmtu(:,:,:,:,:),cpwu(:,:,:,:)
      complex_dp, allocatable  :: uwan(:,:,:,:),irrep_wan(:,:,:,:),irrep_wan1(:,:,:,:)      
      real_dp,    allocatable  :: wancent(:,:,:)
      real_dp                  :: wscale
      integer                  :: spinw
      logical                  :: wbloch
c      logical                  :: l_3rdOrder,l_numint
c      logical                  :: l_GoldstoneModeEnforced
c      character(80)            :: RestartFileLocation
c      integer                  :: FirstOrderInBSE
c      real_dp                  :: eta,deltaEX,wzero

c MPI Parallelisation
# ifdef MPI
      integer                   :: Msplit            ! type for mpi_comm_split_type
      integer                   :: Mrank0,Msize0     ! for mpi_comm_world
      integer                   :: Mcomm,Mrank,Msize ! current communicator (mostly mpi_comm_world)
      integer                   :: Ncomm,Nrank,Nsize ! in-node communicator
      integer                   :: Ocomm,Orank,Osize ! across-node communicator (containing only Nrank=0)
      integer, allocatable      :: Opnt(:)           ! Mrank -> Orank pointer
      integer, allocatable      :: Mcomms(:)         ! hirarchical layers of communicators (for begin/end/change_split)
      integer                   :: Nrank0,Nsize0,Orank0 ! communicator for physical node
#   ifndef noPID
      integer, allocatable      :: Mpid(:)           ! PID for each process in the node
#   endif
#   ifdef noMSYNC
      integer(mpi_address_kind) :: Maddress_dummy
#   endif
# endif

# ifdef TDDFT
      !xMN  new global variable for the exchange-correlation kernel
      !     array. It is stored in packed form (only one index).
      !     Later one might like to have two variables for magnetic and
      !     dielectric kernel.
      MCOMPLEX_dp,allocatable :: xckernel(:)
      real_dp,    allocatable :: fxc_rad(:), fxc_rad_sqrt(:)
# endif

      contains

c     ------------------

c     Replaces the function modulo(kpoint,1d0) for a kpoint in kpt(:,ikpt).
      function modulo1(kpoint)
      implicit none
      real_dp, intent(in) :: kpoint(3)
      real_dp             :: modulo1(3)
      integer             :: help(3)
      modulo1 = kpoint*nkpt3
      help    = nint(modulo1)
      if(any(abs(help-modulo1)>1d-10)) then
        write(0,'(A,F10.8,2('','',F10.8),A)') 'Argument (',kpoint,') is not an element of the k-point set.'
        write(0,'(A)') 'This usually happens when the kpoint set read from the input data is different from the one defined in BZ.'
        Error('argument not an element of k-point set.')
      endif
      modulo1 = modulo(help,nkpt3)*1d0/nkpt3
      end function modulo1

c     Same for shifted k-point set
      function modulo1r(kpoint)
      implicit none
      real_dp, intent(in) :: kpoint(3)
      real_dp             :: modulo1r(3)
      integer             :: i
      modulo1r = modulo ( kpoint+1d-13 , 1d0 ) - 1d-13
      do i = 1,3
        if(modulo1r(i)<0.or.abs(1-modulo1r(i))<1d-13) modulo1r(i) = 0d0 ! avoid signed zero (-0d0) and values close to 1
      enddo
      end function modulo1r

c     ------------------

# if 0
c     ifc seems to have problems transposing integer arrays. this is a fix.
      function transpose_int ( a )
      implicit none
      integer transpose_int(3,3),a(3,3)
      integer i,j
      do i = 1,3
        do j = 1,3
          transpose_int(i,j) = a(j,i)
        enddo
      enddo
      end function transpose_int
# endif

c     ------------------

c     copy_job:  a replacement for the statement:  job(j) = job(j-1)
      subroutine copy_job(j, new_job)
      implicit none
      integer,       intent(in)   ::  j       ! j = job_index
      type(jobtype), intent(out)  ::  new_job
      if(allocated(new_job%band)) deallocate ( new_job%band )
      if(allocated(new_job%kpt))  deallocate ( new_job%kpt  )
      if(allocated(new_job%spin)) deallocate ( new_job%spin )
      if(allocated(job(j)%band)) then
        allocate ( new_job%band(size(job(j)%band)) )
                   new_job%band(:) = job(j)%band(:)
      end if
      allocate ( new_job%kpt(size(job(j)%kpt )) )
                 new_job%kpt(:) = job(j)%kpt(:)
      allocate ( new_job%spin(size(job(j)%spin)) )
                 new_job%spin(:) = job(j)%spin(:)
      new_job%type      = job(j)%type
      new_job%out       = job(j)%out
      new_job%freq      = job(j)%freq
      new_job%label     = job(j)%label
      new_job%freqfile  = job(j)%freqfile
      new_job%full      = job(j)%full
      new_job%kern      = job(j)%kern
      end subroutine copy_job

c     ------------------

      end module global
