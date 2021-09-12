c--------------------------------------------------------------------------------
c Copyright (c) 2021 Peter Grünberg Institut, Forschungszentrum Jülich, Germany
c This file is part of SPEX and available as free software under the conditions
c of the MIT license as expressed in the LICENSE file in more detail.
c--------------------------------------------------------------------------------

c Construction of Wannier functions

c Treatment of groups of degenerate states:
c   0 : none - cut through degenerate states
c   1 : include them completely
c   2 : exclude them completely
# define degstates 1

c The Wannier centers are set at
c   0 : location returned by Wannier90
c   1 : first-guess center
c   2 : first-guess center if difference is smaller than 1 Bohr, location returned by Wannier90 otherwise
# define set_wancent 2

c The maximally localized Wannier functions can be reprojected onto the first-guess functions (for testing)
c# define reproject

# if degstates != 1
#   warning degstates != 1
# endif

# ifdef reproject
#   warning reproject defined
# endif

# include "cppmacro.h"
# include "restype.h"      

# ifdef LOAD
#   define cmt_ cmt
#   define cpw_ cpw
# endif

      subroutine wannier

      use global
      use wrapper
      use readwrite
      use file
      use util
      use key
      use, intrinsic :: iso_fortran_env
# ifdef MPI      
      use Mwrapper, only: Mcast,Mcastl
# endif      
      Load( use readwrite )

      implicit none
      complex_dp,   allocatable :: proj(:,:),cmt1(:,:,:)
      complex_dp,   allocatable :: mat1(:,:),mat2(:,:),decomp(:,:,:)
      real_dp,      allocatable :: eig(:),fwancent(:,:,:),enehlp(:,:,:)
      real_dp                   :: kvec(3),kvec0(3),dev,lat1(3,3)
      real_dp                   :: chksum,rdum,rdum1,dist0,dist1,dist(ncent)
      integer                   :: wanbandi0,wanbandf0
      integer                   :: backf(3,ncent),backf1(3),pnt(ncent),ilat(3,3),ivec(3)
      integer                   :: ikpt,ikpt1,ispin,s,sb,off,koff,soff
      integer                   :: iwan,iwan1
      integer                   :: itype,ieq,ic,ios,vers
      integer                   :: i,j,i1,i2,i3,l,m,lm,n,l1
      real_dp                   :: intgrf
      integer                   :: iunit,ibpt
      integer                   :: nntot,nsmall1,nsmall2
      integer,      allocatable :: nnlist(:,:),nncell(:,:,:),neigh(:,:)
      integer,      allocatable :: proj_l(:),proj_m(:),proj_radial(:),exclude_bands(:)
      real_dp,      allocatable :: cent1(:,:),proj_site(:,:),proj_z(:,:),proj_x(:,:),proj_zona(:)
      complex_dp,   allocatable :: kolap(:,:,:,:),umat(:,:,:),uwan1(:,:,:),fwan(:,:)
      real_dp,      allocatable :: wann_spreads(:)
      logical,      allocatable :: lwindow(:,:)
      logical                   :: ldum,second,lauto
      logical                   :: l_wanmax,l_wanread,l_inter,l_frozen
      real_dp,      parameter   :: mwarn1 = 1d-6, mwarn2 = 0.01d0
      integer,      parameter   :: nlat_backf = 3
      integer,      allocatable :: lmwan(:),centwan(:),lwan(:)
      real_dp,      allocatable :: eulerwan(:,:),efrozen(:)
      real_dp                   :: kpt1(3,nkpt)
      complex_dp                :: cdum
      complex_dp                :: olapmt(maxlmindx,maxlmindx,ncent,nspin)
      complex_dp                :: cmt2(maxlmindx,ncent,2)
# ifndef LOAD
      complex_dp,   allocatable :: cmt_(:,:,:,:,:)
      MCOMPLEX_dp,  allocatable :: cpw_(:,:,:,:)
# endif
      MCOMPLEX_dp,  allocatable :: cpw2(:),olappw(:,:)
      real                      :: cputime
      integer                   :: ch2i
      character(:), allocatable :: charr(:)
      character(10)             :: parse_out(3)
      character(2), allocatable :: zlabel1(:)
      character(2), parameter   :: zlabel(112) = (/
     &  'H ',                                                                                'He',
     &  'Li','Be',                                                  'B ','C ','N ','O ','F ','Ne',
     &  'Na','Mg',                                                  'Al','Si','P ','S ','Cl','Ar',
     &  'K ','Ca','Sc','Ti','V ','Cr','Mn','Fe','Co','Ni','Cu','Zn','Ga','Ge','As','Se','Br','Kr',
     &  'Rb','Sr','Y ','Zr','Nb','Mo','Tc','Ru','Rh','Pd','Ag','Cd','In','Sn','Sb','Te','I ','Xe',
     &  'Cs','Ba','La','Ce','Pr','Nd','Pm','Sm','Eu','Gd','Tb','Dy','Ho','Er','Tm','Yb',
     &            'Lu','Hf','Ta','W ','Re','Os','Ir','Pt','Au','Hg','Tl','Pb','Bi','Po','At','Rn',
     &  'Fr','Ra','Ac','Th','Pa','U ','Np','Pu','Am','Cm','Bk','Cf','Es','Fm','Md','No',
     &            'Lr','Rf','Db','Sg','Bh','Hs','Mt','Ds','Rg','Cn' /)

# include "interface/orbitals.inc"
# ifdef WAN
      interface
        subroutine wannier_setup(seed__name, mp_grid_loc, num_kpts_loc,
     +    real_lattice_loc, recip_lattice_loc, kpt_latt_loc,
     +    num_bands_tot, num_atoms_loc, atom_symbols_loc,
     +    atoms_cart_loc, gamma_only_loc, spinors_loc, nntot_loc,
     +    nnlist_loc, nncell_loc, num_bands_loc, num_wann_loc,
     +    proj_site_loc, proj_l_loc, proj_m_loc, proj_radial_loc,
     +    proj_z_loc, proj_x_loc, proj_zona_loc, exclude_bands_loc
#   ifdef WANv2
     +    ,proj_s,proj_s_qaxis
#   endif
     +    )
         implicit none
         integer, parameter :: dp = selected_real_kind(15,300)
         integer, parameter :: num_nnmax=12
         character(len=*), intent(in) :: seed__name
         integer, dimension(3), intent(in) :: mp_grid_loc
         integer, intent(in) :: num_kpts_loc
         real(kind=dp), dimension(3,3), intent(in) :: real_lattice_loc
         real(kind=dp), dimension(3,3), intent(in) :: recip_lattice_loc
         real(kind=dp), dimension(3,num_kpts_loc), intent(in) ::
     +      kpt_latt_loc
         integer, intent(in) :: num_bands_tot
         integer, intent(in) :: num_atoms_loc
         character(len=*), dimension(num_atoms_loc), intent(in) ::
     +      atom_symbols_loc
         real(kind=dp), dimension(3,num_atoms_loc), intent(in) ::
     +      atoms_cart_loc
         logical, intent(in) :: gamma_only_loc
         logical, intent(in) :: spinors_loc
         integer, intent(out) :: nntot_loc
         integer, dimension(num_kpts_loc,num_nnmax), intent(out) ::
     +      nnlist_loc
         integer,dimension(3,num_kpts_loc,num_nnmax), intent(out) ::
     +      nncell_loc
         integer, intent(out) :: num_bands_loc
         integer, intent(out) :: num_wann_loc
         real(kind=dp), dimension(3,num_bands_tot), intent(out) ::
     +      proj_site_loc
         integer, dimension(num_bands_tot), intent(out) :: proj_l_loc
         integer, dimension(num_bands_tot), intent(out) :: proj_m_loc
         integer, dimension(num_bands_tot), intent(out) ::
     +      proj_radial_loc
         real(kind=dp), dimension(3,num_bands_tot), intent(out) ::
     +      proj_z_loc
         real(kind=dp), dimension(3,num_bands_tot), intent(out) ::
     +      proj_x_loc
         real(kind=dp), dimension(num_bands_tot), intent(out) ::
     +      proj_zona_loc
         integer, dimension(num_bands_tot), intent(out) ::
     +      exclude_bands_loc
#   ifdef WANv2
         integer, dimension(num_bands_tot), optional,intent(out) :: proj_s
         real(kind=dp), dimension(3,num_bands_tot), optional, intent(out) :: proj_s_qaxis
#   endif
        end subroutine wannier_setup
        subroutine wannier_run(seed__name,mp_grid_loc,num_kpts_loc,
     +    real_lattice_loc,recip_lattice_loc,kpt_latt_loc,num_bands_loc,
     +    num_wann_loc,nntot_loc,num_atoms_loc,atom_symbols_loc,
     +    atoms_cart_loc,gamma_only_loc,M_matrix_loc,A_matrix_loc,
     +    eigenvalues_loc,
     +    U_matrix_loc,U_matrix_opt_loc,lwindow_loc,wann_centres_loc,
     +    wann_spreads_loc,spread_loc)
         implicit none
         integer, parameter :: dp = selected_real_kind(15,300)
         character(len=*), intent(in) :: seed__name
         integer, dimension(3), intent(in) :: mp_grid_loc
         integer, intent(in) :: num_kpts_loc
         real(kind=dp), dimension(3,3), intent(in) :: real_lattice_loc
         real(kind=dp), dimension(3,3), intent(in) :: recip_lattice_loc
         real(kind=dp), dimension(3,num_kpts_loc), intent(in) ::
     +         kpt_latt_loc
         integer, intent(in) :: num_bands_loc
         integer, intent(in) :: num_wann_loc
         integer, intent(in) :: nntot_loc
         integer, intent(in) :: num_atoms_loc
         character(len=*), dimension(num_atoms_loc), intent(in) ::
     +         atom_symbols_loc
         real(kind=dp), dimension(3,num_atoms_loc), intent(in) ::
     +         atoms_cart_loc
         logical, intent(in) :: gamma_only_loc
         complex(kind=dp), dimension(num_bands_loc,num_bands_loc,
     +         nntot_loc,num_kpts_loc), intent(in) :: M_matrix_loc
         complex(kind=dp),
     +      dimension(num_bands_loc,num_wann_loc,num_kpts_loc),
     +      intent(in) :: A_matrix_loc
         real(kind=dp), dimension(num_bands_loc,num_kpts_loc),
     +      intent(in) :: eigenvalues_loc
         complex(kind=dp),
     +      dimension(num_wann_loc,num_wann_loc,num_kpts_loc),
     +      intent(out) :: U_matrix_loc
         complex(kind=dp),
     +      dimension(num_bands_loc,num_wann_loc,num_kpts_loc),
     +      optional, intent(out) :: U_matrix_opt_loc
         logical, dimension(num_bands_loc,num_kpts_loc),
     +      optional, intent(out) :: lwindow_loc
         real(kind=dp), dimension(3,num_wann_loc),
     +      optional, intent(out) :: wann_centres_loc
         real(kind=dp), dimension(num_wann_loc), optional,
     +      intent(out) :: wann_spreads_loc
         real(kind=dp), dimension(3), optional,
     +      intent(out) :: spread_loc
       end subroutine wannier_run
      end interface
# endif
      
      beginSingle

      allocate(efrozen(nspin1))
      lauto  = .false.
      second = .false.

      write(6,'(//A)') '### subroutine: wannier ### '

# if degstates != 1
      write(6,'(/A,I2)') 'Treatment of degenerate states (see wannier.f):',degstates
# endif

      ! Input from spex.inp
      call getkey(inp,'WBLOCH',    wbloch, section='WANNIER', default=.false.)
      call getkey(inp,'UREAD',  l_wanread, section='WANNIER', default=.false.)
      call getkey(inp,'INTERPOL', l_inter, section='WANNIER', default=.false.)
      call getkey(inp,'MAXIMIZE',l_wanmax, section='WANNIER', default=.false.)
      call getkey(inp,'FROZEN',   efrozen, section='WANNIER', status=i, writeout=.false., allow_eV=.true.)      
      l_frozen = i/=0
      if(i==1) efrozen = huge(0d0)
      if(l_wanmax.or.l_frozen) then
# ifndef WAN
        Error('MAXIMIZE or FROZEN specified, but code has been compiled without Wannier90 support (-DWAN).')
# endif
# ifdef MPI
        Info('Wannier90 library is not parallelized. Will be run in serial.')
# endif
      endif
      if(l_wanread) then
        write(6,'(/A)') 'U matrices for Wannier functions are read in.'
        write(6,'( A)') '(Make sure that the symmetry-equivalent wave functions are related by the same symmetry operations!)'
        call read_wannier
        call checkmem('uwan',16d0*size(uwan))
        write(6,'(A,F6.1," MB")') 'Size of Wannier U matrix:',size(uwan)*16d0 / megabyte
        goto 1
      endif

      backf = 0
      call getkey(inp,'BACKFOLD',charr,section='WANNIER',status=i)
      if(i==1) then
        Error('Arguments missing after BACKFOLD')
      else if(i==0) then

        ! Automatic: Determine integer "backfolding" vectors (backf) that bring atoms (cent+backf) as close as possible to each other
        call backfold(lat1,lat)
        lat1  = matmul(invert(lat1),lat)
        ilat  = nint( lat1 ) ; if(sum(abs(ilat-lat1))>1d-12) Bug('ilat not integer.')
        dist0 = sum ( [ ((matmul(lat,cent(:,i)-cent(:,ic))**2, i=1,ncent), ic=1,ncent) ] ) ! sum all pairwise distances
        backf = 0
        m     = 0
        do
          dist = [ (sum( [ (matmul(lat,cent(:,i)+backf(:,i)-cent(:,ic)-backf(:,ic))**2, i=1,ncent) ] ), ic=1,ncent) ]
          call rorderp(pnt,dist,ncent)
          ldum = .false.
          do i = ncent,1,-1
            ic     = pnt(i)
            backf1 = 0
            dist1  = dist(ic)
            do i1 = -nlat_backf , nlat_backf
            do i2 = -nlat_backf , nlat_backf
            do i3 = -nlat_backf , nlat_backf
              ivec = matmul( ilat,[i1,i2,i3] )
              rdum = 0
              do j = 1,ncent
                if(j/=ic) rdum = rdum + sum(matmul(lat,cent(:,j)+backf(:,j)-cent(:,ic)-backf(:,ic)-ivec)**2)
              enddo
              if(rdum<dist1-1d-4) then
                backf1 = ivec
                dist1  = rdum
              endif
            enddo
            enddo
            enddo
            ldum        = ldum .or. any(backf1/=0)
            backf(:,ic) = backf(:,ic) + backf1
          enddo
          if(ldum) then ; m = m + 1 ; if(m>2*ncent) Error('Backfolding vectors did not converge.')
          else          ; exit
          endif
        enddo

        if(any(backf/=0)) then ! undo a global shift (by subtracting the average lattice vector)
          do i = 1,3
            backf(i,:) = backf(i,:) - nint ( dble(sum(backf(i,:))-1d-12) / ncent )
          enddo
        endif
        if(any(backf/=0)) then          
          dist1 = sum ( [ ((matmul(lat,cent(:,i)+backf(:,i)-cent(:,ic)-backf(:,ic))**2, i=1,ncent), ic=1,ncent) ] )
          write(6,'(/A)') 'Atom backfolding vectors:'
          write(6,'(3I3)') backf
          write(6,'(2(A,F8.3))') 'Average distance',sqrt(dist0/(ncent*(ncent-1))),' reduced to',sqrt(dist1/(ncent*(ncent-1)))
          if(dist1>dist0) Error('Average distance increased.')
        endif

      else if(i==2) then

        ! Manual backfolding
        if(size(charr)>ncent) Error('Too many arguments after BACKFOLD.')
        do i = 1,size(charr)
          if(charr(i)/='()') then
            call parser(parse_out,charr(i),'(1,2,3)','BACKFOLD')
            do j = 1,3 ; backf(j,i) = ch2i(parse_out(j),'BACKFOLD') ; enddo
          endif
        enddo
        deallocate ( charr )

      endif

      write(6,*)

      ! Define orbitals and energy window
      call getkey(inp,'ORBITALS',   charr, section='WANNIER', status=i)
      if     (i==0) then ; Error('Missing ORBITALS definition.')
      else if(i==1) then ; Error('Missing arguments after ORBITALS.')
      endif
      call orbitals(lwan,lmwan,centwan,eulerwan,nwan,charr(3),size(charr)-2)
      centwan = abs(centwan)
      do i = 1,nwan
        do j = 1,i-1
          if(lmwan(j)==lmwan(i).and.centwan(j)==centwan(i)) Error('Double orbital definition.')
        enddo
      enddo
      if(l_soc) then
        call reallocate(lwan,      2*nwan) ; lwan      (nwan+1:2*nwan) = lwan      (:nwan)
        call reallocate(lmwan,     2*nwan) ; lmwan     (nwan+1:2*nwan) = lmwan     (:nwan)
        call reallocate(centwan,   2*nwan) ; centwan   (nwan+1:2*nwan) = centwan   (:nwan)        
        call reallocate(eulerwan,3,2*nwan) ; eulerwan(:,nwan+1:2*nwan) = eulerwan(:,:nwan)
        nwan = nwan * 2
      endif
      wanbandi = ch2i(charr(1),'ORBITALS')
      lauto    = charr(2)=='*'
      if(lauto) then ; wanbandf = wanbandi + nwan - 1 ; write(6,'(/A/)') 'Upper band index set to '//Chr(wanbandf)
      else           ; wanbandf = ch2i(charr(2),'ORBITALS')
      endif
      deallocate(charr)

 100  nwanband = wanbandf - wanbandi + 1 ; if(nwanband<0) Error('Second argument must be larger than first.')
      if(wanbandi<1)             Error('Lowest band index must be positive.')
      if(wanbandf>minval(nband)) Error('Highest band index for Wannier construction exceeds max. band index.')
      if(nwanband<nwan)          Error('Less bands than Wannier functions.')

      wanbandi0 = wanbandi
      wanbandf0 = wanbandf

# if degstates == 1
c
c     Redefine wanbandi and wanbandf if necessary
      do ispin = 1,nspin1
        do ikpt = 1,nkpt2
          i = deg(wanbandi0,ikpt,ispin) ; if(i<wanbandi0) wanbandi = min(wanbandi,i)
          i = deg(wanbandf0,ikpt,ispin) ; if(i<wanbandf0) i        = deg(i,ikpt,ispin) ; wanbandf = max(wanbandf,i)
        enddo
      enddo
# endif

      nwanband = wanbandf - wanbandi + 1

      ! Read spex.uwan if requested
      ldum = iand(restart,R_uwan)/=0
      if(ldum) inquire(file='spex.uwan',exist=ldum)
      if(ldum) then
        iunit = fopen('spex.uwan',form='unformatted',status='old')
        read(iunit,iostat=ios) vers
        if      (ios/=0) then ; Error('Read error (version number, "spex.uwan")')
        else if(vers>-1) then ; Error('File format of "spex.uwan" too old. Cannot be read.')
        else if(vers<-2) then ; Error('File format of "spex.uwan" too new. Cannot be read.')
        endif
        read(iunit,iostat=ios) i,j,kvec,rdum
        if(ios/=0) Error('Read error (nkpt, nspin, kptadd, checksum, "spex.uwan")')
        if(i/=nkpt)   Error('Number of k points in "spex.uwan" differs from present one: '//Chr(nkpt))
        if(j/=nspin1) Error('Number of spins in "spex.uwan" differs from present one: '//Chr(nspin1))
        if(sum(abs(kvec-kptadd))>1d-8)Error('Additional (+) k point in "spex.uwan" differs from present one: ('//Chn(kvec,',')//')')
        read(iunit,iostat=ios) i1,i2,i3,m,n
        if(ios/=0) Error('Read error (parameters, "spex.uwan")')
        if(any([nwan,wanbandi,wanbandf]/=[i1,i2,i3])) then
          write(6,'(A)')  'Dimensions of Wannier U matrix set to the ones in "spex.uwan":'
          write(6,'(A)')  '- Number of Wannier functions: '//Chr(i1)   //' instead of '//Chr(nwan)
          write(6,'(A/)') '- Wannier band window: '//Chn((/i2,i3/),'-')//' instead of '//Chn((/wanbandi,wanbandf/),'-')
          Warn('Dimensions of Wannier U matrix set to the ones in "spex.uwan".')
          nwan = i1 ; wanbandi = i2 ; wanbandf = i3 ; nwanband = i3 - i2 + 1
        endif
        if(any([wanbandi0,wanbandf0]/=[m,n])) then
          write(6,'(A)') 'Wannier U matrix in "spex.uwan" was generated from ORBITALS '//Chr(m)//' '//Chr(n)
          wanbandi0 = m
          wanbandf0 = n
        endif
        if(vers==-1) then
          call get_checksum(chksum,rdum1,1,maxband,[(i,i=1,nkpti),(i,i=nkpt+1,nkpt+nkpti2)],nkpti+nkpti2,0)
          if(abs(rdum-chksum)>acc_checksum*sqrt(dble(size(cmt)+size(cpw))))
     &      Error('Checksum failure. Current wave functions inconsistent with stored Wannier coefficients. Remove spex.uwan.')
        else if(vers==-2) then
          call get_checksum(chksum,rdum1,wanbandi,wanbandf,[(i,i=1,nkpti),(i,i=nkpt+1,nkpt+nkpti2)],nkpti+nkpti2,0)
          if(abs(rdum-chksum)>acc_checksum*sqrt(dble(nwanband*(nkpti+nkpti2)*nspin2)))
     &      Error('Checksum failure. Current wave functions inconsistent with stored Wannier coefficients. Remove spex.uwan.')          
        endif
        allocate(wancent(3,nwan,nspin1*(nkpt2/nkpt)))
        Allocate_(uwan,(wanbandi:wanbandf,nwan,nkpt2,nspin1) )
        read(iunit,iostat=ios) wancent ; if(ios/=0) Error('Read error (wancent, "spex.uwan").')
        read(iunit,iostat=ios) uwan    ; if(ios/=0) Error('Read error (uwan, "spex.uwan").')
        call fclose(iunit)
        write(6,'(A)') 'Wannier U matrix read from "spex.uwan".'
        goto 10
      endif

      allocate ( mat2(nwan,nwan),eig(nwan) )
      if(lkptadd) then ; allocate(wancent(3,nwan,nspin1*2),fwancent(3,nwan,nspin1*2),decomp(nwan,nwan,nspin1*2))
      else             ; allocate(wancent(3,nwan,nspin1),  fwancent(3,nwan,nspin1),  decomp(nwan,nwan,nspin1)  )
      endif

      Allocate_(uwan,(wanbandi:wanbandf,nwan,nkpt2,nspin1))
      write(6,'(A,F6.1," MB")') 'Size of Wannier U matrix:',size(uwan)*16d0 / megabyte
      write(6,'(A)')            'Wannier band window: '//Chr(wanbandi)//'-'//Chr(wanbandf)
      write(6,'(/A)')           'Wannier construction / SVD eigenvalues'

 1    allocate ( cmt1(maxlmindx,ncent,wanbandi:wanbandf) )
      allocate ( fwan(wanbandi:wanbandf,nwan) )      

      ! Set wancent
      do ispin = 1,nspin1
        do iwan = 1,nwan
          ic = centwan(iwan)
          if(second) then
# if   set_wancent == 1
            wancent(:,iwan,ispin) = cent(:,ic) + backf(:,ic) ! overwrite wancent with first-guess centers
# elif set_wancent == 2
            if(sum(matmul(lat,cent(:,ic) + backf(:,ic) - wancent(:,iwan,ispin))**2)<1)
     &      wancent(:,iwan,ispin) = cent(:,ic) + backf(:,ic) ! overwrite wancent if Wannier90 has moved center less than 1 Bohr away from first-guess center
# endif
          else
            wancent(:,iwan,ispin) = cent(:,ic) + backf(:,ic) ! always set wancent to first-guess center in first round
          endif
          if(lkptadd) wancent(:,iwan,nspin1+ispin) = wancent(:,iwan,ispin)
        enddo
      enddo

      nsmall1 = 0
      nsmall2 = 0

      do ispin = 1,nspin1
        do ikpt = 1,nkpt2

          if(l_wanread) cycle

          if(.not.second) uwan(:,:,ikpt,ispin) = 0

c
c         Calculate first-guess Wannier functions

          i1 = wanbandi0
          i2 = wanbandf0
# if degstates == 2
          ! Cut states
          do while(deg(i1,ikpt,ispin)<i1.and.i1/=nband(ikpt,ispin))
            i1 = i1 + 1
          enddo
          if(i2/=nband(ikpt,ispin)) then
            if(deg(i2+1,ikpt,ispin)<i2+1) then
              i2 = deg(i2+1,ikpt,ispin) - 1
            endif
          endif
          if(i1>wanbandi.or.i2<wanbandf.and.ikpt<=nkpti) then
            if(.not.second) write(6,'(I1,I6,A,2I4)') ispin,ikpt,'  : ',i1,i2
            if(i1>i2)        Error('Zero bands after cutting of degenerate states.')
            if(i2-i1+1<nwan) Error('Less states than Wannier functions.')
          endif
# elif degstates == 1
          i = deg(i1,ikpt,ispin) ; if(i<i1) i1 = i
          i = deg(i2,ikpt,ispin) ; if(i<i2) i  = deg(i,ikpt,ispin) ; i2 = i
          if(i1<wanbandi0.or.i2>wanbandf0.and.ikpt<=nkpti.and..not.second) then
            write(6,'(I1,I6,A,2I4)') ispin,ikpt,'  : ',i1,i2
          endif
# endif

c
c         First guess: Project onto MT functions (->uwan)
          do iwan = 1,nwan
            ! define projection function in harmonics (->proj)
            ic = centwan(iwan)
            lm = abs(lmwan(iwan))
            l  = sqrt(lm*1d0) + 1d-12
            if(lm>100) then ; allocate ( proj(-2:2,0:2) )
            else            ; allocate ( proj(-l:l,l:l) )
            endif
            proj = 0
            do l = lbound(proj,2),ubound(proj,2)
              call add_proj(proj(-l:l,l),l,lmwan(iwan),eulerwan(:,iwan))
            enddo
            ! Get wavefunctions at current k point (->cmt1)
            if(iwan==1.or.(l_soc.and.iwan==nwan/2+1)) then
              if(iwan==1) then ; s = ispin
              else             ; s = 2
              endif
# ifdef LOAD
              allocate(cmt(maxlmindx,ncent,i1:i2,1,nspin3))
              call read_wavef2([(i,i=i1,i2)],i2-i1+1,ikpt,ispin,cmt)
              do i = i1,i2
                cmt1(:,:,i) = cmt(:,:,i,1, min(nspin3,s) )
              enddo
              deallocate(cmt)
# else
              do i = i1,i2
#   ifndef old_trafo
                if(l_soc) then ; call wavefunction_mt(cmt2,       maxlmindx,0,i,ikpt,ispin) ; cmt1(:,:,i) = cmt2(:,:,s)
                else           ; call wavefunction_mt(cmt1(:,:,i),maxlmindx,0,i,ikpt,ispin)
                endif
#   else
                if(storeibz.and.kptp(ikpt)/=ikpt) then
                  do j = 1,ncent
                    call waveftrafo_mt(cmt1(:,j,i),maxlmindx,ikpt,i,s,j,.false.)
                  enddo
                else
                  cmt1(:,:,i) = cmt(:,:,i,kindx(ikpt),s)
                endif
#   endif
              enddo
# endif
            endif
            ! project (->uwan)
            do itype = 1,ntype
              if(sum(neq(:itype))>=ic) exit
            enddo
            sb = min(s,nspin)
            do i = i1,i2
              cdum = 0
              do l = lbound(proj,2),ubound(proj,2)
                lm = sum ( [ ((2*l1+1)*nindx(l1,itype),l1=0,l-1) ] )
                do m = -l,l
                  do n = 1,nindx(l,itype)
                    lm   = lm + 1 ; if(proj(m,l)==0) cycle
                    cdum = cdum + conjg ( cmt1(lm,ic,i) ) * proj(m,l) *
     &                     intgrf ( bas1(:,n,l,itype,sb)*bas1(:,1,l,itype,sb) +
     &                              bas2(:,n,l,itype,sb)*bas2(:,1,l,itype,sb) ,itype)
                  enddo
                enddo
              enddo
              cdum = cdum * exp( -img * 2*pi * dot_product(kpt(:,ikpt),backf(:,ic)) )
              if(second) then ; soff = 0 ; if(ikpt>nkpt) soff = nspin1
                decomp(iwan,:,soff+ispin) = decomp(iwan,:,soff+ispin) + conjg(cdum) * uwan(i,:,ikpt,ispin) / nkpt
              endif
              fwan(i,iwan) = cdum
            enddo
            deallocate ( proj )
          enddo

          ! Make uwan unitary (singular value decomposition)
          allocate ( mat1(i1:i2,i1:i2) )
          call svd(mat1,eig,mat2,fwan(i1:i2,:))
          fwan(i1:i2,:) = matmat(mat1(:,i1:i1+nwan-1),mat2)
          deallocate ( mat1 )

          if(second) then
# ifdef reproject
            uwan(i1:i2,:,ikpt,ispin) = matmul ( fwan(i1:i2,:) , matmul (
     &        transpose(conjg(fwan(i1:i2,:))) , uwan(i1:i2,:,ikpt,ispin) ) )
            allocate ( mat1(i1:i2,i1:i2) )
            call svd(mat1,eig,mat2,uwan(i1:i2,:,ikpt,ispin))
            uwan(i1:i2,:,ikpt,ispin) = matmat(mat1(:,i1:i1+nwan-1),mat2)
            deallocate ( mat1 )
# endif
            cycle
          endif

          uwan(i1:i2,:,ikpt,ispin) = fwan(i1:i2,:)

          if(ikpt<=nkpti) then
            if(nkpti<200.or.modulo(ikpt,5)==0) then
              write(6,'(I1,I6'NoA) ispin,ikpt
              do i = 1,nwan,5
                if(i/=1) write(6,'(A'NoA) '       '
                write(6,'(5F14.6)') eig(i:min(i+4,nwan))
              enddo
            endif
            if(any(abs(eig)<mwarn1)) then
              if(.not.lauto) Warn('Very small SVD eigenvalue(s) (<'//Chf(mwarn1,'ES8.1')//'). Increase upper band index.')
              nsmall1 = nsmall1 + count(abs(eig)<mwarn1)
              write(6,'(A'NoA) 'Found very small eigenvalue(s):'
              do i = 1,nwan
                if(abs(eig(i))<mwarn1) write(6,'(ES8.1'NoA) eig(i)
              enddo
              write(6,*)
            endif
            if(any(abs(eig)<mwarn2.and.abs(eig)>=mwarn1)) then
              if(.not.lauto) Warn('Small SVD eigenvalue(s) (<'//Chf(mwarn2,'ES8.1')//'). Consider increasing upper band index.')
              nsmall2 = nsmall2 + count(abs(eig)<mwarn2.and.abs(eig)>=mwarn1)
              write(6,'(A'NoA) 'Small eigenvalue(s):'
              do i = 1,nwan
                if(abs(eig(i))<mwarn2.and.abs(eig(i))>=mwarn1) write(6,'(ES8.1'NoA) eig(i)
              enddo
              write(6,*)
            endif
            if(any(abs(eig)<mwarn2).and.lauto) then
              wanbandf = wanbandf0 + 1
              Deallocate_(uwan)
              deallocate(fwan,cmt1,mat2,eig,wancent,fwancent,decomp)
              write(6,'(/A/)') 'Upper band index increased to '//Chr(wanbandf)
              goto 100                
            endif
            if(any(eig>1)) Warn('Eigenvalue greater than one. Orbitals may be overcomplete.')
          endif

        enddo
      enddo

      deallocate ( cmt1,fwan )
      if(second) goto 2

      if(lauto) write(6,'(/A)') 'We use: ORBITALS '//Chr(wanbandi0)//' '//Chr(wanbandf0)//' ...'

      if(nsmall1>0.or.nsmall2>0) write(6,*)
      if(nsmall1>0) then
        write(6,'(A)') Chr(nsmall1)//' SVD eigenvalues smaller than '//Chf(mwarn1,'ES8.1')//', see warnings. '//
     &                 'Upper energy window should be increased.'
      endif
      if(nsmall2>0) then
        if(nsmall1>0) write(6,'(A'NoA) 'In addition, '
        write(6,'(A)') Chr(nsmall2)//' SVD eigenvalues smaller than '//Chf(mwarn2,'ES8.1')//', see warnings.'
      endif
      
      ! test if Wannier set is closed (needed for STOREIBZ, see wavefproducts)
      dev = 0
      do ispin = 1,nspin1
        s = min(ispin,nspin1)
        do ikpt = nkpti+1,nkpt
          if(symkpt(ikpt)>nsymt) then
            mat2 = matmul ( transpose(conjg(uwan(:,:,ikpt,ispin))) , conjg(uwan(:,:,kptp(ikpt),ispin)) )
          else
            mat2 = matmul ( transpose(conjg(uwan(:,:,ikpt,ispin))) ,       uwan(:,:,kptp(ikpt),ispin)  ) * phase(ikpt)
          endif
          dev = max(dev,sum(abs(identity(nwan)-matmul(conjg(transpose(mat2)),mat2))))
        enddo
      enddo
      if(dev>1d-8) then
        if(storeibz.and.wbloch) then
          write(0,'(A,F11.8)') 'wannier: Wannier set not closed wrt symmetry operations. Error:',dev
          write(0,'(A)')       '         Check if the Wannier functions are defined for all equivalent atoms.'
          Error('The Wannier set is not closed wrt the symmetry operations. WBLOCH possible only with STOREBZ.')
        else
          Info('The Wannier set is not closed wrt the symmetry operations. Error: '//Chf(dev,'F11.8'))
        endif
      endif

# ifdef WANv2
c        iunit = fopen('wannier.dmn')
c        write(6,'(A)') 'Symmetry file created by SPEX.'
c        write(6,'(4I5)') nwanband,nsym,nkpti,nkpt
c        write(6,*)
c        write(6,'(10I5)') kptp(:nkpt)
c        write(6,*)
c        write(6,'(10I5)') (i,i=1,nkpti)
c        write(6,*)
c        do ikpt = 1,nkpti ; do isym = 1,nsym
c          write nwan
c        enddo ; enddo
c        write(6,*)
c        do ikpt = 1,nkpti ; do isym = 1,nsym
c          write nwanbnad
c        enddo ; enddo
c        call fclose(iunit)
# endif

      ! Determine first-guess Wannier centers (should go here...)

      ! Maximal localization
# ifdef WAN
      if(l_wanmax.or.l_frozen) then
        write(6,'(/A)') 'Construct MLWFs'
        ldum = iand(restart,R_kolap)/=0
        if(ldum) inquire(file='spex.kolap',exist=ldum)
        allocate ( zlabel1(ncent),cent1(3,ncent),nnlist(nkpt,12),nncell(3,nkpt,12) )
        allocate ( proj_l(nwanband),proj_m(nwanband),proj_radial(nwanband),exclude_bands(nwanband) )
        allocate ( proj_site(3,nwanband),proj_z(3,nwanband),proj_x(3,nwanband),proj_zona(nwanband) )
        zlabel1 = [ ((zlabel(ztype(itype)),ieq=1,neq(itype)),itype=1,ntype) ]
        cent1   = matmul(lat,cent) * angstrom
        do ispin = 1,nspin1
          write(6,*)
          if(nspin1/=1) then
            if(ispin==1) write(6,'(A)') 'Spin up'
            if(ispin==2) write(6,'(A)') 'Spin down'
          endif
          call write_wannier_win(ispin)
          off  = 0                       ! set to 1 for add. kpoint
 4        koff = off*nkpt                ! =0    unless add. kpoint
          soff = off*nspin1              ! =0    unless add. kpoint
          kpt1 = kpt(:,koff+1:koff+nkpt) ! =kpt  unless add. kpoint
          if(off==1) then
            write(6,'(A)') 'Shifted k-point set...'
            do ikpt = 1,nkpt ; kpt1(:,ikpt) = kpt1(:,ikpt) - kptadd ; enddo
          endif
          write(6,'(A'NoA) 'Wannier setup... '
          call cpu_time(cputime)
          call wannier_setup('wannier',nkpt3,nkpt,transpose(lat)*angstrom,transpose(rlat)/angstrom,kpt1,nwanband,
     &      ncent,zlabel1,cent1,.false.,l_soc,nntot,nnlist,nncell,i1,i2,proj_site,proj_l,proj_m,proj_radial,
     &      proj_z,proj_x,proj_zona,exclude_bands)
          if(i1/=nwanband) Error('wannier_setup returned different number of bands.')
          if(i2/=nwan)     Error('wannier_setup returned different number of Wannier functions.')
          allocate ( kolap(wanbandi:wanbandf,wanbandi:wanbandf,nntot,nkpt),umat(nwan,nwan,nkpt) )
          allocate ( uwan1(wanbandi:wanbandf,nwan,nkpt),lwindow(nwanband,nkpt),wann_spreads(nwan) )
          call cpu_done(cputime)
          if(.not.ldum) then
            write(6,'(A'NoA) 'Calculate overlap matrices... '
            call checkmem('cmt_ and cpw_',(MBYTES*maxgpt+16*maxlmindx*ncent)*nwanband*nkpt*nspin3)
            allocate ( cmt_(maxlmindx,ncent,wanbandi:wanbandf,nkpt,nspin3) )
            allocate ( cpw_(maxgpt,wanbandi:wanbandf,nkpt,         nspin3) )
# ifdef LOAD
            call read_wavef0([(i,i=wanbandi,wanbandf)],[(koff+ikpt,ikpt=1,nkpt)],ispin,cmt,cpw)
# else
            if(storeibz) then
              if(off==1) then ; m = nkpti ; n = nkpti2
              else            ; m = 0     ; n = nkpti
              endif
              cmt_(:,:,:,:n,:) = cmt(:,:,wanbandi:wanbandf,m+1:m+n,ispin:max(ispin,nspin3))
              cpw_(:,  :,:n,:) = cpw(:,  wanbandi:wanbandf,m+1:m+n,ispin:max(ispin,nspin3))
              do ikpt = off*nkpt+n+1,(1+off)*nkpt ! loop over equivalents
                ikpt1 = ikpt - off*nkpt
                do i = wanbandi,wanbandf
# ifndef old_trafo
                  call wavefunction_mt(cmt_(:,:,i,ikpt1,:),maxlmindx,0,i,ikpt,ispin)
                  call wavefunction_pw(cpw_(:,  i,ikpt1,:),maxgpt,     i,ikpt,ispin)                  
# else
                  do ic = 1,ncent
                    call waveftrafo_mt(cmt_(:,ic,i,ikpt1,1),maxlmindx,ikpt,i,ispin,ic,.false.)
                    if(l_soc) then
                      call waveftrafo_mt(cmt_(:,ic,i,ikpt1,2),maxlmindx,ikpt,i,2,  ic,.false.)
                      call waveftrafo_soc(cmt_(:,ic,i,ikpt1,:),maxlmindx,symkpt(ikpt))
                    endif
                  enddo
                  call waveftrafo_pw(cpw_(:,i,ikpt1,1),ikpt,i,ispin,.false.)
                  if(l_soc) then
                    call waveftrafo_pw(cpw_(:,i,ikpt1,2),ikpt,i,2,  .false.)
                    call waveftrafo_soc(cpw_(:,i,ikpt1,:),maxgpt,symkpt(ikpt))
                  endif
# endif
                enddo
              enddo
            else
              cmt_(:,:,:,:nkpt,:) = cmt(:,:,wanbandi:wanbandf,koff+1:koff+nkpt,ispin:max(ispin,nspin3))
              cpw_(:,  :,:nkpt,:) = cpw(:,  wanbandi:wanbandf,koff+1:koff+nkpt,ispin:max(ispin,nspin3))
            endif
# endif
            allocate ( neigh(nkpt,nkpt) )
            neigh = 0
            do j = 1,nntot
              kvec = kpt(:,koff+nnlist(1,j)) - kpt(:,koff+1) + nncell(:,1,j)
              call wfolap_init_mt(olapmt,kvec)
              do ikpt = 1,nkpt
                allocate ( cpw2(ngpt(koff+ikpt)) )
                do i = 1,nntot
                  ibpt  = nnlist(ikpt,i)
                  kvec0 = kpt(:,koff+ibpt) - kpt(:,koff+ikpt) + nncell(:,ikpt,i)
                  if(any(abs(kvec0-kvec)>1d-12)) cycle
                  if(neigh(ikpt,ibpt)/=0) then
                    kolap(:,:,i,ikpt) = kolap(:,:,neigh(ikpt,ibpt),ikpt)
                    cycle
                  else if(neigh(ibpt,ikpt)/=0) then
                    kolap(:,:,i,ikpt) = conjg(transpose(kolap(:,:,neigh(ibpt,ikpt),ibpt)))
                    cycle
                  endif
                  neigh(ikpt,ibpt) = i
                  allocate ( olappw(ngpt(koff+ikpt),ngpt(koff+ibpt)) )
                  if(any(nncell(:,ikpt,i)/=0)) then
                    do i2 = 1,ngpt(koff+ibpt) ; do i1 = 1,ngpt(koff+ikpt)
                      olappw(i1,i2) = cstep ( gpt(1,pgpt(i1,koff+ikpt))-gpt(1,pgpt(i2,koff+ibpt))+nncell(1,ikpt,i) ,
     &                                        gpt(2,pgpt(i1,koff+ikpt))-gpt(2,pgpt(i2,koff+ibpt))+nncell(2,ikpt,i) ,
     &                                        gpt(3,pgpt(i1,koff+ikpt))-gpt(3,pgpt(i2,koff+ibpt))+nncell(3,ikpt,i) )
                    enddo ; enddo
                  else
                    call olap_gpt(olappw,ngpt(koff+ikpt),koff+ikpt,koff+ibpt)
                  endif
                  do i2 = wanbandi,wanbandf
                    call wfolap_mt2(cmt2,olapmt,cmt_(:,:,i2,ibpt,1),ispin,.false.)
                    cpw2 = matmul(olappw,cpw_(:ngpt(koff+ibpt),i2,ibpt,1))
                    do i1 = wanbandi,wanbandf
                      kolap(i1,i2,i,ikpt) =   dot_product_mt(cmt_(:,:,i1,ikpt,1),cmt2)
     &                                      + dot_product(cpw_(:ngpt(koff+ikpt),i1,ikpt,1),cpw2)
                    enddo
                  enddo
                  if(l_soc) then
                    do i2 = wanbandi,wanbandf
                      call wfolap_mt2(cmt2,olapmt,cmt_(:,:,i2,ibpt,2),nspin,.false.)
                      cpw2 = matmul(olappw,cpw_(:ngpt(koff+ibpt),i2,ibpt,2))
                      do i1 = wanbandi,wanbandf
                        kolap(i1,i2,i,ikpt) = kolap(i1,i2,i,ikpt)
     &                                        + dot_product_mt(cmt_(:,:,i1,ikpt,2),cmt2)
     &                                        + dot_product(cpw_(:ngpt(koff+ikpt),i1,ikpt,2),cpw2)
                      enddo
                    enddo
                  endif
                  deallocate ( olappw )
                enddo
                deallocate ( cpw2 )
              enddo
            enddo
            deallocate ( neigh )
            deallocate ( cmt_,cpw_ ) ; call checkmem('cmt_ and cpw_',-(MBYTES*maxgpt+16*maxlmindx*ncent)*nwanband*nkpt*nspin3)
            call cpu_done(cputime)
            if(iand(restart,W_kolap)/=0) then
              iunit = fopen('spex.kolap',form='unformatted',status='unknown')
              n     = (ispin-1)*2 ; if(lkptadd) n = n * 2 ; if(off==1) n = n + 2 ! skip spin up and unshifted
              do i = 1,n ; read(iunit) ; enddo
              call get_checksum(chksum,rdum1,wanbandi,wanbandf,[(i,i=1,nkpti),(i,i=nkpt+1,nkpt+nkpti2)],nkpti+nkpti2,0)
              write(iunit) wanbandi,wanbandf,nntot,nkpt,ispin,off,kptadd,chksum,2
              write(iunit) kolap
              call fclose(iunit)
              write(6,'(A)') 'Overlap matrix written to spex.kolap.'
            endif
          else
            iunit = fopen('spex.kolap',form='unformatted',status='old')
            n     = (ispin-1)*2 ; if(lkptadd) n = n * 2 ; if(off==1) n = n + 2 ! skip spin up and unshifted
            do i = 1,n ; read(iunit) ; enddo
            read(iunit,iostat=ios) i1,i2,i3,i,j,l,kvec,rdum,vers
            vers = 2
            if(ios/=0) then
              backspace(iunit)
              read(iunit,iostat=ios) i1,i2,i3,i,j,l,kvec,rdum
              if(ios/=0.or.rdum==0) then
                vers = 0 ; Warn('spex.kolap version 0. No checksum. Cannot check consistency of wavefunctions.')
              else
                vers = 1 ; Warn('spex.kolap version 1. Only global wavefunction checksum can be checked.')
              endif
            endif
            if(vers==1) then
              call get_checksum(chksum,rdum1,1,maxband,[(i,i=1,nkpti),(i,i=nkpt+1,nkpt+nkpti2)],nkpti+nkpti2,0)
              if(abs(rdum-chksum)>acc_checksum*sqrt(dble(size(cmt)+size(cpw))))
     &          Error('Checksum failure. Current wave functions inconsistent with stored Wannier coefficients. Remove spex.kolap.')
            else if(vers==2) then
              call get_checksum(chksum,rdum1,wanbandi,wanbandf,[(i,i=1,nkpti),(i,i=nkpt+1,nkpt+nkpti2)],nkpti+nkpti2,0)
              if(abs(rdum-chksum)>acc_checksum*sqrt(dble(nwanband*(nkpti+nkpti2)*nspin2)))
     &          Error('Checksum failure. Current wave functions inconsistent with stored Wannier coefficients. Remove spex.kolap.')
            endif
            if(vers==0) then
              backspace(iunit)
              read(iunit,iostat=ios) i1,i2,i3,i,j,l,kvec
              if(ios/=0) Error('Read error (parameters, spex.kolap).')
            endif 
            if(any([i1,i2,i3,i,j,l]/=[wanbandi,wanbandf,nntot,nkpt,ispin,off]).or.sum(abs(kvec-kptadd))>1d-8) then
              write(0,'(A,3I4,I6,2I2,3F9.5,A/A)') 'wannier: Parameters changed (spex.kolap:',i1,i2,i3,i,j,l,kvec,')',
     &                                            '         Remove spex.kolap.'
              Error('Parameters in spex.kolap inconsistent with present ones.')
            endif
            read(iunit,iostat=ios) kolap ; if(ios/=0) Error('Read error (kolap, spex.kolap)')
            call fclose(iunit)
            write(6,'(A)') 'Overlap matrix read from spex.kolap.'
          endif
          write(6,'(A'NoA) 'Wannierize... '
          call wannier_run('wannier',nkpt3,nkpt,transpose(lat)*angstrom,transpose(rlat)/angstrom,kpt1,nwanband,
     &      nwan,nntot,ncent,zlabel1,cent1,.false.,kolap,uwan(:,:,koff+1:,ispin),ene(wanbandi:wanbandf,koff+1:koff+nkpt,ispin)*
     &      hartree,umat,uwan1,lwindow,wancent(:,:,soff+ispin),wann_spreads,kvec0)
          wancent(:,:,soff+ispin) = matmul( invert(lat) , wancent(:,:,soff+ispin) / angstrom ) ! transform to internal coordinates
          call cpu_done(cputime)
          do ikpt = 1,nkpt
            uwan(:,:,koff+ikpt,ispin) = matmul(uwan1(:,:,ikpt),umat(:,:,ikpt))
          enddo
          deallocate ( umat,kolap,uwan1,lwindow,wann_spreads )
          if(lkptadd.and.off==0) then
            off = 1
            goto 4
          endif
        enddo
        deallocate ( zlabel1,cent1,nnlist,nncell,proj_l,proj_m,proj_radial,exclude_bands,proj_site,proj_z,proj_x,proj_zona )
      endif
# endif

c     Orbital decomposition (-> decomp(iwan,iwan1) = <firstguess(iwan)|wannier(iwan1)> )
 3    if(.not.l_wanread) then
        second = .true.
        decomp = 0
        goto 1
      endif
 2    if(second) then
        write(6,'(//A)') 'Orbital decomposition 10^3 (projection onto first guess)'
        ic = centwan(1)
        do ispin = 1,nspin1
          if(nspin1==2.and.ispin==1) write(6,'(/A)') 'Spin up'
          if(nspin1==2.and.ispin==2) write(6,'(/A)') 'Spin down'
          do off = 0,1
            if(off==1) write(6,'(A)') 'Shifted k-point set'
            do iwan = 1,nwan
              if(centwan(iwan)/=ic.or.l_soc.and.iwan==nwan/2+1) then ; write(6,*) ; ic = centwan(iwan) ; endif
              do iwan1 = 1,nwan
                write(6,'(I4'NoA) nint(abs(decomp(iwan,iwan1,off*nspin1+ispin))*1000)
              enddo
              write(6,*) !sum(abs(decomp(:,iwan,ispin))**2),sum(abs(decomp(iwan,:,ispin))**2)
            enddo
            if(.not.lkptadd) exit
          enddo
        enddo
        write(6,'(/A'NoA) 'Wannier centers in internal coordinates '
        if(l_wanmax.or.l_frozen) then
# if   set_wancent == 0
          write(6,'(A)') '(as returned by Wannier90 library)'
# elif set_wancent == 1
          write(6,'(A)') '(set at first-guess centers)'
# else
          write(6,'(A)') '(set at first-guess center if Wannier90 shift < 1 Bohr)'
# endif
        else
          write(6,*)
        endif
        do ispin = 1,nspin1
          if(nspin1==2.and.ispin==1) write(6,'(/A)') 'Spin up'
          if(nspin1==2.and.ispin==2) write(6,'(/A)') 'Spin down'
          do off = 0,1
            if(off==1) write(6,'(A)') 'Shifted k-point set'
            do iwan = 1,nwan
              write(6,'(I3,3F10.5)') iwan,wancent(:,iwan,off*nspin1+ispin)
            enddo
            if(.not.lkptadd) exit
          enddo
        enddo
        write(6,*)
      endif

      if(iand(restart,W_uwan)/=0) then
        iunit = fopen('spex.uwan',form='unformatted',status='unknown')
        call get_checksum(chksum,rdum1,wanbandi,wanbandf,[(i,i=1,nkpti),(i,i=nkpt+1,nkpt+nkpti2)],nkpti+nkpti2,0)
        write(iunit) -2
        write(iunit) nkpt,nspin1,kptadd,chksum
        write(iunit) nwan,wanbandi,wanbandf,wanbandi0,wanbandf0
        write(iunit) wancent
        write(iunit) uwan
        call fclose(iunit)
        write(6,'(/A)') 'Wannier U matrix written to spex.uwan.'
      endif

# if 0
      write(6,'(/A)') 'Determine Wannier centers'
      allocate(dcprod(nwanband,nwanband,nkpt,3))
      fwancent = 0
      do ispin = 1,nspin1
        call dwavefproducts(dcprod,[(ikpt,ikpt=1,nkpt)],nkpt,ispin,ispin,wanbandi,wanbandf,wanbandi,wanbandf,.false.,.false.)
        if(l_soc) call dwavefproducts(dcprod,[(ikpt,ikpt=1,nkpt)],nkpt,2,2,wanbandi,wanbandf,wanbandi,wanbandf,.false.,.true.)
        do iwan = 1,nwan
          do i = 1,3
            cdum = 0
            do ikpt = 1,nkpt
              cfac = exp(img * 2*pi * 1 * kpt(1,ikpt) ) !; write(*,*) cfac
              cdum = cdum + imag ( log( 1 + dot_product( cfac * uwan(:,iwan,ikpt,ispin), matmul( dcprod(:,:,ikpt,i),
     &                                                   cfac * uwan(:,iwan,ikpt,ispin) ) ) * sqrt(vol) / 1000 ) ) * 1000
            enddo
            if(abs(imag(cdum))>1d-12) Error('Nonzero imaginary part in Wannier center.')
            fwancent(i,iwan,ispin) = cdum / nkpt
          enddo
          fwancent(:,iwan,ispin) = matmul( invert(lat), fwancent(:,iwan,ispin) )
          write(*,*) fwancent(:,iwan,ispin)
        enddo
      enddo
      deallocate(dcprod)
      stop
# endif

      deallocate ( mat2,decomp,fwancent,efrozen )
 10   deallocate ( lwan,lmwan,centwan,eulerwan )

      call getkey(inp,'SUBSET',charr,section='WANNIER',status=i)
      if     (i==1) then ; Error('List of indices missing after keyword SUBSET.')
      else if(i==2) then
        if(size(charr)>1) Error('Expected only one argument after keyword SUBSET.')
        call parser(parse_out,charr(1),'(1)','SUBSET')
        call str2iarr(lwan,parse_out(1),'SUBSET')
        if(any(lwan<1.or.lwan>nwan)) Error('Wannier indices out of range (keyword SUBSET).')
        nwan = size(lwan)
        do i = 1,nwan
          if(count(lwan==lwan(i))>1) Error('Double Wannier indices after keyword SUBSET.')
        enddo
        do ikpt = 1,nkpt2
          uwan(:,:nwan,ikpt,:) = uwan(:,lwan,ikpt,:)
        enddo
        call checkmem('uwan-realloc',16d0*size(uwan))
        call checkmem('uwan',-16d0*size(uwan)*2)
        call reallocate(uwan,wanbandf,nwan,nkpt2,nspin1)
        call checkmem('uwan',16d0*size(uwan))
        deallocate(lwan,charr)
        if(storeibz.and.wbloch) Warn('The Wannier subset might break symmetry. Usage of WBLOCH without STOREBZ critical!')
      endif

      call getkey(inp,'DISENTGL',ldum,section='WANNIER',default=.false.)
      if(ldum) call disentgl

      NoLoad( call prepare_wbloch )

      call getkey(inp,'IRREP',ldum,section='WANNIER',default=.false.)
      if(ldum) call irreps_wannier

      if(l_inter) then
        call getkey(inp,'INTERPOL', charr, section='WANNIER', status=i)
        ldum = .false.
        if(i==2) then
          n = size(charr)
          if(any(charr=='ENERGY')) then
            ldum = n>1
            if(ldum) ldum = charr(n-1)=='ENERGY'
            if(.not.ldum) Error('If ENERGY is specified after INTERPOL, it has to be the last but one argument.')
          endif
        endif
        if(ldum) then
          allocate(enehlp(wanbandi:wanbandf,nkpti,nspin1))
          call get_energies(enehlp,[wanbandi,wanbandf],nkpti,wanbandi,wanbandf,charr(size(charr)),.false.,0)
          call wannier_interpolation('0',(1d0,0d0)*enehlp,0)
          deallocate(enehlp)
        else
          call wannier_interpolation('0',(1d0,0d0)*ene(wanbandi:wanbandf,:nkpti,:),0)
        endif
        if(allocated(charr)) deallocate(charr)
      endif

      endSingle

# ifdef MPI
      call Mcast(nwan)     ; call Mcast(wanbandi)   ; call Mcast(wanbandf)
      call Mcast(nwanband) ; call Mcastl(irrep_wan) ; call Mcastl(irrep_wan1)     
      call Mcast(wscale)   ; call Mcastl(wancent)   ; call Mcast(wbloch)
      MnoR( Allocate_ ( uwan,(wanbandi:wanbandf,nwan,nkpt2,nspin1) ) )
      call Mcast(uwan)
#   ifndef LOAD
      MpiR( n = size(cmtu,4) )
      call Mcast(n)
      MnoR( Allocate_ ( cmtu,(maxlmindx,ncent,nwan,n,nspin2)   ) )
      MnoR( Allocate_ ( cpwu,(maxgpt,nwan,n,nspin2)            ) )
      call Mcast(cmtu)
      call Mcast(cpwu)
#   endif
# endif

      Rcall getkey(inp,'PLOT',ldum,section='WANNIER',default=.false.) ; Mpi( call Mcast(ldum) )      
      if(ldum) call wannier_plot

      contains

c --------------

      function dot_product_mt(cmt1,cmt2)
      implicit none
      complex_dp             :: dot_product_mt
      complex_dp, intent(in) :: cmt1(maxlmindx,ncent),cmt2(maxlmindx,ncent)
      integer                :: itype,ieq,ic,lm,l
      dot_product_mt = 0
      ic             = 0
      do itype = 1,ntype
        lm = sum ( [ (nindx(l,itype)*(2*l+1),l=0,lcut(itype)) ] )
        do ieq = 1,neq(itype)
          ic             = ic + 1
          dot_product_mt = dot_product_mt + dot_product(cmt1(:lm,ic),cmt2(:lm,ic))
        enddo
      enddo
      end function dot_product_mt

c --------------

c Disentangles band structure according to Wannier set:
c - Construction of H = PHP + (1-P)H(1-P) where P=SUM(m)|Wm><Wm| projects onto the Wannier set. (The off-diagonal blocks are zero.)
c - Diagonalization -> New single particle states (energies, wave functions, degs, U matrices)
      subroutine disentgl
      use global
      use wrapper
      use, intrinsic :: iso_fortran_env
      implicit none
      complex_dp  :: proj(nwanband,nwanband)
      complex_dp  :: hamil(nwanband,nwanband)
      MCOMPLEX_dp :: evec(nwanband,nwanband)
      real_dp     :: eval(nwanband),eord(maxband),rdum
      integer     :: ikpt,ispin,ic,i,j,k
      write(6,'(/A)') 'Bands are disentangled.'
      if(l_soc) Error('Spin-orbit coupling not implemented.')
      if(nwanband==nwan) Warn('All transformations unitary. No effect expected.')
      do ispin = 1,nspin
        do ikpt = 1,nkpt2
          if(kptp(ikpt)/=ikpt) cycle ! loop only over parent k points
          ! projector P = SUM(w) |w><w|
          proj = matmul ( uwan(:,:,ikpt,ispin) , transpose(conjg(uwan(:,:,ikpt,ispin))) )
          ! H - PH - HP + 2PHP
          eval  = ene(wanbandi:wanbandf,ikpt,ispin)
          hamil = 0
          do i = 1,nwanband
            hamil(i,i) = eval(i)
            do j = 1,nwanband
              hamil(i,j) = hamil(i,j) - proj(i,j) * (eval(i)+eval(j)) + 2 * sum(proj(i,:)*eval*proj(:,j))
            enddo
          enddo
          ! Check and diagonalize
          rdum = sum ( [ ((abs(hamil(i,j)-conjg(hamil(j,i))),j=1,i),i=1,nwanband) ] )
          if(rdum>1d-10) Error('Projected Hamiltonian not Hermitian.')
# ifdef INV
          rdum = sum ( [ ((abs(hamil(i,j)-hamil(j,i)),j=1,i),i=1,nwanband) ] )
          if(rdum>1d-10) Error('Projected Hamiltonian not symmetric.')
          call diagonalize(evec,eval,real(hamil))
# else
          call diagonalize(evec,eval,hamil)
# endif
          ! Change wave functions
          k = kindx(ikpt)
          do ic = 1,ncent
            do i = 1,size(cmt,1)
              cmt(i,ic,wanbandi:wanbandf,k,ispin) = matmul ( cmt(i,ic,wanbandi:wanbandf,k,ispin) , evec )
            enddo
          enddo
          do i = 1,ngpt(ikpt)
            cpw(i,wanbandi:wanbandf,k,ispin) = matmul ( cpw(i,wanbandi:wanbandf,k,ispin) , evec )
          enddo
          ! Transformation to equivalent k points (energies, wave functions, U matrices: U1=<phi1|w>=<phi1|phi0><phi0|w>=conjg(evec)*U0)
          do k = 1,nkpt2
            if(kptp(k)/=ikpt) cycle
            ! Energies
            ene(wanbandi:wanbandf,k,ispin) = eval
            ! Degeneracies (carbon copy of getinput, we loop over all bands for simplicity)
            i = 1
            do while(i<=nband(k,ispin))
              j = i
              do while(j<nband(k,ispin))
                if(abs(ene(j+1,k,ispin)-ene(j,k,ispin))>edeg) exit
                j = j + 1
              enddo
              deg(i,k,ispin)   = j
              if(i<j) deg(i+1:j,k,ispin) = i
              rdum             = sum(ene(i:j,k,ispin)) / (j-i+1)
              ene(i:j,k,ispin) = rdum
              i                = j + 1
            enddo
            ! Wave functions
            if(k/=ikpt) then
              if(.not.storeibz) then
# ifndef old_trafo
                do i = 1,nband(k,ispin)
                  call waveftrafo(cmt(:,:,:,k,ispin),cpw(:,:,k,ispin),symkpt(k),kptp(k),ispin)                  
                enddo
# else
                call waveftrafo2(cmt(:,:,:,k,ispin),cpw(:,:,k,ispin),kptp(k),symkpt(k),ispin)
# endif
              endif
            endif
            ! Wannier functions
            do i = 1,nwan
              if(symkpt(k)>nsymt) then
                uwan(:,i,k,ispin) = matmul ( uwan(:,i,k,ispin) , evec )
                call teststop('DISENTGL+time reversal')
              else
# ifdef INV
                uwan(:,i,k,ispin) = matmul ( uwan(:,i,k,ispin)/phase(k) , MCONJG(evec) ) / conjg(phase(k))
# else
                uwan(:,i,k,ispin) = matmul ( uwan(:,i,k,ispin)          , MCONJG(evec) )
# endif
              endif
            enddo
          enddo
          ! Check ordering
          eord(:nband(ikpt,ispin)) = ene(:nband(ikpt,ispin),ikpt,ispin)
          call rorder(eord,nband(ikpt,ispin))
          if(any(eord(:nband(ikpt,ispin))/=ene(:nband(ikpt,ispin),ikpt,ispin)))
     &      Error('Band ordering changed (not implemented).')
        enddo
      enddo
      end subroutine disentgl

c --------------

      subroutine write_wannier_win(ispin)
      implicit none
      integer, intent(in) :: ispin
      integer             :: iunit,iunit_s,ind,i,k,ios
      logical             :: lexist,lentry(4),found
      real_dp             :: minproj
      real_dp             :: efrozen1,efrozen0,r
      character(256)      :: line,line1
      minproj = dble(nwan)/nwanband*0.99d0 ! each frozen state must have a total projection into the Wannier basis of minproj or higher
      if(l_frozen) then
        if(efrozen(ispin)==huge(0d0)) then
          efrozen1 = minval(ene(wanbandi+nwan,:,ispin)) - 1d-6
          efrozen0 = efrozen1
          do i = wanbandf,wanbandi,-1
            found = .false.
            do k = 1,nkpti
              if(sum(abs(uwan(i,:,k,ispin))**2)<minproj) then
                efrozen1 = min( efrozen1 , ene(i,k,ispin) - 1d-6 )
                found    = .true.
              endif
            enddo
            if(.not.found) exit
          enddo
          if(found) Warn('Possibly bad choice of Wannier orbitals!')
          r = huge(0d0)
          do i = wanbandi,wanbandf
            do k = 1,nkpti
              r = min(r,sum(abs(uwan(i,:,k,ispin))**2))
            enddo
          enddo
          if(r<minproj*0.8d0) Warn('Some of the bands in the frozen window have low projection ('//Chr(r)//') into Wannier basis.')
          write(6,'(A)') 'Parameter dis_froz_max set to '//Chf(efrozen1*hartree,'F10.6')//' eV'
          if(efrozen0==efrozen1) then
            write(6,'(A)') 'This is the highest possible dis_froz_max with this number of Wannier functions.'
          else
            write(6,'(A)') 'Highest possible dis_froz_max with this number of Wannier functions: '//
     &                     Chf(efrozen0*hartree,'F10.6')//' eV'
          endif
        else
          efrozen1 = efrozen(ispin)
        endif
      endif
      lentry = .false.
      inquire(file='wannier.win',exist=lexist)
      if(lexist) then
        write(6,'(A)') 'File wannier.win exists and will be used.'
        iunit   = fopen('wannier.win',status='old',action='readwrite')
        iunit_s = fopen(status='scratch')
        do
          read(iunit,'(A)',iostat=ios) line ; if(ios/=0) exit
          line1 = line
          ind   = index(line1,'#') ; if(ind/=0) line1 = line1(:ind-1)
          ind   = index(line1,'!') ; if(ind/=0) line1 = line1(:ind-1)
          line1 = adjustl(line1)
          if(line1(:9)=='num_wann ') then
            if(lentry(1)) Error('Double keyword in wannier.win: num_wann.')
            lentry(1) = .true.
            read(line1(9:),*,iostat=ios) i
            if(ios/=0) Error('wannier.win: Error reading num_wann argument.')
            if(i/=nwan) then
              write(iunit_s,'(A)') 'num_wann '//Chr(nwan)
              write(6,'(A)') 'wannier.win: num_wann overwritten'
            else
              write(iunit_s,'(A)') trim(line)
            endif
          else if(line1(:10)=='num_bands ') then
            if(lentry(2)) Error('Double keyword in wannier.win: num_bands.')
            lentry(2) = .true.
            read(line1(10:),*,iostat=ios) i
            if(ios/=0) Error('wannier.win: Error reading num_bands argument.')
            if(i/=nwanband) then
              write(iunit_s,'(A)') 'num_bands '//Chr(nwanband)
              write(6,'(A)') 'wannier.win: num_bands overwritten'
            else
              write(iunit_s,'(A)') trim(line)
            endif
          else if(line1(:13)=='dis_froz_max ') then
            if(lentry(3)) Error('Double keyword in wannier.win: dis_froz_max.')
            lentry(3) = .true.
            read(line1(13:),*,iostat=ios) r
            if(ios/=0) Error('wannier.win: Error reading dis_froz_max argument.')
            if(l_frozen) then
              if(abs(r-efrozen1*hartree)>1d-8) then
                write(iunit_s,'(A)') 'dis_froz_max '//Chf(efrozen1*hartree,'F18.10')
                write(6,'(A)') 'wannier.win: dis_froz_max overwritten'
              else
                write(iunit_s,'(A)') trim(line)
              endif
            else
              write(6,'(A)') 'wannier.win: dis_froz_max removed'
            endif
          else if(line1(:9)=='num_iter ') then
            if(lentry(4)) Error('Double keyword in wannier.win: num_iter.')
            lentry(4) = .true.
            read(line1(9:),*,iostat=ios) i
            if(ios/=0) Error('wannier.win: Error reading num_iter argument.')
            if(l_wanmax) then
              if(i==0) then ; write(6,'(A)') 'wannier.win: num_iter removed'
              else          ; write(iunit_s,'(A)') trim(line)
              endif
            else
              if(i/=0) then
                write(iunit_s,'(A)') 'num_iter 0'
                write(6,'(A)') 'wannier.win: num_iter overwritten'
              else
                write(iunit_s,'(A)') trim(line)
              endif
            endif
          else
            write(iunit_s,'(A)') trim(line)
          endif
        enddo
        rewind(iunit)
        rewind(iunit_s)
        do ! copy scratch -> wannier.win
          read(iunit_s,'(A)',iostat=ios) line ; if(ios/=0) exit ; write(iunit,'(A)') trim(line)
        enddo
        call fclose(iunit_s)        
      else
        iunit = fopen('wannier.win',status='new')
        write(6,'(A)') 'Write default wannier.win.'
      endif
      if(.not.lentry(1)) then
        write(iunit,'(A)') 'num_wann '//Chr(nwan)
        if(lexist) write(6,'(A)') 'wannier.win: num_wann added'
      endif
      if(.not.lentry(2)) then
        write(iunit,'(A)') 'num_bands '//Chr(nwanband)
        if(lexist) write(6,'(A)') 'wannier.win: num_bands added'
      endif
      if(.not.lentry(3).and.l_frozen) then
        write(iunit,'(A)') 'dis_froz_max '//Chf(efrozen1*hartree,'F18.10')
        if(lexist) write(6,'(A)') 'wannier.win: dis_froz_max added'
      endif
      if(.not.lentry(4).and..not.l_wanmax) then
        write(iunit,'(A)') 'num_iter 0'
        if(lexist) write(6,'(A)') 'wannier.win: num_iter added'
      endif
      call fclose(iunit)
      end subroutine write_wannier_win

c --------------

      end

c --------------

      ! Returns nwan, lwan, lmwan, centwan, eulerwan interpreted from ORBITALS line [charr(:ncharr)]
      ! lwan(:)        : 0 - 3     : l-quantum number of full shell
      !                  otherwise : lwan = l + 100*multiplicity / l = l-quantum_number or 10,11,... for hybrids
      !                  examples  : 302 (t2g), 101 (py), 210 (sp hybrids), 502 (d, same as 2, but not treated as full shell!)
      ! lmwan(:)       : real harmonics index (see below), hybrids >= 100
      ! centwan(:)     : > 0 :  centwan is atom index
      !                  < 0 : -centwan is atom index, full atom type defined
      ! eulerwan(:3,:) : Euler angels for orbital rotation
c begin interface
      subroutine orbitals(lwan,lmwan,centwan,eulerwan,nwan,charr,ncharr)
      use global, only: ntype,neq
      use, intrinsic :: iso_fortran_env !inc
      implicit none
      real_dp,      intent(out), allocatable :: eulerwan(:,:)
      integer,      intent(out), allocatable :: lwan(:),lmwan(:),centwan(:)
      integer,      intent(out)              :: nwan
      integer,      intent(in)               :: ncharr
      character(*), intent(in)               :: charr(ncharr)
c end interface
      character(:), allocatable              :: label
      integer                                :: i,ic
      if(ncharr==0) Bug('Orbitals definition missing.')
      allocate(character(len(charr(1))) :: label)
      nwan = 0
      i    = 0
      ic   = 0
      do
        i = i + 1
        call get_orbital_label(label,ic,i,charr,ncharr) ; if(label==' ') exit
        call def_wannierorbital(label,ic)
      enddo
      if(nwan==0) Error('No orbitals defined.')
      deallocate(label)

      contains

      ! for orbitals: updates lwan, lmwan, centwan, eulerwan of current orbital definition in label/ic and counts up nwan
      subroutine def_wannierorbital(label,ic)
      use, intrinsic :: iso_fortran_env
      implicit none
      character(*), intent(in) :: label
      integer                  :: ic
      integer                  :: ind,ind1,nwan0,i,j,l,ios
      real_dp                  :: angl,ch2r
      logical                  :: define
      if(label==' ') return
      define = .false.
      nwan0  = nwan
 1    ind    = index(label,'/')
      if(define) eulerwan(:,nwan+1) = 0
      if(ind==0) then
        ind  = len_trim(label)
      else
        ind1 = index(label(ind+1:),'/')                         ; if(ind1==0) Error('Incorrect definition of Euler angles.')
        angl = ch2r(label(ind+1:ind+ind1-1),'ORBITALS/PROJECT') ; if(define) eulerwan(1,nwan+1) = angl        
        ind  = ind + ind1
        ind1 = index(label(ind+1:),'/')                         ; if(ind1==0) Error('Incorrect definition of Euler angles.')
        angl = ch2r(label(ind+1:ind+ind1-1),'ORBITALS/PROJECT') ; if(define) eulerwan(2,nwan+1) = angl        
        ind  = ind + ind1        
        angl = ch2r(label(ind+1:),          'ORBITALS/PROJECT') ; if(define) eulerwan(3,nwan+1) = angl        
        ind  = index(label,'/') - 1
      endif
      select case(label(:ind))
        case('s')    ; nwan = nwan + 1 ; l =   0 ; if(define) lmwan(nwan)        = 0
        case('S')    ; nwan = nwan + 1 ; l =   0 ; if(define) lmwan(nwan)        = 0
        case('p')    ; nwan = nwan + 3 ; l =   1 ; if(define) lmwan(nwan-2:nwan) =  [ 2,3,1 ]
        case('P')    ; nwan = nwan + 3 ; l =   1 ; if(define) lmwan(nwan-2:nwan) = -[ 1,2,3 ]
        case('d')    ; nwan = nwan + 5 ; l =   2 ; if(define) lmwan(nwan-4:nwan) =  [ 6,7,5,8,4]
        case('D')    ; nwan = nwan + 5 ; l =   2 ; if(define) lmwan(nwan-4:nwan) = -[ 4,5,6,7,8 ]
        case('f')    ; nwan = nwan + 7 ; l =   3 ; if(define) lmwan(nwan-6:nwan) =  [ 12,13,11,14,10,15,9 ]
        case('F')    ; nwan = nwan + 7 ; l =   3 ; if(define) lmwan(nwan-6:nwan) = -[ 9,10,11,12,13,14,15 ]
        case('g')    ; nwan = nwan + 9 ; l =   4 ; if(define) lmwan(nwan-8:nwan) =  [ 16,17,18,19,20,21,22,23,24 ]        
        case('G')    ; nwan = nwan + 9 ; l =   4 ; if(define) lmwan(nwan-8:nwan) = -[ 16,17,18,19,20,21,22,23,24 ]
        case('px')   ; nwan = nwan + 1 ; l = 101 ; if(define) lmwan(nwan)        = 3
        case('py')   ; nwan = nwan + 1 ; l = 101 ; if(define) lmwan(nwan)        = 1
        case('pz')   ; nwan = nwan + 1 ; l = 101 ; if(define) lmwan(nwan)        = 2
        case('dxy')  ; nwan = nwan + 1 ; l = 102 ; if(define) lmwan(nwan)        = 4
        case('dyz')  ; nwan = nwan + 1 ; l = 102 ; if(define) lmwan(nwan)        = 5
        case('dz2')  ; nwan = nwan + 1 ; l = 102 ; if(define) lmwan(nwan)        = 6
        case('dxz')  ; nwan = nwan + 1 ; l = 102 ; if(define) lmwan(nwan)        = 7
        case('dx2y2'); nwan = nwan + 1 ; l = 102 ; if(define) lmwan(nwan)        = 8
        case('t2g')  ; nwan = nwan + 3 ; l = 302 ; if(define) lmwan(nwan-2:nwan) = [4,5,7]
        case('eg')   ; nwan = nwan + 2 ; l = 202 ; if(define) lmwan(nwan-1:nwan) = [6,8]
        case('fz3')  ; nwan = nwan + 1 ; l = 103 ; if(define) lmwan(nwan)        = 12
        case('fxz2') ; nwan = nwan + 1 ; l = 103 ; if(define) lmwan(nwan)        = 13
        case('fyz2') ; nwan = nwan + 1 ; l = 103 ; if(define) lmwan(nwan)        = 11
        case('fzxy') ; nwan = nwan + 1 ; l = 103 ; if(define) lmwan(nwan)        = 14
        case('fxyz') ; nwan = nwan + 1 ; l = 103 ; if(define) lmwan(nwan)        = 10
        case('fxxy') ; nwan = nwan + 1 ; l = 103 ; if(define) lmwan(nwan)        = 15
        case('fyxy') ; nwan = nwan + 1 ; l = 103 ; if(define) lmwan(nwan)        = 9
        case('sp')   ; nwan = nwan + 2 ; l = 210 ; if(define) lmwan(nwan-1:nwan) = [ 101,102 ]
        case('sp2')  ; nwan = nwan + 3 ; l = 311 ; if(define) lmwan(nwan-2:nwan) = [ 103,104,105 ]
        case('sp3')  ; nwan = nwan + 4 ; l = 412 ; if(define) lmwan(nwan-3:nwan) = [ 106,107,108,109 ]
        case('sp3d') ; nwan = nwan + 5 ; l = 513 ; if(define) lmwan(nwan-4:nwan) = [ 110,111,112,113,114 ]
        case('sp3d2'); nwan = nwan + 6 ; l = 614 ; if(define) lmwan(nwan-5:nwan) = [ 115,116,117,118,119,120 ]
        case default ; Error('Missing or wrong orbital label: '//trim(label(:ind)))
      end select
      if(define) then
        centwan(nwan0+1:nwan) = ic
        lwan   (nwan0+1:nwan) = l
        do i = nwan0+2,nwan
          eulerwan(:,i) = eulerwan(:,nwan0+1)
        enddo
      else
        define = .true. ; call realloc_wan(nwan)        
        nwan   = nwan0
        goto 1
      endif
      end subroutine def_wannierorbital

      subroutine realloc_wan(nwan)
      use util, only: reallocate
      implicit none
      integer, intent(in) :: nwan
      if(allocated(lwan)) then
        if(nwan>size(lwan)) then
          call reallocate(lwan,nwan)
          call reallocate(lmwan,nwan)
          call reallocate(centwan,nwan)
          call reallocate(eulerwan,3,nwan)
        endif
      else
        allocate(lwan(nwan),lmwan(nwan),centwan(nwan),eulerwan(3,nwan))
      endif
      end subroutine realloc_wan

      end

c --------------

c     Returns ilabel-th orbital label in "label" and atom index in "ic" from orbitals definition charr(:ncharr).
      subroutine get_orbital_label(label,ic,ilabel,charr,ncharr)
      use global, only: neq,ncent,atype
      implicit none
      integer,      intent(in)  :: ncharr,ilabel
      character(*), intent(in)  :: charr(ncharr)
      character(*), intent(out) :: label
      integer,      intent(out) :: ic
      integer                   :: itype
      integer                   :: i,j,l,p1,p2
      if(ncharr==0) return
      j  = 0
      ic = 0
      do i = 1,ncharr
        l  = len_trim(charr(i))
        ic = ic + 1
        if(ic>ncent)                                           Error('Atom count error in ORBITALS or PROJECT.')
        if(all( charr(i)(1:1)//charr(i)(l:l) /= ['[]','()'] )) Error('Expected (...) or [...] in line ORBITALS or PROJECT.')
        if(charr(i)(1:1)=='[') then
          itype = atype(ic)
          if(ic/=sum(neq(:itype-1))+1) Error('Expected (...) instead of [...].')
        else
          itype = 0
        endif
 1      p2 = 0
        do while(p2<l-1)
          p1    = p2 + 2
          p2    = index(charr(i)(p1:),',') + p1 - 2 ; if(p2==p1-2) p2 = l - 1
          label = charr(i)(p1:p2)                   ; if(label==' ') cycle
          j     = j + 1                             ; if(j==ilabel) then ; if(itype/=0) ic = -ic ; return ; endif
        enddo
        if(itype/=0.and.ic<sum(neq(:itype))) then
          ic = ic + 1
          goto 1
        endif
      enddo
      label = ' '
      end

c --------------

c     Adds harmonics expansion coefficients of l to proj: L(r)=SUM(lm) proj(lm * Y_lm(r)
c     according to lmwan and eulerwan (in deg) from subroutine orbitals.
      subroutine add_proj(proj,l,lmwan,eulerwan)
      use global, only: pi
      use, intrinsic :: iso_fortran_env
      implicit none
      integer,    intent(in)    :: l,lmwan
      real_dp,    intent(in)    :: eulerwan(3)
      complex_dp, intent(inout) :: proj(-l:l)
      complex_dp                :: a(-l:l),b(-l:l),wignerd(-l:l,-l:l,0:l)
      real_dp                   :: rot(3,3)
      integer                   :: lm,m
      real_dp,    parameter     :: sq2=1/sqrt(2d0),sq3=1/sqrt(3d0),sq6=1/sqrt(6d0)
      real_dp,    parameter     :: hybrid(0:8,20) = reshape ( [
c         s      py      pz      px    dxy   dyz    dz2     dxz    dx2y2
     &   sq2 ,   0d0 ,   0d0 ,   sq2 , 0d0 , 0d0 ,  0d0   , 0d0 ,   0d0 , ! sp-1
     &   sq2 ,   0d0 ,   0d0 ,  -sq2 , 0d0 , 0d0 ,  0d0   , 0d0 ,   0d0 , ! sp-2
     &   sq3 ,   sq2 ,   0d0 ,  -sq6 , 0d0 , 0d0 ,  0d0   , 0d0 ,   0d0 , ! sp2-1
     &   sq3 ,  -sq2 ,   0d0 ,  -sq6 , 0d0 , 0d0 ,  0d0   , 0d0 ,   0d0 , ! sp2-2
     &   sq3 ,   0d0 ,   0d0 , 2*sq6 , 0d0 , 0d0 ,  0d0   , 0d0 ,   0d0 , ! sp2-3
     &  .5d0 ,  .5d0 ,  .5d0 ,  .5d0 , 0d0 , 0d0 ,  0d0   , 0d0 ,   0d0 , ! sp3-1
     &  .5d0 , -.5d0 , -.5d0 ,  .5d0 , 0d0 , 0d0 ,  0d0   , 0d0 ,   0d0 , ! sp3-2
     &  .5d0 ,  .5d0 , -.5d0 , -.5d0 , 0d0 , 0d0 ,  0d0   , 0d0 ,   0d0 , ! sp3-3
     &  .5d0 , -.5d0 ,  .5d0 , -.5d0 , 0d0 , 0d0 ,  0d0   , 0d0 ,   0d0 , ! sp3-4
     &   sq3 ,   sq2 ,   0d0 ,  -sq6 , 0d0 , 0d0 ,  0d0   , 0d0 ,   0d0 , ! sp3d-1
     &   sq3 ,  -sq2 ,   0d0 ,  -sq6 , 0d0 , 0d0 ,  0d0   , 0d0 ,   0d0 , ! sp3d-2
     &   sq3 ,   0d0 ,   0d0 , 2*sq6 , 0d0 , 0d0 ,  0d0   , 0d0 ,   0d0 , ! sp3d-3
     &   0d0 ,   0d0 ,   sq2 ,   0d0 , 0d0 , 0d0 ,  sq2   , 0d0 ,   0d0 , ! sp3d-4
     &   0d0 ,   0d0 ,  -sq2 ,   0d0 , 0d0 , 0d0 ,  sq2   , 0d0 ,   0d0 , ! sp3d-5
     &   sq6 ,   0d0 ,   0d0 ,  -sq2 , 0d0 , 0d0 , -sq3/2 , 0d0 ,  .5d0 , ! sp3d2-1
     &   sq6 ,   0d0 ,   0d0 ,   sq2 , 0d0 , 0d0 , -sq3/2 , 0d0 ,  .5d0 , ! sp3d2-2
     &   sq6 ,  -sq2 ,   0d0 ,   0d0 , 0d0 , 0d0 , -sq3/2 , 0d0 , -.5d0 , ! sp3d2-3
     &   sq6 ,   sq2 ,   0d0 ,   0d0 , 0d0 , 0d0 , -sq3/2 , 0d0 , -.5d0 , ! sp3d2-4
     &   sq6 ,   0d0 ,  -sq2 ,   0d0 , 0d0 , 0d0 ,  sq3   , 0d0 ,   0d0 , ! sp3d2-5
     &   sq6 ,   0d0 ,   sq2 ,   0d0 , 0d0 , 0d0 ,  sq3   , 0d0 ,   0d0   ! sp3d2-6
     &   ] , [ 9,20] )
      lm = abs(lmwan)
      if(lm>100) then
        if(l>2) return
        b = hybrid(l**2:(l+1)**2-1,lmwan-100)
        call harmonics_r2c(a,b,l)
      else if(lm<l**2.or.lm>=(l+1)**2) then
        return
      else
        m = lm - l*(l+1)
        if(lmwan<=0) then
          a    = 0
          a(m) = 1d0
        else
          b    = 0
          b(m) = 1d0
          call harmonics_r2c(a,b,l)
        endif
      endif
      if(any(eulerwan/=0)) then
        call rotation_euler(rot,eulerwan*pi/180)
        call dwigner(wignerd,rot,l)
        a = matmul(wignerd(-l:l,-l:l,l),a)
      endif
      proj = proj + a
      end

c --------------

      ! Premultiply cmt/cpw with uwan -> cmtu, cpwu
      subroutine prepare_wbloch
      use global
      Mpi( use Mwrapper )
      Load( use readwrite ) LoadC( only: read_wavef3 )
      use, intrinsic :: iso_fortran_env
      implicit none
      integer :: ispin,ikpt,ic,s,ikpt1,nkpt1,kpt1(nkpt2)
      kpt1 = [(ikpt,ikpt=1,nkpt2)]
      if(storeibz) then ; nkpt1 = nkpti + nkpti2 ; kpt1(nkpti+1:nkpt1) = nkpt + [(ikpt,ikpt=1,nkpti2)]
      else              ; nkpt1 = nkpt2
      endif
      if(associated(cmtu)) Bug('array cmtu already allocated.')
      Allocate_ ( cmtu,(maxlmindx,ncent,nwan,nkpt1,nspin2) )
      Allocate_ ( cpwu,(maxgpt,nwan,nkpt1,nspin2)          )
# ifndef LOAD
      write(6,'(A,F6.1," MB") ') 'Size of array cmtu:',size(cmtu)*16d0 / megabyte
      write(6,'(A,F6.1," MB"/)') 'Size of array cpwu:',size(cpwu)*16d0 / megabyte
# endif
      do ispin = 1,nspin2
        Load( call read_wavef3(wanbandi,wanbandf,kpt1,nkpt1,ispin) )
        s = min(ispin,nspin1)
        do ikpt = 1,nkpt1 ; Mcycle(ikpt)
          ikpt1 = kpt1(ikpt)
          ! Premultiply coefficients -> cmtu and cpwu
          do ic = 1,ncent
            cmtu(:,ic,:,ikpt,ispin) = matmul ( cmt(:,ic,wanbandi:wanbandf,ikpt, ifLoad(1,ispin) ) , uwan(:,:,ikpt1,s) )
          enddo
          cpwu(:,:,ikpt,ispin)      = matmul ( cpw(:,   wanbandi:wanbandf,ikpt, ifLoad(1,ispin) ) , uwan(:,:,ikpt1,s) )
        enddo
        Load( deallocate(cmt,cpw) )
        Mdistr(cmtu(:,:,:,ikpt,ispin),ikpt,1,nkpt1)
        Mdistr(cpwu(:,:,  ikpt,ispin),ikpt,1,nkpt1)
      enddo
      end
      
c --------------

      subroutine wannier_plot
      use global
      use key
      use util, only: chr
      use file
      use wrapper, only: identity
# ifdef MPI      
      use Mwrapper, only: Mcast,Msum,node_allocate
# endif
      implicit none
      complex_dp,   allocatable :: wavef(:,:)
      complex_dp,   pointer_cnt :: wanorb(:,:,:,:)
      real_dp,      allocatable :: rx(:,:)
      character(:), allocatable :: lines(:)
      character(80)             :: parse_out(2,3)
      complex_dp                :: phs
      real_dp                   :: origin(3)
      real_dp                   :: span(3,3)
      real_dp                   :: rdum,sum2
      integer                   :: iunit
      integer                   :: iwan,iband,ikpt,ispin
      integer                   :: nx,ny,nz,ix,iy,iz
      integer                   :: itype,ieq,ic
      integer                   :: ch2i
      real                      :: cputime
      real_dp                   :: ch2r
# ifdef MPI
      type(c_ptr)               :: ptr
      integer                   :: win_wanorb,Merr
# endif      
      Load( Error('Not implemented for LOAD. Inform the developers.') )
      Rbegin
      write(6,'(A'NoA) 'Calculate 3D plots of Wannier orbitals... '
      call cpu_time(cputime)
      call getkey(inp,'PLOT',lines,section='WANNIER')
      origin = 0
      span   = lat
      nx     = 20
      ny     = 20
      nz     = 20
      if(size(lines)>0) then
        ix = 1
        if     (lines(1)=='SC')   then ; do iy = 1,3 ; span(:,iy) = nkpt3(iy) * span(:,iy) ; enddo
        else if(lines(1)=='CART') then ; span = identity(3)
        else if(lines(1)/='UC')   then ; ix = 0
        endif
        if(size(lines)>ix) then
          if(index(lines(ix+1),'{')/=0) then
            ix = ix + 1
            call parser(parse_out,lines(ix),'{1:2,3:4,5:6}','PLOT')
            do iy = 1,3
              rdum       = ch2r(parse_out(1,iy),'PLOT')
              origin     = origin + span(:,iy) * rdum
            enddo
            do iy = 1,3
              rdum       = ch2r(parse_out(2,iy),'PLOT') - ch2r(parse_out(1,iy),'PLOT')
              span(:,iy) = span(:,iy) * rdum
            enddo
          endif
        endif
        if(size(lines)>ix) then
          if(size(lines)<ix+3) Error('Three numbers of data points nx,ny,nz expected after PLOT.')
          ix = ix + 1 ; nx = ch2i(lines(ix),'PLOT')
          ix = ix + 1 ; ny = ch2i(lines(ix),'PLOT')
          ix = ix + 1 ; nz = ch2i(lines(ix),'PLOT')          
        endif
      endif
      Rend
# ifdef MPI
      call Mcast(span) ; call Mcast(origin) ; call Mcast(nx) ; call Mcast(ny) ; call Mcast(nz)
# endif
      allocate(rx(3,nx),wavef(wanbandi:wanbandf,nx))
      Nallocate(wanorb,(S_ nx,ny,nz,nwan S_))
      ifO wanorb = 0
      ic = -1
      do ispin = 1,nspin1
        do iz = 1,nz
        do iy = 1,ny ; McycleP(ic)
        do ix = 1,nx
          rx(:,ix) = origin + span(:,1) * (ix-1)/max(nx-1,1) + span(:,2) * (iy-1)/max(ny-1,1) + span(:,3) * (iz-1)/max(nz-1,1)
          rx(:,ix) = matmul ( rx(:,ix) , rlat ) / (2*pi)
        enddo
        do ikpt = 1,nkpt
          if(nwanband<=3) then ! faster for small nwanband due to memory management I guess; checked only for 2D NiI2!
            do iband = wanbandi,wanbandf
              call wavefunction_r(wavef(iband,:),rx,nx,iband,ikpt,ispin)
            enddo
          else
            call wavefunction_r1(wavef,rx,nx,wanbandi,wanbandf,ikpt,ispin)
          endif
          do iwan = 1,nwan
            do ix = 1,nx
              wanorb(ix,iy,iz,iwan) = wanorb(ix,iy,iz,iwan) + sum ( uwan(:,iwan,ikpt,ispin) * wavef(:,ix) )
            enddo
          enddo
        enddo
        enddo ! iy
        enddo ! iz
        Nfence(wanorb)
        MpiO( call Msum(wanorb,comm=Ocomm) )
        Rbegin
        wanorb = wanorb / nkpt
        do iwan = 1,nwan
          sum2 = sum(abs(wanorb(:,:,:,iwan))**2)              
          phs  = sqrt( sum(wanorb(:,:,:,iwan)**2) )
          phs  = phs / abs(phs) ! average phase
          rdum = sum(real( wanorb(:,:,:,iwan) / phs )**2) / sum2
          ix   = nint(rdum*100)
          if(nspin1==2) then ; iunit = fopen('wannier'//Chr(iwan)//'s'//Chr(ispin)//'.xsf',status='unknown')
          else               ; iunit = fopen('wannier'//Chr(iwan)//'.xsf',status='unknown')
          endif
          write(iunit,'(A)') 'CRYSTAL','PRIMVEC'
          write(iunit,'(3F13.7)') lat * angstrom
          write(iunit,'(A)') 'PRIMCOORD'
          write(iunit,'(2I4)') ncent,1
          ic = 0
          do itype = 1,ntype
            do ieq = 1,neq(itype)
              ic = ic + 1
              write(iunit,'(I4,3F13.7)') ztype(itype), matmul(lat,cent(:,ic)) * angstrom
            enddo
          enddo
          write(iunit,'(/A)') 'BEGIN_BLOCK_DATAGRID_3D'
          write(iunit,'(A)') '  Wannier orbital'
          if(ix<100) Warn('Wannier orbital '//Chr(iwan)//': phase quality '//Chr(ix)//'%')
          if(ix>70) then
            wanorb(:,:,:,iwan) = wanorb(:,:,:,iwan) / phs ! the real part has the component along phs 
            write(iunit,'(A'NoA) '  BEGIN_DATAGRID_3D_'
          else
            wanorb(:,:,:,iwan) = abs(wanorb(:,:,:,iwan))
            write(iunit,'(A'NoA) '  BEGIN_DATAGRID_3D_ABS'
            if(ix<60) then
              Info('Absolute Wannier orbital plotted. Phase quality '//Chr(ix)//' < 60%. Phase information would be meaningless.')
            else
              Info('Absolute Wannier orbital plotted. Phase quality '//Chr(ix)//' < 70% loo low.')
            endif            
          endif
          write(iunit,'(A'NoA) 'orbital='//Chr(iwan)
          if(nspin1==2) then ; write(iunit,'(A)') '_spin='//Chr(ispin)
          else               ; write(iunit,*)
          endif
          write(iunit,'(3I6)') nx,ny,nz
          write(iunit,'(3F13.7)') origin * angstrom
          write(iunit,'(3F13.7)') span   * angstrom
          do iz = 1,nz
          do iy = 1,ny ; write(iunit,'(A'NoA) '    ' 
          do ix = 1,nx
            write(iunit,'(A'NoA) Chf(real(wanorb(ix,iy,iz,iwan)),'F17.7')//' '                
          enddo ; write(iunit,*)
          enddo ; if(iz<nz) write(iunit,*)
          enddo
          write(iunit,'(A)') '  END_DATAGRID_3D'
          Rwrite(iunit,'(A)') 'END_BLOCK_DATAGRID_3D'
          Rcall fclose(iunit)
        enddo        
        Rend
      enddo
      Rcall cpu_done(cputime)
      if(nwan==1) then ; Rwrite(6,'(A)') 'Wannier plot data written to wannier1.xsf'
      else             ; Rwrite(6,'(A)') 'Wannier plot data written to wannier1.xsf, wannier2.xsf, ...'
      endif
      deallocate(rx,wavef)
      Ndeallocate(wanorb)
      end

c --------------

c     Wannier interpolation of a diagonal Hamiltonian matrix (elements given in diag).
c
c     Real-space Hamiltonian is written to "hamiltonian"//suffx.
c     Interpolated bands are written to "bands"//suffx.
      subroutine wannier_interpolation(suffx,diag,spin)
      use global
      use util, only: chr,average_deg
      use file
      use key
      use wrapper
      use, intrinsic :: iso_fortran_env
      implicit none
      character(*), intent(in)  :: suffx
      integer,      intent(in)  :: spin
      complex_dp,   intent(in)  :: diag(nwanband,nkpti,*)
      complex_dp                :: hamiltonian(nwan,nwan),eigc(nwan),eigv(nwan,nwan)
      complex_dp,   allocatable :: hamiltonianr(:,:,:),cmtr(:,:,:),cmtk(:,:,:,:),bst(:,:,:)
      character(256)            :: qfile
      character(:), allocatable :: lines(:)
      real_dp                   :: eig(nwan),qvec0(3)
      real_dp,      allocatable :: qpath(:),qvec(:,:),proj(:,:,:),project(:,:,:,:),wproject(:,:,:,:)
      integer,      allocatable :: lpt(:,:),ldeg(:,:),pnt(:)
      integer                   :: nlpt,ndiff
      integer,      allocatable :: lpt1(:,:),ldeg1(:,:) ! only for bandinfo
      integer                   :: nlpt1,ndiff1         ! only for bandinfo
      integer                   :: ispin,i,iwan,nq,iq,ikpt
      integer                   :: ic
      integer                   :: iunit,iunitr
      logical                   :: lcompl
# include "interface/getqlist.inc"
# include "interface/wannier_spatial.inc"
# include "interface/band_info.inc"

      if(spin==0) then
        lcompl = any(imag(diag(:,:,:nspin1))/=0)
      else
        lcompl = any(imag(diag(:,:,1))/=0)
      endif
      
      call getkey(inp,'INTERPOL', lines, section='WANNIER', status=i)
      qfile = ' '
      if(i==2) then ; if(lines(1)/='ENERGY') qfile = lines(1) ; endif      
      if(allocated(lines)) deallocate(lines)
      call getqlist(qvec,nq,100,qfile)
      if(nq==0) Error('No list of q vectors provided (KPTPATH or file) for Wannier interpolation.')

      if(spin==0) then ; allocate ( bst(nq,nwan,nspin1),    wproject(nspin1,   nwan,nq,nwan) )
      else             ; allocate ( bst(nq,nwan,spin:spin), wproject(spin:spin,nwan,nq,nwan) )
      endif

      iunit  = fopen('bands'      //suffx,status='unknown')
      iunitr = fopen('hamiltonian'//suffx,status='unknown')

      allocate ( qpath(0:nq) )
      
      if(bandinfo) then
        call getkey(inp,'PROJECT',lines,status=i)
        if(i==0)      then ; Bug('PROJECT not found.')
        else if(i==1) then ; allocate(character(80) :: lines(ntype)) ; lines = '[s,p,d,f,g]'
        endif
        allocate(cmtk(maxlmindx,ncent,nwan,nspin3))
      endif

      call bands_header(iunit,'# Wannier-interpolated band structure',qpath,qvec,nq)

      if(lcompl) allocate(pnt(nwan))

      do ispin = 1,nspin1 ; if(spin/=0.and.spin/=ispin) cycle
        if(spin==0) then ; call wannier_spatial(hamiltonianr,diag(:,:,ispin),lpt,nlpt,ldeg,ndiff,ispin,0)
        else             ; call wannier_spatial(hamiltonianr,diag,           lpt,nlpt,ldeg,ndiff,ispin,0)
        endif
        if(bandinfo) call wannier_spatial(cmtr,[(0d0,0d0)],lpt1,nlpt1,ldeg1,ndiff1,ispin,2)
        write(iunitr,'(A)')    'Real-space Wannier Hamiltonian  < W(n0) | H | W(n''R) >'
        write(iunitr,'(A,I3)') 'Number of W functions:',nwan
        write(iunitr,'(A,I5)') 'Number of R vectors:',nlpt
        if(nspin1==2) then
          if(ispin==1) write(iunitr,'(/A)') 'Spin up'
          if(ispin==2) write(iunitr,'(/A)') 'Spin down'
        endif
        do i = 1,nlpt
          write(iunitr,*)
          write(iunitr,'(4I3)')                         lpt(:,i)
          write(iunitr,'('//chr(nwan)//'(2F15.10,3X))') transpose(hamiltonianr(i,:,:))
        enddo
        do iq = 1,nq
          call wannier_interpolate(hamiltonian,hamiltonianr,nwan**2,lpt,nlpt,qvec(:,iq))
          if(lcompl) then
            call diagonalize_gen(eigv,eigc,hamiltonian)
            call rorderpf(pnt,real(eigc),nwan) ! LAPACK routine does not order the eigenvalues
            bst(iq,:,ispin) = eigc(pnt)
            eigv            = eigv(:,pnt)
            eig             = real(eigc)
          else
            call diagonalize    (eigv,eig ,hamiltonian)
            bst(iq,:,ispin) = eig
          endif
          do iwan = 1,nwan
            wproject(ispin,iwan,iq,:) = average_deg( abs(eigv(iwan,:))**2 , eig )
          enddo
          if(bandinfo) then
            call wannier_interpolate(cmtk,cmtr,maxlmindx*ncent*nwan*nspin3,lpt1,nlpt1,qvec(:,iq))
            do ic = 1,ncent
              cmtk(:,ic,:,1) = matmul(cmtk(:,ic,:,1),eigv)
              if(l_soc) cmtk(:,ic,:,2) = matmul(cmtk(:,ic,:,2),eigv)
            enddo
            call band_info(lines,cmtk,spin=ispin,result=proj)
            if(.not.allocated(project)) then
              if(spin==0) then ; allocate(project(nspin1,   size(proj,1),nq,nwan))
              else             ; allocate(project(spin:spin,size(proj,1),nq,nwan))
              endif
              project = 0
            endif
            do i = 1,size(proj,1)
              project(ispin,i,iq,:) = average_deg( proj(i,1,:) , eig )
            enddo
            deallocate(proj)
          endif
        enddo
        deallocate ( lpt,hamiltonianr,ldeg )
        if(bandinfo) deallocate(lpt1,cmtr,ldeg1)
      enddo

      if(lcompl) deallocate(pnt)

      bst = bst * escale ! convert Hartree to energy unit (global.f)

      if(spin==0) then ; ispin = nspin1
      else             ; ispin = 1
      endif
      if(bandinfo) then
        write(iunit,'(A,9X,A,8X,A'NoA) '#','     ','       '
        if(ispin==2) write(iunit,'(8X,A'NoA) '       '
        call write_orbitals_header(iunit,lines,size(lines),ispin,10,1) ; write(iunit,*)
      endif
      if(lcompl) then
        write(iunit,'(A,9X,A,8X,A,8X,A'NoA) '#','kpath','real(e)','imag(e)'
        if(ispin==2) write(iunit,'(16X,A'NoA) '              '
      else
        write(iunit,'(A,9X,A,8X,A'NoA)      '#','kpath',' energy'
        if(ispin==2) write(iunit,'(8X,A'NoA) '       '
      endif
      if(bandinfo) call write_orbitals_header(iunit,lines,size(lines),ispin,10,2)
      write(iunit,'(A)') '    | Wannier projections ...'

      if(bandinfo) deallocate(cmtk,lines)

      ispin = size(bst,3)
      do iwan = 1,nwan
        do iq = 1,nq
          if(lcompl) then ; write(iunit,'('//chr(2*ispin+1)//'F15.9'NoA) qpath(iq),     bst(iq,iwan,:)
          else            ; write(iunit,'('//chr(  ispin+1)//'F15.9'NoA) qpath(iq),real(bst(iq,iwan,:))
          endif
          if(bandinfo) write(iunit,'('//chr(ispin*size(project,2))//'F10.4'NoA) project(:,:,iq,iwan)
          write(iunit,'('//chr(ispin*size(wproject,2))//'F10.4'NoA) wproject(:,:,iq,iwan)
          write(iunit,*)
        enddo
        write(iunit,*)
        write(iunit,*)
      enddo
      deallocate(bst,qpath,qvec)
      deallocate(wproject)
      if(bandinfo) deallocate(project)
      call fclose(iunit)
      call fclose(iunitr)

      write(6,'(A)') 'Band structure written to bands'//suffx

      end

c --------------

c     Calculates spatial Wannier representation of matk
c     Returns in matr spatial representation of matk(i,j) = <phi_ki|A|phi_kj>
c                                            or cmt(i,k) (for mode=2)
c
c     mode = 0: matk contains diagonal elements, e.g., (1d0,0d0)*ene(wanbandi:wanbandf,:nkpti,ispin)
c     mode = 1: matk contains full matrices,     e.g., (1d0,0d0)*matk(nwanband,nwanband,:nkpti)
c     mode = 2: MT coefficients are transformed:       cmt(:,:,wanbandi:wandbandf,:nkpti,ispin)
c               (matk is not referenced in this case)
c
c     WS_wancent:
c     If defined, Wannier centers are taken into account in the construction of the Wigner-Seitz cell.
c     Improves symmetry invariance of interpolated matrix. (Same as "use_ws_distance" of Wannier90.)
c     (Not applied for mode=2, as it only makes sense for interpolations of the type <wan|A|wan>.)
# define WS_wancent
c begin interface
      subroutine wannier_spatial(matr,matk,lpt,nlpt,ldeg,ndiff,ispin,mode)
      use global
      use wrapper
      use, intrinsic :: iso_fortran_env !inc
      implicit none
      integer,                 intent(inout) :: nlpt,ndiff
      integer,    allocatable, intent(inout) :: lpt(:,:),ldeg(:,:)
      integer,                 intent(in)    :: ispin,mode
      complex_dp,              intent(in)    :: matk(*)
      complex_dp, allocatable, intent(inout) :: matr(:,:,:)
c end interface
      complex_dp, allocatable                :: carr(:),cphase(:),cmt1(:,:,:,:)
      complex_dp                             :: cdum
      real_dp                                :: corn(3,8),vec(3),mat(3,3),rdiff(3,0:nwan*(nwan-1)),rdum
      real_dp,    allocatable                :: rdist(:)
      integer,    allocatable                :: neigh(:,:)
      integer                                :: idiff,diff(nwan,nwan)
      integer                                :: i1,i2,i3
      integer                                :: ilpt,ikpt,nneigh,ib
      integer                                :: i,j,iwan,iwan1
      integer                                :: itype,ic,l,lm,lm1,s
      logical                                :: first,compl
      real                                   :: time

      call cpu_time(time)

      if(mode<0.or.mode>2) Bug('Unknown mode.')

      rdiff(:,0) = 0
      ndiff      = 0
# ifdef WS_wancent
      if(mode<2) then
        i1 = 0
        do i = 1,nwan ; diff(i,i) = 0
          jloop: do j = 1,nwan ; if(i==j) cycle
            vec = wancent(:,j,ispin) - wancent(:,i,ispin)
            do i2 = 0,i1
              if(all(abs(vec-rdiff(:,i2))<1d-6)) then
                diff(i,j) = i2
                cycle jloop
              endif
            enddo
            i1          = i1 + 1
            diff(i,j)   = i1
            rdiff(:,i1) = vec
          enddo jloop
        enddo
        ndiff = i1
      endif
# endif

      if(.not.allocated(lpt)) then

      ! Define supercell neighbors (to construct Wigner-Seitz borders)
      ! (a) try to produce a more compact supercell by backfolding
      call backfold(mat,[(lat(:,i)*nkpt3(i),i=1,3)])
      ! (b) determine corners of "primitive WS cell" defined only by the neighbors 100 -100 010 etc. (wrt mat)
      do i = 1,3
        mat(:,i) = mat(:,i) / 2
        mat(:,i) = mat(:,i) / sum(mat(:,i)**2)
      enddo
      mat = transpose(mat)
      i   = 0
      do i1 = -1,1,2 ; do i2 = -1,1,2 ; do i3 = -1,1,2
        i = i + 1
        call solve(corn(:,i),mat,[i1,i2,i3]*1d0)
      enddo ; enddo ; enddo
      ! (c) add all neighbors whose planes cut the "primitive WS cell" and therefore could be relevant
      first  = .true.
 1    compl  = .false.
      nneigh = 6
      ib     = 1 ! extra shells
      i      = 1
      do while(ib>=0)
        i     = i + 1
        compl = .true.
        do i1 = -i,i
          do i2 = -i+abs(i1),i-abs(i1)
            i3  = i - abs(i1) - abs(i2)
            vec = matmul(lat,[i1,i2,i3]*nkpt3) / 2
            vec = vec / (sum(vec**2))
            do j = 1,8
              if(dot_product(vec,corn(:,j))>1-1d-10) then
                compl  = .false.
                nneigh = nneigh + 1 ; if(.not.first) neigh(:,nneigh) = [i1,i2,i3] * nkpt3
                if(i3/=0) then
                  nneigh = nneigh + 1 ; if(.not.first) neigh(:,nneigh) = [-i1,-i2,-i3] * nkpt3
                endif
                exit
              endif
            enddo
          enddo
        enddo
        if(compl) then ; ib = ib - 1
        else           ; ib = 1
        endif
      enddo
      if(first) then
        allocate ( neigh(3,nneigh) )
        neigh(:,:6) = 0
        do i = 1,3
          neigh(i,2*i-1:2*i) = [ nkpt3(i),-nkpt3(i) ] ! 100 -100 010 etc. always included
        enddo
        first = .false.
        goto 1
      endif

      ! Determine lattice points that lie in Wigner-Seitz cell (of supercell) (->lpt)
      allocate ( rdist(0:nneigh) )
      rdist = 0
      first = .true.
 2    nlpt  = 0
      compl = .false.
      i     = -1
      do while(.not.compl)
        i     = i + 1
        compl = .true.
        do i1 = -i,i
          do i2 = -i+abs(i1),i-abs(i1)
            i3 = i-abs(i1)-abs(i2)
            do
              do idiff = 0,ndiff
                rdist(0)  =    sqrt(sum(matmul(lat,[i1,i2,i3] + rdiff(:,idiff)             )**2))
                rdist(1:) = [ (sqrt(sum(matmul(lat,[i1,i2,i3] + rdiff(:,idiff) - neigh(:,j))**2)),j=1,nneigh) ]
                if(all(rdist(0)<rdist(1:)+1d-6)) then
                  compl = .false.
                  nlpt  = nlpt + 1
                  if(.not.first) lpt(:,nlpt) = [ i1,i2,i3 ]
                  exit
                endif
              enddo
              if(i3<=0) exit
              i3 = -i3
            enddo
          enddo
        enddo
      enddo
      if(first) then
        allocate( lpt(3,nlpt),ldeg(nlpt,0:ndiff) )
        first = .false.
        goto 2
      endif
      ! Calculate ldeg
      ldeg = 0
      do idiff = 0,ndiff
        do ilpt = 1,nlpt
          rdist(0)  =    sqrt(sum(matmul(lat,lpt(:,ilpt)+rdiff(:,idiff)           )**2))
          rdist(1:) = [ (sqrt(sum(matmul(lat,lpt(:,ilpt)+rdiff(:,idiff)-neigh(:,j))**2)),j=1,nneigh) ]
          if(all(rdist(1:)>rdist(0)-1d-6)) ldeg(ilpt,idiff) = count(abs(rdist(0)-rdist)<1d-6)
        enddo
        rdum = 0
        do ilpt = 1,nlpt
          if(ldeg(ilpt,idiff)>0) rdum = rdum + 1d0/ldeg(ilpt,idiff)
        enddo
        if(abs(rdum-nkpt)>nkpt*1d-12) Bug('Sum of weighted lattice points differs from nkpt.')
      enddo
      deallocate(rdist,neigh)

      else
        if(.not.allocated(ldeg)) Bug('ldeg not allocated.')
      endif

      if(.not.allocated(matr)) then
        if(mode==2) then ; allocate(matr(nlpt,maxlmindx*ncent*nwan,nspin3))
        else             ; allocate(matr(nlpt,nwan,nwan))
        endif        
      endif

      ! Calculate Hamiltonian at Wigner-Seitz points
      matr = 0

      allocate ( carr(nwanband),cphase(nlpt) )
      if(mode==2) allocate(cmt1(maxlmindx,ncent,nwanband,nspin3))
      do ikpt = 1,nkpt
        do ilpt = 1,nlpt
          cphase(ilpt) = exp( -img * 2*pi * dot_product(kpt(:,ikpt),lpt(:,ilpt)) ) / nkpt
        enddo
        if(mode<2) then ! mode = 0 and 1
          if(mode==0) then ; ib = ( kptp(ikpt) - 1 ) * nwanband    ! diag
          else             ; ib = ( kptp(ikpt) - 1 ) * nwanband**2 ! matrix
          endif
          do iwan1 = 1,nwan
            if(mode==0) then ! diag
              carr = uwan(:,iwan1,ikpt,ispin) * matk(ib+1:ib+nwanband)
            else             ! matrix
              carr = do_matvec ( matk(ib+1:ib+nwanband**2) , uwan(:,iwan1,ikpt,ispin) , symkpt(ikpt)>nsymt )
            endif
            do iwan = 1,nwan
              cdum = dotprod(uwan(:,iwan,ikpt,ispin),carr)
              do i = 1,nlpt
                matr(i,iwan,iwan1) = matr(i,iwan,iwan1) + cphase(i) * cdum
              enddo
            enddo
          enddo
        else ! mode = 2
          do ib = 1,nwanband
            call wavefunction_mt(cmt1(:,:,ib,:),maxlmindx,0,ib-1+wanbandi,ikpt,ispin)
          enddo
          do iwan = 1,nwan
            do s = 1,2
              do ic = 1,ncent
                itype = atype(ic)
                do lm = 1,sum( [ (nindx(l,itype)*(2*l+1) , l=0,lcut(itype)) ] )
                  cdum = sum( cmt1(lm,ic,:,s) * uwan(:,iwan,ikpt,ispin) )
                  lm1  = ((iwan-1)*ncent+(ic-1)) * maxlmindx + lm ! composite index lm1: lm,ic,iwan
                  do i = 1,nlpt
                    matr(i,lm1,s) = matr(i,lm1,s) + cphase(i) * cdum
                  enddo
                enddo
              enddo
              if(.not.l_soc) exit
            enddo            
          enddo
        endif
      enddo
      deallocate( carr,cphase )
      if(mode==2) deallocate(cmt1)

# ifdef WS_wancent
      if(mode<2) then
        do ilpt = 1,nlpt
          do i = 1,nwan
            do j = 1,nwan
              if(     ldeg(ilpt,diff(i,j))==0) then ; matr(ilpt,i,j) = 0
              else if(ldeg(ilpt,diff(i,j))> 1) then ; matr(ilpt,i,j) = matr(ilpt,i,j) / ldeg(ilpt,diff(i,j))
              endif
            enddo
          enddo
        enddo
      else
        do ilpt = 1,nlpt
          matr(ilpt,:,:) = matr(ilpt,:,:) / ldeg(ilpt,0)
        enddo
      endif
# else
      do ilpt = 1,nlpt
        matr(ilpt,:,:) = matr(ilpt,:,:) / ldeg(ilpt,0)
      enddo
# endif

      contains

      function do_matvec(mat,vec,ltrs)
      use, intrinsic :: iso_fortran_env
      implicit none
      complex_dp             :: do_matvec(nwanband)
      complex_dp, intent(in) :: mat(nwanband,nwanband),vec(nwanband)
      logical,    intent(in) :: ltrs
      if(ltrs) then ; do_matvec = matvec(transpose(mat),vec)
      else          ; do_matvec = matvec(mat,vec)
      endif
      end function do_matvec

      end

c --------------

c     Returns in matk the Wannier interpolated matrix (inverse Fourier transformation from matr as obtained from wannier_spatial)
c
      subroutine wannier_interpolate(matk,matr,dim,lpt,nlpt,kvec)
      use global, only: img,pi,nwan
      use, intrinsic :: iso_fortran_env
      implicit none
      integer,    intent(in)  :: nlpt,lpt(3,nlpt),dim
      real_dp,    intent(in)  :: kvec(3)
      complex_dp, intent(out) :: matk(dim)
      complex_dp, intent(in)  :: matr(nlpt,dim)
      complex_dp              :: cphase(nlpt)
      integer                 :: ilpt,i
      do ilpt = 1,nlpt
        cphase(ilpt) = exp( img * 2*pi * dot_product(kvec,lpt(:,ilpt)) )
      enddo
      do i = 1,dim
        matk(i) = sum(cphase*matr(:,i))
      enddo
      end

c --------------





c BACKUP
# if 0
      do ic = 1,ncent
        dist(ic) = sum ( matmul(lat,cent(:,ic)-cent(:,1))**2 )
        write(*,'(I2,4F10.5)') ic,matmul(lat,cent(:,ic)-cent(:,1)) , dist(ic)
      enddo
      dev   = sum(sqrt(dist))/(ncent-1)
      backf = 0
      l     = 5
      do ic = 2,ncent
      do i1 = -l,l
      do i2 = -l,l
      do i3 = -l,l
        rdum = sum(matmul(lat,cent(:,ic)-cent(:,1)+[i1,i2,i3])**2)
        if(rdum<dist(ic)-1d-6) then
          backf(:,ic) = [i1,i2,i3]
          write(*,'(I2,2F10.5)') ic,dist(ic),rdum
          dist(ic)    = rdum
          write(*,'(I2,3F10.5)') ic,matmul(lat,cent(:,ic)-cent(:,1)+[i1,i2,i3])
        endif
      enddo
      enddo
      enddo
      enddo
# endif
# if 0
      subroutine wannier_interpolation
      use global
      use util, only: chr
      use wrapper
      use file
      use, intrinsic :: iso_fortran_env
      implicit none
      complex_dp              :: carr(nkpt),carr1(nwanband),cdum
      complex_dp              :: hamiltonian(nwan,nwan)
      complex_dp, allocatable :: hamiltonian_r(:,:,:),cphase(:)
      real_dp,    allocatable :: bst(:,:,:),rpath(:)
      real_dp                 :: kvec(3),kvec0(3),eig(nwan)
      integer,    allocatable :: lpt(:,:),neigh(:,:),ldeg(:)
      integer                 :: nlpt,ilpt,nneigh
      integer                 :: iunit,iunitr,iunitk,iunitl
      integer                 :: ivec(3)
      integer                 :: ipath,ikpt,ikpt_path,ispin
      integer                 :: i,j,iwan,iwan1
      integer                 :: i1,i2,i3
      logical                 :: first,compl,klist
      character(80)           :: line
      real                    :: cputime

      write(6,'(/A'NoA) 'Wannier interpolation... '
      call cpu_time(cputime)

      first  = .true.
 1    nneigh = 0
      do i1 = -3,3
        do i2 = -3,3
          do i3 = -3,3
            if(any([i1,i2,i3]/=0)) then
              nneigh = nneigh + 1
              if(.not.first) neigh(:,nneigh) = [i1,i2,i3]*nkpt3
            endif
          enddo
        enddo
      enddo
      if(first) then
        allocate ( neigh(3,nneigh) )
        first = .false.
        goto 1
      endif

      ! Define Wigner-Seitz cell (of supercell) (->lpt)
      allocate ( rpath(0:nneigh) )
      rpath = 0
      first = .true.
 2    compl = .false.
      nlpt  =  0
      i     = -1
      do while(.not.compl)
        i     = i + 1
        compl = .true.
        do i1 = -i,i
          do i2 = -i+abs(i1),i-abs(i1)
            i3 = i-abs(i1)-abs(i2)
            do
              rpath(0)  =    sum(matmul(lat,[i1,i2,i3]           )**2)
              rpath(1:) = [ (sum(matmul(lat,[i1,i2,i3]+neigh(:,j))**2)+1d-6,j=1,nneigh) ]
              if(all(rpath(0)<rpath(1:))) then
                compl = .false.
                nlpt  = nlpt + 1
                if(.not.first) then
                  lpt(:,nlpt) = [ i1,i2,i3 ]
                  ldeg(nlpt)  = count(abs(rpath(0)-rpath)<1d-6)
                endif
              endif
              if(i3<=0) exit
              i3 = -i3
            enddo
          enddo
        enddo
      enddo
      if(first) then
        first = .false.
        allocate ( lpt(3,nlpt),ldeg(nlpt) )
        goto 2
      endif
      deallocate ( rpath,neigh )

      ! Prepare band structure plot
      ikpt_path = 0
      do ipath = 1,nkpt_path-1
        do i = 0,100
          if(i==0.and.ipath/=1) cycle
          ikpt_path = ikpt_path + 1
        enddo
      enddo
      allocate ( bst(ikpt_path,nwan,nspin1),rpath(0:ikpt_path) )

      iunit  = fopen('bands0',status='unknown')
      iunitr = fopen('hamiltonian_r',status='unknown')
      write(iunitr,'(A)') '# Real-space Wannier Hamiltonian  < W(n0) | H | W(n''R) >'

      inquire(file='klist',exist=klist)
      if(klist) then
        iunitk = fopen('hamiltonian_k',status='unknown')
        write(iunitk,'(A)') '# k-space Wannier Hamiltonian  < W(nk) | H | W(n''k) >'
      endif

      allocate ( hamiltonian_r(nlpt,nwan,nwan),cphase(nlpt) )

      do ispin = 1,nspin1

      ! Calculate Hamiltonian at Wigner-Seitz points
      hamiltonian_r = 0
      do iwan1 = 1,nwan
        do ikpt = 1,nkpt
          do ilpt = 1,nlpt
            cphase(ilpt) = exp( -img * 2*pi * dot_product(kpt(:,ikpt),lpt(:,ilpt)) )
          enddo
          carr1(:) = uwan(:,iwan1,ikpt,ispin) * ene(wanbandi:wanbandf,ikpt,ispin)
          do iwan = 1,nwan
            cdum                        = dot_product(uwan(:,iwan,ikpt,ispin),carr1) / nkpt
            hamiltonian_r(:,iwan,iwan1) = hamiltonian_r(:,iwan,iwan1) + cdum * cphase
          enddo
        enddo
      enddo
      if(nspin1==2) then
        if(ispin==1) write(iunitr,'(/A)') 'Spin up'
        if(ispin==2) write(iunitr,'(/A)') 'Spin down'
      endif
      do ilpt = 1,nlpt
        write(iunitr,*)
        write(iunitr,'(4I3)')                         lpt(:,ilpt),ldeg(ilpt)
        write(iunitr,'('//chr(nwan)//'(2F15.10,3X))') transpose(hamiltonian_r(ilpt,:,:))
      enddo

      do ilpt = 1,nlpt
        hamiltonian_r(ilpt,:,:) = hamiltonian_r(ilpt,:,:) / ldeg(ilpt)
      enddo

      ! Calculate band structure (->bst)
      ikpt_path = 0
      rpath(0)  = 0
      kvec0     = kpt_path(:,1)
      do ipath = 1,nkpt_path-1
        do i = 0,100
          if(i==0.and.ipath/=1) cycle
          kvec             = ( kpt_path(:,ipath)*(100-i) + kpt_path(:,ipath+1)*i ) / 100
          ikpt_path        = ikpt_path + 1
          rpath(ikpt_path) = rpath(ikpt_path-1) + sqrt(sum(matmul(rlat,kvec-kvec0)**2))

          do ilpt = 1,nlpt
            cphase(ilpt) = exp( img * 2*pi * dot_product(kvec,lpt(:,ilpt)) )
          enddo

          do iwan1 = 1,nwan
            do iwan = 1,nwan
              hamiltonian(iwan,iwan1) = sum ( cphase * hamiltonian_r(:,iwan,iwan1) )
            enddo
          enddo

          call diagonalize(eig,hamiltonian)
          bst(ikpt_path,:,ispin) = eig

          kvec0 = kvec
        enddo
        if(ipath/=nkpt_path-1) write(iunit,'(''#'',F15.10,3F10.5)') rpath(ikpt_path),kvec
      enddo

      ! If klist exists, write Hamiltonians to file hamiltonian_k
      if(klist) then
        iunitl = fopen('klist',status='old')
        do
          read(iunitl,'(A)',iostat=i) line
          if(i/=0) exit
          if(line==' '.or.line(:1)=='#') cycle
          read(line,*,iostat=i) kvec
          if(i/=0) Error('Error in klist.')
          do ilpt = 1,nlpt
            cphase(ilpt) = exp( img * 2*pi * dot_product(kvec,lpt(:,ilpt)) ) / nkpt
          enddo
          do iwan1 = 1,nwan
            do iwan = 1,nwan
              hamiltonian(iwan,iwan1) = sum ( cphase * hamiltonian_r(:,iwan,iwan1) )
            enddo
          enddo
          call diagonalize(eig,hamiltonian)
          write(iunitk,*)
          write(iunitk,'(3F10.5)') kvec
          write(iunitk,'('//chr(nwan)//'(2F15.10,3X))') transpose(hamiltonian)
          write(iunitk,*)
          write(iunitk,'(F15.10)') eig
        enddo
        call fclose(iunitl)
      endif

      enddo ! Spin loop
      deallocate ( lpt,ldeg,hamiltonian_r,cphase )
      call fclose(iunitr)
      if(klist) call fclose(iunitk)

      do iwan = 1,nwan
        if(nspin1==1) then ; write(iunit,'(2F15.10)') (rpath(i),bst(i,iwan,1),i=1,ikpt_path)
        else               ; write(iunit,'(3F15.10)') (rpath(i),bst(i,iwan,:),i=1,ikpt_path)
        endif
        write(iunit,*)
      enddo

      deallocate ( bst )
      call fclose(iunit)

      call cpu_done(cputime)

      write(6,'(A)') 'Band structure written to bands.'
      if(klist) then
        write(6,'(A)') 'Hamiltonians written to hamiltonian_r (real-space) and hamiltonian_k (k-space).'
      else
        write(6,'(A)') 'Hamiltonians (real-space) written to hamiltonian_r.'
      endif

      end
# endif

c --------------

