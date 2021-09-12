c--------------------------------------------------------------------------------
c Copyright (c) 2021 Peter Grünberg Institut, Forschungszentrum Jülich, Germany
c This file is part of SPEX and available as free software under the conditions
c of the MIT license as expressed in the LICENSE file in more detail.
c--------------------------------------------------------------------------------

# include "cppmacro.h"
# include "jobtype.h"
# include "restype.h"

      subroutine correlation_wannier(job1)

      use global
      use arrays
      use key
      Mpi( use Mwrapper )

      use, intrinsic :: iso_fortran_env
      implicit none
      type(jobtype), intent(in) :: job1
      real_dp                   :: contour_par(4)
      real_dp,    allocatable   :: freqt(:)
      complex_dp, allocatable   :: matrix(:,:,:)
      complex_dp                :: screen(nwan**2,nwan**2)
      integer,    parameter     :: spin0 = 3
      integer                   :: nfreqt,n
      integer                   :: ikpt,typ
      logical                   :: lkptc(nkpt),ldum
# ifdef MPI
      type(c_ptr)               :: ptr
      complex_dp,  pointer_cnt  :: Ninit_selfc(:)
      integer                   :: Merr
# endif
# ifdef CHECK_SENERGY
      integer                   :: i
# endif
# include "interface/define_freq.inc"

      if(job1%type/=J_GT) return

      Rwrite(6,'(//A)') '### subroutine: correlation_wannier ###'

c
c     Define frequencies
      Rbegin
      call define_freq(job1,contour_par,spin0-2,.false.,.true.)
      call getkey(inp,'FPADE',ldum,section='SUSCEP',default=.false.)
      if(ldum) then ; call getfreqt(freqt,nfreqt) ; deallocate(freq) ; allocate(freq(nfreq)) ;  freq = freqt ; endif
      Rend
      Mpi( call Mbroadcast )
      if(job1%full) Error('GT FULL not implemented yet.')

      Allocate_( matrix, ( nwan**2,nwan**2,nfreq+nfreqc ) )
      if(oselfc==1) then ; Nallocate0( selfc, (S_ size(job1%band),nfreq  S_) )
      else               ; Nallocate0( selfc, (S_ size(job1%band),nfreqr S_) )
      endif

c
c     Get screened interaction from spex.cou
c     (We calculate W and T for spin0 and transform to the other spin if needed)
      Rcall readcou(ikpt,matrix,screen,.false.,.true.,spin0,[0,0,0],1,(0d0,0d0),1)
      Rcall reorder4(screen)
      Rif(ikpt==0)     Error('spex.cou not found.')
      Rif(ikpt/=nkpti) Warn('W incomplete or calculated with different k-point set.')

      Rbegin
      lkptc = .false.
      if(iand(restart,R_sigt)/=0) then
        call read_dump(2,lkptc,nkpt,'spex.sigt',job1)
        if(allocated(irrep_wan).and.any(lkptc(nkpti+1:))) Error('Restart from calculation without IRREP?')
      endif
      n = nkpt ; if(allocated(irrep_wan)) n = nkpti
      if(iand(restart,W_sigt)/=0 .and. all(.not.lkptc(:n)) ) then
        call write_dump(2,[0],0,'spex.sigt',job1 MpiC(0)MpiC(0) ) ! truncate sig file
      endif
      Rend

      Mpi(  call Mcast(screen) )
      MpiO( call Mcast(selfc, comm=Ocomm) ) ; Nfence(selfc)
      Mpi(  call Mcast(lkptc) )

      do ikpt = 1,nkpt

        if(lkptc(ikpt)) cycle
        if(allocated(irrep_wan).and.ikpt>nkpti) exit

        Rwrite(6,'(/A)')             '========='
        Rwrite(6,'(/A,I4,A,3F10.5)') 'K point (',ikpt,' ):',kpt(:,ikpt)
        if(iand(restart,R_cot)/=0) then ; typ = 5 ; call readcot(typ)
        else                            ; typ = 0
        endif
        if(typ<4) call susceptibility_wannier(ikpt,spin0,allocated(irrep_wan),nfreq+nfreqc,[img*freq,freqc],matrix)
        if(typ<5) call bethesalpeter_t(matrix,screen,[img*freq,freqc],nfreq+nfreqc,ikpt)
        if(typ<4.and.iand(restart,W_cot)/=0) call writecot(5)
        call selfenergy_wannier(job1,ikpt,matrix,spin0,allocated(irrep_wan))
        Rif(iand(restart,W_sigt)/=0) call write_dump(2,[ikpt],1,'spex.sigt',job1 MpiC(0)MpiC(0) )
# ifdef CHECK_SENERGY
        Nfence(selfc)
        Rbegin
        write(6,'(F10.5,2F15.8)') (freqr(i),selfc(1,i),i=1,nfreqr)
        selfc = 0
        Rend
        Nfence(selfc)
#   warning CHECK_SENERGY defined
# endif

      enddo

      Deallocate_(matrix)

      contains

c ------------------------

c     readcot:  Reads  spex.cot
c     writecot: Writes spex.cot

# define ERRSTOP if(Herr/=0) Error('Fatal HDF5 error.')

      subroutine readcot(typ)
      use util, only: chr
      use hdf5
      use Hwrapper
      use, intrinsic :: iso_fortran_env
      implicit none
      integer, intent(inout)  :: typ
      integer(HID_T)          :: Hfile,Hkpt
      integer                 :: Herr
      character(38)           :: Hname
      character(2)            :: cjob
      character(6)            :: label1
      integer                 :: iarr(2),i
      real_dp                 :: rarr(4)
      real_dp,    allocatable :: freq1(:)
      complex_dp, allocatable :: freqc1(:)
      logical                 :: ldum,lnew
      logical, save           :: first = .true.
# ifdef HDF5ser
      MnoR( return ) ! only root writes in case of HDF5ser
# endif
      cjob = ' ' !; if(maxval(job%indx)>1) write(cjob,'(I2.2)') job1%indx
      ifR inquire(file='spex'//trim(cjob)//'.cot',exist=ldum)
      Mpi( call Mcast(ldum) )
      if(.not.ldum) then
        typ = -1
        return
      endif
      ! Open file
      call hdf_fopen(Hfile,'spex'//trim(cjob)//'.cot',0)
      ! Open k-point group (return if it does not exist)
      write(Hname,'(F12.10,'','',F12.10,'','',F12.10)') kpt(:,ikpt)
      call h5lexists_f(Hfile,Hname,ldum,Herr) ; ERRSTOP
      if(.not.ldum) then
        typ = -1
        call hdf_fclose(Hfile)
        return
      endif
      call h5gopen_f(Hfile,Hname,Hkpt,Herr) ; ERRSTOP
      call h5aexists_f(Hkpt,'type',ldum,Herr) ; ERRSTOP
      if(.not.ldum) Error('Type not found in spex.cot')
      ! Read file attributes
      if(first) then
        first = .false.
        call hdf_rdwr_a(Hfile,'irrep',0,iarr(1))
        if(iarr(1)==0.and..not.allocated(irrep_wan)) Warn('Data in spex.cot is IRREP-symmetrized.')
        if(iarr(1)==1.and.     allocated(irrep_wan)) Warn('Data in spex.cot is not IRREP-symmetrized.')
        call h5aexists_f(Hfile,'nfreqt',ldum,Herr) ; ERRSTOP
        if(ldum) then ; call hdf_rdwr_a(Hfile,'nfreqt',0,iarr(1:1)) ; iarr(2) = 0
        else          ; call hdf_rdwr_a(Hfile,'nfreq,nfreqc',0,iarr)
        endif
        lnew = .false.
        if(nfreq/=0) then
          if(iarr(1)==0) Error('No imaginary frequency mesh in spex'//trim(cjob)//'.cot.')
          allocate(freq1(iarr(1)))
          call h5aexists_f(Hfile,'freqt',ldum,Herr) ; ERRSTOP
          if(ldum) then ; call hdf_rdwr_a(Hfile,'freqt',0,freq1)
          else          ; call hdf_rdwr_a(Hfile,'freq',0,freq1)
          endif
          ldum = .false.
          if(iarr(1)/=nfreq) then
            nfreq = iarr(1) ; deallocate(freq) ; allocate(freq(nfreq))
            lnew  = .true.
            ldum  = .true.
          else if(any(freq1/=freq)) then
            ldum  = .true.
          endif
          if(ldum) then
            freq = freq1
            RWarn('Imaginary frequency mesh replaced by the one from spex'//trim(cjob)//'.cot: '//)
     &        Chr(freq(1))//' : '//Chr(freq(nfreq))
            Rif(nfreq==1) Warn('Imaginary frequency mesh contains a single point.')
          endif
          deallocate(freq1)
        endif
        if(nfreqc/=0) then
          if(iarr(2)==0) Error('No real frequency mesh in spex'//trim(cjob)//'.cot.')
          allocate(freqc1(iarr(2)))
          call hdf_rdwr_a(Hfile,'freqc',0,freqc1)
          ldum = .false.
          if(iarr(2)/=nfreqc) then
            nfreqc = iarr(2) ; deallocate(freqc) ; allocate(freqc(nfreqc))
            lnew   = .true.
            ldum   = .true.
          else if(any(freqc1/=freqc)) then
            ldum   = .true.
          endif
          if(ldum) then
            freqc = freqc1
            if(any(abs(imag(freqc))>1d-12)) Error('Frequencies for T interpolation have imaginary part.')
            RWarn('Frequency mesh for T interpolation replaced by the one from spex'//trim(cjob)//'.cot: '//)
     &        Chr(real(freqc(1)))//' : '//Chr(real(freqc(nfreqc)))
          endif
          deallocate(freqc1)
        endif
        if(lnew) then
          Deallocate_(matrix) ; Allocate_(matrix,(nwan**2,nwan**2,nfreq+nfreqc))
        endif
      endif
      ! Read k-point attributes
      call hdf_rdwr_a(Hkpt,'type',   0,i) ; if(i>typ) Error('Matrix type inconsistent.')
      typ = i
      call hdf_rdwr_a(Hkpt,'label',  0,label1)
      call hdf_rdwr_a(Hkpt,'dim',    0,i) ; if(i/=nwan) Error('Number of Wannier functions incorrect.')
      call hdf_rdwr_a(Hkpt,'contour',0,rarr)
      call hdf_rdwr_a(Hkpt,'oselfc', 0,i)
      Rif(i/=oselfc) Warn('Data taken from different contour type.')
      Rif(any(rarr/=contour_par).and.oselfc>1) Warn('Contour parameters different.')
      ! Read data
      call hdf_rdwr(Hkpt,'matrix',0,matrix,0)
      if     (typ==4) then ; Rwrite(6,'(A'NoA) 'Susceptibility '
      else if(typ==5) then ; Rwrite(6,'(A'NoA) 'T '
      else                 ; RBug('Unknown type index: '//trim(chr(typ))//'.')
      endif
      Rwrite(6,'(A'NoA) 'matrices read from spex'//trim(cjob)//'.cot'
      Rif(label1/=job1%label) write(6,'(A'NoA) ' (job '//trim(label1)//')'
      Rwrite(6,*)
      ! Close all
      call h5gclose_f(Hkpt,Herr) ; ERRSTOP
      call hdf_fclose(Hfile)
      end subroutine readcot

      subroutine writecot(typ)
      use util, only: chr
      use hdf5
      use Hwrapper
      use, intrinsic :: iso_fortran_env
      implicit none
      integer, intent(in) :: typ
      integer(HID_T)      :: Hfile,Hkpt
      integer             :: Herr
      character(38)       :: Hname,fname
      character(2)        :: cjob
      logical             :: ldum
      integer             :: i
      cjob  = ' ' !; if(maxval(job%indx)>1) write(cjob,'(I2.2)') job1%indx
      fname = 'spex'//trim(cjob)//'.cot' !; Mpi( if(Mcolor/=0) fname(len_trim(fname)+1:) = '.c'//trim(chr(Mcolor)) )
      ! Open file
      call hdf_fopen(Hfile,trim(fname),1)
      ! Open k-point group
      write(Hname,'(F12.10,'','',F12.10,'','',F12.10)') kpt(:,ikpt)
      call h5lexists_f(Hfile,Hname,ldum,Herr) ; ERRSTOP
      if(ldum) Error('Group '//Hname//' already exists in '//trim(fname)//'.')
      call h5gcreate_f(Hfile,Hname,Hkpt,Herr) ; ERRSTOP
      select case(typ)
        case (4)     ; Rwrite(6,'(A'NoA) 'Susceptibility matrices '
        case (5)     ; Rwrite(6,'(A'NoA) 'T matrices '
        case default ; RBug('Unknown type index: '//trim(chr(typ))//'.')
      end select
      Rwrite(6,'(A)') 'appended to '//trim(fname)
      ! Write file attributes
      i = 1 ; if(allocated(irrep_wan)) i = 0
      call hdf_rdwr_a(Hfile,'irrep',2,i)
      call hdf_rdwr_a(Hfile,'nfreq,nfreqc',2,[nfreq,nfreqc])
      if(nfreq >0) call hdf_rdwr_a(Hfile,'freq',2,freq)
      if(nfreqc>0) call hdf_rdwr_a(Hfile,'freqc',2,freqc)
      ! Write k-point attributes
      call hdf_rdwr_a(Hkpt,'type',   1,typ)
      call hdf_rdwr_a(Hkpt,'label',  1,job1%label)
      call hdf_rdwr_a(Hkpt,'dim',    1,nwan)
      call hdf_rdwr_a(Hkpt,'oselfc', 1,oselfc)
      call hdf_rdwr_a(Hkpt,'contour',1,contour_par)
      ! Write data
      call hdf_rdwr(Hkpt,'matrix',1,matrix)
      ! Close all
      call h5gclose_f(Hkpt,Herr) ; ERRSTOP
      call hdf_fclose(Hfile)
      end subroutine writecot

c ------------------------

      ! Copied from Mathias.
      subroutine getfreqt (freqt,nfreqt)
      use global, only: inp
      use file
      use key
      use, intrinsic :: iso_fortran_env
      implicit none
      real_dp, intent(out), allocatable :: freqt(:)
      integer, intent(out)              :: nfreqt
      character(:),         allocatable :: charr(:)
      real_dp                           :: ch2r
      real_dp                           :: fpade(3)
      real_dp                           :: x,maxfreqt
      integer                           :: i

      allocate(character(80) :: charr(3))
      call getkey(inp,'FPADE', charr, section='SUSCEP')!,  default=['50.','1.02','0.5'] )
      fpade(1) = ch2r(charr(1),'FPADE') ; if(fpade(1)<=0) stop 'getfreqt: First argument to FPADE must be positive.'
      fpade(2) = ch2r(charr(2),'FPADE')
      fpade(3) = ch2r(charr(3),'FPADE')
      deallocate(charr)

      if(fpade(1)<=0) Error('Define frequency mesh for Pade approximation. use keyword PADE' )
      !maxfreqt  = maxval(ene(wanbandf,:,:)) - minval(ene(wanbandi,:,:))
      maxfreqt = fpade(3)
      nfreqt = fpade(1)
      x = fpade(2)
      write(6,'(A,I5,2F12.8)') 'Frequency mesh for Pade approximant:', nfreqt,x,maxfreqt
      allocate(freqt(nfreqt))
      freqt = 0
      do i = 1,nfreqt
        freqt(i) = (x**i - x)/(x**nfreqt-x) * maxfreqt ; write(6,'(I3,F20.10)') i,freqt(i)
      enddo

      end subroutine getfreqt

c ------------------------

# ifdef MPI
      subroutine Mbroadcast
      call Mcast(contour_par); call Mcast(oselfc)
      call Mcastl(freqt); call Mcast(nfreqt)
      call Mcastl(freq) ; call Mcast(nfreq)
      call Mcastl(freqr); call Mcast(nfreqr)
      call Mcastl(freqc); call Mcast(nfreqc)
      end subroutine Mbroadcast
# endif

      end

