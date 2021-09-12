c--------------------------------------------------------------------------------
c Copyright (c) 2021 Peter Grünberg Institut, Forschungszentrum Jülich, Germany
c This file is part of SPEX and available as free software under the conditions
c of the MIT license as expressed in the LICENSE file in more detail.
c--------------------------------------------------------------------------------

      program spex

# include "version.h"
# include "cppmacro.h"
# include "jobtype.h"

      use global
      use util, only: chr
      use readwrite, only: nlo,llo
      use file
      use m_fft
      use hdf5
      use timer_util
      use, intrinsic :: iso_fortran_env
      Mpi( use Mwrapper )

      implicit none
      real                    :: time1,time2
      integer                 :: j,ijob,jtype,err MpiC(Merr)

      mem => mem_
      mem =  0
# if defined(MPI) && defined(WRTMEM)
      mem0     => mem0_
      mem_prv0 => mem_prv0_
      mem0     =  0
      mem_prv0 =  0
# endif
      
      CHKMEM0(start)      

c
c     Catch option -x for launcher script
      call commandline

c
c     Initialize MPI
# ifdef MPI
      call mpi_init(Merr) ; CHKMEM0(MPIinit)
      call mpi_comm_size(mpi_comm_world,Msize,Merr)
      call mpi_comm_rank(mpi_comm_world,Mrank,Merr)
      allocate ( Mcomms(1) )
      Mcomm  = mpi_comm_world
      Mcomms = Mcomm
      Msize0 = Msize
      Mrank0 = Mrank
# endif
      call h5open_f(err) ; if(err<0) Error('Could not initialize HDF5.')
      CHKMEM0(HDF5init)

c
c     Header
      Rbegin
      call header(0)
      call cpu_time(time1)
      Rend

c
c     Read in data and construct mixed basis and Coulomb matrix (and xc kernel)
      call getinput
      call checkinput

      if(njob/=0) then
        if(l_soc.and.nspin==2) Error('Non-collinear magnetism implemented up to here.')
        if(any(job%type/=J_KS)) then
          if(any(job%kern=='ALDA')) then
            if(any(cores)) then ; call mixedbasis(4)
            else                ; call mixedbasis(3)
            endif
          else
            if(any(cores)) then ; call mixedbasis(2)
            else                ; call mixedbasis(1)
            endif
          endif
        endif

# if defined(__GFORTRAN__) && __GNUC__ < 5
#   warning rewritten to avoid gfortran 4.9.1 bug
        do ijob = 1,njob
          if(any(job(ijob)%type==[J_HF,J_GW,J_RPA,J_HFE,J_SX,J_COSX,J_PBE0]).and.coulstore) then
            call coulombmatrix0
            exit
          endif
        enddo
# else
        if(any([(any(job(ijob)%type==[J_HF,J_GW,J_RPA,J_HFE,J_SX,J_COSX,J_PBE0]),ijob=1,njob)]).and.coulstore)
     &                            call coulombmatrix0
# endif
# ifdef TDDFT
        if(any(job%kern=='ALDA')) call xckernelmatrix
# endif
      endif

c
c     Loop over jobs
      do ijob = 1,njob

        Rwrite(6,'(///A,I3,A)') '##### Job index:',job(ijob)%indx,' #####'

        jtype = job(ijob)%type

        if(allocated(job(ijob)%band))        call prepare_offdiag(job(ijob))

        if(any(jtype==J_EXCH).and.coulstore) call exchange(job(ijob))

        if(jtype/=J_KS)                      call correlation(job(ijob))

        if(jtype==J_GT) call correlation_wannier (job(ijob))
        if(any(jtype==[J_SUSR,J_GOLD]))                          call spectra_wannier(job(ijob))
        if(any(jtype==[J_HF,J_GW,J_SX,J_COSX,J_PBE0,J_GT,J_KS])) call quasiparticle  (job(ijob))
        if(any(jtype==[J_HFE,J_RPA]))                            call xc_energy      (job(ijob))

        call dealloc_job

      enddo

      if     (otimer==1) then ; call timer_print([''],select=.true.,title='TIMINGS',keep=.false.,prnt=.true. andR )        
      else if(otimer==2) then ; call timer_print([''],select=.true.,title='TIMINGS',keep=.false.,prnt=.true.,unit= ifMpi(Mrank,0) )
      endif

      call cpu_time(time2) ; Rwrite(6,'(/A,I10)') 'Timing:',nint(time2-time1)

      Ocall fclose(inp)
      Mpi( call mpi_barrier(Mcomm,Merr) )

c
c     Deallocate arrays

      if(allocated(uwan))      tDeallocate(uwan)
      if(associated(cmtu))     tDeallocate(cmtu)
      if(associated(cpwu))     tDeallocate(cpwu)
      if(allocated(cstep))     tDeallocate(cstep)
      if(allocated(cores))     deallocate ( cores )
      if(allocated(omit))      deallocate ( omit )
      if(allocated(kpt_path))  deallocate ( kpt_path )
      if(allocated(constfunc)) deallocate ( constfunc )
      if(allocated(cblock))    deallocate ( cblock )
      if(allocated(gridf))     deallocate ( gridf )
      if(allocated(rgrid0))    deallocate ( rgrid0 )
      if(allocated(nbasm))     deallocate ( nbasm,nindxm,basm,dbasm,gptm,pgptm,ngptm,pntgptm,lcutm )
      if(allocated(rpa_int))   deallocate ( rpa_int )
      deallocate ( job,neq,ztype,cent,pcent,tcent )
      deallocate ( lcutc,nindxc,ecore )
      deallocate ( nband,phase,ene,deg )
      deallocate ( lcut,nindx,lcutp,ubas,dubas,ebas,gpt,pgpt,ngpt,pntgpt )
      deallocate ( kpt,kptp,kindx,nkpt3,pkpt,symkpt,kptsym,gkptsym,wintgr,gauss )
      deallocate ( symtab,sym )
      deallocate ( nlh,llh,nmlh,mlh,clh,symcent )
      deallocate ( dwgn,fac,sfac,facfac,grid,rgrid )
      deallocate ( nlo,llo )
      Deallocate_( bas1  )
      Deallocate_( bas2  )
      Deallocate_( core1 )
      Deallocate_( core2 )      
      Mpi( deallocate ( Opnt ) )

# ifdef INV
      if(allocated(cfft1)) then
        deallocate ( cfft1,cfft2,cffts,cfftp,rfft1,rfft2,rffts,rfftp )
        nullify ( cfft1_w,cfft2_w,cffts_w,cfftp_w )
        call dfftw_destroy_plan(plan_fft1_w)
        call dfftw_destroy_plan(plan_fft2_w)
        call dfftw_destroy_plan(plan_ffts_w)
        call dfftw_destroy_plan(plan_fftp_w)
# else
      if(allocated(fft1)) then
        deallocate ( fft1,fft2,ffts,fftp )
        nullify ( cfft1,cfft2,cffts,cfftp,rfft1,rfft2,rffts,rfftp )
# endif
        call dfftw_destroy_plan(plan_fft1)
        call dfftw_destroy_plan(plan_fft2)
        call dfftw_destroy_plan(plan_ffts)
        call dfftw_destroy_plan(plan_fftp)
      endif

      call h5close_f(err) ; if(err<0) Error('Could not close HDF5.')

      call dealloc_final

      Nfence_( mem )
      CHKMEM0(final)
      Oif(abs(mem)>1d-6) Warn('Final checkmem returned nonzero value: '//Chr(mem/megabyte)//' MB.' Mpi(//' NodeRank: '//Chr(Orank)))
      if(nalloc_/=0) then
        Warn('Checkmem allocation/deallocation calls unmatched: '//Chr(nalloc_))
     &    Mpi(//' NodeRank: '//Chr(Orank)//', Rank: '//Chr(Mrank))
      endif
      Nfence_(mem)
      Nfree(mem)
# ifdef WRTMEM      
      Nfence_(mem0)
      Nfence_(mem_prv0)
      Nfree(mem0)
      Nfree(mem_prv0)
# endif      
      
# ifdef MPI
#   ifndef noPID
      deallocate(Mpid)
#   endif
      if(Nrank==0) call mpi_comm_free(Ocomm,Merr)
      call mpi_comm_free(Ncomm,Merr)
      call mpi_finalize(Merr)
# endif

      end

c -------------

      subroutine commandline
      use global, only: l_soc
      use hdf5
      use readwrite, only: invs,read_param,dftshort,dftlong
      implicit none
      character(10) :: arg
      integer       :: err
      call get_command_argument(1,arg)
      if(arg=='-x') then
        call write_macros
        write(6,'(A)') 'DFT short="'//trim(dftshort)//'", long="'//trim(dftlong)//'"'
        call h5open_f(err) ; if(err<0) Error('Could not initialize HDF5 (-x).')
        call read_param
        write(6,'(2(A,L1))') 'invs=',invs,', l_soc=',l_soc
        call h5close_f(err) ; if(err<0) Error('Could not close HDF5 (-x).')
        stop
      else if(arg=='-v'.or.arg=='--version') then
        call header(1)
        stop
      else if(arg=='-h'.or.arg=='--help') then
        write(6,'(A)') 'Usage: spex [OPTIONS]'
        write(6,'(A)') '       -v | --version  output version information and exit'
        write(6,'(A)') '       -h | --help     display this help and exit'
        write(6,'(A)') '       --inp=file      use ''file'' as input file instead of the default ''spex.inp'''
        write(6,'(A)') '       -o=opt          use special option ''opt'''
        write(6,'(A)') '       -w              set WRTKPT'
        write(6,'(A)') '       -x              write basic DFT flags (used by Spex launcher script)'
        stop
      endif
      end

c -------------

      subroutine write_macros
      write(6,'(A)') 'Macros:' Inv(//' INV') Mpi(//' MPI') Load(//' LOAD') Check(//' CHECK') Wrtmem(//' WRTMEM')
# ifdef WAN
     & //' WAN'
# endif
# ifdef WANv2
     & //' WANv2'
# endif
# ifdef WANnoMPI
     & //' WANnoMPI'
# endif
# ifdef noSHARED
     & //' noSHARED'
# endif
# ifdef noMSYNC
     & //' noMSYNC'
# endif
# ifdef noMSYNC2
     & //' noMSYNC2'
# endif
# ifdef HDF5ser
     & //' HDF5ser'        
# endif
# ifdef CHECK_SENERGY
     & //' CHECK_SENERGY'
# endif
# ifdef CPPTRAD
     & //' CPPTRAD'
# endif
# ifdef noWARN
     & //' noWARN'
# endif
# ifdef noPID
     & //' noPID'
# endif
# ifdef M_SIZEOF
     & //' M_SIZEOF
# endif
# ifdef noDOLLAR
     & //' noDOLLAR'
# endif
# ifdef noARITHM
     & //' noARITHM'
# endif
# ifdef unbug82065
     & //' unbug82065'
# endif
# ifdef imag
     & //' imag'
# endif
# ifdef noHOSTNM
     & //' noHOSTNM'
# endif
# ifdef noCOMPNM
     & //' noCOMPNM'
# endif
# ifdef noEXECMD
     & //' noEXECMD'
# endif      
# ifdef TIMING
     & //' TIMING'
# endif
# ifdef F08STOP
     & //' F08STOP'      
# endif
# ifdef F77STOP
     & //' F77STOP'      
# endif      
      end

c -------------

# include "config.h"
# include "make.h"

      subroutine header(mode)
      use readwrite, only: dftlong
      implicit none
      Mpi( include 'mpif.h' )
      integer,     intent(in) :: mode
      integer                 :: dt0(3),dt1(8)
      character(3), parameter :: month(12) = [ 'Jan','Feb','Mar','Apr','May','Jun','Jul','Aug','Sep','Oct','Nov','Dec' ]
      character(1), parameter :: githash = ' ' ! used if GITHASH undefined
# ifndef noHOSTNM
      character(40)           :: name
      integer                 :: hostnm,j
# endif
# ifdef MPI
      character(MPI_MAX_LIBRARY_VERSION_STRING) :: mpilib
      integer                                   :: n,Merr
# endif
      if(mode==0) then
        write(6,'(A)') '########################'
        write(6,'(A)') '######### SPEX #########'
        write(6,'(A)') '########################'
      else
        write(6,'(A'NoA) 'SPEX '
      endif
      dt0 = [PACKDATE]
      write(6,'(A,I3,'','',I5,A)') 'Version '//VERSION//'  ('//month(dt0(1)),dt0(2),dt0(3),')  '//GITHASH
      call date_and_time(values=dt1)
      write(6,'(A,I3,'','',I5,I4.2,2('':'',I2.2))') 'Execution time: '//month(dt1(2)),dt1(3),dt1(1),dt1(5:7)
# ifndef noCOMPNM
      write(6,'(A)') 'Compiler: '// COMPNM
# endif
# ifdef MPI
      call mpi_get_library_version(mpilib,n,Merr)
      if(index(mpilib,achar(10))/=0) n = index(mpilib,achar(10))-1
      write(6,'(A)') 'MPI libr: '//mpilib(:n)
# endif
# ifndef noHOSTNM
      j = hostnm(name)
      if(j==0) write(6,'(A,A)') 'Hostname: ',trim(name)
# endif
      write(6,'(A)') 'Interfaced to '//trim(dftlong)
      end

c -------------

      subroutine dealloc_job
      use global
      use arrays
      use, intrinsic :: iso_fortran_env
      implicit none
      Mpi( integer :: Merr )
      if(allocated(block))     deallocate ( block     )
      if(allocated(psub))      deallocate ( psub      )
      if(allocated(irrep_sub)) deallocate ( irrep_sub )
      if(allocated(freq))      deallocate ( freq      )
      if(allocated(freqr))     deallocate ( freqr     )
      if(allocated(freqc))     deallocate ( freqc     )      
      if(allocated(selfx))    tDeallocate ( selfx     )
      if(associated(selfc))  tNdeallocate ( selfc     )
      nullify(selfc)
      end
      
c -------------

      subroutine dealloc_final
      use global
      implicit none
      integer :: Merr
      NoLoad( if(associated(cmt)) tNdeallocate(cmt) )
      NoLoad( if(associated(cpw)) tNdeallocate(cpw) )      
      if(associated(coulomb)) tNdeallocate(coulomb)
      if(associated(vpw))     tNdeallocate(vpw)
      if(associated(vpw_xc))  tNdeallocate(vpw_xc)
      if(associated(vmt))     tNdeallocate(vmt)
      if(associated(vmt_xc))  tNdeallocate(vmt_xc)      
      end

c -------------      

