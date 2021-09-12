c--------------------------------------------------------------------------------
c Copyright (c) 2021 Peter Grünberg Institut, Forschungszentrum Jülich, Germany
c This file is part of SPEX and available as free software under the conditions
c of the MIT license as expressed in the LICENSE file in more detail.
c--------------------------------------------------------------------------------

# include "cppmacro.h"

c     ----------------

c     Keeps track of the amount of used memory (mem in bytes).
c     Stops if storage exceeds a predefined maximal amount (maxmem) or falls below 0 MB (bug!).
c
c     (In case of warnings about nonzero final checkmem, type: grep <???@??_>CHECKMEM spex.out | sort -b --key=1.82,1.90
c     to check all allocate-deallocate pairs.)      
      subroutine checkmem(arraystring,addmem)
      use global
      use util,   only: chr
      use, intrinsic :: iso_fortran_env
      Mpi2( use, intrinsic :: iso_c_binding )
      implicit none
      real_dp, intent(in)       :: addmem
      character*(*), intent(in) :: arraystring
      Mpi( integer              :: Merr )
# ifdef WRTMEM
#   warning WRTMEM defined
#   ifdef noPID
#     error WRTMEM needs getpid(), which is not supported by your Fortran compiler (noPID defined)
#   endif
      character(6)              :: chr_nalloc
      real_dp                   :: total,add MpiC(private)
      logical                   :: kernel_old
      integer, parameter        :: unit = 0
      integer                   :: prv,stk,tot
      integer                   :: getpid      
# endif
      if(arraystring==' ') Bug('arraystring missing.')
      if(addmem/=0) then
        Nacc1_r_( mem,,addmem ) ! addmem added to mem (virtual node)
      endif
# ifdef WRTMEM
      chr_nalloc = ' '
      if     (addmem>0) then ; chr_nalloc = ' '//Chr(nalloc_+1) ; nalloc_ = nalloc_ + 1
      else if(addmem<0) then ; chr_nalloc = ' '//Chr(nalloc_)   ; nalloc_ = nalloc_ - 1
      endif
      add = addmem / megabyte
#   ifdef MPI
      call mem_usage(-1,tot,prv,stk,kernel_old)
      total   = dble(tot) / kilobyte + add ! VmRSS updates later
      private = dble(prv) / kilobyte
      if(addmem/=0) then
        Nacc1_r_( mem0,,addmem ) ! addmem added to mem0 (physical node)
      endif
      if(arraystring(:1)/='#') then        
        if(kernel_old) then ; private = private + min(add,0d0) ! VmDat updates immediately
        else                ; private = private + add          ! VmRSS updates later
        endif
        if(addmem/=0) then
          Nacc1_r_( mem_prv0,,addmem ) ! addmem added to mem_prv0 (private data on physical node)
        endif
      endif
      if(arraystring=='start') write(unit,'(/A,6X,A,7X,A,5X,A,5X,A,4X,A,3X,A,2X,A)')
     &  '#CHECKMEM (in MB)','add','mem','total','stack','memprv','private','arrayname #alloc'
      write(unit,'(/A,F14.6,5F10.2,2X,A)') Chf(Nrank0,'I3.3')//'@'//Chf(Orank0,'I2.2')//
     &  '_CHECKMEM',add,mem0/megabyte,total,real(stk)/kilobyte,mem_prv0/megabyte,private,
     &  adjustl(trim(arraystring))//trim(chr_nalloc)
#   else
      call mem_usage(getpid(),tot,prv,stk,kernel_old)
      total = dble(tot) / kilobyte + add ! VmRSS updates later
      if(arraystring=='start') write(unit,'(/A,2X,A,7X,A,5X,A,5X,A,2X,A)')
     &  '#CHECKMEM (in MB)','add','mem','total','stack','arrayname #alloc'
      write(unit,'(/A,F14.6,3F10.2,2X,A)')
     &  'CHECKMEM',add,mem/megabyte,total,real(stk)/kilobyte,
     &  adjustl(trim(arraystring))//trim(chr_nalloc)
#   endif
# endif
      if(addmem/=0) then
        if(mem>maxmem) then
          write(0,'(A)') 'checkmem: Maximal allowed storage exceeded when trying to allocate array '//trim(arraystring)//'.'
          write(0,'(A)') '          Increase MEM to at least '//Chr(ceiling( Mpi(((Nsize0-1)/Nsize+1)*) mem / megabyte))//' MB !'
          Mpi( call write_node_info )
          Error('Maximal allowed storage exceeded.')
        else if(mem<-1d-3) then
          write(0,'(A)')      'checkmem: Current storage fell below 0 MB when trying to deallocate array '//trim(arraystring)//'.'
          Mpi( call write_node_info )
          Error('Storage fell below 0 MB.')
        endif
      endif

# ifdef MPI
      contains
      
      subroutine write_node_info
      implicit none
      write(0,'(A'NoA)   '          Global rank '//Chr(Mrank0)//' on node rank '//Chr(Nrank)//'@'//Chr(Orank)
      if(Nsize/=Nsize0) then ; write(0,'(A)') ' (virtual) '//Chr(Nrank0)//'@'//Chr(Orank0)//' (physical)'
      else                   ; write(0,*)
      endif
      end subroutine write_node_info
# endif
      
      end

c     ----------------

c
c     Returns (all in KB)
c       tot : Total data         = VmRSS - RssFile,                          (stack + heap + shared)
c       prv : Private data,      = VmRSS - RssFile - RssShmem ( = RssAnon )  (stack + heap)
c       stk : Private stack data = VmStk                                     (stack)
c     of process with PID pid.
c     If pid==-1 (MPI), the corresponding sums over the node are returned.
c
c     Old Linux kernel (RssFile and RssShmem not available):
c       tot = VmRSS
c       prv = VmDat
c       stk = VmStk
c
c     Glossary
c     VmRSS : Resident set size (including shared, stack)
c     VmDat : Total address space (shared and stack substracted)
c
c     The call of checkmem is immediately after the allocate and before the deallocate [because of size(array)]
c     VmRSS : Updates when arrays are defined, not when they are allocated, therefore rss+add is the relevant memory
c     VmDat : Updates immediately when arrays are allocated, therefore dat+min(add,0d0) is the relevant memory
c
c# define KERNEL_OLD
# ifdef KERNEL_OLD
#   warning KERNEL_OLD defined
# endif

# ifdef WRTMEM
      recursive subroutine mem_usage(pid,tot,prv,stk,kernel_old)
      use file
      use global, only: kilobyte MpiC(Mpid) MpiC(Nsize0) MpiC(Nrank0)
      use util,   only: chr
      use, intrinsic :: iso_fortran_env
      implicit none
      integer, intent(in)  :: pid
      integer, intent(out) :: tot,prv,stk
      logical, intent(out) :: kernel_old
      character(len=200)   :: filename=' '
      character(len=80)    :: line
      character(len=8)     :: pid_char=' '
      integer              :: iunit
      integer              :: rss,dat,rss_file,rss_shmem
      integer              :: getpid
      logical              :: ifxst
#   ifdef MPI
      integer              :: tot1,prv1,stk1,i
      logical, save        :: mpid_nodef = .true.
      if(pid==-1) then
        prv = 0
        stk = 0        
        tot = 0
        if(.not.allocated(Mpid)) then
          call mem_usage(getpid(),tot,prv,stk,kernel_old)
          Warn('Mpid undefined. Memory information only for current process.')
        else
          if(mpid_nodef) then
            Warn('Mpid undefined so far, now defined.')
            mpid_nodef = .false.
          endif
          do i = 0,Nsize0-1
            call mem_usage(Mpid(i),tot1,prv1,stk1,kernel_old)
            if(i==0) then ; ifxst = kernel_old
            else          ; if(kernel_old.neqv.ifxst) Error('Different cores run on different Linux kernels.')
            endif
            !write(*,*) i,Mpid(i),dble(tot1)/kilobyte,dble(prv1)/kilobyte,dble(stk1)/kilobyte
            tot = tot + tot1
            prv = prv + prv1
            stk = stk + stk1           
          enddo
          !ifR read(*,*)
        endif
        return
      endif
#   endif      
      write(pid_char,'(I8)') pid
      filename = '/proc/'//trim(adjustl(pid_char))//'/status'
      ! Read system file
      inquire (file=filename,exist=ifxst)
      if (.not.ifxst) then
        Error('System file '//trim(filename)// Mpi(' (Nrank='//Chr(minloc(abs(Mpid-pid))-1)//')'//) ' does not exist.')
      endif
      iunit     = fopen(filename,action='read',status='old')
      rss       = -1
      dat       = -1
      stk       = -1
      rss_file  = -1
      rss_shmem = -1
      do
        read(iunit,'(A)',end=120) line
        call getval(rss,'VmRSS:')
        call getval(dat,'VmData:')
        call getval(stk,'VmStk:')
        call getval(rss_file,'RssFile:')
        call getval(rss_shmem,'RssShmem:')
      enddo
 120  call fclose(iunit)
# ifdef KERNEL_OLD
      if( .true. ) then
# else        
      if( rss_file == -1 .and. rss_shmem == -1 ) then ! old Linux kernel
# endif
        kernel_old = .true.
        if( any( [dat,rss,stk] == -1 ) ) Error('VmRSS, VmData, VmStk could not be read from status file.')
        tot = rss
        prv = dat
      else                                            ! new Linux kernel
        kernel_old = .false.
        if( any( [rss_file,rss_shmem,stk] == -1 ) ) Error('RssFile, RssShmem, VmStk could not be read from status file.')
        tot = rss - rss_file
        prv = rss - rss_file - rss_shmem
      endif      

      contains
      subroutine getval(val,key)
      implicit none
      integer                  :: val
      character(*), intent(in) :: key
      integer                  :: ind
      ind = index(line,trim(key))
      if(ind/=0) read(line(len_trim(key)+ind:),*) val
      end subroutine getval
      
      end
# endif
