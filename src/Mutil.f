c--------------------------------------------------------------------------------
c Copyright (c) 2021 Peter Grünberg Institut, Forschungszentrum Jülich, Germany
c This file is part of SPEX and available as free software under the conditions
c of the MIT license as expressed in the LICENSE file in more detail.
c--------------------------------------------------------------------------------

# include "cppmacro.h"

c     --------

c MPI write routines:
c     Mwrite   : Each process in Mcomm writes a line (line0). The actual output is done by root.
c     Mwrite_n : Each process writes arbitrarily many lines. The actual output is done by global root.
c                (The lines are stored in a file unless the process is also the global root process and non-ordered output is requested.)

c     If Mrank/=0, send line0 to root process (Mrank=0).
c     If Mrank==0, write line0, receive line0 from all other processes and write them.
c     If any line0=='skipline', this line is skipped (not written).
      subroutine Mwrite(line0)
      use global, only: Mrank,Msize,Mcomm
      use, intrinsic :: iso_fortran_env
      implicit none
      include 'mpif.h'
      character(*), intent(in) :: line0
      character(80)            :: line
      integer                  :: Merr,Mstat(mpi_status_size)
      integer                  :: i
      line = line0
      if(Mrank/=0) then
        call mpi_ssend(line,80,mpi_character,0,0,Mcomm,Merr)
      else
        if(line/='skipline') write(6,'(A)') trim(line)
        do i = 1,Msize-1
          call mpi_recv(line,80,mpi_character,mpi_any_source,0,Mcomm,Mstat,Merr)
          if(line/='skipline') write(6,'(A)') trim(line)
        enddo
      endif
      end

c     Each process writes output to unit 6. The routine takes care that the blocks of output (i.e., several lines of output
c     called "blocks" below) are written separately by the global root (Mrank0==0) to unit 6 (avoiding unsynchronized writing).
c
c     mode=0 : initialization: processes open files for blocks; called before writing of blocks
c     mode=1 : blocks are collected and written (if ready) by global root; called after writing of blocks
c     mode=2 : finalize: all remaining blocks are written by global root (the other processes wait in mpi_wait); called at the end
c              (exception for ii<0: there are no remaining blocks, finalization only deallocates an array)
c
c     There are several kinds of block orderings defined by the argument "ii". Here, "direct write" means: progress of global root is visible.
c     
c     ii     : if ii/=0, abs(ii) = index of block to be written
c              if ii>0,  blocks are ordered     (no MPI synchronization overhead, no direct write)
c              if ii==0, blocks are not ordered (no MPI synchronization overhead, direct write)
c              if ii<0,  blocks are ordered     (MPI synchronization overhead, direct write)
c                        global root must always write the first block of a group of blocks
c     nn     : global root:     nn = number of remaining blocks to be written (nn is intent(inout)!)
c              other processes: nn = number of blocks to be written (nn unchanged)
c              if ii<0, nn should be identical to the number of processes that write      
c
c     Typical usage: ("do" is used because a do loop is the most common use case, but any kind of loop or no loop at all can be used)
c
c     (a) ii==0: unordered blocks, only the total number of blocks (nn) is needed
c       nn = ...
c       do
c         call Mwrite_n(0,0,nn)
c         ... arbitrarily many processes write to unit 6 ...
c         call Mwrite_n(1,0,nn)
c       enddo
c       call Mwrite_n(2,0,nn)
c
c     (b) ordered blocks (ii>0), the total number of blocks (nn) and an index (ii running from 1 to nn) are needed
c       nn = ...
c       do
c         ii = ...      
c         call Mwrite_n(0,ii,nn)
c         ... arbitrarily many processes write to unit 6 (each block gets a unique index ii) ...
c         call Mwrite_n(1,ii,nn)
c       enddo
c       call Mwrite_n(2,1,nn) ! 2nd argument must be positive, otherwise not referenced
c      
c     (c) ordered blocks (ii<0), the total number of blocks (nn) and an index (ii running from -1 to -nn) are needed
c       do
c         nn = ... ! typically, nn = Msize
c         ii = ...
c         call Mwrite_n(0,-ii,nn)
c         ... all processes write to unit 6 (nn outputs with index ii running from -1 to -nn) ...
c         call Mwrite_n(1,-ii,nn) ! global root writes the other nn-1 blocks (including MPI synchronization)
c         call Mwrite_n(2,-ii,nn) ! must be here if nn is different in next loop
c       enddo
c       call Mwrite_n(2,-1,nn) ! must be here if not given before / 2nd argument must be negative, otherwise not referenced;
c                              ! nn is not referenced
c
c     Note: No checks are made about correct usage. Wrong usage may lead to unpredictable behavior.
c           In terms of ease of implementation, (a) is easiest, then (b), then (c).
c      
      subroutine Mwrite_n(mode,ii,nn)
      use global
      use util,   only: chr
      implicit none
      integer, intent(in)        :: mode,ii
      integer, intent(inout)     :: nn
      integer, save, allocatable :: Mreq(:)      
      integer                    :: Merr,Mstat(mpi_status_size),color,rank,i,j,aii
      logical                    :: lflag,lfile
      aii   = abs(ii)
      lfile = ii>0 .or. Mrank0/=0
      select case(mode)
        case(0) ! initialize: Non-root process opens file where unit-6 output will be stored
          if(lfile) then
            if(.not.allocated(Mreq)) then
              allocate(Mreq(nn))
              Mreq = mpi_request_null
            endif
            i = size(Mreq) - aii + 1 ; if(ii==0) i = next_req()
            open(6,file='spex.w.mwrite.'//Chr(uin)//'r'//Chr(Mrank0)//'q'//Chr(i))
          endif
        case(1) ! write: Non-root process closes file and sends flag to root process; root process writes file when flag is received
          if(lfile) then
            flush(6)
            close(6)
            j = next_req()
            i = size(Mreq) - aii + 1 ; if(ii==0) i = j
            if(ii>=0) then ; call mpi_isend(Mrank0,1,mpi_integer,0,100+i,mpi_comm_world,Mreq(j),Merr)
            else           ; call mpi_send (Mrank0,1,mpi_integer,0,100+i,mpi_comm_world,        Merr)
            endif
          endif
          if(Mrank0==0) then
            if(ii<=0) nn = nn - 1
            do while(nn>0)
              i = 100 + nn ; if(ii==0) i = mpi_any_tag
              if(ii>=0) then ; call mpi_iprobe(mpi_any_source,i,mpi_comm_world,lflag,Mstat,Merr)
              else           ; call mpi_probe (mpi_any_source,i,mpi_comm_world,      Mstat,Merr) ; lflag = .true.
              endif
              if(lflag.and.Mstat(mpi_tag)>100) then
                call write_from_file
                nn = nn - 1
              else if(ii>=0) then
                exit
              endif              
            enddo
          endif
        case(2) ! finalize: Root process writes all remaining files; Non-root processes wait for mpi_recv to be finished (with mpi_wait)
          if(Mrank0==0.and.ii>=0) then
            do while(nn>0)
              i = 100 + nn ; if(ii==0) i = mpi_any_tag
              call mpi_probe(mpi_any_source,i,mpi_comm_world,Mstat,Merr)
              if(Mstat(mpi_tag)>100) then
                call write_from_file
                nn = nn - 1
              endif
            enddo
          endif
          if(lfile.and.allocated(Mreq)) then
            if(ii>=0) then
              do i = 1,size(Mreq)
                if(Mreq(i)/=mpi_request_null) call mpi_wait(Mreq(i),Mstat,Merr)
              enddo
            endif
            deallocate(Mreq)
          endif
        case default
          Bug('Unknown mode')
      end select

      contains

      function next_req()
      implicit none
      integer :: next_req
      do next_req = 1,size(Mreq)
        if(Mreq(next_req)==mpi_request_null) return
      enddo
      Bug('No request handle available.')
      end function next_req

      subroutine write_from_file
      use file
      implicit none
      integer        :: rank,itag,irank,ios,Merr,iunit
      character(256) :: line
      rank = Mstat(mpi_source)
      itag = Mstat(mpi_tag) - 100
      call mpi_recv(irank,1,mpi_integer,rank,100+itag,mpi_comm_world,Mstat,Merr)
      if(irank/=rank) Bug('rank sent deviates from source rank.')
      iunit = fopen('spex.w.mwrite.'//Chr(uin)//'r'//Chr(rank)//'q'//Chr(itag),status='old')
      do
        read(iunit,'(A)',iostat=ios) line ; if(ios/=0) exit
        write(6,'(A)') trim(line)
      enddo
      call fclose(iunit,status='delete')
      end subroutine write_from_file
      
      end

c     --------      

      subroutine begin_split(Mcolor)
      use global, only: Mrank,Msize,Mcomm,Mcomms
      use util,   only: reallocate
      implicit none
      include 'mpif.h'
      integer, intent(in) :: Mcolor
      integer             :: Merr,n
      n = size(Mcomms)
      if(Mcomm/=Mcomms(n)) Bug('Not in lowest communicator.')
      n = n + 1
      call reallocate ( Mcomms,n )
      call mpi_comm_split(Mcomm,Mcolor,Mrank,Mcomms(n),Merr)
      Mcomm = Mcomms(n)
      if(Mcomm/=mpi_comm_null) then
        call mpi_comm_size(Mcomm,Msize,Merr)
        call mpi_comm_rank(Mcomm,Mrank,Merr)
      endif
      end

      subroutine end_split
      use global, only: Mrank,Msize,Mcomm,Mcomms
      use util,   only: reallocate
      implicit none
      include 'mpif.h'
      integer :: Merr,n
      n = size(Mcomms)
      if(n==1)                 Bug('Cannot free mpi_comm_world.')
      if(Mcomm/=Mcomms(n))     Bug('Not in lowest communicator.')
      if(Mcomm/=mpi_comm_null) call mpi_comm_free(Mcomm,Merr)
      n = n - 1
      call reallocate ( Mcomms,n )
      Mcomm = Mcomms(n)
      call mpi_comm_size(Mcomm,Msize,Merr)
      call mpi_comm_rank(Mcomm,Mrank,Merr)
      end

      subroutine change_split(j)
      use global, only: Mrank,Msize,Mcomm,Mcomms
      implicit none
      include 'mpif.h'
      integer, intent(in) :: j
      integer                i,Merr
      do i = 1,size(Mcomms)
        if(Mcomm==Mcomms(i)) exit
      enddo
      i = i + j
      if(i<1.or.i>size(Mcomms)) Bug('Out of range of communicators.')
      Mcomm = Mcomms(i)
      if(Mcomm/=mpi_comm_null) then
        call mpi_comm_size(Mcomm,Msize,Merr)
        call mpi_comm_rank(Mcomm,Mrank,Merr)
      endif
      end

      subroutine begin_split_nodes(Mcolor)
      use global, only: Mcomms,Mcomm
      implicit none
      integer, intent(in) :: Mcolor
      call begin_split(Mcolor)
      call def_nodes
      end

      subroutine end_split_nodes
      use global, only: Mcomms,Mcomm
      implicit none
      call end_split
      call def_nodes
      end

      subroutine def_nodes
      use global,   only: Msplit,Mrank0,Msize0,Mrank,Msize,Mcomm,Nrank,Nsize,Ncomm,Orank,Osize,Ocomm,Opnt,maxmem
      use Mwrapper, only: Mcast,Msum,node_allocate
      use, intrinsic :: iso_fortran_env
      use, intrinsic :: iso_c_binding
      implicit none
      include 'mpif.h'
      integer     :: i,Merr,Mcolor
      type(c_ptr) :: ptr
      if(Mcomm/=mpi_comm_null) then
        ! Define intra-node communicator -> Ncomm
        if(Ncomm/=mpi_comm_null) call mpi_comm_free(Ncomm,Merr)
        if(Msplit<=-huge(0)+Msize0) then
          call mpi_comm_split(Mcomm,Mrank0/(Msplit+huge(0)),0,Ncomm,Merr)
        else
          call mpi_comm_split_type(Mcomm,Msplit,0,mpi_info_null,Ncomm,Merr)
        endif
        call mpi_comm_size(Ncomm,Nsize,Merr)
        call mpi_comm_rank(Ncomm,Nrank,Merr)
        ! Define inter-node communicator -> Ocomm
        if(Nrank==0) then ; Mcolor = 0
        else              ; Mcolor = MPI_UNDEFINED
        endif
        if(Ocomm/=mpi_comm_null) call mpi_comm_free(Ocomm,Merr)
        call mpi_comm_split(Mcomm,Mcolor,0,Ocomm,Merr)
        if(Nrank==0) then
          call mpi_comm_size(Ocomm,Osize,Merr)
          call mpi_comm_rank(Ocomm,Orank,Merr)
          if(Mrank==0.and.Orank/=0) Error('Root must have rank 0 in Ocomm.')
        endif
        call Mcast(Osize,comm=Ncomm)
        call Mcast(Orank,comm=Ncomm)
        ! Define pointer Mrank -> Orank
        if(allocated(Opnt)) deallocate(Opnt)
        allocate ( Opnt(0:Msize-1) )
        Opnt(Mrank) = Orank
        do i = 0,Msize-1
          call Mcast(Opnt(i),i)
        enddo
      endif
      end

c     --------

c     Function pcycle_skip(i,n,...) returns .true. if current process should work on ith loop (with a total of n loops).
c
c     The other arguments are for disabling certain process groupings:
c
c     xmem = mem1 / mem0  with
c     mem0 = available memory per node
c     mem1 = estimated memory requirement per subgroup
c     is used to avoid work divisions that would exhaust the memory. The more subgroups run in parallel, the higher
c     the memory demand. If the expected memory requirement exceeds the available memory for a particular subdivision,
c     this subdivision is discarded. (Obviously, xmem = 0d0 disables this feature.)
c     
c     If lscalapack=.true., the number of processes working on ith loop are restricted to the ones "allowed" by grid_proc.
c     (In short, lscalapack should be set .true. if ScaLAPACK takes significant cost in the loop.)
c
c     In all loops and processes, the arguments xmem and lscalapack must be identical!!      
c
c      
c     Typical usage: do i = 1,n ; if(pcycle_skip(i,n,0d0,.true.) cycle ! this if clause can be abbreviated with macro Pcycle of cppmacro.h
c                      beginSplit(i)                              ! optional communicator splitting (required for ScaLAPACK)
c                      ...
c                      endSplit
c                    enddo

      function pcycle_skip(indx,n,xmem,lscalapack)
      use global,  only: Mrank,Msize,Nsize
      use util,    only: chr
      use wrapper, only: grid_proc
      use, intrinsic :: iso_fortran_env
      implicit none
      logical             :: pcycle_skip
      integer, intent(in) :: indx,n
      real_dp, intent(in) :: xmem
      logical, intent(in) :: lscalapack
      integer             :: prim(ceiling(log(1d0*Msize)/log(2d0))),nprim      
      integer             :: ii,nn,m,i,nproc,mproc,ngrp,ncol,nrow
      real_dp             :: rdum1,rdum2
      call primfac(prim,nprim,Msize)
      if(product(prim(:nprim))/=Msize) Bug('Prime factorization of '//Chr(Msize)//' failed.')
      pcycle_skip = .true.
      ii          = indx
      nn          = n
      do
        mproc = Msize
        do m = 0,2**nprim-2
          nproc = product(prim(:nprim),mask=[(iand(m,2**i)/=0,i=0,nprim-1)])
          if(lscalapack) then
            rdum1 = sqrt(1d0*nproc)          ! for old check
            rdum2 = sqrt(0.25d0+nproc)-0.5d0 !
            call grid_proc(nrow,ncol,nproc,0,0)
            if( min ( abs(rdum1-nint(rdum1)) , abs(rdum2-nint(rdum2)) ) > 1d-10 .neqv. nproc/=nrow*ncol ) ! compare old with new check
     &      Bug('ScaLAPACK conditions are inequivalent. Please check!')
            if(nproc/=nrow*ncol) cycle ! skip nproc's that are "not allowed" by grid_proc (ScaLAPACK)
          endif
          if(nproc<Msize.and.((Nsize-1)/nproc+1)*xmem>1) cycle ! skip subgroups that would exceed available memory
          if(Msize/nproc<=nn) mproc = min(mproc,nproc)
        enddo
        ngrp = Msize/mproc
        if(ngrp>=ii) then
          pcycle_skip = Mrank/mproc /= ii-1
          return
        endif
        nn = nn - ngrp
        ii = ii - ngrp
        if(nn==0) Bug('Loop finished without finding valid subgroup.')
      enddo
      end

c     --------

c     Distributes ranks over sub-blocks of an array.
c     ndim              : Number of dimensions of array.
c     dim(:ndim)        : Absolute values: Dimensions of array, e.g., [20,8,6] for array(20,8,6).
c     scaling(dim,ndim) : External function giving the computational cost of an array block of dimensions dim(:ndim).
c     bounds(:4,:ndim)  : bounds(1:2,i) are the lower and upper bounds of ith dimension of sub-block for current Mrank.
c                         bounds(3:4,i) are the "rank" and "size" of ith dimension of sub-block for current Mrank.
c                         If bounds(:2,:)=[1,0,dim+1,huge(0)], current Mrank is taken idle for this step.
      subroutine Mblocks(bounds,dim,ndim,scaling)
      use global,   only: Msize,Mrank
      use Mwrapper, only: Msum      
      use util,     only: chr
      use, intrinsic :: iso_fortran_env
      implicit none
      integer, intent(in)  :: ndim,dim(ndim)
      integer, intent(out) :: bounds(4,ndim)
      real_dp, external    :: scaling
      logical              :: lidle
      integer              :: rnk,siz(ndim)
      integer              :: i,nidle
      if(any(dim<=0)) Bug('Dimensions must be positive.')
      call get_siz(siz,dim,ndim,Msize)
      i = product(siz)
      if(Mrank>=i) then
        call take_idle
        nidle = 1
      else
        nidle = 0
        do i = 1,ndim
          rnk = mod( Mrank/product(siz(:i-1)) , siz(i) )         
          if(siz(i)<=dim(i)) then
            bounds(1,i) =  1+rnk *dim(i)/siz(i)
            bounds(2,i) = (1+rnk)*dim(i)/siz(i)
            bounds(3,i) = rnk
            bounds(4,i) = siz(i)
          else if(rnk<dim(i)) then
            bounds(:2,i)= rnk + 1
            bounds(3,i) = rnk
            bounds(4,i) = dim(i)
          else
            nidle = 1
            call take_idle
            exit
          endif
        enddo
      endif
      call Msum(nidle)
      Rif(nidle>0) Info(Chr(nidle)//' MPI processes temporarily idle.')

      contains

      subroutine get_siz(siz,dim,ndim,size0)
      use, intrinsic :: iso_fortran_env
      implicit none
      integer, intent(in)  :: ndim,dim(ndim),size0
      integer, intent(out) :: siz(ndim)
      integer, allocatable :: prd(:,:)
      real_dp, allocatable :: cost(:)      
      real_dp              :: mincost
      integer              :: nprd,d(ndim)
      integer              :: icost,isize0,i
      mincost = huge(0d0)
      do isize0 = size0,1,-1
        call iproducts(i,  nprd,-ndim,isize0) ; allocate(prd(ndim,nprd),cost(nprd))
        call iproducts(prd,nprd, ndim,isize0)
        do i = 1,nprd
          where( mod(dim,prd(:,i))/=0 )
            d = 1
          else where
            d = 0
          end where
          cost(i) = scaling(dim/prd(:,i)+d,ndim)
        enddo
        icost = minloc(cost,1)
        if(cost(icost)<mincost) then
          mincost = cost(icost)-1d-12
          siz     = prd(:,icost)
        endif
        deallocate(prd,cost)
      enddo
      end subroutine get_siz

      subroutine take_idle
      implicit none
      bounds(1,:) = 1
      bounds(2,:) = 0
      bounds(3,:) = dim + 1
      bounds(4,:) = huge(0)      
      end subroutine take_idle

      end
      
      function scaling_linear(dim,ndim)
      use, intrinsic :: iso_fortran_env
      implicit none
      real_dp             :: scaling_linear
      integer, intent(in) :: ndim,dim(ndim)
      scaling_linear = product(dim)
      end

c     --------
c
c     Initializes MPI work sharing.
c
c     Root creates a file "spex.w.task.XXX" with MPI I/O containing the numbers 1,2,3,...,number+size.
c     Each number stands for a task. (There are "size" extra numbers to avoid that the file pointer hits EOF.)
c     Afterwards, this file is opened by all processes (in Mcomm).
c
c     number : number of tasks
c     size   : number of "workers" that will carry out the tasks (can be number of processes or subcommunicators etc.)
c     Mfile  : MPI file handle
c     Mrank  : rank of current process (only process with Mrank=0 will write the file)
c     Mcomm  : MPI communicator handle (must at least contain all processes that will later call Mtask_indx)
c
c     Routine "Mtasks_indx" reads the tasks using a shared file pointer.
c     Must be finalized with "call mpi_file_close(Mfile,Merr)".
c
      subroutine Mtask_init(number,size,Mfile,Mrank,Mcomm)
      use global,   only: Mpid,Nrank0,Mrank0
      use util,     only: chr
      use Mwrapper, only: Mcast
      implicit none
      include 'mpif.h'
      integer, intent(out) :: Mfile
      integer, intent(in)  :: number,size,Mrank,Mcomm
      integer              :: Merr,Mstat(mpi_status_size)
      integer              :: i,pid1
      integer, save        :: count = 0 ! pid1 and count create a unique filename for this round
      if(Mrank==0) then
        count = count + 1
# ifdef noPID
        pid1  = Mrank0 ! fallback if PID is not available
# else
        pid1  = Mpid(Nrank0) ! root's PID
# endif        
        call mpi_file_open(mpi_comm_self,'spex.w.mtask.'//Chr(pid1)//Chr(count),ior(mpi_mode_wronly,mpi_mode_create),
     &                     mpi_info_null,Mfile,Merr)
        do i = 1,number+size
          call mpi_file_write(Mfile,i,1,mpi_integer,Mstat,Merr)
        enddo
        call mpi_file_close(Mfile,Merr)
        i = count
      endif
      call Mcast(i,comm=Mcomm)
      call Mcast(pid1,comm=Mcomm)
      call mpi_file_open(Mcomm,'spex.w.mtask.'//Chr(pid1)//Chr(i),ior(mpi_mode_rdonly,mpi_mode_delete_on_close),
     &                   mpi_info_null,Mfile,Merr)
      end

c
c     Reads the next task from the list in "spex.task.UIN".
c
c     indx  : index of next task
c     work  : index of current worker (1..size) (needed only if size>0)
c     size  : if size>0, then size = number of workers (not checked!) and
c                        the tasks 1..size will be assigned to the workers 1..size in the first loop
c             if size=0, no ordering in the first loop
c     Mfile : MPI file handle (as returned by "Mtask_init").
c
      subroutine Mtask_indx(indx,work,size,Mfile)
      implicit none
      include 'mpif.h'
      integer, intent(out) :: indx
      integer, intent(in)  :: Mfile,work,size
      integer              :: Merr,Mstat(mpi_status_size)
      integer              :: i
      call mpi_file_read_shared(Mfile,indx,1,mpi_integer,Mstat,Merr)
      if(indx<=size) indx = work
      end

c     --------

# ifdef noMSYNC2
c     Can be used for mpi_f_sync_reg, which does not exist in some MPI implementations.
      subroutine donothing(buf)
#   ifdef CHECK_SHARED
#     warning CHECK_SHARED defined
      real r
      call random_number(r)
      call sleepqq(int(r*500))
#   endif
      end
# endif

c     --------
