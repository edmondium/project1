c--------------------------------------------------------------------------------
c Copyright (c) 2021 Peter Grünberg Institut, Forschungszentrum Jülich, Germany
c This file is part of SPEX and available as free software under the conditions
c of the MIT license as expressed in the LICENSE file in more detail.
c--------------------------------------------------------------------------------

c
c     Wrapper functions for ScaLAPACK.
c
c     Input matrices are assumed to be present (and identical) in all processes.
c
c     Output matrix: Each process writes only its BLACS subblock to the output matrix using cyclic(mode=1).
c                    Before the subblock is written, the output matrix is initialized with zeros.
c                    In case of shared memory: Note that only those shared arrays are set to zero that are
c                      accessible for the processes entering here. There might be other instances of the 
c                      array residing in other nodes, which have to be set to zero in the calling routine.
c             win_*: window handle of output matrix, must be present if and only if output matrix is in shared memory.
c                    present:     The initialization is done with calls to Nfence and only Nroot initializes.
c                                 To avoid a deadlock, all processes sharing the memory must call the routine
c                                 (but not necessarily with the same portion of the array).
c                    not present: Each process initializes its local output matrix with zeros.
c
c     If the output matrix is in shared memory, the nodes communicator Ncomm has to be defined (e.g., using
c     begin_split_nodes / end_split_nodes). Otherwise, the nodes communicator is not needed.  
c
c     Mdiagonalize: eigenvalues are returned identically in all processes in eval
c
c     Data distribution afterwards in the following way:
c     call Msum(matrix)                                               ! distributed memory (win_* not present)
c     Nfence(matrix) ; Ocall Msum(matrix,comm=Ocomm) ; Nfence(matrix) ! shared memory      (win_* present)
c
c     A local grid is defined and exited after usage.
c
c     Usage of ScaLAPACK can be switched off here:
c# define noSCALAPACK

c     --------

# ifndef DZtyp
c
c   (A) Double precision
#   define DZtyp 1
#   include "w_scalapack.inc"
#   undef DZtyp
c
c   (B) Double complex      
#   define DZtyp 2
c
# endif      

# if DZtyp == 1
#   define DZ    d
#   define DZsub DSY
#   define DZstr 'DSY'      
#   define Datatyp real_dp
#   define Conjg dble
#   define RWORK      
# elif DZtyp == 2
#   define DZ    z
#   define DZsub ZHE
#   define DZstr 'ZHE'      
#   define Datatyp complex_dp      
#   define RWORK rwork,lrwork,      
# else
#   error DZtyp unknown
# endif

# ifndef noSCALAPACK
#   define HighAcc
#   define CLUSTERSIZE max(n/100,8)
#   define ORFAC 1d-8

c     --------

c
c     Mdiagonalize0... : Backend functions
c

c     --------
      
c     Full solution to eigenvalue problem
      subroutine Mdiagonalize0_ DZ (eval,evec,win_evec,mat,matp,lwrite)
      use global,   only: Mrank,Msize,Mcomm,Nrank
      use util,     only: chr
      use Mwrapper, only: Mcast
      implicit none
      include 'mpif.h'
      real_dp, intent(out)           :: eval(:)
      Datatyp, intent(out), optional :: evec(:,:)
      Datatyp, intent(in),  optional :: mat(:,:),matp(:)
      Datatyp, allocatable           :: mat1(:,:),evec1(:,:),work(:)
      Datatyp                        :: cdum
      integer, intent(in),  optional :: win_evec
      logical, intent(in),  optional :: lwrite
      logical                        :: wrt
# if DZtyp == 2
      real_dp, allocatable           :: rwork(:)
      integer                        :: lrwork
# endif      
      character                      :: mod
      integer                        :: n,npr,npc,info,lwork,m1,n1
      integer                        :: ctxt,rp,cp,desc(9)
      integer                        :: Merr
      integer                        :: numroc
# ifdef noMSYNC
      integer(mpi_address_kind)      :: Maddress_dummy
# endif
      wrt = .false. ; if(present(lwrite)) wrt = lwrite
      n   = size(eval)
      if     (present(mat))  then ; call check_dim(shape(mat), [n,n],      'mat', __LINE__,__FILE__)
      else if(present(matp)) then ; call check_dim(shape(matp),[n*(n+1)/2],'matp',__LINE__,__FILE__)
      else                        ; Bug('No input matrix given.')
      endif
      if(present(evec)) then
        call check_dim(shape(evec),[n,n],'evec',__LINE__,__FILE__)
        mod = 'V'
      else
        mod = 'N'
      endif
      ctxt = Mcomm
      call grid_proc(npr,npc,Msize,n,n)
      call blacs_gridinit(ctxt,'R',npr,npc)
      call blacs_gridinfo(ctxt,npr,npc,rp,cp)
      if(Mrank==0.and.wrt) write(6,'(A'NoA) ' ('//trim(chr(npr))//'*'//trim(chr(npc))//') '
      Rif(rp==-1) Bug('Root process not in BLACS.')
      beginSplit(-min(0,rp))
      if(rp/=-1) then
        m1 = numroc(n,blk,rp,0,npr)
        n1 = numroc(n,blk,cp,0,npc)
# if DZtyp == 2
        lrwork = 4*n - 2
        allocate(rwork(lrwork))
# endif        
        allocate ( mat1(m1,n1) )
        if(present(evec)) then ; allocate( evec1(m1,n1) )
        else                   ; allocate( evec1(1,1)   )
        endif
        if(present(mat)) then ; call cyclic_ DZ  (mat1,m1,n1,mat, n,n,rp,cp,npr,npc,blk,0)
        else                  ; call cyclic_ DZ p(mat1,m1,n1,matp,n,  rp,cp,npr,npc,blk,0)
        endif
        call descinit(desc,n,n,blk,blk,0,0,ctxt,m1,info) ; if(info/=0) Error('descinit failed.')
        call p DZsub ev(mod,'U',n,mat1,1,1,desc,eval,evec1,1,1,desc,cdum,-1, RWORK info)
        if(info/=0) Error('P'//DZstr//'EV failed (first call).')
        lwork = nint(dble(cdum))
        allocate ( work(lwork) )
        call p DZsub ev(mod,'U',n,mat1,1,1,desc,eval,evec1,1,1,desc,work,lwork, RWORK info)
        if(info/=0) Error('P'//DZstr//'EV failed with info: '//chr(info))
        call blacs_gridexit(ctxt)
      endif
      if(present(evec)) then
        if(present(win_evec)) then ! evec is in shared memory
          Nfence(evec)
          ifO evec = 0
          Nfence(evec)
        else                       ! evec is in distributed memory
          evec = 0 
        endif
      endif
      if(rp/=-1) then
        if(present(evec)) call cyclic_ DZ (evec1,m1,n1,evec,n,n,rp,cp,npr,npc,blk,1)
        deallocate ( mat1,evec1,work )
# if DZtyp == 2
        deallocate ( rwork )
# endif
      endif
      endSplit
      call Mcast(eval)
      end subroutine Mdiagonalize0_ DZ

c     --------

c     Selected solutions to eigenvalue problem
      subroutine Mdiagonalize0_ DZ x(eval,evec,win_evec,mat,matp,ir1,ir2,r1,r2,lwrite)
      use global,   only: Mrank,Msize,Mcomm,Nrank
      use util,     only: chr
      use Mwrapper, only: Mcast
      implicit none
      include 'mpif.h'
      real_dp, intent(out)             :: eval(:)
      Datatyp, intent(out),   optional :: evec(:,:)
      Datatyp, intent(in),    optional :: mat(:,:),matp(:)
      real_dp, intent(in),    optional :: r1,r2
      integer, intent(in),    optional :: ir1,win_evec
      integer, intent(inout), optional :: ir2
      logical, intent(in),    optional :: lwrite
      logical                          :: wrt
      Datatyp, allocatable             :: mat1(:,:),evec1(:,:),work(:)
      real_dp, allocatable             :: gap(:)
# if DZtyp == 2
      real_dp, allocatable             :: rwork(:)
      integer                          :: lrwork
# endif      
      real_dp                          :: abstol,pdlamch
      real_dp                          :: r1_,r2_,orfac=ORFAC
      character                        :: mod,cmod
      integer, allocatable             :: iwork(:),ifail(:),iclustr(:)
      integer                          :: clsize
      integer                          :: ir1_,ir2_
      integer                          :: n,npr,npc,info,m,nz,lwork,liwork,m1,n1
      integer                          :: ctxt,rp,cp,desc(9),anb,nps,nn,n0,m0
      integer                          :: Merr
      integer                          :: iceil,numroc,pjlaenv
# ifdef noMSYNC
      integer(mpi_address_kind)        :: Maddress_dummy
# endif
      wrt  = .false. ; if(present(lwrite)) wrt = lwrite
      n    = size(eval)
      ir1_ = 1          ; if(present(ir1)) ir1_ = ir1
      ir2_ = n          ; if(present(ir2)) ir2_ = ir2
      r1_  = -huge(0d0) ; if(present(r1))  r1_  = r1
      r2_  =  huge(0d0) ; if(present(r2))  r2_  = r2      
      m    = ir2_ - ir1_ + 1
      mod  = 'N' ; if(present(evec)) mod = 'V'
      cmod = 'A'
      if(present(ir1).or.present(ir2)) cmod = 'I'
      if(present(r1) .or.present(r2)) then
        if(r1<r2) then
          if(cmod=='I'.and.ir1<=ir2) Bug('Valid index and value ranges given at the same time.')
          cmod = 'V'
        endif
      endif
      if(n==0) Error('zero dimension in eigenvalue problem.')
      if(present(mat)) then
        if(present(matp)) Bug('Optional arguments inconsistent.')        
        call check_dim(shape(mat), [n,n],'mat',__LINE__,__FILE__)
      else
        if(present(mat)) Bug('Optional arguments inconsistent.')
        nn = n*(n+1)/2
        call check_dim(shape(matp), [nn],'matp',__LINE__,__FILE__)
      endif
      if(cmod=='I') then
        if(ir1_<1)    Error('non-positive lower index.')
        if(ir1_>ir2_) Error('negative index range.')
        if(ir2_>n)    Error('index of selected eigenvalues exceeds maximal number.')
        if(ir1_==1.and.ir2_==n) cmod = 'A'
        if(present(evec)) call check_dim(shape(evec),[n,-m],'evec',__LINE__,__FILE__)
      else
        if(present(evec)) call check_dim(shape(evec),[n,n],'evec',__LINE__,__FILE__)
      endif
      ctxt = Mcomm
      call grid_proc(npr,npc,Msize,n,n)
      call blacs_gridinit(ctxt,'R',npr,npc)
      call blacs_gridinfo(ctxt,npr,npc,rp,cp)
      if(Mrank==0.and.wrt) write(6,'(A'NoA) ' ('//trim(chr(npr))//'*'//trim(chr(npc))//') '
      Rif(rp==-1) Bug('Root process not in BLACS.')
      beginSplit(-min(0,rp))
      if(rp/=-1) then
        m1     = numroc(n,blk,rp,0,npr)
        n1     = numroc(n,blk,cp,0,npc)
        nn     = max(n,blk,2)
        m0     = numroc(nn,blk,0,0,npr)
        n0     = numroc(nn,blk,0,0,npc)
        liwork = 6 * max(n,npr*npc+1,4)
        abstol = 2 * pdlamch(ctxt,'S')
        allocate ( iwork(liwork),ifail(n),iclustr(2*npr*npc),gap(npr*npc) )
        allocate ( mat1(m1,n1) )
        if(present(evec)) then ; allocate( evec1(m1,n1) )
        else                   ; allocate( evec1(1,1)   )
        endif
        clsize = CLUSTERSIZE
 1      continue
# if DZtyp == 2
        lwork  = n + ( m0 + n0 + blk ) * blk
        lrwork = 4*n + max( 5*nn , m0*n0 ) + iceil(m,npr*npc)*nn + (clsize-1)*n
# else
        lwork  = 5*n + max( 5*nn , m0*n0 + 2*blk*blk ) + iceil(n,npr*npc)*nn + (clsize-1)*n
# endif        
# ifdef HighAcc
        anb    = pjlaenv(ctxt,3,'P'//DZstr//'TTRD','L',0,0,0,0)
        nps    = max( numroc(n,1,0,0,int(sqrt(dble(npr*npc)))) , 2*anb )
        lwork  = max( lwork , n + 2*(anb+1)*(4*nps+2)+(nps+4)*nps , 2*m0*blk + n0*blk + blk*blk )
# endif
# if DZtyp == 2
        if(allocated(rwork)) deallocate(rwork)
        allocate ( rwork(max(3,lrwork)) )
# endif
        if(allocated(work)) deallocate(work)
        allocate ( work(max(3,lwork)) )
        call descinit(desc,n,n,blk,blk,0,0,ctxt,m1,info) ; if(info/=0) Error('descinit failed.')                
        if(present(mat)) then
          call cyclic_ DZ  (mat1,m1,n1, mat,n,n,rp,cp,npr,npc,blk,0)
        else
          call cyclic_ DZ p(mat1,m1,n1, matp,n,rp,cp,npr,npc,blk,0)
        endif
        call p DZsub evx(mod,cmod,'U',n,mat1,1,1,desc,r1_,r2_,ir1_,ir2_,abstol,m,nz,eval,orfac,
     &                   evec1,1,1,desc,work,lwork, RWORK iwork,liwork,ifail,iclustr,gap,info)
        if(info/=0) then
          if(clsize<n-1) then
            clsize = min(2*clsize,n-1)
            RWarn('P'//DZstr//'EVX failed with info: '//Chr(info)//'. Retry with larger CLUSTERSIZE: '//Chr(clsize))
            goto 1
          endif
          Error('P'//DZstr//'EVX failed with info: '//Chr(info))
        endif
        if(cmod=='I') then
          if(m /=ir2_-ir1_+1) Error('Not enough eigenvalues found.')
          if(nz/=ir2_-ir1_+1) Error('Not enough eigenvectors found.')
        endif
        if(present(ir2)) ir2 = m
      endif
      if(present(evec)) then
        if(present(win_evec)) then ! evec is in shared memory
          Nfence(evec)
          ifO evec = 0
          Nfence(evec)
        else                       ! evec is in distributed memory
          evec = 0 
        endif
      endif
      if(rp/=-1) then
        if(present(evec)) call cyclic_ DZ (evec1,m1,n1,evec,n,m,rp,cp,npr,npc,blk,1)
        call blacs_gridexit(ctxt)
        deallocate ( mat1,evec1 )
        deallocate ( work,iwork,ifail,iclustr,gap )
# if DZtyp == 2
        deallocate ( rwork )
# endif        
      endif
      endSplit
      if(present(ir2)) call Mcast(ir2)
      call Mcast(eval)
      end subroutine Mdiagonalize0_ DZ x

c     --------

c     Selected solutions to generalized eigenvalue problem
      subroutine Mdiagonalize0_ DZ ox(eval,evec,win_evec,mat,olap,matp,olapp,ir1,ir2,r1,r2,lwrite)
      use global,   only: Mrank,Msize,Mcomm,Nrank
      use util,     only: chr
      use Mwrapper, only: Mcast
      !use timer_util
      implicit none
      include 'mpif.h'
      real_dp, intent(out)             :: eval(:)
      Datatyp, intent(out),   optional :: evec(:,:)
      Datatyp, intent(in),    optional :: mat(:,:),olap(:,:),matp(:),olapp(:)
      real_dp, intent(in),    optional :: r1,r2
      integer, intent(in),    optional :: ir1,win_evec
      integer, intent(inout), optional :: ir2
      logical, intent(in),    optional :: lwrite
      logical                          :: wrt
# if DZtyp == 2
      real_dp, allocatable             :: rwork(:)
      integer                          :: lrwork,mm
# endif
      Datatyp, allocatable             :: mat1(:,:),olap1(:,:),evec1(:,:),work(:)
      real_dp, allocatable             :: gap(:)
      real_dp                          :: abstol,pdlamch
      real_dp                          :: r1_,r2_,orfac=ORFAC
      character                        :: mod,cmod
      integer, allocatable             :: iwork(:),ifail(:),iclustr(:)
      integer                          :: clsize
      integer                          :: ir1_,ir2_
      integer                          :: n,npr,npc,info,m,nz,lwork,liwork,m1,n1
      integer                          :: ctxt,rp,cp,desc(9),anb,nps,nn,n0,m0
      integer                          :: Merr
      integer                          :: iceil,numroc,pjlaenv
# ifdef noMSYNC
      integer(mpi_address_kind)        :: Maddress_dummy
# endif      
      !call timer_start('ScaLAPACK diag total')
      wrt  = .false. ; if(present(lwrite)) wrt = lwrite
      n    = size(eval)
      ir1_ = 1          ; if(present(ir1)) ir1_ = ir1
      ir2_ = n          ; if(present(ir2)) ir2_ = ir2
      r1_  = -huge(0d0) ; if(present(r1))  r1_  = r1
      r2_  =  huge(0d0) ; if(present(r2))  r2_  = r2      
      m    = ir2_ - ir1_ + 1
      mod  = 'N' ; if(present(evec)) mod = 'V'
      cmod = 'A'
      if(present(ir1).or.present(ir2)) cmod = 'I'
      if(present(r1) .or.present(r2)) then
        if(r1<r2) then
          if(cmod=='I'.and.ir1<=ir2) Bug('Valid index and value ranges given at the same time.')
          cmod = 'V'
        endif
      endif
      if(n==0) Error('zero dimension in eigenvalue problem.')
      if(present(mat)) then
        if(.not.present(olap).or.present(matp).or.present(olapp)) Bug('Optional arguments inconsistent.')        
        call check_dim(shape(mat), [n,n],'mat',__LINE__,__FILE__)
        call check_dim(shape(olap),[n,n],'olap',__LINE__,__FILE__)
      else
        if(.not.present(olapp).or.present(mat).or.present(olap)) Bug('Optional arguments inconsistent.')
        nn = n*(n+1)/2
        call check_dim(shape(matp), [nn],'matp',__LINE__,__FILE__)
        call check_dim(shape(olapp),[nn],'olapp',__LINE__,__FILE__)
      endif
      if(cmod=='I') then
        if(ir1_<1)    Error('non-positive lower index.')
        if(ir1_>ir2_) Error('negative index range.')
        if(ir2_>n)    Error('index of selected eigenvalues exceeds maximal number.')
        if(ir1_==1.and.ir2_==n) cmod = 'A'
        if(present(evec)) call check_dim(shape(evec),[n,-m],'evec',__LINE__,__FILE__)
      else
        if(present(evec)) call check_dim(shape(evec),[n,n],'evec',__LINE__,__FILE__)
      endif
      ctxt = Mcomm
      call grid_proc(npr,npc,Msize,n,n)
      call blacs_gridinit(ctxt,'R',npr,npc)
      call blacs_gridinfo(ctxt,npr,npc,rp,cp)
      if(Mrank==0.and.wrt) write(6,'(A'NoA) ' ('//trim(chr(npr))//'*'//trim(chr(npc))//') '
      Rif(rp==-1) Bug('Root process not in BLACS.')
      beginSplit(-min(0,rp))
      if(rp/=-1) then
        m1     = numroc(n,blk,rp,0,npr)
        n1     = numroc(n,blk,cp,0,npc)
        nn     = max(n,blk,2)
        m0     = numroc(nn,blk,0,0,npr)
        n0     = numroc(nn,blk,0,0,npc)
        liwork = 6 * max(n,npr*npc+1,4)
        abstol = 2 * pdlamch(ctxt,'S')
        allocate ( ifail(n),iclustr(2*npr*npc),gap(npr*npc) )
        Allocate_( iwork,(liwork) ) ; Wrtmem( iwork = 0 )
        Allocate_( mat1,(m1,n1) )   ; Wrtmem( mat1  = 0 )
        Allocate_( olap1,(m1,n1) )  ; Wrtmem( olap1 = 0 )
c        allocate ( mat1(m1,n1),olap1(m1,n1) )
        if(present(evec)) then ; Allocate_( evec1,(m1,n1) ) ; Wrtmem( evec1 = 0 )
        else                   ; Allocate_( evec1,(1,1)   )
        endif
        clsize = CLUSTERSIZE
 1      continue
# if DZtyp == 2
        lwork  = n + ( m0 + n0 + blk ) * blk
# else
        lwork  = 5*n + max( 5*nn , m0*n0 + 2*blk*blk ) + iceil(n,npr*npc)*nn + (clsize-1)*n
# endif
# ifdef HighAcc
        anb    = pjlaenv(ctxt,3,'P'//DZstr//'TTRD','L',0,0,0,0)
        nps    = max( numroc(n,1,0,0,int(sqrt(dble(npr*npc)))) , 2*anb )
        lwork  = max( lwork , n + 2*(anb+1)*(4*nps+2)+(nps+4)*nps , 2*m0*blk + n0*blk + blk*blk )
# endif
# if DZtyp == 2
        m      = n ; if(cmod=='I') m = ir2 - ir1 + 1
        nn     = max(n,blk,2)
        mm     = max(m,blk,2)
        m0     = numroc(nn,blk,0,0,npr)
        n0     = numroc(mm,blk,0,0,npc)
        lrwork = 4*n + max( 5*nn , m0*n0 ) + iceil(m,npr*npc)*nn + (clsize-1)*n
        if(allocated(rwork)) tDeallocate(rwork)
        Allocate_ ( rwork,(max(3,lrwork)) ) ; Wrtmem( rwork = 0 )
# endif
        if(allocated(work)) tDeallocate(work)
        Allocate_ ( work,(max(3,lwork)) )   ; Wrtmem( work = 0 )
        call descinit(desc,n,n,blk,blk,0,0,ctxt,m1,info) ; if(info/=0) Error('descinit failed.')
        if(present(mat)) then
          call cyclic_ DZ ( mat1,m1,n1, mat,n,n,rp,cp,npr,npc,blk,0)
          call cyclic_ DZ (olap1,m1,n1,olap,n,n,rp,cp,npr,npc,blk,0)
        else
          call cyclic_ DZ p( mat1,m1,n1, matp,n,rp,cp,npr,npc,blk,0)
          call cyclic_ DZ p(olap1,m1,n1,olapp,n,rp,cp,npr,npc,blk,0)
        endif
        !call timer_start('ScaLAPACK diag p*gvx call')
        call p DZsub gvx(1,mod,cmod,'U',n,mat1,1,1,desc,olap1,1,1,desc,r1_,r2_,ir1_,ir2_,abstol,m,nz,eval,orfac,
     &                   evec1,1,1,desc,work,lwork, RWORK iwork,liwork,ifail,iclustr,gap,info)
        !call timer_stop('ScaLAPACK diag p*gvx call')
        if(info/=0) then
          if(clsize<n-1) then
            clsize = min(2*clsize,n-1)
            RWarn('P'//DZstr//'GVX failed with info: '//Chr(info)//'. Retry with larger CLUSTERSIZE: '//Chr(clsize))
            goto 1
          endif
          Error('P'//DZstr//'GVX failed with info: '//Chr(info))
        endif
        if(cmod=='I') then
          if(m /=ir2_-ir1_+1) Error('Not enough eigenvalues found.')
          if(nz/=ir2_-ir1_+1) Error('Not enough eigenvectors found.')
        endif
        if(present(ir2)) ir2 = m
      endif
      if(present(evec)) then
        if(present(win_evec)) then ! evec is in shared memory
          Nfence(evec)
          ifO evec = 0
          Nfence(evec)
        else                       ! evec is in distributed memory
          evec = 0 
        endif
      endif
      if(rp/=-1) then
        if(present(evec)) call cyclic_ DZ (evec1,m1,n1,evec,n,m,rp,cp,npr,npc,blk,1)
        call blacs_gridexit(ctxt)
        Deallocate_(mat1)
        Deallocate_(olap1)
        Deallocate_(evec1)
c        deallocate ( mat1,olap1,evec1 )
        Deallocate_ ( iwork )
        Deallocate_ ( work )
        deallocate ( ifail,iclustr,gap )
# if DZtyp == 2
        Deallocate_ ( rwork )
# endif
      endif
      endSplit
      if(present(ir2)) call Mcast(ir2)
      call Mcast(eval)
      !call timer_stop('ScaLAPACK diag total')
      end subroutine Mdiagonalize0_ DZ ox

c     --------

c
c     Mdiagonalize... : Frontend functions
c

c     --------      

      subroutine Mdiagonalize_ DZ (eval,mat,ir1,ir2,r1,r2,lwrite)
      implicit none
      real_dp, intent(out)             :: eval(:)      
      Datatyp, intent(in)              :: mat(:,:)
      real_dp, intent(in),    optional :: r1,r2
      integer, intent(in),    optional :: ir1
      integer, intent(inout), optional :: ir2
      logical, intent(in),    optional :: lwrite
      logical                          :: wrt
      wrt = .false.
      if(present(lwrite)) wrt = lwrite
      select case( which_ranges(present(ir1),present(ir2),present(r1),present(r2)) )
        case(0) ; call Mdiagonalize0_ DZ  (eval,mat=mat,                            lwrite=wrt)
        case(1) ; call Mdiagonalize0_ DZ x(eval,mat=mat,ir1=ir1,ir2=ir2,            lwrite=wrt)
        case(2) ; call Mdiagonalize0_ DZ x(eval,mat=mat,                r1=r1,r2=r2,lwrite=wrt)        
        case(3) ; call Mdiagonalize0_ DZ x(eval,mat=mat,ir1=ir1,ir2=ir2,r1=r1,r2=r2,lwrite=wrt)
      end select
      end subroutine Mdiagonalize_ DZ

c     --------

      subroutine Mdiagonalize_ DZ p(eval,mat,ir1,ir2,r1,r2,lwrite)
      implicit none
      real_dp, intent(out)             :: eval(:)
      Datatyp, intent(in)              :: mat(:)
      real_dp, intent(in),    optional :: r1,r2
      integer, intent(in),    optional :: ir1
      integer, intent(inout), optional :: ir2
      logical, intent(in),    optional :: lwrite
      logical                          :: wrt
      wrt = .false.
      if(present(lwrite)) wrt = lwrite
      select case( which_ranges(present(ir1),present(ir2),present(r1),present(r2)) )
        case(0) ; call Mdiagonalize0_ DZ  (eval,matp=mat,                            lwrite=wrt)
        case(1) ; call Mdiagonalize0_ DZ x(eval,matp=mat,ir1=ir1,ir2=ir2,            lwrite=wrt)
        case(2) ; call Mdiagonalize0_ DZ x(eval,matp=mat,                r1=r1,r2=r2,lwrite=wrt)        
        case(3) ; call Mdiagonalize0_ DZ x(eval,matp=mat,ir1=ir1,ir2=ir2,r1=r1,r2=r2,lwrite=wrt)
      end select
      end subroutine Mdiagonalize_ DZ p

c     --------

      subroutine Mdiagonalize_ DZ v(evec,win_evec,eval,mat,ir1,ir2,r1,r2,lwrite)
      implicit none
      real_dp, intent(out)             :: eval(:)
      Datatyp, intent(out)             :: evec(:,:)
      Datatyp, intent(in)              :: mat(:,:)
      real_dp, intent(in),    optional :: r1,r2
      integer, intent(in),    optional :: ir1,win_evec
      integer, intent(inout), optional :: ir2
      logical, intent(in),    optional :: lwrite
      logical                          :: wrt
      wrt = .false.   
      if(present(lwrite)) wrt = lwrite
      select case( which_ranges(present(ir1),present(ir2),present(r1),present(r2),present(win_evec)) )
        case(0) ; call Mdiagonalize0_ DZ  (eval,evec         ,mat=mat,                            lwrite=wrt)
        case(1) ; call Mdiagonalize0_ DZ x(eval,evec,         mat=mat,ir1=ir1,ir2=ir2,            lwrite=wrt)
        case(2) ; call Mdiagonalize0_ DZ x(eval,evec,         mat=mat,                r1=r1,r2=r2,lwrite=wrt)        
        case(3) ; call Mdiagonalize0_ DZ x(eval,evec,         mat=mat,ir1=ir1,ir2=ir2,r1=r1,r2=r2,lwrite=wrt)
        case(4) ; call Mdiagonalize0_ DZ  (eval,evec,win_evec,mat=mat,                            lwrite=wrt)
        case(5) ; call Mdiagonalize0_ DZ x(eval,evec,win_evec,mat=mat,ir1=ir1,ir2=ir2,            lwrite=wrt)
        case(6) ; call Mdiagonalize0_ DZ x(eval,evec,win_evec,mat=mat,                r1=r1,r2=r2,lwrite=wrt)        
        case(7) ; call Mdiagonalize0_ DZ x(eval,evec,win_evec,mat=mat,ir1=ir1,ir2=ir2,r1=r1,r2=r2,lwrite=wrt)
      end select
      end subroutine Mdiagonalize_ DZ v

c     --------

      subroutine Mdiagonalize_ DZ pv(evec,win_evec,eval,mat,ir1,ir2,r1,r2,lwrite)
      implicit none
      real_dp, intent(out)             :: eval(:)
      Datatyp, intent(out)             :: evec(:,:)
      Datatyp, intent(in)              :: mat(:)
      real_dp, intent(in),    optional :: r1,r2
      integer, intent(in),    optional :: ir1,win_evec
      integer, intent(inout), optional :: ir2
      logical, intent(in),    optional :: lwrite
      logical                          :: wrt
      wrt = .false.
      if(present(lwrite)) wrt = lwrite
      select case( which_ranges(present(ir1),present(ir2),present(r1),present(r2),present(win_evec)) )
        case(0) ; call Mdiagonalize0_ DZ  (eval,evec,         matp=mat,                            lwrite=wrt)
        case(1) ; call Mdiagonalize0_ DZ x(eval,evec,         matp=mat,ir1=ir1,ir2=ir2,            lwrite=wrt)
        case(2) ; call Mdiagonalize0_ DZ x(eval,evec,         matp=mat,                r1=r1,r2=r2,lwrite=wrt)        
        case(3) ; call Mdiagonalize0_ DZ x(eval,evec,         matp=mat,ir1=ir1,ir2=ir2,r1=r1,r2=r2,lwrite=wrt)
        case(4) ; call Mdiagonalize0_ DZ  (eval,evec,win_evec,matp=mat,                            lwrite=wrt)
        case(5) ; call Mdiagonalize0_ DZ x(eval,evec,win_evec,matp=mat,ir1=ir1,ir2=ir2,            lwrite=wrt)
        case(6) ; call Mdiagonalize0_ DZ x(eval,evec,win_evec,matp=mat,                r1=r1,r2=r2,lwrite=wrt)        
        case(7) ; call Mdiagonalize0_ DZ x(eval,evec,win_evec,matp=mat,ir1=ir1,ir2=ir2,r1=r1,r2=r2,lwrite=wrt)
      end select
      end subroutine Mdiagonalize_ DZ pv

c     --------

      subroutine Mdiagonalize_ DZ o(eval,mat,olap,ir1,ir2,r1,r2,lwrite)
      implicit none
      real_dp, intent(out)             :: eval(:)
      Datatyp, intent(in)              :: mat(:,:),olap(:,:)
      real_dp, intent(in),    optional :: r1,r2
      integer, intent(in),    optional :: ir1
      integer, intent(inout), optional :: ir2
      logical, intent(in),    optional :: lwrite
      logical                          :: wrt
      wrt = .false.
      if(present(lwrite)) wrt = lwrite
      select case( which_ranges(present(ir1),present(ir2),present(r1),present(r2)) )
        case(0) ; call Mdiagonalize0_ DZ ox(eval,mat=mat,olap=olap,                            lwrite=wrt)
        case(1) ; call Mdiagonalize0_ DZ ox(eval,mat=mat,olap=olap,ir1=ir1,ir2=ir2,            lwrite=wrt)
        case(2) ; call Mdiagonalize0_ DZ ox(eval,mat=mat,olap=olap,                r1=r1,r2=r2,lwrite=wrt)        
        case(3) ; call Mdiagonalize0_ DZ ox(eval,mat=mat,olap=olap,ir1=ir1,ir2=ir2,r1=r1,r2=r2,lwrite=wrt)
      end select
      end subroutine Mdiagonalize_ DZ o

c     --------

      subroutine Mdiagonalize_ DZ po(eval,mat,olap,ir1,ir2,r1,r2,lwrite)
      implicit none
      real_dp, intent(out)             :: eval(:)
      Datatyp, intent(in)              :: mat(:),olap(:)
      real_dp, intent(in),    optional :: r1,r2
      integer, intent(in),    optional :: ir1
      integer, intent(inout), optional :: ir2
      logical, intent(in),    optional :: lwrite
      logical                          :: wrt
      wrt = .false.
      if(present(lwrite)) wrt = lwrite
      select case( which_ranges(present(ir1),present(ir2),present(r1),present(r2)) )
        case(0) ; call Mdiagonalize0_ DZ ox(eval,matp=mat,olapp=olap,                            lwrite=wrt)
        case(1) ; call Mdiagonalize0_ DZ ox(eval,matp=mat,olapp=olap,ir1=ir1,ir2=ir2,            lwrite=wrt)
        case(2) ; call Mdiagonalize0_ DZ ox(eval,matp=mat,olapp=olap,                r1=r1,r2=r2,lwrite=wrt)        
        case(3) ; call Mdiagonalize0_ DZ ox(eval,matp=mat,olapp=olap,ir1=ir1,ir2=ir2,r1=r1,r2=r2,lwrite=wrt)
      end select
      end subroutine Mdiagonalize_ DZ po

c     --------

      subroutine Mdiagonalize_ DZ vo(evec,win_evec,eval,mat,olap,ir1,ir2,r1,r2,lwrite)
      implicit none
      real_dp, intent(out)             :: eval(:)
      Datatyp, intent(out)             :: evec(:,:)
      Datatyp, intent(in)              :: mat(:,:),olap(:,:)
      real_dp, intent(in),    optional :: r1,r2
      integer, intent(in),    optional :: ir1,win_evec
      integer, intent(inout), optional :: ir2
      logical, intent(in),    optional :: lwrite
      logical                          :: wrt
      wrt = .false.
      if(present(lwrite)) wrt = lwrite
      select case( which_ranges(present(ir1),present(ir2),present(r1),present(r2),present(win_evec)) )
        case(0) ; call Mdiagonalize0_ DZ ox(eval,evec,         mat=mat,olap=olap,                            lwrite=wrt)
        case(1) ; call Mdiagonalize0_ DZ ox(eval,evec,         mat=mat,olap=olap,ir1=ir1,ir2=ir2,            lwrite=wrt)
        case(2) ; call Mdiagonalize0_ DZ ox(eval,evec,         mat=mat,olap=olap,                r1=r1,r2=r2,lwrite=wrt)        
        case(3) ; call Mdiagonalize0_ DZ ox(eval,evec,         mat=mat,olap=olap,ir1=ir1,ir2=ir2,r1=r1,r2=r2,lwrite=wrt)
        case(4) ; call Mdiagonalize0_ DZ ox(eval,evec,win_evec,mat=mat,olap=olap,                            lwrite=wrt)
        case(5) ; call Mdiagonalize0_ DZ ox(eval,evec,win_evec,mat=mat,olap=olap,ir1=ir1,ir2=ir2,            lwrite=wrt)
        case(6) ; call Mdiagonalize0_ DZ ox(eval,evec,win_evec,mat=mat,olap=olap,                r1=r1,r2=r2,lwrite=wrt)        
        case(7) ; call Mdiagonalize0_ DZ ox(eval,evec,win_evec,mat=mat,olap=olap,ir1=ir1,ir2=ir2,r1=r1,r2=r2,lwrite=wrt)
      end select
      end subroutine Mdiagonalize_ DZ vo

c     --------

      subroutine Mdiagonalize_ DZ pvo(evec,win_evec,eval,mat,olap,ir1,ir2,r1,r2,lwrite)
      implicit none
      real_dp, intent(out)             :: eval(:)
      Datatyp, intent(out)             :: evec(:,:)
      Datatyp, intent(in)              :: mat(:),olap(:)
      real_dp, intent(in),    optional :: r1,r2
      integer, intent(in),    optional :: ir1,win_evec
      integer, intent(inout), optional :: ir2
      logical, intent(in),    optional :: lwrite
      logical                          :: wrt
      wrt = .false.
      if(present(lwrite)) wrt = lwrite
      select case( which_ranges(present(ir1),present(ir2),present(r1),present(r2),present(win_evec)) )
        case(0) ; call Mdiagonalize0_ DZ ox(eval,evec,         matp=mat,olapp=olap,                            lwrite=wrt)
        case(1) ; call Mdiagonalize0_ DZ ox(eval,evec,         matp=mat,olapp=olap,ir1=ir1,ir2=ir2,            lwrite=wrt)
        case(2) ; call Mdiagonalize0_ DZ ox(eval,evec,         matp=mat,olapp=olap,                r1=r1,r2=r2,lwrite=wrt)        
        case(3) ; call Mdiagonalize0_ DZ ox(eval,evec,         matp=mat,olapp=olap,ir1=ir1,ir2=ir2,r1=r1,r2=r2,lwrite=wrt)
        case(4) ; call Mdiagonalize0_ DZ ox(eval,evec,win_evec,matp=mat,olapp=olap,                            lwrite=wrt)
        case(5) ; call Mdiagonalize0_ DZ ox(eval,evec,win_evec,matp=mat,olapp=olap,ir1=ir1,ir2=ir2,            lwrite=wrt)
        case(6) ; call Mdiagonalize0_ DZ ox(eval,evec,win_evec,matp=mat,olapp=olap,                r1=r1,r2=r2,lwrite=wrt)        
        case(7) ; call Mdiagonalize0_ DZ ox(eval,evec,win_evec,matp=mat,olapp=olap,ir1=ir1,ir2=ir2,r1=r1,r2=r2,lwrite=wrt)
      end select
      end subroutine Mdiagonalize_ DZ pvo

c     --------

c
c     Minverse0... : Backend functions
c

c     --------

      subroutine Minverse0_ DZ (mat,matp,win_mat)
      use global,   only: Mrank,Msize,Mcomm,Nrank
      use util,     only: chr
      use Mwrapper, only: Mcast
      !use timer_util
      implicit none
      include 'mpif.h'
      integer, intent(in),    optional :: win_mat
      Datatyp, intent(inout), optional :: mat(:,:),matp(:)
      Datatyp, allocatable             :: mat1(:,:),work(:)
      Datatyp                          :: cdum
      integer, allocatable             :: ipiv(:),iwork(:)
      integer                          :: n,nn,npr,npc,info,lwork,liwork,m1,n1
      integer                          :: ctxt,rp,cp,desc(9),Merr,win_matp
      integer                          :: numroc
# ifdef noMSYNC
      integer(mpi_address_kind)        :: Maddress_dummy
# endif
      !call timer_start('ScaLAPACK inv total')
      if(present(mat).eqv.present(matp)) Bug('Optional arguments inconsistent.')
      if(present(mat)) then
        n = size(mat,1) ; call check_dim(shape(mat),[n,n],'mat',__LINE__,__FILE__)
      else
        nn = size(matp) ; n = nint(sqrt(0.25d0+2*nn)-0.5d0) ; if(n*(n+1)-2*nn/=0) Bug('Input matrix has wrong size.')
      endif
      if(n==0) Bug('matrix to be inverted has rank 0')
      ctxt = Mcomm
      !call timer_start('ScaLAPACK inv grid_proc')
      call grid_proc(npr,npc,Msize,n,n)
      call blacs_gridinit(ctxt,'R',npr,npc)
      call blacs_gridinfo(ctxt,npr,npc,rp,cp)
      !call timer_stop('ScaLAPACK inv grid_proc')
      Rif(rp==-1) Bug('Root process not in BLACS.')
      beginSplit(-min(0,rp))
      if(rp/=-1) then
        m1 = numroc(n,blk,rp,0,npr)
        n1 = numroc(n,blk,cp,0,npc)
        allocate ( mat1(m1,n1) )
        !call timer_start('ScaLAPACK inv desinit')
        call descinit(desc,n,n,blk,blk,0,0,ctxt,m1,info)            ; if(info/=0) Error('descinit failed.')
        !call timer_stop ('ScaLAPACK inv desinit')
        if(present(mat)) then
          allocate ( ipiv(m1+blk) )
          !call timer_start('ScaLAPACK inv cyclic_distr')
          call cyclic_ DZ (mat1,m1,n1, mat,n,n,rp,cp,npr,npc,blk,0)
          !call timer_stop ('ScaLAPACK inv cyclic_distr')
          !call timer_start('ScaLAPACK inv p*getrf')
          call p DZ getrf(n,n,mat1,1,1,desc,ipiv,info)                 ; if(info/=0) Error('ScaLAPACK call failed.')
          !call timer_stop ('ScaLAPACK inv p*getrf')
          !call timer_start('ScaLAPACK inv p*getri')
          call p DZ getri(n,mat1,1,1,desc,ipiv,cdum,-1,liwork,-1,info) ; if(info/=0) Error('ScaLAPACK call failed (first call).')
          lwork = nint(dble(cdum))
          allocate ( work(lwork),iwork(liwork) )
          call p DZ getri(n,mat1,1,1,desc,ipiv,work,lwork,iwork,liwork,info)
          !call timer_stop ('ScaLAPACK inv p*getri')
          if(info/=0) Error('ScaLAPACK call failed (second call) with info: '//chr(info))
        else
          !call timer_start('ScaLAPACK inv cyclic')
          call cyclic_ DZ p(mat1,m1,n1, matp,n,rp,cp,npr,npc,blk,0)
          !call timer_stop ('ScaLAPACK inv cyclic')
          !call timer_start('ScaLAPACK inv p*potrf')
          call p DZ potrf('U',n,mat1,1,1,desc,info)                    ; if(info/=0) Error('ScaLAPACK call failed.')
          !call timer_stop ('ScaLAPACK inv p*potrf')
          !call timer_start('ScaLAPACK inv p*potri')
          call p DZ potri('U',n,mat1,1,1,desc,info)
          !call timer_stop ('ScaLAPACK inv p*potri')
          if(info/=0) Error('ScaLAPACK call failed with info: '//chr(info))
        endif
      endif
      !call timer_start('ScaLAPACK inv init')
      if(present(win_mat)) then ! mat (matp) is in shared memory
        if(present(mat)) then
          Nfence(mat)
          ifO mat = 0
          Nfence(mat)
        else
          win_matp = win_mat
          Nfence(matp)
          ifO matp = 0
          Nfence(matp)
        endif
      else                      ! mat (matp) is in distributed memory
        if(present(mat)) then ; mat  = 0
        else                  ; matp = 0
        endif
      endif
      !call timer_stop('ScaLAPACK inv init')      
      if(rp/=-1) then
        if(present(mat)) then
          !call timer_start('ScaLAPACK inv cyclic_compile')
          call cyclic_ DZ (mat1,m1,n1,mat,n,n,rp,cp,npr,npc,blk,1)
          !call timer_stop ('ScaLAPACK inv cyclic_compile')
          deallocate ( mat1,ipiv,work,iwork )
        else
          !call timer_start('ScaLAPACK inv cyclic')
          call cyclic_ DZ p(mat1,m1,n1,matp,n,rp,cp,npr,npc,blk,1)
          !call timer_stop ('ScaLAPACK inv cyclic')
          deallocate ( mat1 )
        endif
        !call timer_start('ScaLAPACK inv gridexit')
        call blacs_gridexit(ctxt)
        !call timer_stop ('ScaLAPACK inv gridexit')
      endif
      endSplit
      !call timer_stop('ScaLAPACK inv total')
      end subroutine Minverse0_ DZ

c     --------

c
c     Minverse... : Frontend functions
c

c     --------

      subroutine Minverse_ DZ (mat,win_mat)
      implicit none
      Datatyp, intent(inout)        :: mat(:,:)
      integer, intent(in), optional :: win_mat
      if(present(win_mat)) then ; call Minverse0_ DZ (mat=mat,win_mat=win_mat)
      else                      ; call Minverse0_ DZ (mat=mat)
      endif
      end subroutine Minverse_ DZ

c     --------

      subroutine Minverse_ DZ p(matp,win_matp)
      implicit none
      Datatyp, intent(inout)        :: matp(:)
      integer, intent(in), optional :: win_matp
      if(present(win_matp)) then ; call Minverse0_ DZ (matp=matp,win_mat=win_matp)
      else                       ; call Minverse0_ DZ (matp=matp)
      endif
      end subroutine Minverse_ DZ p
      
c     --------

c mode=0: Distributes matrix mat(:m,:n) onto processes in block cyclic form into mat1(:m1,:n1)
c mode=1: Compile process-bound data mat1(:m1,:n1) into global matrix mat(:m,:n)
c Inout: mat(:m,:n)    - global matrix
c        mat1(:m1,:n1) - local portion of matrix
c Input: nrow,ncol     - process grid nrow x ncol
c        row,col       - coordinates of local process
c        blk           - block size (e.g., 64)
      subroutine cyclic_ DZ (mat1,m1,n1,mat,m,n,row,col,nrow,ncol,blk,mode)
      implicit none
      integer, intent(in)  :: m,n,m1,n1,row,col,nrow,ncol,blk,mode
      Datatyp              :: mat1(m1,n1)
      Datatyp              :: mat(m,n)
      integer              :: ic,ir,icol,irow,blkr,blkc
      integer              :: i,j
      if(all(mode/=[0,1])) Bug('mode/=0,1')
      ic = 0
      do icol = col*blk,n-1,blk*ncol
        blkc = min(blk,n-icol)
        ir   = 0
        do irow = row*blk,m-1,blk*nrow
          blkr = min(blk,m-irow)
          if(mode==0) then ; mat1(ir+1:ir+blkr,ic+1:ic+blkc) = mat(irow+1:irow+blkr,icol+1:icol+blkc)
          else             ; mat(irow+1:irow+blkr,icol+1:icol+blkc) = mat1(ir+1:ir+blkr,ic+1:ic+blkc)
          endif          
          ir   = ir + blkr
        enddo
        if(ir/=m1) Bug('Wrong number of rows.')
        ic = ic + blkc
      enddo
      if(ic/=n1.and.n==m) Bug('Wrong number of columns.')
      end subroutine cyclic_ DZ

c     --------

c Same for packed symmetric matrix mat(:n*(n+1)/2)
      subroutine cyclic_ DZ p(mat1,m1,n1,mat,n,row,col,nrow,ncol,blk,mode)
      implicit none
      integer, intent(in)  :: n,m1,n1,row,col,nrow,ncol,blk,mode
      Datatyp              :: mat1(m1,n1)
      Datatyp              :: mat(n*(n+1)/2)
      integer              :: ic,ir,icol,irow,blkr,blkc
      integer              :: i,j,ii,jj
      if(all(mode/=[0,1])) Bug('mode/=0,1')
      ic = 0
      do icol = col*blk,n-1,blk*ncol
        blkc = min(blk,n-icol)
        ir   = 0
        do irow = row*blk,n-1,blk*nrow
          blkr = min(blk,n-irow)
          do j = 1,blkc
            do i = 1,blkr
              ii = irow + i
              jj = icol + j
              if(mode==0) then
                if(ii<=jj) then ; mat1(ir+i,ic+j) =        mat(jj*(jj-1)/2+ii)
                else            ; mat1(ir+i,ic+j) = Conjg( mat(ii*(ii-1)/2+jj) )
                endif
              else
                if(ii<=jj)        mat(jj*(jj-1)/2+ii) =    mat1(ir+i,ic+j)
              endif
            enddo
          enddo
          ir = ir + blkr
        enddo
        if(ir/=m1) Bug('Wrong number of rows.')
        ic = ic + blkc
      enddo     
      if(ic/=n1.and.n>0) Bug('Wrong number of columns.')
      end subroutine cyclic_ DZ p

c     --------

# else
#   error macro noSCALAPACK disabled

c     --------      

      subroutine Mdiagonalize_ DZ (eval,mat,ir1,ir2,r1,r2,lwrite)
      use global,   only: Mrank
      use Mwrapper, only: Mcast
      implicit none
      real_dp, intent(out)             :: eval(:)      
      Datatyp, intent(in)              :: mat(:,:)
      real_dp, intent(in),    optional :: r1,r2
      integer, intent(in),    optional :: ir1
      integer, intent(inout), optional :: ir2
      logical, intent(in),    optional :: lwrite
      integer                          :: i,j
      select case( which_ranges(present(ir1),present(ir2),present(r1),present(r2)) )
        case(0) ; Rcall diagonalize(eval,mat                )
        case(1) ; Bug('Not implemented: Rcall diagonalize(eval,mat,ir1,ir2,0d0,0d0)')
        case(2) ; Bug('Not implemented: Rcall diagonalize(eval,mat,i,j,    r1,r2  )')
        case(3) ; Bug('Not implemented: Rcall diagonalize(eval,mat,ir1,ir2,r1,r2  )')
      end select
      call Mcast(eval)
      if(present(ir2)) call Mcast(ir2)
      end subroutine Mdiagonalize_ DZ

c     --------

      subroutine Mdiagonalize_ DZ p(eval,mat,ir1,ir2,r1,r2,lwrite)
      use global,   only: Mrank
      use Mwrapper, only: Mcast
      implicit none
      real_dp, intent(out)             :: eval(:)
      Datatyp, intent(in)              :: mat(:)
      real_dp, intent(in),    optional :: r1,r2
      integer, intent(in),    optional :: ir1
      integer, intent(inout), optional :: ir2
      logical, intent(in),    optional :: lwrite
      integer                          :: i,j
      select case( which_ranges(present(ir1),present(ir2),present(r1),present(r2)) )
        case(0) ; Rcall diagonalize(eval,mat                )
        case(1) ; Bug('Not implemented: Rcall diagonalize(eval,mat,ir1,ir2,0d0,0d0)')
        case(2) ; Bug('Not implemented: Rcall diagonalize(eval,mat,i,j,    r1,r2  )')
        case(3) ; Bug('Not implemented: Rcall diagonalize(eval,mat,ir1,ir2,r1,r2  )')
      end select      
      call Mcast(eval)
      if(present(ir2)) call Mcast(ir2)
      end subroutine Mdiagonalize_ DZ p

c     --------

      subroutine Mdiagonalize_ DZ v(evec,eval,mat,ir1,ir2,r1,r2,lwrite)
      use global,   only: Mrank
      use Mwrapper, only: Mcast
      implicit none
      real_dp, intent(out)             :: eval(:)
      Datatyp, intent(out)             :: evec(:,:)
      Datatyp, intent(in)              :: mat(:,:)
      real_dp, intent(in),    optional :: r1,r2
      integer, intent(in),    optional :: ir1
      integer, intent(inout), optional :: ir2
      logical, intent(in),    optional :: lwrite
      integer                          :: i,j
      select case( which_ranges(present(ir1),present(ir2),present(r1),present(r2)) )
        case(0) ; Rcall diagonalize(evec,eval,mat                )
        case(1) ; Rcall diagonalize(evec,eval,mat,ir1,ir2,0d0,0d0)
        case(2) ; Rcall diagonalize(evec,eval,mat,i,j,    r1,r2  )
        case(3) ; Rcall diagonalize(evec,eval,mat,ir1,ir2,r1,r2  )
      end select
      call Mcast(eval)
      call Mcast(evec)
      if(present(ir2)) call Mcast(ir2)
      end subroutine Mdiagonalize_ DZ v

c     --------

      subroutine Mdiagonalize_ DZ pv(evec,eval,mat,ir1,ir2,r1,r2,lwrite)
      use global,   only: Mrank
      use Mwrapper, only: Mcast
      implicit none
      real_dp, intent(out)             :: eval(:)
      Datatyp, intent(out)             :: evec(:,:)
      Datatyp, intent(in)              :: mat(:)
      real_dp, intent(in),    optional :: r1,r2
      integer, intent(in),    optional :: ir1
      integer, intent(inout), optional :: ir2
      logical, intent(in),    optional :: lwrite
      integer                          :: i,j
      select case( which_ranges(present(ir1),present(ir2),present(r1),present(r2)) )
        case(0) ; Rcall diagonalize(evec,eval,mat                )
        case(1) ; Rcall diagonalize(evec,eval,mat,ir1,ir2,0d0,0d0)
        case(2) ; Rcall diagonalize(evec,eval,mat,i,j,    r1,r2  )
        case(3) ; Rcall diagonalize(evec,eval,mat,ir1,ir2,r1,r2  )
      end select
      call Mcast(eval)
      call Mcast(evec)
      if(present(ir2)) call Mcast(ir2)
      end subroutine Mdiagonalize_ DZ pv

c     --------

      subroutine Mdiagonalize_ DZ o(eval,mat,olap,ir1,ir2,r1,r2,lwrite)
      use global,   only: Mrank
      use Mwrapper, only: Mcast
      implicit none
      real_dp, intent(out)             :: eval(:)
      Datatyp, intent(in)              :: mat(:,:),olap(:,:)
      real_dp, intent(in),    optional :: r1,r2
      integer, intent(in),    optional :: ir1
      integer, intent(inout), optional :: ir2
      logical, intent(in),    optional :: lwrite
      integer                          :: i,j
      select case( which_ranges(present(ir1),present(ir2),present(r1),present(r2)) )
        case(0) ; Bug('Not implemented: Rcall diagonalize(eval,mat,olap                )')
        case(1) ; Bug('Not implemented: Rcall diagonalize(eval,mat,olap,ir1,ir2,0d0,0d0)')
        case(2) ; Bug('Not implemented: Rcall diagonalize(eval,mat,olap,i,j,    r1,r2  )')
        case(3) ; Bug('Not implemented: Rcall diagonalize(eval,mat,olap,ir1,ir2,r1,r2  )')
      end select
      call Mcast(eval)
      if(present(ir2)) call Mcast(ir2)
      end subroutine Mdiagonalize_ DZ o

c     --------

      subroutine Mdiagonalize_ DZ po(eval,mat,olap,ir1,ir2,r1,r2,lwrite)
      use global,   only: Mrank
      use Mwrapper, only: Mcast
      implicit none
      real_dp, intent(out)             :: eval(:)
      Datatyp, intent(in)              :: mat(:),olap(:)
      real_dp, intent(in),    optional :: r1,r2
      integer, intent(in),    optional :: ir1
      integer, intent(inout), optional :: ir2
      logical, intent(in),    optional :: lwrite
      integer                          :: i,j
      select case( which_ranges(present(ir1),present(ir2),present(r1),present(r2)) )
        case(0) ; Bug('Not implemented: Rcall diagonalize(eval,mat,olap                )')
        case(1) ; Bug('Not implemented: Rcall diagonalize(eval,mat,olap,ir1,ir2,0d0,0d0)')
        case(2) ; Bug('Not implemented: Rcall diagonalize(eval,mat,olap,i,j,    r1,r2  )')
        case(3) ; Bug('Not implemented: Rcall diagonalize(eval,mat,olap,ir1,ir2,r1,r2  )')
      end select
      call Mcast(eval)
      if(present(ir2)) call Mcast(ir2)
      end subroutine Mdiagonalize_ DZ po

c     --------

      subroutine Mdiagonalize_ DZ vo(evec,eval,mat,olap,ir1,ir2,r1,r2,lwrite)
      use global,   only: Mrank
      use Mwrapper, only: Mcast
      implicit none
      real_dp, intent(out)             :: eval(:)
      Datatyp, intent(out)             :: evec(:,:)
      Datatyp, intent(in)              :: mat(:,:),olap(:,:)
      real_dp, intent(in),    optional :: r1,r2
      integer, intent(in),    optional :: ir1
      integer, intent(inout), optional :: ir2
      logical, intent(in),    optional :: lwrite
      integer                          :: i,j
      select case( which_ranges(present(ir1),present(ir2),present(r1),present(r2)) )
        case(0) ; Rcall diagonalize(evec,eval,mat,olap                )
        case(1) ; Rcall diagonalize(evec,eval,mat,olap,ir1,ir2,0d0,0d0)
        case(2) ; Rcall diagonalize(evec,eval,mat,olap,i,j,    r1,r2  )
        case(3) ; Rcall diagonalize(evec,eval,mat,olap,ir1,ir2,r1,r2  )
      end select
      call Mcast(eval)
      call Mcast(evec)
      if(present(ir2)) call Mcast(ir2)
      end subroutine Mdiagonalize_ DZ vo

c     --------

      subroutine Mdiagonalize_ DZ pvo(evec,eval,mat,olap,ir1,ir2,r1,r2,lwrite)
      use global,   only: Mrank
      use Mwrapper, only: Mcast
      implicit none
      real_dp, intent(out)             :: eval(:)
      Datatyp, intent(out)             :: evec(:,:)
      Datatyp, intent(in)              :: mat(:),olap(:)
      real_dp, intent(in),    optional :: r1,r2
      integer, intent(in),    optional :: ir1
      integer, intent(inout), optional :: ir2
      logical, intent(in),    optional :: lwrite
      integer                          :: i,j
      select case( which_ranges(present(ir1),present(ir2),present(r1),present(r2)) )
        case(0) ; Rcall diagonalize(evec,eval,mat,olap                )
        case(1) ; Rcall diagonalize(evec,eval,mat,olap,ir1,ir2,0d0,0d0)
        case(2) ; Rcall diagonalize(evec,eval,mat,olap,i,j,    r1,r2  )
        case(3) ; Rcall diagonalize(evec,eval,mat,olap,ir1,ir2,r1,r2  )
      end select
      call Mcast(eval)
      call Mcast(evec)
      if(present(ir2)) call Mcast(ir2)
      end subroutine Mdiagonalize_ DZ pvo

c     --------

      subroutine Minverse_ DZ (mat)
      use global,   only: Mrank
      use Mwrapper, only: Mcast
      implicit none
      Datatyp, intent(inout) :: mat(:,:)
      Rcall inverse(mat)
      call Mcast(mat)
      end subroutine Minverse_ DZ

c     --------

      subroutine Minverse_ DZ p(matp)
      use global,   only: Mrank
      use Mwrapper, only: Mcast
      implicit none
      Datatyp, intent(inout) :: matp(:)
      Rcall inverse(matp)
      call Mcast(matp)
      end subroutine Minverse_ DZ p

# endif

c     --------
      
# undef DZ
# undef DZsub
# undef DZstr      
# undef Datatyp
# undef Conjg
# undef RWORK

# if DZtyp == 2
#   undef DZtyp
# endif      
