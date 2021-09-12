c--------------------------------------------------------------------------------
c Copyright (c) 2021 Peter Grünberg Institut, Forschungszentrum Jülich, Germany
c This file is part of SPEX and available as free software under the conditions
c of the MIT license as expressed in the LICENSE file in more detail.
c--------------------------------------------------------------------------------

# include "cppmacro.h"

c     ----------------

c     Orders rarr(1:n) according to size.
      subroutine rorder(rarr,n)
      use, intrinsic :: iso_fortran_env
      implicit none
      integer, intent(in)    :: n
      real_dp, intent(inout) :: rarr(1:n)
      real_dp                :: rhlp(1:n)
      integer                :: i,j,k
      rhlp = rarr
      do i=1,n
        rarr(i) = rhlp(i)
        do j=1,i-1
          if(rarr(j)>rhlp(i)) then
            do k=i,j+1,-1
              rarr(k) = rarr(k-1)
            enddo
            rarr(j) = rhlp(i)
            exit
          endif
        enddo
      enddo
      end

c     ---------------

c     Orders rarr(1:n) according to size and returns a correspondingly defined pointer in pnt.
      subroutine rorderp(pnt,rarr,n)
      use, intrinsic :: iso_fortran_env
      implicit none
      integer, intent(in)  :: n
      integer, intent(out) :: pnt(n)
      real_dp, intent(in)  :: rarr(n)
      integer              :: i,j,k
      do i = 1,n
        pnt(i) = i
        do j = 1,i-1
          if(rarr(pnt(j))>rarr(i)) then
            do k = i,j+1,-1
              pnt(k) = pnt(k-1)
            enddo
            pnt(j) = i
            exit
          endif
        enddo
      enddo
      end

c     ---------------

c     Same as rorder but uses a divide and conquer algorithm (much faster than rorderf).
      recursive subroutine rorderf(rarr,n)
      use, intrinsic :: iso_fortran_env
      implicit none
      integer, intent(in)    :: n
      real_dp, intent(inout) :: rarr(n)
      real_dp, allocatable   :: rarr0(:)
      real_dp                :: ravg,rtiny
      integer                :: n1,n2,i
      ! call rorder for small arrays
      if(n<=10) then
        call rorder(rarr,n)
        return
      endif
      ! divide
      allocate ( rarr0(n) )
      ravg  = sum(rarr)/n
      rarr0 = rarr
      n1    = 0
      n2    = n+1
      do i = 1,n
        if(rarr0(i)<=ravg) then
          n1       = n1 + 1
          rarr(n1) = rarr0(i)
        else
          n2       = n2 - 1
          rarr(n2) = rarr0(i)
        endif
      enddo
      deallocate ( rarr0 )
      n2 = n - n2 + 1
      ! rarr(:) numerically identical?
      if(n1==0.or.n2==0) then
        rtiny = sqrt(1d0*n) * 1d-13
        do i = 2,n
          if(abs(rarr(i)-rarr(1))>rtiny) Bug('Expected rarr(:) to be numerically identical.')
        enddo
        return
      endif
      ! conquer
      call rorderf(rarr,      n1)
      call rorderf(rarr(n1+1),n2)
      end

c     ---------------

c     Same as rorderp but uses a divide and conquer algorithm (much faster than rorderp).
      recursive subroutine rorderpf(pnt,rarr,m)
      use, intrinsic :: iso_fortran_env
      implicit none
      integer, intent(in)  :: m
      integer, intent(out) :: pnt(abs(m))
      real_dp, intent(in)  :: rarr(abs(m))
      real_dp, allocatable :: rarr1(:)
      real_dp              :: ravg,rtiny
      integer, allocatable :: pnt1(:)
      integer              :: n,n1,n2,i
      n = abs(m)      
      ! initialize on first call (the other calls use a negative m)
      if(m>0) then
        do i = 1,n
          pnt(i) = i
        enddo
      endif
      ! call rorderp for small arrays
      if(n<=10) then
        allocate ( pnt1(n) )
        call rorderp(pnt1,rarr,n)
        pnt = pnt(pnt1)
        deallocate ( pnt1 )
        return
      endif
      ! divide
      allocate ( rarr1(n),pnt1(n) )
      ravg  = sum(rarr)/n
      n1    = 0
      n2    = n+1
      do i = 1,n
        if(rarr(i)<=ravg) then
          n1        = n1 + 1
          rarr1(n1) = rarr(i)
          pnt1(n1)  = pnt(i)
        else
          n2        = n2 - 1
          rarr1(n2) = rarr(i)
          pnt1(n2)  = pnt(i)
        endif
      enddo
      n2 = n - n2 + 1
      ! rarr(:) numerically identical?
      if(n1==0.or.n2==0) then
        deallocate ( pnt1 )
        rtiny = sqrt(1d0*n) * 1d-13
        do i = 2,n
          if(abs(rarr(i)-rarr(1))>rtiny) Bug('Expected rarr(:) to be numerically identical.')
        enddo
        return
      endif
      ! conquer
      pnt = pnt1
      deallocate ( pnt1 )
      call rorderpf(pnt,      rarr1,      -n1)
      call rorderpf(pnt(n1+1),rarr1(n1+1),-n2)
      deallocate ( rarr1 )
      end

c     ---------------

c     Orders iarr(1:n) according to size.
      subroutine iorder(iarr,n)
      implicit none
      integer, intent(in)    :: n
      integer, intent(inout) :: iarr(1:n)
      integer                :: ihlp(1:n)
      integer                :: i,j,k
      ihlp = iarr
      do i=1,n
        iarr(i) = ihlp(i)
        do j=1,i-1
          if(iarr(j)>ihlp(i)) then
            do k=i,j+1,-1
              iarr(k) = iarr(k-1)
            enddo
            iarr(j) = ihlp(i)
            exit
          endif
        enddo
      enddo
      end

c     ---------------

c     Same as iorder but uses a divide and conquer algorithm (much faster than iorder).
      recursive subroutine iorderf(iarr,n)
      use, intrinsic :: iso_fortran_env
      implicit none
      integer, intent(in)    :: n
      integer, intent(inout) :: iarr(n)
      integer, allocatable   :: iarr0(:)
      integer                :: n1,n2,i
      real_dp                :: ravg
      ! call iorder for small arrays
      if(n<=10) then
        call iorder(iarr,n)
        return
      endif
      ! divide
      allocate ( iarr0(n) )
      ravg  = sum(1d0*iarr)/n
      iarr0 = iarr
      n1    = 0
      n2    = n+1
      do i = 1,n
        if(iarr0(i)<=ravg) then
          n1       = n1 + 1
          iarr(n1) = iarr0(i)
        else
          n2       = n2 - 1
          iarr(n2) = iarr0(i)
        endif
      enddo
      deallocate ( iarr0 )
      n2 = n - n2 + 1
      ! iarr all identical?
      if(n1==0.or.n2==0) then
        do i = 2,n          
          if(iarr(i)/=iarr(1)) Bug('Expected iarr to be all identical.')
        enddo
        return
      endif
      ! conquer
      call iorderf(iarr,      n1)
      call iorderf(iarr(n1+1),n2)
      end

c     ---------------

c     Orders rarr(:,:) hirarchically according to size from "right to left" and returns pointer in pnt:
c     rarr(:,n2) is ordered, rarr(:,n2-1) is ordered where rarr(:,n2) is identical, etc.
c     If lswitch = true, the dimensions of rarr(:,:) are switched.
      recursive subroutine rnorderpf(pnt,rarr,dim1,dim2,n1,n2,lswitch)
      use, intrinsic :: iso_fortran_env
      implicit none
      integer, intent(in)  :: dim1,dim2,n1,n2
      integer, intent(out) :: pnt(*)
      logical, intent(in)  :: lswitch
      real_dp, intent(in)  :: rarr(dim1,dim2)
      integer, allocatable :: pnt1(:)
      integer              :: j1,j2
      if(lswitch) then
        if(any([n1,n2]>[dim2,dim1])) Bug('Dimensions too small.')
        call rorderpf(pnt,rarr(n2,:),n1)
      else
        if(any([n1,n2]>[dim1,dim2])) Bug('Dimensions too small.')
        call rorderpf(pnt,rarr(:,n2),n1)
      endif
      if(n2>1) then
        j1 = 1
        do while(j1<n1)
          j2 = j1
          do while(j2<n1)
            if(lswitch) then ; if(abs(rarr(n2,pnt(j2+1))-rarr(n2,pnt(j2)))>1d-12) exit
            else             ; if(abs(rarr(pnt(j2+1),n2)-rarr(pnt(j2),n2))>1d-12) exit
            endif
            j2 = j2 + 1     
          enddo
          if(j2-j1>0) then
            allocate(pnt1(j1:j2))
            if(lswitch) then ; call rnorderpf(pnt1,rarr(:n2-1,pnt(j1:j2)),n2-1,j2-j1+1,j2-j1+1,n2-1,lswitch)
            else             ; call rnorderpf(pnt1,rarr(pnt(j1:j2),:n2-1),j2-j1+1,n2-1,j2-j1+1,n2-1,lswitch)
            endif
            pnt(j1:j2) = pnt( pnt1(j1:j2) - 1 + j1 )
            deallocate(pnt1)
          endif          
          j1 = j2 + 1
        enddo
      endif
      end

c     ---------------

      ! Returns degeneracy limit and spread of val(1:n), i.e.,
      ! it looks for biggest gap in val(1:n) and returns gap center in "limit" and ratio in "spread".
      ! prepare = .true.  : constructs auxiliary array of steps abs(val(i+1)-val(i)) internally
      ! prepare = .false. : input array is assumed to be the auxiliary step array 
      ! Array val(:) must be size-ordered (increasing or decreasing if prepare=.true., else increasing).
      recursive subroutine deg_limit(limit,spread,val,n,prepare)
      use, intrinsic :: iso_fortran_env
      implicit none
      logical, intent(in)  :: prepare
      integer, intent(in)  :: n
      real_dp, intent(in)  :: val(n)
      real_dp, intent(out) :: limit,spread
      real_dp, allocatable :: val1(:)
      real_dp              :: spread1
      integer              :: i

      if(prepare) then
        allocate(val1(n-1))
        do i = 1,n-1
          if( val(i)>=val(i+1) .neqv. val(1)>=val(n) ) Bug('Input array not size-ordered (prepare=.true.).')
          val1(i) = abs( val(i+1) - val(i) )
        enddo
        call rorderf(val1,n-1)
        call deg_limit(limit,spread,val1,n-1,.false.)
        deallocate(val1)
      else
        do i = 1,n-1
          if( val(i)<0        ) Bug('Negative value in input array.')
          if( val(i)>val(i+1) ) Bug('Input array not size-ordered (prepare=.false.).')
        enddo
        spread = 0
        limit  = 0
        do i = 1,n-1
          if(val(i)>0) then
            spread1 = val(i+1) / val(i)
            if(spread1>spread) then
              spread = spread1
              limit  = sqrt(val(i)*val(i+1))
            endif
          endif
        enddo
        if(limit==0) Bug('Automatic determination of degeneracy limit failed.')
      endif
      end

c     ---------------

      subroutine iadd(arr,val,nval,narr) ! not used
      implicit none
      integer, intent(in)    :: val,nval
      integer, intent(inout) :: arr(nval,*),narr
      integer                :: i
      if(all( [ (all(val/=arr(:,i)),i=1,narr) ] )) then
        narr        = narr + 1
        arr(:,narr) = val
      endif
      end

c     ---------------

      subroutine mergeintarr(a,b,n,na,dima) ! not used
      implicit none
      integer, intent(in)    :: n,dima
      integer, intent(inout) :: na
      integer, intent(out)   :: a(n,*)
      integer, intent(in)    :: b(:,:)
      integer                :: nb,i,j
      logical                :: found
      nb  = size(b,2)
      do i=1,nb
        found=.false.
        do j=1,na
          if(all(b(:,i)==a(:,j))) found=.true.
        enddo
        if(.not.found) then
          na=na+1
          if(na>dima) Error('Dimension of destination array exceeded.')
          a(:,na)=b(:,i)
        endif
      enddo
      end

c     ----------------

c     Returns index of integer value ival in integer array iarr(1:narr).
      function position(ival,iarr,narr)
      implicit none
      integer             :: position
      integer, intent(in) :: narr,iarr(narr),ival
      integer             :: i
      do i = 1,narr
        if(iarr(i)==ival) then
          position = i
          return
        endif
      enddo
      position = 0
      end

c     ----------------

c     Returns .true. if string contains a readable integer, otherwise .false. .
      function isinteger(string)
      implicit none
      logical                  :: isinteger
      character(*), intent(in) :: string
      integer                  :: i,i0
      i0        = verify(string,' ')
      isinteger = i0/=0
      if(isinteger) then
        if(string(i0:i0)=='-') i0 = i0 + 1
        isinteger = i0<=len_trim(string)
        isinteger = isinteger .and. verify(string(i0:len_trim(string)),'0123456789')==0
      endif
      end

c     Returns .true. if string does not contain a readable integer, otherwise .false. .
      function noninteger(string)
      use, intrinsic :: iso_fortran_env
      implicit none
      logical                  :: noninteger
      character(*), intent(in) :: string
      logical                  :: isinteger
      noninteger = .not.isinteger(string)
      end

c     Returns .true. if string contains a readable real number, otherwise .false. .
      function isreal(string)
      implicit none
      logical                  :: isreal
      character(*), intent(in) :: string
      character                :: c
      integer                  :: i,i0
      logical                  :: number
      i0     = verify(string,' ')
      isreal = i0/=0
      if(isreal) then
        if(string(i0:i0)=='-') i0 = i0 + 1
        isreal = i0<=len_trim(string)
        isreal = isreal .and. verify(string(i0:len_trim(string)),'0123456789.')==0
        isreal = isreal .and. index(string,'.')==index(string,'.',back=.true.)
      endif
      end

c     ------------------

c     Returns the real number converted from character string 'ch'.
c     Allows eV unit if key starts with '*'
c     (In case of error the program issues an error message and stops.
c     If key is not empty, the keyword key is reported in the error message.)
      function ch2r(ch,key)
      use global, only: hartree
      use, intrinsic :: iso_fortran_env
      implicit none
      real_dp                   :: ch2r
      character(*), intent(in)  :: ch,key
      integer                   :: ios,i
      logical                   :: isreal
      ios = 1
      i   = max(1, len_trim(ch) - 1 )
      if(key(:1)=='*'.and.ch(i:)=='eV') then
        if(isreal(ch(:i-1))) then
          read(ch(:i-1),*,iostat=ios) ch2r
          ch2r = ch2r / hartree
        endif
      else
        if(isreal(ch)) read(ch,*,iostat=ios) ch2r
      endif
      i = 1 ; if(key(:1)=='*') i = 2
      if(ios/=0) then
        write(0,'(A)') 'ch2r: Not a real number: '//trim(ch)
        if(key(i:)/=' ')
     &  write(0,'(A)') '      Check input file after keyword '//trim(key(i:))//'.'
        Error('Not a real number: '//trim(ch))
      endif
      end

c     Returns the integer number converted from character string 'ch'. (See above)
      function ch2i(ch,key)
      use, intrinsic :: iso_fortran_env
      implicit none
      integer                   :: ch2i
      character(*), intent(in)  :: ch,key
      integer                   :: ios
      real_dp                   :: r
      logical                   :: isinteger
      ios = 1
      if(isinteger(ch)) read(ch,*,iostat=ios) r
      if(ios/=0.or.nint(r)/=r.or.index(ch,',')/=0) then
        write(0,'(A)') 'ch2i: Not an integer: '//trim(ch)
        if(key/=' ')
     &  write(0,'(A)') '      Check input file  after keyword '//trim(key)//'.'
        Error('Not an integer: '//trim(ch))
      endif
      ch2i = nint(r)
      end

c     ----------------

      subroutine spex_error(ekind,file,line,msg)
      Mpi(  use global )
      Load( use hdf5 )
      implicit none
      character(*), intent(in) :: ekind,file,msg
      integer,      intent(in) :: line
# if defined(LOAD) || defined(MPI)
      logical                  :: lmpi
      integer                  :: err
# endif
      if(msg/=' ') write(0,'(A,I4.4,A)') 'SPEX-'//ekind//' ('//file//':',line,') '//msg
      Load( call h5close_f(err) )
# ifdef MPI
      call mpi_initialized(lmpi,err)
      if(lmpi) call mpi_abort(mpi_comm_world,max(1,line),err)
# endif
# ifdef F08STOP
      error stop
# elif defined(F77STOP)
      stop 1
# else
      call exit(1)
# endif
      end

c     ----------------

      subroutine finish
      use global
      use file, only: fclose
      Load( use hdf5 )
      implicit none
# if defined(LOAD) || defined(MPI)
      integer      :: err,Merr
# endif
      logical      :: ldum
# ifdef LOAD
      call h5close_f(err) ; if(err<0) Error('Could not close HDF5.')
# endif
      inquire(inp,opened=ldum)
      if(ldum) call fclose(inp)
      call dealloc_job
      call dealloc_final
# ifdef MPI
      Oif(Ocomm/=mpi_comm_null) call mpi_comm_free(Ocomm,Merr)
      if( Ncomm/=mpi_comm_null) then
        Nfence_(mem)
        Nfree(mem)
#   ifdef WRTMEM
        Nfence_(mem0)
        Nfence_(mem_prv0)
        Nfree(mem0)
        Nfree(mem_prv0)
#   endif
        call mpi_comm_free(Ncomm,Merr)
      endif
      call mpi_finalize(Merr)
# endif
      stop
      end

c     ----------------

      subroutine teststop(feature)
      implicit none
      character(*), intent(in) :: feature
      write(0,*)
      write(0,'(A)') 'TESTSTOP. This feature ('//feature//') has never been tested.'
      write(0,'(A)') '          Remove the test stop, rerun, and check carefully!'
      Error('Teststop feature: '//feature)
      end

c     ----------------

      subroutine disable(feature)
      implicit none
      character(*), intent(in) :: feature
      write(0,*)
      write(0,'(A)') 'DISABLED. This feature ('//feature//') has been disabled. Please report.'
      Error('Disabled feature: '//feature)
      end

c     ----------------

c     Replaces normal isnan (because some compilers don't know isnan)
# ifndef noARITHM
      function isnan(r)
      use, intrinsic :: ieee_arithmetic
      use, intrinsic :: iso_fortran_env
      implicit none
      logical             :: isnan
      real_dp, intent(in) :: r
      isnan = ieee_is_nan(r)
      end
# endif

c     ----------------

c     Returns degenerate subspace deg0..deg1 for iband,ikpt,ispin
      subroutine getdeg(deg0,deg1,iband,ikpt,ispin)
      use global, only: deg,nband
      use, intrinsic :: iso_fortran_env
      implicit none
      integer, intent(out) :: deg0,deg1
      integer, intent(in)  :: iband,ikpt,ispin
      if(iband<0.or.iband>nband(ikpt,ispin)) Bug('iband out of range.')
      deg0 = iband
      deg1 = deg(iband,ikpt,ispin)
      if(deg1<deg0) then
        deg0 = deg1
        deg1 = deg(deg0,ikpt,ispin)
      endif
      end

c     ----------------

c     The following routines need a modular environment.

c     ----------------      

      module util

      interface reallocate
      module procedure  reallocate_r1, reallocate_c1, reallocate_i1, reallocate_ch1,
     &                  reallocate_r2, reallocate_c2, reallocate_i2, reallocate_ch2,
     &                  reallocate_r3, reallocate_c3, reallocate_i3, reallocate_ch3,
     &                  reallocate_r4, reallocate_c4, reallocate_i4, reallocate_ch4,
     &                  reallocate_r5, reallocate_c5, reallocate_i5, reallocate_ch5
      end interface

      interface chr
      module procedure chr_s,chr_r,chr_i,chr_rn,chr_in
      end interface

c     ----------------

      contains

# include "util.inc"      
c# define ARGTYPE 1
c# include "reallocate.inc"
c# undef  ARGTYPE
c# define ARGTYPE 2
c# include "reallocate.inc"
c# undef  ARGTYPE
c# define ARGTYPE 3
c# include "reallocate.inc"
c# undef  ARGTYPE
c# define ARGTYPE 4
c# include "reallocate.inc"
c# undef  ARGTYPE
c# define ARGTYPE 5
c# include "reallocate.inc"
c# undef  ARGTYPE
c# define ARGTYPE 6
c# include "reallocate.inc"
c# undef  ARGTYPE
c# define ARGTYPE 7
c# include "reallocate.inc"
c# undef  ARGTYPE
c# define ARGTYPE 8
c# include "reallocate.inc"
c# undef  ARGTYPE
c# define ARGTYPE 9
c# include "reallocate.inc"
c# undef  ARGTYPE
c# define ARGTYPE 10
c# include "reallocate.inc"
c# undef  ARGTYPE
c# define ARGTYPE 11
c# include "reallocate.inc"
c# undef  ARGTYPE
c# define ARGTYPE 12
c# include "reallocate.inc"
c# undef  ARGTYPE
c# define ARGTYPE 13
c# include "reallocate.inc"
c# undef  ARGTYPE
c# define ARGTYPE 14
c# include "reallocate.inc"
c# undef  ARGTYPE
c# define ARGTYPE 15
c# include "reallocate.inc"
c# undef  ARGTYPE
c# define ARGTYPE 16
c# include "reallocate.inc"
c# undef  ARGTYPE
c# define ARGTYPE 17
c# include "reallocate.inc"
c# undef  ARGTYPE
c# define ARGTYPE 18
c# include "reallocate.inc"
c# undef  ARGTYPE
c# define ARGTYPE 19
c# include "reallocate.inc"
c# undef  ARGTYPE
c# define ARGTYPE 20
c# include "reallocate.inc"
c# undef  ARGTYPE

c     ----------------
      
c     Parses string for a series of (positive) integers which are returned in iarr.
c     The string must have a syntax like '1,3-5,7'.
c     (In case of error the program issues an error message and stops.)
c     If key is not empty, the keyword key is reported in the error message.)
      subroutine str2iarr(iarr,string,key)
      implicit none
      integer,      allocatable :: iarr(:)
      character(*), intent(in)  :: string,key
      character                 :: c
      character(10)             :: help
      logical                   :: def
      integer                   :: lenstr,i,j,k,first
      integer                   :: ch2i
      if(allocated(iarr)) Bug('Array iarr already allocated.')
      lenstr = len_trim(string)
      def    = .false.
      do
        help   = ' '
        first  = -1
        k      = 0
        do i = 1,lenstr
          c = string(i:i)
          if(i==lenstr) then ; help = trim(help)//c ; c = ',' ; endif
          if(c==',') then
            if(first/=-1) then
              do j = first,ch2i(help,key)
                k = k + 1
                if(def) iarr(k) = j
              enddo
            else
              k = k + 1
              if(def) iarr(k) = ch2i(help,key)
            endif
            first = -1
            help  = ' '
            cycle
          else if(c=='-') then
            first = ch2i(help,key)
            help  = ' '
            cycle
          endif
          if(len_trim(help)==10) Error('Integer longer than ten digits: '//trim(string))
          help = trim(help)//c
        enddo
        if(first/=-1) then
          write(6,'(A'NoA) 'str2iarr: Could not parse '//trim(string)
          if(key/=' ') write(6,'(A)') ' after keyword '//trim(key)//'.'
          Error('Parse error.')
        endif
        if(def) exit
        def = .true.
        allocate(iarr(k))        
      enddo
      end subroutine str2iarr

c     ----------------

c     Divides string "A,B,..." into list items A, B, etc.
c     Returns end positions of each entry in pos(:): A=string(1:pos(1)), B=string(pos(1)+2:pos(2)), etc.
c     For simplicity, pos(0):=-1, so A=string(pos(0)+2:pos(1)).x
c     Comma is the default separator. Can be changed with the optional parameter sep.
c     If blank (' ') is the separator, multiple blanks between arguments are allowed.
c     The number of returned arguments is ubound(pos,1).      
      subroutine strlist(pos,string,sep0)
      implicit none
      integer, allocatable, intent(out) :: pos(:)
      character(*),         intent(in)  :: string
      character,  optional, intent(in)  :: sep0
      character                         :: sep
      character(:), allocatable         :: string1
      integer                           :: n,i,ind
      logical                           :: lwrd
      if(allocated(pos)) deallocate(pos)
      if(string==' ') then
        allocate(pos(0))
        return
      endif
      sep = ',' ; if(present(sep0)) sep = sep0
      do
        n    = 1
        lwrd = sep/=' '
        do i = 1,len_trim(string)
          if(string(i:i)==sep) then
            if(lwrd) then
              if(allocated(pos)) pos(n) = i - 1
              if(sep==' ') lwrd = .false.
              n = n + 1
            endif
          else
            lwrd = .true.
          endif
        enddo
        if(allocated(pos)) exit
        allocate(pos(0:n))
        pos(0) = -1
        pos(n) = len_trim(string)
      enddo
      end subroutine strlist

c     ----------------

      function chr_s(r,form)
      use, intrinsic :: iso_fortran_env
      implicit none
      character(40)                      :: chr_s
      real,         intent(in)           :: r
      character(*), intent(in), optional :: form
      character(10)                      :: form1
      if(present(form)) then ; form1 = form
      else                   ; form1 = 'F10.5'
      endif
      call enlarge_form(form1,dble(r))
      write(chr_s,'('//form1//')') r
      if(chr_s(:1)=='*') Warn('Output overflow (***) detected with form '//trim(form1))
      if(.not.present(form)) chr_s = adjustl(chr_s)
      end function chr_s

      function chr_r(r,form)
      use, intrinsic :: iso_fortran_env
      implicit none
      character(40)                      :: chr_r
      real_dp,      intent(in)           :: r
      character(*), intent(in), optional :: form
      character(10)                      :: form1      
      if(present(form)) then ; form1 = form
      else                   ; form1 = 'F10.5'
      endif
      call enlarge_form(form1,r)      
      write(chr_r,'('//form1//')') r
      if(.not.present(form)) chr_r = adjustl(chr_r)
      end function chr_r

      subroutine enlarge_form(form,r)
      use, intrinsic :: iso_fortran_env
      implicit none
      real_dp,      intent(in)    :: r
      character(*), intent(inout) :: form
      integer                     :: n,m1,m2,ind,ch2i
      if(r==0) return
      if(form(:1)=='F') then
        ind = index(form,'.')
        m1  = ch2i(form(2:ind-1),' ')
        m2  = ch2i(form( ind+1:),' ')
        n   = int(log(abs(r))/log(10d0)+1e-15) + 1 ; if(r<0) n = n + 1
        if(n>m1-m2-1) form = 'F'//trim(chr(n+m2+1))//form(ind:)
      endif
      end subroutine enlarge_form

      function chr_i(i,form)
      implicit none
      character(40)                      :: chr_i,hlp
      integer,      intent(in)           :: i
      character(*), intent(in), optional :: form
      write(hlp,'(I10)') i
      hlp = adjustl(hlp)
      if(present(form)) then
        write(chr_i,'('//form//')') i
        if(len_trim(hlp)>len_trim(chr_i)) chr_i = hlp
      else
        chr_i = hlp
      endif
      end function chr_i

      function chr_rn(r,sep,form)
      use, intrinsic :: iso_fortran_env
      implicit none
      real_dp,      intent(in)           :: r(:)
      character(*), intent(in), optional :: form,sep
      character(40*size(r))              :: chr_rn
      integer                            :: j,n
      chr_rn = ' '
      do j = 1,size(r)
        n = len_trim(chr_rn)
        if(j>1.and.present(sep)) then
          chr_rn(n+1:n+len(sep)) = sep
          n                      = n + len(sep)
        endif    
        if(present(form)) then ; chr_rn(n+1:) = adjustl(chr_r(r(j),form))
        else                   ; chr_rn(n+1:) = adjustl(chr_r(r(j)))
        endif
      enddo
      end function chr_rn

      function chr_in(i,sep,form)
      implicit none
      integer,      intent(in)           :: i(:)
      character(*), intent(in), optional :: form,sep
      character(40*size(i))              :: chr_in
      integer                            :: j,n
      chr_in = ' '
      do j = 1,size(i)
        n = len_trim(chr_in)
        if(j>1.and.present(sep)) then
          chr_in(n+1:n+len(sep)) = sep
          n                      = n + len(sep)
        endif
        if(present(form)) then ; chr_in(n+1:) = adjustl(chr_i(i(j),form))
        else                   ; chr_in(n+1:) = adjustl(chr_i(i(j)))
        endif
      enddo
      end function chr_in

c     ----------------

c     Returns true if states n1 and n2 (ikpt,ispin) are degenerate.
c     States are degenerate if their energies are different by less than edeg1. [If edeg1 is not present, deg(...) is used.]
      function same_eigenspace(n1,n2,ikpt,ispin,edeg1)
      use global, only: ene,deg
      use, intrinsic :: iso_fortran_env
      implicit none
      logical                       :: same_eigenspace
      real_dp, intent(in), optional :: edeg1
      integer, intent(in)           :: n1,n2,ikpt,ispin
      integer                       :: first1,first2
# ifdef switch_off_symmetry
      same_eigenspace = .false.
# else
      if(present(edeg1)) then
        same_eigenspace = abs(ene(n1,ikpt,ispin)-ene(n2,ikpt,ispin))<=edeg1
      else
        first1 = deg(n1,ikpt,ispin) ; if(first1>n1) first1 = n1
        first2 = deg(n2,ikpt,ispin) ; if(first2>n2) first2 = n2
        same_eigenspace = first1==first2
      endif
# endif
      end function same_eigenspace

c     ----------------

c
c     Average val(:) over degenerate states defined in eig(:).
      function average_deg(val,eig) result(out)
      use global, only: edeg
      use, intrinsic :: iso_fortran_env
      implicit none
      real_dp, intent(in) :: val(:),eig(:)
      real_dp             :: out(size(val))
      integer             :: i,j,n
      n = size(val) ; if(size(eig)/=n) Bug('Dimensions of val(:) and eig(:) differ.')      
      i = 1
      do while(i<=n)
        j = i
        do while(j<n)
          if(eig(j+1)-eig(i)>edeg) exit
          j = j + 1
        enddo
        out(i:j) = sum(val(i:j)) / (j-i+1)
        i        = j + 1
      enddo
      end function average_deg

c     ----------------

      end module util

c     ----------------

      subroutine cpu_done(cputime)
      use util, only: chr
      use, intrinsic :: iso_fortran_env
      implicit none
      real, intent(inout) :: cputime
      real                :: cputime1
      call cpu_time(cputime1)
c      write(6,'(A,F8.2,A)') 'done   ( Timing:',cputime1-cputime,' )'
      write(6,'(A)') 'done   ( Timing: '//Chf(cputime1-cputime,'F16.2')//' )'
      cputime = cputime1
      end

      subroutine cpu_stop(cputime)
      use, intrinsic :: iso_fortran_env
      implicit none
      real, intent(inout) :: cputime
      real                :: cputime1
      call cpu_time(cputime1)
      cputime = cputime1 - cputime
      end subroutine cpu_stop

c     ----------------

      function chr_band(n,k,s)
      use global, only: nspin1
      use util,   only: chr
      implicit none
      character(30)       :: chr_band
      integer, intent(in) :: n,k,s
      if(nspin1==2) then
        if(s==1) then ; chr_band = trim(chr(k))//':u('//trim(chr(n))//')'
        else          ; chr_band = trim(chr(k))//':d('//trim(chr(n))//')'
        endif
      else
        chr_band = trim(chr(k))//':('//trim(chr(n))//')'
      endif
      end

c     ----------------

c
c     Runs shell command (using execute_command_line or system)
      subroutine shell_command(command,ios)
      implicit none
      character(*), intent(in)  :: command
      integer,      intent(out) :: ios
# ifdef noEXECMD
      call system(command)
      ios = 0
# else      
      character(80)             :: err
      integer                   :: ext
      call execute_command_line(command,wait=.true.,exitstat=ext,cmdstat=ios,cmdmsg=err)
      if(ios/=0) then
        write(0,'(A)') 'execute_command_line returned error: '//trim(err)
        write(0,'(A)') 'Offending command: '//command
        if(err=='Invalid command line')
     &    write(0,'(A)') 'Note: This is probably a compiler bug. Try compiling Spex with "-DnoEXECMD"!'
      else
        ios = ext
      endif
# endif
      end

c     ----------------

c
c     Waits "time" seconds (CPU busy waiting).
      subroutine asleep(time)
      use, intrinsic :: iso_fortran_env
      implicit none
      real, intent(in) :: time
      integer_dp       :: walltime0,walltime,walldiff,wallrate,wallmax
      call system_clock(walltime0,wallrate,wallmax)
      if(wallrate*time>wallmax) Bug('asleep time argument exceeds wallmax.')      
      do
        call system_clock(walltime,wallrate,wallmax)
        walldiff = walltime - walltime0
        if(walldiff<0) walldiff = walldiff + wallmax
        if(walldiff>=wallrate*time) return
      enddo
      end
      
