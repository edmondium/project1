c--------------------------------------------------------------------------------
c Copyright (c) 2021 Peter Grünberg Institut, Forschungszentrum Jülich, Germany
c This file is part of SPEX and available as free software under the conditions
c of the MIT license as expressed in the LICENSE file in more detail.
c--------------------------------------------------------------------------------

# ifndef IDatatyp
#   define IDatatyp 1
#   define Datatyp real_dp
#   define SUFX    r
#   include __FILE__
#   define IDatatyp 2
#   define Datatyp complex_dp
#   define SUFX    c
#   include __FILE__
#   define IDatatyp 3
#   define Datatyp integer
#   define SUFX    i
#   include __FILE__
#   define IDatatyp 4
#   define Datatyp character(*)
#   define SUFX    ch
# endif

c     Permutes the dimensions dims=[dims(1),dims(2),...] of an array a(:,:...) according to
c     dims_new=[dims(perm(1)),dims(perm(2)),...] and also redistributes the array elements.

# if IDatatyp != 4
#   define BITS 4*8-1
      subroutine reshape_permute_ SUFX (array,dims,perm)
      use, intrinsic :: iso_fortran_env
      implicit none      
      Datatyp, intent(inout) :: array(*)
      integer, intent(in)    :: dims(:),perm(:)
      integer, allocatable   :: def(:)
      integer, parameter     :: bits = BITS ! number of bits used in integer array def
      Datatyp                :: buf,dum
      integer                :: i,j,m,n
      if(size(dims)/=size(perm)) Bug('Dimensions of dims and perm inconsistent.')
      n = size(dims)
      do i = 1,n
        if(all(perm/=i)) Bug('Incorrect permutation array.')
      enddo
      m   = product(dims) ; allocate(def(0:m/bits))
      def = 0
      i   = 1
      do
        do while(isdef(i,def))
          i = i + 1
          if(i>m) goto 1
        enddo
        buf = array(i)
        j   = i
        do
          call get_index(j)
          dum      = array(j)
          array(j) = buf
          buf      = dum
          call putdef(j,def)
          if(j==i) exit
        enddo
      enddo
 1    continue
      deallocate(def)

      contains
      
      subroutine get_index(index)
      implicit none
      integer, intent(inout) :: index
      integer                :: indices(n),index1,i,prod
      index1 = index
      do i = n,1,-1
        prod       = product(dims(:i-1))
        indices(i) = (index1-1) / prod
        index1     = index1 - indices(i) * prod
      enddo
      index = 1
      do i = 1,n
        index = index + indices(perm(i)) * product(dims(perm(:i-1)))
      enddo
      end subroutine get_index

      end subroutine reshape_permute_ SUFX

c     ----------------

c     Resorts array according to external routine get_index. The call get_index(index,args) turns an index
c     of the old ordering into the one of the new ordering (args(:) are optional integer arguments).
c     m = size of array      

      subroutine array_resort_ SUFX (array,m,get_index,args)
      use, intrinsic :: iso_fortran_env
      implicit none      
      Datatyp,  intent(inout) :: array(*)
      integer,  intent(in)    :: m,args(:)
      external                :: get_index
      integer,  allocatable   :: def(:)
      integer, parameter      :: bits = BITS ! number of bits used in integer array def
      Datatyp                 :: buf,dum
      integer                 :: i,j
      Allocate_(def,(0:m/bits))
      def = 0
      i   = 1
      do
        do while(isdef(i,def))
          i = i + 1
          if(i>m) goto 1
        enddo
        buf = array(i)
        j   = i
        do
          call get_index(j,args)
          dum      = array(j)
          array(j) = buf
          buf      = dum
          call putdef(j,def)
          if(j==i) exit
        enddo
      enddo
 1    continue
      Deallocate_(def)
      end subroutine array_resort_ SUFX

# endif

c     ----------------

# if IDatatyp == 1

c     Needed by reshape_permute and array_resort

      function isdef(index,def)
      implicit none
      logical             :: isdef
      integer, intent(in) :: index,def(0:)
      integer             :: idef,ibit
      integer, parameter  :: bits = BITS ! number of bits used in integer array def      
      idef  = index / bits
      ibit  = index - idef * bits
      isdef = iand(def(idef),2**ibit)/=0
      end function isdef

      subroutine putdef(index,def)
      integer, intent(in)    :: index
      integer, intent(inout) :: def(0:)
      integer                :: idef,ibit
      integer, parameter     :: bits = BITS ! number of bits used in integer array def      
      if(isdef(index,def)) Bug('Matrix element already resorted.')
      idef      = index / bits
      ibit      = index - idef * bits
      def(idef) = def(idef) + 2**ibit
      end subroutine putdef

# endif
      
c     ----------------      

c     Reallocates an allocated array.
c     Usage: call reallocate ( array, first new upper bound, second new upper bound, ... )
c            (Lower bounds cannot be changed.)

# if IDatatyp == 4
#   define ADD *(len(array))
#   define INIT ' '
# else
#   define ADD
#   define INIT 0
# endif

      subroutine reallocate_ SUFX 1(array,upper1B)
      use, intrinsic :: iso_fortran_env
      implicit none
      Datatyp, allocatable :: array(:),help(:)ADD
      integer :: upper1A,upper1B,lower1
      upper1A = ubound(array,1)
      lower1  = lbound(array,1)
      allocate (help(lower1:upper1A))
      help    = array
      deallocate (array)
      allocate (array(lower1:upper1B))
      array                 = INIT
      upper1A               = min(upper1A,upper1B)
      array(lower1:upper1A) = help(lower1:upper1A)
      deallocate(help)
      end subroutine reallocate_ SUFX 1

      subroutine reallocate_ SUFX 2(array,upper1B,upper2B)
      use, intrinsic :: iso_fortran_env
      implicit none
      Datatyp, allocatable :: array(:,:),help(:,:)ADD
      integer              :: upper1A,upper2A,upper1B,upper2B,lower1,lower2
      upper1A = ubound(array,1)
      upper2A = ubound(array,2)
      lower1  = lbound(array,1)
      lower2  = lbound(array,2)
      allocate (help(lower1:upper1A,lower2:upper2A))
      help    = array
      deallocate (array)
      allocate (array(lower1:upper1B,lower2:upper2B))
      array   = INIT
      upper1A = min(upper1A,upper1B)
      upper2A = min(upper2A,upper2B)
      array(lower1:upper1A,lower2:upper2A) = help(lower1:upper1A,lower2:upper2A)
      deallocate(help)
      end subroutine reallocate_ SUFX 2

      subroutine reallocate_ SUFX 3(array,upper1B,upper2B,upper3B)
      use, intrinsic :: iso_fortran_env
      implicit none
      Datatyp, allocatable :: array(:,:,:),help(:,:,:)ADD
      integer              :: upper1A,upper2A,upper3A,upper1B,upper2B,upper3B,lower1,lower2,lower3
      upper1A = ubound(array,1)
      upper2A = ubound(array,2)
      upper3A = ubound(array,3)
      lower1  = lbound(array,1)
      lower2  = lbound(array,2)
      lower3  = lbound(array,3)
      allocate (help(lower1:upper1A,lower2:upper2A,lower3:upper3A))
      help    = array
      deallocate (array)
      allocate (array(lower1:upper1B,lower2:upper2B,lower3:upper3B))
      array   = INIT
      upper1A = min(upper1A,upper1B)
      upper2A = min(upper2A,upper2B)
      upper3A = min(upper3A,upper3B)
      array(lower1:upper1A,lower2:upper2A,lower3:upper3A) = help(lower1:upper1A,lower2:upper2A,lower3:upper3A)
      deallocate(help)
      end subroutine reallocate_ SUFX 3
      
      subroutine reallocate_ SUFX 4(array,upper1B,upper2B,upper3B,upper4B)
      use, intrinsic :: iso_fortran_env
      implicit none
      Datatyp, allocatable :: array(:,:,:,:),help(:,:,:,:)ADD
      integer              :: upper1A,upper2A,upper3A,upper4A,upper1B,upper2B,upper3B,upper4B,lower1,lower2,lower3,lower4
      upper1A = ubound(array,1)
      upper2A = ubound(array,2)
      upper3A = ubound(array,3)
      upper4A = ubound(array,4)
      lower1  = lbound(array,1)
      lower2  = lbound(array,2)
      lower3  = lbound(array,3)
      lower4  = lbound(array,4)
      allocate (help(lower1:upper1A,lower2:upper2A,lower3:upper3A,lower4:upper4A))
      help    = array
      deallocate (array)
      allocate (array(lower1:upper1B,lower2:upper2B,lower3:upper3B,lower4:upper4B))
      array   = INIT
      upper1A = min(upper1A,upper1B)
      upper2A = min(upper2A,upper2B)
      upper3A = min(upper3A,upper3B)
      upper4A = min(upper4A,upper4B)
      array(lower1:upper1A,lower2:upper2A,lower3:upper3A,lower4:upper4A) =
     &  help(lower1:upper1A,lower2:upper2A,lower3:upper3A,lower4:upper4A)
      deallocate(help)
      end subroutine reallocate_ SUFX 4
      
      subroutine reallocate_ SUFX 5(array,upper1B,upper2B,upper3B,upper4B,upper5B)
      use, intrinsic :: iso_fortran_env
      implicit none
      Datatyp, allocatable :: array(:,:,:,:,:),help(:,:,:,:,:)ADD
      integer :: upper1A,upper2A,upper3A,upper4A,upper5A,upper1B,upper2B,upper3B,upper4B,upper5B,lower1,lower2,lower3,lower4,lower5
      upper1A = ubound(array,1)
      upper2A = ubound(array,2)
      upper3A = ubound(array,3)
      upper4A = ubound(array,4)
      upper5A = ubound(array,5)
      lower1  = lbound(array,1)
      lower2  = lbound(array,2)
      lower3  = lbound(array,3)
      lower4  = lbound(array,4)
      lower5  = lbound(array,5)
      allocate (help(lower1:upper1A,lower2:upper2A,lower3:upper3A,lower4:upper4A,lower5:upper5A))
      help    = array
      deallocate (array)
      allocate (array(lower1:upper1B,lower2:upper2B,lower3:upper3B,lower4:upper4B,lower5:upper5B))
      array   = INIT
      upper1A = min(upper1A,upper1B)
      upper2A = min(upper2A,upper2B)
      upper3A = min(upper3A,upper3B)
      upper4A = min(upper4A,upper4B)
      upper5A = min(upper5A,upper5B)
      array(lower1:upper1A,lower2:upper2A,lower3:upper3A,lower4:upper4A,lower5:upper5B) =
     &  help(lower1:upper1A,lower2:upper2A,lower3:upper3A,lower4:upper4A,lower5:upper5B)      
      deallocate(help)
      end subroutine reallocate_ SUFX 5
      
# undef ADD
# undef INIT
# undef SUFX
# undef Datatyp
# undef IDatatyp
