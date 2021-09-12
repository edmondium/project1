c--------------------------------------------------------------------------------
c Copyright (c) 2021 Peter Grünberg Institut, Forschungszentrum Jülich, Germany
c This file is part of SPEX and available as free software under the conditions
c of the MIT license as expressed in the LICENSE file in more detail.
c--------------------------------------------------------------------------------

# include "cppmacro.h"

      module Mwrapper

      use global, only: Mcomm
      use, intrinsic :: iso_c_binding
      use, intrinsic :: iso_fortran_env

# ifdef MaxChunkBytes
#   define Chunk MaxChunk / (storage_size(a(1))/8)
# else
#   define Chunk MaxChunk
# endif

      interface Mcast
      module procedure  Mcast_l,  Mcast_i,  Mcast_r,  Mcast_c,
     &                  Mcast_l1, Mcast_i1, Mcast_r1, Mcast_c1,
     &                  Mcast_l2, Mcast_i2, Mcast_r2, Mcast_c2,
     &                            Mcast_i3, Mcast_r3, Mcast_c3,
     &                                      Mcast_r4, Mcast_c4,
     &                                                Mcast_c5,
     &                  Mcast_string
      end interface

      interface Mcastl
      module procedure  Mcast_ll1,Mcast_il1,Mcast_rl1,Mcast_cl1,
     &                  Mcast_ll2,Mcast_il2,Mcast_rl2,Mcast_cl2,
     &                  Mcast_ll3,Mcast_il3,Mcast_rl3,Mcast_cl3,
     &                            Mcast_il4,Mcast_rl4,Mcast_cl4,
     &                                      Mcast_rl5,Mcast_cl5
      end interface

      interface Msum
      module procedure Msum_i,Msum_r,Msum_c,Msum_i1,Msum_i2,Msum_i3,Msum_r1,Msum_r2,Msum_r3,Msum_r4,Msum_r5,
     &                                      Msum_c1,Msum_c2,Msum_c3,Msum_c4,Msum_c5,Msum_c6
      end interface

      contains

c     --------

      subroutine Mcast1_l(a,n,rank,comm)
      implicit none
      include 'mpif.h'
      logical,           intent(inout) :: a(*)
      integer,           intent(in)    :: rank,comm
      integer(c_size_t), intent(in)    :: n
      integer(c_size_t)                :: i
      integer                          :: Merr
      integer,           parameter     :: chunk = Chunk
      do i = 1,n,chunk
        call mpi_bcast(a(i),int(min(n-i+1,chunk)),mpi_logical,rank,comm,Merr)
      enddo
      end subroutine Mcast1_l

      subroutine Mcast1_i(a,n,rank,comm)
      implicit none
      include 'mpif.h'
      integer,           intent(inout) :: a(*)
      integer,           intent(in)    :: rank,comm
      integer(c_size_t), intent(in)    :: n
      integer(c_size_t)                :: i
      integer                          :: Merr
      integer,           parameter     :: chunk = Chunk      
      do i = 1,n,chunk
        call mpi_bcast(a(i),int(min(n-i+1,chunk)),mpi_integer,rank,comm,Merr)
      enddo
      end subroutine Mcast1_i

      subroutine Mcast1_r(a,n,rank,comm)
      implicit none
      include 'mpif.h'
      real_dp,           intent(inout) :: a(*)
      integer,           intent(in)    :: rank,comm
      integer(c_size_t), intent(in)    :: n
      integer(c_size_t)                :: i
      integer                          :: Merr
      integer,           parameter     :: chunk = Chunk
      do i = 1,n,chunk
        call mpi_bcast(a(i),int(min(n-i+1,chunk)),mpi_double_precision,rank,comm,Merr)
      enddo
      end subroutine Mcast1_r

      subroutine Mcast1_c(a,n,rank,comm)
      implicit none
      include 'mpif.h'
      complex_dp,        intent(inout) :: a(*)
      integer,           intent(in)    :: rank,comm
      integer(c_size_t), intent(in)    :: n
      integer(c_size_t)                :: i
      integer                          :: Merr
      integer(c_size_t), parameter     :: chunk = Chunk
      do i = 1,n,chunk
        call mpi_bcast(a(i),int(min(n-i+1,chunk)),mpi_double_complex,rank,comm,Merr)
      enddo
      end subroutine Mcast1_c

c     --------

      subroutine Mcast_l(a,rank,comm)
      implicit none
      include 'mpif.h'
      logical                       :: a
      integer, intent(in), optional :: rank,comm
      integer                       :: Merr,rnk,cmm
      rnk = 0     ; if(present(rank)) rnk = rank
      cmm = Mcomm ; if(present(comm)) cmm = comm
      call mpi_bcast(a,1,mpi_logical,rnk,cmm,Merr)
      end subroutine Mcast_l

c     --------

      subroutine Mcast_i(a,rank,comm)
      implicit none
      include 'mpif.h'
      integer                       :: a
      integer, intent(in), optional :: rank,comm
      integer                       :: Merr,rnk,cmm
      rnk = 0     ; if(present(rank)) rnk = rank
      cmm = Mcomm ; if(present(comm)) cmm = comm
      call mpi_bcast(a,1,mpi_integer,rnk,cmm,Merr)
      end subroutine Mcast_i

c     --------

      subroutine Mcast_r(a,rank,comm)
      implicit none
      include 'mpif.h'
      real_dp                       :: a
      integer, intent(in), optional :: rank,comm
      integer                       :: Merr,rnk,cmm
      rnk = 0     ; if(present(rank)) rnk = rank
      cmm = Mcomm ; if(present(comm)) cmm = comm
      call mpi_bcast(a,1,mpi_double_precision,rnk,cmm,Merr)
      end subroutine Mcast_r

c     --------

      subroutine Mcast_c(a,rank,comm)
      implicit none
      include 'mpif.h'
      complex_dp                       :: a
      integer,    intent(in), optional :: rank,comm
      integer                          :: Merr,rnk,cmm
      rnk = 0     ; if(present(rank)) rnk = rank
      cmm = Mcomm ; if(present(comm)) cmm = comm
      call mpi_bcast(a,1,mpi_double_complex,rnk,cmm,Merr)
      end subroutine Mcast_c

c     --------

      subroutine Mcast_l1(a,rank,comm)
      use global, only: Mrank
      implicit none
      include 'mpif.h'
      logical                       :: a(:)
      integer, intent(in), optional :: rank,comm
      integer(c_size_t)             :: n
      integer                       :: rnk,cmm
      rnk = 0     ; if(present(rank)) rnk = rank
      cmm = Mcomm ; if(present(comm)) cmm = comm
      n   = size(a,kind=c_size_t)
      call Mcast1_l(a,n,rnk,cmm)
      end subroutine Mcast_l1

c     --------
      subroutine Mcast_l2(a,rank,comm)
      use global, only: Mrank
      implicit none
      include 'mpif.h'
      logical                       :: a(:,:)
      integer, intent(in), optional :: rank,comm
      integer(c_size_t)             :: n
      integer                       :: rnk,cmm
      rnk = 0     ; if(present(rank)) rnk = rank
      cmm = Mcomm ; if(present(comm)) cmm = comm
      n   = size(a,kind=c_size_t)
      call Mcast1_l(a,n,rnk,cmm)
      end subroutine Mcast_l2

c     --------

      subroutine Mcast_i1(a,rank,comm)
      use global, only: Mrank
      implicit none
      include 'mpif.h'
      integer                       :: a(:)
      integer, intent(in), optional :: rank,comm
      integer(c_size_t)             :: n
      integer                       :: rnk,cmm
      rnk = 0     ; if(present(rank)) rnk = rank
      cmm = Mcomm ; if(present(comm)) cmm = comm
      n   = size(a,kind=c_size_t)
      call Mcast1_i(a,n,rnk,cmm)
      end subroutine Mcast_i1

c     --------

      subroutine Mcast_i2(a,rank,comm)
      use global, only: Mrank
      implicit none
      include 'mpif.h'
      integer                       :: a(:,:)
      integer, intent(in), optional :: rank,comm
      integer(c_size_t)             :: n
      integer                       :: rnk,cmm
      rnk = 0     ; if(present(rank)) rnk = rank
      cmm = Mcomm ; if(present(comm)) cmm = comm
      n   = size(a,kind=c_size_t)
      call Mcast1_i(a,n,rnk,cmm)
      end subroutine Mcast_i2

c     --------

      subroutine Mcast_i3(a,rank,comm)
      use global, only: Mrank
      implicit none
      include 'mpif.h'
      integer                       :: a(:,:,:)
      integer, intent(in), optional :: rank,comm
      integer(c_size_t)             :: n
      integer                       :: rnk,cmm
      rnk = 0     ; if(present(rank)) rnk = rank
      cmm = Mcomm ; if(present(comm)) cmm = comm
      n   = size(a,kind=c_size_t)
      call Mcast1_i(a,n,rnk,cmm)
      end subroutine Mcast_i3

c     --------

      subroutine Mcast_r1(a,rank,comm)
      use global, only: Mrank
      implicit none
      include 'mpif.h'
      real_dp                       :: a(:)
      integer, intent(in), optional :: rank,comm
      integer(c_size_t)             :: n
      integer                       :: rnk,cmm
      rnk = 0     ; if(present(rank)) rnk = rank
      cmm = Mcomm ; if(present(comm)) cmm = comm
      n   = size(a,kind=c_size_t)
      call Mcast1_r(a,n,rnk,cmm)
      end subroutine Mcast_r1

c     --------

      subroutine Mcast_r2(a,rank,comm)
      use global, only: Mrank
      implicit none
      include 'mpif.h'
      real_dp                       :: a(:,:)
      integer, intent(in), optional :: rank,comm
      integer(c_size_t)             :: n
      integer                       :: rnk,cmm
      rnk = 0     ; if(present(rank)) rnk = rank
      cmm = Mcomm ; if(present(comm)) cmm = comm
      n   = size(a,kind=c_size_t)
      call Mcast1_r(a,n,rnk,cmm)
      end subroutine Mcast_r2

c     --------

      subroutine Mcast_r3(a,rank,comm)
      use global, only: Mrank
      implicit none
      include 'mpif.h'
      real_dp                       :: a(:,:,:)
      integer, intent(in), optional :: rank,comm
      integer(c_size_t)             :: n
      integer                       :: rnk,cmm
      rnk = 0     ; if(present(rank)) rnk = rank
      cmm = Mcomm ; if(present(comm)) cmm = comm
      n   = size(a,kind=c_size_t)
      call Mcast1_r(a,n,rnk,cmm)
      end subroutine Mcast_r3

c     --------

      subroutine Mcast_r4(a,rank,comm)
      use global, only: Mrank
      implicit none
      include 'mpif.h'
      real_dp                       :: a(:,:,:,:)
      integer, intent(in), optional :: rank,comm
      integer(c_size_t)             :: n
      integer                       :: rnk,cmm
      rnk = 0     ; if(present(rank)) rnk = rank
      cmm = Mcomm ; if(present(comm)) cmm = comm
      n   = size(a,kind=c_size_t)
      call Mcast1_r(a,n,rnk,cmm)
      end subroutine Mcast_r4

c     --------

      subroutine Mcast_c1(a,rank,comm)
      use global, only: Mrank
      implicit none
      include 'mpif.h'
      complex_dp                       :: a(:)
      integer,    intent(in), optional :: rank,comm
      integer(c_size_t)                :: n
      integer                          :: rnk,cmm
      rnk = 0     ; if(present(rank)) rnk = rank
      cmm = Mcomm ; if(present(comm)) cmm = comm
      n   = size(a,kind=c_size_t)
      call Mcast1_c(a,n,rnk,cmm)
      end subroutine Mcast_c1

c     --------

      subroutine Mcast_c2(a,rank,comm)
      use global, only: Mrank
      implicit none
      include 'mpif.h'
      complex_dp                       :: a(:,:)
      integer,    intent(in), optional :: rank,comm
      integer(c_size_t)                :: n
      integer                          :: rnk,cmm
      rnk = 0     ; if(present(rank)) rnk = rank
      cmm = Mcomm ; if(present(comm)) cmm = comm
      n   = size(a,kind=c_size_t)
      call Mcast1_c(a,n,rnk,cmm)
      end subroutine Mcast_c2

c     --------

      subroutine Mcast_c3(a,rank,comm)
      use global, only: Mrank
      implicit none
      include 'mpif.h'
      complex_dp                       :: a(:,:,:)
      integer,    intent(in), optional :: rank,comm
      integer(c_size_t)                :: n
      integer                          :: rnk,cmm
      rnk = 0     ; if(present(rank)) rnk = rank
      cmm = Mcomm ; if(present(comm)) cmm = comm
      n   = size(a,kind=c_size_t)
      call Mcast1_c(a,n,rnk,cmm)
      end subroutine Mcast_c3

c     --------

      subroutine Mcast_c4(a,rank,comm)
      use global, only: Mrank
      implicit none
      include 'mpif.h'
      complex_dp                       :: a(:,:,:,:)
      integer,    intent(in), optional :: rank,comm
      integer(c_size_t)                :: n
      integer                          :: rnk,cmm
      rnk = 0     ; if(present(rank)) rnk = rank
      cmm = Mcomm ; if(present(comm)) cmm = comm
      n   = size(a,kind=c_size_t)
      call Mcast1_c(a,n,rnk,cmm)
      end subroutine Mcast_c4

c     --------

      subroutine Mcast_c5(a,rank,comm)
      use global, only: Mrank
      implicit none
      include 'mpif.h'
      complex_dp                       :: a(:,:,:,:,:)
      integer,    intent(in), optional :: rank,comm
      integer(c_size_t)                :: n
      integer                          :: rnk,cmm
      rnk = 0     ; if(present(rank)) rnk = rank
      cmm = Mcomm ; if(present(comm)) cmm = comm
      n   = size(a,kind=c_size_t)
      call Mcast1_c(a,n,rnk,cmm)
      end subroutine Mcast_c5

c     --------

      subroutine Mcast_string(string,rank,comm)
      implicit none
      include 'mpif.h'
      character(*)                       :: string
      integer,      intent(in), optional :: rank,comm
      integer                            :: Merr,rnk,cmm
      rnk = 0     ; if(present(rank)) rnk = rank
      cmm = Mcomm ; if(present(comm)) cmm = comm
      call mpi_bcast(string,len(string),mpi_character,rnk,cmm,Merr)
      end subroutine Mcast_string

c     --------

      subroutine Mcast_ll1(a,rank,comm)
      use global, only: Mrank
      implicit none
      include 'mpif.h'
      logical, allocatable          :: a(:)
      integer, intent(in), optional :: rank,comm
      integer(c_size_t)             :: n
      integer                       :: b(2),Merr,rnk,cmm
      rnk = 0     ; if(present(rank)) rnk = rank
      cmm = Mcomm ; if(present(comm)) cmm = comm
      if(Mrank==rnk) then
        b(1) = huge(0)
        if(allocated(a)) then
          b(1) = lbound(a,1) ; b(2) = ubound(a,1)
        endif
      else
        if(allocated(a)) Bug('Array already allocated.')
      endif
      call mpi_bcast(b,2,mpi_integer,rnk,cmm,Merr)
      if(b(1)==huge(0)) return
      if(Mrank/=rnk) allocate ( a(b(1):b(2)) )
      n = size(a,kind=c_size_t)
      call Mcast1_l(a,n,rnk,cmm)
      end subroutine Mcast_ll1

c     --------

      subroutine Mcast_ll2(a,rank,comm)
      use global, only: Mrank
      implicit none
      include 'mpif.h'
      logical, allocatable          :: a(:,:)
      integer, intent(in), optional :: rank,comm
      integer(c_size_t)             :: n
      integer                       :: b(2,2),Merr,rnk,cmm
      rnk = 0     ; if(present(rank)) rnk = rank
      cmm = Mcomm ; if(present(comm)) cmm = comm
      if(Mrank/=rnk) then
        b(1,1) = huge(0)
        if(allocated(a)) then
          b(1,:) = lbound(a) ; b(2,:) = ubound(a)
        endif
      else
        if(allocated(a)) Bug('Array already allocated.')
      endif
      call mpi_bcast(b,4,mpi_integer,rnk,cmm,Merr)
      if(b(1,1)==huge(0)) return
      if(Mrank/=rnk) allocate ( a(b(1,1):b(2,1),b(1,2):b(2,2)) )
      n = size(a,kind=c_size_t)
      call Mcast1_l(a,n,rnk,cmm)
      end subroutine Mcast_ll2

c     --------

      subroutine Mcast_ll3(a,rank,comm)
      use global, only: Mrank
      implicit none
      include 'mpif.h'
      logical, allocatable          :: a(:,:,:)
      integer, intent(in), optional :: rank,comm
      integer(c_size_t)             :: n
      integer                       :: b(2,3),Merr,rnk,cmm
      rnk = 0     ; if(present(rank)) rnk = rank
      cmm = Mcomm ; if(present(comm)) cmm = comm
      if(Mrank==rnk) then
        b(1,1) = huge(0)
        if(allocated(a)) then
          b(1,:) = lbound(a) ; b(2,:) = ubound(a)
        endif
      else
        if(allocated(a)) Bug('Array already allocated.')
      endif
      call mpi_bcast(b,6,mpi_integer,rnk,cmm,Merr)
      if(b(1,1)==huge(0)) return
      if(Mrank/=rnk) allocate ( a(b(1,1):b(2,1),b(1,2):b(2,2),b(1,3):b(2,3)) )
      n = size(a,kind=c_size_t)
      call Mcast1_l(a,n,rnk,cmm)
      end subroutine Mcast_ll3

c     --------

      subroutine Mcast_il1(a,rank,comm)
      use global, only: Mrank
      implicit none
      include 'mpif.h'
      integer, allocatable          :: a(:)
      integer, intent(in), optional :: rank,comm
      integer(c_size_t)             :: n
      integer                       :: b(2),Merr,rnk,cmm
      rnk = 0     ; if(present(rank)) rnk = rank
      cmm = Mcomm ; if(present(comm)) cmm = comm
      if(Mrank==rnk) then
        b(1) = huge(0)
        if(allocated(a)) then
          b(1) = lbound(a,1) ; b(2) = ubound(a,1)
        endif
      else
        if(allocated(a)) Bug('Array already allocated.')
      endif
      call mpi_bcast(b,2,mpi_integer,rnk,cmm,Merr)
      if(b(1)==huge(0)) return
      if(Mrank/=rnk) allocate ( a(b(1):b(2)) )
      n = size(a,kind=c_size_t)
      call Mcast1_i(a,n,rnk,cmm)
      end subroutine Mcast_il1

c     --------

      subroutine Mcast_il2(a,rank,comm)
      use global, only: Mrank
      implicit none
      include 'mpif.h'
      integer, allocatable          :: a(:,:)
      integer, intent(in), optional :: rank,comm
      integer(c_size_t)             :: n
      integer                       :: b(2,2),Merr,rnk,cmm
      rnk = 0     ; if(present(rank)) rnk = rank
      cmm = Mcomm ; if(present(comm)) cmm = comm
      if(Mrank==rnk) then
        b(1,1) = huge(0)
        if(allocated(a)) then
          b(1,:) = lbound(a) ; b(2,:) = ubound(a)
        endif
      else
        if(allocated(a)) Bug('Array already allocated.')
      endif
      call mpi_bcast(b,4,mpi_integer,rnk,cmm,Merr)
      if(b(1,1)==huge(0)) return
      if(Mrank/=rnk) allocate ( a(b(1,1):b(2,1),b(1,2):b(2,2)) )
      n = size(a,kind=c_size_t)
      call Mcast1_i(a,n,rnk,cmm)
      end subroutine Mcast_il2

c     --------

      subroutine Mcast_il3(a,rank,comm)
      use global, only: Mrank
      implicit none
      include 'mpif.h'
      integer, allocatable          :: a(:,:,:)
      integer, intent(in), optional :: rank,comm
      integer(c_size_t)             :: n
      integer                       :: b(2,3),Merr,rnk,cmm
      rnk = 0     ; if(present(rank)) rnk = rank
      cmm = Mcomm ; if(present(comm)) cmm = comm
      if(Mrank==rnk) then
        b(1,1) = huge(0)
        if(allocated(a)) then
          b(1,:) = lbound(a) ; b(2,:) = ubound(a)
        endif
      else
        if(allocated(a)) Bug('Array already allocated.')
      endif
      call mpi_bcast(b,6,mpi_integer,rnk,cmm,Merr)
      if(b(1,1)==huge(0)) return
      if(Mrank/=rnk) allocate ( a(b(1,1):b(2,1),b(1,2):b(2,2),b(1,3):b(2,3)) )
      n = size(a,kind=c_size_t)
      call Mcast1_i(a,n,rnk,cmm)
      end subroutine Mcast_il3

c     --------

      subroutine Mcast_il4(a,rank,comm)
      use global, only: Mrank
      implicit none
      include 'mpif.h'
      integer, allocatable          :: a(:,:,:,:)
      integer, intent(in), optional :: rank,comm
      integer(c_size_t)             :: n
      integer                       :: b(2,4),Merr,rnk,cmm
      rnk = 0     ; if(present(rank)) rnk = rank
      cmm = Mcomm ; if(present(comm)) cmm = comm
      if(Mrank==rnk) then
        b(1,1) = huge(0)
        if(allocated(a)) then
          b(1,:) = lbound(a) ; b(2,:) = ubound(a)
        endif
      else
        if(allocated(a)) Bug('Array already allocated.')
      endif
      call mpi_bcast(b,8,mpi_integer,rnk,cmm,Merr)
      if(b(1,1)==huge(0)) return
      if(Mrank/=rnk) allocate ( a(b(1,1):b(2,1),b(1,2):b(2,2),b(1,3):b(2,3),b(1,4):b(2,4)) )
      n = size(a,kind=c_size_t)
      call Mcast1_i(a,n,rnk,cmm)
      end subroutine Mcast_il4

c     --------

      subroutine Mcast_rl1(a,rank,comm)
      use global, only: Mrank
      implicit none
      include 'mpif.h'
      real_dp, allocatable          :: a(:)
      integer, intent(in), optional :: rank,comm
      integer(c_size_t)             :: n
      integer                       :: b(2),Merr,rnk,cmm
      rnk = 0     ; if(present(rank)) rnk = rank
      cmm = Mcomm ; if(present(comm)) cmm = comm
      if(Mrank==rnk) then
        b(1) = huge(0)
        if(allocated(a)) then
          b(1) = lbound(a,1) ; b(2) = ubound(a,1)
        endif
      else
        if(allocated(a)) Bug('Array already allocated.')
      endif
      call mpi_bcast(b,2,mpi_integer,rnk,cmm,Merr)
      if(b(1)==huge(0)) return
      if(Mrank/=rnk) allocate ( a(b(1):b(2)) )
      n = size(a,kind=c_size_t)
      call Mcast1_r(a,n,rnk,cmm)
      end subroutine Mcast_rl1

c     --------

      subroutine Mcast_rl2(a,rank,comm)
      use global, only: Mrank
      implicit none
      include 'mpif.h'
      real_dp, allocatable          :: a(:,:)
      integer, intent(in), optional :: rank,comm
      integer(c_size_t)             :: n
      integer                       :: b(2,2),Merr,rnk,cmm
      rnk = 0     ; if(present(rank)) rnk = rank
      cmm = Mcomm ; if(present(comm)) cmm = comm
      if(Mrank==rnk) then
        b(1,1) = huge(0)
        if(allocated(a)) then
          b(1,:) = lbound(a) ; b(2,:) = ubound(a)
        endif
      else
        if(allocated(a)) Bug('Array already allocated.')
      endif
      call mpi_bcast(b,4,mpi_integer,rnk,cmm,Merr)
      if(b(1,1)==huge(0)) return
      if(Mrank/=rnk) allocate ( a(b(1,1):b(2,1),b(1,2):b(2,2)) )
      n = size(a,kind=c_size_t)
      call Mcast1_r(a,n,rnk,cmm)
      end subroutine Mcast_rl2

c     --------

      subroutine Mcast_rl3(a,rank,comm)
      use global, only: Mrank
      implicit none
      include 'mpif.h'
      real_dp, allocatable          :: a(:,:,:)
      integer, intent(in), optional :: rank,comm
      integer(c_size_t)             :: n
      integer                       :: b(2,3),Merr,rnk,cmm
      rnk = 0     ; if(present(rank)) rnk = rank
      cmm = Mcomm ; if(present(comm)) cmm = comm
      if(Mrank==rnk) then
        b(1,1) = huge(0)
        if(allocated(a)) then
          b(1,:) = lbound(a) ; b(2,:) = ubound(a)
        endif
      else
        if(allocated(a)) Bug('Array already allocated.')
      endif
      call mpi_bcast(b,6,mpi_integer,rnk,cmm,Merr)
      if(b(1,1)==huge(0)) return
      if(Mrank/=rnk) allocate ( a(b(1,1):b(2,1),b(1,2):b(2,2),b(1,3):b(2,3)) )
      n = size(a,kind=c_size_t)
      call Mcast1_r(a,n,rnk,cmm)
      end subroutine Mcast_rl3

c     --------

      subroutine Mcast_rl4(a,rank,comm)
      use global, only: Mrank
      implicit none
      include 'mpif.h'
      real_dp, allocatable          :: a(:,:,:,:)
      integer, intent(in), optional :: rank,comm
      integer(c_size_t)             :: n
      integer                       :: b(2,4),Merr,rnk,cmm
      rnk = 0     ; if(present(rank)) rnk = rank
      cmm = Mcomm ; if(present(comm)) cmm = comm
      if(Mrank==rnk) then
        b(1,1) = huge(0)
        if(allocated(a)) then
          b(1,:) = lbound(a) ; b(2,:) = ubound(a)
        endif
      else
        if(allocated(a)) Bug('Array already allocated.')
      endif
      call mpi_bcast(b,8,mpi_integer,rnk,cmm,Merr)
      if(b(1,1)==huge(0)) return
      if(Mrank/=rnk) allocate ( a(b(1,1):b(2,1),b(1,2):b(2,2),b(1,3):b(2,3),b(1,4):b(2,4)) )
      n = size(a,kind=c_size_t)
      call Mcast1_r(a,n,rnk,cmm)
      end subroutine Mcast_rl4

c     --------

      subroutine Mcast_rl5(a,rank,comm)
      use global, only: Mrank
      implicit none
      include 'mpif.h'
      real_dp, allocatable          :: a(:,:,:,:,:)
      integer, intent(in), optional :: rank,comm
      integer(c_size_t)             :: n
      integer                       :: b(2,5),Merr,rnk,cmm
      rnk = 0     ; if(present(rank)) rnk = rank
      cmm = Mcomm ; if(present(comm)) cmm = comm
      if(Mrank==rnk) then
        b(1,1) = huge(0)
        if(allocated(a)) then
          b(1,:) = lbound(a) ; b(2,:) = ubound(a)
        endif
      else
        if(allocated(a)) Bug('Array already allocated.')
      endif
      call mpi_bcast(b,10,mpi_integer,rnk,cmm,Merr)
      if(b(1,1)==huge(0)) return
      if(Mrank/=rnk) allocate ( a(b(1,1):b(2,1),b(1,2):b(2,2),b(1,3):b(2,3),b(1,4):b(2,4),b(1,5):b(2,5)) )
      n = size(a,kind=c_size_t)
      call Mcast1_r(a,n,rnk,cmm)
      end subroutine Mcast_rl5

c     --------

      subroutine Mcast_cl1(a,rank,comm)
      use global, only: Mrank
      implicit none
      include 'mpif.h'
      complex_dp, allocatable       :: a(:)
      integer, intent(in), optional :: rank,comm
      integer(c_size_t)             :: n
      integer                       :: b(2),Merr,rnk,cmm
      rnk = 0     ; if(present(rank)) rnk = rank
      cmm = Mcomm ; if(present(comm)) cmm = comm
      if(Mrank==rnk) then
        b(1) = huge(0)
        if(allocated(a)) then
          b(1) = lbound(a,1) ; b(2) = ubound(a,1)
        endif
      else
        if(allocated(a)) Bug('Array already allocated.')
      endif
      call mpi_bcast(b,2,mpi_integer,rnk,cmm,Merr)
      if(b(1)==huge(0)) return
      if(Mrank/=rnk) allocate ( a(b(1):b(2)) )
      n = size(a,kind=c_size_t)
      call Mcast1_c(a,n,rnk,cmm)
      end subroutine Mcast_cl1

c     --------

      subroutine Mcast_cl2(a,rank,comm)
      use global, only: Mrank
      implicit none
      include 'mpif.h'
      complex_dp, allocatable       :: a(:,:)
      integer, intent(in), optional :: rank,comm
      integer(c_size_t)             :: n
      integer                       :: b(2,2),Merr,rnk,cmm
      rnk = 0     ; if(present(rank)) rnk = rank
      cmm = Mcomm ; if(present(comm)) cmm = comm
      if(Mrank==rnk) then
        b(1,1) = huge(0)
        if(allocated(a)) then
          b(1,:) = lbound(a) ; b(2,:) = ubound(a)
        endif
      else
        if(allocated(a)) Bug('Array already allocated.')
      endif
      call mpi_bcast(b,4,mpi_integer,rnk,cmm,Merr)
      if(b(1,1)==huge(0)) return
      if(Mrank/=rnk) allocate ( a(b(1,1):b(2,1),b(1,2):b(2,2)) )
      n = size(a,kind=c_size_t)
      call Mcast1_c(a,n,rnk,cmm)
      end subroutine Mcast_cl2

c     --------

      subroutine Mcast_cl3(a,rank,comm)
      use global, only: Mrank
      implicit none
      include 'mpif.h'
      complex_dp, allocatable       :: a(:,:,:)
      integer, intent(in), optional :: rank,comm
      integer(c_size_t)             :: n
      integer                       :: b(2,3),Merr,rnk,cmm
      rnk = 0     ; if(present(rank)) rnk = rank
      cmm = Mcomm ; if(present(comm)) cmm = comm
      if(Mrank==rnk) then
        b(1,1) = huge(0)
        if(allocated(a)) then
          b(1,:) = lbound(a) ; b(2,:) = ubound(a)
        endif
      else
        if(allocated(a)) Bug('Array already allocated.')
      endif
      call mpi_bcast(b,6,mpi_integer,rnk,cmm,Merr)
      if(b(1,1)==huge(0)) return
      if(Mrank/=rnk) allocate ( a(b(1,1):b(2,1),b(1,2):b(2,2),b(1,3):b(2,3)) )
      n = size(a,kind=c_size_t)
      call Mcast1_c(a,n,rnk,cmm)
      end subroutine Mcast_cl3

c     --------

      subroutine Mcast_cl4(a,rank,comm)
      use global, only: Mrank
      implicit none
      include 'mpif.h'
      complex_dp, allocatable       :: a(:,:,:,:)
      integer, intent(in), optional :: rank,comm
      integer(c_size_t)             :: n
      integer                       :: b(2,4),Merr,rnk,cmm
      rnk = 0     ; if(present(rank)) rnk = rank
      cmm = Mcomm ; if(present(comm)) cmm = comm
      if(Mrank==rnk) then
        b(1,1) = huge(0)
        if(allocated(a)) then
          b(1,:) = lbound(a) ; b(2,:) = ubound(a)
        endif
      else
        if(allocated(a)) Bug('Array already allocated.')
      endif
      call mpi_bcast(b,8,mpi_integer,rnk,cmm,Merr)
      if(b(1,1)==huge(0)) return
      if(Mrank/=rnk) allocate ( a(b(1,1):b(2,1),b(1,2):b(2,2),b(1,3):b(2,3),b(1,4):b(2,4)) )
      n = size(a,kind=c_size_t)
      call Mcast1_c(a,n,rnk,cmm)
      end subroutine Mcast_cl4

c     --------

      subroutine Mcast_cl5(a,rank,comm)
      use global, only: Mrank
      implicit none
      include 'mpif.h'
      complex_dp, allocatable       :: a(:,:,:,:,:)
      integer, intent(in), optional :: rank,comm
      integer(c_size_t)             :: n
      integer                       :: b(2,5),Merr,rnk,cmm
      rnk = 0     ; if(present(rank)) rnk = rank
      cmm = Mcomm ; if(present(comm)) cmm = comm
      if(Mrank==rnk) then
        b(1,1) = huge(0)
        if(allocated(a)) then
          b(1,:) = lbound(a) ; b(2,:) = ubound(a)
        endif
      else
        if(allocated(a)) Bug('Array already allocated.')
      endif
      call mpi_bcast(b,10,mpi_integer,rnk,cmm,Merr)
      if(b(1,1)==huge(0)) return
      if(Mrank/=rnk) allocate ( a(b(1,1):b(2,1),b(1,2):b(2,2),b(1,3):b(2,3),b(1,4):b(2,4),b(1,5):b(2,5)) )
      n = size(a,kind=c_size_t)
      call Mcast1_c(a,n,rnk,cmm)
      end subroutine Mcast_cl5

c     --------

      subroutine Msum1_i(a,n,rank,comm)
      use global, only: Mrank
      implicit none
      include 'mpif.h'
      integer,           intent(inout) :: a(*)
      integer,           intent(in)    :: rank,comm
      integer(c_size_t), intent(in)    :: n
      integer(c_size_t)                :: i
      integer                          :: Merr
      integer,           parameter     :: chunk = Chunk
      if(Mrank==rank) then
        do i = 1,n,chunk
          call mpi_reduce(mpi_in_place,a(i),int(min(n-i+1,chunk)),mpi_integer,mpi_sum,rank,comm,Merr)
        enddo
      else if(rank>=0) then
        do i = 1,n,chunk
          call mpi_reduce(a(i),           0,int(min(n-i+1,chunk)),mpi_integer,mpi_sum,rank,comm,Merr)
        enddo
      else
        do i = 1,n,chunk
          call mpi_allreduce(mpi_in_place,a(i),int(min(n-i+1,chunk)),mpi_integer,mpi_sum,comm,Merr)
        enddo
      endif
      end subroutine Msum1_i

c     --------

      subroutine Msum1_r(a,n,rank,comm)
      use global, only: Mrank
      implicit none
      include 'mpif.h'
      real_dp,           intent(inout) :: a(*)
      integer,           intent(in)    :: rank,comm
      integer(c_size_t), intent(in)    :: n
      integer(c_size_t)                :: i
      integer                          :: Merr
      integer,           parameter     :: chunk = Chunk
      if(Mrank==rank) then
        do i = 1,n,chunk
          call mpi_reduce(mpi_in_place,a(i),int(min(n-i+1,chunk)),mpi_double_precision,mpi_sum,rank,comm,Merr)
        enddo
      else if(rank>=0) then
        do i = 1,n,chunk
          call mpi_reduce(a(i),           0,int(min(n-i+1,chunk)),mpi_double_precision,mpi_sum,rank,comm,Merr)
        enddo
      else
        do i = 1,n,chunk
          call mpi_allreduce(mpi_in_place,a(i),int(min(n-i+1,chunk)),mpi_double_precision,mpi_sum,comm,Merr)
        enddo
      endif
      end subroutine Msum1_r

c     --------

      subroutine Msum1_c(a,n,rank,comm)
      use global, only: Mrank
      implicit none
      include 'mpif.h'
      complex_dp,        intent(inout) :: a(*)
      integer,           intent(in)    :: rank,comm
      integer(c_size_t), intent(in)    :: n
      integer(c_size_t)                :: i
      integer                          :: Merr
      integer,           parameter     :: chunk = Chunk
      if(Mrank==rank) then
        do i = 1,n,chunk
          call mpi_reduce(mpi_in_place,a(i),int(min(n-i+1,chunk)),mpi_double_complex,mpi_sum,rank,comm,Merr)
        enddo
      else if(rank>=0) then
        do i = 1,n,chunk
          call mpi_reduce(a(i),           0,int(min(n-i+1,chunk)),mpi_double_complex,mpi_sum,rank,comm,Merr)
        enddo
      else
        do i = 1,n,chunk
          call mpi_allreduce(mpi_in_place,a(i),int(min(n-i+1,chunk)),mpi_double_complex,mpi_sum,comm,Merr)
        enddo
      endif
      end subroutine Msum1_c

c     --------

      subroutine Msum_i(a,rank,comm)
      use global, only: Mrank
      implicit none
      include 'mpif.h'
      integer                       :: a,b
      integer                       :: Merr,cmm
      integer, intent(in), optional :: rank,comm
      cmm = Mcomm ; if(present(comm)) cmm = comm
      if(present(rank)) then
        call mpi_reduce(a,b,1,mpi_integer,mpi_sum,rank,cmm,Merr)
        if(Mrank==rank) a = b
      else
        call mpi_allreduce(a,b,1,mpi_integer,mpi_sum,cmm,Merr)
        a = b
      endif
      end subroutine Msum_i

c     --------

      subroutine Msum_r(a,rank,comm)
      use global, only: Mrank
      implicit none
      include 'mpif.h'
      real_dp                       :: a,b
      integer                       :: Merr,cmm
      integer, intent(in), optional :: rank,comm
      cmm = Mcomm ; if(present(comm)) cmm = comm
      if(present(rank)) then
        call mpi_reduce(a,b,1,mpi_double_precision,mpi_sum,rank,cmm,Merr)
        if(Mrank==rank) a = b
      else
        call mpi_allreduce(a,b,1,mpi_double_precision,mpi_sum,cmm,Merr)
        a = b
      endif
      end subroutine Msum_r

c     --------

      subroutine Msum_c(a,rank,comm)
      use global, only: Mrank
      implicit none
      include 'mpif.h'
      complex_dp                    :: a,b
      integer                       :: Merr,cmm
      integer, intent(in), optional :: rank,comm
      cmm = Mcomm ; if(present(comm)) cmm = comm
      if(present(rank)) then
        call mpi_reduce(a,b,1,mpi_double_complex,mpi_sum,rank,cmm,Merr)
        if(Mrank==rank) a = b
      else
        call mpi_allreduce(a,b,1,mpi_double_complex,mpi_sum,cmm,Merr)
        a = b
      endif
      end subroutine Msum_c

c     --------

      subroutine Msum_i1(a,rank,comm)
      use global, only: Mrank
      implicit none
      integer                       :: a(:)
      integer, intent(in), optional :: rank,comm
      integer                       :: rnk,cmm
      integer(c_size_t)             :: n
      cmm = Mcomm ; if(present(comm)) cmm = comm
      rnk = -1    ; if(present(rank)) rnk = rank
      n   = size(a,kind=c_size_t)
      call Msum1_i(a,n,rnk,cmm)
      end subroutine Msum_i1

c     --------

      subroutine Msum_i2(a,rank,comm)
      use global, only: Mrank
      implicit none
      integer                       :: a(:,:)
      integer, intent(in), optional :: rank,comm
      integer                       :: rnk,cmm
      integer(c_size_t)             :: n
      cmm = Mcomm ; if(present(comm)) cmm = comm
      rnk = -1    ; if(present(rank)) rnk = rank
      n   = size(a,kind=c_size_t)
      call Msum1_i(a,n,rnk,cmm)
      end subroutine Msum_i2

c     --------

      subroutine Msum_i3(a,rank,comm)
      use global, only: Mrank
      implicit none
      integer                       :: a(:,:,:)
      integer, intent(in), optional :: rank,comm
      integer                       :: rnk,cmm
      integer(c_size_t)             :: n
      cmm = Mcomm ; if(present(comm)) cmm = comm
      rnk = -1    ; if(present(rank)) rnk = rank
      n   = size(a,kind=c_size_t)
      call Msum1_i(a,n,rnk,cmm)
      end subroutine Msum_i3

c     --------

      subroutine Msum_r1(a,rank,comm)
      use global, only: Mrank
      implicit none
      real_dp                       :: a(:)
      integer, intent(in), optional :: rank,comm
      integer                       :: rnk,cmm
      integer(c_size_t)             :: n
      cmm = Mcomm ; if(present(comm)) cmm = comm
      rnk = -1    ; if(present(rank)) rnk = rank
      n   = size(a,kind=c_size_t)
      call Msum1_r(a,n,rnk,cmm)
      end subroutine Msum_r1

c     --------

      subroutine Msum_r2(a,rank,comm)
      use global, only: Mrank
      implicit none
      real_dp                       :: a(:,:)
      integer, intent(in), optional :: rank,comm
      integer                       :: rnk,cmm
      integer(c_size_t)             :: n
      cmm = Mcomm ; if(present(comm)) cmm = comm
      rnk = -1    ; if(present(rank)) rnk = rank
      n   = size(a,kind=c_size_t)
      call Msum1_r(a,n,rnk,cmm)
      end subroutine Msum_r2

c     --------

      subroutine Msum_r3(a,rank,comm)
      use global, only: Mrank
      implicit none
      real_dp                       :: a(:,:,:)
      integer, intent(in), optional :: rank,comm
      integer                       :: rnk,cmm
      integer(c_size_t)             :: n
      cmm = Mcomm ; if(present(comm)) cmm = comm
      rnk = -1    ; if(present(rank)) rnk = rank
      n   = size(a,kind=c_size_t)
      call Msum1_r(a,n,rnk,cmm)
      end subroutine Msum_r3

c     --------

      subroutine Msum_r4(a,rank,comm)
      use global, only: Mrank
      implicit none
      real_dp                       :: a(:,:,:,:)
      integer, intent(in), optional :: rank,comm
      integer                       :: rnk,cmm
      integer(c_size_t)             :: n
      cmm = Mcomm ; if(present(comm)) cmm = comm
      rnk = -1    ; if(present(rank)) rnk = rank
      n   = size(a,kind=c_size_t)
      call Msum1_r(a,n,rnk,cmm)
      end subroutine Msum_r4

c     --------

      subroutine Msum_r5(a,rank,comm)
      use global, only: Mrank
      implicit none
      real_dp                       :: a(:,:,:,:,:)
      integer, intent(in), optional :: rank,comm
      integer                       :: rnk,cmm
      integer(c_size_t)             :: n
      cmm = Mcomm ; if(present(comm)) cmm = comm
      rnk = -1    ; if(present(rank)) rnk = rank
      n   = size(a,kind=c_size_t)
      call Msum1_r(a,n,rnk,cmm)
      end subroutine Msum_r5

c     --------

      subroutine Msum_c1(a,rank,comm)
      use global, only: Mrank
      implicit none
      complex_dp                    :: a(:)
      integer, intent(in), optional :: rank,comm
      integer                       :: rnk,cmm
      integer(c_size_t)             :: n
      cmm = Mcomm ; if(present(comm)) cmm = comm
      rnk = -1    ; if(present(rank)) rnk = rank
      n   = size(a,kind=c_size_t)
      call Msum1_c(a,n,rnk,cmm)
      end subroutine Msum_c1

c     --------

      subroutine Msum_c2(a,rank,comm)
      use global, only: Mrank
      implicit none
      complex_dp                    :: a(:,:)
      integer, intent(in), optional :: rank,comm
      integer                       :: rnk,cmm
      integer(c_size_t)             :: n
      cmm = Mcomm ; if(present(comm)) cmm = comm
      rnk = -1    ; if(present(rank)) rnk = rank
      n   = size(a,kind=c_size_t)
      call Msum1_c(a,n,rnk,cmm)
      end subroutine Msum_c2

c     --------

      subroutine Msum_c3(a,rank,comm)
      use global, only: Mrank
      implicit none
      complex_dp                    :: a(:,:,:)
      integer, intent(in), optional :: rank,comm
      integer                       :: rnk,cmm
      integer(c_size_t)             :: n
      cmm = Mcomm ; if(present(comm)) cmm = comm
      rnk = -1    ; if(present(rank)) rnk = rank
      n   = size(a,kind=c_size_t)
      call Msum1_c(a,n,rnk,cmm)
      end subroutine Msum_c3

c     --------

      subroutine Msum_c4(a,rank,comm)
      use global, only: Mrank
      implicit none
      complex_dp                    :: a(:,:,:,:)
      integer, intent(in), optional :: rank,comm
      integer                       :: rnk,cmm
      integer(c_size_t)             :: n
      cmm = Mcomm ; if(present(comm)) cmm = comm
      rnk = -1    ; if(present(rank)) rnk = rank
      n   = size(a,kind=c_size_t)
      call Msum1_c(a,n,rnk,cmm)
      end subroutine Msum_c4

c     --------

      subroutine Msum_c5(a,rank,comm)
      use global, only: Mrank
      implicit none
      complex_dp                    :: a(:,:,:,:,:)
      integer, intent(in), optional :: rank,comm
      integer                       :: rnk,cmm
      integer(c_size_t)             :: n
      cmm = Mcomm ; if(present(comm)) cmm = comm
      rnk = -1    ; if(present(rank)) rnk = rank
      n   = size(a,kind=c_size_t)
      call Msum1_c(a,n,rnk,cmm)
      end subroutine Msum_c5

c     --------

      subroutine Msum_c6(a,rank,comm)
      use global, only: Mrank
      implicit none
      complex_dp                    :: a(:,:,:,:,:,:)
      integer, intent(in), optional :: rank,comm
      integer                       :: rnk,cmm
      integer(c_size_t)             :: n
      cmm = Mcomm ; if(present(comm)) cmm = comm
      rnk = -1    ; if(present(rank)) rnk = rank
      n   = size(a,kind=c_size_t)
      call Msum1_c(a,n,rnk,cmm)
      end subroutine Msum_c6

c     --------

      subroutine node_allocate(win,ptr,byte,dims)
      use global, only: Ncomm,Nrank
      use, intrinsic :: iso_c_binding
      implicit none
      include 'mpif.h'
      type(c_ptr), intent(out)       :: ptr
      integer,     intent(out)       :: win
      integer,     intent(in)        :: byte,dims(:)
      integer                        :: disp,Merr,dim,i
      integer(kind=mpi_address_kind) :: msize
      do i = 1,size(dims)
        MpiO( dim = dims(i) )
        call Mcast(dim,comm=Ncomm)
        MnoO( if(dim/=dims(i)) Bug('Dimensions differ from those of Orank.') )
      enddo
      if(Nrank==0) then ; msize = product(dims*1_mpi_address_kind) * byte
      else              ; msize = 0_mpi_address_kind
      endif
      if(msize<0) Bug('Negative array size.')
      disp = 1
      call mpi_win_allocate_shared(msize,disp,mpi_info_null,Ncomm,ptr,win,Merr)
      if(Nrank/=0) call mpi_win_shared_query(win,0,msize,disp,ptr,Merr)
      end subroutine node_allocate

# if 0
      subroutine Mgather_r1(a,na,rank)
      use global, only: Mrank,Msize
      implicit none
      include 'mpif.h'
      real_dp                       :: a(:),b(size(a))
      integer, intent(in)           :: na
      integer                       :: Merr,nb(Msize),d(Msize)
      integer, intent(in), optional :: rank
      if(present(rank)) then
        call mpi_gather(na,1,mpi_integer,nb,1,mpi_integer,rank,Mcomm,Merr)
        if(Mrank==rank) then
          do i = 1,Msize
            d(i) = sum(nb(:i-1)) + 1
          enddo
        endif
        call mpi_gatherv(a(d(Mrange)),nb(Mrange),mpi_double_precision,b,nb,d,mpi_double_precision,rank,Mcomm,Merr)
        if(Mrank==rank) a = b
      else
        Error('not implemented')
        call mpi_allgatherv(a,b,size(a),mpi_double_precision,mpi_sum,Mcomm,Merr)
        a = b
      endif
      end subroutine Mgather_r1
# endif

c     --------

      end module Mwrapper

c     --------




