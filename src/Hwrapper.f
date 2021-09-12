c--------------------------------------------------------------------------------
c Copyright (c) 2021 Peter Grünberg Institut, Forschungszentrum Jülich, Germany
c This file is part of SPEX and available as free software under the conditions
c of the MIT license as expressed in the LICENSE file in more detail.
c--------------------------------------------------------------------------------

c Wrapper routines for reading from and writing to HDF5 files.
c      
c hdf_fopen  : Opens HDF5 file.
c hdf_fclose : Closes HDF5 file.
c hdf_rdwr   : Reads or writes data set.
c hdf_rdwr_a : Reads or writes data as attribute.
c hdf_dim    : Queries dimensions of data set.
c
c See below for more details.

# if defined(MPI) && defined(HDF5ser)
#   ifndef M_PART   /* only give warning once */
#     warning MPI but HDF5 library only serial; this might lead to extra I/O overhead and is not well tested yet.
#   endif
#   undef MPI
#   define MPI_
# endif

# include "cppmacro.h"

# ifdef MaxChunkBytes
#   define Chunk MaxChunk / (storage_size(buf)/8)      
# else
#   define Chunk MaxChunk
# endif

# define ERRSTOP(arg) if(Herr<0) Error('Fatal HDF5 error: '//trim(arg))
# define LARGEATT(name,sub) if(storage_size(buf)/8*size(buf)>64*1024) then ; if(mode>0) Info('Large attribute '//name//' (>64kB) written as dataset.') ; call sub ; return ; endif

# ifndef M_PART           /* subdivide file into parts (necessary for preprocessor loops) */
#   define M_PART 1
#   include "Hwrapper.f"
#   undef M_PART
#   define M_PART 2
#   include "Hwrapper.f"
#   undef M_PART
#   define M_PART 3
#   include "Hwrapper.f"
#   undef M_PART
#   define M_PART 4
# endif

c --------------
# if M_PART == 1
c --------------

      module Hwrapper

      use, intrinsic :: iso_fortran_env

      integer :: Hpos ! array dimension that offset refers to (rightmost is zeroth)

      interface hdf_rdwr
      module procedure hdf_rdwr_i1,hdf_rdwr_i2,hdf_rdwr_i3,hdf_rdwr_i4,hdf_rdwr_i5,
     &                 hdf_rdwr_r1,hdf_rdwr_r2,hdf_rdwr_r3,hdf_rdwr_r4,hdf_rdwr_r5,
     &                 hdf_rdwr_c1,hdf_rdwr_c2,hdf_rdwr_c3,hdf_rdwr_c4
      end interface

      interface hdf_rdwr_a
      module procedure hdf_rdwr_a_i, hdf_rdwr_a_r, hdf_rdwr_a_str,
     &                 hdf_rdwr_a_i1,hdf_rdwr_a_r1,hdf_rdwr_a_c1
      end interface

      private :: array_divide

      contains

c     --------

c     Opens a HDF5 file "Hname" and returns file handle in Hfile.
c       mode = 0 : readonly           (file opened for reading only)
c              1 : readwrite          (file opened for reading and writing; create file if it does not exist)
c              2 : truncate&readwrite (file opened for reading and writing; always create new file)

      subroutine hdf_fopen(Hfile,Hname,mode)
      use hdf5
# ifdef MPI
      use Mwrapper
      use global, only: Mrank,Mcomm
# endif
      implicit none
      character(*),     intent(in)  :: Hname
      integer,          intent(in)  :: mode
      integer(HID_T),   intent(out) :: Hfile
      integer                       :: Herr,Hmode
      logical                       :: exist,Hpar
      Mpi( integer(HID_T)           :: Hplist )
      Mpi( include 'mpif.h' )
      if(mode==0) then ; Hmode = H5F_ACC_RDONLY_F
      else             ; Hmode = H5F_ACC_RDWR_F
      endif
      ifR inquire(file=Hname,exist=exist)
      Rif(mode==0.and..not.exist) Error('File '//trim(Hname)//' not found.')
# ifdef MPI
      call mpi_initialized(Hpar,Herr)
      if(Hpar) then
        call Mcast(exist)
        call h5pcreate_f(H5P_FILE_ACCESS_F,Hplist,Herr)          ; ERRSTOP(Hname)
        call h5pset_fapl_mpio_f(Hplist,Mcomm,mpi_info_null,Herr) ; ERRSTOP(Hname)
        if(mode==2.or..not.exist) then
          call h5fcreate_f(Hname,H5F_ACC_TRUNC_F,Hfile,Herr,
     &                                        access_prp=Hplist) ; ERRSTOP(Hname)
        else
          call h5fopen_f(Hname,Hmode,Hfile,Herr,Hplist)          ; ERRSTOP(Hname)
        endif
        call h5pclose_f(Hplist,Herr)                             ; ERRSTOP(Hname)
        return
      endif
# endif
      if(mode==2.or..not.exist) then
        call h5fcreate_f(Hname,H5F_ACC_TRUNC_F,Hfile,Herr)       ; ERRSTOP(Hname)
      else
        call h5fopen_f(Hname,Hmode,Hfile,Herr)                   ; ERRSTOP(Hname)
      endif
      end subroutine hdf_fopen

c     --------

c     Closes HDF5 file (with handle Hfile)

      subroutine hdf_fclose(Hfile)
      use hdf5
      implicit none
      integer(HID_T), intent(in) :: Hfile
      integer                    :: Herr
      call h5fclose_f(Hfile,Herr) ; ERRSTOP('file close')
      end subroutine hdf_fclose

c     --------

c     Queries dimensions of data set "Hname" at Hloc

      subroutine hdf_dim(Hloc,Hname,Hdim)
      use hdf5
      implicit none
      character(*),     intent(in)  :: Hname
      integer(HID_T),   intent(in)  :: Hloc
      integer(HID_T)                :: Hset,Hspc
      integer(HSIZE_T), intent(out) :: Hdim(:)
      integer(HSIZE_T)              :: Hmaxdim(size(Hdim))
      integer                       :: Herr
      call h5dopen_f(Hloc,Hname,Hset,Herr)                     ; ERRSTOP(Hname)
      call h5dget_space_f(Hset,Hspc,Herr)                      ; ERRSTOP(Hname)
      call h5sget_simple_extent_dims_f(Hspc,Hdim,Hmaxdim,Herr) ; ERRSTOP(Hname)
      if(Herr/=size(Hdim)) Error('Wrong number of dimensions.')
      call h5sclose_f(Hspc,Herr)                               ; ERRSTOP(Hname)
      call h5dclose_f(Hset,Herr)                               ; ERRSTOP(Hname)
      end subroutine hdf_dim

      subroutine hdf_rdwr_a_str(Hloc0,Hname0,mode,str)
      use hdf5
      implicit none
      integer(HID_T), intent(in) :: Hloc0
      integer,        intent(in) :: mode
      character(*),   intent(in) :: Hname0
      character(len(Hname0))     :: Hname
      character(*)               :: str
      integer(HID_T)             :: Hspc,Hatt,Htyp,Hloc
      integer(HSIZE_T)           :: Hdim(1)
      integer(SIZE_T)            :: Hlen
      integer                    :: Herr
      call hdf_gopen(Hloc,Hname,Hloc0,Hname0)
      if(mode==2) Error('mode==2 not implemented for this routine.')
      Hlen    = len(str)
      Hdim(1) = 1
      call h5tcopy_f(H5T_NATIVE_CHARACTER,Htyp,Herr)     ; ERRSTOP(Hname)
      call h5tset_size_f(Htyp,Hlen,Herr)                 ; ERRSTOP(Hname)
      if(mode==0) then
        call h5aopen_f(Hloc,Hname,Hatt,Herr)             ; ERRSTOP(Hname)
        call h5aread_f(Hatt,Htyp,str,Hdim,Herr)          ; ERRSTOP(Hname)
      else
        call h5screate_simple_f(1,Hdim,Hspc,Herr)        ; ERRSTOP(Hname)
        call h5acreate_f(Hloc,Hname,Htyp,Hspc,Hatt,Herr) ; ERRSTOP(Hname)
        call h5awrite_f(Hatt,Htyp,str,Hdim,Herr)         ; ERRSTOP(Hname)
        call h5sclose_f(Hspc,Herr)                       ; ERRSTOP(Hname)
      endif
      call h5aclose_f(Hatt,Herr)                         ; ERRSTOP(Hname)
      if(Hloc/=Hloc0) call h5gclose_f(Hloc,Herr)         ; ERRSTOP(Hname)
      end subroutine hdf_rdwr_a_str

c     Read/Write data attribute "Hname" from/to Hloc
c       mode = 0 : read
c              1 : write
c              2 : write/check (write if attibute does not exist, otherwise check if it is identical to buf)

      recursive subroutine hdf_rdwr_a_i1(Hloc0,Hname0,mode,buf)
      use hdf5
      implicit none
      integer(HID_T), intent(in) :: Hloc0
      integer,        intent(in) :: mode
      character(*),   intent(in) :: Hname0
      character(len(Hname0))     :: Hname
      logical                    :: ldum
      integer                    :: buf(:),buf1(size(buf))
      integer(HID_T)             :: Hspc,Hatt,Hloc
      integer(HSIZE_T)           :: Hdim(1)
      integer                    :: Herr
      LARGEATT(Hname0,hdf_rdwr_i1(Hloc0,Hname0,mode,buf))
      call hdf_gopen(Hloc,Hname,Hloc0,Hname0)
      Hdim(1) = size(buf)
      if(mode==0) then
        call h5aopen_f(Hloc,Hname,Hatt,Herr)                           ; ERRSTOP(Hname)
        call h5aread_f(Hatt,H5T_NATIVE_INTEGER,buf,Hdim,Herr)          ; ERRSTOP(Hname)
      else
        if(mode==2) then
          call h5aexists_f(Hloc,Hname,ldum,Herr)                       ; ERRSTOP(Hname)
          if(ldum) then
            call hdf_rdwr_a(Hloc,Hname,0,buf1)
            if(any(buf/=buf1)) Error('Attribute(s) '//trim(Hname)//' incorrect.')
            return
          endif
        endif
        call h5screate_simple_f(1,Hdim,Hspc,Herr)                      ; ERRSTOP(Hname)
        call h5acreate_f(Hloc,Hname,H5T_NATIVE_INTEGER,Hspc,Hatt,Herr) ; ERRSTOP(Hname)
        call h5awrite_f(Hatt,H5T_NATIVE_INTEGER,buf,Hdim,Herr)         ; ERRSTOP(Hname)
        call h5sclose_f(Hspc,Herr)                                     ; ERRSTOP(Hname)
      endif
      call h5aclose_f(Hatt,Herr)                                       ; ERRSTOP(Hname)
      if(Hloc/=Hloc0) call h5gclose_f(Hloc,Herr)                       ; ERRSTOP(Hname)
      end subroutine hdf_rdwr_a_i1

      recursive subroutine hdf_rdwr_a_r1(Hloc0,Hname0,mode,buf)
      use hdf5
      implicit none
      integer(HID_T), intent(in) :: Hloc0
      integer,        intent(in) :: mode
      character(*),   intent(in) :: Hname0
      character(len(Hname0))     :: Hname
      real_dp                    :: buf(:),buf1(size(buf))
      logical                    :: ldum
      integer(HID_T)             :: Hspc,Hatt,Hloc
      integer(HSIZE_T)           :: Hdim(1)
      integer                    :: Herr
      LARGEATT(Hname0,hdf_rdwr_r1(Hloc0,Hname0,mode,buf))
      call hdf_gopen(Hloc,Hname,Hloc0,Hname0)
      Hdim(1) = size(buf)
      if(mode==0) then
        call h5aopen_f(Hloc,Hname,Hatt,Herr)                          ; ERRSTOP(Hname)
        call h5aread_f(Hatt,H5T_NATIVE_DOUBLE,buf,Hdim,Herr)          ; ERRSTOP(Hname)
      else
        if(mode==2) then
          call h5aexists_f(Hloc,Hname,ldum,Herr)                      ; ERRSTOP(Hname)
          if(ldum) then
            call hdf_rdwr_a(Hloc,Hname,0,buf1)
            if(any(abs(buf-buf1)>1d-10)) Error('Attribute(s) '//trim(Hname)//' incorrect.')
            return
          endif
        endif
        call h5screate_simple_f(1,Hdim,Hspc,Herr)                     ; ERRSTOP(Hname)
        call h5acreate_f(Hloc,Hname,H5T_NATIVE_DOUBLE,Hspc,Hatt,Herr) ; ERRSTOP(Hname)
        call h5awrite_f(Hatt,H5T_NATIVE_DOUBLE,buf,Hdim,Herr)         ; ERRSTOP(Hname)
        call h5sclose_f(Hspc,Herr)                                    ; ERRSTOP(Hname)
      endif
      call h5aclose_f(Hatt,Herr)                                      ; ERRSTOP(Hname)
      if(Hloc/=Hloc0) call h5gclose_f(Hloc,Herr)                      ; ERRSTOP(Hname)
      end subroutine hdf_rdwr_a_r1

      recursive subroutine hdf_rdwr_a_c1(Hloc0,Hname0,mode,buf)
      use hdf5
      implicit none
      integer(HID_T), intent(in) :: Hloc0
      integer,        intent(in) :: mode
      character(*),   intent(in) :: Hname0
      character(len(Hname0))     :: Hname
      complex_dp                 :: buf(:),buf1(size(buf))
      real_dp                    :: rbuf(2,size(buf))
      logical                    :: ldum
      integer(HID_T)             :: Hspc,Hatt,Hloc
      integer(HSIZE_T)           :: Hdim(2)
      integer                    :: Herr
      LARGEATT(Hname0,hdf_rdwr_c1(Hloc0,Hname0,mode,buf))
      call hdf_gopen(Hloc,Hname,Hloc0,Hname0)
      Hdim(1) = 2
      Hdim(2) = size(buf)
      if(mode==0) then
        call h5aopen_f(Hloc,Hname,Hatt,Herr)                          ; ERRSTOP(Hname)
        call h5aread_f(Hatt,H5T_NATIVE_DOUBLE,rbuf,Hdim,Herr)         ; ERRSTOP(Hname)
        buf = rbuf(1,:) + (0d0,1d0) * rbuf(2,:)
      else
        if(mode==2) then
          call h5aexists_f(Hloc,Hname,ldum,Herr)                      ; ERRSTOP(Hname)
          if(ldum) then
            call hdf_rdwr_a(Hloc,Hname,0,buf1)
            if(any(abs(buf-buf1)>1d-10)) Error('Attribute(s) '//trim(Hname)//' incorrect.')
            return
          endif
        endif
        rbuf(1,:) = real(buf)
        rbuf(2,:) = imag(buf)
        call h5screate_simple_f(2,Hdim,Hspc,Herr)                     ; ERRSTOP(Hname)
        call h5acreate_f(Hloc,Hname,H5T_NATIVE_DOUBLE,Hspc,Hatt,Herr) ; ERRSTOP(Hname)
        call h5awrite_f(Hatt,H5T_NATIVE_DOUBLE,rbuf,Hdim,Herr)        ; ERRSTOP(Hname)
        call h5sclose_f(Hspc,Herr)                                    ; ERRSTOP(Hname)
      endif
      call h5aclose_f(Hatt,Herr)                                      ; ERRSTOP(Hname)
      if(Hloc/=Hloc0) call h5gclose_f(Hloc,Herr)                      ; ERRSTOP(Hname)
      end subroutine hdf_rdwr_a_c1

      subroutine hdf_rdwr_a_i(Hloc,Hname,mode,buf)
      use hdf5
      implicit none
      integer(HID_T), intent(in) :: Hloc
      integer,        intent(in) :: mode
      character(*),   intent(in) :: Hname
      integer                    :: buf
      integer                    :: buf1(1)
      if(mode/=0) buf1(1) = buf
      call hdf_rdwr_a_i1(Hloc,Hname,mode,buf1)
      if(mode==0) buf = buf1(1)
      end subroutine hdf_rdwr_a_i

      subroutine hdf_rdwr_a_r(Hloc,Hname,mode,buf)
      use hdf5
      implicit none
      integer(HID_T), intent(in) :: Hloc
      integer,        intent(in) :: mode
      character(*),   intent(in) :: Hname
      real_dp                    :: buf
      real_dp                    :: buf1(1)
      if(mode/=0) buf1(1) = buf
      call hdf_rdwr_a_r1(Hloc,Hname,mode,buf1)
      if(mode==0) buf = buf1(1)
      end subroutine hdf_rdwr_a_r

c     --------

c ----------------
# elif M_PART == 2
c ----------------

c     Read/Write data set "Hname" from/to Hloc
c
c     mode = 0 : read
c            1 : write
c            2 : write/check (write data set if it does not exist, otherwise just return)
c            3 : do not write, only create dataset (only dimensions of buf are used, not the contents)
c
c     off   = offset for dimension Hpos (see Hpos above)
c     str   = stride for dimension Hpos
c     indep = true:  independent parallel access (H5FD_MPIO_INDEPENDENT)
c             false: collective parallel access (H5FD_MPIO_COLLECTIVE)  (default)
c
c     mode = 3 and indep are not fully implemented yet.

c     --------

# ifndef TYP             /* outer loop over data types */
#   define TYP  integer
#   define TYP_ i
#   define TYPH H5T_NATIVE_INTEGER
#   include "Hwrapper.f"
#   define TYP  real_dp
#   define TYP_ r
#   define TYPH H5T_NATIVE_DOUBLE
# endif

# ifndef DIM             /* inner loop over dimensions */
#   define DIM 1
#   define DIMS :
#   define DIMB lb(1):ub(1)
#   include "Hwrapper.f"
#   define DIM 2
#   define DIMS :,:
#   define DIMB lb(1):ub(1),lb(2):ub(2)
#   include "Hwrapper.f"
#   define DIM 3
#   define DIMS :,:,:
#   define DIMB lb(1):ub(1),lb(2):ub(2),lb(3):ub(3)
#   include "Hwrapper.f"
#   define DIM 4
#   define DIMS :,:,:,:
#   define DIMB lb(1):ub(1),lb(2):ub(2),lb(3):ub(3),lb(4):ub(4)
#   include "Hwrapper.f"
#   define DIM 5
#   define DIMS :,:,:,:,:
#   define DIMB lb(1):ub(1),lb(2):ub(2),lb(3):ub(3),lb(4):ub(4),lb(5):ub(5)
#   include "Hwrapper.f"
#   undef TYP           /* undefine macros at end of outer loop */
#   undef TYP_
#   undef TYPH
# else

      recursive subroutine hdf_rdwr_ TYP_ DIM (Hloc,Hname,mode,buf,off,str,Hset0  MpiC(indep) )
      use hdf5
      implicit none
      integer(HID_T),    intent(in)    :: Hloc
      integer,           intent(in)    :: mode
      character(*),      intent(in)    :: Hname
      TYP                              :: buf( DIMS )
      integer, optional, intent(in)    :: off,str
      integer, optional, intent(inout) :: Hset0
      logical                          :: ldum,Hpar
      integer                          :: off1,str1
      integer(HID_T)                   :: Hset,Hspc,Hmem MpiC(Hplist)
      integer(HSIZE_T)                 :: Hdim( DIM ),Hoff( DIM ),Hstr( DIM )
      integer                          :: Herr
      integer                          :: lb( DIM ),ub( DIM ),d,i
#   ifdef MPI
      logical, optional, intent(in) :: indep
#   endif
      if(present(off)) then ; off1 = off ; else ; off1 = -1 ; endif
      if(present(str)) then ; str1 = str ; else ; str1 = -1 ; endif
      if(present(Hset0).and.any(mode==[0,2])) Bug('Dataset identifier Hset0 present but mode=0 or 1.')
      if(.not.present(Hset0).and.mode==3)     Bug('Dataset identifier Hset0 not present but mode=3.')
      if(present(Hset0).and.mode==1) Hset = Hset0
      if(mode==3) Error('Not yet implemented: mode==3.')
      if(size(buf,kind=int_dp)>Chunk) then
        if(present(Hset0)) Error('Dataset division not implemented for Hset0.')
        Info('Dataset '//trim(Hname)//' divided as its size exceeds an internal HDF5 limit.')
        if(off1>0.or.abs(str1)/=1) Error('off/=0 or str/=1. Dataset division not implemented for this case.')
        lb = lbound(buf)
        ub = ubound(buf)
        call array_divide(d,i,lb,ub)
        ub(d) = i
        call hdf_rdwr(Hloc,trim(Hname)//'_0',mode,buf( DIMB ),off1,str1)
        ub(d) = ubound(buf,d)
        lb(d) = i+1
        call hdf_rdwr(Hloc,trim(Hname)//'_1',mode,buf( DIMB ),off1,str1)
        return
      endif
      if(mode==2) then
        call h5lexists_f(Hloc,Hname,ldum,Herr) ; ERRSTOP(Hname)
        if(ldum) return
      endif
      Hdim = shape(buf)
#   ifdef MPI
      Hplist = h5p_default_f
      call mpi_initialized(Hpar,Herr)
      if(Hpar) then
        call h5pcreate_f(H5P_DATASET_XFER_F,Hplist,Herr)                                 ; ERRSTOP(Hname)
        ldum = present(indep)
        if(ldum) ldum = indep
        if(ldum) then
          call h5pset_dxpl_mpio_f(Hplist,H5FD_MPIO_INDEPENDENT_F,Herr)                   ; ERRSTOP(Hname)
        else
          call h5pset_dxpl_mpio_f(Hplist,H5FD_MPIO_COLLECTIVE_F,Herr)                    ; ERRSTOP(Hname)          
        endif
      endif
#   endif
      if(off1>=0) then ! read/write dataset hyperslabs (write: dataset must exist already)
        if(Hpos<0.or.Hpos>1) Bug('Wrong position.')
        Hoff = 0 ;             Hoff(max(1, DIM -Hpos)) = off1
        Hstr = 1 ; if(str1>=0) Hstr(max(1, DIM -Hpos)) = str1
        call h5dopen_f(Hloc,Hname,Hset,Herr)                                             ; ERRSTOP(Hname)
        call h5dget_space_f(Hset,Hspc,Herr)                                              ; ERRSTOP(Hname)
        call h5sselect_hyperslab_f(Hspc,H5S_SELECT_SET_F,Hoff,Hdim,Herr,Hstr)            ; ERRSTOP(Hname)
        call h5screate_simple_f( DIM ,Hdim,Hmem,Herr)                                    ; ERRSTOP(Hname)
        if(mode==0) then
          call h5dread_f(Hset, TYPH ,buf,Hdim,Herr,Hmem,Hspc MpiC(Hplist) )              ; ERRSTOP(Hname)
        else
          call h5dwrite_f(Hset, TYPH ,buf,Hdim,Herr,Hmem,Hspc MpiC(Hplist) )             ; ERRSTOP(Hname)
        endif
        call h5sclose_f(Hmem,Herr)                                                       ; ERRSTOP(Hname)
        call h5sclose_f(Hspc,Herr)                                                       ; ERRSTOP(Hname)
      else             ! read/write full datasets (write: dataset is created)
        if(mode==0) then
          call h5dopen_f(Hloc,Hname,Hset,Herr)                                           ; ERRSTOP(Hname)
          call h5dread_f(Hset, TYPH ,buf,Hdim,Herr MpiC(xfer_prp=Hplist) )               ; ERRSTOP(Hname)
        else
          call h5screate_simple_f( DIM ,Hdim,Hspc,Herr)                                  ; ERRSTOP(Hname)
          call h5dcreate_f(Hloc,Hname, TYPH ,Hspc,Hset,Herr)                             ; ERRSTOP(Hname)
          if(mode<=2) call h5dwrite_f(Hset, TYPH ,buf,Hdim,Herr MpiC(xfer_prp=Hplist) )  ; ERRSTOP(Hname)
          call h5sclose_f(Hspc,Herr)                                                     ; ERRSTOP(Hname)
        endif
      endif
#   ifdef MPI
      if(Hpar) call h5pclose_f(Hplist,Herr)                                              ; ERRSTOP(Hname)
#   endif
      call h5dclose_f(Hset,Herr)                                                         ; ERRSTOP(Hname)
      end subroutine hdf_rdwr_ TYP_ DIM

#   undef DIM           /* undefine macros at end of inner loop */
#   undef DIMS
#   undef DIMB
# endif
      
c     --------

c ----------------
# elif M_PART == 3
c ----------------

# ifndef DIM           /* loop over dimensions */
#   define DIM 1
#   define DIMS :
#   include "Hwrapper.f"
#   define DIM 2
#   define DIMS :,:
#   include "Hwrapper.f"
#   define DIM 3
#   define DIMS :,:,:
#   include "Hwrapper.f"
#   define DIM 4
#   define DIMS :,:,:,:
# endif

c The pointer solution leads to a segfault in the HDF5 library. Therefore, it is disabled:
# define NOPTR
      subroutine hdf_rdwr_c DIM (Hloc,Hname,mode,buf,off,str)
      use hdf5
      use iso_c_binding
      implicit none
      integer(HID_T),    intent(in) :: Hloc
      integer,           intent(in) :: mode
      character(*),      intent(in) :: Hname
      integer, optional, intent(in) :: off,str
# ifdef NOPTR
      complex_dp                    :: buf( DIMS )
      if(present(off)) then
        if(present(str)) then
          call real_hdf_rdwr_c DIM (Hloc,Hname,mode,buf,shape(buf),off,str)
        else
          call real_hdf_rdwr_c DIM (Hloc,Hname,mode,buf,shape(buf),off,-1)
        endif
      else
        call real_hdf_rdwr_c DIM (Hloc,Hname,mode,buf,shape(buf),-1,-1)
      endif
# else
      complex_dp,        target     :: buf( DIMS )
      real_dp,           pointer    :: rbuf(:, DIMS ) => null()      
      type(c_ptr)                   :: ptr
      ptr = c_loc(buf)
      call c_f_pointer(ptr,rbuf,shape=[2,shape(buf)])
      if(present(off)) then
        if(present(str)) then
          call hdf_rdwr(Hloc,Hname,mode,rbuf,off,str)
        else
          call hdf_rdwr(Hloc,Hname,mode,rbuf,off)
        endif
      else
        call hdf_rdwr(Hloc,Hname,mode,rbuf)
      endif
      nullify(rbuf)
# endif
      end subroutine hdf_rdwr_c DIM

#   undef DIM           /* undefine macros at end of loop */
#   undef DIMS

c     --------

c ----------------
# elif M_PART == 4
c ----------------

      subroutine hdf_gopen(Hloc,Hname,Hloc0,Hname0)
      use hdf5
      implicit none
      integer(HID_T), intent(out) :: Hloc
      integer(HID_T), intent(in)  :: Hloc0
      character(*),   intent(out) :: Hname
      character(*),   intent(in)  :: Hname0
      integer                     :: Herr,ind
      ind = index(Hname0,'/',back=.true.)
      if(ind/=0) then
        if(Hname0(:1)=='/') Bug('Leading "/".')
        call h5gopen_f(Hloc0,Hname0(:ind-1),Hloc,Herr) ; ERRSTOP(Hname)
        Hname = Hname0(ind+1:)
      else
        Hname = Hname0
        Hloc  = Hloc0
      endif
      end subroutine hdf_gopen

c     --------

      subroutine array_divide(d,i,lb,ub)
      implicit none
      integer, intent(out) :: d,i
      integer, intent(in)  :: lb(:),ub(:)
      if(size(lb)/=size(ub)) Bug('Different sizes of lb and ub.')
      do d = size(lb),1,-1
        if(ub(d)-lb(d)>0) then
          i = ( ub(d) - lb(d) ) / 2 + lb(d)
          return
        endif
      enddo
      Bug('Subroutine array_divide failed. Array cannot be divided.')
      end subroutine array_divide

c     --------      

      end module Hwrapper

c     --------

c     Real calls for complex datatype

c     --------

      subroutine real_hdf_rdwr_c1(Hloc,Hname,mode,buf,b,off,str)
      use hdf5
      use Hwrapper
      use, intrinsic :: iso_fortran_env
      implicit none
      integer(HID_T), intent(in) :: Hloc
      integer,        intent(in) :: mode,b(1)
      character(*),   intent(in) :: Hname
      real_dp                    :: buf(2,b(1))
      integer,        intent(in) :: off,str
      if(off==-1) then
        call hdf_rdwr(Hloc,Hname,mode,buf)
      else
        if(str==-1) then
          call hdf_rdwr(Hloc,Hname,mode,buf,off)
        else
          call hdf_rdwr(Hloc,Hname,mode,buf,off,str)
        endif
      endif
      end

c     --------

      subroutine real_hdf_rdwr_c2(Hloc,Hname,mode,buf,b,off,str)
      use hdf5
      use Hwrapper
      use, intrinsic :: iso_fortran_env
      implicit none
      integer(HID_T), intent(in) :: Hloc
      integer,        intent(in) :: mode,b(2)
      character(*),   intent(in) :: Hname
      real_dp                    :: buf(2,b(1),b(2))
      integer,        intent(in) :: off,str
      if(off==-1) then
        call hdf_rdwr(Hloc,Hname,mode,buf)
      else
        if(str==-1) then
          call hdf_rdwr(Hloc,Hname,mode,buf,off)
        else
          call hdf_rdwr(Hloc,Hname,mode,buf,off,str)
        endif
      endif
      end

c     --------

      subroutine real_hdf_rdwr_c3(Hloc,Hname,mode,buf,b,off,str)
      use hdf5
      use Hwrapper
      use, intrinsic :: iso_fortran_env
      implicit none
      integer(HID_T), intent(in) :: Hloc
      integer,        intent(in) :: mode,b(3)
      character(*),   intent(in) :: Hname
      real_dp                    :: buf(2,b(1),b(2),b(3))
      integer,        intent(in) :: off,str
      if(off==-1) then
        call hdf_rdwr(Hloc,Hname,mode,buf)
      else
        if(str==-1) then
          call hdf_rdwr(Hloc,Hname,mode,buf,off)
        else
          call hdf_rdwr(Hloc,Hname,mode,buf,off,str)
        endif
      endif
      end

c     --------

      subroutine real_hdf_rdwr_c4(Hloc,Hname,mode,buf,b,off,str)
      use hdf5
      use Hwrapper
      use, intrinsic :: iso_fortran_env
      implicit none
      integer(HID_T), intent(in) :: Hloc
      integer,        intent(in) :: mode,b(4)
      character(*),   intent(in) :: Hname
      real_dp                    :: buf(2,b(1),b(2),b(3),b(4))
      integer,        intent(in) :: off,str
      if(off==-1) then
        call hdf_rdwr(Hloc,Hname,mode,buf)
      else
        if(str==-1) then
          call hdf_rdwr(Hloc,Hname,mode,buf,off)
        else
          call hdf_rdwr(Hloc,Hname,mode,buf,off,str)
        endif
      endif
      end

# endif

# if defined(MPI_) && defined(HDF5ser)
#   undef MPI_
#   define MPI
# endif
