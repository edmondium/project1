c--------------------------------------------------------------------------------
c Copyright (c) 2021 Peter Grünberg Institut, Forschungszentrum Jülich, Germany
c This file is part of SPEX and available as free software under the conditions
c of the MIT license as expressed in the LICENSE file in more detail.
c--------------------------------------------------------------------------------

c     Initialization routine for Fast Fourier Transformation (FFT)
c     (fftw library)
c
# include "cppmacro.h"

      module m_fft

      use, intrinsic :: iso_fortran_env
      implicit none
# include <fftw3.f>
      integer                         :: nnfft,nfft(3)
      integer                         :: imat(3,3),imati(3,3)
      integer_dp                      :: plan_fft1,plan_fft2,plan_ffts,plan_fftp
# ifdef INV
      real_dp,    allocatable, target :: cfft1(:,:,:),cfft2(:,:,:),cffts(:,:,:),cfftp(:,:,:)
      complex_dp, allocatable, target :: rfft1(:,:,:),rfft2(:,:,:),rffts(:,:,:),rfftp(:,:,:)
      complex_dp, pointer_cnt         :: cfft1_w(:,:,:),cfft2_w(:,:,:),cffts_w(:,:,:),cfftp_w(:,:,:)
      integer_dp                      :: plan_fft1_w,plan_fft2_w,plan_ffts_w,plan_fftp_w
# else
      complex_dp, allocatable, target :: fft1(:,:,:),fft2(:,:,:),ffts(:,:,:),fftp(:,:,:)
      complex_dp, pointer_cnt         :: cfft1(:,:,:),cfft2(:,:,:),cffts(:,:,:),cfftp(:,:,:)
      complex_dp, pointer_cnt         :: rfft1(:,:,:),rfft2(:,:,:),rffts(:,:,:),rfftp(:,:,:)
# endif

      contains

c     --------------------------

c     Initializes Fast Fourier transform: defines box of G vectors that contains all points in the sphere of radius gcutf
      subroutine fft_init(gcutf)
      use global, only: rlat MpiC(Mrank)
      use util
      use, intrinsic :: iso_fortran_env
      implicit none
      real_dp, intent(in)  :: gcutf
      real_dp              :: rlat1(3,3)
      integer, parameter   :: ident(3,3) = reshape ( [ 1,0,0, 0,1,0, 0,0,1 ] , [3,3] )

      Rwrite(6,'(/A)') 'Initialization of FFT'

c     Redefine reciprocal lattice vectors according to imat (from determine_imat below), might be replaced by backfold once!
      rlat1 = matmul(rlat,imat)

c     Define box of G vectors by powers of 2, 3, and 5.
      call gcount_sphere(rlat1,gcutf)

# ifdef INV
      allocate ( cfft1(0:nfft(1)-1,0:nfft(2)-1,0:nfft(3)-1) )
      allocate ( cfft2(0:nfft(1)-1,0:nfft(2)-1,0:nfft(3)-1) )
      allocate ( cffts(0:nfft(1)-1,0:nfft(2)-1,0:nfft(3)-1) )
      allocate ( cfftp(0:nfft(1)-1,0:nfft(2)-1,0:nfft(3)-1) )
c#   ifndef WAN / this is needed only for wavefproducts3 I suppose. Could be changed for the slimmer def. below once wavefproducts3 is not used anymore.
c      allocate ( rfft1(0:nfft(1)/2,0:nfft(2)-1,0:nfft(3)-1) )
c      allocate ( rfft2(0:nfft(1)/2,0:nfft(2)-1,0:nfft(3)-1) )
c      allocate ( rffts(0:nfft(1)/2,0:nfft(2)-1,0:nfft(3)-1) )
c      allocate ( rfftp(0:nfft(1)/2,0:nfft(2)-1,0:nfft(3)-1) )
c#   else
      allocate ( rfft1(0:nfft(1)-1,0:nfft(2)-1,0:nfft(3)-1) )
      allocate ( rfft2(0:nfft(1)-1,0:nfft(2)-1,0:nfft(3)-1) )
      allocate ( rffts(0:nfft(1)-1,0:nfft(2)-1,0:nfft(3)-1) )
      allocate ( rfftp(0:nfft(1)-1,0:nfft(2)-1,0:nfft(3)-1) )
c#   endif
      call dfftw_plan_dft_r2c_3d(plan_fft1, nfft(1),nfft(2),nfft(3),cfft1,rfft1,FFTW_MEASURE)
      call dfftw_plan_dft_r2c_3d(plan_fft2, nfft(1),nfft(2),nfft(3),cfft2,rfft2,FFTW_MEASURE)
      call dfftw_plan_dft_r2c_3d(plan_ffts, nfft(1),nfft(2),nfft(3),cffts,rffts,FFTW_MEASURE)
      call dfftw_plan_dft_c2r_3d(plan_fftp, nfft(1),nfft(2),nfft(3),rfftp,cfftp,FFTW_MEASURE)
      cfft1_w => rfft1
      cfft2_w => rfft2
      cffts_w => rffts
      cfftp_w => rfftp
      call dfftw_plan_dft_3d(plan_fft1_w, nfft(1),nfft(2),nfft(3),cfft1_w,rfft1,FFTW_BACKWARD,FFTW_MEASURE)
      call dfftw_plan_dft_3d(plan_fft2_w, nfft(1),nfft(2),nfft(3),cfft2_w,rfft2,FFTW_BACKWARD,FFTW_MEASURE)
      call dfftw_plan_dft_3d(plan_ffts_w, nfft(1),nfft(2),nfft(3),cffts_w,rffts,FFTW_BACKWARD,FFTW_MEASURE)
      call dfftw_plan_dft_3d(plan_fftp_w, nfft(1),nfft(2),nfft(3),rfftp,cfftp_w,FFTW_FORWARD, FFTW_MEASURE)
# else
      allocate ( fft1(0:nfft(1)-1,0:nfft(2)-1,0:nfft(3)-1) )
      allocate ( fft2(0:nfft(1)-1,0:nfft(2)-1,0:nfft(3)-1) )
      allocate ( ffts(0:nfft(1)-1,0:nfft(2)-1,0:nfft(3)-1) )
      allocate ( fftp(0:nfft(1)-1,0:nfft(2)-1,0:nfft(3)-1) )
      cfft1 => fft1 ; rfft1 => fft1
      cfft2 => fft2 ; rfft2 => fft2
      cffts => ffts ; rffts => ffts
      cfftp => fftp ; rfftp => fftp
      call dfftw_plan_dft_3d(plan_fft1, nfft(1),nfft(2),nfft(3),cfft1,rfft1,FFTW_BACKWARD,FFTW_MEASURE)
      call dfftw_plan_dft_3d(plan_fft2, nfft(1),nfft(2),nfft(3),cfft2,rfft2,FFTW_BACKWARD,FFTW_MEASURE)
      call dfftw_plan_dft_3d(plan_ffts, nfft(1),nfft(2),nfft(3),cffts,rffts,FFTW_BACKWARD,FFTW_MEASURE)
      call dfftw_plan_dft_3d(plan_fftp, nfft(1),nfft(2),nfft(3),rfftp,cfftp,FFTW_FORWARD, FFTW_MEASURE)
# endif

      end subroutine fft_init

c     --------------------------

c     Returns the number of G points in a sphere of radius gcutf (according to lattice vectors defined in rlat1)
      subroutine gcount_sphere(rlat1,gcutf)
      Mpi2(use global, only: Mrank)
      use, intrinsic :: iso_fortran_env
      implicit none
      real_dp, intent(in) :: gcutf,rlat1(3,3)
      real_dp             :: vec(3),rdum
      logical             :: found
      integer             :: n,nn,i,i2,i3,i5,i7,ix,iy,iz
      n    = 0
      i    = 0
      nfft = -1000
      do
        found = .false.
        do ix = -n,n
          do iy = -(n-abs(ix)),n-abs(ix)
            iz   = n - abs(ix) - abs(iy)
 1          vec  = ix*rlat1(:,1) + iy*rlat1(:,2) + iz*rlat1(:,3)
            rdum = sum(vec**2)
            if(rdum<=gcutf**2) then
              i       = i + 1
              found   = .true.
              nfft(1) = max(nfft(1),ix)
              nfft(2) = max(nfft(2),iy)
              nfft(3) = max(nfft(3),iz)
            endif
            if(iz>0) then
              iz = -iz
              goto 1
            endif
          enddo
        enddo
        if(.not.found) exit
        n = n + 1
      enddo
      nfft = 2*nfft + 1
      Rwrite(6,'(2X,A,I7)') 'Number of G points in sphere:',i
      Rwrite(6,'(2X,A,''('',I3.3,''x'',I3.3,''x'',I3.3,'')'')') 'Minimal mesh: ',nfft

      if(any(nfft>2**16)) Error('more than 2^16 points in one direction.')
      do i = 1,3
        nn = huge(0)
        do i2 = 0,16
          do i3 = 0,11
            do i5 = 0,7
              do i7 = 0,6
                n = 2**i2 * 3**i3 * 5**i5 * 7**i7
                if(n>=nfft(i).and.n<nn) then ; nn = n ; endif
              enddo
            enddo
          enddo
        enddo
        if(nn<nfft(i)) Error('could not determine power decomposition.')
        nfft(i) = nn
      enddo
      nnfft = product(nfft)
      Rwrite(6,'(2X,A,''('',I3.3,''x'',I3.3,''x'',I3.3,'')'')') 'Mesh defined: ',nfft
      Rwrite(6,'(2X,A,I10)') 'Number of G points in FFT:',nnfft
      end subroutine gcount_sphere

c     --------------------------

c     Search for a set of reciprocal basis vectors (rlat1) that produces (or is likely to produce) an FFT box of minimal size (->imat)
      subroutine determine_imat
      use global, only: rlat MpiC(Mrank)
      use, intrinsic :: iso_fortran_env
      implicit none
      integer            :: imat1(3,3)
      integer            :: b
      integer, parameter :: ident(3,3) = reshape ( [ 1,0,0, 0,1,0, 0,0,1 ] , [3,3] )
      real_dp            :: rlat1(3,3),len2(3),rdum,rdum0
      integer            :: i11,i12,i13,i21,i22,i23,i31,i32,i33,i
      b     = 2
      imat  = ident
      len2  = [ (sum(rlat(:,i)**2),i=1,3) ]
      rdum0 = (dot_product(rlat(:,1),rlat(:,2)))**2/(len2(1)*len2(2)) +
     &        (dot_product(rlat(:,2),rlat(:,3)))**2/(len2(2)*len2(3)) +
     &        (dot_product(rlat(:,3),rlat(:,1)))**2/(len2(3)*len2(1))
      do i11 = -b,b
      do i12 = -b,b
      do i13 = -b,b
      do i21 = -b,b
      do i22 = -b,b
      do i23 = -b,b
      do i31 = -b,b
      do i32 = -b,b
      do i33 = -b,b
        i = i11*i22*i33 + i12*i23*i31 + i13*i21*i32 - i31*i22*i13 - i32*i23*i11 - i33*i21*i12
        if(i/=1) cycle
        imat1 = reshape ( [ i11,i21,i31,i12,i22,i32,i13,i23,i33 ] , [ 3,3 ] )
        rlat1 = matmul ( rlat , imat1 )
        len2  = [ (sum(rlat1(:,i)**2),i=1,3) ]
        rdum  = (dot_product(rlat1(:,1),rlat1(:,2)))**2/(len2(1)*len2(2)) +
     &          (dot_product(rlat1(:,2),rlat1(:,3)))**2/(len2(2)*len2(3)) +
     &          (dot_product(rlat1(:,3),rlat1(:,1)))**2/(len2(3)*len2(1))
        if(rdum<rdum0-1d-12) then
          rdum0 = rdum
          imat  = imat1
        endif
      enddo
      enddo
      enddo
      enddo
      enddo
      enddo
      enddo
      enddo
      enddo
      imati = reshape ( (/ determ(1,1) , -determ(1,2) ,  determ(1,3) ,
     &                    -determ(2,1) ,  determ(2,2) , -determ(2,3) ,
     &                     determ(3,1) , -determ(3,2) ,  determ(3,3) /) , [ 3,3 ] )
      if(any(matmul(imat,imati)/=ident)) Bug('calculation of inverse failed.')
      if(any(imat/=ident)) then
        Rwrite(6,'(/A)') 'BZ transformation matrix:'
        Rwrite(6,'(3I3)') transpose(imat)
      endif

      contains

      function determ(i,j)
      use, intrinsic :: iso_fortran_env
      implicit none
      integer             :: determ
      integer, intent(in) :: i,j
      integer             :: i1,i2,j1,j2
      if     (i==1) then ; i1 = 2 ; i2 = 3
      else if(i==2) then ; i1 = 1 ; i2 = 3
      else if(i==3) then ; i1 = 1 ; i2 = 2
      endif
      if     (j==1) then ; j1 = 2 ; j2 = 3
      else if(j==2) then ; j1 = 1 ; j2 = 3
      else if(j==3) then ; j1 = 1 ; j2 = 2
      endif
      determ = imat(i1,j1)*imat(i2,j2)-imat(i1,j2)*imat(i2,j1)
      end function determ

c     --------------------------

      end subroutine determine_imat

c     --------------------------

      end module m_fft


