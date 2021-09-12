c--------------------------------------------------------------------------------
c Copyright (c) 2021 Peter Grünberg Institut, Forschungszentrum Jülich, Germany
c This file is part of SPEX and available as free software under the conditions
c of the MIT license as expressed in the LICENSE file in more detail.
c--------------------------------------------------------------------------------

# include "cppmacro.h"

      subroutine vprod(a,b,c)
      use, intrinsic :: iso_fortran_env
      implicit none

      real_dp, intent(out) :: a(3)
      real_dp, intent(in)  :: b(3),c(3)

      a(1) = b(2)*c(3) - b(3)*c(2)
      a(2) = b(3)*c(1) - b(1)*c(3)
      a(3) = b(1)*c(2) - b(2)*c(1)

      end

c     ---------------

      function gptnorm(gpt)
      use global, only: rlat

      use, intrinsic :: iso_fortran_env
      implicit none
      real_dp             :: gptnorm
      integer, intent(in) :: gpt(3)

      gptnorm = sqrt ( sum ( matmul ( rlat(:,:),gpt(:) ) **2 ) )

      end

c     ---------------

      subroutine getangles(theta,phi,r1)
      use, intrinsic :: iso_fortran_env
      implicit none
      real_dp theta,phi,r1(3),r(3),rdum,pi
      pi = 4*atan(1d0)
      r = r1
      rdum = sqrt(sum(r**2))
      if(rdum<1d-12) then
        theta = 0
        phi   = 0
      else
        r = r/rdum
        theta = acos(r(3))
        if(theta<1d-12) then
          phi = 0
        else
          rdum = r(1)/sin(theta)
          if(abs(rdum)>1) rdum = 1d0 ! obviously necessary for G95 !?
          phi = acos(rdum)
          if(r(2)/sin(theta)<0) phi = 2*pi-phi
        endif
      endif
      end

c     ---------------

      function determinant(mat)
      use, intrinsic :: iso_fortran_env
      implicit none
      real_dp             :: determinant
      real_dp, intent(in) :: mat(3,3)
      determinant = mat(1,1) * (mat(2,2)*mat(3,3)-mat(3,2)*mat(2,3)) +
     &              mat(2,1) * (mat(3,2)*mat(1,3)-mat(1,2)*mat(3,3)) +
     &              mat(3,1) * (mat(1,2)*mat(2,3)-mat(2,2)*mat(1,3))
      end

      function cdeterminant(mat)
      use, intrinsic :: iso_fortran_env
      implicit none
      complex_dp             :: cdeterminant
      complex_dp, intent(in) :: mat(3,3)
      cdeterminant = mat(1,1) * (mat(2,2)*mat(3,3)-mat(3,2)*mat(2,3)) +
     &               mat(2,1) * (mat(3,2)*mat(1,3)-mat(1,2)*mat(3,3)) +
     &               mat(3,1) * (mat(1,2)*mat(2,3)-mat(2,2)*mat(1,3))
      end

c     ---------------

      function inverse3(mat)
      use, intrinsic :: iso_fortran_env
      implicit none
      real_dp             :: inverse3(3,3)
      real_dp, intent(in) :: mat(3,3)
      real_dp             :: determinant
      inverse3(1,1) = mat(2,2)*mat(3,3) - mat(3,2)*mat(2,3)
      inverse3(1,2) = mat(3,2)*mat(1,3) - mat(1,2)*mat(3,3)
      inverse3(1,3) = mat(1,2)*mat(2,3) - mat(2,2)*mat(1,3)
      inverse3(2,1) = mat(2,3)*mat(3,1) - mat(3,3)*mat(2,1)
      inverse3(2,2) = mat(3,3)*mat(1,1) - mat(1,3)*mat(3,1)
      inverse3(2,3) = mat(1,3)*mat(2,1) - mat(2,3)*mat(1,1)
      inverse3(3,1) = mat(2,1)*mat(3,2) - mat(3,1)*mat(2,2)
      inverse3(3,2) = mat(3,1)*mat(1,2) - mat(1,1)*mat(3,2)
      inverse3(3,3) = mat(1,1)*mat(2,2) - mat(2,1)*mat(1,2)
      inverse3      = inverse3 / determinant(mat)
      end

c     ---------------

c     Returns the k-point index of the sum kpt(:,ikpt1)+kpt(:,ikpt2)+G where G is a reciprocal lattice vector.
      function kptsum(ikpt1,ikpt2)
      use global, only: kpt,nkpt3,pkpt,pkpt1,modulo1,modulo1r,nkpt,kptadd
      use, intrinsic :: iso_fortran_env
      implicit none
      integer              :: kptsum
      integer, intent(in)  :: ikpt1,ikpt2
      integer              :: g(3)
      real_dp              :: kvec(3)
      if(ikpt1>nkpt.and.ikpt2>nkpt) Bug('both k-point indices of shifted set.')
      if(ikpt1<=nkpt.and.ikpt2<=nkpt) then
        kvec   = modulo1(kpt(:,ikpt1)+kpt(:,ikpt2)) ! = kpt(:,ikpt1)+kpt(:,ikpt2)+G
        g      = nint( kvec * nkpt3 )
        kptsum = pkpt(g(1),g(2),g(3))
      else
        kvec   = modulo1r(kpt(:,ikpt1)+kpt(:,ikpt2)-kptadd)
        g      = nint( kvec * nkpt3 )
        kptsum = pkpt1(g(1),g(2),g(3))
      endif
      end

c     Returns the k-point index of the difference kpt(:,ikpt1)-kpt(:,ikpt2)+G where G is a reciprocal lattice vector.
      function kptdiff(ikpt1,ikpt2)
      use global, only: kpt,nkpt3,pkpt,pkpt1,modulo1,modulo1r,nkpt,kptadd
      use, intrinsic :: iso_fortran_env
      implicit none
      integer              :: kptdiff
      integer, intent(in)  :: ikpt1,ikpt2
      integer              :: g(3)
      real_dp              :: kvec(3)
      if(ikpt1<=nkpt.and.ikpt2>nkpt) Bug('First/second k-point index of normal/shifted set.')
      if(ikpt1>nkpt.and.ikpt2<=nkpt) then
        kvec    = modulo1(kpt(:,ikpt1)-kpt(:,ikpt2)-kptadd)
        g       = nint( kvec * nkpt3 )
        kptdiff = pkpt1(g(1),g(2),g(3))
      else
        kvec    = modulo1r(kpt(:,ikpt1)-kpt(:,ikpt2)) ! = kpt(:,ikpt1)-kpt(:,ikpt2)+G
        g       = nint( kvec * nkpt3 )
        kptdiff = pkpt(g(1),g(2),g(3))
      endif
      end

c     ---------------

c     Returns coefficients {c_lm} of
c
c     vec * r = SUM c   Y  (r)  ;  |r| = 1
c                 m  1m  1m
c
      subroutine vector2harmonics(coeff,vec)
      use global, only: pi,img
      use, intrinsic :: iso_fortran_env
      implicit none
      complex_dp, intent(out) :: coeff(3)
      complex_dp, intent(in)  :: vec(3)
      coeff(1)  = sqrt(2*pi/3)  * (   vec(1) + img*vec(2) )
      coeff(2)  = sqrt(4*pi/3)  *     vec(3)
      coeff(3)  = sqrt(2*pi/3)  * ( - vec(1) + img*vec(2) )
      end

c     Returns coefficients {c_lm} of
c
c      T
c     r * mat * r = c   Y  (r) + SUM c   Y  (r)  ;  |r| = 1
c                    00  00        m  2m  2m
c
c     Note that mat is symmetrized (->mat1) because r^T*mat*r = r^T*mat1*r.
      subroutine matrix2harmonics(coeff0,coeff2,mat)
      use global, only: pi,img
      use, intrinsic :: iso_fortran_env
      implicit none
      complex_dp, intent(out) :: coeff0,coeff2(5)
      complex_dp, intent(in)  :: mat(3,3)
      complex_dp              :: mat1(3,3)
c      if(abs(mat(1,2)-mat(2,1))+
c     &   abs(mat(1,3)-mat(3,1))+
c     &   abs(mat(2,3)-mat(3,2))>1d-8) then
c        write(6,'(A,3(/3F20.10))') 'matrix2harmonics: input matrix not symmetric:',mat
c        Error('input matrix not symmetric.')
c      endif
      mat1      = ( mat + transpose(mat) ) / 2
      coeff0    = sqrt(4*pi/9)  * (   mat1(1,1) + mat1(2,2) + mat1(3,3) )
      coeff2(1) = sqrt(2*pi/15) * (   mat1(1,1) - mat1(2,2) + 2*img * mat1(1,2) )
      coeff2(2) = sqrt(8*pi/15) * (   mat1(1,3) + img*mat1(2,3) )
      coeff2(3) = sqrt(4*pi/45) * ( 2*mat1(3,3) - mat1(1,1) - mat1(2,2) )
      coeff2(4) = sqrt(8*pi/15) * ( - mat1(1,3) + img*mat1(2,3) )
      coeff2(5) = sqrt(2*pi/15) * (   mat1(1,1) - mat1(2,2) - 2*img * mat1(1,2) )
      end

c     Transforms coefficients coeff to coeff' such that
c                                       *
c     SUM coeff   Y  (r) = SUM coeff'  Y  (r) .
c      lm      lm  lm       lm      lm  lm
      subroutine conjgharmonics(coeff,ll)
      use, intrinsic :: iso_fortran_env
      implicit none
      integer,    intent(in)    :: ll
      complex_dp, intent(inout) :: coeff((ll+1)**2)
      complex_dp                :: coeff1(-ll:ll)
      integer                   :: l,m,lm
      lm = 0
      do l = 0,ll
        coeff1(-l:l) = coeff(lm+1:lm+2*l+1)
        do m = -l,l
          lm        = lm + 1
          coeff(lm) = coeff1(-m) * (-1)**m
        enddo
      enddo
      end

c     Inverse operation of matrix2harmonics (see above). The returned matrix is symmetric.
      subroutine harmonics2matrix(mat,coeff0,coeff2)
      use global, only: pi,img
      use, intrinsic :: iso_fortran_env
      implicit none
      complex_dp, intent(out) :: mat(3,3)
      complex_dp, intent(in)  :: coeff0,coeff2(5)
      real_dp                 :: f00,f20,f21,f22
      f00      = sqrt( 9/(4*pi))
      f20      = sqrt(45/(4*pi))
      f21      = sqrt(15/(8*pi))
      f22      = sqrt(15/(2*pi))
      mat(1,1) = f00 * coeff0    !   Mxx + Myy +  Mzz
      mat(2,1) = f20 * coeff2(3) ! - Mxx - Myy + 2Mzz
      mat(3,3) = ( mat(1,1) + mat(2,1) ) / 3       ! Mzz
      mat(2,2) =   mat(1,1) - mat(3,3)             ! Mxx + Myy
      mat(2,1) = f22 * coeff2(1) !   Mxx - Myy + 2iMxy
      mat(3,1) = f22 * coeff2(5) !   Mxx - Myy - 2iMxy
      mat(3,1) = ( mat(2,1) + mat(3,1) ) / 2       ! Mxx - Myy
      mat(1,1) = ( mat(2,2) + mat(3,1) ) / 2       ! Mxx
      mat(2,2) = ( mat(2,2) - mat(3,1) ) / 2       ! Myy
      mat(1,2) = ( mat(2,1) - mat(3,1) ) / (2*img) ! Mxy
      mat(2,1) = f21 * coeff2(2) !   Mxz + iMyz
      mat(3,1) = f21 * coeff2(4) ! - Mxz + iMyz
      mat(1,3) = ( mat(2,1) - mat(3,1) ) / 2       ! Mxz
      mat(2,3) = ( mat(2,1) + mat(3,1) ) / (2*img) ! Myz
      mat(2,1) = mat(1,2)
      mat(3,1) = mat(1,3)
      mat(3,2) = mat(2,3)
      end

c     ---------------

c     The routines convert SUM(m) a(m) complexY(lm) = SUM(m) b(m) realY(lm)
c
      subroutine harmonics_r2c(a,b,l)
      use global, only: img
      use, intrinsic :: iso_fortran_env
      implicit none
      integer,    intent(in)  :: l
      complex_dp, intent(out) :: a(-l:l)
      complex_dp, intent(in)  :: b(-l:l)
      real_dp,    parameter   :: sq = 1/sqrt(2d0)
      integer                 :: m
      a = 0
      do m = -l,-1 ! m < 0
        a( m) = a( m) - sq * b(m) * (-1)**m * img
        a(-m) = a(-m) + sq * b(m)           * img
      enddo
      a(0) = b(0)  ! m = 0
      do m = 1,l   ! m > 0
        a( m) = a( m) + sq * b(m) * (-1)**m
        a(-m) = a(-m) + sq * b(m)
      enddo
      end
      subroutine harmonics_c2r(b,a,l)
      use global, only: img
      use, intrinsic :: iso_fortran_env
      implicit none
      integer,    intent(in)  :: l
      complex_dp, intent(in)  :: a(-l:l)
      complex_dp, intent(out) :: b(-l:l)
      real_dp,    parameter   :: sq = 1/sqrt(2d0)
      integer                 :: m
      b = 0
      do m = -l,-1 ! m < 0
        b( m) = b( m) + sq * a(m) * (-1)**m * img
        b(-m) = b(-m) + sq * a(m)
      enddo
      b(0) = a(0)  ! m = 0
      do m = 1,l   ! m > 0
        b( m) = b( m) + sq * a(m) * (-1)**m
        b(-m) = b(-m) - sq * a(m)           * img
      enddo
      end

c     ---------------

c     Returns the rotation matrix (theta around axis)
      subroutine rotation_axis(rot,axis,theta)
      use, intrinsic :: iso_fortran_env
      implicit none
      real_dp, intent(out) :: rot(3,3)
      real_dp, intent(in)  :: axis(3),theta
      real_dp              :: u(3),rdum,ct,st
      rdum     = sqrt(sum(axis**2)) ; if(rdum<1d-15) Error('axis is null vector.')
      u        = axis / rdum
      ct       = cos(theta)
      st       = sin(theta)
      rot(1,1) = u(1)*u(1)*(1-ct) + ct
      rot(2,1) = u(1)*u(2)*(1-ct) + u(3)*st
      rot(3,1) = u(1)*u(3)*(1-ct) - u(2)*st
      rot(1,2) = u(1)*u(2)*(1-ct) - u(3)*st
      rot(2,2) = u(2)*u(2)*(1-ct) + ct
      rot(3,2) = u(2)*u(3)*(1-ct) + u(1)*st
      rot(1,3) = u(1)*u(3)*(1-ct) + u(2)*st
      rot(2,3) = u(2)*u(3)*(1-ct) - u(1)*st
      rot(3,3) = u(3)*u(3)*(1-ct) + ct
      end

c     ---------------

c     Returns the rotation matrix (Euler angles a(1:3), x convention)
      subroutine rotation_euler(rot,a)
      use, intrinsic :: iso_fortran_env
      implicit none
      real_dp, intent(out) :: rot(3,3)
      real_dp, intent(in)  :: a(3)
      real_dp              :: ca,sa,cb,sb,cc,sc
      ca       = cos(a(1)) ; cb = cos(a(2)) ; cc = cos(a(3))
      sa       = sin(a(1)) ; sb = sin(a(2)) ; sc = sin(a(3))
      rot(1,1) =   ca*cc - sa*cb*sc
      rot(2,1) = - ca*sc - sa*cb*cc
      rot(3,1) =           sa*sb
      rot(1,2) =   sa*cc + ca*cb*sc
      rot(2,2) = - sa*sc + ca*cb*cc
      rot(3,2) =         - ca*sb
      rot(1,3) =   sb*sc
      rot(2,3) =   sb*cc
      rot(3,3) =      cb
      end

c     ---------------

c     Returns the rotation matrix given by the polar angle theta and azimuthal angle phi
      subroutine rotation_sph(rot,theta,phi)
      use, intrinsic :: iso_fortran_env
      implicit none
      real_dp, intent(out) :: rot(3,3)
      real_dp, intent(in)  :: theta,phi
      real_dp              :: ct,st,cp,sp
      ct       = cos(theta) ; cp = cos(phi)
      st       = sin(theta) ; sp = sin(phi)
      rot(1,1) =  ct*cp
      rot(2,1) =  ct*sp
      rot(3,1) = -st
      rot(1,2) = -sp
      rot(2,2) =  cp
      rot(3,2) =  0
      rot(1,3) =  st*cp
      rot(2,3) =  st*sp
      rot(3,3) =  ct
      end

c     ---------------

c     Returns backfolded vectors in mat in terms of lat (minimization of length).
c     lat and mat are in cartesian coordinates.
      subroutine backfold(mat,lat)
      use, intrinsic :: iso_fortran_env
      implicit none
      real_dp, intent(out) :: mat(3,3)
      real_dp, intent(in)  :: lat(3,3)
      real_dp              :: vec(3),vol
      integer              :: pnt(3),i,j
      call vprod(vec,lat(:,1),lat(:,2))
      vol = dot_product(vec,lat(:,3))
      mat = lat
      vec = [(sum(mat(:,i)**2),i=1,3)] ; call rorderp(pnt,vec,3) ; pnt(3:1:-1) = pnt
      mat = mat(:,pnt)
      do i = 1,2
        do j = i+1,3
          do
            vec = mat(:,i) + mat(:,j) * sign(1d0,-dot_product(mat(:,i),lat(:,j)))
            if(sum(vec**2)>sum(mat(:,i)**2)-1d-12) exit
            mat(:,i) = vec
          enddo
        enddo
      enddo
      mat(:,pnt) = mat
      call vprod(vec,mat(:,1),mat(:,2))
      if(abs(dot_product(vec,mat(:,3))-vol)/vol>1d-12) Bug('Volume changed.')
      end

c     ---------------

      subroutine lattice_vectors(gpt,ngpt,kvec,gcut,rlat)
      use, intrinsic :: iso_fortran_env
      implicit none
      integer, intent(out) :: ngpt,gpt(3,*)
      real_dp, intent(in)  :: kvec(3),gcut,rlat(3,3)
      real_dp              :: kg(3),kgn
      integer              :: g(3),ix,iy,iz,i,n
      logical              :: found
      i     = 0
      n     = 0
      found = .true.
      do while(found)
        found = .false.
        do ix = -n,n
          do iy = -(n-abs(ix)),n-abs(ix)
            iz = n - abs(ix) - abs(iy)
            do
              g   = [ ix,iy,iz ] - nint(kvec)
              kg  = matmul(rlat,kvec+g)
              kgn = sum(kg**2)
              if(kgn<=gcut**2) then
                i     = i + 1
                found = .true.
                if(ngpt>0) gpt(:,i) = g
              endif
              if(iz<=0) exit
              iz = -iz
            enddo
          enddo
        enddo
        n = n + 1
      enddo      
      if(ngpt==0)     then ; ngpt = i
      else if(i>ngpt) then ; Bug('Dimension of gpt too small.')
      endif
      end

      function n_lattice_vectors(kvec,gcut,rlat) result(ngpt)      
      use, intrinsic :: iso_fortran_env
      implicit none
      integer             :: ngpt
      real_dp, intent(in) :: kvec(3),gcut,rlat(3,3)
      integer             :: idum(3)
      ngpt = 0
      call lattice_vectors(idum,ngpt,kvec,gcut,rlat)
      end
      
