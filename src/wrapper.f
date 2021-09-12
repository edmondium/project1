c--------------------------------------------------------------------------------
c Copyright (c) 2021 Peter Grünberg Institut, Forschungszentrum Jülich, Germany
c This file is part of SPEX and available as free software under the conditions
c of the MIT license as expressed in the LICENSE file in more detail.
c--------------------------------------------------------------------------------

# include "cppmacro.h"

# ifndef M_PART
#   define M_PART 1
#   include __FILE__
#   undef M_PART
#   define M_PART 2
#   include __FILE__
#   undef M_PART
#   define M_PART 3
# endif

c Matrix multiplications (matmat) can be additionally compiled as subroutines (p_matmat):
c# define P_MATMAT

c -----------------------------
# if M_PART == 1
c -----------------------------
      
      module wrapper

      use, intrinsic :: iso_fortran_env

# ifdef MPI
      integer, parameter :: blk = 64 ! ScaLAPACK block size
# endif

      interface diagonal
      module procedure diagonal_d,diagonal_z
      end interface

      interface blockmat
      module procedure  blockmat_d,blockmat_z
      end interface

      interface packmat
      module procedure  packmat_d,packmat_z
      end interface

      interface unpackmat
      module procedure  unpackmat_d,unpackmat_z
      end interface

      interface p_packmat
      module procedure  p_packmat_d,p_packmat_z
      end interface

      interface p_unpackmat
      module procedure  p_unpackmat_d,p_unpackmat_z
      end interface

      interface unpackrow
      module procedure  unpackrow_z,unpackrow_zm
      end interface

      interface dotprod
      module procedure  dotprod_dd, dotprod_dz, dotprod_zd, dotprod_zz
      end interface
      interface matvec
      module procedure  matvec_dd,  matvec_zd,  matvec_zz,
     &                  matvec_dpd, matvec_dpz, matvec_zpd, matvec_zpz
      end interface
      interface macvec
      module procedure  macvec_dd,  macvec_zd,  macvec_zz
      end interface

      interface matmat
      module procedure  matmat_dpdp, matmat_dpzp, matmat_zpdp, matmat_zpzp,
     &                  matmat_dpdm, matmat_dmdp, matmat_dmdm,
     &                  matmat_dmzp, matmat_dpzm, matmat_dmzm,
     &                  matmat_zmdp, matmat_zpdm, matmat_zmdm,
     &                  matmat_zmzp, matmat_zpzm, matmat_zmzm
      end interface

      interface macmat
      module procedure  macmat_dmdm, macmat_dmzm, macmat_zmdm, macmat_zmzm
      end interface

      interface matmac
      module procedure  matmac_dmdm, matmac_dmzm, matmac_zmdm, matmac_zmzm
      end interface

# ifdef P_MATMAT

      interface p_matmat
      module procedure  p_matmat_dpdp, p_matmat_dpzp, p_matmat_zpdp, p_matmat_zpzp,
     &                  p_matmat_dpdm, p_matmat_dmdp, p_matmat_dmdm,
     &                  p_matmat_dmzp, p_matmat_dpzm, p_matmat_dmzm,
     &                  p_matmat_zmdp, p_matmat_zpdm, p_matmat_zmdm,
     &                  p_matmat_zmzp, p_matmat_zpzm, p_matmat_zmzm
      end interface

      interface p_macmat
      module procedure  p_macmat_dmdm, p_macmat_dmzm, p_macmat_zmdm, p_macmat_zmzm
      end interface

      interface p_matmac
      module procedure  p_matmac_dmdm, p_matmac_dmzm, p_matmac_zmdm, p_matmac_zmzm
      end interface

# endif

      interface diagmat
      module procedure  diagmat_dd,  diagmat_dz,  diagmat_zd,  diagmat_zz
      end interface diagmat

      interface matdiag
      module procedure  matdiag_dd,  matdiag_dz,  matdiag_zd,  matdiag_zz      
      end interface matdiag

      interface diagmatdiag
      module procedure  diagmatdiag_ddd, diagmatdiag_dzd, diagmatdiag_ddpd, diagmatdiag_dzpd
      end interface
      
      interface p_diagmatdiag
      module procedure  p_diagmatdiag_ddd,  p_diagmatdiag_dzd,  p_diagmatdiag_ddpd,  p_diagmatdiag_dzpd,
     &                  p_diagmatdiag_ddd1, p_diagmatdiag_dzd1, p_diagmatdiag_ddpd1, p_diagmatdiag_dzpd1
      end interface

      interface p_plusdiag
      module procedure  p_plusdiag_d,  p_plusdiag_z,  p_plusdiag_dp,  p_plusdiag_zp,
     &                  p_plusdiag1_d, p_plusdiag1_z, p_plusdiag1_dp, p_plusdiag1_zp
      end interface

      interface unitarytrafo
      module procedure unitarytrafo_d,  unitarytrafo_z,   unitarytrafo_dp,  unitarytrafo_zp,
     &                 unitarytrafo_d1, unitarytrafo_z1,  unitarytrafo_dp1, unitarytrafo_zp1,
     &                 unitarytrafo_zd, unitarytrafo_zd1, unitarytrafo_dp2, unitarytrafo_zp2
      end interface

      interface diagonalize ! the routines ending on 's' should be reformulated using '...x' once
      module procedure  diagonalize_de,  diagonalize_dv,  diagonalize_dpe,  diagonalize_dpv,
     &                  diagonalize_ze,  diagonalize_zv,  diagonalize_zpe,  diagonalize_zpv,
     &                  diagonalize_deo, diagonalize_dvo, diagonalize_dpeo, diagonalize_dpvo,
     &                  diagonalize_zeo, diagonalize_zvo, diagonalize_zpeo, diagonalize_zpvo,
     &                  diagonalize_dvs, diagonalize_dvos,diagonalize_dpvs, diagonalize_dpvos,
     &                  diagonalize_zvs, diagonalize_zvos,diagonalize_zpvs, diagonalize_zpvos,
     &                  diagonalize_dvx, diagonalize_dvox,diagonalize_dpvx, diagonalize_dpvox,
     &                  diagonalize_zvx, diagonalize_zvox,diagonalize_zpvx, diagonalize_zpvox
      end interface

      interface diagonalize_gen
      module procedure  diagonalize_gen_e, diagonalize_gen_v
      end interface

      interface geteigen
      module procedure  geteigen_zpvo
      end interface

      interface inverse
      module procedure  inverse_d,  inverse_dp,  inverse_z,  inverse_zp,
     &                  inverse_d1, inverse_dp1, inverse_z1, inverse_zp1
      end interface

      interface invert
      module procedure  invert_d, invert_dp, invert_z, invert_zp
      end interface

      interface sqrtmat
      module procedure  sqrtmat_d,  sqrtmat_dp,  sqrtmat_z,  sqrtmat_zp,
     &                  sqrtmat_d1, sqrtmat_dp1, sqrtmat_z1, sqrtmat_zp1
      end interface

      interface sqrtm
      module procedure  sqrtm_d,  sqrtm_dp,  sqrtm_z,  sqrtm_zp
      end interface

      interface solve
      module procedure  solve_d, solve_z
      end interface

      interface trace
      module procedure  trace_d, trace_z, trace_dp, trace_zp
      end interface

# ifdef MPI
      interface Mdiagonalize
      module procedure  Mdiagonalize_d,    Mdiagonalize_z,
     &                  Mdiagonalize_dp,   Mdiagonalize_zp,
     &                  Mdiagonalize_dv,   Mdiagonalize_zv,
     &                  Mdiagonalize_dpv,  Mdiagonalize_zpv,
     &                  Mdiagonalize_do,   Mdiagonalize_zo,
     &                  Mdiagonalize_dpo,  Mdiagonalize_zpo,
     &                  Mdiagonalize_dvo,  Mdiagonalize_zvo,
     &                  Mdiagonalize_dpvo, Mdiagonalize_zpvo
      end interface

      private :: Mdiagonalize0_d, Mdiagonalize0_dx, Mdiagonalize0_dox,
     &           Mdiagonalize0_z, Mdiagonalize0_zx, Mdiagonalize0_zox

      interface Minverse
      module procedure  Minverse_d, Minverse_dp, Minverse_z, Minverse_zp
      end interface

      private :: Minverse0_d, Minverse0_z
# endif

      private :: check_dim, which_ranges

      contains

c     --------

      function identity(n)
      implicit none
      integer, intent(in) :: n
      integer             :: identity(n,n)
      integer             :: i
      identity = 0
      do i = 1,n
        identity(i,i) = 1
      enddo
      end function identity

      function identityp(n)
      implicit none
      integer, intent(in) :: n
      integer             :: identityp(n*(n+1)/2)
      integer             :: i,j
      identityp = 0
      j         = 0
      do i = 1,n
        j            = j + i
        identityp(j) = 1
      enddo
      end function identityp

c     --------

      function diagonal_d(d)
      implicit none
      real_dp, intent(in) :: d(:)
      real_dp             :: diagonal_d(size(d),size(d))
      integer             :: i
      diagonal_d = 0
      do i = 1,size(d)
        diagonal_d(i,i) = d(i)
      enddo
      end function diagonal_d

      function diagonal_z(d)
      implicit none
      complex_dp, intent(in) :: d(:)
      complex_dp             :: diagonal_z(size(d),size(d))
      integer                :: i
      diagonal_z = 0
      do i = 1,size(d)
        diagonal_z(i,i) = d(i)
      enddo
      end function diagonal_z

      function diagonalp(d)
      implicit none
      real_dp, intent(in) :: d(:)
      real_dp             :: diagonalp(size(d)*(size(d)+1)/2)
      integer             :: i,j
      diagonalp = 0
      j         = 0
      do i = 1,size(d)
        j            = j + i
        diagonalp(j) = d(i)
      enddo
      end function diagonalp

c     --------

      function blockmat_d(a,b)
      implicit none
      real_dp, intent(in) :: a(:,:),b(:,:)
      real_dp             :: blockmat_d(size(a,1)+size(b,1),size(a,1)+size(b,1))
      integer             :: na,nb
      na                      = size(a,1) ; call check_dim(shape(a),[na,na],'a',__LINE__)
      nb                      = size(b,1) ; call check_dim(shape(b),[nb,nb],'b',__LINE__)
      blockmat_d              = 0d0
      blockmat_d(  :na,  :na) = a
      blockmat_d(na+1:,na+1:) = b
      end function blockmat_d

      function blockmat_z(a,b)
      implicit none
      complex_dp, intent(in) :: a(:,:),b(:,:)
      complex_dp             :: blockmat_z(size(a,1)+size(b,1),size(a,1)+size(b,1))
      integer                :: na,nb
      na                      = size(a,1) ; call check_dim(shape(a),[na,na],'a',__LINE__)
      nb                      = size(b,1) ; call check_dim(shape(b),[nb,nb],'b',__LINE__)
      blockmat_z              = 0d0
      blockmat_z(  :na,  :na) = a
      blockmat_z(na+1:,na+1:) = b
      end function blockmat_z

c     --------

# include "w_packmat.inc"

c     --------
      
      function unpackrow_z(mat,m)
      implicit none
      integer,    intent(in) :: m
      complex_dp, intent(in) :: mat(:)
      complex_dp             :: unpackrow_z(nint(sqrt(0.25d0+2*size(mat))-0.5d0))
      integer                :: i,j,n
      n               = size(unpackrow_z)
      i               = m*(m-1)/2
      unpackrow_z(:m) = conjg(mat(i+1:i+m))
      i               = i + m
      do j = m+1,n
        i              = i + j - 1
        unpackrow_z(j) = mat(i)
      enddo
      end function unpackrow_z

      function unpackrow_zm(mat,m)
      implicit none
      integer,    intent(in) :: m
      complex_dp, intent(in) :: mat(:,:)
      complex_dp             :: unpackrow_zm(nint(sqrt(0.25d0+2*size(mat,1))-0.5d0),size(mat,2))
      integer                :: i,j,n
      n                  = size(unpackrow_zm,1)
      i                  = m*(m-1)/2
      unpackrow_zm(:m,:) = conjg(mat(i+1:i+m,:))
      i                  = i + m
      do j = m+1,n
        i                 = i + j - 1
        unpackrow_zm(j,:) = mat(i,:)
      enddo
      end function unpackrow_zm

c     --------

      function dotprod_dd(vec1,vec2)
      implicit none
      real_dp, intent(in) :: vec1(:),vec2(:)
      real_dp             :: dotprod_dd
      integer             :: n
      real_dp             :: ddot
      n          = size(vec1) ; call check_dim(shape(vec2),[n],'vec2',__LINE__)
      dotprod_dd = ddot(n,vec1,1,vec2,1)
      end function dotprod_dd

      function dotprod_dz(vec1,vec2)
      implicit none
      real_dp,    intent(in) :: vec1(:)
      complex_dp, intent(in) :: vec2(:)
      complex_dp             :: dotprod_dz
      integer                :: n
      real_dp                :: ddot
      n          = size(vec1) ; call check_dim(shape(vec2),[n],'vec2',__LINE__)
      dotprod_dz = ddot(n,vec1,1,real(vec2),1)+(0d0,1d0)*ddot(n,vec1,1,imag(vec2),1)
      end function dotprod_dz

      function dotprod_zd(vec1,vec2)
      implicit none
      complex_dp, intent(in) :: vec1(:)
      real_dp,    intent(in) :: vec2(:)
      complex_dp             :: dotprod_zd
      integer                :: n
      real_dp                :: ddot
      n          = size(vec1) ; call check_dim(shape(vec2),[n],'vec2',__LINE__)
      dotprod_zd = ddot(n,real(vec1),1,vec2,1)-(0d0,1d0)*ddot(n,imag(vec1),1,vec2,1)
      end function dotprod_zd

      function dotprod_zz(vec1,vec2)
      implicit none
      complex_dp, intent(in) :: vec1(:),vec2(:)
      complex_dp             :: dotprod_zz
      integer                :: n
      complex_dp             :: zdotc
      n          = size(vec1) ; call check_dim(shape(vec2),[n],'vec2',__LINE__)
      dotprod_zz = zdotc(n,vec1,1,vec2,1)
      end function dotprod_zz

c     --------

      function matvec_dd(mat,vec)
      implicit none
      real_dp, intent(in) :: vec(:)
      real_dp, intent(in) :: mat(:,:)
      real_dp             :: matvec_dd(size(mat,1))
      integer             :: n1,n2
      n1 = size(mat,1)
      n2 = size(mat,2) ; call check_dim(shape(vec),[n2],'vec',__LINE__) 
      call dgemv('N',n1,n2,1d0,mat,n1,vec,1,0d0,matvec_dd,1)
      end function matvec_dd

      function matvec_zd(mat,vec)
      implicit none
      real_dp,    intent(in) :: vec(:)
      complex_dp, intent(in) :: mat(:,:)
      complex_dp             :: matvec_zd(size(mat,1))
      integer                :: n1,n2
      n1 = size(mat,1)
      n2 = size(mat,2) ; call check_dim(shape(vec),[n2],'vec',__LINE__)
      call dgemv('N',2*n1,n2,1d0,mat,2*n1,vec,1,0d0,matvec_zd,1)
      end function matvec_zd

      function matvec_zz(mat,vec)
      implicit none
      complex_dp, intent(in) :: vec(:)
      complex_dp, intent(in) :: mat(:,:)
      complex_dp             :: matvec_zz(size(mat,1))
      integer                :: n1,n2
      n1 = size(mat,1)
      n2 = size(mat,2) ; call check_dim(shape(vec),[n2],'vec',__LINE__)
      call zgemv('N',n1,n2,(1d0,0d0),mat,n1,vec,1,(0d0,0d0),matvec_zz,1)
      end function matvec_zz

      function matvec_dpd(mat,vec)
      implicit none
      real_dp, intent(in)  :: mat(:),vec(:)
      real_dp              :: matvec_dpd(size(vec))
      integer              :: nn,n
      n  = size(vec)
      nn = n*(n+1)/2
      call check_dim(shape(mat),[nn],'mat',__LINE__)
      call dspmv('U',n,1d0,mat,vec,1,0d0,matvec_dpd,1)
      end function matvec_dpd

      function matvec_dpz(mat,vec)
      implicit none
      real_dp,    intent(in) :: mat(:)
      complex_dp, intent(in) :: vec(:)
      complex_dp             :: matvec_dpz(size(vec))
      integer                :: nn,n
      n  = size(vec)
      nn = n*(n+1)/2 ; call check_dim(shape(mat),[nn],'mat',__LINE__)
      call dspmv('U',n,1d0,mat,imag(vec),1,0d0,matvec_dpz,2) ; matvec_dpz = (0d0,1d0)*real(matvec_dpz)
      call dspmv('U',n,1d0,mat,real(vec),1,1d0,matvec_dpz,2)
      end function matvec_dpz

      function matvec_zpd(mat,vec)
      implicit none
      complex_dp, intent(in) :: mat(:)
      real_dp,    intent(in) :: vec(:)
      complex_dp             :: matvec_zpd(size(vec))
      matvec_zpd = matvec_zpz(mat,(1d0,0d0)*vec)
      end function matvec_zpd

      function matvec_zpz(mat,vec)
      implicit none
      complex_dp, intent(in)  :: mat(:),vec(:)
      complex_dp              :: matvec_zpz(size(vec))
      integer                 :: nn,n
      n  = size(vec)
      nn = n*(n+1)/2 ; call check_dim(shape(mat),[nn],'mat',__LINE__)
      call zhpmv('U',n,(1d0,0d0),mat,vec,1,(0d0,0d0),matvec_zpz,1)
      end function matvec_zpz

c     --------

      function triavec(mat,vec)
      implicit none
      real_dp, intent(in) :: mat(:),vec(:)
      real_dp             :: triavec(size(vec))
      integer             :: nn,n
      n       = size(vec)
      nn      = n*(n+1)/2 ; call check_dim(shape(mat),[nn],'mat',__LINE__)
      triavec = vec
      call dtpmv('U','N','N',n,mat,triavec,1)
      end function triavec

c     --------

c ------------------------------
# elif M_PART == 2
c ------------------------------

# ifndef INLOOP
#   define INLOOP
#   define SUBROUT_FUNC function
#   define MATOUT
#   define INTENTOUT
#   define RESULTMAT result(matout)
#   ifdef P_MATMAT
#     include __FILE__
#     define SUBROUT_FUNC subroutine p_
#     define MATOUT       matout,
#     define INTENTOUT   ,intent(out)
#     define RESULTMAT
#   endif
# endif
      
      SUBROUT_FUNC matmat_dpdp( MATOUT mat1,mat2) RESULTMAT
      implicit none
      real_dp, intent(in)  :: mat1(:),mat2(:)
      real_dp  INTENTOUT   :: matout(nint(sqrt(0.25d0+2*size(mat1))-0.5d0),
     &                               nint(sqrt(0.25d0+2*size(mat1))-0.5d0))
      real_dp, allocatable :: vec(:),vec2(:)
      integer              :: nn,n,k1,i,j,k
      nn = size(mat1)
      n  = nint(sqrt(0.25d0+2*nn)-0.5d0) ; allocate ( vec(n),vec2(n) )
      if(n*(n+1)-2*nn/=0) Bug('Input matrix has wrong size.')
      call check_dim(shape(mat2),[nn],'mat2',__LINE__)
      k  = 0
      do i = 1,n
        vec2(:i) = mat2(k+1:k+i)
        k1       = k+2*i
        do j = i+1,n
          vec2(j) = mat2(k1)
          k1      = k1 + j
        enddo
        call dspmv('U',n,1d0,mat1,vec2,1,0d0,vec,1)
        matout(:,i) = vec
        k = k + i
      enddo
      deallocate ( vec,vec2 )
      end SUBROUT_FUNC matmat_dpdp

      SUBROUT_FUNC matmat_dpzp( MATOUT mat1,mat2) RESULTMAT
      implicit none
      real_dp,    intent(in)  :: mat1(:)
      complex_dp, intent(in)  :: mat2(:)
      complex_dp  INTENTOUT   :: matout(nint(sqrt(0.25d0+2*size(mat1))-0.5d0),
     &                                  nint(sqrt(0.25d0+2*size(mat1))-0.5d0))
      real_dp,    allocatable :: vecr(:),veci(:)
      complex_dp, allocatable :: vec2(:)
      integer                 :: nn,n,k1,i,j,k
      nn = size(mat1)
      n  = nint(sqrt(0.25d0+2*nn)-0.5d0) ; allocate ( vecr(n),veci(n),vec2(n) )
      if(n*(n+1)-2*nn/=0) Bug('Input matrix has wrong size.')
      call check_dim(shape(mat2),[nn],'mat2',__LINE__)
      k  = 0
      do i = 1,n
        vec2(:i) = mat2(k+1:k+i)
        k1       = k+2*i
        do j = i+1,n
          vec2(j) = conjg(mat2(k1))
          k1      = k1 + j
        enddo
        call dspmv('U',n,1d0,mat1,real(vec2),1,0d0,vecr,1)
        call dspmv('U',n,1d0,mat1,imag(vec2),1,0d0,veci,1)
        matout(:,i) = vecr + (0d0,1d0) * veci
        k = k + i
      enddo
      deallocate ( vecr,veci,vec2 )
      end SUBROUT_FUNC matmat_dpzp

      SUBROUT_FUNC matmat_zpdp( MATOUT mat1,mat2) RESULTMAT
      implicit none
      complex_dp, intent(in)  :: mat1(:)
      real_dp,    intent(in)  :: mat2(:)
      complex_dp  INTENTOUT   :: matout(nint(sqrt(0.25d0+2*size(mat1))-0.5d0),
     &                                  nint(sqrt(0.25d0+2*size(mat1))-0.5d0))
      real_dp,    allocatable :: vecr(:),veci(:)
      complex_dp, allocatable :: vec1(:)
      integer                 :: nn,n,k1,i,j,k
      nn = size(mat1)
      n  = nint(sqrt(0.25d0+2*nn)-0.5d0) ; allocate ( vecr(n),veci(n),vec1(n) )
      if(n*(n+1)-2*nn/=0) Bug('Input matrix has wrong size.')
      call check_dim(shape(mat2),[nn],'mat2',__LINE__)
      k  = 0
      do i = 1,n
        vec1(:i) = conjg(mat1(k+1:k+i))
        k1       = k+2*i
        do j = i+1,n
          vec1(j) = mat1(k1)
          k1      = k1 + j
        enddo
        call dspmv('U',n,1d0,mat2,real(vec1),1,0d0,vecr,1)
        call dspmv('U',n,1d0,mat2,imag(vec1),1,0d0,veci,1)
        matout(i,:) = vecr + (0d0,1d0) * veci
        k = k + i
      enddo
      deallocate ( vecr,veci,vec1 )
      end SUBROUT_FUNC matmat_zpdp

      SUBROUT_FUNC matmat_zpzp( MATOUT mat1,mat2) RESULTMAT
      implicit none
      complex_dp, intent(in)  :: mat1(:),mat2(:)
      complex_dp  INTENTOUT   :: matout(nint(sqrt(0.25d0+2*size(mat1))-0.5d0),
     &                                  nint(sqrt(0.25d0+2*size(mat1))-0.5d0))
      complex_dp, allocatable :: vec(:),vec2(:)
      integer                 :: nn,n,k1,i,j,k
      nn = size(mat1)
      n  = nint(sqrt(0.25d0+2*nn)-0.5d0) ; allocate ( vec(n),vec2(n) )
      if(n*(n+1)-2*nn/=0) Bug('Input matrix has wrong size.')
      call check_dim(shape(mat2),[nn],'mat2',__LINE__)
      k  = 0
      do i = 1,n
        vec2(:i) = mat2(k+1:k+i)
        k1       = k+2*i
        do j = i+1,n
          vec2(j) = conjg(mat2(k1))
          k1      = k1 + j
        enddo
        call zhpmv('U',n,(1d0,0d0),mat1,vec2,1,(0d0,0d0),vec,1)
        matout(:,i) = vec
        k = k + i
      enddo
      deallocate ( vec,vec2 )
      end SUBROUT_FUNC matmat_zpzp

      SUBROUT_FUNC matmat_dpdm( MATOUT mat1,mat2) RESULTMAT
      implicit none
      real_dp, intent(in)  :: mat1(:),mat2(:,:)
      real_dp  INTENTOUT   :: matout(size(mat2,1),size(mat2,2))
      real_dp, allocatable :: vec(:),vec2(:)
      integer              :: nn,n,i,n2
      n  = size(mat2,1) ; nn = n*(n+1)/2 ; allocate ( vec(n),vec2(n) )
      n2 = size(mat2,2)
      call check_dim(shape(mat1),[nn],'mat1',__LINE__)
      do i = 1,n2
        vec2 = mat2(:,i)
        call dspmv('U',n,1d0,mat1,vec2,1,0d0,vec,1)
        matout(:,i) = vec
      enddo
      deallocate ( vec,vec2 )
      end SUBROUT_FUNC matmat_dpdm

      SUBROUT_FUNC matmat_dmdp( MATOUT mat1,mat2) RESULTMAT
      implicit none
      real_dp, intent(in)  :: mat1(:,:),mat2(:)
      real_dp  INTENTOUT   :: matout(size(mat1,1),size(mat1,2))
      real_dp, allocatable :: vec(:),vec2(:)
      integer              :: nn,n,i,n1
      n  = size(mat1,2) ; nn = n*(n+1)/2 ; allocate ( vec(n),vec2(n) )
      n1 = size(mat1,1)
      call check_dim(shape(mat2),[nn],'mat2',__LINE__)
      do i = 1,n1
        vec2 = mat1(i,:)
        call dspmv('U',n,1d0,mat2,vec2,1,0d0,vec,1)
        matout(i,:) = vec
      enddo
      deallocate ( vec,vec2 )
      end SUBROUT_FUNC matmat_dmdp

      SUBROUT_FUNC matmat_dmdm( MATOUT mat1,mat2) RESULTMAT
      implicit none
      real_dp, intent(in) :: mat1(:,:),mat2(:,:)
      real_dp  INTENTOUT  :: matout(size(mat1,1),size(mat2,2))
      integer             :: n1,n,n2
      n1 = size(mat1,1)
      n  = size(mat1,2)
      n2 = size(mat2,2)
      call check_dim(shape(mat2),[n,n2],'mat2',__LINE__)
      call dgemm('N','N',n1,n2,n,1d0,mat1,n1,mat2,n,0d0,matout,n1)
      end SUBROUT_FUNC matmat_dmdm

c      subroutine p_matmat_dmdm(matout,mat1,mat2)
c      real_dp, intent(out) :: matout(:,:)
c      real_dp, intent(in)  :: mat1(:,:),mat2(:,:)
c      integer              :: n1,n,n2
c      n1 = size(mat1,1)
c      n  = size(mat1,2)
c      n2 = size(mat2,2)
c      call check_dim(shape(matout),[n1,n2],'mat2',__LINE__)
c      call check_dim(shape(mat2),  [n, n2],'mat2',__LINE__)
c      call dgemm('N','N',n1,n2,n,1d0,mat1,n1,mat2,n,0d0,matout,n1)
c      end subroutine p_matmat_dmdm

      SUBROUT_FUNC matmat_dpzm( MATOUT mat1,mat2) RESULTMAT
      implicit none
      real_dp,    intent(in)  :: mat1(:)
      complex_dp, intent(in)  :: mat2(:,:)
      complex_dp  INTENTOUT   :: matout(size(mat2,1),size(mat2,2))
      real_dp,    allocatable :: vecr(:),veci(:)
      complex_dp, allocatable :: vec2(:)
      integer                 :: nn,n,n2,i
      n  = size(mat2,1) ; nn = n*(n+1)/2 ; allocate ( vecr(n),veci(n),vec2(n) )
      n2 = size(mat2,2)
      call check_dim(shape(mat1),[nn],'mat1',__LINE__)
      do i = 1,n2
        vec2 = mat2(:,i)
        call dspmv('U',n,1d0,mat1,real(vec2),1,0d0,vecr,1)
        call dspmv('U',n,1d0,mat1,imag(vec2),1,0d0,veci,1)
        matout(:,i) = vecr + (0d0,1d0) * veci
      enddo
      deallocate ( vecr,veci,vec2 )
      end SUBROUT_FUNC matmat_dpzm

      SUBROUT_FUNC matmat_dmzp( MATOUT mat1,mat2) RESULTMAT
      implicit none
      real_dp,    intent(in)  :: mat1(:,:)
      complex_dp, intent(in)  :: mat2(:)
      complex_dp  INTENTOUT   :: matout(size(mat1,1),size(mat1,2))
      complex_dp, allocatable :: vec1(:),vec(:)
      integer                 :: nn,n,n1,i
      n  = size(mat1,2) ; nn = n*(n+1)/2 ; allocate ( vec(n),vec1(n) )
      n1 = size(mat1,1)
      call check_dim(shape(mat2),[nn],'mat2',__LINE__)
      do i = 1,n1
        vec1 = mat1(i,:)
        call zhpmv('U',n,(1d0,0d0),mat2,vec1,1,(0d0,0d0),vec,1)
        matout(i,:) = conjg(vec)
      enddo
      deallocate ( vec,vec1 )
      end SUBROUT_FUNC matmat_dmzp

      SUBROUT_FUNC matmat_dmzm( MATOUT mat1,mat2) RESULTMAT
      implicit none
      real_dp,    intent(in) :: mat1(:,:)
      complex_dp, intent(in) :: mat2(:,:)
      complex_dp  INTENTOUT  :: matout(size(mat1,1),size(mat2,2))
      real_dp                :: matr(size(mat1,1),size(mat2,2)),mati(size(mat1,1),size(mat2,2))
      integer                :: n,n1,n2
      n1 = size(mat1,1)
      n  = size(mat1,2)
      n2 = size(mat2,2)
      call check_dim(shape(mat2),[n,n2],'mat2',__LINE__)
      call dgemm('N','N',n1,n2,n,1d0,mat1,n1,real(mat2),n,0d0,matr,n1)
      call dgemm('N','N',n1,n2,n,1d0,mat1,n1,imag(mat2),n,0d0,mati,n1)
      matout = matr + (0d0,1d0) * mati
      end SUBROUT_FUNC matmat_dmzm

      SUBROUT_FUNC matmat_zpdm( MATOUT mat1,mat2) RESULTMAT
      implicit none
      complex_dp, intent(in)  :: mat1(:)
      real_dp,    intent(in)  :: mat2(:,:)
      complex_dp  INTENTOUT   :: matout(size(mat2,1),size(mat2,2))
      complex_dp, allocatable :: vec(:),vec2(:)
      integer                 :: nn,n,i,n2
      n  = size(mat2,1) ; nn = n*(n+1)/2 ; allocate ( vec(n),vec2(n) )
      n2 = size(mat2,2)
      call check_dim(shape(mat1),[nn],'mat1',__LINE__)
      do i = 1,n2
        vec2 = mat2(:,i)
        call zhpmv('U',n,(1d0,0d0),mat1,vec2,1,(0d0,0d0),vec,1)
        matout(:,i) = vec
      enddo
      deallocate ( vec,vec2 )
      end SUBROUT_FUNC matmat_zpdm

      SUBROUT_FUNC matmat_zmdp( MATOUT mat1,mat2) RESULTMAT
      implicit none
      complex_dp, intent(in)  :: mat1(:,:)
      real_dp,    intent(in)  :: mat2(:)
      complex_dp  INTENTOUT   :: matout(size(mat1,1),size(mat1,2))
      complex_dp, allocatable :: vec1(:)
      real_dp,    allocatable :: vecr(:),veci(:)
      integer                 :: nn,n,n1,i
      n  = size(mat1,2) ; nn = n*(n+1)/2 ; allocate ( vecr(n),veci(n),vec1(n) )
      n1 = size(mat1,1)
      call check_dim(shape(mat2),[nn],'mat2',__LINE__)
      do i = 1,n1
        vec1 = conjg(mat1(i,:))
        call dspmv('U',n,1d0,mat2,real(vec1),1,0d0,vecr,1)
        call dspmv('U',n,1d0,mat2,imag(vec1),1,0d0,veci,1)
        matout(i,:) = vecr - (0d0,1d0) * veci
      enddo
      deallocate ( vecr,veci,vec1 )
      end SUBROUT_FUNC matmat_zmdp

      SUBROUT_FUNC matmat_zmdm( MATOUT mat1,mat2) RESULTMAT
      implicit none
      complex_dp, intent(in) :: mat1(:,:)
      real_dp,    intent(in) :: mat2(:,:)
      complex_dp  INTENTOUT  :: matout(size(mat1,1),size(mat2,2))
c      real_dp                :: matr(size(mat1,1),size(mat2,2)),mati(size(mat1,1),size(mat2,2))
      integer                :: n,n1,n2
      n1 = size(mat1,1)
      n  = size(mat1,2)
      n2 = size(mat2,2)
      call check_dim(shape(mat2),[n,n2],'mat2',__LINE__)
      call dgemm('N','N',2*n1,n2,n,1d0,mat1,2*n1,mat2,n,0d0,matout,2*n1)
c      call dgemm('N','N',n1,n2,n,1d0,real(mat1),n1,mat2,n,0d0,matr,n1)
c      call dgemm('N','N',n1,n2,n,1d0,imag(mat1),n1,mat2,n,0d0,mati,n1)
c      matmat_zmdm = matr + (0d0,1d0) * mati
      end SUBROUT_FUNC matmat_zmdm

      SUBROUT_FUNC matmat_zpzm( MATOUT mat1,mat2) RESULTMAT
      implicit none
      complex_dp, intent(in)  :: mat1(:),mat2(:,:)
      complex_dp  INTENTOUT   :: matout(size(mat2,1),size(mat2,2))
      complex_dp, allocatable :: vec(:),vec2(:)
      integer                 :: nn,n,i,n2
      n  = size(mat2,1) ; nn = n*(n+1)/2 ; allocate ( vec(n),vec2(n) )
      n2 = size(mat2,2)
      call check_dim(shape(mat1),[nn],'mat1',__LINE__)
      do i = 1,n2
        vec2 = mat2(:,i)
        call zhpmv('U',n,(1d0,0d0),mat1,vec2,1,(0d0,0d0),vec,1)
        matout(:,i) = vec
      enddo
      deallocate ( vec,vec2 )
      end SUBROUT_FUNC matmat_zpzm

      SUBROUT_FUNC matmat_zmzp( MATOUT mat1,mat2) RESULTMAT
      implicit none
      complex_dp, intent(in)  :: mat1(:,:),mat2(:)
      complex_dp  INTENTOUT   :: matout(size(mat1,1),size(mat1,2))
      complex_dp, allocatable :: vec(:),vec2(:)
      integer                 :: nn,n,i,n1
      n  = size(mat1,2) ; nn = n*(n+1)/2 ; allocate ( vec(n),vec2(n) )
      n1 = size(mat1,1)
      call check_dim(shape(mat2),[nn],'mat2',__LINE__)
      do i = 1,n1
        vec2 = conjg(mat1(i,:))
        call zhpmv('U',n,(1d0,0d0),mat2,vec2,1,(0d0,0d0),vec,1)
        matout(i,:) = conjg(vec)
      enddo
      deallocate ( vec,vec2 )
      end SUBROUT_FUNC matmat_zmzp

      SUBROUT_FUNC matmat_zmzm( MATOUT mat1,mat2) RESULTMAT
      implicit none
      complex_dp, intent(in) :: mat1(:,:),mat2(:,:)
      complex_dp  INTENTOUT  :: matout(size(mat1,1),size(mat2,2))
      integer                :: n1,n,n2
      n1 = size(mat1,1)
      n  = size(mat1,2)
      n2 = size(mat2,2)
      call check_dim(shape(mat2),[n,n2],'mat2',__LINE__)
      call zgemm('N','N',n1,n2,n,(1d0,0d0),mat1,n1,mat2,n,(0d0,0d0),matout,n1)
      end SUBROUT_FUNC matmat_zmzm

c     --------

      SUBROUT_FUNC macmat_dmdm( MATOUT mat1,mat2) RESULTMAT
      implicit none
      real_dp, intent(in) :: mat1(:,:),mat2(:,:)
      real_dp  INTENTOUT  :: matout(size(mat1,2),size(mat2,2))
      integer             :: n1,n,n2
      n1 = size(mat1,2)
      n  = size(mat1,1)
      n2 = size(mat2,2)
      call check_dim(shape(mat2),[n,n2],'mat2',__LINE__)
      call dgemm('C','N',n1,n2,n,1d0,mat1,n,mat2,n,0d0,matout,n1)
      end SUBROUT_FUNC macmat_dmdm

      SUBROUT_FUNC macmat_dmzm( MATOUT mat1,mat2) RESULTMAT
      implicit none
      real_dp,    intent(in) :: mat1(:,:)
      complex_dp, intent(in) :: mat2(:,:)
      complex_dp  INTENTOUT  :: matout(size(mat1,2),size(mat2,2))
      real_dp                :: matr(size(mat1,2),size(mat2,2)),mati(size(mat1,2),size(mat2,2))
      integer                :: n,n1,n2
      n1 = size(mat1,2)
      n  = size(mat1,1)
      n2 = size(mat2,2)
      call check_dim(shape(mat2),[n,n2],'mat2',__LINE__)
      call dgemm('C','N',n1,n2,n,1d0,mat1,n,real(mat2),n,0d0,matr,n1)
      call dgemm('C','N',n1,n2,n,1d0,mat1,n,imag(mat2),n,0d0,mati,n1)
      matout = matr + (0d0,1d0) * mati
      end SUBROUT_FUNC macmat_dmzm

      SUBROUT_FUNC macmat_zmdm( MATOUT mat1,mat2) RESULTMAT
      implicit none
      complex_dp, intent(in) :: mat1(:,:)
      real_dp,    intent(in) :: mat2(:,:)
      complex_dp  INTENTOUT  :: matout(size(mat1,2),size(mat2,2))
      real_dp                :: matr(size(mat1,2),size(mat2,2)),mati(size(mat1,2),size(mat2,2))
      integer                :: n,n1,n2
      n1 = size(mat1,2)
      n  = size(mat1,1)
      n2 = size(mat2,2)
      call check_dim(shape(mat2),[n,n2],'mat2',__LINE__)
      call dgemm('C','N',n1,n2,n,1d0,real(mat1),n,mat2,n,0d0,matr,n1)
      call dgemm('C','N',n1,n2,n,1d0,imag(mat1),n,mat2,n,0d0,mati,n1)
      matout = matr - (0d0,1d0) * mati
      end SUBROUT_FUNC macmat_zmdm

      SUBROUT_FUNC macmat_zmzm( MATOUT mat1,mat2) RESULTMAT
      implicit none
      complex_dp, intent(in) :: mat1(:,:),mat2(:,:)
      complex_dp  INTENTOUT  :: matout(size(mat1,2),size(mat2,2))
      integer                :: n,n1,n2
      n1 = size(mat1,2)
      n  = size(mat1,1)
      n2 = size(mat2,2)
      call check_dim(shape(mat2),[n,n2],'mat2',__LINE__)
      call zgemm('C','N',n1,n2,n,(1d0,0d0),mat1,n,mat2,n,(0d0,0d0),matout,n1)
      end SUBROUT_FUNC macmat_zmzm

c     --------

      SUBROUT_FUNC matmac_dmdm( MATOUT mat1,mat2) RESULTMAT
      implicit none
      real_dp, intent(in) :: mat1(:,:),mat2(:,:)
      real_dp  INTENTOUT  :: matout(size(mat1,1),size(mat2,1))
      integer             :: n1,n,n2
      n1 = size(mat1,1)
      n  = size(mat1,2)
      n2 = size(mat2,1)
      call check_dim(shape(mat2),[n2,n],'mat2',__LINE__)
      call dgemm('N','C',n1,n2,n,1d0,mat1,n1,mat2,n2,0d0,matout,n1)
      end SUBROUT_FUNC matmac_dmdm

      SUBROUT_FUNC matmac_dmzm( MATOUT mat1,mat2) RESULTMAT
      implicit none
      real_dp,    intent(in) :: mat1(:,:)
      complex_dp, intent(in) :: mat2(:,:)
      complex_dp  INTENTOUT  :: matout(size(mat1,1),size(mat2,1))
      integer                :: n,n1,n2
      n1 = size(mat1,1)
      n  = size(mat1,2)
      n2 = size(mat2,1)
      call check_dim(shape(mat2),[n2,n],'mat2',__LINE__)
      call dgemm('N','C',n1,2*n2,n,1d0,mat1,n1,mat2,2*n2,0d0,matout,n1)
      do n = 1,n2
        call redistribute_cvec(matout(:,n),n1)
      enddo
      end SUBROUT_FUNC matmac_dmzm

      SUBROUT_FUNC matmac_zmdm( MATOUT mat1,mat2) RESULTMAT
      implicit none
      complex_dp, intent(in) :: mat1(:,:)
      real_dp,    intent(in) :: mat2(:,:)
      complex_dp  INTENTOUT  :: matout(size(mat1,1),size(mat2,1))
      integer                :: n,n1,n2
      n1 = size(mat1,1)
      n  = size(mat1,2)
      n2 = size(mat2,1)
      call check_dim(shape(mat2),[n,n2],'mat2',__LINE__)
      call dgemm('N','C',2*n1,n2,n,1d0,mat1,2*n1,mat2,n2,0d0,matout,2*n1)
      end SUBROUT_FUNC matmac_zmdm

      SUBROUT_FUNC matmac_zmzm( MATOUT mat1,mat2) RESULTMAT
      implicit none
      complex_dp, intent(in) :: mat1(:,:),mat2(:,:)
      complex_dp  INTENTOUT  :: matout(size(mat1,1),size(mat2,1))
      integer                :: n1,n,n2
      n1 = size(mat1,1)
      n  = size(mat1,2)
      n2 = size(mat2,1)
      call check_dim(shape(mat2),[n2,n],'mat2',__LINE__)
      call zgemm('N','C',n1,n2,n,(1d0,0d0),mat1,n1,mat2,n2,(0d0,0d0),matout,n1)
      end SUBROUT_FUNC matmac_zmzm

# undef MATOUT
# undef INTENTOUT
# undef RESULTMAT
# undef SUBROUT_FUNC

c ------------------------------
# elif M_PART == 3
c ------------------------------

c     --------

      function macvec_dd(mat,vec)
      implicit none
      real_dp, intent(in) :: vec(:)
      real_dp, intent(in) :: mat(:,:)
      real_dp             :: macvec_dd(size(mat,2))
      integer             :: n1,n2
      n1 = size(mat,1)
      n2 = size(mat,2)
      call check_dim(shape(vec),[n1],'vec',__LINE__)
      call dgemv('C',n1,n2,1d0,mat,n1,vec,1,0d0,macvec_dd,1)
      end function macvec_dd

      function macvec_zd(mat,vec)
      implicit none
      real_dp,    intent(in) :: vec(:)
      complex_dp, intent(in) :: mat(:,:)
      complex_dp             :: macvec_zd(size(mat,2))
      integer                :: n1,n2
      n1 = size(mat,1)
      n2 = size(mat,2)
      call check_dim(shape(vec),[n1],'vec',__LINE__)
      call dgemv('C',n1,n2,1d0,imag(mat),n1,vec,1,0d0,macvec_zd,2) ; macvec_zd = (0d0,-1d0)*real(macvec_zd)
      call dgemv('C',n1,n2,1d0,real(mat),n1,vec,1,1d0,macvec_zd,2)
      end function macvec_zd

      function macvec_zz(mat,vec)
      implicit none
      complex_dp, intent(in) :: vec(:)
      complex_dp, intent(in) :: mat(:,:)
      complex_dp             :: macvec_zz(size(mat,2))
      integer                :: n1,n2
      n1 = size(mat,1)
      n2 = size(mat,2)
      call check_dim(shape(vec),[n1],'vec',__LINE__)
      call zgemv('C',n1,n2,(1d0,0d0),mat,n1,vec,1,(0d0,0d0),macvec_zz,1)
      end function macvec_zz

c     --------

      function diagmat_dd(diag,mat)
      implicit none
      real_dp, intent(in) :: diag(:),mat(:,:)
      real_dp             :: diagmat_dd(size(mat,1),size(mat,2))
      integer             :: n,i
      n = size(diag)
      call check_dim(shape(mat(:,1)),[n],'mat',__LINE__)
      do i = 1,n
        diagmat_dd(i,:) = diag(i) * mat(i,:)
      enddo
      end function diagmat_dd

      function diagmat_dz(diag,mat)
      implicit none
      real_dp,    intent(in) :: diag(:)
      complex_dp, intent(in) :: mat(:,:)
      complex_dp             :: diagmat_dz(size(mat,1),size(mat,2))
      integer                :: n,i
      n = size(diag)
      call check_dim(shape(mat(:,1)),[n],'mat',__LINE__)
      do i = 1,n
        diagmat_dz(i,:) = diag(i) * mat(i,:)
      enddo
      end function diagmat_dz

      function diagmat_zd(diag,mat)
      implicit none
      real_dp,    intent(in) :: mat(:,:)
      complex_dp, intent(in) :: diag(:)      
      complex_dp             :: diagmat_zd(size(mat,1),size(mat,2))
      integer                :: n,i
      n = size(diag)
      call check_dim(shape(mat(:,1)),[n],'mat',__LINE__)
      do i = 1,n
        diagmat_zd(i,:) = diag(i) * mat(i,:)
      enddo
      end function diagmat_zd

      function diagmat_zz(diag,mat)
      implicit none
      complex_dp, intent(in) :: diag(:),mat(:,:)
      complex_dp             :: diagmat_zz(size(mat,1),size(mat,2))
      integer                :: n,i
      n = size(diag)
      call check_dim(shape(mat(:,1)),[n],'mat',__LINE__)
      do i = 1,n
        diagmat_zz(i,:) = diag(i) * mat(i,:)
      enddo
      end function diagmat_zz

      ! had to rename diag -> dia due to a strange Intel compiler bug
      function matdiag_dd(mat,dia)
      implicit none
      real_dp, intent(in) :: dia(:),mat(:,:)
      real_dp             :: matdiag_dd(size(mat,1),size(mat,2))
      integer             :: n,i
      n = size(dia)
      call check_dim(shape(mat(1,:)),[n],'mat',__LINE__)
      do i = 1,n
        matdiag_dd(:,i) = dia(i) * mat(:,i)
      enddo
      end function matdiag_dd

      function matdiag_dz(mat,dia)
      implicit none
      real_dp,    intent(in) :: dia(:)
      complex_dp, intent(in) :: mat(:,:)
      complex_dp             :: matdiag_dz(size(mat,1),size(mat,2))
      integer                :: n,i
      n = size(dia)
      call check_dim(shape(mat(1,:)),[n],'mat',__LINE__)
      do i = 1,n
        matdiag_dz(:,i) = dia(i) * mat(:,i)
      enddo
      end function matdiag_dz

      function matdiag_zd(mat,dia)
      implicit none
      real_dp,    intent(in) :: mat(:,:)
      complex_dp, intent(in) :: dia(:)
      complex_dp             :: matdiag_zd(size(mat,1),size(mat,2))
      integer                :: n,i
      n = size(dia)
      call check_dim(shape(mat(1,:)),[n],'mat',__LINE__)
      do i = 1,n
        matdiag_zd(:,i) = dia(i) * mat(:,i)
      enddo
      end function matdiag_zd

      function matdiag_zz(mat,dia)
      implicit none
      complex_dp, intent(in) :: dia(:),mat(:,:)
      complex_dp             :: matdiag_zz(size(mat,1),size(mat,2))
      integer                :: n,i
      n = size(dia)
      call check_dim(shape(mat(1,:)),[n],'mat',__LINE__)
      do i = 1,n
        matdiag_zz(:,i) = dia(i) * mat(:,i)
      enddo
      end function matdiag_zz

c     --------

      function diagmatdiag_ddd(diag1,mat,diag2)
      implicit none
      real_dp, intent(in) :: diag1(:),mat(:,:),diag2(:)
      real_dp             :: diagmatdiag_ddd(size(diag1),size(diag2))
      integer             :: n1,n2,i
      n1 = size(diag1)
      n2 = size(diag2)
      call check_dim(shape(mat),[n1,n2],'mat',__LINE__)
      do i = 1,n2
        diagmatdiag_ddd(:,i) = diag2(i) * mat(:,i)
      enddo
      do i = 1,n1
        diagmatdiag_ddd(i,:) = diag1(i) * diagmatdiag_ddd(i,:)
      enddo
      end function diagmatdiag_ddd

      function diagmatdiag_dzd(diag1,mat,diag2)
      implicit none
      real_dp,    intent(in) :: diag1(:),diag2(:)
      complex_dp, intent(in) :: mat(:,:)
      complex_dp             :: diagmatdiag_dzd(size(diag1),size(diag2))
      integer                :: n1,n2,i
      n1 = size(diag1)
      n2 = size(diag2)      
      call check_dim(shape(mat),[n1,n2],'mat',__LINE__)
      do i = 1,n2
        diagmatdiag_dzd(:,i) = diag2(i) * mat(:,i)
      enddo
      do i = 1,n1
        diagmatdiag_dzd(i,:) = diag1(i) * diagmatdiag_dzd(i,:)
      enddo
      end function diagmatdiag_dzd

      function diagmatdiag_ddpd(diag,mat)
      implicit none
      real_dp, intent(in) :: diag(:),mat(:)
      real_dp             :: diagmatdiag_ddpd(size(mat))
      integer             :: n,i,j,k
      n = size(diag)
      call check_dim(shape(mat),[n*(n+1)/2],'mat',__LINE__)
      k = 0
      do j = 1,n
        do i = 1,j
          k                   = k + 1
          diagmatdiag_ddpd(k) = diag(i) * diag(j) * mat(k)
        enddo
      enddo
      end function diagmatdiag_ddpd

      function diagmatdiag_dzpd(diag,mat)
      implicit none
      real_dp,    intent(in) :: diag(:)
      complex_dp, intent(in) :: mat(:)
      complex_dp             :: diagmatdiag_dzpd(size(mat))
      integer                :: n,i,j,k
      n = size(diag)
      call check_dim(shape(mat),[n*(n+1)/2],'mat',__LINE__)
      k = 0
      do j = 1,n
        do i = 1,j
          k                   = k + 1
          diagmatdiag_dzpd(k) = diag(i) * diag(j) * mat(k)
        enddo
      enddo
      end function diagmatdiag_dzpd

c     --------

      subroutine p_diagmatdiag_ddd(out,diag1,mat,diag2)
      implicit none
      real_dp, intent(in)  :: diag1(:),mat(:,:),diag2(:)      
      real_dp, intent(out) :: out(:,:)
      integer              :: n1,n2,i
      n1 = size(diag1)
      n2 = size(diag2)      
      call check_dim(shape(mat),[n1,n2],'mat',__LINE__)
      call check_dim(shape(mat),[n1,n2],'out',__LINE__)
      do i = 1,n2
        out(:,i) = diag2(i) * mat(:,i)
      enddo
      do i = 1,n1
        out(i,:) = diag1(i) * out(i,:)
      enddo
      end subroutine p_diagmatdiag_ddd

      subroutine p_diagmatdiag_dzd(out,diag1,mat,diag2)
      implicit none
      complex_dp, intent(out) :: out(:,:)      
      complex_dp, intent(in)  :: mat(:,:)
      real_dp,    intent(in)  :: diag1(:),diag2(:)      
      integer                 :: n1,n2,i
      n1 = size(diag1)
      n2 = size(diag2)      
      call check_dim(shape(mat),[n1,n2],'mat',__LINE__)
      call check_dim(shape(mat),[n1,n2],'out',__LINE__)
      do i = 1,n2
        out(:,i) = diag2(i) * mat(:,i)
      enddo
      do i = 1,n1
        out(i,:) = diag1(i) * out(i,:)
      enddo
      end subroutine p_diagmatdiag_dzd

      subroutine p_diagmatdiag_ddpd(out,diag,mat)
      implicit none
      real_dp, intent(in)  :: diag(:),mat(:)
      real_dp, intent(out) :: out(:)
      integer              :: n,i,j,k
      n = size(diag)
      call check_dim(shape(mat),[n*(n+1)/2],'mat',__LINE__)
      call check_dim(shape(mat),[n*(n+1)/2],'out',__LINE__)      
      k = 0
      do j = 1,n
        do i = 1,j
          k      = k + 1
          out(k) = diag(i) * diag(j) * mat(k)
        enddo
      enddo
      end subroutine p_diagmatdiag_ddpd

      subroutine p_diagmatdiag_dzpd(out,diag,mat)
      implicit none
      complex_dp, intent(out) :: out(:)
      complex_dp, intent(in)  :: mat(:)
      real_dp,    intent(in)  :: diag(:)
      integer                 :: n,i,j,k
      n = size(diag)
      call check_dim(shape(mat),[n*(n+1)/2],'mat',__LINE__)
      call check_dim(shape(mat),[n*(n+1)/2],'out',__LINE__)      
      k = 0
      do j = 1,n
        do i = 1,j
          k      = k + 1
          out(k) = diag(i) * diag(j) * mat(k)
        enddo
      enddo
      end subroutine p_diagmatdiag_dzpd

      subroutine p_diagmatdiag_ddd1(diag1,mat,diag2)
      implicit none
      real_dp, intent(in)    :: diag1(:),diag2(:)      
      real_dp, intent(inout) :: mat(:,:)
      integer                :: n1,n2,i
      n1 = size(diag1)
      n2 = size(diag2)      
      call check_dim(shape(mat),[n1,n2],'mat',__LINE__)
      do i = 1,n2
        mat(:,i) = diag2(i) * mat(:,i)
      enddo
      do i = 1,n1
        mat(i,:) = diag1(i) * mat(i,:)
      enddo
      end subroutine p_diagmatdiag_ddd1

      subroutine p_diagmatdiag_dzd1(diag1,mat,diag2)
      implicit none
      complex_dp, intent(inout) :: mat(:,:)
      real_dp,    intent(in)    :: diag1(:),diag2(:)      
      integer                   :: n1,n2,i
      n1 = size(diag1)
      n2 = size(diag2)
      call check_dim(shape(mat),[n1,n2],'mat',__LINE__)
      do i = 1,n2
        mat(:,i) = diag2(i) * mat(:,i)
      enddo
      do i = 1,n1
        mat(i,:) = diag1(i) * mat(i,:)
      enddo
      end subroutine p_diagmatdiag_dzd1

      subroutine p_diagmatdiag_ddpd1(diag,mat)
      implicit none
      real_dp, intent(in)    :: diag(:)
      real_dp, intent(inout) :: mat(:)
      integer                :: n,i,j,k
      n = size(diag)
      call check_dim(shape(mat),[n*(n+1)/2],'mat',__LINE__)
      k = 0
      do j = 1,n
        do i = 1,j
          k      = k + 1
          mat(k) = diag(i) * diag(j) * mat(k)
        enddo
      enddo
      end subroutine p_diagmatdiag_ddpd1

      subroutine p_diagmatdiag_dzpd1(diag,mat)
      implicit none
      complex_dp, intent(inout) :: mat(:)
      real_dp,    intent(in)    :: diag(:)
      integer                   :: n,i,j,k
      n = size(diag)
      call check_dim(shape(mat),[n*(n+1)/2],'mat',__LINE__)
      k = 0
      do j = 1,n
        do i = 1,j
          k      = k + 1
          mat(k) = diag(i) * diag(j) * mat(k)
        enddo
      enddo
      end subroutine p_diagmatdiag_dzpd1

c     --------

      subroutine p_plusdiag1_d(d,mat)
      implicit none
      real_dp, intent(in)    :: d
      real_dp, intent(inout) :: mat(:,:)
      integer                :: n,i
      n = min(size(mat,1),size(mat,2))
      do i = 1,n
        mat(i,i) = d + mat(i,i)
      enddo      
      end subroutine p_plusdiag1_d
 
      subroutine p_plusdiag1_z(z,mat)
      implicit none
      real_dp,    intent(in)    :: z
      complex_dp, intent(inout) :: mat(:,:)
      integer                   :: n,i
      n = min(size(mat,1),size(mat,2))
      do i = 1,n
        mat(i,i) = z + mat(i,i)
      enddo      
      end subroutine p_plusdiag1_z

      subroutine p_plusdiag1_dp(d,mat)
      implicit none
      real_dp, intent(in)    :: d
      real_dp, intent(inout) :: mat(:)
      integer                :: n,i,j
      n = nint(sqrt(0.25d0+2*size(mat))-0.5d0) ; call check_dim(shape(mat),[n*(n+1)/2],'mat',__LINE__)
      j = 0
      do i = 1,n
        j      = j + i
        mat(j) = d + mat(j)
      enddo
      end subroutine p_plusdiag1_dp

      subroutine p_plusdiag1_zp(z,mat)
      implicit none
      real_dp,    intent(in)    :: z
      complex_dp, intent(inout) :: mat(:)
      integer                   :: n,i,j
      n = nint(sqrt(0.25d0+2*size(mat))-0.5d0) ; call check_dim(shape(mat),[n*(n+1)/2],'mat',__LINE__)
      j = 0
      do i = 1,n
        j      = j + i
        mat(j) = z + mat(j)
      enddo
      end subroutine p_plusdiag1_zp

      subroutine p_plusdiag_d(diag,mat)
      implicit none
      real_dp, intent(in)    :: diag(:)
      real_dp, intent(inout) :: mat(:,:)
      integer                :: n,i
      n = size(diag)
      call check_dim(shape(mat),[n,n],'mat',__LINE__)      
      do i = 1,n
        mat(i,i) = diag(i) + mat(i,i)
      enddo      
      end subroutine p_plusdiag_d

      subroutine p_plusdiag_z(diag,mat)
      implicit none
      real_dp,    intent(in)    :: diag(:)
      complex_dp, intent(inout) :: mat(:,:)
      integer                   :: n,i
      n = size(diag)
      call check_dim(shape(mat),[n,n],'mat',__LINE__)      
      do i = 1,n
        mat(i,i) = diag(i) + mat(i,i)
      enddo      
      end subroutine p_plusdiag_z

      subroutine p_plusdiag_dp(diag,mat)
      implicit none
      real_dp, intent(in)    :: diag(:)
      real_dp, intent(inout) :: mat(:)
      integer                :: n,i,j
      n = size(diag) ; call check_dim(shape(mat),[n*(n+1)/2],'mat',__LINE__)
      j = 0
      do i = 1,n
        j      = j + i
        mat(j) = diag(i) + mat(j)
      enddo
      end subroutine p_plusdiag_dp

      subroutine p_plusdiag_zp(diag,mat)
      implicit none
      real_dp,    intent(in)    :: diag(:)
      complex_dp, intent(inout) :: mat(:)
      integer                   :: n,i,j
      n = size(diag) ; call check_dim(shape(mat),[n*(n+1)/2],'mat',__LINE__)
      j = 0
      do i = 1,n
        j      = j + i
        mat(j) = diag(i) + mat(j)
      enddo
      end subroutine p_plusdiag_zp

c     --------

# include "w_unitarytrafo.inc"

c     --------

      function vecmatvec(vec,mat1,mat2) ! by 40% faster than dotprod(vec,matvec(mat,vec)) ( mat1=real(mat), mat2=imag(mat) )
      implicit none
      real_dp                :: vecmatvec
      real_dp,    intent(in) :: mat1(:),mat2(:)
      complex_dp, intent(in) :: vec(:)
      real_dp                :: vec1(size(vec)),vec2(size(vec)),r1,r2
      integer                :: n,i,j,k
      n         = size(vec)
      k         = 0
      vec1      = real(vec)
      vec2      = imag(vec)
      vecmatvec = 0
      do j = 1,n
        r1 = 0
        r2 = 0
        do i = 1,j-1
          k  = k + 1
          r1 = r1 + vec1(i) * mat1(k) + vec2(i) * mat2(k)
          r2 = r2 + vec2(i) * mat1(k) - vec1(i) * mat2(k)
        enddo
        k  = k + 1
        r1 = 2*r1 + vec1(j) * mat1(k)
        r2 = 2*r2 + vec2(j) * mat1(k)
        vecmatvec = vecmatvec + vec1(j) * r1 + vec2(j) * r2
      enddo
      end function vecmatvec

c     --------

      subroutine diagonalize_de(eval,mat)
      implicit none
      real_dp, intent(out) :: eval(:)
      real_dp, intent(in)  :: mat(:,:)
      real_dp, allocatable :: mat1(:,:),work(:)
      integer              :: n,info
      n = size(eval)
      if(n==0) Error('zero dimension in eigenvalue problem.')
      call check_dim(shape(mat),[n,n],'mat',__LINE__)
      allocate ( mat1(n,n),work(3*n) ) ; mat1 = mat
      call dsyev('N','U',n,mat1,n,eval,work,3*n,info) ; if(info/=0) Error('dsyev failed.')
      deallocate ( mat1,work )
      end subroutine diagonalize_de

      subroutine diagonalize_dv(evec,eval,mat)
      implicit none
      real_dp, intent(out) :: eval(:),evec(:,:)
      real_dp, intent(in)  :: mat(:,:)
      real_dp, allocatable :: work(:)
      integer              :: n,info
      n = size(eval)
      if(n==0) Error('zero dimension in eigenvalue problem.')
      call check_dim(shape(mat), [n,n],'mat',__LINE__)
      call check_dim(shape(evec),[n,n],'evec',__LINE__)
      allocate ( work(3*n) ) ; evec = mat
      call dsyev('V','U',n,evec,n,eval,work,3*n,info) ; if(info/=0) Error('dsyev failed.')
      deallocate ( work )
      end subroutine diagonalize_dv

      subroutine diagonalize_dpe(eval,mat)
      implicit none
      real_dp, intent(out) :: eval(:)
      real_dp, intent(in)  :: mat(:)
      real_dp, allocatable :: mat1(:),work(:)
      integer              :: n,nn,info
      n = size(eval) ; nn = n*(n+1)/2
      if(n==0) Error('zero dimension in eigenvalue problem.')
      call check_dim(shape(mat),[nn],'mat',__LINE__)
      allocate ( mat1(nn),work(3*n) ) ; mat1 = mat
      call dspev('N','U',n,mat1,eval,work,n,work,info) ; if(info/=0) Error('dspev failed.')
      deallocate ( mat1,work )
      end subroutine diagonalize_dpe

      subroutine diagonalize_dpv(evec,eval,mat)
      implicit none
      real_dp, intent(out) :: eval(:),evec(:,:)
      real_dp, intent(in)  :: mat(:)
      real_dp, allocatable :: mat1(:),work(:)
      integer              :: n,nn,info
      n = size(eval) ; nn = n*(n+1)/2
      if(n==0) Error('zero dimension in eigenvalue problem.')
      call check_dim(shape(mat),[nn],'mat',__LINE__)
      call check_dim(shape(evec),[n,n],'evec',__LINE__)
      allocate ( mat1(nn),work(3*n) ) ; mat1 = mat
      call dspev('V','U',n,mat1,eval,evec,n,work,info) ; if(info/=0) Error('dspev failed.')
      deallocate ( mat1,work )
      end subroutine diagonalize_dpv

      subroutine diagonalize_ze(eval,mat)
      implicit none
      real_dp,    intent(out) :: eval(:)
      complex_dp, intent(in)  :: mat(:,:)
      complex_dp, allocatable :: mat1(:,:),work(:)
      real_dp,    allocatable :: rwork(:)
      integer                 :: n,info
      n = size(eval)
      if(n==0) Error('zero dimension in eigenvalue problem.')
      call check_dim(shape(mat),[n,n],'mat',__LINE__)
      allocate ( mat1(n,n),work(3*n),rwork(3*n) ) ; mat1 = mat
      call zheev('N','U',n,mat1,n,eval,work,3*n,rwork,info) ; if(info/=0) Error('zheev failed.')
      deallocate ( mat1,work,rwork )
      end subroutine diagonalize_ze

      subroutine diagonalize_zv(evec,eval,mat)
      implicit none
      real_dp,    intent(out) :: eval(:)
      complex_dp, intent(out) :: evec(:,:)
      complex_dp, intent(in)  :: mat(:,:)
      complex_dp, allocatable :: work(:)
      real_dp,    allocatable :: rwork(:)
      integer                 :: n,info
      n = size(eval)
      if(n==0) Error('zero dimension in eigenvalue problem.')
      call check_dim(shape(mat), [n,n],'mat',__LINE__)
      call check_dim(shape(evec),[n,n],'evec',__LINE__)
      allocate ( work(3*n),rwork(3*n) ) ; evec = mat
      call zheev('V','U',n,evec,n,eval,work,3*n,rwork,info) ; if(info/=0) Error('zheev failed.')
      deallocate ( work,rwork )
      end subroutine diagonalize_zv

      subroutine diagonalize_zpe(eval,mat)
      implicit none
      real_dp,    intent(out) :: eval(:)
      complex_dp, intent(in)  :: mat(:)
      complex_dp, allocatable :: mat1(:),work(:)
      real_dp,    allocatable :: rwork(:)
      integer                 :: n,nn,info
      n = size(eval) ; nn = n*(n+1)/2
      if(n==0) Error('zero dimension in eigenvalue problem.')
      call check_dim(shape(mat),[nn],'mat',__LINE__)
      allocate ( mat1(nn),work(3*n),rwork(3*n) ) ; mat1 = mat
      call zhpev('N','U',n,mat1,eval,work,n,work,rwork,info) ; if(info/=0) Error('zhpev failed.')
      deallocate ( mat1,work,rwork )
      end subroutine diagonalize_zpe

      subroutine diagonalize_zpv(evec,eval,mat)
      implicit none
      real_dp,    intent(out) :: eval(:)
      complex_dp, intent(out) :: evec(:,:)
      complex_dp, intent(in)  :: mat(:)
      complex_dp, allocatable :: mat1(:),work(:)
      real_dp,    allocatable :: rwork(:)
      integer                 :: n,nn,info
      n = size(eval) ; nn = n*(n+1)/2
      if(n==0) Error('zero dimension in eigenvalue problem.')
      call check_dim(shape(mat),[nn],'mat',__LINE__)
      call check_dim(shape(evec),[n,n],'evec',__LINE__)
      allocate ( mat1(nn),work(3*n),rwork(3*n) ) ; mat1 = mat
      call zhpev('V','U',n,mat1,eval,evec,n,work,rwork,info) ; if(info/=0) Error('zhpev failed.')
      deallocate ( mat1,work,rwork )
      end subroutine diagonalize_zpv

      subroutine diagonalize_deo(eval,mat,olap)
      implicit none
      real_dp, intent(out) :: eval(:)
      real_dp, intent(in)  :: mat(:,:),olap(:,:)
      real_dp, allocatable :: mat1(:,:),olap1(:,:),work(:)
      integer              :: n,info
      n = size(eval)
      if(n==0) Error('zero dimension in eigenvalue problem.')
      call check_dim(shape(mat), [n,n],'mat',__LINE__)
      call check_dim(shape(olap),[n,n],'olap',__LINE__)
      allocate ( mat1(n,n),olap1(n,n),work(3*n) ) ; mat1 = mat ; olap1 = olap
      call dsygv(1,'N','U',n,mat1,n,olap1,n,eval,work,3*n,info) ; if(info/=0) Error('dsygv failed.')
      deallocate ( mat1,olap1,work )
      end subroutine diagonalize_deo

      subroutine diagonalize_dvo(evec,eval,mat,olap)
      implicit none
      real_dp, intent(out) :: eval(:),evec(:,:)
      real_dp, intent(in)  :: mat(:,:),olap(:,:)
      real_dp, allocatable :: olap1(:,:),work(:)
      integer              :: n,info
      n = size(eval)
      if(n==0) Error('zero dimension in eigenvalue problem.')
      call check_dim(shape(mat), [n,n],'mat',__LINE__)
      call check_dim(shape(evec),[n,n],'evec',__LINE__)
      call check_dim(shape(olap),[n,n],'olap',__LINE__)      
      allocate ( olap1(n,n),work(3*n) ) ; evec = mat ; olap1 = olap
      call dsygv(1,'V','U',n,evec,n,olap1,n,eval,work,3*n,info) ; if(info/=0) Error('dsygv failed.')
      deallocate ( olap1,work )
      end subroutine diagonalize_dvo

      subroutine diagonalize_dpeo(eval,mat,olap)
      implicit none
      real_dp, intent(out) :: eval(:)
      real_dp, intent(in)  :: mat(:),olap(:)
      real_dp, allocatable :: mat1(:),olap1(:),work(:)
      integer              :: n,nn,info
      n = size(eval) ; nn = n*(n+1)/2
      if(n==0) Error('zero dimension in eigenvalue problem.')
      call check_dim(shape(mat), [nn],'mat',__LINE__)
      call check_dim(shape(olap),[nn],'olap',__LINE__)
      allocate ( mat1(nn),olap1(nn),work(3*n) ) ; mat1 = mat ; olap1 = olap
      call dspgv(1,'N','U',n,mat1,olap1,eval,work,n,work,info) ; if(info/=0) Error('dspgv failed.')
      deallocate ( mat1,olap1,work )
      end subroutine diagonalize_dpeo

      subroutine diagonalize_dpvo(evec,eval,mat,olap)
      implicit none
      real_dp, intent(out) :: eval(:),evec(:,:)
      real_dp, intent(in)  :: mat(:),olap(:)
      real_dp, allocatable :: mat1(:),olap1(:),work(:)
      integer              :: n,nn,info
      n = size(eval) ; nn = n*(n+1)/2
      if(n==0) Error('zero dimension in eigenvalue problem.')
      call check_dim(shape(mat), [nn],'mat',__LINE__)
      call check_dim(shape(olap),[nn],'olap',__LINE__)
      call check_dim(shape(evec),[n,n],'evec',__LINE__)
      allocate ( mat1(nn),olap1(nn),work(3*n) ) ; mat1 = mat ; olap1 = olap
      call dspgv(1,'V','U',n,mat1,olap1,eval,evec,n,work,info) ; if(info/=0) Error('dspgv failed.')
      deallocate ( mat1,olap1,work )
      end subroutine diagonalize_dpvo

      subroutine diagonalize_zeo(eval,mat,olap)
      implicit none
      real_dp,    intent(out) :: eval(:)
      complex_dp, intent(in)  :: mat(:,:),olap(:,:)
      complex_dp, allocatable :: mat1(:,:),olap1(:,:),work(:)
      real_dp,    allocatable :: rwork(:)
      integer                 :: n,info
      n = size(eval)
      if(n==0) Error('zero dimension in eigenvalue problem.')
      call check_dim(shape(mat), [n,n],'mat',__LINE__)
      call check_dim(shape(olap),[n,n],'olap',__LINE__)
      allocate ( mat1(n,n),olap1(n,n),work(3*n),rwork(3*n) ) ; mat1 = mat ; olap1 = olap
      call zhegv(1,'N','U',n,mat1,n,olap1,n,eval,work,3*n,rwork,info) ; if(info/=0) Error('zhegv failed.')
      deallocate ( mat1,olap1,work,rwork )
      end subroutine diagonalize_zeo

      subroutine diagonalize_zvo(evec,eval,mat,olap)
      implicit none
      real_dp,    intent(out) :: eval(:)
      complex_dp, intent(out) :: evec(:,:)
      complex_dp, intent(in)  :: mat(:,:),olap(:,:)
      complex_dp, allocatable :: olap1(:,:),work(:)
      real_dp,    allocatable :: rwork(:)
      integer                 :: n,info
      n = size(eval)
      if(n==0) Error('zero dimension in eigenvalue problem.')
      call check_dim(shape(mat), [n,n],'mat',__LINE__)
      call check_dim(shape(olap),[n,n],'olap',__LINE__)
      call check_dim(shape(evec),[n,n],'evec',__LINE__)
      allocate ( olap1(n,n),work(3*n),rwork(3*n) ) ; evec = mat ; olap1 = olap
      call zhegv(1,'V','U',n,evec,n,olap1,n,eval,work,3*n,rwork,info) ; if(info/=0) Error('zhegv failed.')
      deallocate ( olap1,work,rwork )
      end subroutine diagonalize_zvo

      subroutine diagonalize_zpeo(eval,mat,olap)
      implicit none
      real_dp,    intent(out) :: eval(:)
      complex_dp, intent(in)  :: mat(:),olap(:)
      complex_dp, allocatable :: mat1(:),olap1(:),work(:)
      real_dp,    allocatable :: rwork(:)
      integer                 :: n,nn,info
      n = size(eval) ; nn = n*(n+1)/2
      if(n==0) Error('zero dimension in eigenvalue problem.')
      call check_dim(shape(mat), [nn],'mat',__LINE__)
      call check_dim(shape(olap),[nn],'olap',__LINE__)
      allocate ( mat1(nn),olap1(nn),work(3*n),rwork(3*n) ) ; mat1 = mat ; olap1 = olap
      call zhpgv(1,'N','U',n,mat1,olap1,eval,work,n,work,rwork,info) ; if(info/=0) Error('zhpev failed.')
      deallocate ( mat1,olap1,work,rwork )
      end subroutine diagonalize_zpeo

      subroutine diagonalize_zpvo(evec,eval,mat,olap)
      implicit none
      real_dp,    intent(out) :: eval(:)
      complex_dp, intent(out) :: evec(:,:)
      complex_dp, intent(in)  :: mat(:),olap(:)
      complex_dp, allocatable :: mat1(:),olap1(:),work(:)
      real_dp,    allocatable :: rwork(:)
      integer                 :: n,nn,info
      n = size(eval) ; nn = n*(n+1)/2
      if(n==0) Error('zero dimension in eigenvalue problem.')
      call check_dim(shape(mat), [nn],'mat',__LINE__)
      call check_dim(shape(olap),[nn],'olap',__LINE__)
      call check_dim(shape(evec),[n,n],'evec',__LINE__)
      allocate ( mat1(nn),olap1(nn),work(3*n),rwork(3*n) ) ; mat1 = mat ; olap1 = olap
      call zhpgv(1,'V','U',n,mat1,olap1,eval,evec,n,work,rwork,info) ; if(info/=0) Error('zhpgv failed.')
      deallocate ( mat1,olap1,work,rwork )
      end subroutine diagonalize_zpvo

      subroutine diagonalize_dvs(evec,eval,mat,m)
      implicit none
      real_dp, intent(out) :: eval(:),evec(:,:)
      real_dp, intent(in)  :: mat(:,:)
      real_dp, allocatable :: work(:),mat1(:,:)
      real_dp              :: abstol,dlamch
      integer, intent(in)  :: m
      integer, allocatable :: iwork(:),ifail(:)
      integer              :: n,ma,idum,info
      ma = abs(m)
      n  = size(eval)
      if(n==0) Error('zero dimension in eigenvalue problem.')
      if(ma>n) Error('number of selected eigenvalues exceeds maximal number.')
      call check_dim(shape(mat),[n,n],'mat',__LINE__)
      call check_dim(shape(evec),[n,-ma],'evec',__LINE__)
      allocate ( work(8*n),iwork(5*n),mat1(n,n),ifail(n) ) ; mat1 = mat
      abstol = 2 * dlamch('S')
      if(m>0) then
        call dsyevx('V','I','U',n,mat1,n,0d0,0d0,n-ma+1,n,abstol,idum,eval,evec,n,work,8*n,iwork,ifail,info)
      else
        call dsyevx('V','I','U',n,mat1,n,0d0,0d0,    1,ma,abstol,idum,eval,evec,n,work,8*n,iwork,ifail,info)
      endif
      if(info/=0) Error('dsyevx failed.')
      deallocate ( work,iwork,mat1,ifail )
      end subroutine diagonalize_dvs

      subroutine diagonalize_dvos(evec,eval,mat,olap,m)
      implicit none
      real_dp, intent(out) :: eval(:),evec(:,:)
      real_dp, intent(in)  :: mat(:,:),olap(:,:)
      real_dp, allocatable :: work(:),mat1(:,:),olap1(:,:)
      real_dp              :: abstol,dlamch
      integer, intent(in)  :: m
      integer, allocatable :: iwork(:),ifail(:)
      integer              :: n,ma,idum,info
      ma = abs(m)
      n  = size(eval)
      if(n==0) Error('zero dimension in eigenvalue problem.')
      if(ma>n) Error('number of selected eigenvalues exceeds maximal number.')
      call check_dim(shape(mat), [n,n],'mat',__LINE__)
      call check_dim(shape(olap),[n,n],'olap',__LINE__)
      call check_dim(shape(evec),[n,-ma],'evec',__LINE__)      
      allocate ( work(8*n),iwork(5*n),mat1(n,n),olap1(n,n),ifail(n) ) ; mat1 = mat ; olap1 = olap
      abstol = 2 * dlamch('S')
      if(m>0) then
        call dsygvx(1,'V','I','U',n,mat1,n,olap1,n,0d0,0d0,n-ma+1,n,abstol,idum,eval,evec,n,work,8*n,iwork,ifail,info)
      else
        call dsygvx(1,'V','I','U',n,mat1,n,olap1,n,0d0,0d0,    1,ma,abstol,idum,eval,evec,n,work,8*n,iwork,ifail,info)
      endif
      if(info/=0) Error('dsygvx failed.')
      deallocate ( work,iwork,mat1,olap1,ifail )
      end subroutine diagonalize_dvos

      subroutine diagonalize_dpvs(evec,eval,mat,m)
      implicit none
      real_dp, intent(out) :: eval(:),evec(:,:)
      real_dp, intent(in)  :: mat(:)
      real_dp, allocatable :: work(:),mat1(:)
      real_dp              :: abstol,dlamch
      integer, intent(in)  :: m
      integer, allocatable :: iwork(:),ifail(:)
      integer              :: n,nn,ma,idum,info
      ma = abs(m)
      n  = size(eval)
      nn = n*(n+1)/2
      if(n==0) Error('zero dimension in eigenvalue problem.')
      if(ma>n) Error('number of selected eigenvalues exceeds maximal number.')
      call check_dim(shape(mat),[nn],'mat',__LINE__)
      call check_dim(shape(evec),[n,-ma],'evec',__LINE__)
      allocate ( work(8*n),iwork(5*n),mat1(nn),ifail(n) ) ; mat1 = mat
      abstol = 2 * dlamch('S')
      if(m>0) then
        call dspevx('V','I','U',n,mat1,0d0,0d0,n-ma+1,n,abstol,idum,eval,evec,n,work,iwork,ifail,info)
      else
        call dspevx('V','I','U',n,mat1,0d0,0d0,    1,ma,abstol,idum,eval,evec,n,work,iwork,ifail,info)
      endif
      if(info/=0) Error('dspevx failed.')
      deallocate ( work,iwork,mat1,ifail )
      end subroutine diagonalize_dpvs

      subroutine diagonalize_dpvos(evec,eval,mat,olap,m)
      implicit none
      real_dp, intent(out) :: eval(:),evec(:,:)
      real_dp, intent(in)  :: mat(:),olap(:)
      real_dp, allocatable :: work(:),mat1(:),olap1(:)
      real_dp              :: abstol,dlamch
      integer, intent(in)  :: m
      integer, allocatable :: iwork(:),ifail(:)
      integer              :: n,nn,ma,idum,info
      ma = abs(m)
      n  = size(eval)
      nn = n*(n+1)/2
      if(n==0) Error('zero dimension in eigenvalue problem.')
      if(ma>n) Error('number of selected eigenvalues exceeds maximal number.')
      call check_dim(shape(mat), [nn],'mat',__LINE__)
      call check_dim(shape(olap),[nn],'olap',__LINE__)
      call check_dim(shape(evec),[n,-ma],'evec',__LINE__)
      allocate ( work(8*n),iwork(5*n),mat1(nn),olap1(nn),ifail(n) ) ; mat1 = mat ; olap1 = olap
      abstol = 2 * dlamch('S')
      if(m>0) then
        call dspgvx(1,'V','I','U',n,mat1,olap1,0d0,0d0,n-ma+1,n,abstol,idum,eval,evec,n,work,iwork,ifail,info)
      else
        call dspgvx(1,'V','I','U',n,mat1,olap1,0d0,0d0,    1,ma,abstol,idum,eval,evec,n,work,iwork,ifail,info)
      endif
      if(info/=0) Error('dspgvx failed.')
      deallocate ( work,iwork,mat1,olap1,ifail )
      end subroutine diagonalize_dpvos

      subroutine diagonalize_zvs(evec,eval,mat,m)
      implicit none
      real_dp,    intent(out) :: eval(:)
      complex_dp, intent(out) :: evec(:,:)
      complex_dp, intent(in)  :: mat(:,:)
      complex_dp, allocatable :: work(:),mat1(:,:)
      real_dp,    allocatable :: rwork(:)
      real_dp                 :: abstol,dlamch
      integer,    intent(in)  :: m
      integer,    allocatable :: iwork(:),ifail(:)
      integer                 :: n,ma,idum,info
      ma = abs(m)
      n  = size(eval)
      if(n==0) Error('zero dimension in eigenvalue problem.')
      if(ma>n) Error('number of selected eigenvalues exceeds maximal number.')
      call check_dim(shape(mat),[n,n],'mat',__LINE__)
      call check_dim(shape(evec),[n,-ma],'evec',__LINE__)
      allocate ( work(2*n),rwork(7*n),iwork(5*n),mat1(n,n),ifail(n) ) ; mat1 = mat
      abstol = 2 * dlamch('S')
      ! ZHEEVX calls ZLANHE. In older Intel MKL libraries, this function has a potential bug which can lead to a segfault.
      if(m>0) then
        call zheevx('V','I','U',n,mat1,n,0d0,0d0,n-ma+1,n,abstol,idum,eval,evec,n,work,2*n,rwork,iwork,ifail,info)
      else
        call zheevx('V','I','U',n,mat1,n,0d0,0d0,    1,ma,abstol,idum,eval,evec,n,work,2*n,rwork,iwork,ifail,info)
      endif
      if(info/=0) Error('zheevx failed.')
      deallocate ( work,rwork,iwork,mat1,ifail )
      end subroutine diagonalize_zvs

      subroutine diagonalize_zvos(evec,eval,mat,olap,m)
      implicit none
      real_dp,    intent(out) :: eval(:)
      complex_dp, intent(out) :: evec(:,:)
      complex_dp, intent(in)  :: mat(:,:),olap(:,:)
      complex_dp, allocatable :: work(:),mat1(:,:),olap1(:,:)
      real_dp,    allocatable :: rwork(:)
      real_dp                 :: abstol,dlamch
      integer,    intent(in)  :: m
      integer,    allocatable :: iwork(:),ifail(:)
      integer                 :: n,ma,idum,info
      ma = abs(m)
      n  = size(eval)
      if(n==0) Error('zero dimension in eigenvalue problem.')
      if(ma>n) Error('number of selected eigenvalues exceeds maximal number.')
      call check_dim(shape(mat), [n,n],'mat',__LINE__)
      call check_dim(shape(olap),[n,n],'olap',__LINE__)
      call check_dim(shape(evec),[n,-ma],'evec',__LINE__)
      allocate ( work(2*n),rwork(7*n),iwork(5*n),mat1(n,n),olap1(n,n),ifail(n) ) ; mat1 = mat ; olap1 = olap
      abstol = 2 * dlamch('S')
      if(m>0) then
        call zhegvx(1,'V','I','U',n,mat1,n,olap1,n,0d0,0d0,n-ma+1,n,abstol,idum,eval,evec,n,work,2*n,rwork,iwork,ifail,info)
      else
        call zhegvx(1,'V','I','U',n,mat1,n,olap1,n,0d0,0d0,    1,ma,abstol,idum,eval,evec,n,work,2*n,rwork,iwork,ifail,info)
      endif
      if(info/=0) Error('zhegvx failed.')
      deallocate ( work,rwork,iwork,mat1,olap1,ifail )
      end subroutine diagonalize_zvos

      subroutine diagonalize_zpvs(evec,eval,mat,m)
      implicit none
      real_dp,    intent(out) :: eval(:)
      complex_dp, intent(out) :: evec(:,:)
      complex_dp, intent(in)  :: mat(:)
      complex_dp, allocatable :: work(:),mat1(:)
      real_dp,    allocatable :: rwork(:)
      real_dp                 :: abstol,dlamch
      integer,    intent(in)  :: m
      integer,    allocatable :: iwork(:),ifail(:)
      integer                 :: n,nn,ma,idum,info
      ma = abs(m)
      n  = size(eval)
      nn = n*(n+1)/2
      if(n==0) Error('zero dimension in eigenvalue problem.')
      if(ma>n) Error('number of selected eigenvalues exceeds maximal number.')
      call check_dim(shape(mat), [nn],'mat',__LINE__)
      call check_dim(shape(evec),[n,-ma],'evec',__LINE__)
      allocate ( work(2*n),rwork(7*n),iwork(5*n),mat1(nn),ifail(n) ) ; mat1 = mat
      abstol = 2 * dlamch('S')
      if(m>0) then
        call zhpevx('V','I','U',n,mat1,0d0,0d0,n-ma+1,n,abstol,idum,eval,evec,n,work,iwork,ifail,info)
      else
        call zhpevx('V','I','U',n,mat1,0d0,0d0,    1,ma,abstol,idum,eval,evec,n,work,iwork,ifail,info)
      endif
      if(info/=0) Error('zhpevx failed.')
      deallocate ( work,rwork,iwork,mat1,ifail )
      end subroutine diagonalize_zpvs

      subroutine diagonalize_zpvos(evec,eval,mat,olap,m)
      implicit none
      real_dp,    intent(out) :: eval(:)
      complex_dp, intent(out) :: evec(:,:)
      complex_dp, intent(in)  :: mat(:),olap(:)
      complex_dp, allocatable :: work(:),mat1(:),olap1(:)
      real_dp,    allocatable :: rwork(:)
      real_dp                 :: abstol,dlamch
      integer,    intent(in)  :: m
      integer,    allocatable :: iwork(:),ifail(:)
      integer                 :: n,nn,ma,idum,info
      ma = abs(m)
      n  = size(eval)
      nn = n*(n+1)/2
      if(n==0) Error('zero dimension in eigenvalue problem.')
      if(ma>n) Error('number of selected eigenvalues exceeds maximal number.')
      call check_dim(shape(mat), [nn],'mat',__LINE__)
      call check_dim(shape(olap),[nn],'olap',__LINE__)
      call check_dim(shape(evec),[n,-ma],'evec',__LINE__)
      allocate ( work(2*n),rwork(7*n),iwork(5*n),mat1(nn),olap1(nn),ifail(n) ) ; mat1 = mat ; olap1 = olap
      abstol = 2 * dlamch('S')
      if(m>0) then
        call zhpgvx(1,'V','I','U',n,mat1,olap1,0d0,0d0,n-ma+1,n,abstol,idum,eval,evec,n,work,rwork,iwork,ifail,info)
      else
        call zhpgvx(1,'V','I','U',n,mat1,olap1,0d0,0d0,    1,ma,abstol,idum,eval,evec,n,work,rwork,iwork,ifail,info)
      endif
      if(info/=0) Error('zhpgvx failed.')
      deallocate ( work,rwork,iwork,mat1,olap1,ifail )
      end subroutine diagonalize_zpvos

      ! routines for diagonalization: eigenvalue range [r1,r2) or index range [ir1,ir2].
      ! the number of actually found eigenvectors is returned in ir2.

      subroutine diagonalize_dvx(evec,eval,mat,ir1,ir2,r1,r2)
      implicit none
      real_dp, intent(out)   :: eval(:),evec(:,:)
      real_dp, intent(in)    :: mat(:,:)
      real_dp, allocatable   :: work(:),mat1(:,:)
      real_dp                :: abstol,dlamch
      real_dp, intent(in)    :: r1,r2
      integer, intent(in)    :: ir1
      integer, intent(inout) :: ir2
      integer, allocatable   :: iwork(:),ifail(:)
      integer                :: n,m,idum,info
      n = size(eval)
      m = ir2 - ir1 + 1
      if(n==0) Error('zero dimension in eigenvalue problem.')
      call check_dim(shape(mat),[n,n],'mat',__LINE__)
      if(r1>=r2) then
        if(m<0) Error('negative index range.')
        if(m>n) Error('number of selected eigenvalues exceeds maximal number.')
        call check_dim(shape(evec),[n,-m],'evec',__LINE__)
      else
        call check_dim(shape(evec),[n,n],'evec',__LINE__)
      endif
      allocate ( work(8*n),iwork(5*n),mat1(n,n),ifail(n) ) ; mat1 = mat
      abstol = 2 * dlamch('S')
      if(r1<r2) then
        call dsyevx('V','V','U',n,mat1,n,r1,r2,0,0,      abstol,idum,eval,evec,n,work,8*n,iwork,ifail,info)
      else
        call dsyevx('V','I','U',n,mat1,n,0d0,0d0,ir1,ir2,abstol,idum,eval,evec,n,work,8*n,iwork,ifail,info)
      endif
      ir2 = idum
      if(info/=0) Error('dsyevx failed.')
      deallocate ( work,iwork,mat1,ifail )
      end subroutine diagonalize_dvx

      subroutine diagonalize_dvox(evec,eval,mat,olap,ir1,ir2,r1,r2)
      implicit none
      real_dp, intent(out)   :: eval(:),evec(:,:)
      real_dp, intent(in)    :: mat(:,:),olap(:,:)
      real_dp, allocatable   :: work(:),mat1(:,:),olap1(:,:)
      real_dp                :: abstol,dlamch
      real_dp, intent(in)    :: r1,r2
      integer, intent(in)    :: ir1
      integer, intent(inout) :: ir2
      integer, allocatable   :: iwork(:),ifail(:)
      integer                :: n,m,idum,info
      n = size(eval)
      m = ir2 - ir1 + 1
      if(n==0) Error('zero dimension in eigenvalue problem.')
      call check_dim(shape(mat), [n,n],'mat',__LINE__)
      call check_dim(shape(olap),[n,n],'olap',__LINE__)
      if(r1>=r2) then
        if(m<0) Error('negative index range.')
        if(m>n) Error('number of selected eigenvalues exceeds maximal number.')
        call check_dim(shape(evec),[n,-m],'evec',__LINE__)
      else
        call check_dim(shape(evec),[n,n],'evec',__LINE__)
      endif
      allocate ( work(8*n),iwork(5*n),mat1(n,n),olap1(n,n),ifail(n) ) ; mat1 = mat ; olap1 = olap
      abstol = 2 * dlamch('S')
      if(r1<r2) then
        call dsygvx(1,'V','V','U',n,mat1,n,olap1,n,r1,r2,0,0,      abstol,idum,eval,evec,n,work,8*n,iwork,ifail,info)
      else
        call dsygvx(1,'V','I','U',n,mat1,n,olap1,n,0d0,0d0,ir1,ir2,abstol,idum,eval,evec,n,work,8*n,iwork,ifail,info)
      endif
      ir2 = idum
      if(info/=0) Error('dsygvx failed.')
      deallocate ( work,iwork,mat1,olap1,ifail )
      end subroutine diagonalize_dvox

      subroutine diagonalize_dpvx(evec,eval,mat,ir1,ir2,r1,r2)
      implicit none
      real_dp, intent(out)   :: eval(:),evec(:,:)
      real_dp, intent(in)    :: mat(:)
      real_dp, allocatable   :: work(:),mat1(:)
      real_dp                :: abstol,dlamch
      real_dp, intent(in)    :: r1,r2
      integer, intent(in)    :: ir1
      integer, intent(inout) :: ir2
      integer, allocatable   :: iwork(:),ifail(:)
      integer                :: n,m,nn,idum,info
      n  = size(eval)
      nn = n*(n+1)/2
      m  = ir2 - ir1 + 1
      if(n==0) Error('zero dimension in eigenvalue problem.')
      call check_dim(shape(mat),[nn],'mat',__LINE__)
      if(r1>=r2) then
        if(m<0) Error('negative index range.')
        if(m>n) Error('number of selected eigenvalues exceeds maximal number.')
        call check_dim(shape(evec),[n,-m],'evec',__LINE__)
      else        
        call check_dim(shape(evec),[n,n],'evec',__LINE__)
      endif
      allocate ( work(8*n),iwork(5*n),mat1(nn),ifail(n) ) ; mat1 = mat
      abstol = 2 * dlamch('S')
      if(r1<r2) then
        call dspevx('V','V','U',n,mat1,r1,r2,0,0,      abstol,idum,eval,evec,n,work,iwork,ifail,info)
      else
        call dspevx('V','I','U',n,mat1,0d0,0d0,ir1,ir2,abstol,idum,eval,evec,n,work,iwork,ifail,info)
      endif
      ir2 = idum
      if(info/=0) Error('dspevx failed.')
      deallocate ( work,iwork,mat1,ifail )
      end subroutine diagonalize_dpvx

      subroutine diagonalize_dpvox(evec,eval,mat,olap,ir1,ir2,r1,r2)
      implicit none
      real_dp, intent(out)   :: eval(:),evec(:,:)
      real_dp, intent(in)    :: mat(:),olap(:)
      real_dp, allocatable   :: work(:),mat1(:),olap1(:)
      real_dp                :: abstol,dlamch
      real_dp, intent(in)    :: r1,r2
      integer, intent(in)    :: ir1
      integer, intent(inout) :: ir2
      integer, allocatable   :: iwork(:),ifail(:)
      integer                :: n,nn,m,idum,info
      n  = size(eval)
      nn = n*(n+1)/2
      m  = ir2 - ir1 + 1
      if(n==0) Error('zero dimension in eigenvalue problem.')
      call check_dim(shape(mat), [nn],'mat',__LINE__)
      call check_dim(shape(olap),[nn],'olap',__LINE__)
      if(r1>=r2) then
        if(m<0) Error('negative index range.')
        if(m>n) Error('number of selected eigenvalues exceeds maximal number.')
        call check_dim(shape(evec),[n,-m],'evec',__LINE__)
      else
        call check_dim(shape(evec),[n,n],'evec',__LINE__)
      endif
      allocate ( work(8*n),iwork(5*n),mat1(nn),olap1(nn),ifail(n) ) ; mat1 = mat ; olap1 = olap
      abstol = 2 * dlamch('S')
      if(r1<r2) then
        call dspgvx(1,'V','V','U',n,mat1,olap1,r1,r2,0,0,      abstol,idum,eval,evec,n,work,iwork,ifail,info)
      else
        call dspgvx(1,'V','I','U',n,mat1,olap1,0d0,0d0,ir1,ir2,abstol,idum,eval,evec,n,work,iwork,ifail,info)
      endif
      ir2 = idum
      if(info/=0) Error('dspgvx failed.')
      deallocate ( work,iwork,mat1,olap1,ifail )
      end subroutine diagonalize_dpvox

      subroutine diagonalize_zvx(evec,eval,mat,ir1,ir2,r1,r2)
      implicit none
      real_dp,    intent(out)   :: eval(:)
      complex_dp, intent(out)   :: evec(:,:)
      complex_dp, intent(in)    :: mat(:,:)
      complex_dp, allocatable   :: work(:),mat1(:,:)
      real_dp,    allocatable   :: rwork(:)
      real_dp                   :: abstol,dlamch
      real_dp,    intent(in)    :: r1,r2
      integer,    intent(in)    :: ir1
      integer,    intent(inout) :: ir2
      integer,    allocatable   :: iwork(:),ifail(:)
      integer                   :: n,m,idum,info
      n = size(eval)
      m = ir2 - ir1 + 1
      if(n==0) Error('zero dimension in eigenvalue problem.')
      call check_dim(shape(mat),[n,n],'mat',__LINE__)
      if(r1>=r2) then
        if(m<0) Error('negative index range.')
        if(m>n) Error('number of selected eigenvalues exceeds maximal number.')
        call check_dim(shape(evec),[n,-m],'evec',__LINE__)
      else
        call check_dim(shape(evec),[n,n],'evec',__LINE__)
      endif
      allocate ( work(2*n),rwork(7*n),iwork(5*n),mat1(n,n),ifail(n) ) ; mat1 = mat
      abstol = 2 * dlamch('S')
      ! ZHEEVX calls ZLANHE. In older Intel MKL libraries, this function has a potential bug which can lead to a segfault.
      if(r1<r2) then
        call zheevx('V','V','U',n,mat1,n,r1,r2,0,0,      abstol,idum,eval,evec,n,work,2*n,rwork,iwork,ifail,info)
      else
        call zheevx('V','I','U',n,mat1,n,0d0,0d0,ir1,ir2,abstol,idum,eval,evec,n,work,2*n,rwork,iwork,ifail,info)
      endif
      ir2 = idum
      if(info/=0) Error('zheevx failed.')
      deallocate ( work,rwork,iwork,mat1,ifail )
      end subroutine diagonalize_zvx

      subroutine diagonalize_zvox(evec,eval,mat,olap,ir1,ir2,r1,r2)
      implicit none
      real_dp,    intent(out)   :: eval(:)
      complex_dp, intent(out)   :: evec(:,:)
      complex_dp, intent(in)    :: mat(:,:),olap(:,:)
      complex_dp, allocatable   :: work(:),mat1(:,:),olap1(:,:)
      real_dp,    allocatable   :: rwork(:)
      real_dp                   :: abstol,dlamch
      real_dp,    intent(in)    :: r1,r2
      integer,    intent(in)    :: ir1
      integer,    intent(inout) :: ir2
      integer,    allocatable   :: iwork(:),ifail(:)
      integer                   :: n,m,idum,info
      n = size(eval)
      m = ir2 - ir1 + 1
      if(n==0) Error('zero dimension in eigenvalue problem.')
      call check_dim(shape(mat), [n,n],'mat',__LINE__)
      call check_dim(shape(olap),[n,n],'olap',__LINE__)
      if(r1>=r2) then
        if(m<0) Error('negative index range.')
        if(m>n) Error('number of selected eigenvalues exceeds maximal number.')
        call check_dim(shape(evec),[n,-m],'evec',__LINE__)
      else
        call check_dim(shape(evec),[n,n],'evec',__LINE__)
      endif
      allocate ( work(2*n),rwork(7*n),iwork(5*n),mat1(n,n),olap1(n,n),ifail(n) ) ; mat1 = mat ; olap1 = olap
      abstol = 2 * dlamch('S')
      if(r1<r2) then
        call zhegvx(1,'V','V','U',n,mat1,n,olap1,n,r1,r2,0,0,      abstol,idum,eval,evec,n,work,2*n,rwork,iwork,ifail,info)
      else
        call zhegvx(1,'V','I','U',n,mat1,n,olap1,n,0d0,0d0,ir1,ir2,abstol,idum,eval,evec,n,work,2*n,rwork,iwork,ifail,info)
      endif
      ir2 = idum
      if(info/=0) Error('zhegvx failed.')
      deallocate ( work,rwork,iwork,mat1,olap1,ifail )
      end subroutine diagonalize_zvox

      subroutine diagonalize_zpvx(evec,eval,mat,ir1,ir2,r1,r2)
      implicit none
      real_dp,    intent(out)   :: eval(:)
      complex_dp, intent(out)   :: evec(:,:)
      complex_dp, intent(in)    :: mat(:)
      complex_dp, allocatable   :: work(:),mat1(:)
      real_dp,    allocatable   :: rwork(:)
      real_dp                   :: abstol,dlamch
      real_dp,    intent(in)    :: r1,r2
      integer,    intent(in)    :: ir1
      integer,    intent(inout) :: ir2
      integer,    allocatable   :: iwork(:),ifail(:)
      integer                   :: n,nn,m,idum,info
      n  = size(eval)
      nn = n*(n+1)/2
      m  = ir2 - ir1 + 1
      if(n==0) Error('zero dimension in eigenvalue problem.')
      call check_dim(shape(mat),[nn],'mat',__LINE__)
      if(r1>=r2) then
        if(m<0) Error('negative index range.')
        if(m>n) Error('number of selected eigenvalues exceeds maximal number.')
        call check_dim(shape(evec),[n,-m],'evec',__LINE__)
      else
        call check_dim(shape(evec),[n,n],'evec',__LINE__)
      endif
      allocate ( work(2*n),rwork(7*n),iwork(5*n),mat1(nn),ifail(n) ) ; mat1 = mat
      abstol = 2 * dlamch('S')
      if(r1<r2) then
        call zhpevx('V','V','U',n,mat1,r1,r2,0,0,      abstol,idum,eval,evec,n,work,iwork,ifail,info)
      else
        call zhpevx('V','I','U',n,mat1,0d0,0d0,ir1,ir2,abstol,idum,eval,evec,n,work,iwork,ifail,info)
      endif
      ir2 = idum
      if(info/=0) Error('zhpevx failed.')
      deallocate ( work,rwork,iwork,mat1,ifail )
      end subroutine diagonalize_zpvx

      subroutine diagonalize_zpvox(evec,eval,mat,olap,ir1,ir2,r1,r2)
      implicit none
      real_dp,    intent(out)   :: eval(:)
      complex_dp, intent(out)   :: evec(:,:)
      complex_dp, intent(in)    :: mat(:),olap(:)
      complex_dp, allocatable   :: work(:),mat1(:),olap1(:)
      real_dp,    allocatable   :: rwork(:)
      real_dp                   :: abstol,dlamch
      real_dp,    intent(in)    :: r1,r2
      integer,    intent(in)    :: ir1
      integer,    intent(inout) :: ir2
      integer,    allocatable   :: iwork(:),ifail(:)
      integer                   :: n,nn,m,idum,info
      n  = size(eval)
      nn = n*(n+1)/2
      m  = ir2 - ir1 + 1
      if(n==0) Error('zero dimension in eigenvalue problem.')
      call check_dim(shape(mat), [nn],'mat',__LINE__)
      call check_dim(shape(olap),[nn],'olap',__LINE__)
      if(r1>=r2) then
        if(m>n) Error('number of selected eigenvalues exceeds maximal number.')
        if(m<0) Error('negative index range.')
        call check_dim(shape(evec),[n,-m],'evec',__LINE__)
      else
        call check_dim(shape(evec),[n,n],'evec',__LINE__)
      endif
      allocate ( work(2*n),rwork(7*n),iwork(5*n),mat1(nn),olap1(nn),ifail(n) ) ; mat1 = mat ; olap1 = olap
      abstol = 2 * dlamch('S')
      if(r1<r2) then
        call zhpgvx(1,'V','V','U',n,mat1,olap1,r1,r2,0,0,      abstol,idum,eval,evec,n,work,rwork,iwork,ifail,info)
      else
        call zhpgvx(1,'V','I','U',n,mat1,olap1,0d0,0d0,ir1,ir2,abstol,idum,eval,evec,n,work,rwork,iwork,ifail,info)
      endif
      ir2 = idum
      if(info/=0) Error('zhpgvx failed.')
      deallocate ( work,rwork,iwork,mat1,olap1,ifail )
      end subroutine diagonalize_zpvox

c     --------

      subroutine diagonalize_gen_e(eval,mat)
      implicit none
      complex_dp, intent(in)  :: mat(:,:)
      complex_dp, intent(out) :: eval(:)
      complex_dp, allocatable :: help(:,:),work(:)
      real_dp,    allocatable :: rwork(:)
      integer                 :: n,info,pnt(size(eval))
      n = size(eval)
      call check_dim(shape(mat),[n,n],'mat',__LINE__)
      allocate ( help(n,n),work(4*n),rwork(2*n) )
      help = mat
      call zgeev('N','N',n,help,n,eval,mat,n,mat,n,work,4*n,rwork,info)
      deallocate(help,work,rwork)
      if(info/=0) Error('zgeev failed.')
      call rorderp(pnt,real(eval),n)
      eval = eval(pnt)
      end subroutine diagonalize_gen_e

      subroutine diagonalize_gen_v(evec,eval,mat)
      implicit none
      complex_dp, intent(in)  :: mat(:,:)
      complex_dp, intent(out) :: evec(:,:),eval(:)
      complex_dp, allocatable :: help(:,:),work(:)
      real_dp,    allocatable :: rwork(:)
      integer                 :: n,info
      n = size(eval)
      call check_dim(shape(mat), [n,n],'mat',__LINE__)
      call check_dim(shape(evec),[n,n],'evec',__LINE__)
      allocate ( help(n,n),work(4*n),rwork(2*n) )
      help = mat
      call zgeev('N','V',n,help,n,eval,evec,n,evec,n,work,4*n,rwork,info)
      deallocate(help,work,rwork)
      if(info/=0) Error('zgeev failed.')
      end subroutine diagonalize_gen_v

# if 0
      subroutine diagonalize_gen2(evecl,evecr,eval,mat)
      implicit none
      complex_dp, intent(in)  :: mat(:,:)
      complex_dp, intent(out) :: evecl(:,:),evecr(:,:),eval(:)
      complex_dp, allocatable :: help(:,:),work(:)
      real_dp,    allocatable :: rwork(:)
      integer                 :: n,info
      n = size(eval)
      call check_dim(shape(mat),  [n,n],'mat',__LINE__)
      call check_dim(shape(evecl),[n,n],'evecl',__LINE__)
      call check_dim(shape(evecr),[n,n],'evecr',__LINE__)
      allocate ( help(n,n),work(4*n),rwork(2*n) )
      help = mat
      call zgeev('V','V',n,help,n,eval,evecl,n,evecr,n,work,4*n,rwork,info)
      if(info/=0) Error('zgeev failed.')
      end subroutine diagonalize_gen2
# endif

c     --------

      subroutine geteigen_zpvo(evec,eval,mat,olap)
      implicit none
      real_dp,    intent(out) :: eval
      complex_dp, intent(out) :: evec(:)
      complex_dp, intent(in)  :: mat(:),olap(:)
      complex_dp, allocatable :: mat1(:),olap1(:),work(:),evec1(:,:)
      real_dp,    allocatable :: eval1(:),rwork(:)
      integer,    allocatable :: iwork(:),ifail(:)
      integer                 :: n,nn,info
      real_dp                 :: dlamch
      Error('disabled!')
      n = size(evec) ; nn = n*(n+1)/2
      call check_dim(shape(mat), [nn],'mat',__LINE__)
      call check_dim(shape(olap),[nn],'olap',__LINE__)
      allocate ( mat1(nn),olap1(nn),eval1(n),evec1(n,n),work(2*n),rwork(7*n),iwork(5*n),ifail(n) ) ; mat1 = mat ; olap1 = olap
      call zhpgvx(1,'V','I','U',n,mat1,olap1,0d0,0d0,n,n,2*dlamch('S'),1,eval1,evec1,n,work,rwork,iwork,ifail,info)
      if(info/=0) Error('zhpgvx failed.')
      evec = evec1(:,1)
      eval = eval1(1)
      deallocate ( mat1,olap1,eval1,evec1,work,rwork,iwork,ifail )
      end subroutine geteigen_zpvo

c     --------

      subroutine inverse_d(mati,mat)
      implicit none
      real_dp, intent(out) :: mati(:,:)
      real_dp, intent(in)  :: mat(:,:)
      real_dp, allocatable :: work(:)
      integer, allocatable :: ipiv(:)
      integer              :: n,info
      n = size(mat,1) ; if(n==0) return
      call check_dim(shape(mat), [n,n],'mat',__LINE__)
      call check_dim(shape(mati),[n,n],'mati',__LINE__)      
      mati = mat ; allocate(ipiv(n),work(n))
      call dgetrf(n,n,mati,n,ipiv,info)      ; if(info/=0) Error('dgetrf failed.')
      call dgetri(n,mati,n,ipiv,work,n,info) ; if(info/=0) Error('dgetri failed.')
      deallocate(ipiv,work)
      end subroutine inverse_d

      subroutine inverse_dp(mati,mat)
      implicit none
      real_dp, intent(out) :: mati(:)
      real_dp, intent(in)  :: mat(:)
      integer              :: n,nn,info
      nn = size(mat,1) ; n = nint(sqrt(0.25d0+2*nn)-0.5d0) ; if(n==0) return
      if(n*(n+1)-2*nn/=0) Bug('Input matrix has wrong size.')
      call check_dim(shape(mati),[nn],'mati',__LINE__)
      mati = mat
      call dpptrf('U',n,mati,info) ; if(info/=0) Error('dpptrf failed.')
      call dpptri('U',n,mati,info) ; if(info/=0) Error('dpptri failed.')
      end subroutine inverse_dp

      subroutine inverse_z(mati,mat)
      implicit none
      complex_dp, intent(out) :: mati(:,:)
      complex_dp, intent(in)  :: mat(:,:)
      complex_dp, allocatable :: work(:)
      integer,    allocatable :: ipiv(:)
      integer                 :: n,info
      n = size(mat,1) ; if(n==0) return
      call check_dim(shape(mat), [n,n],'mat',__LINE__)
      call check_dim(shape(mati),[n,n],'mati',__LINE__)
      mati = mat ; allocate(ipiv(n),work(n))
      call zgetrf(n,n,mati,n,ipiv,info)      ; if(info/=0) Error('zgetrf failed.')
      call zgetri(n,mati,n,ipiv,work,n,info) ; if(info/=0) Error('zgetri failed.')
      deallocate(ipiv,work)
      end subroutine inverse_z

      subroutine inverse_zp(mati,mat)
      implicit none
      complex_dp, intent(out) :: mati(:)
      complex_dp, intent(in)  :: mat(:)
      integer                 :: n,nn,info
      nn = size(mat,1) ; n = nint(sqrt(0.25d0+2*nn)-0.5d0) ; if(n==0) return
      if(n*(n+1)-2*nn/=0) Bug('Input matrix has wrong size.')
      call check_dim(shape(mati),[nn],'mati',__LINE__)
      mati = mat
      call zpptrf('U',n,mati,info) ; if(info/=0) Error('zpptrf failed.')
      call zpptri('U',n,mati,info) ; if(info/=0) Error('zpptri failed.')
      end subroutine inverse_zp

      subroutine inverse_d1(mat)
      implicit none
      real_dp, intent(inout) :: mat(:,:)
      real_dp, allocatable   :: work(:)
      integer, allocatable   :: ipiv(:)
      integer                :: n,info
      n = size(mat,1) ; if(n==0) return
      call check_dim(shape(mat),[n,n],'mat',__LINE__)
      allocate ( ipiv(n),work(n) )
      call dgetrf(n,n,mat,n,ipiv,info)      ; if(info/=0) Error('dgetrf failed.')
      call dgetri(n,mat,n,ipiv,work,n,info) ; if(info/=0) Error('dgetri failed.')
      deallocate(ipiv,work)
      end subroutine inverse_d1

      subroutine inverse_dp1(mat)
      implicit none
      real_dp, intent(inout) :: mat(:)
      integer                :: n,nn,info
      nn = size(mat,1) ; n = nint(sqrt(0.25d0+2*nn)-0.5d0) ; if(n==0) return
      if(n*(n+1)-2*nn/=0) Bug('Input matrix has wrong size.')
      call dpptrf('U',n,mat,info) ; if(info/=0) Error('dpptrf failed.')
      call dpptri('U',n,mat,info) ; if(info/=0) Error('dpptri failed.')
      end subroutine inverse_dp1

      subroutine inverse_z1(mat)
      implicit none
      complex_dp, intent(inout) :: mat(:,:)
      complex_dp, allocatable   :: work(:)
      integer,    allocatable   :: ipiv(:)
      integer                   :: n,info
      n = size(mat,1) ; if(n==0) return
      call check_dim(shape(mat),[n,n],'mat',__LINE__)
      allocate(ipiv(n),work(n))
      call zgetrf(n,n,mat,n,ipiv,info)      ; if(info/=0) Error('zgetrf failed.')
      call zgetri(n,mat,n,ipiv,work,n,info) ; if(info/=0) Error('zgetri failed.')
      deallocate(ipiv,work)
      end subroutine inverse_z1

      subroutine inverse_zp1(mat)
      implicit none
      complex_dp, intent(inout) :: mat(:)
      complex_dp, allocatable   :: work(:)
      integer,    allocatable   :: ipiv(:)
      integer                   :: n,nn,info
      nn = size(mat,1) ; n = nint(sqrt(0.25d0+2*nn)-0.5d0) ; if(n==0) return
      if(n*(n+1)-2*nn/=0) Bug('Input matrix has wrong size.')
      allocate(ipiv(n),work(n))
      call zhptrf('U',n,mat,ipiv,info)      ; if(info/=0) Error('zpptrf failed.')
      call zhptri('U',n,mat,ipiv,work,info) ; if(info/=0) Error('zpptri failed.')
      deallocate(ipiv,work)
      end subroutine inverse_zp1

c     --------

      function invert_d(mat)
      implicit none
      real_dp, intent(in)  :: mat(:,:)
      real_dp                 invert_d(size(mat,1),size(mat,1))
      call inverse(invert_d,mat)
      end function invert_d

      function invert_dp(mat)
      implicit none
      real_dp, intent(in)  :: mat(:)
      real_dp                 invert_dp(size(mat))
      call inverse(invert_dp,mat)
      end function invert_dp

      function invert_z(mat)
      implicit none
      complex_dp, intent(in)  :: mat(:,:)
      complex_dp                 invert_z(size(mat,1),size(mat,1))
      call inverse(invert_z,mat)
      end function invert_z

      function invert_zp(mat)
      implicit none
      complex_dp, intent(in)  :: mat(:)
      complex_dp                 invert_zp(size(mat))
      call inverse(invert_zp,mat)
      end function invert_zp

c     --------

      ! Return square root of mat (assumed to be hermitian). The result is hermitian if mat is real, otherwise complex.

      subroutine sqrtmat_d(matout,matin)
      implicit none
      real_dp, intent(out) :: matout(:,:)
      real_dp, intent(in)  :: matin(:,:)
      real_dp, allocatable :: eval(:),evec(:,:)
      integer              :: n,i
      n = size(matin,1)
      call check_dim(shape(matin), [n,n],'matin',__LINE__)
      call check_dim(shape(matout),[n,n],'matout',__LINE__)
      allocate ( evec(n,n),eval(n) )
      call diagonalize(evec,eval,matin) ; if(any(eval<0d0)) Error('negative eigenvalue.')
      do i = 1,n
        evec(:,i) = sqrt(sqrt(eval(i))) * evec(:,i)
      enddo
      matout = matmul(evec,transpose(evec))
      deallocate ( evec,eval )
      end subroutine sqrtmat_d

      subroutine sqrtmat_dp(matout,matin)
      implicit none
      real_dp, intent(out) :: matout(:)
      real_dp, intent(in)  :: matin(:)
      real_dp, allocatable :: eval(:),evec(:,:)
      integer              :: nn,n,i
      nn = size(matin,1) ; n = nint(sqrt(0.25d0+2*nn)-0.5d0) ; if(n*(n+1)-2*nn/=0) Bug('Input matrix has wrong size.')
      call check_dim(shape(matout),[nn],'matout',__LINE__)
      allocate ( evec(n,n),eval(n) )
      call diagonalize(evec,eval,matin) ; if(any(eval<0d0)) Error('negative eigenvalue.')
      do i = 1,n
        evec(:,i) = sqrt(sqrt(eval(i))) * evec(:,i)
      enddo
      matout = packmat ( matmul(evec,transpose(evec)) )
      deallocate ( evec,eval )
      end subroutine sqrtmat_dp

      subroutine sqrtmat_z(matout,matin)
      implicit none
      complex_dp, intent(out) :: matout(:,:)
      complex_dp, intent(in)  :: matin(:,:)
      complex_dp, allocatable :: evec(:,:)
      real_dp,    allocatable :: eval(:)
      integer                 :: n
      n = size(matin,1)
      call check_dim(shape(matin), [n,n],'matin',__LINE__)
      call check_dim(shape(matout),[n,n],'matout',__LINE__)
      allocate ( evec(n,n),eval(n) )
      call diagonalize(evec,eval,matin)
      matout = matmul(evec,diagmat(sqrt((1d0,0d0)*eval),conjg(transpose(evec))))
      deallocate ( evec,eval )
      end subroutine sqrtmat_z

      subroutine sqrtmat_zp(matout,matin)
      implicit none
      complex_dp, intent(out) :: matout(:)
      complex_dp, intent(in)  :: matin(:)
      complex_dp, allocatable :: evec(:,:)
      real_dp,    allocatable :: eval(:)
      integer                 :: nn,n
      nn = size(matin,1) ; n = nint(sqrt(0.25d0+2*nn)-0.5d0) ; if(n*(n+1)-2*nn/=0) Bug('Input matrix has wrong size.')
      call check_dim(shape(matout),[nn],'matout',__LINE__)
      allocate ( eval(n),evec(n,n) )
      call diagonalize(evec,eval,matin)
      matout = packmat( matmul(evec,diagmat(sqrt((1d0,0d0)*eval),conjg(transpose(evec)))) )
      deallocate ( evec,eval )
      end subroutine sqrtmat_zp

      subroutine sqrtmat_d1(mat)
      implicit none
      real_dp, intent(inout) :: mat(:,:)
      real_dp, allocatable   :: eval(:),evec(:,:)
      integer                :: n,i
      n = size(mat,1)
      call check_dim(shape(mat),[n,n],'mat',__LINE__)
      allocate ( evec(n,n),eval(n) )
      call diagonalize(evec,eval,mat) ; if(any(eval<0d0)) Error('negative eigenvalue.')
      do i = 1,n
        evec(:,i) = sqrt(sqrt(eval(i))) * evec(:,i)
      enddo
      mat = matmul(evec,transpose(evec))
      deallocate ( evec,eval )
      end subroutine sqrtmat_d1

      subroutine sqrtmat_dp1(mat)
      implicit none
      real_dp, intent(inout) :: mat(:)
      real_dp, allocatable   :: eval(:),evec(:,:)
      integer                :: nn,n,i
      nn = size(mat,1) ; n = nint(sqrt(0.25d0+2*nn)-0.5d0) ; if(n*(n+1)-2*nn/=0) Bug('Input matrix has wrong size.')
      allocate ( evec(n,n),eval(n) )
      call diagonalize(evec,eval,mat) ; if(any(eval<0d0)) Error('negative eigenvalue.')
      do i = 1,n
        evec(:,i) = sqrt(sqrt(eval(i))) * evec(:,i)
      enddo
      mat = packmat( matmul(evec,transpose(evec)) )
      deallocate ( eval,evec )
      end subroutine sqrtmat_dp1

      subroutine sqrtmat_z1(mat)
      implicit none
      complex_dp, intent(inout) :: mat(:,:)
      complex_dp, allocatable   :: evec(:,:)
      real_dp,    allocatable   :: eval(:)
      integer                   :: n
      n = size(mat,1)
      call check_dim(shape(mat),[n,n],'mat',__LINE__)
      allocate ( evec(n,n),eval(n) )
      call diagonalize(evec,eval,mat)
      mat = matmul(evec,diagmat(sqrt((1d0,0d0)*eval),conjg(transpose(evec))))
      deallocate ( eval,evec )
      end subroutine sqrtmat_z1

      subroutine sqrtmat_zp1(mat)
      implicit none
      complex_dp, intent(inout) :: mat(:)
      complex_dp, allocatable   :: evec(:,:)
      real_dp,    allocatable   :: eval(:)
      integer                   :: nn,n
      nn = size(mat,1) ; n = nint(sqrt(0.25d0+2*nn)-0.5d0) ; if(n*(n+1)-2*nn/=0) Bug('Input matrix has wrong size.')
      allocate ( eval(n),evec(n,n) )
      call diagonalize(evec,eval,mat) ; if(any(eval<0d0)) Error('negative eigenvalue.')
      mat = packmat( matmul(evec,diagmat(sqrt((1d0,0d0)*eval),conjg(transpose(evec)))) )
      deallocate ( evec,eval )
      end subroutine sqrtmat_zp1

c     --------

# if 0
      subroutine inverse_sqrtmat(matout,matin)
      implicit none
      complex_dp, intent(out) :: matout(:,:)
      complex_dp, intent(in)  :: matin(:,:)
      complex_dp, allocatable :: evecl(:,:),evecr(:,:)
      complex_dp, allocatable :: eval(:)
      integer                 :: n,i
      n = size(matin,1)
      call check_dim(shape(matin), [n,n],'matin',__LINE__)
      call check_dim(shape(matout),[n,n],'matout',__LINE__)
      allocate ( evecl(n,n),evecr(n,n),eval(n) )
      call diagonalize_gen2(evecl,evecr,eval,matin)
      do i = 1,n
        evecr(:,i) = evecr(:,i) /  sqrt(eval(i))
      enddo
      matout = matmul(evecr,conjg(transpose(evecl)))
      deallocate ( evecr,evecl,eval )
      end subroutine inverse_sqrtmat
# endif

c     --------

      function sqrtm_d(mat)
      implicit none
      real_dp, intent(in)  :: mat(:,:)
      real_dp              :: sqrtm_d(size(mat,1),size(mat,1))
      call sqrtmat(sqrtm_d,mat)
      end function sqrtm_d

      function sqrtm_dp(mat)
      implicit none
      real_dp, intent(in)  :: mat(:)
      real_dp              :: sqrtm_dp(size(mat))
      call sqrtmat(sqrtm_dp,mat)
      end function sqrtm_dp

      function sqrtm_z(mat)
      implicit none
      complex_dp, intent(in)  :: mat(:,:)
      complex_dp              :: sqrtm_z(size(mat,1),size(mat,2))
      call sqrtmat(sqrtm_z,mat)
      end function sqrtm_z

      function sqrtm_zp(mat)
      implicit none
      complex_dp, intent(in)  :: mat(:)
      complex_dp              :: sqrtm_zp(size(mat))
      call sqrtmat(sqrtm_zp,mat)
      end function sqrtm_zp

c     --------

      subroutine solve_d(vecout,mat,vecin)
      implicit none
      real_dp, intent(out) :: vecout(:)
      real_dp, intent(in)  :: mat(:,:),vecin(:)
      real_dp              :: mat1(size(mat,1),size(mat,2))
      integer              :: n,m,info,ipiv(size(vecout))
      n = size(vecout)
      m = min(n,size(vecin))
      call check_dim(shape(mat),[n,n],'mat',__LINE__)
      vecout(:m)   = vecin(:m)
      vecout(m+1:) = 0d0
      mat1         = mat
      call dgesv(n,1,mat1,n,ipiv,vecout,n,info)
      if(info/=0) Error('dgesv failed.')
      end subroutine solve_d

      subroutine solve_z(vecout,mat,vecin)
      implicit none
      complex_dp, intent(out) :: vecout(:)
      complex_dp, intent(in)  :: mat(:,:),vecin(:)
      complex_dp              :: mat1(size(mat,1),size(mat,2))
      integer                 :: n,m,info,ipiv(size(vecout))
      n = size(vecout)
      m = min(n,size(vecin))
      call check_dim(shape(mat),[n,n],'mat',__LINE__)
      vecout(:m)   = vecin(:m)
      vecout(m+1:) = 0d0
      mat1         = mat
      call zgesv(n,1,mat1,n,ipiv,vecout,n,info)
      if(info/=0) Error('zgesv failed.')
      end subroutine solve_z

c     --------

      subroutine svd(mat1,diag,mat2,mat)
      implicit none
      real_dp,    intent(out) :: diag(:)
      complex_dp, intent(out) :: mat1(:,:),mat2(:,:)
      complex_dp, intent(in)  :: mat(:,:)
      complex_dp              :: hlp(size(mat,1),size(mat,2)),work(3*(size(mat,1)+size(mat,2)))
      real_dp                 :: rwork(5*min(size(mat,1),size(mat,2)))
      integer                 :: m,n,info
      m   = size(mat,1)
      n   = size(mat,2)
      hlp = mat
      call check_dim(shape(mat1),[m,m],'mat1',__LINE__)
      call check_dim(shape(mat2),[n,n],'mat2',__LINE__)
      if(size(diag)<min(m,n)) Error('array for diagonal elements too small.')
      call zgesvd('A','A',m,n,hlp,m,diag,mat1,m,mat2,n,work,3*(m+n),rwork,info)
      if(info/=0) Error('zgesvd failed.')
      end subroutine svd

c     --------

      function trace_d(mat)
      implicit none
      real_dp             :: trace_d
      real_dp, intent(in) :: mat(:,:)
      integer             :: i,n
      n       = size(mat,1) ; call check_dim(shape(mat),[n,n],'mat',__LINE__)
      trace_d = 0
      do i = 1,n
        trace_d = trace_d + mat(i,i)
      enddo
      end function trace_d

      function trace_z(mat)
      implicit none
      complex_dp             :: trace_z
      complex_dp, intent(in) :: mat(:,:)
      integer                :: i,n
      n       = size(mat,1) ; call check_dim(shape(mat),[n,n],'mat',__LINE__)
      trace_z = 0
      do i = 1,n
        trace_z = trace_z + mat(i,i)
      enddo
      end function trace_z

      function trace_dp(mat)
      implicit none
      real_dp             :: trace_dp
      real_dp, intent(in) :: mat(:)
      integer             :: i,j,n
      n        = size(mat) ; n = nint(sqrt(0.25d0+2*n)-0.5d0) ; if(n*(n+1)-2*size(mat)/=0) Bug('Input matrix has wrong size.')
      j        = 0
      trace_dp = 0
      do i = 1,n
        j        = j + i
        trace_dp = trace_dp + mat(j)
      enddo
      if(j/=size(mat)) Error('matrix has wrong size.')
      end function trace_dp

      function trace_zp(mat)
      implicit none
      complex_dp             :: trace_zp
      complex_dp, intent(in) :: mat(:)
      integer                :: i,j,n
      n        = size(mat) ; n = nint(sqrt(0.25d0+2*n)-0.5d0) ; if(n*(n+1)-2*size(mat)/=0) Bug('Input matrix has wrong size.')
      j        = 0
      trace_zp = 0
      do i = 1,n
        j        = j + i
        trace_zp = trace_zp + mat(j)
      enddo
      if(j/=size(mat)) Error('matrix has wrong size.')
      end function trace_zp

c     --------

      subroutine check_dim(dims1,dims0,name,line,file)
      use util
      implicit none
      integer,      intent(in)           :: dims0(:),dims1(:),line
      character(*), intent(in)           :: name
      character(*), intent(in), optional :: file
      character(80)                      :: file1
      integer                            :: i
      if(present(file)) then ; file1 = ' File: '//trim(file)
      else                   ; file1 = ' File: wrapper.f'
      endif
      if(size(dims0)/=size(dims1)) Bug('Array '//name//' has wrong number of dimensions.'//trim(file1))
      do i = 1,size(dims0)
        if(dims0(i)<0) then
          if(dims1(i)<-dims0(i)) Bug('Dimension '//Chr(i)//' of array '//name//' too small. Line: '//Chr(line)//trim(file1))
        else
          if(dims1(i)/=dims0(i)) Bug('Dimension '//Chr(i)//' of array '//name//' is wrong. Line: '//Chr(line)//trim(file1))
        endif          
      enddo
      end subroutine check_dim

c     --------

      function which_ranges(l1,l2,r1,r2,win)
      implicit none
      integer                       :: which_ranges
      logical, intent(in)           :: l1,l2,r1,r2
      logical, intent(in), optional :: win
      if((l1.neqv.l2).or.(r1.neqv.r2)) Bug('Optional ranges not given pairwise.')
      if(l1) then
        if(r1) then ; which_ranges = 3
        else        ; which_ranges = 1
        endif
      else
        if(r1) then ; which_ranges = 2
        else        ; which_ranges = 0
        endif
      endif
      if(present(win)) then
        if(win) which_ranges = which_ranges + 4
      endif
      end function which_ranges

c     --------

# ifdef MPI
      
c     Wrapper functions for ScaLAPACK

c     use Mwrapper   ! force wrapper to depend on Mwrapper
c     use timer_util ! force wrapper to depend on timer_util

#   include "w_scalapack.inc"

c Returns process grid nrow x ncol <= nprocs for array (m,n). (nrow=ncol or nrow+1=ncol)
c (If m==0, the array size is ignored.)
      subroutine grid_proc(nrow,ncol,nprocs,m,n)
      use, intrinsic :: iso_fortran_env
      implicit none
      integer, intent(out) :: nrow,ncol
      integer, intent(in)  :: nprocs,m,n
      integer              :: i
# ifdef scalapack_test
#   warning small loc for testing
      real_dp, parameter   :: loc=1000d0
# else
      real_dp, parameter   :: loc=500000d0 ! minimal local portion (e.g., 500000)
# endif
      i = 2
      do
        i    = i + 1
        nrow = i/2
        ncol = (i+1)/2
        if(nrow*ncol>nprocs)                 exit
        if(1d0*m*n/(nrow*ncol)<loc.and.m/=0) exit
        if(m>0.and.m<=blk.and.nrow>1)        exit ! avoid descinit failure
      enddo
      nrow = (i-1)/2
      ncol = i/2
      end subroutine grid_proc

# endif

c
c     --------

      end module wrapper

c     --------

      subroutine redistribute_cvec(vec,n)
      use, intrinsic :: iso_fortran_env
      implicit none
      integer, intent(in)    :: n
      real_dp, intent(inout) :: vec(2*n)
      real_dp                :: re(n),im(n)
      integer                :: i
      re = vec(:n)
      im = vec(n+1:)
      do i = 1,n
        vec(2*i-1) =  re(i)
        vec(2*i)   = -im(i)
      enddo
      end

c ------------------------------
# endif
c ------------------------------
