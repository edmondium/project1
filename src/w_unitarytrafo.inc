c--------------------------------------------------------------------------------
c Copyright (c) 2021 Peter Grünberg Institut, Forschungszentrum Jülich, Germany
c This file is part of SPEX and available as free software under the conditions
c of the MIT license as expressed in the LICENSE file in more detail.
c--------------------------------------------------------------------------------

c *** Unitary transformations ***
c
c call unitarytrafo( out , mat , ctrafo , m ) : returns unitary transformation of mat in out
c                                               (m=1: ctrafo**T * mat * ctrafo, m=2: ctrafo * mat * ctrafo**T)
c call unitarytrafo( mat , ctrafo , m )       : returns result in mat (in place)
c
c suffix  out  mat  ctrafo   (m=full matrix, p=packed matrix)
c -       m    m    m
c p       p    p    m
c 1       -    m    m
c p1      -    p    m
c p2      m    p    m
c zd      m    m    m
c zd1     -    m    m
c
c All matrices are assumed to have the same datatype (real_dp or complex_dp), except in routines ending
c on "zd" and "zd1", where the datatypes are complex_dp and real_dp for mat/out and ctrafo, respectively.     
c
c The result matrix out (or mat) can have dimensions larger than the actual result but not smaller.
c Routines ending on "1" return the result in mat. The leading dimension of mat must then not be smaller
c than the result matrix.      
c
c More stack-friendly and (usually faster) than the explicit matrix operations.
c Could be implemented more efficiently if BLAS (level 3) routines were complete
c (e.g., dsyrk2 and zherk2 with result in packed storage or operation A*B + B**T*A**T).
c
c
c Performance of routines p and p2 (valid for about n>50 and for both real and complex versions):
c
c p(m=1) < p2(m=1)
c p2(m=2) < p(m=2)
c p2(m=1) < p(m=1,2) < p2(m=2)
c

# ifndef DZtyp
#   define DZtyp 1
#   include "w_unitarytrafo.inc"
#   define INPLACE
#   include "w_unitarytrafo.inc"
#   undef DZtyp
#   define DZtyp 2
#   include "w_unitarytrafo.inc"
#   define INPLACE
#   include "w_unitarytrafo.inc"
#   undef DZtyp
#   define DZtyp 3
#   include "w_unitarytrafo.inc"
#   define INPLACE      
# endif

# if DZtyp == 1
#   define DZ d
#   define SH s
#   define SH2 sy
#   define Im0      
#   define Datatyp real_dp
#   define Conjg dble
# elif DZtyp == 2
#   define DZ z
#   define SH h
#   define SH2 he
#   define Im0 ,0d0      
#   define Datatyp complex_dp
#   define Conjg conjg
# elif DZtyp != 3
#   error DZtyp unknown      
# endif

# ifdef INPLACE
#   define Output mat
# else
#   define Output out
# endif

# if DZtyp != 3

# ifdef INPLACE
      subroutine unitarytrafo_ DZ 1 (mat,ctrafo,m)
      implicit none
      Datatyp, intent(inout) :: mat(:,:)
# else
      subroutine unitarytrafo_ DZ (out,mat,ctrafo,m)      
      implicit none
      Datatyp, intent(out) :: out(:,:)
      Datatyp, intent(in)  :: mat(:,:)
# endif
      Datatyp, intent(in)  :: ctrafo(:,:)
      integer, intent(in)  :: m
      Datatyp, allocatable :: help(:,:)
      integer              :: n,n1
      if(all(m/=[1,2])) Bug('Last argument must 1 or 2.')
      n  = size(mat,1)
      n1 = size(ctrafo,3-m)
# ifdef INPLACE
      if(n1>n) Bug('Output matrix must not be larger than input matrix for inplace routines.')
# else
      call check_dim(shape(out),[-n1,n1],'out',__LINE__,__FILE__)
# endif      
      call check_dim(shape(mat),[n,n],  'mat',__LINE__,__FILE__)
      if(m==1) then
        call check_dim(shape(ctrafo),[n,n1],'ctrafo',__LINE__,__FILE__)
        allocate ( help(n,n1) )
        call DZ gemm('N','N',n,n1,n,(1d0 Im0),mat,n,ctrafo,n,(0d0 Im0),help,n)
        call DZ gemm('C','N',n1,n1,n,(1d0 Im0),ctrafo,n,help,n,(0d0 Im0), Output ,size( Output ,1))
      else
        call check_dim(shape(ctrafo),[n1,n],'ctrafo',__LINE__,__FILE__)
        allocate ( help(n1,n) )
        call DZ gemm('N','N',n1,n,n,(1d0 Im0),ctrafo,n1,mat,n,(0d0 Im0),help,n1)
        call DZ gemm('N','C',n1,n1,n,(1d0 Im0),help,n1,ctrafo,n1,(0d0 Im0), Output ,size( Output ,1))
      endif
      deallocate(help)
# ifdef INPLACE
      end subroutine unitarytrafo_ DZ 1
# else
      end subroutine unitarytrafo_ DZ
# endif

c     -------

# ifdef INPLACE
      subroutine unitarytrafo_ DZ p1 (mat,ctrafo,m)
      implicit none
      Datatyp, intent(inout) :: mat(:)
# else
      subroutine unitarytrafo_ DZ p(out,mat,ctrafo,m)      
      implicit none      
      Datatyp, intent(out) :: out(:)
      Datatyp, intent(in)  :: mat(:)
# endif      
      Datatyp, intent(in)  :: ctrafo(:,:)
      integer, intent(in)  :: m
      Datatyp, allocatable :: help(:,:),out1(:,:)
      integer, parameter   :: step = 16, minpack = 16
      integer              :: n,n1,n2,pack
      integer              :: i,j,k,l
      if(all(m/=[1,2])) Bug('Last argument must 1 or 2.')
      n  = nint(sqrt(0.25d0+2*size(mat))-0.5d0)
      n1 = size(ctrafo,3-m)
# ifdef INPLACE
      if(n1>n) Bug('Output matrix must not be larger than input matrix.')
# else
      call check_dim(shape(out),[-n1*(n1+1)/2],'out',__LINE__,__FILE__)     
# endif
      call check_dim(shape(mat),[n*(n+1)/2],'mat',__LINE__,__FILE__)
      pack = max(minpack,n1/step+1)
      allocate( out1(pack,n1) )
      if(m==1) then
        call check_dim(shape(ctrafo),[n,n1],'ctrafo',__LINE__,__FILE__)        
        allocate( help(n,n1) )
        do i = 1,n1
          call DZ SH pmv('U',n,(1d0 Im0),mat,ctrafo(:,i),1,(0d0 Im0),help(1,i),1)
        enddo
        do i = 1,n1,pack
          n2 = min(pack,n1-i+1)
          call DZ gemm('C','N',n2,n1-i+1,n,(1d0 Im0),ctrafo(1,i),n,help(1,i),n,(0d0 Im0),out1,pack)
          do j = 1,n1-i+1
            l = (j+i-1)*(j+i-2)/2 + i
            do k = 1,min(j,pack)
              Output(l) = out1(k,j)
              l         = l + 1                            
            enddo
          enddo
        enddo
      else
        call check_dim(shape(ctrafo),[n1,n],'ctrafo',__LINE__,__FILE__)
        allocate( help(n,n1) )
        do i = 1,n1
          call DZ SH pmv('U',n,(1d0 Im0),mat,Conjg(ctrafo(i,:)),1,(0d0 Im0),help(1,i),1)
        enddo
        do i = 1,n1,pack
          n2 = min(pack,n1-i+1)
          call DZ gemm('N','N',n2,n1-i+1,n,(1d0 Im0),ctrafo(i,1),n1,help(1,i),n,(0d0 Im0),out1,pack)
          do j = 1,n1-i+1
            l = (j+i-1)*(j+i-2)/2 + i
            do k = 1,min(j,pack)
              Output(l) = out1(k,j)
              l         = l + 1                            
            enddo
          enddo
        enddo
      endif
      deallocate(help,out1)
# ifdef INPLACE
      end subroutine unitarytrafo_ DZ p1
# else
      end subroutine unitarytrafo_ DZ p
# endif

c     -------

# ifndef INPLACE

      subroutine unitarytrafo_ DZ p2(out,mat,ctrafo,m)      
      implicit none      
      Datatyp, intent(out) :: out(:,:)
      Datatyp, intent(in)  :: mat(:)
      Datatyp, intent(in)  :: ctrafo(:,:)
      integer, intent(in)  :: m
      Datatyp, allocatable :: help(:,:)
      integer              :: n,n1,n2,pack
      integer              :: i,j,k,l
      if(all(m/=[1,2])) Bug('Last argument must 1 or 2.')
      n  = nint(sqrt(0.25d0+2*size(mat))-0.5d0)
      n1 = size(ctrafo,3-m)
      call check_dim(shape(out),[-n1,n1],'out',__LINE__,__FILE__)
      call check_dim(shape(mat),[n*(n+1)/2],'mat',__LINE__,__FILE__)
      if(m==1) then
        allocate( help(n,n1) )
        call check_dim(shape(ctrafo),[n,n1],'ctrafo',__LINE__,__FILE__)
        do i = 1,n1
          call DZ SH pmv('U',n,(1d0 Im0),mat,ctrafo(:,i),1,(0d0 Im0),help(1,i),1)
        enddo
        call DZ gemm('C','N',n1,n1,n,(1d0 Im0),ctrafo,n,help,n,(0d0 Im0),out,n1)
      else
        allocate( help(n1,n) )
        call check_dim(shape(ctrafo),[n1,n],'ctrafo',__LINE__,__FILE__)
        j = 1
        k = 1
        do i = 1,n
          help(:,i) = ctrafo(:,i)
          call DZ gemv('N',n1,k-j+1,(1d0 Im0),ctrafo,n1,mat(j),1,-mat(k)/2,help(1,i),1)
          j = k + 1
          k = j + i
        enddo
        call DZ SH2 r2k('U','N',n1,n,(1d0 Im0),help,n1,ctrafo,n1,(0d0 Im0),out,n1)
        do j = 1,n1
          do i = j+1,n1
            out(i,j) = Conjg( out(j,i) )
          enddo
        enddo
      endif
      deallocate(help)
      end subroutine unitarytrafo_ DZ p2

# endif

# undef DZ
# undef SH
# undef SH2
# undef Im0
# undef Datatyp
# undef Conjg

# else

# ifdef INPLACE
      subroutine unitarytrafo_zd1 (mat,ctrafo,m)
      implicit none
      complex_dp, intent(inout) :: mat(:,:)
# else
      subroutine unitarytrafo_zd (out,mat,ctrafo,m)      
      implicit none
      complex_dp, intent(out) :: out(:,:)
      complex_dp, intent(in)  :: mat(:,:)
# endif
      real_dp,    intent(in)  :: ctrafo(:,:)
      integer,    intent(in)  :: m
      real_dp,    allocatable :: help(:,:),vec(:)
      integer                 :: n,n1,n2
      integer                 :: i,j,j1
      if(all(m/=[1,2])) Bug('Last argument must 1 or 2.')
      n  = size(mat,1)
      n1 = size(ctrafo,3-m)
      n2 = size( Output ,1)
# ifdef INPLACE
      if(n1>n) Bug('Output matrix must not be larger than input matrix')
# else
      call check_dim(shape(out),[-n1,n1],'out',__LINE__,__FILE__)
# endif
      call check_dim(shape(mat),[n,n],  'mat',__LINE__,__FILE__)
      allocate ( help(2*n,n1),vec(2*max(n,n1)) )
      if(m==1) then
        call check_dim(shape(ctrafo),[n,n1],'ctrafo',__LINE__,__FILE__)
        call dgemm('N','N',2*n,n1,n,1d0,mat,2*n,ctrafo,n,0d0,help,2*n)
      else
        call check_dim(shape(ctrafo),[n1,n],'ctrafo',__LINE__,__FILE__)
        call dgemm('N','C',2*n,n1,n,1d0,mat,2*n,ctrafo,n1,0d0,help,2*n)
      endif
      do i = 1,n1
        vec(1) = help(2,i)
        do j = 2,n
          help(j,i)   = help(2*j-1,i)
          vec(j)      = help(2*j,i)
        enddo
        do j = 1,n
          help(n+j,i) = vec(j)
        enddo
      enddo
      if(m==1) then ; call dgemm('C','N',n1,2*n1,n,1d0,ctrafo,n, help,n,0d0, Output ,n2)
      else          ; call dgemm('N','N',n1,2*n1,n,1d0,ctrafo,n1,help,n,0d0, Output ,n2)
      endif
      do i = 1,n1
        if(n2/=n1) then          
          do j = 1,n1/2
            vec(2*j-1) = real(Output(j,i))
            vec(2*j)   = imag(Output(j,i))
          enddo
          if(mod(n1,2)==1) vec(n1) = real(Output(n1/2+1,i))
          if(mod(n2,2)==0) then
            j1 = n2/2
            do j = 1,n1/2
              vec(n1+2*j-1) = real(Output(j+j1,i))
              vec(n1+2*j)   = imag(Output(j+j1,i))
            enddo
            if(mod(n1,2)==1) vec(2*n1) = real(Output(n1/2+1+j1,i))
          else
            j1        = n2/2+1
            vec(n1+1) = imag(Output(j1,i))
            do j = 1,(n1-1)/2
              vec(n1+2*j)   = real(Output(j+j1,i))
              vec(n1+2*j+1) = imag(Output(j+j1,i))
            enddo
            if(mod(n1,2)==0) vec(2*n1) = real(Output(n1/2+j1,i))
          endif
        else
          do j = 1,n1
            vec(2*j-1) = real(Output(j,i))
            vec(2*j)   = imag(Output(j,i))
          enddo
        endif
        do j = 1,n1
          Output(j,i) = cmplx( vec(j) , vec(n1+j) , kind(1d0) )
        enddo
      enddo
      deallocate ( help,vec )
# ifdef INPLACE
      end subroutine unitarytrafo_zd1
# else
      end subroutine unitarytrafo_zd
# endif

# endif

# undef Output

# ifdef INPLACE
#   if DZtyp == 3
#     undef DZtyp
#   endif
#   undef INPLACE
# endif

