c--------------------------------------------------------------------------------
c Copyright (c) 2021 Peter Grünberg Institut, Forschungszentrum Jülich, Germany
c This file is part of SPEX and available as free software under the conditions
c of the MIT license as expressed in the LICENSE file in more detail.
c--------------------------------------------------------------------------------

c        packmat ( mat ) : returns matrix mat in packed storage (upper triangle)
c call p_packmat ( mat ) : --- " --- (more stack-economical)

# ifndef DZtyp
#   define DZtyp 1
#   include "w_packmat.inc"
#   define SUB
#   include "w_packmat.inc"
#   undef DZtyp      
#   define DZtyp 2
#   include "w_packmat.inc"
#   define SUB
# endif

# if DZtyp == 1      
#   define DZ d
#   define Datatyp real_dp
# elif DZtyp == 2
#   define DZ z
#   define Datatyp complex_dp
# else
#   error DZtyp unknown     
# endif
       
# ifdef SUB
      subroutine p_packmat_ DZ (out,mat,ltest)
# else      
      function packmat_ DZ (mat,ltest) result(out)
# endif
      use util, only: chr
      implicit none
      logical, intent(in), optional :: ltest
      Datatyp, intent(in)           :: mat(:,:)
      Datatyp                       :: out(size(mat,1)*(size(mat,1)+1)/2)
      real_dp                       :: rdum
      real_dp                       :: r1,r2,i1,i2
      logical                       :: ltest1
      integer                       :: n,i,j,k
      ltest1 = .true.
      if(present(ltest)) ltest1 = ltest
      n = size(mat,1)
      call check_dim(shape(mat),[n,n],'mat',__LINE__)
      k = 0
      do j = 1,n
        do i = 1,j
          k      = k + 1
          out(k) = mat(i,j)
          if(ltest1) then
# if DZtyp == 1
            r1   = mat(i,j)
            r2   = mat(j,i)
            rdum = abs(r1) + abs(r2)
            if(rdum>1) then ; rdum = ( abs(r1-r2) ) / rdum
            else            ; rdum = ( abs(r1-r2) )
            endif
            if     (rdum>1d-6) then ; Error('Matrix not symmetric. Deviation: '//Chf(rdum,'ES8.1'))
            else if(rdum>1d-8) then ; Warn('Matrix deviates from symmetry by more than 10^-8.')
# else
            r1   = real(mat(i,j)) ; i1 =  imag(mat(i,j))
            r2   = real(mat(j,i)) ; i2 = -imag(mat(j,i))
            rdum = abs(r1) + abs(r2) + abs(i1) + abs(i2)
            if(rdum>1) then ; rdum = ( abs(r1-r2) + abs(i1-i2) ) / rdum
            else            ; rdum = ( abs(r1-r2) + abs(i1-i2) )
            endif
            if     (rdum>1d-6) then ; Error('Matrix not Hermitian. Deviation: '//Chf(rdum,'ES8.1'))
            else if(rdum>1d-8) then ; Warn('Matrix deviates from hermiticity by more than 10^-8.')
# endif
            endif
          endif
        enddo
      enddo

# ifdef SUB
      end subroutine p_packmat_ DZ      
# else      
      end function packmat_ DZ
# endif

c     ------

# ifdef SUB
      subroutine p_unpackmat_ DZ (out,mat,ltest)
# else      
      function unpackmat_ DZ (mat,ltest) result(out)
# endif      
      implicit none      
      logical, intent(in), optional :: ltest
      Datatyp, intent(in)           :: mat(:)
      Datatyp                       :: out(nint(sqrt(0.25d0+2*size(mat))-0.5d0),
     &                                     nint(sqrt(0.25d0+2*size(mat))-0.5d0))
      integer                       :: n,nn,i,j,k
      nn = size(mat) ; n = nint(sqrt(0.25d0+2*nn)-0.5d0) ; if(n*(n+1)-2*nn/=0) Bug('Input matrix has wrong size.')
      k  = 0
      do j = 1,n
        do i = 1,j
          k        = k + 1
          out(i,j) = mat(k)
# if DZtyp == 1
          out(j,i) = mat(k)
        enddo
# else          
          out(j,i) = conjg(mat(k))
        enddo
        if(present(ltest)) then
          if(ltest) then
            if(abs(imag(mat(k)))>1d-10) Error('Diagonal element of hermitian matrix not real.')
          endif
        endif
# endif        
      enddo
# ifdef SUB
      end subroutine p_unpackmat_ DZ
# else
      end function unpackmat_ DZ
# endif

# undef DZ
# undef Datatyp

# ifdef SUB
#   if DZtyp == 2      
#     undef DZtyp
#   endif
#   undef SUB
# endif
