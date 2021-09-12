c--------------------------------------------------------------------------------
c Copyright (c) 2021 Peter Grünberg Institut, Forschungszentrum Jülich, Germany
c This file is part of SPEX and available as free software under the conditions
c of the MIT license as expressed in the LICENSE file in more detail.
c--------------------------------------------------------------------------------

c     "Rotates" a mixed-basis quantity from k0 to k1 (or from k1 to k0)
c     matin(:dim1,:dim2) -> matout(:dim1,:dim2)
c
c     Routines
c       mtrafo      : Input and output arrays are complex-valued.
c       mtrafo_r    : Input and output arrays are real-valued. (Only -DINV, otherwise not compiled)
c       mtrafo_pw   : Same as mtrafo for PW part.
c       mtrafo_pw_r : Same as mtrafo_r for PW part. (Only -DINV, otherwise not compiled)
c       mtrafo1*    : "inplace" routines, i.e., input/output matrices are the same array.
c       mtrafo*_r   : Apart from rotation, result matrix is multiplied by phase factor to make it real (see KPT1ADD).
c      
c     Two types of rotation depending on the sign of isym0 [isym=abs(isym0), ikpt1=kptsym(ikpt0,isym)]
c     isym0 > 0 : matin at k0 (ikpt0) is        rotated to k1 (ikpt1)
c     isym0 < 0 : matin at k1 (ikpt1) is (back-)rotated to k0 (ikpt0)
c
c     Notation for symmetry isym: P = (A,a,C), R = (A,a) = spatial part
c                                 A = spatial rotation
c                                 a = translation
c                                 C = 1 (without time reversal) or conjg() (time reversal)
c
c     P acting on M:
c     MT and C=1:     P M_k0,almp(r) =   M_k0,almp[A^-1(r-a)] = exp(i g r_a') exp[-i (k'+g) a] SUM(m') M_k1,a'lm'p(r) D_m'm(A)  | D(A)  = Wigner matrix
c     MT and C=conjg: P M_k0,almp(r) = C M_k0,almp[A^-1(r-a)] = exp(i g r_a') exp[-i (k'+g) a] SUM(m') M_k1,a'lm'p(r) D_m'm(CA) | D(CA) = Wigner matrix with complex conjugation
c     PW and C=1:     P M_k0,G(r)    =   M_k0,G   [A^-1(r-a)] = exp[-i (k'+g+AG) a] M_k1,(g+AG)
c     PW and C=conjg: P M_k0,G(r)    = C M_k0,G   [A^-1(r-a)] = exp[-i (k'+g-AG) a] M_k1,(g-AG)
c     with k' = Ak - g and k' = -Ak - g for C=1/conjg, respectively. (dwgn and sym%rrot contain C and +/-1, respectively, already.)
c
c     We can generally write  P M_k0,I  = SUM(I') S(II')       M_k1,I'  | S is unitary
c                               M_k1,I' = SUM(I)  S*(II')    P M_k0,I
c                               M_k0,I  = SUM(I') CS(II') P^-1 M_k1,I'  | CS (without space): C acts on S only, C S (with space): C would act on S and the rest.
c                          P^-1 M_k1,I' = SUM(I)  CS*(II')     M_k0,I
c     with [the nonzero elements S(II')]
c     MT: S(almp,a'lm'p) = exp(i g r_a') exp[-i (k'+g) a] D_m'm(CA) | C=1/conjg with/without time reversal. Wigner matrix is unitary.
c     PW: S(G,g+-AG)     = exp[-i (k'+g+-AG) a]                     |        +- with/without time reversal
c
c
c     Trafos for isym0 > 0:
c
c     mode = 1 : Bra vector is transformed
c                matout(I') =           <      M_k1,I'|      B >                | Variable transformation with spatial part R=(A,a)
c                           =           < R^-1 M_k1,I'| R^-1 B >
c                           = C         < P^-1 M_k1,I'| P^-1 B >                | B assumed invariant wrt P
c                           = C         < P^-1 M_k1,I'|      B >
c                           = C SUM(I)  <      M_k0,I |      B > CS(II')
c                           =   SUM(I) C<      M_k0,I |      B >  S(II')
c                           =   SUM(I) Cmatin(I)                  S(II')
c
c     mode = 2 : Ket vector is transformed:
c                matout(I') =           < B | M_k1,I'>                          | Analogous derivation
c                           =   SUM(I) Cmatin(I) S*(II')
c
c     mode = 3 : Matrix is transformed:
c                matout(I'J') =         <       M_k1,I'| B |       M_k1,J'>     | B assumed invariant wrt spatial part R
c                             =         < CP^-1 M_k1,I'| B | CP^-1 M_k1,J'>     | assuming symmetry B(rr')=B(r'r) in case of time reversal
c                             = SUM(IJ) <       M_k0,I | B |       M_k0,J > * /  S(II') *  S*(JJ')  for C=1  (no time reversal)
c                                                                             \ CS(IJ') * CS*(JI')  for C=conjg (time reversal)
c                             = SUM(IJ) Tmatin(IJ) * S(II') * S*(JJ')           | T=transpose if C=conjg (time reversal), otherwise T=1
c
c     Trafos for isym0 < 0:
c
c     mode = 1 : Bra vector is transformed
c                matout(I) =           <   M_k0,I  |   B >                      | Variable transformation with spatial part R=(A,a)
c                          =           < R M_k0,I  | R B >
c                          = C         < P M_k0,I  | P B >                      | B is invariant wrt P
c                          = C         < P M_k0,I  |   B >
c                          = C SUM(I') <   M_k1,I' |   B > S*(II')
c                          = C SUM(I') matin(I')           S*(II')
c
c     mode = 2 : Ket vector is transformed:
c                matout(I) =           < B | M_k0,I  >                          | Analogous derivation
c                          = C SUM(I') matin(I') S(II')
c
c     mode = 3 : Matrix is transformed:
c                matout(IJ) =             <    M_k0,I  | B |    M_k0,J  >       | B assumed invariant wrt spatial part R
c                           =             < CP M_k0,I  | B | CP M_k0,J  >       | assuming symmetry B(rr')=B(r'r) in case of time reversal
c                           =   SUM(I'J') <    M_k1,I' | B |    M_k1,J' > * / S*(II') * S(JJ')  for C=1  (no time reversal)
c                                                                           \ S*(JI') * S(IJ')  for C=conjg (time reversal)
c                           = T SUM(I'J') Tmatin(I'J') * S*(II') * S(JJ')       | note T definition
c
c     mode =  1,  2,  3 : normal matrix dimensions (rows=bras, columns=kets)
c     mode = -1, -2, -3 : identical to mode = 1, 2, 3 but with flipped matrix dimensions (rows=kets, columns=bras)
c      
c     ladd = .true.     : result is added to matout (not available in "inplace" routines)
c

# ifndef SUBTYPE
#   define SUBTYPE 1
#   include "trafo.inc"
#   define SUBTYPE 3
#   ifdef INV
#     include "trafo.inc"      
#     define SUBTYPE 2
#     include "trafo.inc"
#     define SUBTYPE 4
#   endif
# endif

# if SUBTYPE == 1 || SUBTYPE == 3
#   define KPT1ADD kpt(:,ikpt1)+
# else
#   define KPT1ADD
# endif

# if SUBTYPE <= 2
#   define SUB12(arg) arg
#   define NBASPADD nbasp+
# else
#   define SUB12(arg)
#   define NBASPADD
# endif

# ifdef INPLACE
#   define SUF 1
#   define ARGS matout,dim1,dim2,ikpt0,isym0,mode0
#   define MATIN matout
#   define LADD .false.
# else
#   define SUF
#   define ARGS matout,matin,dim1,dim2,ikpt0,isym0,mode0,ladd
# endif

# if   SUBTYPE == 1
      subroutine mtrafo SUF       ( ARGS )
# elif SUBTYPE == 2
      subroutine mtrafo SUF _r    ( ARGS )
# elif SUBTYPE == 3
      subroutine mtrafo SUF _pw   ( ARGS )
# elif SUBTYPE == 4
      subroutine mtrafo SUF _pw_r ( ARGS )
# else
#   error "Unknown subtype"
# endif

      use global
      use, intrinsic :: iso_fortran_env

      implicit none
      integer,    intent(in)    :: ikpt0,isym0,dim1,dim2,mode0
# ifdef INPLACE
#   if SUBTYPE == 1 || SUBTYPE == 3
      complex_dp, intent(inout) :: matout(dim1,dim2)
#   else
      real_dp,    intent(inout) :: matout(dim1,dim2)
#   endif
# else
      logical,    intent(in)    :: ladd
#   if SUBTYPE == 1 || SUBTYPE == 3
      complex_dp, intent(in)    :: matin(dim1,dim2)
      complex_dp, intent(inout) :: matout(dim1,dim2)
#   else
      real_dp,    intent(in)    :: matin(dim1,dim2)
      real_dp,    intent(inout) :: matout(dim1,dim2)
#   endif
# endif
# if SUBTYPE <= 2
      complex_dp                :: dwgnt(-maxlcutm:maxlcutm,-maxlcutm:maxlcutm,0:maxlcutm),cexp
      integer                   :: itype,ieq,ic,ic1,l,nn,nlm,lm,lm1,m,m1
# endif
      complex_dp                :: matin1(dim1,dim2)
      complex_dp                :: matout1(dim1,dim2)
      complex_dp                :: cdum,carr(ngptm(ikpt0))
      integer                   :: mode,smode
      integer                   :: isym,ssym,ikpt1,g(3),g1(3),iarr(ngptm(ikpt0))
      integer                   :: i,i1,i2
      integer                   :: igptm,igptm2,igptp
      real                      :: time1

      call cpu_time(time1)

      if(all(mode0/=[-3,-2,-1,1,2,3])) Bug('mode unknown.')
      isym  = abs(isym0)
      mode  = abs(mode0)
      ssym  = sign(1,isym0)
      smode = sign(1,mode0) ; if(smode==-1.and.mode<3) mode = 3 - mode
      ikpt1 = kptsym(ikpt0,isym)
      g     = nint ( matmul(sym(isym)%rrot,kpt(:,ikpt0))-kpt(:,ikpt1) )
# if SUBTYPE <= 2
      do l = 0,maxlcutm
        if(ssym==1) then ; dwgnt(-l:l,-l:l,l) = conjg(transpose(dwgn(-l:l,-l:l,l,isym)))
        else             ; dwgnt(-l:l,-l:l,l) =                 dwgn(-l:l,-l:l,l,isym)
        endif
        if(smode==-1)      dwgnt(-l:l,-l:l,l) = conjg( dwgnt(-l:l,-l:l,l) )
      enddo
# endif

c     Check dims
      i = NBASPADD ngptm(ikpt0)
      if((mode==1.or.mode==3).and.i>dim1) Bug('First dimension too small.')
      if((mode==2.or.mode==3).and.i>dim2) Bug('Second dimension too small.')

      i1 = dim1 ; if(mode/=2) i1 = i
      i2 = dim2 ; if(mode/=1) i2 = i
# ifdef INV
      matin1(:i1,:i2) = MATIN(:i1,:i2)
      SUB12( call desymmetrize(matin1,dim1,dim2,mode0) )
# else
      if(isym0>nsymt) then
        if(mode==3) then ; matin1(:i,:i)   = transpose( MATIN (:i,:i))
        else             ; matin1(:i1,:i2) = conjg( MATIN (:i1,:i2))
        endif
      else               ; matin1(:i1,:i2) = MATIN (:i1,:i2)
      endif
# endif

# if SUBTYPE <= 2
      cexp = exp(ssym * smode * img * 2*pi * dot_product( KPT1ADD g,sym(isym)%transl))
# endif

c     Right-multiplication
      if(mode==2.or.mode==3) then
#       if SUBTYPE <= 2
        ! MT
        ic = 0
        lm = 0
        do itype = 1,ntype
          nlm = sum( [((2*l+1)*nindxm(l,itype),l=0,lcutm(itype))] )
          do ieq = 1,neq(itype)
            ic                       = ic + 1
            ic1                      = pcent(ic,isym)
            lm1                      = lm ; if(ssym==1) lm1 = lm1 + (ic1-ic) * nlm
            matout1(:,lm1+1:lm1+nlm) = 0
            i1                       = lm1
            do l = 0,lcutm(itype)
              nn = nindxm(l,itype)
              do m1 = -l,l
                i = lm ; if(ssym==-1) i = i + (ic1-ic) * nlm
                do m  = -l,l
                  cdum = dwgnt(m,m1,l) ; if(cdum.ne.0) matout1(:,i1+1:i1+nn) = matout1(:,i1+1:i1+nn) + matin1(:,i+1:i+nn) * cdum
                  i    = i + nn
                enddo
                i1 = i1 + nn
              enddo
              lm = lm + (2*l+1) * nn
            enddo
            cdum                     = cexp * exp(-ssym * smode * img * 2*pi * dot_product(g,cent(:,ic1)))
            matout1(:,lm1+1:lm1+nlm) = matout1(:,lm1+1:lm1+nlm) * cdum
          enddo
        enddo
#       endif
        ! PW
        do igptm = 1,ngptm(ikpt0)
          igptp  = pgptm(igptm,ikpt0)
          g1     = matmul(sym(isym)%rrot,gptm(:,igptp)) + g
          igptm2 = pntgptm(g1(1),g1(2),g1(3),ikpt1)
          cdum   = exp(ssym * smode * img * 2*pi * dot_product( KPT1ADD g1,sym(isym)%transl))
          if(igptm2==0) Error('G point not found in G-point set (right-multiplication).')
          if(ssym==1) then ; matout1(:,NBASPADD igptm2) = cdum * matin1(:,NBASPADD igptm)
          else             ; matout1(:,NBASPADD igptm)  = cdum * matin1(:,NBASPADD igptm2)
          endif
        enddo
        matin1 = matout1
      endif

c     Left-multiplication
      if(mode==1.or.mode==3) then
#       if SUBTYPE <= 2
        ! MT
        cexp  = conjg(cexp)
        dwgnt = conjg(dwgnt)
        ic    = 0
        lm    = 0
        do itype = 1,ntype
          nlm = sum( [((2*l+1)*nindxm(l,itype),l=0,lcutm(itype))] )
          do ieq = 1,neq(itype)
            ic                       = ic + 1
            ic1                      = pcent(ic,isym)
            lm1                      = lm ; if(ssym==1) lm1 = lm1 + (ic1-ic) * nlm
            matout1(lm1+1:lm1+nlm,:) = 0
            i1                       = lm1
            do l = 0,lcutm(itype)
              nn = nindxm(l,itype)
              do m1 = -l,l
                i = lm ; if(ssym==-1) i = i + (ic1-ic) * nlm
                do m  = -l,l
                  cdum = dwgnt(m,m1,l) ; if(cdum.ne.0) matout1(i1+1:i1+nn,:) = matout1(i1+1:i1+nn,:) + matin1(i+1:i+nn,:) * cdum
                  i    = i + nn
                enddo
                i1 = i1 + nn
              enddo
              lm = lm + (2*l+1) * nn
            enddo
            cdum                     = cexp * exp(ssym * smode * img * 2*pi * dot_product(g,cent(:,ic1)))
            matout1(lm1+1:lm1+nlm,:) = matout1(lm1+1:lm1+nlm,:) * cdum
          enddo
        enddo
#       endif
        ! PW
        do igptm = 1,ngptm(ikpt0)
          igptp  = pgptm(igptm,ikpt0)
          g1     = matmul(sym(isym)%rrot,gptm(:,igptp)) + g
          igptm2 = pntgptm(g1(1),g1(2),g1(3),ikpt1)
          if(igptm2==0) Error('G point not found in G-point set (left-multiplication).')
          iarr(igptm) = igptm2
          carr(igptm) = exp(-ssym * smode * img * 2*pi * dot_product( KPT1ADD g1,sym(isym)%transl))
        enddo
        do i2 = 1,dim2
          do i1 = 1,ngptm(ikpt0)
            if(ssym==1) then ; matout1(NBASPADD iarr(i1),i2) = carr(i1) * matin1(NBASPADD      i1, i2)
            else             ; matout1(NBASPADD      i1, i2) = carr(i1) * matin1(NBASPADD iarr(i1),i2)
            endif
          enddo
        enddo
      endif

      i  = NBASPADD ngptm(ikpt0)
      i1 = dim1 ; if(mode/=2) i1 = i
      i2 = dim2 ; if(mode/=1) i2 = i
# ifdef INV
      SUB12( call symmetrize(matout1,dim1,dim2,mode0) )
      if( LADD ) then ; matout(:i1,:i2) = matout(:i1,:i2) + matout1(:i1,:i2)
      else            ; matout(:i1,:i2) = matout1(:i1,:i2)
      endif
# else
      if(isym0<-nsymt) then
        if(mode==3) then
          if( LADD ) then ; matout(:i,:i) = matout(:i,:i) + transpose(matout1(:i,:i))
          else            ; matout(:i,:i) =                 transpose(matout1(:i,:i))
          endif
        else
          if( LADD ) then ; matout(:i1,:i2) = matout(:i1,:i2) + conjg(matout1(:i1,:i2))
          else            ; matout(:i1,:i2) =                   conjg(matout1(:i1,:i2))
          endif
        endif
      else
        if( LADD ) then ; matout(:i1,:i2) = matout(:i1,:i2) + matout1(:i1,:i2)
        else            ; matout(:i1,:i2) =                   matout1(:i1,:i2)
        endif
      endif
# endif

      end

c     ------------------

# undef KPT1ADD
# undef NBASPADD
# undef SUB12

# undef SUF
# undef ARGS

# ifdef INPLACE
#   undef MATIN
#   undef LADD
# endif

# undef SUBTYPE
