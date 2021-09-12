c--------------------------------------------------------------------------------
c Copyright (c) 2021 Peter Grünberg Institut, Forschungszentrum Jülich, Germany
c This file is part of SPEX and available as free software under the conditions
c of the MIT license as expressed in the LICENSE file in more detail.
c--------------------------------------------------------------------------------

c These routines transform quantities (wave functions, matrices, etc.) at k-point ikpt0 through a symmetry operation isym
c to obtain the corresponding quantities at the rotated k-point ikpt1=kptsym(ikpt0,isym).

c --------------

# include "cppmacro.h"

# ifndef old_trafo

c --------------

c Overview for wavefunction transformation routines:
c
c                   coeff.  SOC  ic    description
c waveftrafo_mt_io  in/out  no   yes   transforms MT coefficients in real space
c waveftrafo_pw_io  in/out  no   -     transforms PW coefficients in real space
c waveftrafo2       in/out  no   no    transforms    coefficients of nbnd bands (only called by getinput)
c waveftrafo_soc    in/out  yes  -     transforms    coefficients in spin space (isym)
c waveftrafo1_soc   in/out  yes  -     transforms    coefficients in spin space (SO2 matrix "esoc" provided)
c waveftrafo        out     yes  no    transforms    coefficients of (iband,ikpt0,ispin)
c waveftrafo_mt     out     yes  yes   transforms    coefficients of (iband,ikpt0,ispin), only MT coefficients
c waveftrafo_pw     out     yes  -     transforms    coefficients of (iband,ikpt0,ispin), only PW coefficients
c wavefunction      out     yes  no          gets    coefficients of (iband,ikpt,ispin)
c wavefunction_mt   out     yes  yes         gets MT coefficients of (iband,ikpt,ispin)
c wavefunction_pw   out     yes  -           gets PW coefficients of (iband,ikpt,ispin)
c
c "SOC" : rotation in spin space included
c "ic"  : atom center (type) selectable with ic>0 (ic<0)

c --------------

c Backend routines
c
c waveftrafo_mt_io: Transforms cmtin from ikpt0 to ikpt = kptsym(ikpt0,isym) by symmetry operation isym => cmtout
c waveftrafo_pw_io:            cpwin                                                                    => cpwout
c
c ikpt00 < 0 : (in-place routine) cmtout is also the input matrix. (kpoint index = abs(ikpt00))
c      
c (No rotations in spin space in case of l_soc.)
c
c --------------
c
c waveftrafo_mt_io:
c ic : if ic=0, transform all atom centers
c      if ic>0, transform only atom center ic to pcent(ic,isym) (coefficients then in cmtio(:,1,:))
c      if ic<0, transform only atom type  -ic                   (coefficients then in cmtio(:,1:neq(-ic),:))
c
c --------------
c Note: Explicit interfaces for waveftrafo_mt_io (and waveftrafo_pw_io) included only for GFortran. Arguably not required by Fortran standard.
c --------------
c      
c begin interface
      subroutine waveftrafo_mt_io(cmtout,cmtin,dimout,dimin,ic,isym,ikpt00)

      use global

      use, intrinsic :: iso_fortran_env !inc
      implicit none
      integer,    intent(in)         :: ikpt00,isym,ic,dimout,dimin
      complex_dp, intent(in), target :: cmtin(dimin,*)
      complex_dp, intent(inout)      :: cmtout(dimout,*)
c end interface
      complex_dp, pointer_cnt        :: cmtin1(:,:)
      complex_dp                     :: cmthlp(2*maxlcut+1)
      complex_dp                     :: cexp
      ifInv(real_dp,integer)         :: tr(3)
      real_dp                        :: dp
      integer                        :: ikpt0,l,n,nn,lm0,lm1,lm2
      integer                        :: itype,ieq,icent,icent0,icent1,ikpt
      logical                        :: trl NoInvC(trs)

      ikpt0 = abs(ikpt00)
      if(ic> ncent) Bug('ic>ncent')
      if(ic<-ntype) Bug('ic<-ntype')
      if(ikpt0<0.or.ikpt0>nkpt2) Bug('ikpt0 out of range.')
      if(isym <0.or.isym>nsym)   Bug('isym out of range.')
      if     (ic==0) then ; if(ikpt00<0) then ; allocate(cmtin1(dimout,ncent))    ; else ; cmtin1 => cmtin(:,:ncent)    ; endif
      else if(ic> 0) then ; if(ikpt00<0) then ; allocate(cmtin1(dimout,1))        ; else ; cmtin1 => cmtin(:,:1)        ; endif
      else                ; if(ikpt00<0) then ; allocate(cmtin1(dimout,neq(-ic))) ; else ; cmtin1 => cmtin(:,:neq(-ic)) ; endif
      endif
      if(ikpt00<0) cmtin1 = cmtout(:,:size(cmtin1,2))
      NoInv( trs = isym>nsymt )
      ikpt  = kptsym(ikpt0,isym)
      icent = 0
      do itype = 1,ntype
        do ieq = 1,neq(itype)
          icent  = icent + 1
          icent1 = pcent(icent,isym)
          if(ic==0.or.ic==icent1.or.ic==-itype) then
            if     (ic==0) then ; icent0 = icent
            else if(ic >0) then ; icent0 = 1                          ; icent1 = 1
            else if(ic <0) then ; icent0 = icent - sum(neq(:itype-1)) ; icent1 = icent1 - sum(neq(:itype-1))
            endif
            tr   = tcent(:,icent,isym) Inv(-sym(isym)%transl)
            dp   = dot_product(kpt(:,ikpt),tr)
            trl  = dp/=0 ; if(trl) cexp = exp( -img * (2*pi) * dp )
            lm0  = 0
            do l = 0,lcut(itype)
              nn = nindx(l,itype) ; if(lm0+nn*(2*l+1)>dimout) exit
              do n = 1,nn
                lm1 = lm0 + n
                lm2 = lm0 + n + 2*l*nn
# ifdef INV
                cmthlp(:2*l+1) = cmtin1(lm1:lm2:nn,icent0)
# else
                if(trs) then ; cmthlp(:2*l+1) = conjg ( cmtin1(lm1:lm2:nn,icent0) )
                else         ; cmthlp(:2*l+1) =         cmtin1(lm1:lm2:nn,icent0)
                endif
# endif
                if(trl) then ; cmtout(lm1:lm2:nn,icent1) = cexp * matmul(dwgn(-l:l,-l:l,l,isym),cmthlp(:2*l+1))
                else         ; cmtout(lm1:lm2:nn,icent1) =        matmul(dwgn(-l:l,-l:l,l,isym),cmthlp(:2*l+1))
                endif
              enddo
              lm0 = lm2
            enddo
            if(lm2<dimout) cmtout(lm2+1:,icent1) = 0
          endif          
        enddo
      enddo
      
      if(ikpt00<0) deallocate(cmtin1)

      end

c --------------

c begin interface
      subroutine waveftrafo_pw_io(cpwout,cpwin,isym,ikpt00)

      use global

      use, intrinsic :: iso_fortran_env !inc
      implicit none
      integer,       intent(in)         :: ikpt00,isym
      MCOMPLEX_dp,   intent(in), target :: cpwin(*)
      MCOMPLEX_dp,   intent(inout)      :: cpwout(*)
c end interface
      MCOMPLEX_dp,   pointer_cnt        :: cpwin1(:)
      logical                           :: trl NoInvC(trs)
      integer                           :: g(3),g1(3) InvC(t(3))
      integer                           :: ikpt0,igpt,igpt0,igpt1,isymi,ikpt
      integer, save, allocatable        :: pgptinv(:)
      integer, save                     :: ikpt0_ = 0
      NoInv( complex_dp                 :: cexp )

      if(.not.allocated(pgptinv)) allocate(pgptinv(ngptall))

      ikpt0 = abs(ikpt00)
      if(ikpt0<0.or.ikpt0>nkpt2) Bug('ikpt0 out of range.')
      if(isym <0.or.isym>nsym)   Bug('isym out of range.')
      if(ikpt00<0) then ; allocate(cpwin1(ngpt(ikpt0))) ; cpwin1 = cpwout(:ngpt(ikpt0))
      else              ; cpwin1 => cpwin(:ngpt(ikpt0))
      endif
      ikpt  = kptsym(ikpt0,isym)
      isymi = sym(isym)%inv
      trl   = any(sym(isym)%transl/=0)
# ifdef INV
      t     = nint(2*sym(isym)%transl) ; if(any(abs(2*sym(isym)%transl-t)>1d-12)) Bug('2*transl is not a lattice vector.')
# else
      trs   = isym>nsymt
# endif

      if(ikpt0/=ikpt0_) then
        pgptinv                           = 0
        pgptinv(pgpt(:ngpt(ikpt0),ikpt0)) = [ (igpt,igpt=1,ngpt(ikpt0)) ]
        ikpt0_                            = ikpt0
      endif

      g = nint ( matmul(sym(isym)%rrot,kpt(:,ikpt0))-kpt(:,ikpt) ) ! ROT*kp = k + g,  g = rec. latt. vector
      do igpt1 = 1,ngpt(ikpt)
        igpt  = pgpt(igpt1,ikpt) !; if(sum(matmul(rlat,kpt(:,ikpt)+gpt(:,igpt))**2)>gcut**2) Bug('|k+G|>gcut.')
        g1    = matmul( sym(isymi)%rrot,gpt(:,igpt)-g )
        igpt0 = pgptinv ( pntgpt(g1(1),g1(2),g1(3)) ) ; if(igpt0==0) Bug('zero pointer')
# ifdef INV
        if(trl) then
          if(mod(dot_product(gpt(:,igpt),t),2)/=0) then ; cpwout(igpt1) = -cpwin1(igpt0) ! factor exp(-iGa)=+-1 [a=sym(i)%transl]
          else                                          ; cpwout(igpt1) =  cpwin1(igpt0)
          endif
        else
          cpwout(igpt1) = cpwin1(igpt0)
        endif
# else
        if(trl) then
          cexp  = exp(-img * 2*pi * dot_product (kpt(:,ikpt)+gpt(:,igpt), sym(isym)%transl))
          if(trs) then ; cpwout(igpt1) = cexp * conjg ( cpwin1(igpt0) )
          else         ; cpwout(igpt1) = cexp *         cpwin1(igpt0)
          endif
        else
          if(trs) then ; cpwout(igpt1) = conjg ( cpwin1(igpt0) )
          else         ; cpwout(igpt1) =         cpwin1(igpt0)
          endif
        endif
# endif
      enddo
      if(ikpt00<0) deallocate(cpwin1)

      end

c --------------

c     Carbon copy of waveftrafo_pw_io for complex arrays (needed by wavefproducts3_pw)
# ifdef INV
      subroutine waveftrafo_pw_io_w(cpwout,cpwin,isym,ikpt0)

      use global

      use, intrinsic :: iso_fortran_env
      implicit none
      integer,    intent(in)  :: ikpt0,isym
      complex_dp, intent(in)  :: cpwin(*)
      complex_dp, intent(out) :: cpwout(*)
      logical                 :: trs,trl
      integer                 :: g(3),g1(3),pgptinv(ngptall)
      integer                 :: igpt,igpt0,igpt1,isymi,ikpt
      complex_dp              :: cexp

      if(isym==1) return
      if(ikpt0<0.or.ikpt0>nkpt2) Bug('ikpt0 out of range.')
      if(isym <0.or.isym>nsym)   Bug('isym out of range.')
      ikpt  = kptsym(ikpt0,isym)      
      isymi = sym(isym)%inv
      trl   = any(sym(isym)%transl/=0)
      trs   = isym>nsymt

      pgptinv                           = 0
      pgptinv(pgpt(:ngpt(ikpt0),ikpt0)) = [ (igpt,igpt=1,ngpt(ikpt0)) ]

      g = nint ( matmul(sym(isym)%rrot,kpt(:,ikpt0))-kpt(:,ikpt) ) ! ROT*kp = k + g,  g = rec. latt. vector
      do igpt1 = 1,ngpt(ikpt)
        igpt  = pgpt(igpt1,ikpt) !; if(sum(matmul(rlat,kpt(:,ikpt)+gpt(:,igpt))**2)>gcut**2) Bug('|k+G|>gcut.')
        g1    = matmul( sym(isymi)%rrot,gpt(:,igpt)-g )
        igpt0 = pgptinv ( pntgpt(g1(1),g1(2),g1(3)) ) ; if(igpt0==0) Bug('zero pointer')
        if(trl) then
          cexp  = exp(-img * 2*pi * dot_product (kpt(:,ikpt)+gpt(:,igpt), sym(isym)%transl))
          if(trs) then ; cpwout(igpt1) = cexp * conjg ( cpwin(igpt0) )
          else         ; cpwout(igpt1) = cexp *         cpwin(igpt0)
          endif
        else
          if(trs) then ; cpwout(igpt1) = conjg ( cpwin(igpt0) )
          else         ; cpwout(igpt1) =         cpwin(igpt0)
          endif
        endif
      enddo

      end
# endif

c --------------

c Special routine called by getinput only: wavefunctions are transformed for nbnd bands.
c (Faster than calling waveftrafo for each band.)
      subroutine waveftrafo2(cmtout,cmtin,cpwout,cpwin,nbnd,isym,ikpt0)

      use global

      use, intrinsic :: iso_fortran_env
      implicit none
      integer,     intent(in)  :: ikpt0,isym,nbnd
      complex_dp,  intent(in)  :: cmtin(maxlmindx,ncent,nbnd)
      complex_dp,  intent(out) :: cmtout(maxlmindx,ncent,nbnd)
      MCOMPLEX_dp, intent(in)  :: cpwin(maxgpt,nbnd)
      MCOMPLEX_dp, intent(out) :: cpwout(maxgpt,nbnd)
      complex_dp               :: cmthlp(2*maxlcut+1)
      complex_dp               :: cexp
      logical                  :: trl NoInvC(trs)
      ifInv(real_dp,integer)   :: tr(3)
      integer                  :: pgptinv(ngptall)
      integer                  :: l,n,nn,lm0,lm1,lm2 InvC(t(3))
      integer                  :: itype,ieq,icent,icent1,ikpt,ibnd
      integer                  :: g(3),g1(3)
      integer                  :: igpt,igpt0,igpt1,isymi

      if(ikpt0<0.or.ikpt0>nkpt2) Bug('ikpt0 out of range.')
      if(isym <0.or.isym>nsym)   Bug('isym out of range.')
      NoInv( trs = isym>nsymt )
      ikpt  = kptsym(ikpt0,isym)
      isymi = sym(isym)%inv

      icent = 0
      do itype = 1,ntype
        do ieq = 1,neq(itype)
          icent  = icent + 1
          icent1 = pcent(icent,isym)
          tr     = tcent(:,icent,isym) Inv(-sym(isym)%transl)
          trl    = any(tr/=0) ; if(trl) cexp = exp( -img * (2*pi) * dot_product(kpt(:,ikpt),tr) )
          lm0    = 0
          do l = 0,lcut(itype)
            nn = nindx(l,itype)
            do n = 1,nn
              lm1 = lm0 + n
              lm2 = lm0 + n + 2*l*nn
              do ibnd = 1,nbnd
# ifdef INV
                cmthlp(:2*l+1) = cmtin(lm1:lm2:nn,icent,ibnd)
# else
                if(trs) then ; cmthlp(:2*l+1) = conjg ( cmtin(lm1:lm2:nn,icent,ibnd) )
                else         ; cmthlp(:2*l+1) =         cmtin(lm1:lm2:nn,icent,ibnd)
                endif
# endif
                if(trl) then ; cmtout(lm1:lm2:nn,icent1,ibnd) = cexp * matmul(dwgn(-l:l,-l:l,l,isym),cmthlp(:2*l+1))
                else         ; cmtout(lm1:lm2:nn,icent1,ibnd) =        matmul(dwgn(-l:l,-l:l,l,isym),cmthlp(:2*l+1))
                endif
              enddo
            enddo
            lm0 = lm2
          enddo
        enddo
      enddo

# ifdef INV
      t = nint(2*sym(isym)%transl) ; if(any(abs(2*sym(isym)%transl-t)>1d-12)) Bug('2*transl is not a lattice vector.')
# else
      trl = any(sym(isym)%transl/=0)
# endif

      pgptinv                           = 0
      pgptinv(pgpt(:ngpt(ikpt0),ikpt0)) = [ (igpt,igpt=1,ngpt(ikpt0)) ]

      g = nint ( matmul(sym(isym)%rrot,kpt(:,ikpt0))-kpt(:,ikpt) ) ! ROT*kp = k + g,  g = rec. latt. vector
      do igpt1 = 1,ngpt(ikpt)
        igpt  = pgpt(igpt1,ikpt)
        g1    = matmul( sym(isymi)%rrot,gpt(:,igpt)-g )
        igpt0 = pgptinv ( pntgpt(g1(1),g1(2),g1(3)) ) ; if(igpt0==0) Bug('zero pointer')
# ifdef INV
        trl   = mod(dot_product(gpt(:,igpt),t),2)/=0
# else
        if(trl) cexp = exp(-img * 2*pi * dot_product (kpt(:,ikpt)+gpt(:,igpt), sym(isym)%transl))
# endif
        do ibnd = 1,nbnd
# ifdef INV
          if(trl) then ; cpwout(igpt1,ibnd) = -cpwin(igpt0,ibnd) ! factor exp(-iGa)=+-1 [a=sym(i)%transl]
          else         ; cpwout(igpt1,ibnd) =  cpwin(igpt0,ibnd)
          endif
# else
          if(trl) then            
            if(trs) then ; cpwout(igpt1,ibnd) = cexp * conjg ( cpwin(igpt0,ibnd) )
            else         ; cpwout(igpt1,ibnd) = cexp *         cpwin(igpt0,ibnd)
            endif
          else
            if(trs) then ; cpwout(igpt1,ibnd) = conjg ( cpwin(igpt0,ibnd) )
            else         ; cpwout(igpt1,ibnd) =         cpwin(igpt0,ibnd)
            endif
          endif
# endif
        enddo
      enddo

      end

c --------------      

c Rotation in spin space
c (waveftrafo_soc:)  through symmetry operation isym
c (waveftrafo1_soc:) through SO2 rotation matrix esoc.
c
      subroutine waveftrafo_soc(c,n,isym)
      use global, only: sym,l_soc
      use, intrinsic :: iso_fortran_env
      implicit none
      integer,    intent(in)    :: n,isym
      complex_dp, intent(inout) :: c(n,2)
      integer                   :: i
      if(.not.l_soc) Bug('waveftrafo_soc called but lsoc=.false.')
      do i = 1,n
        c(i,:) = matmul ( c(i,:) , sym(isym)%esoc )
      enddo
      end

      subroutine waveftrafo1_soc(c,n,esoc)
      use global, only: l_soc
      use, intrinsic :: iso_fortran_env
      implicit none
      integer,    intent(in)    :: n
      complex_dp, intent(in)    :: esoc(2,2)
      complex_dp, intent(inout) :: c(n,2)
      integer                   :: i
      if(.not.l_soc) Bug('waveftrafo1_soc called but l_soc=.false.')
      do i = 1,n
        c(i,:) = matmul ( c(i,:) , esoc )
      enddo
      end

c     ------------------

c Rotate wavefunction from ikpt0 through symmetry operation to kptsym(ikpt0,isym).
c
c ispin : if ispin == 0, transform all spins (from 1 to nspin2),
c         otherwise only transform spin ispin, the resulting coefficients are then in cmtout(:,:,1) and/or cpwout(:,1)
c l_soc : transform always all spins (for ispin==0 and 1, otherwise error message),
c         spin rotation is included
c
c waveftrafo_mt : only MT coefficients (atom center can be selected)
c waveftrafo_pw : only PW coefficients
c waveftrafo    : all coefficients
c      
      subroutine waveftrafo(cmtout,cpwout,isym,iband,ikpt0,ispin)

      use global, only: maxlmindx,ncent,maxgpt

      use, intrinsic :: iso_fortran_env
      implicit none
      integer,     intent(in)  :: ikpt0,isym,iband,ispin
      complex_dp,  intent(out) :: cmtout(maxlmindx,*)
      MCOMPLEX_dp, intent(out) :: cpwout(maxgpt,*)

      call waveftrafo_mt(cmtout,maxlmindx,0,isym,iband,ikpt0,ispin)
      call waveftrafo_pw(cpwout,maxgpt,     isym,iband,ikpt0,ispin)

      end

c     ------------------

      recursive subroutine waveftrafo_mt(cmtout,dim,ic,isym,iband,ikpt0,ispin)

      use global

      use, intrinsic :: iso_fortran_env
      implicit none
      integer,    intent(in)  :: ic,ikpt0,isym,iband,ispin,dim
      complex_dp, intent(out) :: cmtout(dim,*)
      complex_dp, pointer_cnt :: cmtin(:,:)
      integer                 :: ncent0,icent,ic0,s
# include "interface/waveftrafo_mt_io.inc"

      if(l_soc.and.ispin==2) Bug('Spin index must not be 2 for SOC.')
      ic0 = ic
      if     (ic==0) then ; ncent0 = ncent                             ; icent = 1
      else if(ic> 0) then ; ncent0 = 1 ; ic0 = pcent(ic,sym(isym)%inv) ; icent = ic0
      else                ; ncent0 = neq(-ic)                          ; icent = 1 + sum(neq(:-ic-1))
      endif

      if(storeibz.and.ikpt0/=kptp(ikpt0)) then
        allocate(cmtin(dim,ncent0*nspin3))
        call wavefunction_mt(cmtin,dim,ic0,iband,ikpt0,ispin)
        call waveftrafo_mt_io(cmtout,cmtin,dim,dim,ic,isym,ikpt0)
        if(l_soc.or.(ispin==0.and.nspin==2)) call waveftrafo_mt_io(cmtout(1,ncent0+1),cmtin(:,ncent0+1),dim,dim,ic,isym,ikpt0)
        deallocate(cmtin)
      else
        s = max(1,ispin)
        cmtin => cmt(:,:,iband,kindx(ikpt0),s)
        call waveftrafo_mt_io(cmtout,cmtin(:,icent),dim,maxlmindx,ic,isym,ikpt0)
        if(l_soc.or.(ispin==0.and.nspin==2)) then
          cmtin => cmt(:,:,iband,kindx(ikpt0),2)
          call waveftrafo_mt_io(cmtout(1,ncent0+1),cmtin(:,icent),dim,maxlmindx,ic,isym,ikpt0)
        endif
        nullify(cmtin)
      endif

      NoInv( if(l_soc) call waveftrafo_soc(cmtout,dim*ncent0,isym) )

      end

c     ------------------

      recursive subroutine waveftrafo_pw(cpwout,dim,isym,iband,ikpt0,ispin)

      use global

      use, intrinsic :: iso_fortran_env
      implicit none
      integer,     intent(in)  :: ikpt0,isym,iband,ispin,dim
      MCOMPLEX_dp, intent(out) :: cpwout(dim,*)
      MCOMPLEX_dp, pointer_cnt :: cpwin(:)
      integer                  :: s
# include "interface/waveftrafo_pw_io.inc"

      if(l_soc.and.ispin==2) Bug('Spin index must not be 2 for SOC.')
      if(storeibz.and.ikpt0/=kptp(ikpt0)) then
        allocate(cpwin(dim*nspin3))
        call wavefunction_pw(cpwin,dim,iband,ikpt0,ispin)
        call waveftrafo_pw_io(cpwout,cpwin,isym,ikpt0)
        if(l_soc.or.(ispin==0.and.nspin==2)) call waveftrafo_pw_io(cpwout(1,2),cpwin(dim+1:),isym,ikpt0)
        deallocate(cpwin)
      else
        s = max(1,ispin)
        cpwin => cpw(:,iband,kindx(ikpt0),s)
        call waveftrafo_pw_io(cpwout,cpwin,isym,ikpt0)
        if(l_soc.or.(ispin==0.and.nspin==2)) then
          cpwin => cpw(:,iband,kindx(ikpt0),2)
          call waveftrafo_pw_io(cpwout(1,2),cpwin,isym,ikpt0)
        endif
        nullify(cpwin)
      endif

      NoInv( if(l_soc) call waveftrafo_soc(cpwout,dim,isym) )

      end

c     ------------------

c Get wavefunction (iband,ikpt,ispin) -> cmtout, cpwout
c
      subroutine wavefunction(cmtout,cpwout,iband,ikpt,ispin)
      use global, only: maxlmindx,maxgpt
      use, intrinsic :: iso_fortran_env
      implicit none
      integer, intent(in) :: iband,ikpt,ispin
      complex_dp          :: cmtout(maxlmindx,*)
      MCOMPLEX_dp         :: cpwout(maxgpt,*)

      call wavefunction_mt(cmtout,maxlmindx,0,iband,ikpt,ispin)
      call wavefunction_pw(cpwout,maxgpt,     iband,ikpt,ispin)
      
      end

c     ------------------

      subroutine wavefunction_mt(cmtout,dim,ic,iband,ikpt,ispin)

      use global

      use, intrinsic :: iso_fortran_env
      implicit none
      integer,    intent(in)  :: ic,ikpt,iband,ispin,dim
      complex_dp, intent(out) :: cmtout(dim,*)
      integer                 :: itype,ieq,icent,i,s,ic0

      if(l_soc.and.ispin==2) Bug('Spin index must not be 2 for SOC.')
      if(dim>maxlmindx) Bug('Leading dimension too large.')

      if(storeibz.and.kptp(ikpt)/=ikpt) then
        call waveftrafo_mt(cmtout,dim,ic,symkpt(ikpt),iband,kptp(ikpt),ispin)
      else
        i = 0
        do s = 1,nspin2 ; if(.not.l_soc.and.ispin/=0.and.ispin/=s) cycle
          icent = 0
          do itype = 1,ntype
            do ieq = 1,neq(itype)
              icent = icent + 1
              if(ic==0.or.icent==ic.or.ic==-itype) then
                i              = i + 1
                cmtout(:dim,i) = cmt(:dim,icent,iband,kindx(ikpt),s)
              endif
            enddo
          enddo
        enddo
      endif

      end

c     ------------------

      subroutine wavefunction_pw(cpwout,dim,iband,ikpt,ispin)

      use global

      use, intrinsic :: iso_fortran_env
      implicit none
      integer,     intent(in)  :: ikpt,iband,ispin,dim
      MCOMPLEX_dp, intent(out) :: cpwout(dim,*)

      if(l_soc.and.ispin==2) Bug('Spin index must not be 2 for SOC.')
      if(dim>maxgpt) Bug('Leading dimension too large.')

      if(storeibz.and.kptp(ikpt)/=ikpt) then
        call waveftrafo_pw(cpwout,dim,symkpt(ikpt),iband,kptp(ikpt),ispin)
      else
        if(ispin==0.or.l_soc) then ; cpwout(:dim,:nspin3) = cpw(:dim,iband,kindx(ikpt),:nspin3)
        else                       ; cpwout(:dim,1)       = cpw(:dim,iband,kindx(ikpt),ispin)
        endif
      endif

      end

c     ------------------

c
c Get wavefunction value at positions vec(:,:nvec) (iband,ikpt,ispin)
c (faster than wavefunction_r1 for one band or, in a loop, for few bands)
c
      subroutine wavefunction_r(wavef,vec,nvec,iband,ikpt,ispin)
      use global
      use, intrinsic :: iso_fortran_env
      implicit none
      integer,    intent(in)  :: nvec,iband,ikpt,ispin
      real_dp,    intent(in)  :: vec(3,nvec)
      complex_dp, intent(out) :: wavef(nvec)
      complex_dp              :: cmt1(maxlmindx,ncent),harm((maxlcut+1)**2)
      MCOMPLEX_dp             :: cpw1(maxgpt)
      real_dp                 :: r(3),rvec(3),rr,rad(maxindx),a(maxgrid-1),b(maxgrid-1),dy
      real                    :: kr,kgr,vec1(3)
      integer                 :: gpt1(3,ngpt(ikpt))
      integer                 :: t(3),ivec,l,m,n,lm,lmn,igpt
      integer                 :: itype,ieq,ic
      if(l_soc) Error('Not implemented: wavefunction_r for SOC.')      
      call wavefunction(cmt1,cpw1,iband,ikpt,ispin)
      gpt1  = gpt(:,pgpt(:ngpt(ikpt),ikpt))
      cpw1  = cpw1 / sqrt(vol)
      wavef = 0
      do ivec = 1,nvec
        t  = int(vec(:,ivec))
        r  = vec(:,ivec) - t
        ic = 0
        loop:
     &  do l = -1,1
        do m = -1,1
        do n = -1,1
          do itype = 1,ntype
            do ieq = 1,neq(itype)
              ic   = ic + 1
              rvec = matmul(lat,r-cent(:,ic)-[l,m,n])
              rr   = sum(rvec**2)
              if(rr<grid(itype)%radius**2) then
                t = t + [l,m,n]
                r = r - [l,m,n]
                exit loop
              endif
            enddo
          enddo
          ic = 0
        enddo
        enddo
        enddo loop
        if(ic>0) then ! in MT
          rr  = sqrt(rr)
          lm  = 0
          lmn = 0
          call harmonicsr(harm,rvec,lcut(itype))
          do l = 0,lcut(itype) ; m = max(2-l,1)
            do n = 1,nindx(l,itype)
              if(rr>=grid(itype)%radius) then
                rad(n) = bas1(grid(itype)%number,n,l,itype,ispin) / rr
              else if(rr<=rgrid(1,itype)) then
                rad(n) = bas1(1,n,l,itype,ispin) / rgrid(1,itype)
              else
                m = 1
                do while(rgrid(m,itype)<rr)
                  m = m + 1
                enddo
                rad(n) = ( ( rgrid(m,itype) - rr   ) * bas1(m-1,n,l,itype,ispin) +
     &                     ( rr - rgrid(m-1,itype) ) * bas1(m,  n,l,itype,ispin) ) / ( ( rgrid(m,itype) - rgrid(m-1,itype) ) * rr )
              endif
            enddo
            do m = -l,l
              lm = lm + 1
              do n = 1,nindx(l,itype)
                lmn         = lmn + 1
                wavef(ivec) = wavef(ivec) + rad(n) * harm(lm) * cmt1(lmn,ic)
              enddo
            enddo
          enddo
          wavef(ivec) = wavef(ivec) * exp( img * 2*pi * dot_product(kpt(:,ikpt),t) )
        else          ! in IR
          kr   = (2*pi) * dot_product( kpt(:,ikpt) , vec(:,ivec) )
          vec1 = vec(:,ivec)
          do igpt = 1,ngpt(ikpt)
            kgr         = kr + (2*pi) * dot_product( gpt1(:,igpt) , vec1 )
            wavef(ivec) = wavef(ivec) + cpw1(igpt) * cmplx( cos(kgr) , sin(kgr) )
          enddo
        endif
      enddo
      end

c     ------------------

c
c Get wavefunction value at positions vec(:,:nvec) (band1:band2,ikpt,ispin)
c (faster than wavefunction_r for many bands)
c
      subroutine wavefunction_r1(wavef,vec,nvec,band1,band2,ikpt,ispin)
      use global
      use, intrinsic :: iso_fortran_env
      implicit none
      integer,    intent(in)  :: nvec,band1,band2,ikpt,ispin
      real_dp,    intent(in)  :: vec(3,nvec)
      complex_dp, intent(out) :: wavef(band1:band2,nvec)
      complex_dp              :: cmt1(band1:band2,maxlmindx,ncent),harm((maxlcut+1)**2),cdum
      MCOMPLEX_dp             :: cpw1(band1:band2,maxgpt)
      real_dp                 :: r(3),rvec(3),rr,rad(maxindx),a(maxgrid-1),b(maxgrid-1),dy
      real                    :: kr,kgr,vec1(3)
      complex                 :: cexp
      integer                 :: gpt1(3,ngpt(ikpt))
      integer                 :: t(3),ivec,l,m,n,lm,lmn,igpt,iband
      integer                 :: itype,ieq,ic
      if(l_soc) Error('Not implemented: wavefunction_r1 for SOC.')
      do iband = band1,band2
        call wavefunction(cmt1(iband,:,:),cpw1(iband,:),iband,ikpt,ispin)
      enddo
      gpt1  = gpt(:,pgpt(:ngpt(ikpt),ikpt))
      cpw1  = cpw1 / sqrt(vol)
      wavef = 0
      do ivec = 1,nvec
        t  = int(vec(:,ivec))
        r  = vec(:,ivec) - t
        ic = 0
        loop:
     &  do l = -1,1
        do m = -1,1
        do n = -1,1
          do itype = 1,ntype
            do ieq = 1,neq(itype)
              ic   = ic + 1
              rvec = matmul(lat,r-cent(:,ic)-[l,m,n])
              rr   = sum(rvec**2)
              if(rr<grid(itype)%radius**2) then
                t = t + [l,m,n]
                r = r - [l,m,n]
                exit loop
              endif
            enddo
          enddo
          ic = 0
        enddo
        enddo
        enddo loop
        if(ic>0) then ! in MT
          rr  = sqrt(rr)
          lm  = 0
          lmn = 0
          call harmonicsr(harm,rvec,lcut(itype))
          do l = 0,lcut(itype) ; m = max(2-l,1)
            do n = 1,nindx(l,itype)
              if(rr>=grid(itype)%radius) then
                rad(n) = bas1(grid(itype)%number,n,l,itype,ispin) / rr
              else if(rr<=rgrid(1,itype)) then
                rad(n) = bas1(1,n,l,itype,ispin) / rgrid(1,itype)
              else
                m = 1
                do while(rgrid(m,itype)<rr)
                  m = m + 1
                enddo
                rad(n) = ( ( rgrid(m,itype) - rr   ) * bas1(m-1,n,l,itype,ispin) +
     &                     ( rr - rgrid(m-1,itype) ) * bas1(m,  n,l,itype,ispin) ) / ( ( rgrid(m,itype) - rgrid(m-1,itype) ) * rr )
              endif
            enddo
            do m = -l,l
              lm = lm + 1
              do n = 1,nindx(l,itype)
                lmn           = lmn + 1
                cdum          = rad(n) * harm(lm)
                wavef(:,ivec) = wavef(:,ivec) + cdum * cmt1(:,lmn,ic)
              enddo
            enddo
          enddo
          wavef(:,ivec) = wavef(:,ivec) * exp( img * 2*pi * dot_product(kpt(:,ikpt),t) )
        else          ! in IR
          kr   = (2*pi) * dot_product( kpt(:,ikpt) , vec(:,ivec) )
          vec1 = vec(:,ivec)
          do igpt = 1,ngpt(ikpt)
            kgr           = kr + (2*pi) * dot_product( gpt1(:,igpt) , vec1 )
            cexp          = cmplx( cos(kgr) , sin(kgr) )
            wavef(:,ivec) = wavef(:,ivec) + cpw1(:,igpt) * cexp
          enddo
        endif
      enddo
      end

# else

c --------------

      recursive subroutine waveftrafo(cmtout,cpwout,ikpt0,isym,iband,ispin)

      use global

      use, intrinsic :: iso_fortran_env
      implicit none
      integer,    intent(in)  :: ikpt0,isym,iband,ispin
      complex_dp, intent(out) :: cmtout(maxlmindx,ncent)
      MCOMPLEX_dp,intent(out) :: cpwout(ngpt(ikpt0))
      complex_dp              :: cmtin(maxlmindx,ncent)
      MCOMPLEX_dp             :: cpwin(ngpt(ikpt0))
      complex_dp              :: cdum,cmthlp(2*maxlcut+1) InvC(cphase)
      logical                 :: trs
      integer                 :: g(3),g1(3),pgptinv(ngptall)
      integer                 :: l,n,nn,lm0,lm1,lm2
      integer                 :: itype,ieq,icent,icent1,igpt,igpt0,igpt1,ikpt1,isymi

      if(storeibz.and.kptp(ikpt0)/=ikpt0) then
        if(l_soc) Error('recursive call to wavetrafo not implemented for SOC.')
        call waveftrafo(cmtin,cpwin,kptp(ikpt0),symkpt(ikpt0),iband,ispin)
      else
        cmtin = cmt(:,:,         iband,kindx(ikpt0),ispin)
        cpwin = cpw(:ngpt(ikpt0),iband,kindx(ikpt0),ispin)
      endif
      ikpt1 = kptsym(ikpt0,isym)
      isymi = sym(isym)%inv
      trs   = isym>nsymt
      Inv( cphase = exp(-img * 2*pi * dot_product(kpt(:,ikpt1),sym(isym)%transl) ) )

      ! MT coefficients
      cmtout = 0
      icent  = 0
      do itype = 1,ntype
        do ieq = 1,neq(itype)
          icent  = icent + 1
          icent1 = pcent(icent,isym)
          cdum   = exp( -img * (2*pi) * dot_product(kpt(:,ikpt1),tcent(:,icent,isym)) ) Inv( / cphase )
          lm0 = 0
          do l = 0,lcut(itype)
            nn = nindx(l,itype)
            do n = 1,nn
              lm1            = lm0 + n
              lm2            = lm0 + n + 2*l*nn
              if(trs) then ; cmthlp(:2*l+1) = conjg ( cmtin(lm1:lm2:nn,icent) )
              else         ; cmthlp(:2*l+1) =         cmtin(lm1:lm2:nn,icent)
              endif
              cmtout(lm1:lm2:nn,icent1) = cdum * matmul(dwgn(-l:l,-l:l,l,isym),cmthlp(:2*l+1))
            enddo
            lm0 = lm2
          enddo
        enddo
      enddo

      ! PW coefficients
      pgptinv                           = 0
      pgptinv(pgpt(:ngpt(ikpt0),ikpt0)) = [ (igpt,igpt=1,ngpt(ikpt0)) ]
      cpwout = 0
      g      = nint ( matmul(sym(isym)%rrot,kpt(:,ikpt0))-kpt(:,ikpt1) ) ! ROT*kp = k + g,  g = rec. latt. vector
      do igpt1 = 1,ngpt(ikpt1)
        igpt  = pgpt(igpt1,ikpt1)
        if(sum(matmul(rlat,kpt(:,ikpt1)+gpt(:,igpt))**2)>gcut**2) Bug('|k+G|>gcut.')
        g1    = matmul( sym(isymi)%rrot,gpt(:,igpt)-g )
        igpt0 = pgptinv ( pntgpt(g1(1),g1(2),g1(3)) ) ; if(igpt0==0) Bug('zero pointer')
# ifdef INV
        cpwout(igpt1) = cpwin(igpt0)
# else
        cdum  = exp(-img * 2*pi * dot_product (kpt(:,ikpt1)+gpt(:,igpt), sym(isym)%transl))
        if(trs) then ; cpwout(igpt1) = cdum * MCONJG ( cpwin(igpt0) )
        else         ; cpwout(igpt1) = cdum *          cpwin(igpt0)
        endif
# endif
      enddo

# ifdef INV
      if(any(sym(isym)%transl/=0)) then ! factor exp(-iGa)=+-1 [a=sym(i)%transl]
        g = nint(2*sym(isym)%transl)    ! 2a must be integer
        do igpt1 = 1,ngpt(ikpt1)
          igpt = pgpt(igpt1,ikpt1)
          if(mod(dot_product(gpt(:,igpt),g),2)/=0) cpwout(igpt1) = -cpwout(igpt1)
        enddo
      endif
# endif

      end

c --------------

c Input/Output version: Transforms cmtio,cpwio from kptp(ikpt) to ikpt (mode=1:mt,2:pw,3:mt&pw)
      subroutine waveftrafo_io(cmtio,cpwio,ikpt,mode)

      use global

      use, intrinsic :: iso_fortran_env
      implicit none
      complex_dp, intent(inout) :: cmtio(maxlmindx,ncent,*)
      MCOMPLEX_dp,intent(inout) :: cpwio(maxgpt,*)
      integer,    intent(in)    :: ikpt,mode
      complex_dp                :: cmtin(maxlmindx,ncent)
      MCOMPLEX_dp               :: cpwin(maxgpt)
      complex_dp                :: cdum,cmthlp(2*maxlcut+1)
      logical                   :: trs
      integer                   :: g(3),g1(3),pgptinv(ngptall),isym,s
      integer                   :: l,n,nn,lm0,lm1,lm2
      integer                   :: itype,ieq,icent,icent1,igpt,igpt0,igpt1,isymi,ikpt0

      if(mode<1.or.mode>3) Bug('Wrong mode.')

      isym  = symkpt(ikpt) ; if(isym==1) return
      isymi = sym(isym)%inv
      ikpt0 = kptp(ikpt)
      trs   = isym>nsymt

      do s = 1,nspin3

      ! MT coefficients
      if(iand(mode,1)/=0) then
        cmtin        = cmtio(:,:,s)
        cmtio(:,:,s) = 0
        icent        = 0
        do itype = 1,ntype
          do ieq = 1,neq(itype)
            icent  = icent + 1
            icent1 = pcent(icent,isym)
            cdum   = exp( -img * (2*pi) * dot_product(kpt(:,ikpt),tcent(:,icent,isym)) ) Inv( / phase(ikpt) )
            lm0 = 0
            do l = 0,lcut(itype)
              nn = nindx(l,itype)
              do n = 1,nn
                lm1            = lm0 + n
                lm2            = lm0 + n + 2*l*nn
                if(trs) then ; cmthlp(:2*l+1) = conjg ( cmtin(lm1:lm2:nn,icent) )
                else         ; cmthlp(:2*l+1) =         cmtin(lm1:lm2:nn,icent)
                endif
                cmtio(lm1:lm2:nn,icent1,s) = cdum * matmul(dwgn(-l:l,-l:l,l,isym),cmthlp(:2*l+1))
              enddo
              lm0 = lm2
            enddo
          enddo
        enddo
      endif

      ! PW coefficients
      if(iand(mode,2)/=0) then
        cpwin                             = cpwio(:,s)
        pgptinv                           = 0
        pgptinv(pgpt(:ngpt(ikpt0),ikpt0)) = [ (igpt,igpt=1,ngpt(ikpt0)) ]
        cpwio(:,s) = 0
        g          = nint ( matmul(sym(isym)%rrot,kpt(:,ikpt0))-kpt(:,ikpt) ) ! ROT*kp = k + g,  g = rec. latt. vector
        do igpt1 = 1,ngpt(ikpt)
          igpt  = pgpt(igpt1,ikpt)
          if(sum(matmul(rlat,kpt(:,ikpt)+gpt(:,igpt))**2)>gcut**2) Bug('|k+G|>gcut.')
          g1    = matmul( sym(isymi)%rrot,gpt(:,igpt)-g )
          igpt0 = pgptinv ( pntgpt(g1(1),g1(2),g1(3)) ) ; if(igpt0==0) Bug('zero pointer')
# ifdef INV
          cpwio(igpt1,s) = cpwin(igpt0)
# else
          cdum  = exp(-img * 2*pi * dot_product (kpt(:,ikpt)+gpt(:,igpt), sym(isym)%transl))
          if(trs) then ; cpwio(igpt1,s) = cdum * MCONJG ( cpwin(igpt0) )
          else         ; cpwio(igpt1,s) = cdum *          cpwin(igpt0)
          endif
# endif
        enddo

# ifdef INV
        if(any(sym(isym)%transl/=0)) then ! factor exp(-iGa)=+-1 [a=sym(i)%transl]
          g = nint(2*sym(isym)%transl)    ! 2a must be integer
          do igpt1 = 1,ngpt(ikpt)
            igpt = pgpt(igpt1,ikpt)
            if(mod(dot_product(gpt(:,igpt),g),2)/=0) cpwio(igpt1,s) = -cpwio(igpt1,s)
          enddo
        endif
# endif
      endif

      enddo ! spin loop

# ifndef INV
      if(l_soc) then
        if(iand(mode,1)/=0) call waveftrafo_soc(cmtio(:,:,:nspin3),maxlmindx*ncent,symkpt(ikpt))
        if(iand(mode,2)/=0) call waveftrafo_soc(cpwio(:,:nspin3),  maxgpt,         symkpt(ikpt))
      endif
# endif

      end

c --------------

c     Only for the muffin-tin part (rotates to ikpt1)
      subroutine waveftrafo_mt(cmtout,maxlm,ikpt1,iband,ispin,ic,lwan)

      use global
      use, intrinsic :: iso_fortran_env

      implicit none
      integer,    intent(in)  :: maxlm
      complex_dp, intent(out) :: cmtout(maxlm)
      integer,    intent(in)  :: ikpt1,iband,ispin,ic
      logical,    intent(in)  :: lwan
      complex_dp              :: cmtin(maxlm)
      complex_dp              :: cdum,cmthlp(2*maxlcut+1)
      logical                 :: trs
      integer                 :: isym,ikindx
      integer                 :: l,n,nn,lm0,lm1,lm2,lcut_
      integer                 :: itype,ic1,isymi

c      if(.not.storeibz)      Bug('storeibz=.false.')
      if(kptp(ikpt1)==ikpt1) Bug('k point in IBZ.')
      if(maxlm<=0)           Bug('Dimension maxlm must be positive.')

      ikindx = kindx(ikpt1)
      isym   = symkpt(ikpt1)
      isymi  = sym(isym)%inv
      ic1    = pcent(ic,isymi)
      trs    = isym>nsymt
      if(lwan) then ; cmtin = cmtu(:maxlm,ic1,iband,ikindx,ispin)
      else          ; cmtin = cmt (:maxlm,ic1,iband,ikindx,ispin)
      endif

      ! Determine itype from ic
      do itype = 1,ntype
        if(sum(neq(:itype))>=ic) exit
      enddo
      ! Determine lcut_ from maxlm
      lm1 = 0
      do l = 0,lcut(itype)
        lm1   = lm1 + (2*l+1)*nindx(l,itype) ; if(lm1>maxlm) exit
        lcut_ = l
      enddo
      ! MT coefficients
      cmtout = 0
      cdum   = exp( -img * (2*pi) * dot_product(kpt(:,ikpt1),tcent(:,ic1,isym)) ) ; Inv( if(.not.lwan) cdum = cdum / phase(ikpt1) )
      lm0    = 0
      do l = 0,lcut_
        nn = nindx(l,itype)
        do n = 1,nn
          lm1            = lm0 + n
          lm2            = lm0 + n + 2*l*nn
          if(trs) then ; cmthlp(:2*l+1) = conjg ( cmtin(lm1:lm2:nn) )
          else         ; cmthlp(:2*l+1) =         cmtin(lm1:lm2:nn)
          endif
          cmtout(lm1:lm2:nn) = cdum * matmul(dwgn(-l:l,-l:l,l,isym),cmthlp(:2*l+1))
        enddo
        lm0 = lm2
      enddo

      end

c --------------

c     Only for the plane-wave part (IBZ->BZ)
      subroutine waveftrafo_pw(cpwout,ikpt1,iband,ispin,lwan)

      use global
      use, intrinsic :: iso_fortran_env

      implicit none
      integer,     intent(in)  :: ikpt1,iband,ispin
      logical,     intent(in)  :: lwan
      MCOMPLEX_dp, intent(out) :: cpwout(ngpt(ikpt1))
      MCOMPLEX_dp              :: cpwin(ngpt(ikpt1))
      integer                  :: g(3),g1(3),pgptinv(ngptall)
      integer                  :: ikpt0,igpt,igpt0,igpt1,isym,isymi,ikindx
      NoInv( logical           :: trs  )
      NoInv( complex_dp        :: cdum )

c      if(.not.storeibz) Bug('storeibz=.false.')

# ifdef INV
      if(lwan) then
        call waveftrafo_pw_w(cpwout,ikpt1,iband,ispin)
        return
      endif
# endif
      ikpt0  = kptp(ikpt1)
      ikindx = kindx(ikpt1)
      isym   = symkpt(ikpt1)
      isymi  = sym(isym)%inv
      NoInv( trs = isym>nsymt )
      if(lwan) then ; cpwin = cpwu(:ngpt(ikpt1),iband,ikindx,ispin)
      else          ; cpwin = cpw (:ngpt(ikpt1),iband,ikindx,ispin)
      endif

      ! PW coefficients
      pgptinv                           = 0
      pgptinv(pgpt(:ngpt(ikpt0),ikpt0)) = [ (igpt,igpt=1,ngpt(ikpt0)) ]
      cpwout = 0
      g      = nint ( matmul(sym(isym)%rrot,kpt(:,ikpt0))-kpt(:,ikpt1) ) ! ROT*kp = k + g,  g = rec. latt. vector
      do igpt1 = 1,ngpt(ikpt1)
        igpt  = pgpt(igpt1,ikpt1)
        if(sum(matmul(rlat,kpt(:,ikpt1)+gpt(:,igpt))**2)>gcut**2) Bug('|k+G|>gcut.')
        g1    = matmul( sym(isymi)%rrot,gpt(:,igpt)-g )
        igpt0 = pgptinv ( pntgpt(g1(1),g1(2),g1(3)) ) ; if(igpt0==0) Bug('zero pointer')
# ifdef INV
        cpwout(igpt1) = cpwin(igpt0)
# else
        cdum  = exp(-img * 2*pi * dot_product (kpt(:,ikpt1)+gpt(:,igpt), sym(isym)%transl))
        if(trs) then ; cpwout(igpt1) = cdum * MCONJG ( cpwin(igpt0) )
        else         ; cpwout(igpt1) = cdum *          cpwin(igpt0)
        endif
# endif
      enddo

# ifdef INV
      if(any(sym(isym)%transl/=0)) then ! factor exp(-iGa)=+-1 [a=sym(i)%transl]
        g = nint(2*sym(isym)%transl)    ! 2a must be integer
        do igpt1 = 1,ngpt(ikpt1)
          igpt = pgpt(igpt1,ikpt1)
          if(mod(dot_product(gpt(:,igpt),g),2)/=0) cpwout(igpt1) = -cpwout(igpt1)
        enddo
      endif
# endif

      end

c --------------

c Carbon copy of above for Wannier and inversion (arrays have to be complex)
# ifdef INV
      subroutine waveftrafo_pw_w(cpwout,ikpt1,iband,ispin)
      use global
      use, intrinsic :: iso_fortran_env
      implicit none
      integer,     intent(in)  :: ikpt1,iband,ispin
      complex_dp,  intent(out) :: cpwout(ngpt(ikpt1))
      complex_dp               :: cdum
      integer                  :: g(3),g1(3),pgptinv(ngptall)
      integer                  :: igpt,igpt0,igpt1,isym,isymi,ikindx,ikpt0
c      Error('With the new phases it is not clear how to change waveftrafo_pw_w')
      ikpt0  = kptp(ikpt1)
      ikindx = kindx(ikpt1)
      isym   = symkpt(ikpt1)
      isymi  = sym(isym)%inv
      ! PW coefficients
      pgptinv                           = 0
      pgptinv(pgpt(:ngpt(ikpt0),ikpt0)) = [ (igpt,igpt=1,ngpt(ikpt0)) ]
      cpwout = 0
      g      = nint ( matmul(sym(isym)%rrot,kpt(:,ikpt0))-kpt(:,ikpt1) ) ! ROT*kp = k + g,  g = rec. latt. vector
      do igpt1 = 1,ngpt(ikpt1)
        igpt  = pgpt(igpt1,ikpt1)
        if(sum(matmul(rlat,kpt(:,ikpt1)+gpt(:,igpt))**2)>gcut**2) Bug('|k+G|>gcut.')
        g1    = matmul( sym(isymi)%rrot,gpt(:,igpt)-g )
        igpt0 = pgptinv ( pntgpt(g1(1),g1(2),g1(3)) ) ; if(igpt0==0) Bug('zero pointer')
        cdum  = exp(-img * 2*pi * dot_product (kpt(:,ikpt1)+gpt(:,igpt), sym(isym)%transl))
        cpwout(igpt1) = cdum * cpwu(igpt0,iband,ikindx,ispin)
      enddo
      end
# endif

c --------------

c Same as waveftrafo for all bands (only called by getinput).
      subroutine waveftrafo2(cmtout,cpwout,ikpt0,isym,ispin)

      use global
      use, intrinsic :: iso_fortran_env

      implicit none
      complex_dp, intent(out) :: cmtout(maxlmindx,ncent,maxband)
      MCOMPLEX_dp,intent(out) :: cpwout(maxgpt,maxband)
      integer,    intent(in)  :: ikpt0,isym,ispin
      complex_dp              :: cdum,cmthlp(2*maxlcut+1) InvC(cphase)
      logical                 :: trs
      integer                 :: g(3),g1(3),pgptinv(ngptall) InvC(idum)
      integer                 :: l,n,nn,lm0,lm1,lm2,itype,ieq,icent,icent1,igpt,igpt0,igpt1,iband,ikpt1,isymi,ikindx

      ikindx = kindx(kptp(ikpt0))
      ikpt1  = kptsym(ikpt0,isym)
      isymi  = sym(isym)%inv
      trs    = isym>nsymt
      Inv( cphase = exp(-img * 2*pi * dot_product(kpt(:,ikpt1),sym(isym)%transl) ) )

      ! MT coefficients
      cmtout = 0
      icent  = 0
      do itype = 1,ntype
        do ieq = 1,neq(itype)
          icent  = icent + 1
          icent1 = pcent(icent,isym)
          cdum   = exp( -img * (2*pi) * dot_product(kpt(:,ikpt1),tcent(:,icent,isym)) ) Inv( / cphase )
          lm0 = 0
          do l = 0,lcut(itype)
            nn = nindx(l,itype)
            do n = 1,nn
              lm1 = lm0 + n
              lm2 = lm0 + n + 2*l*nn
              do iband = 1,nband(ikpt1,min(ispin,nspin))
                if(trs) then ; cmthlp(:2*l+1) = conjg ( cmt(lm1:lm2:nn,icent,iband,ikindx,ispin) )
                else         ; cmthlp(:2*l+1) =         cmt(lm1:lm2:nn,icent,iband,ikindx,ispin)
                endif
                cmtout(lm1:lm2:nn,icent1,iband) = cdum * matmul(dwgn(-l:l,-l:l,l,isym),cmthlp(:2*l+1))
              enddo
            enddo
            lm0 = lm2
          enddo
        enddo
      enddo

      ! PW coefficients
      pgptinv                           = 0
      pgptinv(pgpt(:ngpt(ikpt0),ikpt0)) = [ (igpt,igpt=1,ngpt(ikpt0)) ]
      cpwout = 0
      g      = nint ( matmul(sym(isym)%rrot,kpt(:,ikpt0))-kpt(:,ikpt1) ) ! ROT*kp = k + g,  g = rec. latt. vector
      do igpt1 = 1,ngpt(ikpt1)
        igpt  = pgpt(igpt1,ikpt1)
        if(sum(matmul(rlat,kpt(:,ikpt1)+gpt(:,igpt))**2)>gcut**2) Bug('|k+G|>gcut.')
        g1    = matmul( sym(isymi)%rrot,gpt(:,igpt)-g )
        igpt0 = pgptinv ( pntgpt(g1(1),g1(2),g1(3)) ) ; if(igpt0==0) Bug('zero pointer')
# ifdef INV
        cdum  = exp(-img * 2*pi * dot_product (gpt(:,igpt), sym(isym)%transl))
        idum  = nint(real(cdum))
        cpwout(igpt1,:) = idum * cpw(igpt0,:,ikindx,ispin)
# else
        cdum  = exp(-img * 2*pi * dot_product (kpt(:,ikpt1)+gpt(:,igpt), sym(isym)%transl))
        if(trs) then ; cpwout(igpt1,:) = cdum * MCONJG ( cpw(igpt0,:,ikindx,ispin) )
        else         ; cpwout(igpt1,:) = cdum *          cpw(igpt0,:,ikindx,ispin)
        endif
# endif
      enddo

      end

c     ------------------

      subroutine waveftrafo_soc(c,n,isym)
      use global, only: sym
      use, intrinsic :: iso_fortran_env
      implicit none
      integer,    intent(in)    :: n,isym
      complex_dp, intent(inout) :: c(n,2)
      integer                   :: i
      do i = 1,n
        c(i,:) = matmul ( c(i,:) , sym(isym)%esoc )
      enddo
      end

      subroutine waveftrafo1_soc(c,n,esoc)
      use, intrinsic :: iso_fortran_env
      implicit none
      integer,    intent(in)    :: n
      complex_dp, intent(in)    :: esoc(2,2)
      complex_dp, intent(inout) :: c(n,2)
      integer                   :: i
      do i = 1,n
        c(i,:) = matmul ( c(i,:) , esoc )
      enddo
      end

c     ------------------

      subroutine wavefunction(cmtout,cpwout,iband,ikpt,ispin)
      use global, only: storeibz,kptp,symkpt,cmt,cpw,maxlmindx,ncent,ngpt,l_soc,kindx
      use, intrinsic :: iso_fortran_env
      implicit none
      integer, intent(in) :: iband,ikpt,ispin
      complex_dp          :: cmtout(maxlmindx,ncent)
      MCOMPLEX_dp         :: cpwout(ngpt(ikpt))
      if(storeibz.and.kptp(ikpt)/=ikpt) then
        call waveftrafo(cmtout,cpwout,kptp(ikpt),symkpt(ikpt),iband,ispin)
        if(l_soc) Error('Not implemented for SOC. Inform the developers.')        
      else
        cmtout = cmt(:,:,        iband,kindx(ikpt),ispin)
        cpwout = cpw(:ngpt(ikpt),iband,kindx(ikpt),ispin)
      endif
      end

c     ------------------

      subroutine wavefunction_r(wavef,vec,nvec,iband,ikpt,ispin)
      use global
      use, intrinsic :: iso_fortran_env
      implicit none
      integer,    intent(in)  :: nvec,iband,ikpt,ispin
      real_dp,    intent(in)  :: vec(3,nvec)
      complex_dp, intent(out) :: wavef(nvec)
      complex_dp              :: cmt1(maxlmindx,ncent),harm((maxlcut+1)**2)
      MCOMPLEX_dp             :: cpw1(ngpt(ikpt))
      real_dp                 :: r(3),rvec(3),rr,rad(maxindx),a(maxgrid-1),b(maxgrid-1),dy
      real                    :: kr,kgr,vec1(3)
      integer                 :: gpt1(3,ngpt(ikpt))
      integer                 :: t(3),ivec,l,m,n,lm,lmn,igpt
      integer                 :: itype,ieq,ic
      call wavefunction(cmt1,cpw1,iband,ikpt,ispin)
      gpt1  = gpt(:,pgpt(:ngpt(ikpt),ikpt))
      cpw1  = cpw1 / sqrt(vol)
      wavef = 0
      do ivec = 1,nvec
        t  = int(vec(:,ivec))
        r  = vec(:,ivec) - t
        ic = 0
        loop:
     &  do l = -1,1
        do m = -1,1
        do n = -1,1
          do itype = 1,ntype
            do ieq = 1,neq(itype)
              ic   = ic + 1
              rvec = matmul(lat,r-cent(:,ic)-[l,m,n])
              rr   = sum(rvec**2)
              if(rr<grid(itype)%radius**2) then
                t = t + [l,m,n]
                r = r - [l,m,n]
                exit loop
              endif
            enddo
          enddo
          ic = 0
        enddo
        enddo
        enddo loop
        if(ic>0) then ! in MT
          rr  = sqrt(rr)
          lm  = 0
          lmn = 0
          call harmonicsr(harm,rvec,lcut(itype))
          do l = 0,lcut(itype) ; m = max(2-l,1)
            do n = 1,nindx(l,itype)
              if(rr>=grid(itype)%radius) then
                rad(n) = bas1(grid(itype)%number,n,l,itype,ispin) / rr
              else if(rr<=rgrid(1,itype)) then
                rad(n) = bas1(1,n,l,itype,ispin) / rgrid(1,itype)
              else
                m = 1
                do while(rgrid(m,itype)<rr)
                  m = m + 1
                enddo
                rad(n) = ( ( rgrid(m,itype) - rr   ) * bas1(m-1,n,l,itype,ispin) +
     &                     ( rr - rgrid(m-1,itype) ) * bas1(m,  n,l,itype,ispin) ) / ( ( rgrid(m,itype) - rgrid(m-1,itype) ) * rr )
              endif
            enddo
            do m = -l,l
              lm = lm + 1
              do n = 1,nindx(l,itype)
                lmn         = lmn + 1
                wavef(ivec) = wavef(ivec) + rad(n) * harm(lm) * cmt1(lmn,ic)
              enddo
            enddo
          enddo
          wavef(ivec) = wavef(ivec) * exp( img * 2*pi * dot_product(kpt(:,ikpt),t) )
        else          ! in IR
          kr   = (2*pi) * dot_product( kpt(:,ikpt) , vec(:,ivec) )
          vec1 = vec(:,ivec)
          do igpt = 1,ngpt(ikpt)
            kgr         = kr + (2*pi) * dot_product( gpt1(:,igpt) , vec1 )
            wavef(ivec) = wavef(ivec) + cpw1(igpt) * cmplx( cos(kgr) , sin(kgr) )
          enddo
        endif
      enddo
      end

      subroutine wavefunction_r1(wavef,vec,nvec,band1,band2,ikpt,ispin)
      use global
      use, intrinsic :: iso_fortran_env
      implicit none
      integer,    intent(in)  :: nvec,band1,band2,ikpt,ispin
      real_dp,    intent(in)  :: vec(3,nvec)
      complex_dp, intent(out) :: wavef(band1:band2,nvec)
      complex_dp              :: cmt1(band1:band2,maxlmindx,ncent),harm((maxlcut+1)**2),cdum
      MCOMPLEX_dp             :: cpw1(band1:band2,ngpt(ikpt))
      real_dp                 :: r(3),rvec(3),rr,rad(maxindx),a(maxgrid-1),b(maxgrid-1),dy
      real                    :: kr,kgr,vec1(3)
      complex                 :: cexp
      integer                 :: gpt1(3,ngpt(ikpt))
      integer                 :: t(3),ivec,l,m,n,lm,lmn,igpt,iband
      integer                 :: itype,ieq,ic
      do iband = band1,band2
        call wavefunction(cmt1(iband,:,:),cpw1(iband,:),iband,ikpt,ispin)
      enddo
      gpt1  = gpt(:,pgpt(:ngpt(ikpt),ikpt))
      cpw1  = cpw1 / sqrt(vol)
      wavef = 0
      do ivec = 1,nvec
        t  = int(vec(:,ivec))
        r  = vec(:,ivec) - t
        ic = 0
        loop:
     &  do l = -1,1
        do m = -1,1
        do n = -1,1
          do itype = 1,ntype
            do ieq = 1,neq(itype)
              ic   = ic + 1
              rvec = matmul(lat,r-cent(:,ic)-[l,m,n])
              rr   = sum(rvec**2)
              if(rr<grid(itype)%radius**2) then
                t = t + [l,m,n]
                r = r - [l,m,n]
                exit loop
              endif
            enddo
          enddo
          ic = 0
        enddo
        enddo
        enddo loop
        if(ic>0) then ! in MT
          rr  = sqrt(rr)
          lm  = 0
          lmn = 0
          call harmonicsr(harm,rvec,lcut(itype))
          do l = 0,lcut(itype) ; m = max(2-l,1)
            do n = 1,nindx(l,itype)
              if(rr>=grid(itype)%radius) then
                rad(n) = bas1(grid(itype)%number,n,l,itype,ispin) / rr
              else if(rr<=rgrid(1,itype)) then
                rad(n) = bas1(1,n,l,itype,ispin) / rgrid(1,itype)
              else
                m = 1
                do while(rgrid(m,itype)<rr)
                  m = m + 1
                enddo
                rad(n) = ( ( rgrid(m,itype) - rr   ) * bas1(m-1,n,l,itype,ispin) +
     &                     ( rr - rgrid(m-1,itype) ) * bas1(m,  n,l,itype,ispin) ) / ( ( rgrid(m,itype) - rgrid(m-1,itype) ) * rr )
              endif
            enddo
            do m = -l,l
              lm = lm + 1
              do n = 1,nindx(l,itype)
                lmn           = lmn + 1
                cdum          = rad(n) * harm(lm)
                wavef(:,ivec) = wavef(:,ivec) + cdum * cmt1(:,lmn,ic)
              enddo
            enddo
          enddo
          wavef(:,ivec) = wavef(:,ivec) * exp( img * 2*pi * dot_product(kpt(:,ikpt),t) )
        else          ! in IR
          kr   = (2*pi) * dot_product( kpt(:,ikpt) , vec(:,ivec) )
          vec1 = vec(:,ivec)
          do igpt = 1,ngpt(ikpt)
            kgr           = kr + (2*pi) * dot_product( gpt1(:,igpt) , vec1 )
            cexp          = cmplx( cos(kgr) , sin(kgr) )
            wavef(:,ivec) = wavef(:,ivec) + cpw1(:,igpt) * cexp
          enddo
        endif
      enddo
      end

c     ------------------

      subroutine wavefunction_mt
      Error('wavefunction_mt not available for -Dold_trafo.')
      end

# endif

c     ------------------

      ! Symmetrizes MT part of input matrix according to inversion symmetry.
      ! This is achieved by a transformation to
      !
      !        1/sqrt(2) * ( exp(ikR) Y_lm(r-R) + (-1)**(l+m) exp(-ikR) Y_l,-m(r+R) )
      ! and                                                                               if R/=0 or m<0
      !        i/sqrt(2) * ( exp(ikR) Y_lm(r-R) - (-1)**(l+m) exp(-ikR) Y_l,-m(r+R) )
      ! or
      !        i * Y_l,0(r)                                                               if R=0, m=0, and l odd .
      !
      ! These functions have the property f(-r)=conjg(f(r)) which makes the output matrix real symmetric.
      ! Note that the transformation remains valid if tcent(ic,invsym)/=0 (i.e., pcent(ic,invsym) is not at -R).
      ! (Array MAT is overwritten.)
      !
      ! imode0 = sgn * ( imode + 3*typ + 3*(ntype+1)*indx + 3*(ntype+1)*symm_maxindx*dim )
      !
      ! imode = 1       : symmetrize bra of matrix
      ! imode = 2       : symmetrize ket of matrix
      ! imode = 3       : symmetrize bra and ket of matrix
      ! sgn   = 1       : bra = rows    , ket = columns (normal order)
      ! sgn   = -1      : bra = columns , ket = rows    (order flipped)
      ! typ   > 0       : symmetrize only part of atom type typ.
      ! typ   = 0       : symmetrize everything.
      ! indx  = 0       : employ product basis (nindxm, lcutm)
      ! indx  = 1       : employ LAPW basis (nindx, lcut)
      ! indx  = 2       : employ LO part of LAPW basis (3:nindx, lcut)
      ! indx  > 2       : employ LAPW basis and set nindx = indx-2 for all atom types (indx<symm_maxindx)
      ! dim             : "dimension" of (lmn)-indx [for example, if mat has the form mat(:dim1,:dim,:ncent)]
      ! symm_maxindx is declared in "global.f" and set to a default value. (Must be increased to the maximal "indx" + 1 if needed.)
      !
      ! The matrix "mat" is assumed to have the dimensions mat(:dim1,:dim2) unless dim2=-1.
      ! Then, "mat" is assumed atom-diagonal, e.g., with the dimensions mat(:dim1,:,:ncent)
      ! (see subroutine "symm_desymm_atomdiag").
      !
      subroutine symmetrize(mat,dim1,dim2,imode0)
      use global
      use, intrinsic :: iso_fortran_env
      implicit none
      integer,    intent(in)    :: imode0,dim1,dim2
      complex_dp, intent(inout) :: mat(dim1,*)
      complex_dp                :: carr(max(dim1,dim2)),cfac,img1
      real_dp                   :: rfac
      integer                   :: imode,typ,indx,dim,sgn
      integer                   :: lcut1,nindx1(0:max(maxlcut,maxlcutm))
      integer                   :: i,j,i0,itype,ieq,ic,ic1,l,m,n,nn,ifac,ishift
      if(dim2==-1) then
        call symm_desymm_atomdiag(mat,dim1,imode0,1)
        return
      endif
      i = 0
      call imode_def(imode,typ,indx,dim,sgn,imode0,dim1,dim2,i)
      if(imode<3.and.sgn==-1) imode = 3 - imode
      img1 = sgn * img
      rfac = sqrt(0.5d0)
      cfac = sqrt(0.5d0) * img1
      ic   = 0
      i0   = 0
      do itype = 1,ntype
        if     (indx==0) then ; lcut1 = lcutm(itype) ; nindx1(:lcut1) = nindxm(:lcut1,itype)
        else if(indx==1) then ; lcut1 = lcut(itype)  ; nindx1(:lcut1) = nindx(:lcut1,itype)
        else if(indx==2) then ; lcut1 = lcut(itype)  ; nindx1(:lcut1) = nindx(:lcut1,itype) - 2
        else                  ; lcut1 = lcut(itype)  ; nindx1(:lcut1) = indx - 2
        endif
        if(dim>0) then ; nn = dim
        else           ; nn = sum( [ ((2*l+1)*nindx1(l),l=0,lcut1) ] )
        endif
        do ieq = 1,neq(itype)
          ic  = ic + 1 ; if(typ>0.and.itype/=typ) cycle
          ic1 = pcent(ic,invsym)
          if(ic1>=ic) then
            i = i0
            do l = 0,lcut1
              ifac  = -1
              do m = -l,l
                ifac   = -ifac
                ishift = (ic1-ic)*nn - 2*m*nindx1(l)
                do n = 1,nindx1(l)
                  i = i + 1
                  j = i + ishift
                  if(ic1/=ic.or.m<0) then
                    if(iand(imode,1)/=0) then
                      carr(:dim2)  = mat(i,:dim2)
                      mat(i,:dim2) = ( carr(:dim2) + ifac * mat(j,:dim2) ) * rfac
                      mat(j,:dim2) = ( carr(:dim2) - ifac * mat(j,:dim2) ) * (-cfac)
                    endif
                    if(iand(imode,2)/=0) then
                      carr(:dim1)  = mat(:,i)
                      mat(:,i)     = ( carr(:dim1) + ifac * mat(:,j) ) * rfac
                      mat(:,j)     = ( carr(:dim1) - ifac * mat(:,j) ) * cfac
                    endif
                  else if(m==0.and.ifac==-1) then
                    if(iand(imode,1)/=0) then
                      mat(i,:dim2) = -img1 * mat(i,:dim2)
                    endif
                    if(iand(imode,2)/=0) then
                      mat(:,i)     =  img1 * mat(:,i)
                    endif
                  endif
                enddo
              enddo
            enddo
          endif
          i0 = i0 + nn
        enddo
      enddo
      end

      ! Undoes symmetrization with routine symmetrize.
      ! 
      ! imode0 parameter: see symmetrize
      !
      subroutine desymmetrize(mat,dim1,dim2,imode0)
      use global
      use, intrinsic :: iso_fortran_env
      implicit none
      integer,    intent(in)    :: imode0,dim1,dim2
      complex_dp, intent(inout) :: mat(dim1,dim2)
      complex_dp                :: carr(max(dim1,dim2)),img1
      real_dp                   :: rfac1,rfac2
      integer                   :: imode,typ,indx,dim,sgn
      integer                   :: lcut1,nindx1(0:max(maxlcut,maxlcutm))
      integer                   :: ifac,i,j,i0,itype,ieq,ic,ic1,l,m,n,nn,ishift
      if(dim2==-1) then
        call symm_desymm_atomdiag(mat,dim1,imode0,2)
        return
      endif      
      i = 0
      call imode_def(imode,typ,indx,dim,sgn,imode0,dim1,dim2,i)
      if(imode<3.and.sgn==-1) imode = 3 - imode ! dimensions are flipped if sgn==-1
      img1  = img * sgn                         ! 
      rfac1 = sqrt(0.5d0)
      ic    = 0
      i0    = 0
      do itype = 1,ntype
        if     (indx==0) then ; lcut1 = lcutm(itype) ; nindx1(:lcut1) = nindxm(:lcut1,itype)
        else if(indx==1) then ; lcut1 = lcut(itype)  ; nindx1(:lcut1) = nindx(:lcut1,itype)
        else if(indx==2) then ; lcut1 = lcut(itype)  ; nindx1(:lcut1) = nindx(:lcut1,itype) - 2
        else                  ; lcut1 = lcut(itype)  ; nindx1(:lcut1) = indx - 2
        endif
        if(dim>0) then ; nn = dim
        else           ; nn = sum( [ ((2*l+1)*nindx1(l),l=0,lcut1) ] )
        endif
        do ieq = 1,neq(itype)
          ic  = ic + 1 ; if(typ>0.and.itype/=typ) cycle
          ic1 = pcent(ic,invsym)
          if(ic1>=ic) then
            i = i0
            do l = 0,lcut1
              ifac = -1
              do m = -l,l
                ifac   = -ifac
                rfac2  = rfac1 * ifac
                ishift = (ic1-ic)*nn - 2*m*nindx1(l)
                do n = 1,nindx1(l)
                  i = i + 1
                  j = i + ishift
                  if(ic1/=ic.or.m<0) then
                    if(iand(imode,1)/=0) then
                      carr(:dim2) = mat(i,:)
                      mat(i,:)    = ( carr(:dim2) + img1 * mat(j,:) ) * rfac1
                      mat(j,:)    = ( carr(:dim2) - img1 * mat(j,:) ) * rfac2
                    endif
                    if(iand(imode,2)/=0) then
                      carr(:dim1) = mat(:,i)
                      mat(:,i)    = ( carr(:dim1) - img1 * mat(:,j) ) * rfac1
                      mat(:,j)    = ( carr(:dim1) + img1 * mat(:,j) ) * rfac2
                    endif
                  else if(m==0.and.ifac==-1) then
                    if(iand(imode,1)/=0) then
                      mat(i,:) =  img1 * mat(i,:)
                    endif
                    if(iand(imode,2)/=0) then
                      mat(:,i) = -img1 * mat(:,i)
                    endif
                  endif
                enddo
              enddo
            enddo
          endif
          i0 = i0 + nn
        enddo
      enddo
      end

      ! Interprets imode0 parameter (see symmetrize).
      subroutine imode_def(imode,typ,indx,dim,sgn,imode0,dim1,dim2,typdim)
      use global, only: ntype,neq,nindxm,nindx,lcutm,lcut,nbasp,symm_maxindx
      implicit none
      integer, intent(out)   :: imode,typ,indx,dim,sgn
      integer, intent(inout) :: typdim
      integer, intent(in)    :: imode0,dim1,dim2
      integer                :: itype,n,nn,l
      if(imode0==0) Bug('imode0 is zero.')
      sgn   = sign(1,imode0)
      imode = abs(imode0) - 1
      dim   = imode / (3*(ntype+1)*symm_maxindx)
      imode = imode - (3*(ntype+1)*symm_maxindx) * dim
      indx  = imode / (3*(ntype+1))
      imode = imode - (3*(ntype+1)) * indx      
      typ   = imode / 3
      imode = imode - 3 * typ + 1
      nn    = 0
      do itype = 1,ntype ; if(typ>0.and.itype/=typ) cycle
        if     (indx==0) then ; n = sum([ ((2*l+1)* nindxm(l,itype),  l=0,lcutm(itype)) ])
        else if(indx==1) then ; n = sum([ ((2*l+1)* nindx(l,itype),   l=0,lcut (itype)) ])
        else if(indx==2) then ; n = sum([ ((2*l+1)*(nindx(l,itype)-2),l=0,lcut (itype)) ])
        else                  ; n = sum([ ((2*l+1)*(indx-2),          l=0,lcut (itype)) ])
        endif
        if(dim>0) then
          if(n>dim) Bug('Dimension "dim" too small.')
          n = dim
        endif
        if(typdim==-itype) typdim = n
        if(dim2==-1) then ; nn = max(nn,n)
        else              ; nn = nn + neq(itype) * n
        endif
      enddo
      if(dim2==-1.and.imode/=3) Bug('dim2==-1 but imode/=3.')
      if(sgn==1) then
        if(iand(imode,1)/=0.and.dim1<nn)              Bug('too few rows to (de)symmetrize (imode=1)')
        if(iand(imode,2)/=0.and.dim2<nn.and.dim2/=-1) Bug('too few columns to (de)symmetrize (imode=2)')
      else
        if(iand(imode,1)/=0.and.dim2<nn.and.dim2/=-1) Bug('too few columns to (de)symmetrize (imode=1)')
        if(iand(imode,2)/=0.and.dim1<nn)              Bug('too few rows to (de)symmetrize (imode=2)')
      endif
      end

      ! (De)symmetrization (isub=1 or 2) of matrix "mat", which is atom-diagonal:
      ! From the calling routine, "mat" can have the dimensions
      ! mat(:dim1,:,:ncent) or mat(:dim1,:,:neq(itype))
      ! In this case, the second argument must have the "minimal" dimension. ("dim" in "imode0" must be defined properly.)
      ! The matrix can also have the dimensions mat(:dim1,:) with the second argument a composite variable.
      !
      ! (Called from symmetrize/desymmetrize and calls symmetrize/desymmetrize.)
      subroutine symm_desymm_atomdiag(mat,dim1,imode0,isub)
      use global, only: ntype,neq,symm_maxindx
      use util
      use, intrinsic :: iso_fortran_env      
      implicit none      
      integer,    intent(in)    :: dim1,imode0,isub
      complex_dp, intent(inout) :: mat(dim1,*)
      complex_dp, allocatable   :: mat1(:,:)
      real_dp                   :: err
      integer                   :: imode0_,imode,typ,indx,dim,sgn
      integer                   :: itype,nerr
      integer                   :: i,n,m
      err  = 0
      nerr = 0
      m    = 0
      do itype = 1,ntype
        n = -itype
        call imode_def(imode,typ,indx,dim,sgn,imode0,dim1,-1,n)
        if(typ>0.and.typ/=itype) cycle
        if(imode/=3) Bug('imode/=0 in symm_dsymm_atomdiag.')
        allocate ( mat1(neq(itype)*n,neq(itype)*n) )
        mat1    = 0
        imode0_ = imode0 + 3*(itype-typ) + 3*(ntype+1)*symm_maxindx*(n-dim)
        do i = 0,neq(itype)*n-1,n
          mat1(i+1:i+n,i+1:i+n) = mat(:n,m+i+1:m+i+n)
        enddo
        if(isub==1) then ; call   symmetrize(mat1,neq(itype)*n,neq(itype)*n,imode0_)
        else             ; call desymmetrize(mat1,neq(itype)*n,neq(itype)*n,imode0_)
        endif
        do i = 0,neq(itype)*n-1,n
          mat(:n,m+i+1:m+i+n) = mat1(i+1:i+n,i+1:i+n)
          err                 = err + sum( real ( mat1(    :i,i+1:i+n) )**2 + imag ( mat1(    :i,i+1:i+n) )**2 )
     &                              + sum( real ( mat1(i+n+1:,i+1:i+n) )**2 + imag ( mat1(i+n+1:,i+1:i+n) )**2 )
          nerr                = nerr + (size(mat1,1)-n)*n
        enddo
        write(*,*) 'typeerr',err
        deallocate ( mat1 )
        m = m + n*neq(itype)
      enddo
      err = sqrt(err/nerr) ;       write(*,*) 'err',err
      if(err>1d-10) then
        if(isub==1) then ; Error('Large offsite elements after basis symmetrization: '//Chf(err,'F20.10'))
        else             ; Error('Large offsite elements after basis desymmetrization: '//Chf(err,'F20.10'))
        endif
      endif
      end

c     ------------------

      subroutine matrixtrafo_wan(a,isym,ispin1,ispin2)
      use global, only: nwan,irrep_wan
      use, intrinsic :: iso_fortran_env
      implicit none
      integer,    intent(in)    :: isym,ispin1,ispin2
      complex_dp, intent(inout) :: a(nwan,nwan,nwan,nwan)
      complex_dp                :: b(nwan,nwan,nwan,nwan)
      integer                   :: n1,n2,n3,n4
      b = a
      do n4 = 1,nwan ; do n3 = 1,nwan ; do n2 = 1,nwan ; do n1 = 1,nwan
        a(n1,n2,n3,n4) = sum(      irrep_wan(n1,:,isym,ispin1) *b(:,n2,n3,n4))
      enddo ; enddo ; enddo ; enddo
      b = a
      do n4 = 1,nwan ; do n3 = 1,nwan ; do n2 = 1,nwan ; do n1 = 1,nwan
        a(n1,n2,n3,n4) = sum(conjg(irrep_wan(n2,:,isym,ispin2))*b(n1,:,n3,n4))
      enddo ; enddo ; enddo ; enddo
      b = a
      do n4 = 1,nwan ; do n3 = 1,nwan ; do n2 = 1,nwan ; do n1 = 1,nwan
        a(n1,n2,n3,n4) = sum(conjg(irrep_wan(n3,:,isym,ispin1))*b(n1,n2,:,n4))
      enddo ; enddo ; enddo ; enddo
      b = a
      do n4 = 1,nwan ; do n3 = 1,nwan ; do n2 = 1,nwan ; do n1 = 1,nwan
        a(n1,n2,n3,n4) = sum(      irrep_wan(n4,:,isym,ispin2) *b(n1,n2,n3,:))
      enddo ; enddo ; enddo ; enddo
      end

c     ------------------

# include "trafo.inc"
# define INPLACE
# include "trafo.inc"
# undef INPLACE

c     ------------------

c     "Rotate" for MT part and atom-diagonal array structure
c
c     (De)Symmetrize is not used. (Matrices are complex.)
c
c     typ = 0 : Matrices have dimensions (:dim1,:dim2,:ncent)
c     typ > 0 :         - " -            (:dim1,:dim2,:neq(typ))
c
      subroutine mtrafo_mt(matout,matin,dim1,dim2,ikpt0,isym0,typ,mode,ladd)

      use global
      use wrapper
      use, intrinsic :: iso_fortran_env

      implicit none
      integer,    intent(in)    :: ikpt0,isym0,dim1,dim2,typ,mode
      logical,    intent(in)    :: ladd
      complex_dp, intent(in)    :: matin(dim1,dim2,*)
      complex_dp, intent(inout) :: matout(dim1,dim2,*)
      complex_dp, allocatable   :: matin1(:,:,:)
      complex_dp, allocatable   :: matout1(:,:,:)
      complex_dp                :: dwgnt(-maxlcutm:maxlcutm,-maxlcutm:maxlcutm,0:maxlcutm)
      complex_dp                :: cexp,cdum
      integer                   :: isym,ssym,ikpt1,g(3)
      integer                   :: itype,ieq,ic,ic1,ic0,ic2,ic3,l,nn,i,i1,i2,typ1,typ2,m,m1,lm,nlm

      if(mode<1.or.mode>3)   Bug('mode unknown.')
      if(typ<0.or.typ>ntype) Bug('typ out of range.')
      isym  = abs(isym0)
      ssym  = sign(1,isym0)
      ikpt1 = kptsym(ikpt0,isym)
      g     = nint ( matmul(sym(isym)%rrot,kpt(:,ikpt0))-kpt(:,ikpt1) )
      do l = 0,maxlcutm
        if(ssym==1) then ; dwgnt(-l:l,-l:l,l) = conjg(transpose(dwgn(-l:l,-l:l,l,isym)))
        else             ; dwgnt(-l:l,-l:l,l) =                 dwgn(-l:l,-l:l,l,isym)
        endif
      enddo

c     Input typ
      if(typ==0) then ; typ1 = 1   ; typ2 = ntype ; allocate( matin1(dim1,dim2,ncent),   matout1(dim1,dim2,ncent)    )
      else            ; typ1 = typ ; typ2 = typ   ; allocate( matin1(dim1,dim2,neq(typ)),matout1(dim1,dim2,neq(typ)) )
      endif

c     Check dims
      i = maxval( [ (sum([((2*l+1)*nindxm(l,itype),l=0,lcutm(itype))]) , itype=typ1,typ2) ] )
      if((mode==1.or.mode==3).and.i>dim1) Bug('First dimension too small.')
      if((mode==2.or.mode==3).and.i>dim2) Bug('Second dimension too small.')

      i1 = dim1 ; if(mode/=2) i1 = i
      i2 = dim2 ; if(mode/=1) i2 = i
      nn = size(matin1,3)
# ifdef INV
      matin1(:i1,:i2,:) = matin(:i1,:i2,:nn)
# else
      m = min(dim1,dim2)
      if(isym0>nsymt) then
        if(mode==3) then ; do ic = 1,nn ; matin1(:i,:i,ic) = transpose(matin(:i,:i,ic)) ; enddo
        else             ; matin1(:i1,:i2,:) = conjg(matin(:i1,:i2,:nn))
        endif
      else               ; matin1(:i1,:i2,:) = matin(:i1,:i2,:nn)
      endif
# endif

      cexp = exp(ssym * img * 2*pi * dot_product(kpt(:,ikpt1)+g,sym(isym)%transl))
      ic0  = sum(neq(:typ1-1)) ! atom offset (=0 if typ1=1)

c     Right-multiplication
      if(mode==2.or.mode==3) then
        ! MT
        ic = ic0
        do itype = typ1,typ2
          nlm = sum( [((2*l+1)*nindxm(l,itype),l=0,lcutm(itype))] )
          do ieq = 1,neq(itype)
            ic  = ic + 1
            ic1 = pcent(ic,isym)
            lm  = 0
            if(ssym==1) then ; ic2 = ic -ic0 ; ic3 = ic1-ic0
            else             ; ic2 = ic1-ic0 ; ic3 = ic -ic0
            endif
            if(mode==3) ic3 = ic2 ! Stored in proper atom later
            matout1(:,:nlm,ic3) = 0
            do l = 0,lcutm(itype)
              nn = nindxm(l,itype)
              i1 = lm
              do m1 = -l,l
                i = lm
                do m = -l,l
                  cdum = dwgnt(m,m1,l)
                  if(cdum/=0) matout1(:,i1+1:i1+nn,ic3) = matout1(:,i1+1:i1+nn,ic3) + matin1(:,i+1:i+nn,ic2) * cdum
                  i = i + nn
                enddo
                i1 = i1 + nn
              enddo
              lm = lm + (2*l+1) * nn
            enddo
            cdum                = cexp * exp(-ssym * img * 2*pi * dot_product(g,cent(:,ic1)))
            matout1(:,:nlm,ic3) = matout1(:,:nlm,ic3) * cdum
          enddo
        enddo
        if(mode==3) matin1 = matout1
      endif

c     Left-multiplication
      if(mode==1.or.mode==3) then
        ! MT
        cexp  = conjg(cexp)
        dwgnt = conjg(dwgnt)
        ic    = ic0
        do itype = typ1,typ2
          nlm = sum( [((2*l+1)*nindxm(l,itype),l=0,lcutm(itype))] )
          do ieq = 1,neq(itype)
            ic  = ic + 1
            ic1 = pcent(ic,isym)
            lm  = 0
            if(ssym==1) then ; ic2 = ic -ic0 ; ic3 = ic1-ic0
            else             ; ic2 = ic1-ic0 ; ic3 = ic -ic0
            endif
            matout1(:nlm,:,ic3) = 0
            do l = 0,lcutm(itype)
              nn = nindxm(l,itype)
              i1 = lm
              do m1 = -l,l
                i = lm
                do m = -l,l
                  cdum = dwgnt(m,m1,l)
                  if(cdum/=0) matout1(i1+1:i1+nn,:,ic3) = matout1(i1+1:i1+nn,:,ic3) + matin1(i+1:i+nn,:,ic2) * cdum
                  i = i + nn
                enddo
                i1 = i1 + nn
              enddo
              lm = lm + (2*l+1) * nn
            enddo
            cdum                = cexp * exp(ssym * img * 2*pi * dot_product(g,cent(:,ic1)))
            matout1(:nlm,:,ic3) = matout1(:nlm,:,ic3) * cdum
          enddo
        enddo
      endif

      i  = maxval( [ (sum([((2*l+1)*nindxm(l,itype),l=0,lcutm(itype))]) , itype=typ1,typ2) ] )
      i1 = dim1 ; if(mode/=2) i1 = i
      i2 = dim2 ; if(mode/=1) i2 = i
      nn = size(matin1,3)
# ifdef INV
      if(ladd) then ; matout(:i1,:i2,:nn) = matout(:i1,:i2,:nn) + matout1(:i1,:i2,:)
      else          ; matout(:i1,:i2,:nn) =                       matout1(:i1,:i2,:)
      endif
# else
      if(isym0<-nsymt) then
        if(mode==3) then
          if(ladd) then ; do ic = 1,nn ; matout(:i,:i,ic) = matout(:i,:i,ic) + transpose(matout1(:i,:i,ic)) ; enddo
          else          ; do ic = 1,nn ; matout(:i,:i,ic) =                    transpose(matout1(:i,:i,ic)) ; enddo
          endif
        else
          if(ladd) then ; matout(:i1,:i2,:nn) = matout(:i1,:i2,:nn) + conjg(matout1(:i1,:i2,:))
          else          ; matout(:i1,:i2,:nn) =                       conjg(matout1(:i1,:i2,:))
          endif
        endif
      else
        if(ladd) then ; matout(:i1,:i2,:nn) = matout(:i1,:i2,:nn) + matout1(:i1,:i2,:)
        else          ; matout(:i1,:i2,:nn) =                       matout1(:i1,:i2,:)
        endif
      endif
# endif

      deallocate(matin1,matout1)

      end
