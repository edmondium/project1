c WAVE-FUNCTION PRODUCTS  <M phi | phi >
c
c wavefproducts1_*  : for selfenergy
c wavefproducts2_*  : for susceptibility
c wavefproducts3_*  : for Wannier (phi=Wannier functions)
c wavefproducts4_*  : same as wavefproducts2_* but different output array structure
c wavefproducts*_mt : Muffin-tin part
c wavefproducts*_pw : Plane-wave part
c
c --------------------------------------------------------------------------------

c The parallelization assumes ik1,b1 (respectively k1,b1low,b1up) to be identical over the Mcomm processes (but not b2 etc.).
c It requires a relatively large communication overhead. It can be enabled by the keywords MPIMT and MPIPW (WFPROD).

# ifndef SUBTYPE
#   define SUBTYPE 1
#   include "wavefproducts.f"
#   undef SUBTYPE
#   define SUBTYPE 2
#   include "wavefproducts.f"
#   define SUBWANNIER
#   include "wavefproducts.f"
#   undef SUBWANNIER
#   define SUBORDER      
# endif

# if SUBTYPE != 1 && SUBTYPE != 2
#   error "Unknown subtype"
# endif

c LOAD: cmt and cpw are loaded (different array structure)
c (0) Standard definitions
# define ADJUST_KINDX_IF_NEEDED
# define UNDO_ADJUST_KINDX
# ifdef LOAD
c   (1) cmt,cpw spin index = 1
#   ifndef SUBWANNIER
#     define ISPIN1 1
#     define ISPIN2 1
#   endif
c   (2) macros for redefining kindx to adjust to different kindx ordering of ik2
#   if SUBTYPE == 2 && !defined(SUBWANNIER)
#     undef ADJUST_KINDX_IF_NEEDED
#     undef UNDO_ADJUST_KINDX
#     define ADJUST_KINDX_IF_NEEDED if(.not.storeibz) then ; ikx = kindx(ik2) ; kindx(ik2) = kindx(ik1) ; endif
#     define UNDO_ADJUST_KINDX      if(.not.storeibz) kindx(ik2) = ikx
#   endif
# endif

c Wannier Bloch functions
c (1) Replace cmt->cmtu and cpw->cpwu
# ifdef SUBWANNIER
#   define WSUB _w
# else
#   define WSUB
# endif
c (2) PW coefficients cannot be defined real for Wannier Bloch functions
# if defined SUBWANNIER && defined INV
#   define WSUFFIX _w
#   define INVW
#   undef INV
# else
#   define WSUFFIX
# endif

# ifdef TIMING
#   define Otimer(arg) Rif(otimer/=0) arg
# else
#   define Otimer(arg)
# endif
      
# include "cppmacro.h"

c --------------------------------------------------------------------------------
c
c Muffin-tin part of wave-function products
c
c cprod = < M(k0) phi(b1,k1,s1) | phi(b2,k2,s2) >            ! k2 = k0+k1
c
c SUBTYPE = 1 : k1 fixed, for selfenergy
c SUBTYPE = 2 : k0 fixed, for susceptibility
c
c SUBTYPE     1                 2
c b1          b1(:nb1)          b1low:b1up
c b2          (b2up-b2low)/M+1  b2low:b2up     (striped storage M=Msize if MPI and not LOAD, otherwise M=1)
c k0          k0(:nk)           ik0
c k1          ik1               k1(:nk)
c s1          ispin             ispin1
c s2          ispin             ispin2
c
# if SUBTYPE == 1
      subroutine wavefproducts1_mt(cprod,dim,b1,nb1,ik1,ispin,k0,nk,b2low,b2up,b2inc)
# else
#   if !defined(SUBWANNIER) && !defined(SUBORDER)
      subroutine wavefproducts2_mt(cprod,k1,nk,ik0,ispin1,ispin2,b1low,b1up,b2low,b2up)
#   elif defined(SUBWANNIER)
      subroutine wavefproducts3_mt(cprod,k1,nk,ik0,ispin1,ispin2,b1low,b1up,b2low,b2up)
#   else
      subroutine wavefproducts4_mt(cprod,dim,k1,nk,ik0,ispin1,ispin2,b1low,b1up,b2low,b2up)
#   endif      
# endif
      use global
      use wrapper, only: macmat
      Timing( use timer_util )
      use, intrinsic :: iso_fortran_env
      Mpi( use Mwrapper )
      implicit none
# if SUBTYPE == 1
      integer,     intent(in)    :: nb1,b1(nb1),ik1,ispin,nk,k0(nk),b2low,b2up,b2inc,dim
      MCOMPLEX_dp, intent(inout) :: cprod(dim,nb1, (b2up-b2low) / b2inc + 1 ,nk)
      integer                    :: ispin1,ispin2,ik0
      complex_dp                 :: cexp(nk)      
      complex_dp,  allocatable   :: cmt2(:,:,:,:)
# else
      integer,     intent(in)    :: nk,k1(nk),ik0,ispin1,ispin2,b1low,b1up,b2low,b2up
#   ifdef SUBORDER
      integer,     intent(in)    :: dim
      MCOMPLEX_dp, intent(inout) :: cprod(dim,b2low:b2up,b1low:b1up,nk)      
#   else
      MCOMPLEX_dp, intent(inout) :: cprod(b2low:b2up,b1low:b1up,nk,nbasp)
#   endif
      integer                    :: ik1
      complex_dp,  allocatable   :: cmt2(:,:,:)
      complex_dp                 :: cexp
# endif
      complex_dp,  allocatable   :: cmt1(:,:),mat(:,:)
      integer                    :: lmstart(-maxlcut:maxlcut,0:maxlcut),pnt(maxlmindx)
      integer                    :: ik,ik2,k2(nk),ib,ib1,ib2,s,ikx
      integer                    :: itype,ieq,ic
      integer                    :: l,l1,ll,m,m1,mm,n,n1,nn,lln,ln,ln1,lm,lm1,llm,llm0,nx,nx1,maxlm,maxllm
      logical                    :: def MpiC(mpimt)
      real_dp,     allocatable   :: integral(:,:,:)
      real_dp                    :: gaunt1(0:maxlcut,0:maxlcut,0:maxlcutm,-maxlcut:maxlcut,-maxlcutm:maxlcutm),gnt
      real_dp                    :: gaunt,intgrf
      integer                    :: kptsum
      Mpi( integer              :: Mcnt )

      if(nk==0.or.b2low>b2up) return

# ifdef MPI
#   ifdef SUBWANNIER
      mpimt = .false.
#   else
      mpimt = mpiprod(1)
#   endif
# endif

# if SUBTYPE == 1
      ispin1 = ispin
      ispin2 = ispin
# else
      if(l_soc.and.(ispin1/=1.or.ispin2/=1)) Bug('Wrong spin indices (SOC).')
# endif

      Otimer( call timer_start('WFP mt prep') )

      ! Precalculate Gaunt coefficients
      gaunt1 = 0
      do ll = 0,maxlcutm
        do mm = -ll,ll
          do l = 0,maxlcut
            do m = -l,l
              do l1 = abs(l-ll),min(l+ll,maxlcut),2
                m1                   = mm + m
                gaunt1(l1,l,ll,m,mm) = gaunt(ll,l1,l,mm,m1,m)
              enddo
            enddo
          enddo
        enddo
      enddo

      ! k2 = q+k
      do ik = 1,nk
# if SUBTYPE == 1
        ik0    = k0(ik)
# else
        ik1    = k1(ik)
# endif
        k2(ik) = kptsum(ik0,ik1)
      enddo

      Otimer( call timer_stop('WFP mt prep') )

# if SUBTYPE == 1 || defined(SUBORDER)
      if(dim<nbasp) Bug('Dimension dim too small.')
      cprod(:nbasp,:,:,:) = 0
# else
      cprod = 0
# endif
      llm0  = 0
      ic    = 0
      do itype = 1,ntype
        maxlm  = sum ( [ ((2*l+1)*nindx(l,itype), l=0,lcutp(itype)) ] )
        maxllm = sum ( [ ((2*l+1)*nindxm(l,itype),l=0,lcutm(itype)) ] )

        Otimer( call timer_start('WFP mt prep') )

        allocate ( mat(maxlm,maxllm) )

        ! precalculate radial integrals (->integral1)
        def = .false.
 1      lln = 0 ; llm = 0 ; lm = 0 ; lm1 = 0
        do ll = 0,lcutm(itype)                         ; do nn = 1,nindxm(ll,itype) ; lln = lln + 1 ; ln  = 0
          do l = 0,lcutp(itype)                        ; do n  = 1,nindx(l,itype)   ; ln  = ln  + 1 ; ln1 = 0
            do l1 = abs(ll-l),min(ll+l,lcutp(itype)),2 ; do n1 = 1,nindx(l1,itype)  ; ln1 = ln1 + 1
              if(def) then
                integral(ln,ln1,lln) = intgrf( basm(:,nn,ll,itype) / rgrid(:,itype) *
     &                               ( bas1(:,n,l,itype,ispin1)*bas1(:,n1,l1,itype,ispin2) +
     &                                 bas2(:,n,l,itype,ispin1)*bas2(:,n1,l1,itype,ispin2) ) , itype )
              else
                llm = max(llm,lln)
                lm  = max(lm,ln)
                lm1 = max(lm1,ln1)
              endif
            enddo ; enddo
          enddo ; enddo
        enddo ; enddo
        if(.not.def) then
          def = .true.
          allocate ( integral(lm,lm1,llm) )
          goto 1
        endif

        ! start index
        lm = 0
        do l = 0,lcutp(itype)
          do m = -l,l
            lmstart(m,l) = lm
            lm           = lm + nindx(l,itype)
          enddo
        enddo

        Otimer( call timer_stop('WFP mt prep') )

        if(l_soc) then ; allocate ( cmt1(maxlm,2) )
        else           ; allocate ( cmt1(maxlm,ispin1:ispin1) )
        endif
# if SUBTYPE == 1
        if(l_soc) then ; allocate ( cmt2(maxlm,(b2up-b2low)/b2inc+1,nk,2) )
        else           ; allocate ( cmt2(maxlm,(b2up-b2low)/b2inc+1,nk,ispin2:ispin2) )
        endif
# else        
        if(l_soc) then ; allocate ( cmt2(maxlm,b2low:b2up,2) )
        else           ; allocate ( cmt2(maxlm,b2low:b2up,ispin2:ispin2) )
        endif
# endif

        do ieq = 1,neq(itype)
          ic = ic + 1

# if SUBTYPE == 2
          cexp = exp ( -img * 2*pi * dot_product(kpt(:,ik0),cent(:,ic)) )
# else
#   define IB2  (ib2-b2low)/b2inc+1,ik
#   define B2UP b2up,b2inc
# endif

          do ik = 1,nk
            ik2 = k2(ik)
            ADJUST_KINDX_IF_NEEDED
            do ib2 = b2low, B2UP
              call wavefunction_mt WSUB (cmt2(:, IB2 ,:),maxlm,ic,ib2,ik2,ISPIN2)
            enddo
            UNDO_ADJUST_KINDX

# if SUBTYPE == 2
            ik1 = k1(ik)
            do ib = b1low,b1up
              ib1 = ib
# else
#   undef IB2
#   undef B2UP              
            ik0      = k0(ik)
            cexp(ik) = exp ( -img * 2*pi * dot_product(kpt(:,ik0),cent(:,ic)) )
          enddo
          do ib = 1,nb1
            ib1 = b1(ib)  
# endif

# if defined(LOAD) && SUBTYPE == 1
c            call teststop('LOAD in new wavefproducts.')
            if(storeibz.and.kptp(ik1)/=ik1) Error('Not implemented: LOAD & ik1 not in IBZ.')
            cmt1 = cmtq(:maxlm,ic,ib1,ISPIN1:ISPIN1+nspin3-1)
# else            
            call wavefunction_mt WSUB (cmt1,maxlm,ic,ib1,ik1,ISPIN1)
# endif

c
c           Loop over MPB function

            s = ispin1
            do ! SOC spin loop

            mat = 0 ; Mpi(Mcnt=-1)
            llm = 0
            lln = 0
            do ll = 0,lcutm(itype)
              do mm = -ll,ll

c
c               For each MPB function, calculate  <M phi(ib1) | basis >  ( -> mat, mat1 for SOC)
                lm  = 0
                ln  = 0
                do l = 0,lcutp(itype)
                  nx = nindx(l,itype)
                  do m = -l,l                    

                    m1 = m + mm

                    ln1 = 0
                    do l1 = abs(l-ll),min(l+ll,lcutp(itype)),2                      
                      nx1 = nindx(l1,itype)                      

                      if(l1>=abs(m1)) then                        

                        gnt = gaunt1(l1,l,ll,m,mm)

                        Mpi( Mcnt = Mcnt + 1 )
                        if(gnt/=0 Mpi(.and.(.not.mpimt.or.mod(Mcnt,Msize)==Mrank)) ) then
                          lm1 = lmstart(m1,l1)
                          do nn = 1,nindxm(ll,itype)
                            mat(lm1+1:lm1+nx1,llm+nn) = mat(lm1+1:lm1+nx1,llm+nn) +
     &                                                  matmul(cmt1(lm+1:lm+nx,s),integral(ln+1:ln+nx,ln1+1:ln1+nx1,lln+nn)) * gnt
                          enddo
                        endif

                      endif

                      ln1 = ln1 + nx1
                    enddo ! l1

                    lm = lm + nx
                  enddo ! m
                  ln = ln + nx
                enddo ! l
                llm = llm + nindxm(ll,itype)
              enddo ! mm
              lln = lln + nindxm(ll,itype)
            enddo ! ll
              
            Mpi( if(mpimt) call Msum(mat) )

c
c           Multiply  cprod = <M phi(ib1) | phi(ib2)> = <M phi(ib1) | basis> <basis | phi(ib2)>            

            if(s==ispin1) s = ispin2
            if(llm/=maxllm) Bug('MPB index count error.')
            llm = llm0 + llm
            
# if SUBTYPE == 1
            do ik = 1,nk
              NoInv( cprod(llm0+1:llm,ib,:,ik) = cprod(llm0+1:llm,ib,:,ik) + cexp(ik) * macmat(mat,cmt2(:,:,ik,s) ) )
              Inv(   call symmetrize_to_cprod                              ( cexp(ik) * macmat(mat,cmt2(:,:,ik,s) ) ) )
            enddo
# elif defined(SUBORDER)
            NoInv( cprod(llm0+1:llm,:,ib,ik)   = cprod(llm0+1:llm,:,ib,ik) + cexp * macmat(mat,cmt2(:,:,s) ) )
            Inv(   call symmetrize_to_cprod                                ( cexp * macmat(mat,cmt2(:,:,s) ) ) )
# else                  
            NoInv( cprod(:,ib,ik,llm0+1:llm)   = cprod(:,ib,ik,llm0+1:llm) + cexp * conjg(macmat(cmt2(:,:,s),mat))   )
            Inv(   call symmetrize_to_cprod                                ( cexp * conjg(macmat(cmt2(:,:,s),mat)) ) )
# endif

            if(.not.l_soc.or.s==2) exit
            s = 2
            enddo ! SOC spin loop

          enddo     ! ib(1)
# if SUBTYPE == 2
          enddo     ! ik(1)
# endif
          llm0 = llm
        enddo ! ieq

        deallocate(integral,cmt1,cmt2,mat)

      enddo ! itype

# ifdef SUBWANNIER
      if(storeibz) call transform_wan(cprod,k1,k2,nk,ispin1,ispin2,b2up-b2low+1,nbasp)
# endif

      contains

# ifdef INV
#   if SUBTYPE == 1
#     define CPROD(llm) cprod(llm+1:llm+nn,ib,:,ik)
#     define MAT mat(lm+1:lm+nn,:)
#   elif defined(SUBORDER)
#     define CPROD(llm) cprod(llm+1:llm+nn,:,ib,ik)
#     define MAT mat(lm+1:lm+nn,:)
#   else
#     define CPROD(llm) cprod(:,ib,ik,llm+1:llm+nn)      
#     define MAT mat(:,lm+1:lm+nn)
#   endif
      subroutine symmetrize_to_cprod(mat) 
      implicit none ! uses ic1, llm0, ic, maxllm, itype, ib, ik from parent routine
      complex_dp, intent(in) :: mat(:,:)
      integer                :: llm,llm1,ic1,sgn,nn,ll,mm,lm
      real_dp,    parameter  :: sq = sqrt(0.5d0)
      ic1 = pcent(ic,invsym)
      lm  = 0
      do ll = 0,lcutm(itype)
        do mm = -ll,ll
          nn   = nindxm(ll,itype)
          llm  = llm0 + lm
          llm1 = llm - 2*mm*nn + (ic1-ic)*maxllm          
          if(llm1/=llm) then
            if(llm1>llm) then
              CPROD(llm)  = CPROD(llm)  +       MAT * sq
              CPROD(llm1) = CPROD(llm1) - img * MAT * sq
            else
              sgn         = (-1)**(mm+ll)
              CPROD(llm1) = CPROD(llm1) +       MAT * (sq*sgn)
              CPROD(llm)  = CPROD(llm)  + img * MAT * (sq*sgn)
            endif
          else
            if(mod(ll,2)==0) then
              CPROD(llm)  = CPROD(llm)  +       MAT
            else
              CPROD(llm)  = CPROD(llm)  - img * MAT
            endif
          endif
          lm = lm + nn
        enddo
      enddo
      end subroutine symmetrize_to_cprod
#   undef CPROD
#   undef MAT
# endif

      end

c --------------------------------------------------------------------------------
c
c Plane-wave part of wave-function products
c
c cprod = < k0+G0 phi(b1,k1) | phi(b2,k2) >
c
c
c (1) Fast Fourier transformation (FFT)
c
c FULLPW: Multiply with step function phi(b1,k1)*theta -> phi(b1,k1)
c
c Transform phi(b1,k1,G) -> phi(b1,k1,R)
c           phi(b2,k2,G) -> phi(b2,k2,R)
c Multiply  phi(b1,k1,R)*phi(b2,k2,R) -> prod(b1,b2,k1,k2,R)
c Backtrafo prod(b1,b2,k1,k2,R) -> cprod(b1,b2,k1,k2,G)
c
c
c (2) Explicit Convolutions (default)
c
c FULLPW: Multiply with step function
c         phi~(b1,k1,G) = SUM(G') phi(b1,k1,G') * theta(G-G')
c
c Convolute   < k0+G0     phi~(b1,k1,G1) | phi(b2,k2,G2)      >  ! k2 = k1+k0+g (g=backfolding vector)
c           = < k0+G0     phi~(b1,k1,G1) | phi(b2,k0+k1+g,G2) >  ! k0+g+k1 on the left balances with k0+k1+g on the right
c           = < k0+g+G0-g phi~(b1,k1,G1) | phi(b2,k0+k1+g,G2) >  ! -> G1 = G2-G0+g
c           = SUM(G2) conjg[phi~(b1,k1,G2-G0+g)] * phi(b2,k2,G2) -> cprod(b1,b2,k1,k2,G0)
c

# if SUBTYPE == 1
      subroutine wavefproducts1_pw(cprod,dim,b1,nb1,ik1,ispin,k0,nk,b2low,b2up,b2inc)
# else
#   if !defined(SUBWANNIER) && !defined(SUBORDER)
      subroutine wavefproducts2_pw(cprod,k1,nk,ik0,ispin1,ispin2,b1low,b1up,b2low,b2up)
#   elif defined(SUBWANNIER)
      subroutine wavefproducts3_pw(cprod,k1,nk,ik0,ispin1,ispin2,b1low,b1up,b2low,b2up)
#   else
      subroutine wavefproducts4_pw(cprod,dim,k1,nk,ik0,ispin1,ispin2,b1low,b1up,b2low,b2up)
#   endif
# endif
      use global
      use wrapper, only: dotprod,macmat
      use m_fft
      use, intrinsic :: iso_fortran_env
      Timing( use timer_util )
      Mpi( use Mwrapper )
      implicit none
# if SUBTYPE == 1
      integer,     intent(in)    :: nb1,b1(nb1),ik1,ispin,nk,k0(nk),b2low,b2up,b2inc,dim
      MCOMPLEX_dp, intent(inout) :: cprod(dim,nb1, (b2up-b2low) / b2inc + 1 ,nk)
      integer                    :: ispin1,ispin2
      integer                    :: ik0,ib2_
# else
      integer,     intent(in)    :: nk,k1(nk),ik0,ispin1,ispin2,b1low,b1up,b2low,b2up
#   ifdef SUBORDER
      integer,     intent(in)    :: dim
      MCOMPLEX_dp, intent(inout) :: cprod(dim,b2low:b2up,b1low:b1up,nk)
#   else      
      MCOMPLEX_dp, intent(inout) :: cprod(b2low:b2up,b1low:b1up,nk,ngptm(ik0))
#   endif
      integer                    :: ik1,nb1,nb2
# endif
      complex_dp,  allocatable   :: rfft(:,:,:,:)
      MCOMPLEX_dp, allocatable   :: cpw1(:),cpw1s(:,:),cpw_1(:,:,:),cpw0(:,:,:),cpw3(:,:),cpw2(:,:,:)
      real_dp                    :: sqrtvoli
      integer,     allocatable   :: point(:,:,:),gpt0(:,:),pgptinv(:)
      integer                    :: k2(nk),ik,ik2,s,ib1,ib2,ib,ikx
      integer                    :: gpt2(3,nk),g(3),g1(3),ig,ig1,ig2
      integer                    :: ngpt0,ig0
      integer                    :: ngpt1,ngpt2,pgpt1(maxgpt),pgpt2(maxgpt),pg(maxgpt)
      integer                    :: i,j,k,lower,upper
      integer                    :: def
      integer                    :: kptsum      
# if defined INV || defined INVW
      real_dp,     allocatable   :: step(:)
# else
      complex_dp,  allocatable   :: step(:)
# endif
      Mpi( logical               :: mpipw )

      real cputime

# ifdef MPI
#   ifdef SUBWANNIER
      mpipw = .false.
#   else
      mpipw = mpiprod(2)
#   endif
# endif

      if(nk==0.or.b2low>b2up) return

      call cpu_time(cputime)

      sqrtvoli = 1/sqrt(vol)

# if SUBTYPE == 1
      ispin1 = ispin
      ispin2 = ispin
# else
      if(l_soc.and.(ispin1/=1.or.ispin2/=1)) Bug('Wrong spin indices (SOC).')
      nb1 = b1up - b1low + 1
      nb2 = b2up - b2low + 1
# endif

      if(gcutf/=0) then
        allocate ( cpw1(maxgpt) )
        if(l_soc) allocate ( cpw1s(maxgpt,2) )
      endif

      ! k2 = q+k
      do ik = 1,nk
# if SUBTYPE == 1
        ik0 = k0(ik)
# else
        ik1 = k1(ik)
# endif
        k2(ik)     = kptsum(ik0,ik1)
        gpt2(:,ik) = nint ( kpt(:,k2(ik)) - kpt(:,ik0) - kpt(:,ik1) )
        if(any(abs(kpt(:,ik1)+kpt(:,ik0)+gpt2(:,ik)-kpt(:,k2(ik)))>1d-12)) Bug('error')
      enddo

      Otimer( call timer_start('WFP pw init') )

# if SUBTYPE == 1
      if(dim<minval(nbasm(k0))) Bug('Dimension dim too small.')
      if(dim==nbasp) Error('No IPWs?')
      cprod(nbasp+1:,:,:,:) = 0
# elif defined(SUBORDER)
      if(dim<nbasm(ik0)) Bug('Dimension dim too small.')
      if(dim==nbasp) Error('No IPWs?')
      cprod(nbasp+1:,:,:,:) = 0      
# else
      cprod = 0
# endif

      Otimer( call timer_stop('WFP pw init') )

      if(.not.fullpw) allocate ( pgptinv(ngptall) )

c
c (1) Fast Fourier Transformation
c

      if(gcutf/=0) then
        ! Calculate step function with Fourier components |G|<=2*gcut+gcutm
        if(fullpw) then
          cffts WSUFFIX = 0
          do k = 0,nfft(3)-1     ; g(3) = k ; if(g(3)>nfft(3)/2) g(3) = g(3) - nfft(3)
            do j = 0,nfft(2)-1   ; g(2) = j ; if(g(2)>nfft(2)/2) g(2) = g(2) - nfft(2)
              do i = 0,nfft(1)-1 ; g(1) = i ; if(g(1)>nfft(1)/2) g(1) = g(1) - nfft(1)
                g1                    = matmul(imat,g) ; if(sum(matmul(rlat,g1)**2)>(2*gcut+gcutm)**2) cycle
                cffts WSUFFIX (i,j,k) = cstep(g1(1),g1(2),g1(3))
              enddo
            enddo
          enddo
          call dfftw_execute_dft Inv(_r2c) (plan_ffts WSUFFIX,cffts WSUFFIX,rffts)
        endif

# if SUBTYPE == 2
c       Either ib1 or ib2 as the inner loop (two blocks of code)
c       --- First block ---
        if(b2up-b2low>b1up-b1low) then
          i = 1 ; if(l_soc) i = 2
          allocate ( rfft(size(rfft1,1),size(rfft1,2),size(rfft1,3),nb1*i) )
          rfft = 0
          do ik = 1,nk
            ik1 = k1(ik)
            ik2 = k2(ik)
            do ib = 1,nb1
              ib1 = b1low + ib - 1
# else
        i = 1 ; if(l_soc) i = 2
        allocate ( rfft(size(rfft1,1),size(rfft1,2),size(rfft1,3),nb1*i) )
        rfft = 0
        do ib = 1,nb1
          ib1 = b1(ib)
# endif

# if defined(LOAD) && SUBTYPE == 1
          if(storeibz.and.kptp(ik1)/=ik1) Error('Not implemented: LOAD & ik1 not in IBZ.')
          if(l_soc) then ; cpw1s = cpwq(:,ib1,:)
          else           ; cpw1  = cpwq(:,ib1,ISPIN1)
          endif
# else
          if(l_soc) then ; call wavefunction_pw WSUB (cpw1s,maxgpt,ib1,ik1,ISPIN1)
          else           ; call wavefunction_pw WSUB (cpw1, maxgpt,ib1,ik1,ISPIN1)
          endif
# endif

          do s = 1,nspin3 ; if(l_soc) cpw1 = cpw1s(:,s) ! SOC spin loop
          cfft1 WSUFFIX = 0
          do ig1 = 1,ngpt(ik1)
            g                              = matmul(imati,gpt(:,pgpt(ig1,ik1)))
            g                              = modulo(g,nfft)
            cfft1 WSUFFIX (g(1),g(2),g(3)) = cpw1(ig1)
          enddo
          call dfftw_execute_dft Inv(_r2c) (plan_fft1 WSUFFIX,cfft1 WSUFFIX,rfft1)
          rfft1 = conjg ( rfft1 * sqrtvoli )
          if(fullpw) then
            rfft1 = rfft1 * rffts
          endif
          if(l_soc) then ; rfft(:,:,:,ib+(s-1)*nb1) = rfft1
          else           ; rfft(:,:,:,ib)           = rfft1
          endif
          enddo ! SOC spin loop
        enddo ! ib1
        
# if SUBTYPE == 1
        do ik = 1,nk
          ik0  = k0(ik)
          ik2  = k2(ik)
          ib2_ = 0
          do ib2 = b2low,b2up,b2inc
            ib2_ = ib2_ + 1
# else
          do ib2 = b2low,b2up
# endif

            ADJUST_KINDX_IF_NEEDED
            if(l_soc) then ; call wavefunction_pw WSUB (cpw1s,maxgpt,ib2,ik2,ISPIN2)
            else           ; call wavefunction_pw WSUB (cpw1, maxgpt,ib2,ik2,ISPIN2)
            endif
            UNDO_ADJUST_KINDX

            do s = 1,nspin3 ; if(l_soc) cpw1 = cpw1s(:,s) ! SOC spin loop
            cfft2 WSUFFIX = 0
            do ig2 = 1,ngpt(ik2)
              g                              = matmul(imati,gpt(:,pgpt(ig2,ik2)))
              g                              = modulo(g,nfft)
              cfft2 WSUFFIX (g(1),g(2),g(3)) = cpw1(ig2)
            enddo
            call dfftw_execute_dft Inv(_r2c) (plan_fft2 WSUFFIX,cfft2 WSUFFIX,rfft2)
            rfft2 = rfft2 / nnfft
            do ib = 1,nb1
# if SUBTYPE == 2
              ib1 = b1low + ib - 1
# endif
              if(l_soc) then ; rfftp = rfft(:,:,:,ib+(s-1)*nb1) * rfft2
              else           ; rfftp = rfft(:,:,:,ib)           * rfft2
              endif
              call dfftw_execute_dft Inv(_c2r) (plan_fftp WSUFFIX,rfftp,cfftp WSUFFIX)
              do ig = 1,ngptm(ik0)
                g                          = matmul(imati,gptm(:,pgptm(ig,ik0))-gpt2(:,ik))
                g                          = modulo(g,nfft)
# if SUBTYPE == 1
                cprod(nbasp+ig,ib,ib2_,ik) = cprod(nbasp+ig,ib,ib2_,ik) + cfftp WSUFFIX (g(1),g(2),g(3))
# elif defined(SUBORDER)
                cprod(nbasp+ig,ib2,ib1,ik) = cprod(nbasp+ig,ib2,ib1,ik) + cfftp WSUFFIX (g(1),g(2),g(3))              
# else
                cprod(ib2,ib1,ik,ig)       = cprod(ib2,ib1,ik,ig)       + cfftp WSUFFIX (g(1),g(2),g(3))
# endif
              enddo
            enddo ! ib(1)
            enddo ! SOC spin loop
          enddo   ! ib2
        enddo     ! ik0 / ik

# if SUBTYPE == 2
c       --- Second block ---
        else
          i = 1 ; if(l_soc) i = 2
          allocate ( rfft(size(rfft2,1),size(rfft2,2),size(rfft2,3),nb2*i) )
          rfft = 0
          do ik = 1,nk
            ik1 = k1(ik)
            ik2 = k2(ik)
            do ib2 = b2low,b2up
              ib = ib2 - b2low + 1
              ADJUST_KINDX_IF_NEEDED
              if(l_soc) then ; call wavefunction_pw WSUB (cpw1s,maxgpt,ib2,ik2,ISPIN2)
              else           ; call wavefunction_pw WSUB (cpw1, maxgpt,ib2,ik2,ISPIN2)
              endif
              UNDO_ADJUST_KINDX
              do s = 1,nspin3 ; if(l_soc) cpw1 = cpw1s(:,s) ! SOC spin loop
              cfft2 WSUFFIX = 0
              do ig2 = 1,ngpt(ik2)
                g                              = matmul(imati,gpt(:,pgpt(ig2,ik2)))
                g                              = modulo(g,nfft)
                cfft2 WSUFFIX (g(1),g(2),g(3)) = cpw1(ig2)
              enddo
              call dfftw_execute_dft Inv(_r2c) (plan_fft2 WSUFFIX,cfft2 WSUFFIX,rfft2)
              if(fullpw) then
                rfft2 = rfft2 * rffts
              endif
              if(l_soc) then ; rfft(:,:,:,ib+(s-1)*nb2) = rfft2
              else           ; rfft(:,:,:,ib)           = rfft2
              endif
              enddo ! SOC spin loop
            enddo ! ib2
            do ib1 = b1low,b1up
              if(l_soc) then ; call wavefunction_pw WSUB (cpw1s,maxgpt,ib1,ik1,ISPIN1)
              else           ; call wavefunction_pw WSUB (cpw1, maxgpt,ib1,ik1,ISPIN1)
              endif
              do s = 1,nspin3 ; if(l_soc) cpw1 = cpw1s(:,s) ! SOC spin loop
              cfft1 WSUFFIX = 0
              do ig1 = 1,ngpt(ik1)
                g                              = matmul(imati,gpt(:,pgpt(ig1,ik1)))
                g                              = modulo(g,nfft)
                cfft1 WSUFFIX (g(1),g(2),g(3)) = cpw1(ig1)
              enddo
              call dfftw_execute_dft Inv(_r2c) (plan_fft1 WSUFFIX,cfft1 WSUFFIX,rfft1)
              rfft1 = conjg ( rfft1 * sqrtvoli ) / nnfft
              do ib = 1,nb2
                ib2 = b2low + ib - 1
                if(l_soc) then ; rfftp = rfft(:,:,:,ib+(s-1)*nb2) * rfft1
                else           ; rfftp = rfft(:,:,:,ib)           * rfft1
                endif
                call dfftw_execute_dft Inv(_c2r) (plan_fftp WSUFFIX,rfftp,cfftp WSUFFIX)
                do ig = 1,ngptm(ik0)
                  g                    = matmul(imati,gptm(:,pgptm(ig,ik0))-gpt2(:,ik))
                  g                    = modulo(g,nfft)
#   ifdef SUBORDER
                  cprod(nbasp+ig,ib2,ib1,ik) = cprod(nbasp+ig,ib2,ib1,ik) + cfftp WSUFFIX (g(1),g(2),g(3))
#   else
                  cprod(ib2,ib1,ik,ig)       = cprod(ib2,ib1,ik,ig)       + cfftp WSUFFIX (g(1),g(2),g(3))
#   endif
                enddo
              enddo ! ib(2)
              enddo ! SOC spin loop
            enddo   ! ib1
          enddo     ! ik
        endif
# endif
        deallocate ( rfft )
        deallocate ( cpw1 )
        if(allocated(cpw1s)) deallocate ( cpw1s )

      else

c
c (2) Explicit Convolutions
c

# if SUBTYPE == 2
        do ik = 1,nk
          ik1 = k1(ik) ; ngpt1 = ngpt(ik1) ; pgpt1 = pgpt(:,ik1)
          ik2 = k2(ik) ; ngpt2 = ngpt(ik2) ; pgpt2 = pgpt(:,ik2)
          allocate ( cpw2(ngpt2,b2low:b2up,nspin3),cpw3(ngpt2,ngptm(ik0)) )

          Otimer( call timer_start('WFP pw wavefunction') )
          ADJUST_KINDX_IF_NEEDED
          do ib2 = b2low,b2up
            call wavefunction_pw WSUB(cpw2(:,ib2,:),ngpt2,ib2,ik2,ISPIN2)
          enddo
          UNDO_ADJUST_KINDX
          Otimer( call timer_stop('WFP pw wavefunction') )
# else
        ngpt1 = ngpt(ik1)
        pgpt1 = pgpt(:,ik1)
# endif

        if(.not.fullpw) then
          pgptinv                       = 0
          pgptinv(pgpt(:ngpt(ik1),ik1)) = [ (i,i=1,ngpt(ik1)) ]
        endif

c        write(*,'(A'NoA) '1';call cpu_done(cputime)

c
c       Prepare multiplication with step function (->point,gpw0)
        Otimer( call timer_start('WFP pw gpt0') )
        def   =  0
        g1    = -1
 1      ngpt0 =  0
# if SUBTYPE == 1
        do ik = 1,nk
          ik0 = k0(ik)
          ik2 = k2(ik) ; ngpt2 = ngpt(ik2) ; pgpt2 = pgpt(:,ik2)
# endif
        do ig2 = 1,ngpt2
          do ig = 1,ngptm(ik0)
            g = gpt(:,pgpt2(ig2)) - gptm(:,pgptm(ig,ik0)) + gpt2(:,ik)
            if(def==0) then
              g1(1) = max(g1(1),abs(g(1)))
              g1(2) = max(g1(2),abs(g(2)))
              g1(3) = max(g1(3),abs(g(3)))
            else
              if(point(g(1),g(2),g(3))==0) then
                if(fullpw) then
                  ngpt0                 = ngpt0 + 1
                  point(g(1),g(2),g(3)) = ngpt0
                else if(sum(matmul(rlat,kpt(:,ik1)+g)**2)<=gcut**2) then
                  ngpt0                 = ngpt0 + 1
                  point(g(1),g(2),g(3)) = pgptinv ( pntgpt(g(1),g(2),g(3)) )
                endif
                if(def==2) then
                  gpt0(:,ngpt0) = g
                endif
              endif
            endif
          enddo
        enddo
# if SUBTYPE == 1
        enddo
# endif
        if(def==0) then
          allocate ( point(-g1(1):g1(1),-g1(2):g1(2),-g1(3):g1(3)) )
          point = 0
          def   = 1
          goto 1
        else if(def==1.and.fullpw) then
          allocate ( gpt0(3,ngpt0) )
          point = 0
          def   = 2
          goto 1
        endif

        if(fullpw) allocate ( cpw0(ngpt0,nspin3,nb1) )
        allocate ( cpw_1(ngpt1,nspin3,nb1) )

        Otimer( call timer_stop('WFP pw gpt0') )

c
c       Get coefficients for states b1 (->cpw_1)

        Otimer( call timer_start('WFP pw wavefunction') )
        do ib = 1,nb1
# if SUBTYPE == 1
          ib1 = b1(ib)
# else
          ib1 = b1low + ib - 1
# endif

#   if defined(LOAD) && SUBTYPE == 1
          if(storeibz.and.kptp(ik1)/=ik1) Error('Not implemented: LOAD & ik1 not in IBZ.')
          if(l_soc) then ; cpw_1(:,:,ib) = cpwq(:ngpt1,ib1,:)
          else           ; cpw_1(:,1,ib) = cpwq(:ngpt1,ib1,ISPIN1)
          endif
#   else
          call wavefunction_pw WSUB (cpw_1(:,:,ib),ngpt1,ib1,ik1,ISPIN1)
#   endif

        enddo
        Otimer( call timer_stop('WFP pw wavefunction') )

c
c       fullpw: Multiply cpw_1 with step function (-> cpw0)
        if(fullpw) then
          Otimer ( call timer_start('WFP pw cpw0') )
          allocate(step(ngpt1))
          lower = 1
          upper = ngpt0
          Mpi( if(mpipw) then ; MrangeDef1(lower,upper,ngpt0) ; endif )
          do ig0 = lower,upper
            do ig1 = 1,ngpt1
              g         = gpt0(:,ig0) - gpt(:,pgpt1(ig1))
              step(ig1) = cstep(g(1),g(2),g(3))
            enddo
            cpw0(ig0,:,:) = reshape ( matmul(step,reshape(cpw_1,(/ngpt1,nspin3*nb1/))) , (/nspin3,nb1/) )
          enddo
          deallocate(step)
# ifdef MPI
          if(mpipw) then ; MrangeDistr(cpw0(McolD1(ngpt0,i),:,:),i) ; endif
# endif
          Otimer ( call timer_stop('WFP pw cpw0') )
        endif
        
c
c       Loop over ib1

        do ib = 1,nb1
# if SUBTYPE == 1
          ib1 = b1(ib)
# else
          ib1 = b1low + ib - 1
# endif        

# if SUBTYPE == 1
          do ik = 1,nk
            Otimer ( call timer_start('WFP pw wavefunction') )
            ik0 = k0(ik)
            ik2 = k2(ik) ; ngpt2 = ngpt(ik2) ; pgpt2 = pgpt(:,ik2)
            allocate ( cpw2(ngpt2,size(cprod,3),nspin3),cpw3(ngpt2,maxgptm) )

            ib2_ = 0
            do ib2 = b2low,b2up,b2inc
              ib2_ = ib2_ + 1
              call wavefunction_pw WSUB(cpw2(:,ib2_,:),ngpt2,ib2,ik2,ISPIN2)
            enddo
            Otimer ( call timer_stop('WFP pw wavefunction') )
# endif         
          do s = 1,nspin3 ! SOC spin loop

          ! define cpw3 with which we multiply phi(ib2,ik2)
          if(fullpw) then
            Otimer ( call timer_start('WFP pw cpw3') )
            lower = 1
            upper = ngptm(ik0)
            Mpi( if(mpipw) then ; MrangeDef1(lower,upper,ngptm(ik0)) ; endif )
            do ig = lower,upper
              do ig2 = 1,ngpt2
                g            = gpt(:,pgpt2(ig2)) - gptm(:,pgptm(ig,ik0)) + gpt2(:,ik)
                cpw3(ig2,ig) = cpw0(point(g(1),g(2),g(3)),s,ib)
              enddo
            enddo
# ifdef MPI            
            if(mpipw) then ; MrangeDistr(cpw3(:,McolD1(ngptm(ik0),i)),i) ; endif
# endif
            Otimer ( call timer_stop('WFP pw cpw3') )
          endif

c
c         Loop over mixed-basis IPWs and define cprod
          Otimer ( call timer_start('WFP pw cprod') )
          if(fullpw) then
c           (1) FULLPW
# if SUBTYPE == 1
            cprod(nbasp+1:nbasp+ngptm(ik0),ib,:,ik)  = cprod(nbasp+1:nbasp+ngptm(ik0),ib,:,ik)  + sqrtvoli*macmat(cpw3,cpw2(:,:,s))
# else
#   ifdef SUBORDER            
            cprod(nbasp+1:nbasp+ngptm(ik0),:,ib1,ik) = cprod(nbasp+1:nbasp+ngptm(ik0),:,ib1,ik) + sqrtvoli*macmat(cpw3,cpw2(:,:,s))
#   else
            cprod(:,ib1,ik,:) = cprod(:,ib1,ik,:) + sqrtvoli * MCONJG( macmat(cpw2(:,:,s),cpw3) )
#   endif
# endif
          else
c           (2) APPROXPW (efficiency can be improved but probably not worth it)
            ! define cpw3 with which we multiply phi(ib2,ik2)
            do ig = 1,ngptm(ik0)
              k = 0
              do ig2 = 1,ngpt2
                g         = gpt(:,pgpt2(ig2)) - gptm(:,pgptm(ig,ik0)) + gpt2(:,ik)
                ig1       = point(g(1),g(2),g(3))
                if(ig1/=0) then
                  k         = k + 1
                  cpw3(k,1) = cpw_1(ig1,s,ib)
                  pg(k)     = ig2
                endif
              enddo
              ! multiply with phi(ib2,ik2)
# if SUBTYPE == 1
              do ib2 = 1,size(cprod,3)
                cprod(nbasp+ig,ib,ib2,ik)  = cprod(nbasp+ig,ib,ib2,ik)  + dotprod(cpw3(:k,1),cpw2(pg(:k),ib2,s)) * sqrtvoli
# else
              do ib2 = b2low,b2up
#   ifdef SUBORDER
                cprod(nbasp+ig,ib2,ib1,ik) = cprod(nbasp+ig,ib2,ib1,ik) + dotprod(cpw3(:k,1),cpw2(pg(:k),ib2,s)) * sqrtvoli
#   else
                cprod(ib2,ib1,ik,ig)       = cprod(ib2,ib1,ik,ig)       + dotprod(cpw3(:k,1),cpw2(pg(:k),ib2,s)) * sqrtvoli
#   endif
# endif
              enddo
            enddo
          endif
          Otimer ( call timer_stop('WFP pw cprod') )


          enddo ! SOC spin loop
# if SUBTYPE == 1
            deallocate ( cpw2,cpw3 )
          enddo ! ik(0)
# endif

        enddo ! ib(1)

        deallocate ( cpw_1 )

        if(fullpw) deallocate ( cpw0,gpt0 )
        deallocate ( point )
# if SUBTYPE == 2
          deallocate ( cpw2,cpw3 )
        enddo ! ik(1)
# endif

      endif

      if(.not.fullpw) deallocate ( pgptinv )

# ifdef SUBWANNIER
      if(storeibz) call transform_wan(cprod,k1,k2,nk,ispin1,ispin2,b2up-b2low+1,ngptm(ik0))
# endif

      end

c --------------------------------------------------------------------------------

# ifdef SUBWANNIER

c --------------------------------------------------------------------------------

c If storeibz=.true., wavefproducts returns < M(k0) P(k1) phi(b1,k1p) | P(k2) phi(b2,k2p) >,
c where P(k1) [P(k2)] rotates from k1p (k2p), the parent of k1 (k2), to k1 (k2).
c
c In the case of Wannier functions, this is < M(k0) P(k1) wan(n1,k1p) | P(k2) wan(n2,k2p) >,
c but we need < M(k0) wan(n1,k1) | wan(n2,k2) >.
c
c Here, we transform both Wannier functions using:      
c | wan(n,k) > = sum(b) U(b,n,k) | phi(b,k) > = sum(b) U(b,n,k) | P(k) phi(b,kp) >
c              = sum(b) U(b,n,k) [ sum(n') U*(b,n',kp)   | P(k) wan(n',kp) > ]
c              = sum(n') [ sum(b) U(b,n,k) U*(b,n',kp) ] | P(k) wan(n',kp) >
c                \-----       = trafo(n,n')       -----/
c
      subroutine transform_wan(cprod,k1,k2,nk,s1,s2,dim1,dim4)
      use global
      use wrapper, only: matmat
      use, intrinsic :: iso_fortran_env
      implicit none
      integer,    intent(in)    :: nk,k1(nk),k2(nk),dim1,dim4,s1,s2
      complex_dp, intent(inout) :: cprod(dim1,nwan,nk,dim4)
      complex_dp                :: trafo(nwan,nwan)
      integer                   :: ik,ik1,ik2,ik1p,ik2p,ibas
      do ik = 1,nk
        ik1 = k1(ik) ; ik1p = kptp(ik1)
        ik2 = k2(ik) ; ik2p = kptp(ik2)
        if(ik1>nkpti) then
          if(symkpt(ik1)>nsymt) then
            trafo = matmat ( conjg(transpose(uwan(:,:,ik1,s1))) ,    conjg(uwan(:,:,ik1p,s1)) )
          else
            trafo = matmat ( conjg(transpose(uwan(:,:,ik1,s1))) ,          uwan(:,:,ik1p,s1)  ) * phase(ik1)
          endif
          do ibas = 1,dim4
            cprod(:,:,ik,ibas) = matmat ( cprod(:,:,ik,ibas) , transpose(trafo) )
          enddo
        endif
        if(ik2>nkpti) then
          if(symkpt(ik2)>nsymt) then
            trafo = matmat ( transpose(uwan(:,:,ik2,s2)) ,       uwan(:,:,ik2p,s2) )
          else
            trafo = matmat ( transpose(uwan(:,:,ik2,s2)) , conjg(uwan(:,:,ik2p,s2) ) ) * conjg(phase(ik2))
          endif
          do ibas = 1,dim4
            cprod(:,:,ik,ibas) = matmat ( trafo , cprod(:,:,ik,ibas) )
          enddo
        endif
      enddo
      end

c --------------------------------------------------------------------------------

c   Special wavefunction routines for Wannier Bloch functions (SUBWANNIER)

      subroutine wavefunction_mt_w(cmtout,dim,ic,iband,ikpt,ispin)

      use global

      use, intrinsic :: iso_fortran_env
      implicit none
      integer,    intent(in)  :: ic,ikpt,iband,ispin,dim
      complex_dp, intent(out) :: cmtout(dim,*)
      complex_dp, pointer_cnt :: cmtin(:,:)
      integer                 :: itype,ieq,icent,i,ic0

      if(ispin<1.or.ispin>nspin2) Bug('ispin out of range')
      if(l_soc) Error('Not implemented for SOC.')

      if(storeibz.and.kptp(ikpt)/=ikpt) then
        if     (ic==0) then ; ic0 = ic                              ; icent = 1
        else if(ic> 0) then ; ic0 = pcent(ic,sym(symkpt(ikpt))%inv) ; icent = ic0
        else                ; ic0 = ic                              ; icent = 1 + sum(neq(:-ic-1))
        endif
        cmtin => cmtu(:,:,iband,kindx(ikpt),ispin)
        call waveftrafo_mt_io(cmtout,cmtin(1,icent),dim,maxlmindx,ic,symkpt(ikpt),kptp(ikpt))
        nullify(cmtin)
        if     (ic==0) then ; cmtout(:,:ncent)    = phase(ikpt) * cmtout(:,:ncent) ! undo division by phase factor
        else if(ic> 0) then ; cmtout(:,1     )    = phase(ikpt) * cmtout(:,1)
        else                ; cmtout(:,:neq(-ic)) = phase(ikpt) * cmtout(:,:neq(-ic))
        endif
      else
        i     = 0
        icent = 0
        do itype = 1,ntype
          do ieq = 1,neq(itype)
            icent       = icent + 1 ; if((ic>0.and.icent/=ic).or.(ic<0.and.itype/=-ic)) cycle
            i           = i + 1
            cmtout(:,i) = cmtu(:dim,icent,iband,kindx(ikpt),ispin)
          enddo
        enddo
      endif

      end

c     ------------------

      subroutine wavefunction_pw_w(cpwout,dim,iband,ikpt,ispin)

      use global

      use, intrinsic :: iso_fortran_env
      implicit none
      integer,    intent(in)  :: ikpt,iband,ispin,dim
      complex_dp, intent(out) :: cpwout(dim)
      complex_dp, pointer_cnt :: cpwin(:)

      if(ispin<1.or.ispin>nspin2) Bug('ispin out of range')
      if(l_soc) Error('Not implemented for SOC.')

      if(storeibz.and.kptp(ikpt)/=ikpt) then
        cpwin => cpwu(:,iband,kindx(ikpt),ispin)
        call waveftrafo_pw_io WSUFFIX (cpwout,cpwin,symkpt(ikpt),kptp(ikpt),ispin)
        nullify(cpwin)
      else
        cpwout = cpwu(:dim,iband,kindx(ikpt),ispin)
      endif

      end

c --------------------------------------------------------------------------------

# endif

c --------------------------------------------------------------------------------

# undef WSUB
# undef WSUFFIX
# ifdef INVW
#   undef INVW
#   define INV
# endif
# undef ADJUST_KINDX_IF_NEEDED
# undef UNDO_ADJUST_KINDX
# ifdef ISPIN1
#   undef ISPIN1
#   undef ISPIN2
# endif
