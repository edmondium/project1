c--------------------------------------------------------------------------------
c Copyright (c) 2021 Peter Grünberg Institut, Forschungszentrum Jülich, Germany
c This file is part of SPEX and available as free software under the conditions
c of the MIT license as expressed in the LICENSE file in more detail.
c--------------------------------------------------------------------------------

c Updates the GT self-energy (only diagonal elements) for the k-point ikpt.
c
c
c Generically, GT(w) = 1/(2pi) INT G(w+w') T(-w') dw'. (Note negative argument of T!)
c      
c
c Analytic continuation (oselfc==1)
c
c      s     s          s       1                 inf                ss'         s     s'*      s*    s'               s'    -1
c < phi  | GT (iw) | phi   > = --- SUM(k) SUM(m') INT dw' SUM  tmatrix (k,-iw') U     U        U     U        (iw+iw'-E     )
c      mq               mq     2pi               -inf     1234       1234        qm,1  q-km',2  qm,3  q-km',4          q-km'
c
c Contour integration (oselfc>=2)
c
c      s     s         s       1                 inf                ss'         s     s'*      s*    s'              s'    -1
c < phi  | GT (w) | phi   > = --- SUM(k) SUM(m') INT dw' SUM  tmatrix (k,-iw') U     U        U     U        (w+iw'-E     )
c      mq              mq     2pi               -inf     1234       1234        qm,1  q-km',2  qm,3  q-km',4         q-km'
c
c                                        u/o                  s'                        s'     s     s'*      s*    s'
c                             -/+ SUM(k) SUM(m') theta[+/-(w-E    )] SUM  tmatrix (k,w-E    ) U     U        U     U
c                                                             q-km'  1234               q-km'  qm,1  q-km',2  qm,3  q-km',4
c for w>=0 (w<0). (Energies relative to Fermi energy.)
c
c Using the tetrahedron weight function P(e) (as returned by from tetrahedron_nonlin), the residue part can be reformulated as
c
c     u/o                                                    u/o
c -/+ SUM(m') SUM(k) INT(0,w)  P   (e) Tmatrix(k,w-e) de = - SUM(m') SUM(k) INT(0,w) P   (e) Tmatrix(k,w-e),
c                    INT(w,0)   q-km'                                                 q-km'
c
c where Tmatrix abbreviates SUM(1234) tmatrix ... (omitting several indices for simplicity).
c In the "normal" case, one has instead a weighted sum with tetrahedron weights that takes into account the u/o restriction and the theta function
c
c -/+ SUM(m') SUM(k) p      Tmatrix(k,w-e).
c                     q-km'
c
c k    : ikpt
c s    : job1%spin
c s'   : opposite to s
c ss'  : spint (if s's is required, tmatrix is transformed)
c w',w : freq(:nfreq) (arrays.f)
c lsym : If set, all equivalents to ikpt are summed too by using irrep_wan (IRREP).
c
c If symmetry is used (lsym), the k-point sum is restricted to the EIBZ(q), which can be disabled by
c compiling with "-Dswitch_off_symmetry".
c

# include "cppmacro.h"

      subroutine selfenergy_wannier(job1,ikpt,tmatrix,spint,lsym)

      use global
      use arrays
      use freq_integral
      Mpi( use Mwrapper )

      use, intrinsic :: iso_fortran_env
      implicit none
      type(jobtype), intent(in) :: job1
      integer,       intent(in) :: ikpt,spint
      logical,       intent(in) :: lsym
      complex_dp,    intent(in) :: tmatrix(nwan,nwan,nwan,nwan,*)
      complex_dp                :: t(nwan,nwan,nwan,nwan),uwan1(nwan)
      complex_dp                :: tuuuu(nfreq+nfreqc),tuu(nwan,nwan,nfreq+nfreqc)
      complex_dp                :: cpade(nfreq+1),pole((nfreq+1)/2*smooth(2)),resid((nfreq+1)/2*smooth(2)),cdum,carg
      real_dp                   :: enek,eneq,weight,freqr1(nfreqr),coeff(nfreqc),frqc(nfreqc),enediff
      logical                   :: ltrafo,ldef(nfreq+nfreqc),ldum
      integer                   :: npole1
      integer                   :: kpt1(nkpt),nkpt1,nkpts(nkpt),sym1(nsym),nsym1,deg0,deg1,ideg
      integer                   :: ib,ibandq,ikptq,ispinq,ibandk,ikptk,ikptp,ispink,ifreq,ifreqc,ikpt1
      integer                   :: n1,n2,k,sgn
      integer                   :: kptdiff,kptsum
      integer                   :: nfrq,i
      real_dp                   :: frq(27),pfrq(0:3,26),emin,emax,f0,f1,f2,ene0
      real_dp                   :: tetrahedron_weight
      real                      :: time
      complex_dp                :: pade_func,pade_inter,residues,residues_pade
      integer, save             :: ichk=-1
      Mpi2( integer             :: Merr,mm )

      if(ichk/=-1) read(*,*) ichk

      nfrq   = 1 ! Default if tetrahedron_nonlin not used
      frq(1) = 0 !

      if(lsym.and..not.allocated(irrep_wan)) Bug('Cannot symmetrize without irrep_wan.')
      if(lsym.and.ikpt>nkpti)                Bug('input kpoint not in IBZ.')
      if(nfreqc>0) then
        if(any(imag(freqc)/=0))                               Bug('Frequency for T matrix has imaginary part.')
        if(any([(real(freqc(i+1)-freqc(i))<0,i=1,nfreqc-1)])) Bug('Frequency mesh for T matrix not size-ordered.')
      endif
      if(ikpt==1) then
        if(option==102) then ; RWarn('Gamma point included (by option 102).')
        else                 ; RInfo('Gamma point excluded. Override with -o=102') ; return
        endif
      endif

      Rwrite(6,'(/A'NoA) 'Add contribution to self-energy... '
      call cpu_time(time)

      Nfence(selfc)
      Mpi( mm=-1 )

      if(oselfc==1.or.oselfc==4) freqr1 = freqr

c
c     Loop over bands
      do ib = 1,size(job1%band)
        ibandq = job1%band(ib)
        ikptq  = job1%kpt (ib)
        ispinq = job1%spin(ib)
        ispink = 3 - ispinq
        ltrafo = spint-ispinq/=2

        ! For residue part: Invert freqc (->frqc), unless ltrafo (in which case transform_spinpair implies w->-w); tuuuu(:) will later be ordered accordingly
        if(oselfc/=1.and.nfreqc>0) then
          if(ltrafo) then ; frqc              =  freqc
          else            ; frqc(nfreqc:1:-1) = -freqc
          endif
        endif

        deg0 = ibandq
        deg1 = ibandq
        if(lsym) then ! kpt1 = members of EIBZ(ikptq) that are equivalent to ikpt (all equivalents if "-Dswitch_off_symmetry" is set)
          call getkpt1(kpt1,nkpt1,nkpts,sym1,nsym1,ikptq,ikpt,.false.)
          if(nkpt1<count(kptp==ikpt)) call getdeg(deg0,deg1,ibandq,ikptq,ispinq)
        else          ! kpt1 = only ikpt
          kpt1  = ikpt
          nkpt1 = 1
          nkpts = 1
        endif

        if(oselfc==2.or.oselfc==3) then !; call teststop('oselfc==2,3 for GT')
          eneq = ene(ibandq,ikptq,ispinq) - efermi
          if(eneq<=0) then ; freqr1 = eneq - freqr
          else             ; freqr1 = eneq + freqr
          endif
        endif

c
c       Loop over symmetry equivalent k points (if lsym)
        do k = 1,nkpt1 ; McycleP(mm)
          ikpt1 = kpt1(k)

          if(ltrafo) then ; ikptk = kptsum( ikptq,ikpt1) ! we will use the spinpair-transformed T, which implies k -> -k (and w -> -w)
          else            ; ikptk = kptdiff(ikptq,ikpt1)
          endif
          ikptp = kptp(ikptk) ; if(ikptp>nkpt) ikptp = ikptp - nkpt + nkpti

          ldef = .false.
          tuu  = 0

c
c         Loop over Green-function poles
          do ibandk = wanbandi,wanbandf

c            if(ikpt==1.and.same_eigenspace(ibandq,ibandk,ikptk,ispink).and.option/=102) cycle ! Exclude Goldstone delta peak. Not successful. Now k=0 fully excluded.

            enek = ene(ibandk,ikptk,ispink) - efermi
            ene0 = enek ! origin of energy scale (e.g., for polynomials)

            ! Tetrahedron weight function (FREQINT NONLIN)
            if(any(freqint==[0,3])) then
              call tetrahedron_nonlin(nfrq,frq,pfrq,ene0,ibandk,ikptk,ispink,ichk==1.or.ichk==3)
              if(ichk>=2) call checktetra(nfrq,frq,pfrq)
            endif

c           Multiply T*U*U*U*U -> tuuuu
c           Auxiliary array tuu=T*U*U is calculated if needed (ldum) unless already calculated (ldef)
            tuuuu = 0
            ! determine ifreq for which tuuuu has to be calculated (ldum=.true.) ...
            if(nfreqc>0) then 
              f1 = minval(freqr1) - (ene0+frq(nfrq)) ! determine real-frequency range [f1,f2] needed for tuuuu(:)
              f2 = maxval(freqr1) - (ene0+frq(1)) ; if(ltrafo) then ; f0 = f2 ; f2 = -f1 ; f1 = -f0 ; endif
            endif
            do ifreq = 1,nfreq+nfreqc
              ldum = ifreq<=nfreq                         ! ... has to be calculated for all imaginary frequencies anyway
              if(.not.ldum) then ; ifreqc = ifreq - nfreq ! ... but not for all real frequencies
                ldum = .false.
                if(ifreqc>1)      ldum =           real(freqc(ifreqc-1))<f2 .and. real(freqc(ifreqc))  >f1 !  left interval overlaps with [f1,f2]
                if(ifreqc<nfreqc) ldum = ldum .or. real(freqc(ifreqc))  <f2 .and. real(freqc(ifreqc+1))>f1 ! right interval overlaps with [f1,f2]
              endif
c# warning ldum always true              
c              ldum = .true.
              if(ldum) then
                if(.not.ldef(ifreq)) then
                  t = tmatrix(:,:,:,:,ifreq)
                  if(ltrafo) call transform_spinpair(t)
                  do ideg = deg0,deg1
                    if(ikpt1/=ikpt) then ; uwan1 = matmul ( uwan(ideg,:,ikptq,ispinq) , irrep_wan(:,:,symkpt(ikpt1),ispinq) )
                    else                 ; uwan1 =          uwan(ideg,:,ikptq,ispinq)
                    endif
                    do n2 = 1,nwan ; do n1 = 1,nwan
                      tuu(n1,n2,ifreq) = tuu(n1,n2,ifreq) + dot_product( uwan1 , matmul ( transpose(t(:,n1,:,n2)), uwan1 ) )
                    enddo ; enddo
                  enddo
                  if(ikpt1/=ikpt) call transform_sym(tuu(:,:,ifreq),ikpt1)
                  ldef(ifreq) = .true.
                endif
                tuuuu(ifreq) = tuuuu(ifreq) + dot_product( uwan(ibandk,:,ikptk,ispink), matmul ( tuu(:,:,ifreq),
     &                                                     uwan(ibandk,:,ikptk,ispink) ) )
              endif
            enddo
            tuuuu = tuuuu * nkpts(k) / (deg1-deg0+1)

            ! We need tuuuu(w) = Tmatrix(-w) because Tmatrix has negative argument in convolution. (Already fulfilled if ltrafo, otherwise we transform accordingly.)
            if(.not.ltrafo) then
              tuuuu(:nfreq) = conjg(tuuuu(:nfreq)) ! using T_1234(q,-iw') = conjg[T_3412(q,iw')], which is valid for imag. freq. axis
              if(nfreqc>0) tuuuu(nfreq+1:) = tuuuu(nfreq+nfreqc:nfreq+1:-1)
            endif

c# warning Tmatrix write
c            rewind(100)
c            write(100,'(3F16.8)') (frqc(i),tuuuu(nfreq+i),i=1,nfreqc)
c            write(*,'(A,5I4,L4,F15.6)') 'current (k,k,q+k,m,n):',ikpt,ikpt1,ikptk,ibandq,ibandk,ltrafo,
c     &        maxval(abs(real(tuuuu(nfreq+1:nfreq+nfreqc))))
c            read(*,*)

c
c           Convolute G*tuuuu
            if(freqint/=1.or.nfreqc==0) then
              if(mod(nfreq,2)/=0) then
                RWarn('Number of frequencies odd. Check constraints of pade_init.')
                call pade_init(cpade,img*freq(nfreq:1:-1),tuuuu(nfreq:1:-1),nfreq,1)
                call pade_poles(pole,resid,npole1,img*freq(nfreq:1:-1),tuuuu(nfreq:1:-1),cpade,nfreq,smooth(2),0,.false.)
              else
                call pade_init(cpade,img*freq,tuuuu,nfreq,0)
                call pade_poles(pole,resid,npole1,img*freq,tuuuu,cpade,nfreq,smooth(2),0,.false.)
              endif              
            endif

            do ifreq = 1,size(selfc,2)

              ! Integral

              if(any(freqint==[0,3])) then ! FREQINT NONLIN
                if(npole1>0) then
                  if(oselfc==1) then ; carg = freq(ifreq) * img
                  else               ; carg = freqr1(ifreq)
                  endif                
                  cdum = freqintegral_poles_nonlin(pole,resid,npole1,2,carg,0,ene0,nfrq,frq,pfrq)
                  Nacc1_c( selfc,(ib,ifreq) , cdum )
                endif
              else                         ! FREQINT SPLINE/PADE
                weight = 0
                if(oselfc==1) then ; if(ibandk<=bando) weight = wintgr(ikptk,ibandk,ispink)
                else               ;                   weight = tetrahedron_weight(ibandk,ikptk,ispink,freqr1(ifreq))
                endif
                do sgn = -1,1,2 ! -1 : occ, 1 : unocc
                  if(weight>1d-10) then
                    if(oselfc==1) then ; carg = freq(ifreq) + sgn * img * abs(enek) ! abs(..) avoids sign change of imag. part
                    else               ; carg =               sgn * img * abs(enek-freqr1(ifreq))
                    endif
                    if(freqint/=1) then ; cdum = freqintegral_poles(pole,resid,npole1,2,carg,0)
                    else                ; cdum = freqintegral(tuuuu,freq,nfreq,2,carg,0)
                    endif
                    Nacc1_c( selfc,(ib,ifreq) , cdum * weight )
                  endif
                  weight = 1d0/nkpt - weight
                enddo
              endif

              ! Residues

              if(oselfc/=1) then

                if(any(freqint==[0,3])) then
                  emin = min(0d0,freqr1(ifreq)) ! Bounds for integration over states (relative to Fermi energy)
                  emax = max(0d0,freqr1(ifreq)) ! Since tuuuu(e-w) is the argument, tuuuu is integrated from -emax to -emin.
                  if(emax-emin>1d-8.and.frq(1)+ene0<emax.and.frq(nfrq)+ene0>emin) then
                    if(nfreqc>0) then
                      cdum = residues(emin,emax,ene0,nfrq,frq,pfrq,nfreqc,frqc+freqr1(ifreq),tuuuu(nfreq+1:),.false.) ! INT(0,w) P(e) tuuuu(e-w) de (shifting e mesh up by w gives  INT(0,w) P(e) tuuuu(e) de  as required)
                    else
                      cdum = residues_pade(emin,emax,ene0,nfrq,frq,pfrq,npole1,pole+freqr1(ifreq),resid) ! shifting all poles up by w gives: INT(0,w) P(e) tuuuu(e) de (tuuuu(e) = Pade approximant)
                      if(emin+emax>0) cdum = conjg(cdum) ! tuuuu is integrated from -emax to -emin; so, if emin and emax >= 0, then the argument of tuuuu is always negative, hence the complex conjugation
                    endif
                    if(emin==0d0) cdum = -cdum
                    Nacc1_c( selfc,(ib,ifreq) , cdum )
                  endif
                else
                  weight = -weight ; if(ibandk<=bando) weight = weight + wintgr(ikptk,ibandk,ispink)
                  if(abs(weight)>1d-10) then
                    enediff = enek - freqr1(ifreq)
                    if(nfreqc==0) then
                      if(smooth(2)>1) then ; cdum = pade_inter((1d0,0d0)*enediff,pole,resid,npole1)    ! note that tuuuu(-iw')
                      else                 ; cdum = pade_func ((1d0,0d0)*enediff,img*freq,cpade,nfreq) !
                      endif
                      if(enediff<0) cdum = conjg(cdum)
                    else
                      call getcoeff(coeff,enediff,frqc,nfreqc) ! note that tuuuu(w) = Tmatrix(-w)
                      cdum = sum(coeff*tuuuu(nfreq+1:))
                    endif
                    Nacc1_c( selfc,(ib,ifreq) , cdum * weight )
                  endif
                endif

              endif

            enddo 
          enddo

        enddo ! ikpt1 loop

      enddo ! ib loop

      Mpi( call mpi_barrier(Mcomm,Merr) )
      Rcall cpu_done(time)

# ifdef MPI
      Nfence(selfc)
      Ocall Msum(selfc,0,comm=Ocomm)
      MnoR( ifO selfc = 0 )
      Nfence(selfc)
# endif

      contains

c ------------

c
c     Switches spins in T: Tud_1234(w,q) -> Tdu_4321(-w,-q) (and vice versa)
      subroutine transform_spinpair(mat)
      use, intrinsic :: iso_fortran_env
      implicit none
      complex_dp, intent(inout) :: mat(nwan,nwan,nwan,nwan)
      integer                   :: i,j
      do i = 1,nwan ; do j = 1,nwan
        mat(i,:,:,j) = transpose(mat(i,:,:,j))
      enddo ; enddo
      do i = 1,nwan ; do j = 1,nwan
        mat(:,i,j,:) = transpose(mat(:,i,j,:))
      enddo ; enddo
      end subroutine transform_spinpair

c ------------

c
c     Transforms tuu to symmetry-equivalent k-point ikpt
      subroutine transform_sym(tuu,ikpt)
      use, intrinsic :: iso_fortran_env
      implicit none
      complex_dp, intent(inout) :: tuu(nwan,nwan)
      integer,    intent(in)    :: ikpt
      integer                   :: isym
      isym = symkpt(ikpt)
      if(isym==1) return
      if(isym>nsymt) Error('time-reversal symmetry not implemented.')
      tuu = matmul ( conjg(irrep_wan(:,:,isym,ispink)) , matmul ( tuu , transpose(irrep_wan(:,:,isym,ispink)) ) )
      end subroutine transform_sym

# if 0
c
c     Transforms T to symmetry-equivalent k-point ikpt (backup)
      subroutine transform_sym(t,ikpt)
      use, intrinsic :: iso_fortran_env
      implicit none
      complex_dp, intent(inout) :: t(nwan,nwan,nwan,nwan)
      integer,    intent(in)    :: ikpt
      integer                   :: isym,n1,n2,n3
      isym = symkpt(ikpt) ; if(isym==1) return
      if(isym>nsymt) Error('time-reversal symmetry not implemented.')
      do n3 = 1,nwan
      do n1 = 1,nwan ; do n2 = 1,nwan ; t(:,n1,n2,n3) = matmul(      irrep_wan(:,:,isym,ispinq) ,t(:,n1,n2,n3)) ; enddo ; enddo
      enddo
      do n3 = 1,nwan
      do n1 = 1,nwan ; do n2 = 1,nwan ; t(n1,:,n2,n3) = matmul(conjg(irrep_wan(:,:,isym,ispink)),t(n1,:,n2,n3)) ; enddo ; enddo
      enddo
      do n3 = 1,nwan
      do n1 = 1,nwan ; do n2 = 1,nwan ; t(n1,n2,:,n3) = matmul(conjg(irrep_wan(:,:,isym,ispinq)),t(n1,n2,:,n3)) ; enddo ; enddo
      enddo
      do n3 = 1,nwan
      do n1 = 1,nwan ; do n2 = 1,nwan ; t(n1,n2,n3,:) = matmul(      irrep_wan(:,:,isym,ispink) ,t(n1,n2,n3,:)) ; enddo ; enddo
      enddo
      end subroutine transform_sym
# endif

c ------------

      end

c ------------

      subroutine checktetra(nfrq,frq,pfrq)
      use global, only: nkpt
      use file
      use, intrinsic :: iso_fortran_env
      implicit none
      integer, intent(in) :: nfrq
      real_dp, intent(in) :: frq(nfrq),pfrq(0:3,nfrq-1)
      real_dp             :: e,int
      integer             :: i,j,iunit
      iunit = fopen('checktetra',status='unknown')
      rewind iunit
      write(iunit,'(''#'',I3'NoA) nfrq
      do i = 1,nfrq
        write(iunit,'(F10.5'NoA) frq(i)
      enddo
      write(iunit,*)
      int = 0
      do i = 1,nfrq-1
        do j = 0,100
          e = (frq(i)*(100-j)+frq(i+1)*j)/100
          write(iunit,'(2F20.10)') e,
     &      pfrq(0,i)+e*(pfrq(1,i)+e*(pfrq(2,i)+e*pfrq(3,i)))
        enddo
        e   = frq(i+1)
        int = int + e*(pfrq(0,i)+e*(pfrq(1,i)/2+e*(pfrq(2,i)/3+e*pfrq(3,i)/4)))
        e   = frq(i)
        int = int - e*(pfrq(0,i)+e*(pfrq(1,i)/2+e*(pfrq(2,i)/3+e*pfrq(3,i)/4)))
      enddo
      call fclose(iunit)
      write(*,'(I4,F15.10)') nfrq,abs(int-(1d0*nkpt)**(-1))
      read(*,*)
      end 

c ------------
