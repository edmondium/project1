c--------------------------------------------------------------------------------
c Copyright (c) 2021 Peter Grünberg Institut, Forschungszentrum Jülich, Germany
c This file is part of SPEX and available as free software under the conditions
c of the MIT license as expressed in the LICENSE file in more detail.
c--------------------------------------------------------------------------------

c     Set of routines that calculate the contribution of the Coulomb divergence to a BZ integration.
c
c     gamma_divergence   : divergence   = vol/(8pi^3) INT   1/k**2  d^3k - 1/nkpt SUM(k/=0)   1/k**2  ; divergence declared in global module
c     gamma_divergence_h : divergence_h = vol/(8pi^3) INT 4pi/(kHk) d^3k - 1/nkpt SUM(k/=0) 4pi/(kHk) ; H = packed hermitian 3x3 matrix
c     gamma_divergence_c : divergence_h = vol/(8pi^3) INT 4pi/(kHk) d^3k - 1/nkpt SUM(k/=0) 4pi/(kHk) ; H = complex 3x3 matrix
c
c     Input
c       hmat (H) : only symmetrized matrix relevant
c       hinv     : angular average of 4pi/(kHk), i.e., 4*pi*trace(hh)/3 with hh from invert_angular(hh)
c
c     Note: In coulomb_matrix.f, there is a another similar routine (gamma_divergence2) which takes into account an additional exponential factor.

# include "cppmacro.h"

      subroutine gamma_divergence(lwrite)
      use global
# ifdef MPI
      use Mwrapper, only: Msum
# endif
      use, intrinsic :: iso_fortran_env
      implicit none
      logical, intent(in) :: lwrite
      integer             :: ix,iy,iz,n
      logical             :: found
      real_dp, parameter  :: expo = 5d-3
      real_dp             :: div,rrad,drad,rrad2,kv1(3),kv2(3),kv3(3),k(3),k2,k1,fsmooth
      real                :: cputime
      complex_dp          :: cerf
      call cpu_time(cputime)
      Rif(lwrite) write(6,'(A'NoA) 'Calculate periodic function...'
      divergence = 0
      rrad       = sqrt(-log(expo)/expo)
      drad       = (rvol/nkpt)**(1d0/3)
      rrad2      = (rrad+drad)**2
      kv1        = rlat(:,1)/nkpt3(1)
      kv2        = rlat(:,2)/nkpt3(2)
      kv3        = rlat(:,3)/nkpt3(3)
      n          = 1 Mpi(+Mrank)
      found      = .true.
      do while(found)
        found = .false.
        div   = 0
        do ix = -n,n
          do iy = -(n-abs(ix)),n-abs(ix)
            iz   = n - abs(ix) - abs(iy)
            k(1) = ix*kv1(1) + iy*kv2(1) + iz*kv3(1)
            k(2) = ix*kv1(2) + iy*kv2(2) + iz*kv3(2)
            k(3) = ix*kv1(3) + iy*kv2(3) + iz*kv3(3)
            k2   = k(1)**2   + k(2)**2   + k(3)**2
            if(k2<rrad2) then
              k1 = sqrt(k2)
              k1 = k1 - rrad
              if(iz/=0) then ; fsmooth = 2
              else           ; fsmooth = 1
              endif
              if(abs(k1)<drad) fsmooth = fsmooth * ( (k1/drad)**3 - 3*(k1/drad) + 2 ) / 4
              found = .true.
              div   = div - exp(-expo*k2)/k2 / nkpt * fsmooth
            endif
          enddo
        enddo
        divergence = divergence + div
        n          = n + ifMpi(Msize,1)
      enddo
      Mpi( call Msum(divergence) )
      divergence = divergence + vol / (4*pi**2) * sqrt(pi/expo) * cerf((1d0,0d0)*sqrt(expo)*rrad)
      Rif(lwrite) call cpu_done(cputime)
      end

c     ----------------------------

      subroutine gamma_divergence_h(divergence_h,hmat,hinv)
      use global
      use wrapper, only: unpackmat
      Mpi2( use Mwrapper, only: Msum )
      use, intrinsic :: iso_fortran_env
      implicit none
      real_dp,     intent(out) :: divergence_h
      real_dp,     intent(in)  :: hinv
      MCOMPLEX_dp, intent(in)  :: hmat(6)
      real_dp                  :: hmatr(3,3)
      integer                  :: ix,iy,iz,n,opt
      logical                  :: found
      real_dp, parameter       :: expo = 5d-3
      real_dp                  :: div,rrad,drad,rrad2,kv1(3),kv2(3),kv3(3),k(3),k2,k1,khk,fsmooth
      complex_dp               :: cerf

      if(abs(hmat(1))+abs(hmat(3))+abs(hmat(6))>huge(0d0)) then
        if(hinv/=0) Bug('hinv expected to be zero.')
        divergence_h = 0
        return
      endif

      if        (abs(hmat(2))+abs(hmat(4))+abs(hmat(5))>1d-10) then ; opt = 0 ; hmatr = unpackmat(hmat) ! general matrix (only real part hmatr is relevant)
      else if(abs(hmat(1)-hmat(3))+abs(hmat(1)-hmat(6))>1d-10) then ; opt = 1                           ! diagonal matrix
      else ! scalar
        if(divergence==0) call gamma_divergence(.false.)
        divergence_h = divergence * hinv
        return
      endif

      divergence_h = 0
      rrad         = sqrt(-log(expo)/expo)
      drad         = (rvol/nkpt)**(1d0/3)
      rrad2        = (rrad+drad)**2
      kv1          = rlat(:,1)/nkpt3(1)
      kv2          = rlat(:,2)/nkpt3(2)
      kv3          = rlat(:,3)/nkpt3(3)
      n            = 1 Mpi(+Mrank)
      found        = .true.
      do while(found)
        found = .false.
        div   = 0
        do ix = -n,n
          do iy = -(n-abs(ix)),n-abs(ix)
            iz   = n - abs(ix) - abs(iy)
            k(1) = ix*kv1(1) + iy*kv2(1) + iz*kv3(1)
            k(2) = ix*kv1(2) + iy*kv2(2) + iz*kv3(2)
            k(3) = ix*kv1(3) + iy*kv2(3) + iz*kv3(3)
            k2   = k(1)**2   + k(2)**2   + k(3)**2
            if(k2<rrad2) then
              k1 = sqrt(k2)
              k1 = k1 - rrad
              if     (opt==0) then ; khk = dot_product( k , matmul(hmatr,k) )
              else if(opt==1) then ; khk = hmat(1)*k(1)**2 + hmat(3)*k(2)**2 + hmat(6)*k(3)**2
              endif
              if(iz/=0) then ;  fsmooth = 2
              else           ;  fsmooth = 1
              endif
              if(abs(k1)<drad) fsmooth = fsmooth * ( (k1/drad)**3 - 3*(k1/drad) + 2 ) / 4
              found = .true.
              div   = div - fsmooth * exp(-expo*k2) / khk
            endif
          enddo
        enddo
        divergence_h = divergence_h + div / nkpt
        n            = n + ifMpi(Msize,1)
      enddo
      Mpi( call Msum(divergence_h) )

      divergence_h = 4*pi * divergence_h + vol / (4*pi**2) * sqrt(pi/expo) * cerf((1d0,0d0)*sqrt(expo)*rrad) * hinv

      end

c     ----------------------------

      subroutine gamma_divergence_c(divergence_c,hmat,hinv)
      use global
      use wrapper, only: identity
      Mpi2( use Mwrapper, only: Msum )
      use, intrinsic :: iso_fortran_env
      implicit none
      complex_dp, intent(out) :: divergence_c
      complex_dp, intent(in)  :: hinv
      complex_dp, intent(in)  :: hmat(3,3)
      complex_dp              :: div,khk
      integer                 :: ix,iy,iz,n,opt
      logical                 :: found
      real_dp, parameter      :: expo = 5d-3
      real_dp                 :: rrad,drad,rrad2,kv1(3),kv2(3),kv3(3),k(3),k2,k1,fsmooth
      complex_dp              :: cerf

      if(abs(hmat(1,1))+abs(hmat(2,2))+abs(hmat(3,3))>huge(0d0)) then
        if(hinv/=0) Bug('hinv expected to be zero.')
        divergence_c = 0
        return
      endif

      if(sum(abs(hmat-hmat(1,1)*identity(3)))<1d-10) then ! scalar
        if(divergence==0) call gamma_divergence(.false.)
        divergence_c = divergence * hinv
        return
      else if(sum(abs(hmat))-abs(hmat(1,1))-abs(hmat(2,2))-abs(hmat(3,3))<1d-10) then ; opt = 1 ! diagonal matrix
      else                                                                            ; opt = 0 ! general matrix
      endif

      divergence_c = 0
      rrad         = sqrt(-log(expo)/expo)
      drad         = (rvol/nkpt)**(1d0/3)
      rrad2        = (rrad+drad)**2
      kv1          = rlat(:,1)/nkpt3(1)
      kv2          = rlat(:,2)/nkpt3(2)
      kv3          = rlat(:,3)/nkpt3(3)
      n            = 1 Mpi(+Mrank)
      found        = .true.
      do while(found)
        found = .false.
        div   = 0
        do ix = -n,n
          do iy = -(n-abs(ix)),n-abs(ix)
            iz   = n - abs(ix) - abs(iy)
            k(1) = ix*kv1(1) + iy*kv2(1) + iz*kv3(1)
            k(2) = ix*kv1(2) + iy*kv2(2) + iz*kv3(2)
            k(3) = ix*kv1(3) + iy*kv2(3) + iz*kv3(3)
            k2   = k(1)**2   + k(2)**2   + k(3)**2
            if(k2<rrad2) then
              k1 = sqrt(k2)
              k1 = k1 - rrad
              if     (opt==0) then ; khk = dot_product( k , matmul(hmat,k) )
              else if(opt==1) then ; khk = hmat(1,1)*k(1)**2 + hmat(2,2)*k(2)**2 + hmat(3,3)*k(3)**2
              endif
              if(iz/=0) then ; fsmooth = 2
              else           ; fsmooth = 1
              endif
              if(abs(k1)<drad) fsmooth = fsmooth * ( (k1/drad)**3 - 3*(k1/drad) + 2 ) / 4
              found = .true.
              div   = div - fsmooth * exp(-expo*k2) / khk
            endif
          enddo
        enddo
        divergence_c = divergence_c + div / nkpt
        n            = n + ifMpi(Msize,1)
      enddo
      Mpi( call Msum(divergence_c) )

      divergence_c = 4*pi * divergence_c + vol / (4*pi**2) * sqrt(pi/expo) * cerf((1d0,0d0)*sqrt(expo)*rrad) * hinv

      end

c ---------------------

c
c Returns factor describing the logarithmic divergence anomaly in metallic HF bands and also in the GW self-energy.
c
c This divergence originates from the long-range behavior of the interaction, 1/k**2 for k->0 in Fourier space.
c The calculation of the HF exchange potential (also of the GW self-energy for energies different from efermi)
c requires 1/k**2 to be integrated over the reciprocal space. The largest contribution comes from the region around
c k=0 with k being the momentum of the interaction. This corresponds in the HF exchange potential (the GW self-energy)
c for a state to the region around the kpoint of that state. If this region is cut by the Fermi surface, the logdiv
c factor is different from zero or one.
c
c The integration area can be pictured as a sphere of radius krad with the kpoint of the state in the middle
c (krad = radius of sphere with volume rvol/nkpt). For the GW self-energy, the interaction is not 1/k**2 but
c 1/(k**T*head*k). However, this cannot be integrated analytically. Therefore, we approximate this by
c k**T*hinv*k (hinv obtained from invert_angular). The integration gives
c
c 4*pi*k * trace(hinv)/3 * { 3/(y+3) * [ L(x) + y/3 * A(x) ] }
c with    x = kdist / krad          (kdist = distance to Fermi surface, -1<x<1)
c         y = 2*h33 / (h11+h22) - 1 (anisotropy factor y>-1, matrix h is unitary-transformed from hinv, see code)
c      L(x) = x*log(|x|) - x + 1
c      A(x) = x**3/4 - 3*x/4 + 1/2
c The returned logdiv factor is ONLY the expression in the curly brackets! (In the calling routine, 4*pi*k can be
c replaced by "divergence", which is also linear in k and shows a better k convergence.)
c
c For hinv=identity(3) (e.g., HF), this equation simplifies to 4*pi*k * L(x), then logdiv = L(x).
c Currently, the anisotropy correction is NOT used.
c
c energy       : band energy of state whose HF exchange (GW self-energy) matrix element is calculated
c                measured from efermi, i.e., ene-efermi.      
c moment(:3)   : band energy gradient of state (=expectation value of momentum operator <phi|-i\nabla|phi>).
c hinv(:3,:3)  : returned from invert_angular for head tensor of interaction (see above).
c                [hinv=0 is short for hinv being a multiple of identity(3).]
c metal = true : The integrated area may be cut by the Fermi surface (only energies below the Fermi energy contribute),
c                which is assumed to be a plane at a distance kdist=energy/amoment from sphere center (the kpoint of
c                the state), where amoment is the absolute value of the band gradient (at the state): |<phi|-i\nabla|phi>|**2
c
c If the Fermi surface does not cut the sphere (or if metal = false), then trivially
c   logdiv = 1 if energy < 0
c          = 0 if energy > 0
c Otherwise,
c   logdiv = 3/(y+3) * [ L(x) + y/3 * A(x) ]  (0<logdiv<1)
c
c This gives the logarithmic divergence anomaly of metallic HF bands [see definition of L(x)]. Note that the HF logarithmic
c divergence calculated for jellium follows a slightly different (but similar) function. The logarithmic divergence also
c affects the GW self-energy (and spectral function) but not the quasiparticle bands.
c
      function logdiv_c(energy,moment,hinv)
      use global, only: rvol,nkpt,pi,metal
      use util
      use, intrinsic :: iso_fortran_env
      implicit none
      complex_dp             :: logdiv_c
      real_dp,    intent(in) :: energy,moment(3)
      complex_dp, intent(in) :: hinv(3,3)
      complex_dp             :: h(3,3)
      real_dp                :: trafo(3,3),krad,kdist,x,y,amoment
      if(energy<=0) then ; logdiv_c = 1
      else               ; logdiv_c = 0
      endif
      amoment  = sqrt(sum(moment**2))          ; if(.not.metal.or.amoment<1d-6) return
      krad     = (3*rvol/(4*pi*nkpt))**(1d0/3)
      kdist    = energy / amoment              ; if(kdist==0) then ; logdiv_c = 0.5d0 ; return ; endif
      x        = kdist / krad                  ; if(abs(x)>=1) return
      logdiv_c = ( x*log(abs(x)) - x + 1 ) / 2 ; if(all(hinv==0)) return
      ! transform hinv so that moment points along z and add anisotropy contribution
      trafo(:,3) = moment / amoment
      if(any(abs(trafo(:2,3))>0.5d0)) then ; trafo(:,1) = [ trafo(2,3),-trafo(1,3),0d0]
      else                                 ; trafo(:,1) = [ trafo(3,3),0d0,-trafo(1,3)]
      endif
      trafo(:,1) = trafo(:,1) / sqrt(sum(trafo(:,1)**2))
      call vprod(trafo(:,2),trafo(:,3),trafo(:,1))
      h        = matmul(trafo,matmul(hinv,transpose(trafo))) ! rotated hinv
      y        = 2 * h(3,3) / (h(1,1) + h(2,2)) - 1          ! anisotropy factor
      logdiv_c = 3/(3+y) * ( logdiv_c + y/3 * ( x**3/4 - 3*x/4 + 1d0/2 ) )
      end
      
c Special functions. They call logdiv_c.
c (1) Real version      
      function logdiv_r(energy,moment,hinv)
      use, intrinsic :: iso_fortran_env
      implicit none
      real_dp             :: logdiv_r
      real_dp, intent(in) :: energy,moment(3),hinv(3,3)
      complex_dp          :: logdiv_c
      logdiv_r = logdiv_c(energy,moment,(1d0,0d0)*hinv)
      end
c (2) Isotropic version
      function logdiv(energy,moment)
      use, intrinsic :: iso_fortran_env
      implicit none
      real_dp             :: logdiv
      real_dp, intent(in) :: energy,moment(3)
      complex_dp          :: logdiv_c
      logdiv = logdiv_c(energy,moment,(1d0,0d0)*[0,0,0,0,0,0,0,0,0])
      end
      
