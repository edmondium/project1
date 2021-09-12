c--------------------------------------------------------------------------------
c Copyright (c) 2021 Peter Grünberg Institut, Forschungszentrum Jülich, Germany
c This file is part of SPEX and available as free software under the conditions
c of the MIT license as expressed in the LICENSE file in more detail.
c--------------------------------------------------------------------------------

c This file contains routines that calculate weights for integrations over
c the reciprocal space using the tetrahedron method.
c
c   In this method "reciprocal cubes" defined by 8 k points as cube corners are
c   decomposed into six tetrahedra inside which the electron bands and the
c   Fermi surface is linearized.
c
c Symmetry considerations:
c   Unfortunately, the definition of the tetrahedra
c   is not unique. Any group of six tetrahedra out of the 24 tetrahedra
c   ( 1,2,3,6, 5,7,3,6, 1,5,3,6, 2,4,3,6, 4,8,3,6, 7,8,3,6,
c     5,6,2,7, 1,5,2,7, 1,3,2,7, 8,6,2,7, 4,3,2,7, 8,4,2,7,
c     2,6,1,8, 2,4,1,8, 3,4,1,8, 3,7,1,8, 5,7,1,8, 5,6,1,8,
c     2,6,4,5, 1,2,4,5, 1,3,4,5, 3,7,4,5, 7,8,4,5, 6,8,4,5 )
c   that fills up a given cube completely can be chosen. This non-uniqueness
c   breaks the symmetry of the system. As a result, quantities of symmetry-equivalent
c   k points might not be identical (but similar). However, even if we sum the
c   contribution of all 24 tetrahedra and divide by 4, there can still be
c   symmetry-equivalent k points (e.g. (0,0,1/4) and (1/4,1/4,1/4) in fcc)
c   where quantities are not identical. This comes from the definition of the
c   cubes which also break the symmetry.
c   (For testing, the Gauss method should be used. It is not as accurate but does
c   not suffer from these symmetry-related problems.)
c
c   Starting from version 05.00, the "optimal" tetrahedron decomposition is
c   chosen in get_tetra with the condition that the common diagonal is shortest.
c   To enforce the old behavior (a fixed decomposition), uncomment the following
c# define old_version
c
# ifdef CMP_0405adj
#   define old_version
# endif
c
c
# include "cppmacro.h"
c --------------------
c
c     Initializes the tetrahedron integration for occupied states,
c     i.e. calculates (for any f) the weights w(k,n,s) in the sum
c
c      2                                  2  occ               3
c     SUM SUM wintgr(k,n,s) * f(k,n,s) = SUM SUM INT f(k,n,s) d k .
c     s=1 k,n                            s=1  n
c
c     SUM wintgr(k,n,s)  gives the number of states up to the Fermi energy (efermi).
c     This can be used to determine the Fermi energy.
c
c     dim           : Leading dimension of wintgr_out (if dim>=nkpt, calculation for whole BZ, otherwise, calculation only for IBZ)
c     minb, maxb    : minimal and maximal band index
c     spin          : spin index (spin=0: both spins)
c     lshift=.true. : Weights calculated for shifted k set
      subroutine tetrahedron_init(wintgr_out,dim,minb,maxb,efermi,spin,lshift)

      use global, only: ene,nkpt3,kptp,pkpt,pkpt1,nkpt,nkpti,nkpti2,nspin1,obloechl
      use, intrinsic :: iso_fortran_env

      implicit none
# include "interface/cuttetra.inc"
      integer, intent(in)    :: minb,maxb,spin,dim
      real_dp, intent(in)    :: efermi
      real_dp, intent(out)   :: wintgr_out(dim,minb:maxb,*)
      logical, intent(in)    :: lshift
      integer                :: tetra(4,6),pnt(4),alladd(8)
      integer, pointer       :: pkpt0(:,:,:)
      integer                :: s,s1,n,itetra,k1,k2,k3,n1,i,kk,ntetra
      integer                :: spin1,spin2
      integer, allocatable   :: kcorn(:)
      real_dp                :: ecube(8),etetra(4)
      real_dp                :: wtetra(4,3),vtetra(3),rdum,atria

      if(minb>maxb) return

      if(lshift) then ; pkpt0 => pkpt1 ; if(dim<nkpti2) Bug('First dimension of wintgr too small.')
      else            ; pkpt0 => pkpt  ; if(dim<nkpti)  Bug('First dimension of wintgr too small.')
      endif

      call get_tetra(tetra,n,1)
      do i = 1,8
        alladd(i) = count(tetra==i)
      enddo

      allocate ( kcorn(8) )

      if(spin==0) then ; spin1 = 1    ; spin2 = nspin1
      else             ; spin1 = spin ; spin2 =  spin
      endif
      wintgr_out(:,:,:spin2-spin1+1) = 0

      do n = minb,maxb

        if(lshift) then ; if(minval(ene(n,nkpt+1:,spin1:spin2))>efermi) exit
        else            ; if(minval(ene(n,:nkpt,  spin1:spin2))>efermi) exit
        endif
        do s = spin1,spin2
          s1 = s - spin1 + 1
          do k3 = 0,nkpt3(3)-1     !
            do k2 = 0,nkpt3(2)-1   ! Loop over cubes in BZ
              do k1 = 0,nkpt3(1)-1 !

                kcorn(1) = kptp( pkpt0(k1  ,k2  ,k3  ) ) !
                kcorn(2) = kptp( pkpt0(k1+1,k2  ,k3  ) ) !
                kcorn(3) = kptp( pkpt0(k1  ,k2+1,k3  ) ) !
                kcorn(4) = kptp( pkpt0(k1+1,k2+1,k3  ) ) ! Definition of the cube corners
                kcorn(5) = kptp( pkpt0(k1  ,k2  ,k3+1) ) !
                kcorn(6) = kptp( pkpt0(k1+1,k2  ,k3+1) ) !
                kcorn(7) = kptp( pkpt0(k1  ,k2+1,k3+1) ) !
                kcorn(8) = kptp( pkpt0(k1+1,k2+1,k3+1) ) !_____
                ecube    = ene(n,kcorn,s)  ! Get energies at cube corners
                if(lshift) kcorn = kcorn - nkpt

                if(efermi>=maxval(ecube)) then                                        !
                  do i = 1,8                                                          !
                    wintgr_out(kcorn(i),n,s1) = wintgr_out(kcorn(i),n,s1) + alladd(i) ! treat the trivial case of
                  enddo                                                               ! efermi >= all energies
                  cycle                                                               !
                endif                                                                 !

                do itetra = 1,6

                  etetra = ecube(tetra(:,itetra)) ! Get energies at tetrahedron corners
                  if(all(etetra>=efermi)) cycle ! Nothing to be done if tetrahedron-corner energies larger than efermi

                  call cuttetra(wtetra,vtetra,ntetra,etetra,efermi)

                  do i = 1,4
                    kk                  = kcorn(tetra(i,itetra))
                    wintgr_out(kk,n,s1) = wintgr_out(kk,n,s1) + sum ( vtetra(:ntetra)*wtetra(i,:ntetra) )
                  enddo

                  ! Bloechl correction
                  if(obloechl) then
                    ! Order according to energy
                    call rorderp(pnt,etetra,4)
                    etetra = etetra(pnt)
                    ! Fermi surface inside tetrahedron (composed of one or two triangles)
                    if     (efermi>etetra(4)) then ; cycle
                    else if(efermi>etetra(3)) then
                      atria = 3 * (efermi-etetra(4))**2 / ( (etetra(4)-etetra(1)) * (etetra(4)-etetra(2)) * (etetra(4)-etetra(3)) )
                    else if(efermi>etetra(2)) then
                      atria = 3 * ( (efermi-etetra(2)) * (etetra(3)-efermi) / (etetra(3)-etetra(2))
     &                            + (efermi-etetra(1)) * (etetra(4)-efermi) / (etetra(4)-etetra(1)) ) /
     &                            ( (etetra(4)-etetra(2)) * (etetra(3)-etetra(1)) )
                    else
                      atria = 3 * (efermi-etetra(1))**2 / ( (etetra(2)-etetra(1)) * (etetra(3)-etetra(1)) * (etetra(4)-etetra(1)) )
                    endif
                    ! Add correction to weights
                    do i = 1,4
                      kk                  = kcorn(tetra(pnt(i),itetra))
                      wintgr_out(kk,n,s1) = wintgr_out(kk,n,s1) + atria / 10 * sum(etetra-etetra(i)) ! (another factor 1/4 below)
                    enddo
                  endif

                enddo

              enddo
            enddo
          enddo
        enddo
      enddo

      s1                  = spin2 - spin1 + 1
      wintgr_out(:,:,:s1) = wintgr_out(:,:,:s1) / 4 / 6 / nkpt ! take average over four tetrahedron corners,
                                                               ! six tetrahedra per cube,
                                                               ! and
                                                               ! nkpt cubes

      deallocate ( kcorn )

      if(lshift) then ; k2 = nkpti2
      else            ; k2 = nkpti
      endif
      do s = spin1,spin2
        s1 = s - spin1 + 1
        ! Average over degenerate eigenstates and divide by k multiplicity
        do k1 = 1,k2          
          if(lshift) then ; k3 = k1 + nkpt ; kk = count(kptp(nkpt+1:)==k3)
          else            ; k3 = k1        ; kk = count(kptp(:nkpt)  ==k3)
          endif
          n = minb
          do while(n<=maxb)
            n1 = n
            do while(n1<maxb)
              if(abs(ene(n1+1,k3,s)-ene(n,k3,s))>1d-8) exit
              n1 = n1 + 1
            enddo
            rdum                   = sum(wintgr_out(k1,n:n1,s1))
            wintgr_out(k1,n:n1,s1) = rdum / ((n1-n+1)*kk)
            n                      = n1 + 1
          enddo
        enddo
        ! Distribute symmetry-equivalent k points (if needed)
        if(dim>=nkpt) then
          do k1 = k2+1,nkpt
            if(lshift) then ; kk = kptp(k1+nkpt)-nkpt
            else            ; kk = kptp(k1)
            endif
            wintgr_out(k1,:,s1) = wintgr_out(kk,:,s1)
          enddo
        endif
      enddo      

      end

c -----------------

c Returns tetrahedron corners, 1-8 are the cube-corner indices.
c mode = 1 or 2 : Shortest diagonal is chosen.
c                 If more than one diagonal is shortest:
c      mode = 1 : Returns corners of only one "representative" set of six tetrahedra (ntetra=6)
c      mode = 2 : Returns corners of all tetrahedra (ntetra = 6, 12, 18, or 24)
c mode = 3      : Returns corners of all 24 tetrahedra (ntetra = 24)
c
      subroutine get_tetra(tetra,ntetra,mode)
      use global, only: rlat,nkpt3 MpiC(Mrank)
      use, intrinsic :: iso_fortran_env
      implicit none
      integer, intent(out) :: tetra(4,*),ntetra
      integer, intent(in)  :: mode
      real_dp              :: dmin,diag(4),k(3),rl(3,3)
      integer              :: idmin
# ifdef old_version
#   warning old_version defined
      logical, save :: first=.true.
      ntetra      = 6
      tetra(:,:6) = reshape( [ 1,2,3,6, 5,3,6,7, 1,5,3,6,   ! fixed tetrahedron decomposition
     &                         8,6,7,2, 4,7,2,3, 8,4,7,2 ], ! used until version 04.04
     &                       [ 4,6 ] )
      Rif(first) Info('old_version defined in tetrahedron.f.')
      first = .false.
      return
# endif
      rl(:,1) = rlat(:,1) / nkpt3
      rl(:,2) = rlat(:,2) / nkpt3
      rl(:,3) = rlat(:,3) / nkpt3
      k       = rl(:,1) + rl(:,3) - rl(:,2) ! diagonal 3-6
      diag(1) = sum(k*k)
      k       = rl(:,2) + rl(:,3) - rl(:,1) ! diagonal 2-7
      diag(2) = sum(k*k)
      k       = rl(:,1) + rl(:,2) + rl(:,3) ! diagonal 1-8
      diag(3) = sum(k*k)
      k       = rl(:,1) + rl(:,2) - rl(:,3) ! diagonal 4-5
      diag(4) = sum(k*k)
      idmin   = minloc(diag,1)
      dmin    = diag(idmin)
      ntetra  = 0
      if(mode==3.or.mode==2.and.abs(diag(1)-dmin)<1d-8.or.idmin==1) then
        tetra(:,ntetra+1:ntetra+6) = reshape ( [ 1,2,3,6, 5,7,3,6, 1,5,3,6, 2,4,3,6, 4,8,3,6, 7,8,3,6 ], [ 4,6 ] )
        ntetra                     = ntetra + 6
      endif
      if(mode==3.or.mode==2.and.abs(diag(2)-dmin)<1d-8.or.idmin==2) then
        tetra(:,ntetra+1:ntetra+6) = reshape ( [ 5,6,2,7, 1,5,2,7, 1,3,2,7, 8,6,2,7, 4,3,2,7, 8,4,2,7 ], [ 4,6 ] )
        ntetra                     = ntetra + 6
      endif
      if(mode==3.or.mode==2.and.abs(diag(3)-dmin)<1d-8.or.idmin==3) then
        tetra(:,ntetra+1:ntetra+6) = reshape ( [ 2,6,1,8, 2,4,1,8, 3,4,1,8, 3,7,1,8, 5,7,1,8, 5,6,1,8 ], [ 4,6 ] )
        ntetra                     = ntetra + 6
      endif
      if(mode==3.or.mode==2.and.abs(diag(4)-dmin)<1d-8.or.idmin==4) then
        tetra(:,ntetra+1:ntetra+6) = reshape ( [ 2,6,4,5, 1,2,4,5, 1,3,4,5, 3,7,4,5, 7,8,4,5, 6,8,4,5 ], [ 4,6 ] )
        ntetra                     = ntetra + 6
      endif
      end

c -----------------

c     Cuts the tetrahedron at the Fermi energy efermi.
c     The volume fraction below efermi is given by ntetra (1 or 3) tetrahedra defined by
c     wtetra(1:4,1:ntetra)    : weights of the four tetrahedron corners
c     vtetra(1:ntetra)        : volumes of the ntetra new tetrahedra
c     tetra(1:4,1:4,1:ntetra) : corner points of the new tetrahedra (optional)
c                               [ wtetra = sum(i) tetra(:,i,:) ]
c
c     If the volume fraction above efermi is needed, provide the negative of etetra_in and efermi instead.
c
c begin interface
      subroutine cuttetra(wtetra,vtetra,ntetra,etetra_in,efermi,tetra)
      use, intrinsic :: iso_fortran_env !inc
      implicit none
      real_dp,           intent(in)  :: efermi,etetra_in(4)
      real_dp,           intent(out) :: wtetra(4,3),vtetra(3)
      real_dp, optional, intent(out) :: tetra(4,4,3)
      integer,           intent(out) :: ntetra
c end interface
      real_dp                        :: etetra(4),f(4)
      integer                        :: pnt(4)

      call rorderp(pnt,etetra_in,4) !
      etetra = etetra_in(pnt)       ! etetra is size-ordered

      if(efermi>=etetra(4)) then
        ntetra      = 1
        wtetra(:,1) = 1
        vtetra(1)   = 1
        if(present(tetra)) tetra(pnt,:,1) = reshape ( [
     &    1d0,0d0,0d0,0d0, 0d0,1d0,0d0,0d0, 0d0,0d0,1d0,0d0, 0d0,0d0,0d0,1d0 ],[4,4] )
      else if(efermi>=etetra(3)) then
        ntetra      = 3
        f(1:3)      = (efermi-etetra(1:3)) / (etetra(4)-etetra(1:3))
        wtetra(pnt,1) = [ 1d0,    2-f(2), 1d0,    f(2)           ]
        wtetra(pnt,2) = [ 2-f(1), 1-f(2), 1-f(3), f(1)+f(2)+f(3) ]
        wtetra(pnt,3) = [ 1d0,    1-f(2), 2-f(3), f(2)+f(3)      ]
        vtetra(1)   = f(2)
        vtetra(2)   = f(1)*(1-f(2))*(1-f(3))
        vtetra(3)   = f(3)*(1-f(2))
        if(present(tetra)) tetra(pnt,:,:) = reshape ( (/
     &    1d0,0d0,0d0,0d0, 0d0,1d0,0d0,0d0,     0d0,0d0,1d0,0d0,     0d0,1-f(2),0d0,f(2),
     &    1d0,0d0,0d0,0d0, 1-f(1),0d0,0d0,f(1), 0d0,1-f(2),0d0,f(2), 0d0,0d0,1-f(3),f(3),
     &    1d0,0d0,0d0,0d0, 0d0,0d0,1d0,0d0,     0d0,1-f(2),0d0,f(2), 0d0,0d0,1-f(3),f(3) /),[4,4,3] )
      else if(efermi>=etetra(2)) then
        ntetra      = 3
        f(1:2)      = (efermi-etetra(1)) / (etetra(3:4)-etetra(1))
        f(3:4)      = (efermi-etetra(2)) / (etetra(3:4)-etetra(2))
        wtetra(pnt,1) = [ 3-f(1)-f(2), 1d0,         f(1),      f(2)      ]
        wtetra(pnt,2) = [ 2-f(1)-f(2), 2-f(4),      f(1),      f(2)+f(4) ]
        wtetra(pnt,3) = [ 1-f(1),      3-f(3)-f(4), f(1)+f(3), f(4)      ]
        vtetra(1)   = f(1)*f(2)
        vtetra(2)   = f(1)*(1-f(2))*f(4)
        vtetra(3)   = f(3)*(1-f(1))*f(4)
        if(present(tetra)) tetra(pnt,:,:) = reshape ( [
     &    1d0,0d0,0d0,0d0, 0d0,1d0,0d0,0d0,     1-f(1),0d0,f(1),0d0, 1-f(2),0d0,0d0,f(2),
     &    0d0,1d0,0d0,0d0, 1-f(1),0d0,f(1),0d0, 1-f(2),0d0,0d0,f(2), 0d0,1-f(4),0d0,f(4),
     &    0d0,1d0,0d0,0d0, 1-f(1),0d0,f(1),0d0, 0d0,1-f(4),0d0,f(4), 0d0,1-f(3),f(3),0d0 ],[4,4,3] )
      else if(efermi>=etetra(1)) then
        ntetra      = 1
        f(1:3)      = (efermi-etetra(1)) / (etetra(2:4)-etetra(1))
        wtetra(pnt,1) = [ 4-f(1)-f(2)-f(3), f(1), f(2), f(3) ]
        vtetra(1)   = f(1)*f(2)*f(3)
        if(present(tetra)) tetra(pnt,:,1) = reshape ( [
     &    1d0,0d0,0d0,0d0, 1-f(1),f(1),0d0,0d0, 1-f(2),0d0,f(2),0d0, 1-f(3),0d0,0d0,f(3) ],[4,4] )
      else
        ntetra = 0
      endif
      end

c     This version is a rewrite and returns a different wtetra:
c     wtetra(1:4)             : weights of the old tetrahedron corners: wtetra = sum(i,j) vtetra(j) * tetra(:,i,j)
c     vtetra(1:ntetra)        : volumes of the ntetra new tetrahedra
c     tetra(1:4,1:4,1:ntetra) : corner points of the new tetrahedra (optional)
c
c begin interface
      subroutine cuttetraN(wtetra,vtetra,ntetra,etetra_in,efermi,tetra)
      use, intrinsic :: iso_fortran_env !inc
      implicit none
      real_dp,           intent(in)  :: efermi,etetra_in(4)
      real_dp,           intent(out) :: wtetra(4),vtetra(3)
      real_dp, optional, intent(out) :: tetra(4,4,3)
      integer,           intent(out) :: ntetra
c end interface
      real_dp,           parameter   :: eta = 1d-14
      real_dp                        :: etetra(4),f(4)
      integer                        :: pnt(4)

      ! size-order etetra (could be a bit slower than rorderp, see cuttetra)
      pnt                            = 4
      pnt(minloc(etetra_in,1))       = 1
      pnt(minloc(etetra_in,1,pnt>1)) = 2
      pnt(minloc(etetra_in,1,pnt>2)) = 3
      pnt(pnt)                       = [ 1,2,3,4 ]
      etetra                         = etetra_in(pnt)

      ! define tetrahedra
      if(efermi>=etetra(4)-eta) then
        ntetra         = 1
        wtetra         = 1
        vtetra(1)      = 1
        if(present(tetra)) tetra(pnt,:,1) = reshape ( [
     &    1d0,0d0,0d0,0d0, 0d0,1d0,0d0,0d0, 0d0,0d0,1d0,0d0, 0d0,0d0,0d0,1d0 ],[4,4] )
      else if(efermi>=etetra(3)) then
        ntetra         = 3
        f(1:3)         = (efermi-etetra(1:3)) / (etetra(4)-etetra(1:3))
        vtetra(1)      = f(2)
        vtetra(2)      = f(1)*(1-f(2))*(1-f(3))
        vtetra(3)      = f(3)*(1-f(2))
        wtetra(pnt(1)) = vtetra(1)      + vtetra(2)*(2-f(1)) + vtetra(3)
        wtetra(pnt(2)) = vtetra(1)      + (vtetra(1)+vtetra(2)+vtetra(3))*(1-f(2))
        wtetra(pnt(3)) = vtetra(1)      + vtetra(2)*(1-f(3)) + vtetra(3)*(2-f(3))
        wtetra(pnt(4)) = vtetra(1)*f(2) + vtetra(2)*f(1)     + (vtetra(2)+vtetra(3))*(f(2)+f(3))
        if(present(tetra)) tetra(pnt,:,:) = reshape ( [
     &    1d0,0d0,0d0,0d0, 0d0,1d0,0d0,0d0,     0d0,0d0,1d0,0d0,     0d0,1-f(2),0d0,f(2),
     &    1d0,0d0,0d0,0d0, 1-f(1),0d0,0d0,f(1), 0d0,1-f(2),0d0,f(2), 0d0,0d0,1-f(3),f(3),
     &    1d0,0d0,0d0,0d0, 0d0,0d0,1d0,0d0,     0d0,1-f(2),0d0,f(2), 0d0,0d0,1-f(3),f(3) ],[4,4,3] )
      else if(efermi>etetra(2)) then
        ntetra         = 3
        f(1:2)         = (efermi-etetra(1)) / (etetra(3:4)-etetra(1))
        f(3:4)         = (efermi-etetra(2)) / (etetra(3:4)-etetra(2))
        vtetra(1)      = f(1)*f(2)
        vtetra(2)      = f(1)*(1-f(2))*f(4)
        vtetra(3)      = f(3)*(1-f(1))*f(4)
        wtetra(pnt(1)) = vtetra(1)*(3-f(1)-f(2)) + vtetra(2)*(2-f(1)-f(2)) + vtetra(3)*(1-f(1))
        wtetra(pnt(2)) = vtetra(1)               + vtetra(2)*(2-f(4))      + vtetra(3)*(3-f(3)-f(4))
        wtetra(pnt(3)) = (vtetra(1)+vtetra(2)+vtetra(3))*f(1)              + vtetra(3)*f(3)
        wtetra(pnt(4)) = (vtetra(1)+vtetra(2))*f(2) + (vtetra(2)+vtetra(3))*f(4)
        if(present(tetra)) tetra(pnt,:,:) = reshape ( [
     &    1d0,0d0,0d0,0d0, 0d0,1d0,0d0,0d0,     1-f(1),0d0,f(1),0d0, 1-f(2),0d0,0d0,f(2),
     &    0d0,1d0,0d0,0d0, 1-f(1),0d0,f(1),0d0, 1-f(2),0d0,0d0,f(2), 0d0,1-f(4),0d0,f(4),
     &    0d0,1d0,0d0,0d0, 1-f(1),0d0,f(1),0d0, 0d0,1-f(4),0d0,f(4), 0d0,1-f(3),f(3),0d0 ],[4,4,3] )
      else if(efermi>=etetra(1)+eta) then
        ntetra         = 1
        f(1:3)         = (efermi-etetra(1)) / (etetra(2:4)-etetra(1))
        vtetra(1)      = f(1)*f(2)*f(3)
        wtetra(pnt(1)) = vtetra(1) * (4-f(1)-f(2)-f(3))
        wtetra(pnt(2)) = vtetra(1) * f(1)
        wtetra(pnt(3)) = vtetra(1) * f(2)
        wtetra(pnt(4)) = vtetra(1) * f(3)
        if(present(tetra)) tetra(pnt,:,1) = reshape ( [
     &    1d0,0d0,0d0,0d0, 1-f(1),f(1),0d0,0d0, 1-f(2),0d0,f(2),0d0, 1-f(3),0d0,0d0,f(3) ],[4,4] )
      else
        ntetra = 0
        wtetra = 0
      endif
      end

c -----------------

c     Initializes the tetrahedron integration on a surface given by
c
c      s'           s
c     e  (q+k)  -  e (q)  =  h0     for occupied (n) and unoccupied states (n'),
c      n'           n
c
c     i.e. calculates (for any f) the weights wintgr3 in the sum
c
c                                         occ unocc                  2
c      SUM  wintgr3(n',n,k) * f(n',n,k) = SUM  SUM  INT   f(n',n,k) d k     (spin indices omitted) .
c     k,n,n'                               n    n'     S
c
      subroutine tetrahedron3_init(wintgr3,h0,ikpt,s1,s2,bandi,bandf)

      use global
      use, intrinsic :: iso_fortran_env
      implicit none
      integer, intent(in)    :: s1,s2,bandi,bandf,ikpt
      real_dp, intent(in)    :: h0
      real_dp, intent(out)   :: wintgr3(bandi:bandf,bando,nkpt)
      integer                :: tetra(4,6)
      integer                :: n,itetra,k1,k2,k3,n1,n2,i,itria,itria1
      integer                :: pnt(4),kk
      integer                :: ntria,ntria1,ntria2,bandi1,bandf1,nn
      integer                :: kcorn(8),kpt2(nkpt)
      real_dp                :: ene1(nkpt),ene2(nkpt)
      real_dp                :: ecube1(8),ecube2(8),etetra1(4),etetra2(4),f(4),hcube(8),htetra(4)
      real_dp                :: atria(2),etria1(3),etria2(3),etria2_0(3),rarr(3)
      real_dp                :: wtria1(3,2),atria1(2),wtria2(3,2),atria2(2),tria(4,4),tria1(3,3,2)
      real_dp                :: rdum
      logical                :: ldum1,ldum2,ldum3,ldum4
      integer                :: kptsum

      if(nkpt==1) Error('Tetrahedron integration not supported for single-k-point meshes.')

      wintgr3 = 0 ; if(h0<=0) return
      bandi1  = bandf
      bandf1  = bandi

      call get_tetra(tetra,n,1)

      do i = 1,nkpt
        kpt2(i) = kptsum(i,ikpt)
      enddo

      do n1 = 1,bando
        ene1 = ene(n1,:,s1)
        do n2 = bandi,bandf
c          if(n1/=n2) cycle
          ene2 = ene(n2,kpt2,s2)
          if(minval(ene2-ene1)>=h0) cycle
          if(maxval(ene2-ene1)<=h0) cycle
          do k3 = 0,nkpt3(3)-1     !
            do k2 = 0,nkpt3(2)-1   ! Loop over cubes in BZ
              do k1 = 0,nkpt3(1)-1 !

                kcorn(1) = pkpt(k1  ,k2  ,k3  ) !
                kcorn(2) = pkpt(k1+1,k2  ,k3  ) !
                kcorn(3) = pkpt(k1  ,k2+1,k3  ) !
                kcorn(4) = pkpt(k1+1,k2+1,k3  ) ! Definition of the cube corners
                kcorn(5) = pkpt(k1  ,k2  ,k3+1) !
                kcorn(6) = pkpt(k1+1,k2  ,k3+1) !
                kcorn(7) = pkpt(k1  ,k2+1,k3+1) !
                kcorn(8) = pkpt(k1+1,k2+1,k3+1) !_____

                ldum1 = .true.
                ldum2 = .true.
                ldum3 = .true.
                ldum4 = .true.
                do i = 1,8
                  n         = kcorn(i)              !
                  ecube1(i) = ene1(n)               ! Get energies at cube corners
                  ecube2(i) = ene2(n)               !     for bands n1 and n2
                  hcube (i) = ecube2(i) - ecube1(i) ! Get surface function at cube corners
                  if(ldum1) then ; if(ecube1(i)<efermi) ldum1 = .false. ; endif ! Nothing to be done if cube-corner energies of band n1 larger than efermi
                  if(ldum2) then ; if(ecube2(i)>efermi) ldum2 = .false. ; endif ! Nothing to be done if cube-corner energies of band n2 smaller than efermi
                  if(ldum3) then ; if(hcube (i)<h0)     ldum3 = .false. ; endif ! Nothing to be done if h0 does not fall inside cube
                  if(ldum4) then ; if(hcube (i)>h0)     ldum4 = .false. ; endif !
                enddo
                if(ldum1.or.ldum2.or.ldum3.or.ldum4) cycle

                do itetra = 1,6

                  ldum1 = .true.
                  ldum2 = .true.
                  ldum3 = .true.
                  ldum4 = .true.
                  do i = 1,4
                    n          = tetra(i,itetra) !
                    etetra1(i) = ecube1(n)       ! Get energies at tetrahedron corners for band n1
                    etetra2(i) = ecube2(n)       !                                     for band n2
                    htetra (i) = hcube (n)       ! Get surface function at tetrahedron corners
                    if(ldum1) then ; if(etetra1(i)<efermi) ldum1 = .false. ; endif ! Nothing to be done if tetrahedron-corner energies of band n1 larger than efermi
                    if(ldum2) then ; if(etetra2(i)>efermi) ldum2 = .false. ; endif ! Nothing to be done if tetrahedron-corner energies of band n2 smaller than efermi
                    if(ldum3) then ; if(htetra (i)<h0)     ldum3 = .false. ; endif ! Nothing to be done if h0 does not fall inside tetrahedron
                    if(ldum4) then ; if(htetra (i)>h0)     ldum4 = .false. ; endif !
                  enddo
                  if(ldum1.or.ldum2.or.ldum3.or.ldum4) cycle

                  call rorderp(pnt,htetra,4) !
                  htetra = htetra(pnt)       ! htetra is size-ordered

                  ! Intersecting plane is a triangle (ntria=1) or a quadrangle (ntria=2).
                  ! tria contains the triangle points and atria the triangle area(s).
                  tria = 0
                  if(h0>htetra(3)) then
                    ntria     = 1
                    f(1:3)    = (h0-htetra(1:3))/(htetra(4)-htetra(1:3))
                    atria(1)  = 3 * (h0-htetra(4))**2 /
     &                              ( (htetra(4)-htetra(1)) * (htetra(4)-htetra(2)) * (htetra(4)-htetra(3)) )
                    tria(1,1) = 1-f(1) ; tria(4,1) = f(1)
                    tria(2,2) = 1-f(2) ; tria(4,2) = f(2)
                    tria(3,3) = 1-f(3) ; tria(4,3) = f(3)
                  else if(h0>htetra(2)) then
                    ntria     = 2
                    f(1:2)    = (h0-htetra(1))/(htetra(3:4)-htetra(1))
                    f(3:4)    = (h0-htetra(2))/(htetra(3:4)-htetra(2))
                    rdum      = (htetra(4)-htetra(2)) * (htetra(3)-htetra(1))
                    atria(1)  = 3*f(3) * (htetra(3)-h0) / rdum
                    atria(2)  = 3*f(2) * (htetra(4)-h0) / rdum
                    tria(2,1) = 1-f(3) ; tria(3,1) = f(3) !
                    tria(1,2) = 1-f(1) ; tria(3,2) = f(1) ! The triangles are
                    tria(2,3) = 1-f(4) ; tria(4,3) = f(4) ! 1,2,3 and 2,3,4
                    tria(1,4) = 1-f(2) ; tria(4,4) = f(2) !
                  else ! if(h0>htetra(1)) then
                    ntria     = 1
                    f(1:3)    = (h0-htetra(1))/(htetra(2:4)-htetra(1))
                    atria(1)  = 3*f(1)*f(2) / (htetra(4)-htetra(1))
                    tria(1,1) = 1-f(1) ; tria(2,1) = f(1)
                    tria(1,2) = 1-f(2) ; tria(3,2) = f(2)
                    tria(1,3) = 1-f(3) ; tria(4,3) = f(3)
                  endif

                  tria(pnt,:) = tria

                  do itria = 1,ntria
                    etria1 = matmul(etetra1,tria(:,itria:itria+2))
                    call cuttria(wtria1,atria1,ntria1,etria1,efermi,tria1)
                    if(ntria1>0) then
                      etria2_0 = matmul(etetra2,tria(:,itria:itria+2))
                      if(n2<bandi1) bandi1 = n2 ! lower and upper index of nonzero weights
                      if(n2>bandf1) bandf1 = n2 !
                    endif
                    do itria1 = 1,ntria1
                      etria2 = matmul( etria2_0, tria1(:,:,itria1) )
                      call cuttria(wtria2,atria2,ntria2,-etria2,-efermi)
                      if     (ntria2==0) then ; cycle
                      else if(ntria2==1) then ; rarr =   wtria2(:,1) * atria2(1)
                      else                    ; rarr = ( wtria2(:,1) * atria2(1) + wtria2(:,2) * atria2(2) )
                      endif
                      rdum = atria(itria) * atria1(itria1)
                      do i = 1,4
                        kk                = kcorn(tetra(i,itetra))
                        wintgr3(n2,n1,kk) = wintgr3(n2,n1,kk) + rdum * dot_product ( tria(i,itria:itria+2), rarr )
                      enddo
                    enddo
                  enddo

                enddo
              enddo
            enddo
          enddo
        enddo
      enddo

      rdum = 3d0*6*nkpt                                              !
      wintgr3(bandi1:bandf1,:,:) = wintgr3(bandi1:bandf1,:,:) / rdum ! take average over three triangle corners,
                                                                     ! six tetrahedra per cube,
                                                                     ! nkpt cubes

      if(alignbd) return ! We have to return here in order to make the residual wings small for metals. Calculation with
                         ! a small kptadd is not improved, though. It seems that disentanglement of states does not work
                         ! well here.

c      return

      ! Average over degenerate eigenstates
      ! - occupied states
      do k1 = 1,nkpt
        n = 1
        do while(n<=bando)
          n1 = deg(n,k1,s1) ; if(n1>bando) exit
          do n2 = bandi1,bandf1
            rdum = sum(wintgr3(n2,n:n1,k1)) / (n1-n+1) ; do i = n,n1 ; wintgr3(n2,i,k1) = rdum ; enddo
          enddo
          n  = n1 + 1
        enddo
      enddo
      ! - unoccupied states
      do k1 = 1,nkpt
        k2 = kpt2(k1)
        nn = min(nband(k2,s2),bandf)
        n  = bandi
        do while(deg(n,k2,s2)<n.and.n<nn)
          n = n + 1
        enddo
        if(deg(n,k2,s2)>=n) then
          do while(n<=nn)
            n1 = min(deg(n,k2,s2),nn) !; if(n>bandf1.or.n1<bandi1) then ; n = n1 + 1 ; cycle ; endif
            do n2 = 1,bando
              rdum = sum(wintgr3(n:n1,n2,k1)) / (n1-n+1) ; do i = n,n1 ; wintgr3(i,n2,k1) = rdum ; enddo
            enddo
            n  = n1 + 1
          enddo
        endif
      enddo

      contains

c -----------------

c     Cuts the triangle at the Fermi energy efermi.
c     The area fraction below efermi is given by ntria (1 or 2) triangles defined by
c     wtria(1:4,1:ntria)    : weights of the three triangle corners
c     atria(1:ntria)        : areas of the ntria new triangles
c     tria(1:4,1:4,1:ntria) : corner points of the new triangles (optional)
c
c     If the area fraction above efermi is needed, provide the negative of etria and efermi instead.
c
      subroutine cuttria(wtria,atria,ntria,etria_in,efermi,tria)

      use, intrinsic :: iso_fortran_env
      implicit none
      real_dp,           intent(in)  :: efermi,etria_in(3)
      real_dp,           intent(out) :: wtria(3,2),atria(2)
      real_dp, optional, intent(out) :: tria(3,3,2)
      integer,           intent(out) :: ntria
      real_dp                        :: etria(3),f(2)
      integer                        :: pnt(3)

      if(etria_in(1)>=efermi) then
        if(etria_in(2)>=efermi) then
          if(etria_in(3)>=efermi) then ; ntria = 0 ; return
          endif
        endif
      endif

      etria = etria_in

      call rorderp(pnt,etria,3) !
      etria = etria(pnt)        ! etria is size-ordered

      if(efermi>=etria(3)) then
        ntria        = 1
        wtria(:,1)   = 1
        atria(1)     = 1
        if(present(tria)) then
          tria(pnt,:,1) = reshape ( (/
     &      1d0,0d0,0d0, 0d0,1d0,0d0, 0d0,0d0,1d0 /),[3,3] )
        endif
      else if(efermi>=etria(2)) then
        ntria        = 2
        f            = (efermi-etria(1:2))/(etria(3)-etria(1:2))
        wtria(pnt,1) = [ 1d0,    2-f(2), f(2)      ]
        wtria(pnt,2) = [ 2-f(1), 1-f(2), f(1)+f(2) ]
        atria(1)     = f(2)
        atria(2)     = f(1)*(1-f(2))
        if(present(tria)) then
          tria(pnt,:,:) = reshape ( (/
     &      1d0,0d0,0d0, 0d0,1d0,0d0,     0d0,1-f(2),f(2),
     &      1d0,0d0,0d0, 0d0,1-f(2),f(2), 1-f(1),0d0,f(1) /),[3,3,2] )
        endif
      else ! if(efermi>etria(1)) then
        ntria        = 1
        f            = (efermi-etria(1))/(etria(2:3)-etria(1))
        wtria(pnt,1) = [ 3-f(1)-f(2), f(1), f(2) ]
        atria(1)     = f(1)*f(2)
        if(present(tria)) then
          tria(pnt,:,1) = reshape ( (/
     &      1d0,0d0,0d0, 1-f(1),f(1),0d0, 1-f(2),0d0,f(2) /),[3,3] )
        endif
      endif

      end subroutine cuttria

      end

c -----------------

c     Initializes the tetrahedron integration on the Fermi surface,
c     i.e. calculates (for any f) the weights wintgr4 in the sum
c
c                                   occ                     2
c      SUM  wintgr4(n,k) * f(n,k) = SUM     INT     f(n,k) d k     (spin index omitted) .
c      k,n                           n   e (k) = E
c                                         n       F
c
      subroutine tetrahedron4_init(wintgr4,s)

      use global
      use, intrinsic :: iso_fortran_env
      implicit none
      integer, intent(in)    :: s
      real_dp, intent(out)   :: wintgr4(bando,nkpt)
      integer                :: tetra(4,6)
      integer                :: n,itetra,k1,k2,k3,i,itria,n1
      integer                :: pnt(4),kk
      integer                :: ntria
      integer, allocatable   :: kcorn(:)
      real_dp                :: ene1(nkpt),ecube(8),etetra(4),f(4),atria(2),wtria(4,2),rdum
      real_dp                :: whlp(bando)

      wintgr4 = 0

      call get_tetra(tetra,n,1)

      allocate ( kcorn(8) )

      do n = 1,bando
        ene1 = ene(n,:nkpt,s)
        if(minval(ene1)>=efermi) cycle
        if(maxval(ene1)<=efermi) cycle
        do k3 = 0,nkpt3(3)-1     !
          do k2 = 0,nkpt3(2)-1   ! Loop over cubes in BZ
            do k1 = 0,nkpt3(1)-1 !

              kcorn(1) = pkpt(k1  ,k2  ,k3  ) !
              kcorn(2) = pkpt(k1+1,k2  ,k3  ) !
              kcorn(3) = pkpt(k1  ,k2+1,k3  ) !
              kcorn(4) = pkpt(k1+1,k2+1,k3  ) ! Definition of the cube corners
              kcorn(5) = pkpt(k1  ,k2  ,k3+1) !
              kcorn(6) = pkpt(k1+1,k2  ,k3+1) !
              kcorn(7) = pkpt(k1  ,k2+1,k3+1) !
              kcorn(8) = pkpt(k1+1,k2+1,k3+1) !_____
              ecube    = ene1(kcorn) ! Get energies at cube corners

              do itetra = 1,6

                etetra = ecube(tetra(:,itetra)) ! Get energies at tetrahedron corners
                if(all(etetra>=efermi)) cycle ! Nothing to be done if tetrahedron-corner energies larger than efermi
                if(all(etetra<=efermi)) cycle ! Nothing to be done if tetrahedron-corner energies smaller than efermi

                call rorderp(pnt,etetra,4) !
                etetra = etetra(pnt)       ! etetra is size-ordered

                ! Intersecting plane is a triangle (ntria=1) or a quadrangle (ntria=2).
                ! tria contains the triangle points and atria the triangle area(s).
                if(efermi>etetra(3)) then
                  ntria       = 1
                  f(1:3)      = (efermi-etetra(1:3))/(etetra(4)-etetra(1:3))
                  atria(1)    = 3 * (efermi-etetra(4))**2 /
     &                              ( (etetra(4)-etetra(1)) * (etetra(4)-etetra(2)) * (etetra(4)-etetra(3)) )
                  wtria(:3,1) = 1-f(:3)
                  wtria(4,1)  = f(1)+f(2)+f(3)
                else if(efermi>etetra(2)) then
                  ntria      = 2
                  f(1:2)     = (efermi-etetra(1))/(etetra(3:4)-etetra(1))
                  f(3:4)     = (efermi-etetra(2))/(etetra(3:4)-etetra(2))
                  rdum       = (etetra(4)-etetra(2)) * (etetra(3)-etetra(1))
                  atria(1)   = 3*f(3) * (etetra(3)-efermi) / rdum
                  atria(2)   = 3*f(2) * (etetra(4)-efermi) / rdum
                  wtria(1,1) = 1-f(1)      ; wtria(1,2) = 2-f(1)-f(2)
                  wtria(2,1) = 2-f(3)-f(4) ; wtria(2,2) = 1-f(4)
                  wtria(3,1) = f(1)+f(3)   ; wtria(3,2) = f(1)
                  wtria(4,1) = f(4)        ; wtria(4,2) = f(2)+f(4)
                else ! if(h0>htetra(1)) then
                  ntria       = 1
                  f(1:3)      = (efermi-etetra(1))/(etetra(2:4)-etetra(1))
                  atria(1)    = 3*f(1)*f(2) / (etetra(3)-etetra(1))
                  wtria(1,1)  = 3-f(1)-f(2)-f(3)
                  wtria(2:,1) = f(:3)
                endif

                wtria(pnt,:) = wtria

                do itria = 1,ntria
                  do i = 1,4
                    kk            = kcorn(tetra(i,itetra))
                    wintgr4(n,kk) = wintgr4(n,kk) + atria(itria) * wtria(i,itria)
                  enddo
                enddo

              enddo
            enddo
          enddo
        enddo
      enddo

      wintgr4 = wintgr4 / 3 / 6 / nkpt ! take average over three triangle corners,
     &                                 ! six tetrahedra per cube,
     &                                 ! nkpt cubes
      deallocate(kcorn)

      ! Average over symmetry-equivalent k-points
      allocate( kcorn(nsym) )
      do k1 = 1,nkpti
        n     = 0
        kcorn = 0
        do k2 = 1,nkpt
          if(kptp(k2)==k1) then
            n        = n + 1
            kcorn(n) = k2
          endif
        enddo
        whlp = 0
        do n1 = 1,n
          whlp = whlp + wintgr4(:,kcorn(n1))
        enddo
        whlp = whlp / n
        do n1 = 1,n
          wintgr4(:,kcorn(n1)) = whlp
        enddo
      enddo

      ! Average over degenerate eigenstates (should be rewritten using array deg)
      do k1 = 1,nkpt
        n  = 1
        do while(n<bando)
          n1 = n
          do while(abs(ene(n1+1,k1,s)-ene(n,k1,s))<1d-8)
            n1 = n1 + 1
            if(n1==bando) exit
          enddo
          rdum             = sum(wintgr4(n:n1,k1))
          wintgr4(n:n1,k1) = rdum / (n1-n+1)
          n                = n1 + 1
        enddo
      enddo

      deallocate ( kcorn )

      end

c -----------------

c     Initializes the tetrahedron integration on a surface given by
c
c      s'           s
c     e  (q+k)  -  e (q)  =  w(i)     for occupied (n) and unoccupied states (n') as well as w(i) = frq(i),
c      n'           n
c
c     i.e., calculates (for any f) the weights wintgr5 in the sum
c
c                                         occ unocc                  2
c      SUM  wintgr5(n',n,k) * f(n',n,k) = SUM  SUM  INT   f(n',n,k) d k     (spin indices omitted) .
c     k,n,n'                               n    n'     S(w)
c
c     To be more precise, in order to catch resonances in the vicinity of w(i), the corresponding weight is integrated
c     over the symmetric interval [ w(i) - dw(i) , w(i) + dw(i) ] and then divided by 2*dw(i), where
c
c     dw(i) = a * (w(i+1)-w(i)) = (1-a) * (w(i)-w(i-1)) .
c
c     The last equality ensures that the integration regions touch. Consequently, one has
c
c     a = (w(i)-w(i-1)) / (w(i+1)-w(i-1))
c
c     which is independent of the index i if the set of frequencies has the form w(i) = b*c**i + d.
c     This is tested at the beginning. Note that for the first frequency (i=1), it is assumed that
c     there is no contribution from frequencies below the lower bound w(1)-dw(1). (So, the lower bound
c     of the integration is effectively set at -infinity.)
c
c     If nfrq=0, minimum and maximum frequencies are returned in the first two values of wintgr5. (Can be called by only one rank.)
c
c     The weights of degenerate states are averaged. If b1low..b1up (or b2low..b2up) cuts degenerate subspaces, the range
c     is extended to include the whole subspace(s) (-> d1low,d1up,d2low,d2up). Can be disabled:
c # define disable_cut_average
c
c     frq(:nfrq)   - frequencies w(i)
c     ikpt         - k index
c     kpt1(:nkpt1) - q indices
c     kpt1p(:nkpt) - parents of kpoints in kpt1 (can be 0 if parent is not in kpt1), see below for MPI!
c     kflag        - 1 (normal), 2 (q and q+k switched)
c     s1/2         - s/s'
c     b1/2|low/up  - n=(b1low..b1up), n'=(b2low..b2up)
c
c     MPI: No parallelization done in this routine. So, each process runs independently and
c          only on its share of b1low..b1up, b2low..b2up, kpt1.
c          The array kpt1p of parent k points must be defined in terms of the "local" kpt1.
c          For k points that do not have a parent in kpt1, kpt1p(k) must be zero.
c          (This is identical to the related routine in gauss.f.)
c
c     If b1low=b1up=0, a dispersion-less occupied state at efermi is assumed (for the treatment of core states).
c
c     Array wintgr5 declared as intent(inout) because only part of the array might be defined. The rest should remain
c     unchanged, which is not guaranteed by the Fortran standard if intent(out) is used.
c
# ifdef old_version
#   define New(arg)
# else
#   define New(arg) arg
# endif

      subroutine tetrahedron5_init(wintgr5,frq,nfrq,ikpt,kpt1,kpt1p,nkpt1,kflag,s1,s2,b1low,b1up,b2low,b2up,
     &                             dim1low,dim1up,dim2low,dim2up)

      use global
      use, intrinsic :: iso_fortran_env
      implicit none
# include "interface/cuttetraN.inc"
      integer, intent(in)    :: s1,s2,b1low,b1up,b2low,b2up,dim1low,dim1up,dim2low,dim2up
      integer, intent(in)    :: ikpt,nfrq,nkpt1,kflag
      integer, intent(in)    :: kpt1(nkpt1),kpt1p(nkpt)
      real_dp, intent(in)    :: frq(nfrq)
      real_dp, intent(inout) :: wintgr5(dim2low:dim2up,dim1low:dim1up,nkpt1,max(1,nfrq))
      integer, allocatable   :: ntetra1(:,:)
      integer                :: ntetra2(6),ntetra,ntetra_(6)
      integer                :: tetra(4,6)
      real_dp, allocatable   :: vtetra1(:,:,:),wtetra1(:,:,:),tetra1(:,:,:,:,:)
      real_dp                :: vtetra2(3,6),wtetra2(4,6),wtetra_(4,6),tetra2(4,4,3,6),tetra1_(4,4,3)
      real_dp                :: vtetra(3),wtetra(4),tetra_(4,4,9,6),vtetra_(9,6)
      real_dp                :: ecube1(8),ecube2(8),dcube(8),wadd(8),wadd0(8)
      real_dp                :: etetra1(4),etetra1_(4),etetra2(4),dtetra(4),dtetra_(4)
      real_dp                :: frq1(nfrq)
      real_dp                :: dmin,dmax,emax1,emin2
      real_dp                :: rdum,a,f
      logical                :: kskip(8),kdone(nkpt1),ldum
      integer                :: kcube(8),kcube1(8),kcube2(8),kcubep(8)
      integer                :: kpt2(nkpt),isym
      integer                :: i,ifrq,k,k1,k2,k3,kk,n,n1,n2,nn
      integer                :: itetra,itetra_,itetra1,itetra2
      integer                :: d1low(nkpt),d1up(nkpt),d2low(nkpt),d2up(nkpt)
      integer                :: mind1,maxd1,mind2,maxd2
      integer                :: kptsum
      real                   :: cputime

c      if(option==2) then
c        if(b1low==0) Error('Option=2 disabled for b1low=0')
c        if(any([dim1low,dim1up,dim2low,dim2up]/=[b1low,b1up,b2low,b2up]))
c     &    Error('Option=2 disabled for dim1low<>b1low etc.')
c        call tetrahedron5_init_old(wintgr5,frq,nfrq,ikpt,kpt1,kpt1p,nkpt1,kflag,s1,s2,b1low,b1up,b2low,b2up)
c        return
c      endif
      if(b1low>b1up.or.b2low>b2up.or.nkpt1==0) return
      if((b1low<dim1low.or.b1up>dim1up.or.b2low<dim2low.or.b2up>dim2up).and.nfrq>0)
     &                                    Bug('Band indices out of bounds.')
      if(b1low*b1up==0.and.b1low+b1up/=0) Bug('One of b1low/b1up is zero.')
      if(max(b1up,b2up)>maxeband)         Bug('Upper band index larger than maxeband.')

      call cpu_time(cputime)
      call get_tetra(tetra,n,1)

      ! Check whether frq(:) is in ascending order and calculate factor a (see above).
      if(nfrq==1) Error('Number of Hilbert frequencies must be larger than one.')
      if(nfrq==2) Warn('Only two Hilbert frequencies. This has never been tested.')
      if(nfrq>=2) then
        rdum = 0
        do ifrq = 2,nfrq
          if(frq(ifrq)<=frq(ifrq-1)) Error('input frequency mesh not properly ordered.')
          if(ifrq/=nfrq) then
            a    = (frq(ifrq)-frq(ifrq-1)) / (frq(ifrq+1)-frq(ifrq-1))
            rdum = rdum  + a
            if(abs(rdum/(ifrq-1)-a)>1d-10) Bug('Hilbert mesh not exponential.')
          endif
        enddo
        a = 0.5d0
        if(nfrq>2) a = rdum / (nfrq-2)
      endif

      ! Calculate frequency range midpoints w(i) -> frq1
      if(nfrq>0) then
        do ifrq = 1,nfrq-1
          frq1(ifrq) = frq(ifrq) + a * ( frq(ifrq+1)-frq(ifrq) )
        enddo
        frq1(nfrq) = frq(nfrq) New( + (1-a) * ( frq(nfrq)-frq(nfrq-1) ) )
      endif

      ! Define kpt2 : kpt(:,ikpt) + kpt(:,i)
      do i = 1,nkpt
        kpt2(i) = kptsum(i,ikpt)
      enddo

      ! Add bands below and above to average over degenerate states
      d1low = b1low ; d1up = b1up
      d2low = b2low ; d2up = b2up
# ifndef disable_cut_average
      if(nfrq/=0) then
        do k = 1,nkpt
          if(kflag==1) then ; k1 = k ; k2 = kpt2(k)
          else              ; k2 = k ; k1 = kpt2(k)
          endif
          if(b1low/=0) then
            d1low(k) = min(b1low,deg(b1low,k1,s1))
            d1up(k)  = deg(b1up,k1,s1)
          endif
          d2low(k) = min(b2low,deg(b2low,k2,s2))
          d2up(k)  = deg(b2up,k2,s2)
          if(d1up(k)<b1up) then
            if(d1up(k)==0) then ; d1up(k) = b2up ! deg can be zero if b1up > nband(k1,s1)
            else                ; d1up(k) = deg(d1up(k),k1,s1)
            endif
          endif
          if(d2up(k)<b2up) then
            if(d2up(k)==0) then ; d2up(k) = b2up ! deg can be zero if b2up > nband(k2,s2)
            else                ; d2up(k) = deg(d2up(k),k2,s2)
            endif
          endif
          d1up(k) = min(d1up(k),maxeband)
          d2up(k) = min(d2up(k),maxeband)
        enddo
      endif
# endif

      ! Determine weights with tetrahedron method.
      if(nfrq<=0) then ! return maximum and minimum frequency in wintgr5
        if(size(wintgr5)<2) Error('need at least size(wintgr5)=2.')
        call twoval(wintgr5,huge(0d0),-huge(0d0),0)
      else
        wintgr5(b2low:b2up,b1low:b1up,:,:) = 0
      endif

c      call cpu_done(cputime)

      do k = 1,nkpt

        k1       = nint( kpt(1,k)*nkpt3(1) )
        k2       = nint( kpt(2,k)*nkpt3(2) )
        k3       = nint( kpt(3,k)*nkpt3(3) )
        kcube(1) = pkpt(k1  ,k2  ,k3  ) !
        kcube(2) = pkpt(k1+1,k2  ,k3  ) !
        kcube(3) = pkpt(k1  ,k2+1,k3  ) !
        kcube(4) = pkpt(k1+1,k2+1,k3  ) ! Definition of the cube corners
        kcube(5) = pkpt(k1  ,k2  ,k3+1) !
        kcube(6) = pkpt(k1+1,k2  ,k3+1) !
        kcube(7) = pkpt(k1  ,k2+1,k3+1) !
        kcube(8) = pkpt(k1+1,k2+1,k3+1) !_____
        kcubep   = kpt1p(kcube)
# ifdef MPI
        if( all(kcubep==0) ) cycle
# endif

        if(kflag==1) then ; kcube1 = kcube ; kcube2 = kpt2(kcube) ! kcube1/2 relevant for ene
        else              ; kcube2 = kcube ; kcube1 = kpt2(kcube) !
        endif

        mind1 = minval(d1low(kcube))
        maxd1 = maxval(d1up (kcube))
        mind2 = minval(d2low(kcube))
        maxd2 = maxval(d2up (kcube))

        ! cut occupied bands
        allocate(ntetra1(6,mind1:maxd1),vtetra1(3,6,mind1:maxd1),wtetra1(4,6,mind1:maxd1),tetra1(4,4,3,6,mind1:maxd1))
        ntetra1 = -1
        if(b1low/=0) then
          do n1 = mind1,maxd1
            ecube1 = ene(n1,kcube1,s1)
            emax1  = maxval(ecube1)
            if(emax1>efermi) then
              do itetra = 1,6
                etetra1 = ecube1(tetra(:,itetra))
                call cuttetraN(wtetra1(:,itetra,n1),vtetra1(:,itetra,n1),ntetra1(itetra,n1),etetra1,efermi,tetra1(:,:,:,itetra,n1))
              enddo
            endif
          enddo
        endif

c
c       Loop over virtual transitions n1->n2
        do n2 = mind2,maxd2                 ! n2 can extend beyond (b2low,b2up) (to avoid cutting of degenerate states)
          nn     = min(b2up,max(b2low,n2))  ! nn is restricted to (b2low,b2up)
          ecube2 = ene(n2,kcube2,s2) ; if(all(ecube2<efermi)) cycle
          emin2  = minval(ecube2)

          ! cut unoccupied bands
          ntetra2 = -1
          if(emin2<efermi) then
            do itetra = 1,6
              etetra2 = ecube2(tetra(:,itetra))
              call cuttetraN(wtetra2(:,itetra),vtetra2(:,itetra),ntetra2(itetra),-etetra2,-efermi,tetra2(:,:,:,itetra))
            enddo
          endif

          do n1 = mind1,maxd1                ! n1 can extend beyond (b1low,b1up) (to avoid cutting of degenerate states)
            if(n1/=0) then ; ecube1 = ene(n1,kcube1,s1) ; if(all(ecube1>efermi)) cycle ! normal case
            else           ; ecube1 = efermi-1d-8                                      ! special case dispersion-less core state
            endif
            n      = min(b1up,max(b1low,n1)) ! n is restricted to (b1low,b1up)
            dcube  = ecube2 - ecube1
            dmin   = minval(dcube)
            dmax   = maxval(dcube)
            kskip  = n1<d1low(kcube).or.n1>d1up(kcube).or.
     &               n2<d2low(kcube).or.n2>d2up(kcube).or.kcubep==0
            if(tetraf) then
              do i = 1,8
                kskip(i) = kskip(i).or.all(kcube(i)/=kpt1)
              enddo
            endif

            ! prepare "default" wadd0 (used if dcube is not cut), and cut ecube1 into occupied and ecube2 into unoccupied tetrahedra (->tetra_,ntetra_)
            if(ntetra1(1,n1)==-1.and.ntetra2(1)==-1) then ! ecube1 and ecube2 are not cut by efermi
              wadd0 = [ (count(tetra==i),i=1,8) ]
            else
              do itetra = 1,6
                ntetra_  (itetra) = 0
                wtetra_(:,itetra) = 0 ; if(all(kskip(tetra(:,itetra)))) cycle
                if(ntetra2(1)==-1) then          ! only etetra1 is cut by efermi (just copy tetra1)
                  ntetra_(itetra)       = ntetra1(itetra,n1)
                  wtetra_(:,itetra)     = wtetra1(:,itetra,n1)
                  vtetra_(:3,itetra)    = vtetra1(:,itetra,n1)
                  tetra_(:,:,:3,itetra) = tetra1(:,:,:,itetra,n1)
                else if(ntetra1(1,n1)==-1) then  ! only etetra2 is cut by efermi (just copy tetra2)
                  ntetra_(itetra)       = ntetra2(itetra)
                  wtetra_(:,itetra)     = wtetra2(:,itetra)
                  vtetra_(:3,itetra)    = vtetra2(:3,itetra)
                  tetra_(:,:,:3,itetra) = tetra2(:,:,:3,itetra)
                else                             ! etetra1 and etetra2 are cut by efermi
                  etetra1           = ecube1(tetra(:,itetra))
                  do itetra2 = 1,ntetra2(itetra)
                    etetra1_ = matmul ( etetra1 , tetra2(:,:,itetra2,itetra) )
                    call cuttetraN(wtetra,vtetra,ntetra,etetra1_,efermi,tetra1_ )
                    wtetra_(:,itetra) = wtetra_(:,itetra) + vtetra2(itetra2,itetra) * matmul( tetra2(:,:,itetra2,itetra) , wtetra )
                    do itetra1 = 1,ntetra
                      ntetra_(                   itetra) = ntetra_(itetra) + 1
                      vtetra_(   ntetra_(itetra),itetra) = vtetra2(itetra2,itetra) * vtetra(itetra1)
                      tetra_(:,:,ntetra_(itetra),itetra) = matmul( tetra2(:,:,itetra2,itetra) , tetra1_(:,:,itetra1) )
                    enddo
                  enddo
                endif
              enddo
              wadd0 = 0
              do itetra = 1,6
                wadd0(tetra(:,itetra)) = wadd0(tetra(:,itetra)) + wtetra_(:,itetra)
              enddo
            endif

            ! nfrq=0 : return minimum and maximum frequency
            if(nfrq==0) then
              if(ntetra1(1,n1)==-1.and.ntetra2(1)==-1) then
                call twoval(wintgr5,dmin,dmax,1)
              else
                do itetra = 1,6                  ; dtetra  = dcube(tetra(:,itetra))
                  do itetra_ = 1,ntetra_(itetra) ; dtetra_ = matmul(dtetra,tetra_(:,:,itetra_,itetra))
                    call twoval(wintgr5,minval(dtetra_),maxval(dtetra_),1)
                  enddo
                enddo
              endif
            endif

            ! cut dcube (band-energy difference)
            do ifrq = 1,nfrq
              f = frq1(ifrq)

              if(f<dmin) then      ! (a) f below dmin -> dcube does not contribute
                cycle
              else if(f>dmax) then ! (b) f above dmax -> whole dcube contributes
                wadd = wadd0
              else                 ! (c) f betwen dmin and dmax -> dcube is cut into tetrahedra by f
                if(ntetra1(1,n1)==-1.and.ntetra2(1)==-1) then ! (c0) ecube1 and ecube2 are not cut by efermi
                  wadd = 0
                  do itetra = 1,6
                    dtetra = dcube(tetra(:,itetra))
                    if     (all(dtetra>=f)) then ; cycle
                    else if(all(dtetra< f)) then ; wadd(tetra(:,itetra)) = wadd(tetra(:,itetra)) + 1
                    else
                      call cuttetraN(wtetra,vtetra,ntetra,dtetra,f)
                      wadd(tetra(:,itetra)) = wadd(tetra(:,itetra)) + wtetra
                    endif
                  enddo
                else                                          ! (c1) ecube1 or ecube2 are cut by efermi (use tetra_,ntetra_ defined above
                  wadd = 0
                  do itetra = 1,6 ; if(ntetra_(itetra)==0) cycle
                    dtetra = dcube(tetra(:,itetra))
                    if     (all(dtetra>=f)) then ; cycle
                    else if(all(dtetra< f)) then ; wadd(tetra(:,itetra)) = wadd(tetra(:,itetra)) + wtetra_(:,itetra)
                    else
                      do itetra_ = 1,ntetra_(itetra)
                        dtetra_ = matmul( dtetra , tetra_(:,:,itetra_,itetra) )
                        if     (all(dtetra_>=f)) then ; cycle
                        else if(all(dtetra_< f)) then ; wadd(tetra(:,itetra)) = wadd(tetra(:,itetra)) + vtetra_(itetra_,itetra) *
     &                                                                            sum(tetra_(:,:,itetra_,itetra),2)
                        else
                          call cuttetraN(wtetra,vtetra,ntetra,dtetra_,f)
                          wadd(tetra(:,itetra)) = wadd(tetra(:,itetra)) + vtetra_(itetra_,itetra) *
     &                                            matmul( tetra_(:,:,itetra_,itetra) , wtetra )
                        endif
                      enddo
                    endif
                  enddo
                endif
              endif

              do i = 1,8 ; if(kskip(i)) cycle
                kk                     = kcubep(i)
                wintgr5(nn,n,kk,ifrq)  = wintgr5(nn,n,kk,ifrq) + wadd(i)
              enddo

            enddo

          enddo
        enddo

        deallocate(ntetra1,vtetra1,wtetra1,tetra1)
      enddo

c      call cpu_done(cputime)

      if(nfrq==0) then
        if(wintgr5(dim2low,dim1low,1,1)<-1d-12) Error('Frequency minimum negative.')
        wintgr5(dim2low,dim1low,1,1) = max(0d0,wintgr5(dim2low,dim1low,1,1))
        return
      endif

      if(tetraf) then
        do k = 1,nkpt1
          wintgr5(b2low:b2up,b1low:b1up,k,:) = wintgr5(b2low:b2up,b1low:b1up,k,:) * count(kpt1p==k)
        enddo
      endif

# ifdef old_version
      frq1(nfrq) = frq1(nfrq) + (1-a) * ( frq(nfrq)-frq(nfrq-1) ) / 2
# endif
      rdum = 4d0*6*nkpt
      do ifrq = nfrq,2,-1
        f                                     = 2 * (frq1(ifrq)-frq(ifrq))
        wintgr5(b2low:b2up,b1low:b1up,:,ifrq) = ( wintgr5(b2low:b2up,b1low:b1up,:,ifrq) -
     &                                            wintgr5(b2low:b2up,b1low:b1up,:,ifrq-1) ) / (f*rdum)
      enddo
      wintgr5(b2low:b2up,b1low:b1up,:,1) = wintgr5(b2low:b2up,b1low:b1up,:,1) / (2*(frq1(1)-frq(1))*rdum)

      if(frq(1)==0) wintgr5(b2low:b2up,b1low:b1up,:,1) = 0 ! for FSPEC

      if(alignbd) return ! We have to return here in order to make the residual wings small for metals. Calculation with
                         ! a small kptadd is improved, too, provided that the bands were aligned before with a very very
                         ! small kptadd.

c      call cpu_done(cputime)

      ! Average over degenerate eigenstates
      ! - occupied states
      if(b1low/=0) then
        do k = 1,nkpt1
          k1 = kpt1(k) ; if(kflag==2) k1 = kpt2(k1)
          n  = min(b1low,deg(b1low,k1,s1))
          do while(n<=min(b1up,nband(k1,s1)))
            nn    = min(deg(n,k1,s1),maxeband)
            maxd1 = max(b1low,n)
            mind1 = min(b1up,nn)
            if(maxd1<=mind1) then
              do ifrq = 1,nfrq
                do n2 = b2low,b2up
                  rdum                           = sum(wintgr5(n2,maxd1:mind1,k,ifrq)) / (nn-n+1)
                  wintgr5(n2,maxd1:mind1,k,ifrq) = rdum
                enddo
              enddo
            endif
            n = nn + 1
          enddo
        enddo
      endif
c      call cpu_done(cputime)
      ! - unoccupied states
      do k = 1,nkpt1
        k2 = kpt1(k) ; if(kflag==1) k2 = kpt2(k2)
        n  = min(b2low,deg(b2low,k2,s2))
        do while(n<=min(b2up,nband(k2,s2)))
          nn    = min(deg(n,k2,s2),maxeband)
          maxd2 = max(b2low,n)
          mind2 = min(b2up,nn)
          if(maxd2<=mind2) then
            do ifrq = 1,nfrq
              do n1 = b1low,b1up
                rdum                           = sum(wintgr5(maxd2:mind2,n1,k,ifrq)) / (nn-n+1)
                wintgr5(maxd2:mind2,n1,k,ifrq) = rdum
              enddo
            enddo
          endif
          n = nn + 1
        enddo
      enddo

c      call cpu_done(cputime)

c     Average over the k points explicitly if there are equivalents in kpt1 (eg., if "switch_off_symmetry" is defined in susceptibility)
      kdone = .false.
      do k = 1,nkpt1
        if(kdone(k)) cycle
        ldum = .false.
        kpt2 = 0
        do isym = 1,nsym
          if(kptsym(ikpt,isym)==ikpt) then
            k1 = kpt1p( kptsym(kpt1(k),isym) )
            if(k1/=0) then
              kpt2(k1) = kpt2(k1) + 1
              ldum     = ldum.or.k1/=k
            endif
          endif
        enddo
        if(ldum) then
          n = 0
          do k1 = 1,nkpt1
            if(kpt2(k1)>0) then
              n         = n + 1
              d1up(n)   = kpt2(k1)
              kpt2(n)   = k1
              kdone(k1) = .true.
            endif
          enddo
          nn = sum(d1up(:n))
          do ifrq = 1,nfrq
            do n1 = b1low,b1up
              do n2 = b2low,b2up
                rdum = 0         ; do k1 = 1,n ; rdum = rdum + d1up(k1) * wintgr5(n2,n1,kpt2(k1),ifrq) ; enddo
                rdum = rdum / nn ; do k1 = 1,n ; wintgr5(n2,n1,kpt2(k1),ifrq) = rdum                   ; enddo
              enddo
            enddo
          enddo
        endif
      enddo
c      call cpu_done(cputime)

      contains
      subroutine twoval(arr,r1,r2,mode)
      use, intrinsic :: iso_fortran_env
      implicit none
      integer, intent(in)  :: mode
      real_dp, intent(in)  :: r1,r2
      real_dp, intent(out) :: arr(2)
      if(mode==0) then ; arr = [ r1,r2 ]
      else             ; arr = [ min(arr(1),r1),max(arr(2),r2) ]
      endif
      end subroutine twoval

      end

# if 0

c -----------------

c
c     Old version: Kept for testing (can be run with -o=2 if enabled above)
c
c     MPI: It is assumed that there are one or more processes with the same parameters b1low, b1up, b2low, b2up
c          but different, non-overlapping kpt1. In other words, the range of kpoints kpt1 has been distributed
c          over processes, and the index array kpt1 only contains the portion for each process. Furthermore,
c          it is assumed that the portions are continuous, i.e., 1-4 and 5-8 and 9-11 etc.
c          BUT: The array kpt1p of parent kpoints is defined for the global kpt1 (with respective 'offsets')!
c          This is different from the related routine in gauss.f.
c          Note that all processes of the current communicator must call this routine
c          because it contains a split-communicator call!
c
      subroutine tetrahedron5_init_old(wintgr5,frq,nfrq,ikpt,kpt1,kpt1p,nkpt1,kflag,s1,s2,b1low,b1up,b2low,b2up)

      use global
      Mpi( use Mwrapper )
      use, intrinsic :: iso_fortran_env
      implicit none
# include "interface/cuttetra.inc"
      integer, intent(in)  :: s1,s2,b1low,b1up,b2low,b2up,ikpt,nfrq,nkpt1,kflag
      integer, intent(in)  :: kpt1(nkpt1),kpt1p(nkpt)
      real_dp, intent(in)  :: frq(nfrq)
      real_dp, intent(out) :: wintgr5(b2low:b2up,b1low:b1up,nkpt1,max(1,nfrq))
      integer              :: tetra(4,6)
      real_dp              :: etetra1(4),etetra2(4),dtetra(4)
      real_dp              :: etetra1_(4,3),etetra2_(4),etetra1__(4),etetra2__(4)
      real_dp              :: w0(4),w1(4),frq1,emax1,emin2,dmin,dmax,rdum,a
      real_dp              :: wtetra0(4,3),wtetra1(4,3),wtetra2(4,3),tetra2(4,4,3)
      real_dp              :: vtetra0(3),vtetra1(3),vtetra2(3),tetra1(4,4,3),mat(4,4)
      logical              :: allocc,allunocc,kdone(nkpt)
      integer              :: kpt2(nkpt),kcube(8),ktetra(4),ktetra1(4)
      integer              :: i,ifrq,k,k1,k2,k3,kk,n,n1,n2,nn!,isym
      integer              :: itetra,itetra1,itetra2,ntetra0,ntetra1,ntetra2
      integer              :: d1low(nkpt),d1up(nkpt),d2low(nkpt),d2up(nkpt)
      integer              :: mind1,maxd1,mind2,maxd2
      integer              :: kptsum,isym
# ifdef MPI
      integer              :: ipar
      integer, allocatable :: offset(:)
      real_dp, allocatable :: Mwintgr5(:,:,:)
# endif
      real cputime

      if(max(b1up,b2up)>maxeband) Bug('Upper band index larger than maxeband.')

      call get_tetra(tetra,n,1)

      ! Check whether frq(:) is in ascending order and calculate factor a (see above).
      rdum = 0
      do ifrq = 2,nfrq
        if(frq(ifrq)<=frq(ifrq-1)) Error('input frequency mesh not properly ordered.')
        if(ifrq/=nfrq) then
          rdum = rdum  + (frq(ifrq)-frq(ifrq-1)) / (frq(ifrq+1)-frq(ifrq-1))
          if(ifrq>2.and.abs(rdum/(ifrq-1)-a)>1d-10) Bug('Hilbert mesh not exponential.')
          a    = rdum / (ifrq-1)
        endif
      enddo

      ! Define kpt2 : kpt(:,ikpt) + kpt(:,i)
      do i = 1,nkpt
        kpt2(i) = kptsum(i,ikpt)
      enddo

      ! Add bands below and above to average over degenerate states
      d1low = b1low ; d1up = b1up
      d2low = b2low ; d2up = b2up
# ifndef disable_cut_average
      if(nfrq/=0) then
        do k = 1,nkpt
          if(kflag==1) then ; k1 = k ; k2 = kpt2(k)
          else              ; k2 = k ; k1 = kpt2(k)
          endif
          d1low(k) = min(b1low,deg(b1low,k1,s1))
          d2low(k) = min(b2low,deg(b2low,k2,s2))
          d1up(k)  = deg(b1up,k1,s1)
          d2up(k)  = deg(b2up,k2,s2)
          if(d1up(k)<b1up) then
            if(d1up(k)==0) then ; d1up(k) = b2up ! deg can be zero if b1up > nband(k1,s1)
            else                ; d1up(k) = deg(d1up(k),k1,s1)
            endif
          endif
          if(d2up(k)<b2up) then
            if(d2up(k)==0) then ; d2up(k) = b2up ! deg can be zero if b2up > nband(k2,s2)
            else                ; d2up(k) = deg(d2up(k),k2,s2)
            endif
          endif
          d1up(k) = min(d1up(k),maxeband)
          d2up(k) = min(d2up(k),maxeband)
        enddo
      endif
# endif

      ! Determine weights with tetrahedron method.
      if(nfrq<=0) then ! return maximum and minimum frequency in wintgr5
        if(b2up==b2low) Error('not implemented for b2low=b2up.')
        wintgr5(b2low,  b1low,1,1) =  huge(0d0)
        wintgr5(b2low+1,b1low,1,1) = -huge(0d0)
      else
        wintgr5 = 0      ! normal run
# ifdef MPI
        allocate(offset(0:Msize-1))
        Allocate_(Mwintgr5,(maxval(kpt1p),nfrq,b2low:b2up))
        beginSplit((b1low-1)*maxband+b2low)
        i             = b1up ; call Mcast(i) ; if(i/=b1up) Bug('b1up differs over processes.')
        i             = b2up ; call Mcast(i) ; if(i/=b2up) Bug('b2up differs over processes.')
        Mwintgr5      = 0
        offset        = 0
        offset(Mrank) = kpt1p(kpt1(1)) - 1
        ipar          = -1
        call Msum(offset)
# endif
      endif

      mind1 = minval(d1low)
      maxd1 = maxval(d1up)

      call cpu_time(cputime)

      do n1 = mind1,maxd1
        n = min(b1up,max(b1low,n1))

        do k3 = 0,nkpt3(3)-1
          do k2 = 0,nkpt3(2)-1
            do k1 = 0,nkpt3(1)-1

              kcube(1) = pkpt(k1  ,k2  ,k3  ) !
              kcube(2) = pkpt(k1+1,k2  ,k3  ) !
              kcube(3) = pkpt(k1  ,k2+1,k3  ) !
              kcube(4) = pkpt(k1+1,k2+1,k3  ) ! Definition of the cube corners
              kcube(5) = pkpt(k1  ,k2  ,k3+1) !
              kcube(6) = pkpt(k1+1,k2  ,k3+1) !
              kcube(7) = pkpt(k1  ,k2+1,k3+1) !
              kcube(8) = pkpt(k1+1,k2+1,k3+1) !_____

              do itetra = 1,6 ; Mpi( if(nfrq/=0) then ; McycleP(ipar) ; endif )
                ktetra = kcube(tetra(:,itetra))
                if(tetraf) then
                  do i = 1,4
                    if(all(kpt1p(ktetra(i))/=kpt1p(:ktetra(i)-1))) then ; ktetra1(i) = kpt1p(ktetra(i)) ! restrict to the first parent kpoint
                    else                                                  ; ktetra1(i) = 0
                    endif
                  enddo
                  if(all(ktetra1==0)) cycle
                else
                  ktetra1 = kpt1p(ktetra)
                endif

                if(n1<minval(d1low(ktetra)).or.n1>maxval(d1up(ktetra))) cycle
                mind2 = minval(d2low(ktetra))
                maxd2 = maxval(d2up (ktetra))

                if(kflag==1) then ; etetra1 = ene(n1,ktetra,s1)
                else              ; etetra1 = ene(n1,kpt2(ktetra),s1)
                endif
                rdum    = minval(etetra1) ; if(rdum>efermi) cycle
                emax1   = maxval(etetra1)
                if(emax1<efermi) then ; allocc = .true. ; ntetra1 = 1 ; etetra1_(:,1) = etetra1
                else                  ; allocc = .false.
                                        call cuttetra(wtetra1,vtetra1,ntetra1,etetra1,efermi,tetra1)
                                        do itetra1 = 1,ntetra1
                                          etetra1_(:,itetra1) = matmul( etetra1,tetra1(:,:,itetra1) )
                                        enddo
                endif

                do n2 = mind2,maxd2
                  nn = min(b2up,max(b2low,n2))
                  if(kflag==1) then ; etetra2 = ene(n2,kpt2(ktetra),s2)
                  else              ; etetra2 = ene(n2,ktetra,s2)
                  endif

                  do itetra1 = 1,ntetra1
                    if(allocc) then ; etetra2_ = etetra2
                    else            ; etetra2_ = matmul ( etetra2,tetra1(:,:,itetra1) )
                    endif
                    rdum  = maxval(etetra2_) ; if(rdum<efermi) cycle
                    emin2 = minval(etetra2_)
                    if(emin2>efermi) then ; allunocc = .true. ; ntetra2 = 1
                    else                  ; allunocc = .false.
                                            call cuttetra(wtetra2,vtetra2,ntetra2,-etetra2_,-efermi,tetra2)
                    endif

                    do itetra2 = 1,ntetra2

                      if(allunocc) then ; etetra1__ = etetra1_(:,itetra1)
                                          etetra2__ = etetra2_
                      else              ; etetra1__ = matmul ( etetra1_(:,itetra1),tetra2(:,:,itetra2) )
                                          etetra2__ = matmul ( etetra2_,           tetra2(:,:,itetra2) )
                      endif
                      dtetra = etetra2__ - etetra1__
                      dmin   = minval(dtetra)
                      dmax   = maxval(dtetra)

                      if(nfrq==0) then
                        wintgr5(b2low,  b1low,1,1) = min(wintgr5(b2low,  b1low,1,1),dmin)
                        wintgr5(b2low+1,b1low,1,1) = max(wintgr5(b2low+1,b1low,1,1),dmax)
                        cycle
                      endif

                      if(allocc) then
                        if(allunocc) then ; w1  = 1
                        else              ; mat = vtetra2(itetra2) * tetra2(:,:,itetra2)
                                            w1  = vtetra2(itetra2) * wtetra2(:,itetra2)
                        endif
                      else
                        if(allunocc) then ; mat = vtetra1(itetra1) * tetra1(:,:,itetra1)
                                            w1  = vtetra1(itetra1) * wtetra1(:,itetra1)
                        else              ; mat = vtetra1(itetra1) * vtetra2(itetra2) *
     &                                            matmul ( tetra1(:,:,itetra1) , tetra2(:,:,itetra2) )
                                            w1  = mat(:,1) + mat(:,2) + mat(:,3) + mat(:,4)
                        endif
                      endif

                      do ifrq = 1,nfrq
                        if(ifrq==nfrq) then ; frq1 = frq(nfrq) New( + (1-a) * ( frq(ifrq)-frq(ifrq-1) ) )
                        else                ; frq1 = frq(ifrq)      +     a * ( frq(ifrq+1)-frq(ifrq) )
                        endif
                        if     (frq1<dmin) then ; cycle
                        else if(frq1>dmax) then ; w0 = w1
                        else                    ; call cuttetra(wtetra0,vtetra0,ntetra0,dtetra,frq1)
                                                  if(allocc.and.allunocc) then
                                                    w0 = matmul ( wtetra0(:,:ntetra0),vtetra0(:ntetra0) )
                                                  else
                                                    w0 = matmul ( mat ,
     &                                                   matmul ( wtetra0(:,:ntetra0),vtetra0(:ntetra0) ) )
                                                  endif
                        endif
                        do i = 1,4
                          k  = ktetra1(i) ; if(k==0) cycle
                          kk = ktetra(i)  ; if(n1<d1low(kk).or.n1>d1up(kk).or.n2<d2low(kk).or.n2>d2up(kk)) cycle
# ifndef MPI
                          wintgr5(nn,n,k,ifrq) = wintgr5(nn,n,k,ifrq) + w0(i)
# else
                          Mwintgr5(k,ifrq,nn)  = Mwintgr5(k,ifrq,nn)  + w0(i)
# endif
                        enddo
                      enddo

                    enddo

                  enddo

                enddo
              enddo

            enddo

          enddo
        enddo

# ifdef MPI
        if(nfrq/=0. and. (n==n1.and.n/=b1up .or. n1==maxd1)) then
          call Msum(Mwintgr5)
          do nn = b2low,b2up
            wintgr5(nn,n,:,:) = Mwintgr5(offset(Mrank)+1:offset(Mrank)+nkpt1,:,nn)
          enddo
          Mwintgr5 = 0
        endif
# endif

      enddo

c      call cpu_done(cputime)

      if(nfrq==0) then
        if(wintgr5(b2low,b1low,1,1)<-1d-12) Error('Frequency minimum negative.')
        wintgr5(b2low,b1low,1,1) = max(0d0,wintgr5(b2low,b1low,1,1))
        return
      endif

      if(tetraf) then
        do k = 1,nkpt1
          wintgr5(:,:,k,:) = wintgr5(:,:,k,:) * count(kpt1p==k Mpi(+offset(Mrank)) )
        enddo
      endif

      endSplit ; Mpi(deallocate(offset);Deallocate_(Mwintgr5))

      rdum = 4d0*6*nkpt
      do ifrq = nfrq,2,-1
        if(ifrq==nfrq) then ; frq1 = (1-a) * (frq(ifrq)-frq(ifrq-1)) New(*2)
        else                ; frq1 = 2 * a * (frq(ifrq+1)-frq(ifrq))
        endif
        wintgr5(:,:,:,ifrq) = ( wintgr5(:,:,:,ifrq) - wintgr5(:,:,:,ifrq-1) ) / (frq1*rdum)
      enddo
      wintgr5(:,:,:,1) = wintgr5(:,:,:,1) / (2*a*(frq(2)-frq(1))*rdum)

      if(frq(1)==0) wintgr5(:,:,:,1) = 0 ! for FSPEC

      if(alignbd) return ! We have to return here in order to make the residual wings small for metals. Calculation with
                         ! a small kptadd is improved, too, provided that the bands were aligned before with a very very
                         ! small kptadd.
c                         write(*,*) sum(wintgr5(2,1,:,11))
c                         return
      ! Average over degenerate eigenstates
      ! - occupied states
      do k = 1,nkpt1
        k1 = kpt1(k) ; if(kflag==2) k1 = kpt2(k1)
        n  = min(b1low,deg(b1low,k1,s1))
        do while(n<=min(b1up,nband(k1,s1)))
          nn = min(deg(n,k1,s1),maxeband)
          do ifrq = 1,nfrq
            do n2 = b2low,b2up
              rdum = sum(wintgr5(n2,max(b1low,n):min(b1up,nn),k,ifrq)) / (nn-n+1)
              wintgr5(n2,max(b1low,n):min(b1up,nn),k,ifrq) = rdum
            enddo
          enddo
          n = nn + 1
        enddo
      enddo
      ! - unoccupied states
      do k = 1,nkpt1
        k2 = kpt1(k) ; if(kflag==1) k2 = kpt2(k2)
        n  = min(b2low,deg(b2low,k2,s2))
        do while(n<=min(b2up,nband(k2,s2)))
          nn = min(deg(n,k2,s2),maxeband)
          do ifrq = 1,nfrq
            do n1 = b1low,b1up
              rdum = sum(wintgr5(max(b2low,n):min(b2up,nn),n1,k,ifrq)) / (nn-n+1)
              wintgr5(max(b2low,n):min(b2up,nn),n1,k,ifrq) = rdum
            enddo
          enddo
          n = nn + 1
        enddo
      enddo

c     Here we can average over the k points explicitly if there are equivalents in kpt1 (eg., if "switch_off_symmetry" is defined in susceptibility)
      kdone = .false.
      do k = 1,nkpt1
        if(kdone(k)) cycle
        kpt2 = 0
        n    = 0
        do isym = 1,nsym
          if(kptsym(ikpt,isym)==ikpt) then
            k1 = kptsym(kpt1(k),isym)
            do k2 = 1,nkpt1
              if(kpt1(k2)==k1) then
                n         = n + 1
                kpt2(n)   = k2
                kdone(k2) = .true.
              endif
            enddo
          endif
        enddo
        if(n>1) then
          do ifrq = 1,nfrq
            do n1 = b1low,b1up
              do n2 = b2low,b2up
                rdum = 0        ; do k1 = 1,n ; rdum = rdum + wintgr5(n2,n1,kpt2(k1),ifrq) ; enddo
                rdum = rdum / n ; do k1 = 1,n ; wintgr5(n2,n1,kpt2(k1),ifrq) = rdum        ; enddo
              enddo
            enddo
          enddo
        endif
      enddo

      end

# endif

c -----------------

c     Initializes the tetrahedron integration with a phase factor,
c     i.e., it calculates (for any f) the weights w(k,n,s) in the sum
c
c                                                 3
c     SUM wintgr(k,R) * f(k) = INT f(k) exp(ikR) d k .
c      k
c
c     A band index is not considered.
c
c     The whole idea is a bit stupid I must admit. The weights can be written as exp(ikR)*something, where "something" could be understood
c     as exp(ik(r-R)) integrated over a certain region in reciprocal space. So, the different weights differ only in a phase factor.
c
      subroutine tetrahedron6_init(wintgr_out,site)

      use global
      use, intrinsic :: iso_fortran_env
      implicit none
      complex_dp, intent(out) :: wintgr_out(nkpt)
      integer,    intent(in)  :: site(3)
      integer                 :: tetra(4,6),pnt(4)
      integer                 :: itetra,k1,k2,k3,isym,i,kk
      integer,    allocatable :: kcorn(:)
      complex_dp              :: e(4),i1,i2,i3,j1,j2,j3,cdum
      real_dp                 :: gcube(8),g(4),g21,g31,g32,g41,g42,g43

c This was an attempt to improve the k integration without much success...
# if 0
      complex_dp :: wintgr1,cdum
      real_dp    :: dot(3,3),dot1
      integer    :: ix,iy,iz,n
      integer    :: jx,jy,jz

      real_dp x,smooth,fsmooth,q
      q(x) = - 2*x**3 + 3*x**2
      smooth(x) = q(q(q(q(q(q(q(q(x))))))))

      n = 100

      do k1 = 1,3
        do k2 = 1,3
          dot(k1,k2) = dot_product ( rlat(:,k1)/nkpt3(k1) , rlat(:,k2)/nkpt3(k2) )
        enddo
      enddo

      wintgr1 = 0
      do ix = -n,n-1
        do iy = -n,n-1
          do iz = -n,n-1
            fsmooth = 1d0/nkpt/n**3
            fsmooth = fsmooth * (1-smooth(abs(ix+0.5d0)/n))
            fsmooth = fsmooth * (1-smooth(abs(iy+0.5d0)/n))
            fsmooth = fsmooth * (1-smooth(abs(iz+0.5d0)/n))
            cexp    = exp( img * 2*pi * dot_product( 1d0*[ix,iy,iz]/nkpt3 , site ) / n )
            wintgr1 = wintgr1 + cexp * fsmooth
          enddo
        enddo
      enddo
      write(*,*) site(1) , abs(wintgr1)
c      read(*,*)
      return
# endif

      call get_tetra(tetra,i,1)

      allocate ( kcorn(8) )

      wintgr_out = 0

      do k3 = 0,nkpt3(3)-1     !
        do k2 = 0,nkpt3(2)-1   ! Loop over cubes in BZ
          do k1 = 0,nkpt3(1)-1 !

            kcorn(1) = pkpt(k1  ,k2  ,k3  ) !
            kcorn(2) = pkpt(k1+1,k2  ,k3  ) !
            kcorn(3) = pkpt(k1  ,k2+1,k3  ) !
            kcorn(4) = pkpt(k1+1,k2+1,k3  ) ! Definition of the cube corners
            kcorn(5) = pkpt(k1  ,k2  ,k3+1) !
            kcorn(6) = pkpt(k1+1,k2  ,k3+1) !
            kcorn(7) = pkpt(k1  ,k2+1,k3+1) !
            kcorn(8) = pkpt(k1+1,k2+1,k3+1) !

            gcube(1) = dot_product( 2*pi/nkpt3 * [k1  ,k2  ,k3  ] , site ) !
            gcube(2) = dot_product( 2*pi/nkpt3 * [k1+1,k2  ,k3  ] , site ) !
            gcube(3) = dot_product( 2*pi/nkpt3 * [k1  ,k2+1,k3  ] , site ) !
            gcube(4) = dot_product( 2*pi/nkpt3 * [k1+1,k2+1,k3  ] , site ) ! Get phases at cube corners
            gcube(5) = dot_product( 2*pi/nkpt3 * [k1  ,k2  ,k3+1] , site ) !
            gcube(6) = dot_product( 2*pi/nkpt3 * [k1+1,k2  ,k3+1] , site ) !
            gcube(7) = dot_product( 2*pi/nkpt3 * [k1  ,k2+1,k3+1] , site ) !
            gcube(8) = dot_product( 2*pi/nkpt3 * [k1+1,k2+1,k3+1] , site ) !

            do itetra = 1,6

              ! Phases at tetrahedron corners
              g   = gcube(tetra(:,itetra))

              call rorderp(pnt,g,4) ; g = g(pnt)
              ! Exponentials at tetrahedron corners
              e   = exp( img * g )

              g21 = g(2)-g(1) ; if(g21<1d-8) then ; pnt = pnt(4:1:-1) ; g = g(4:1:-1) ; e = e(4:1:-1) ; g21 = g(2)-g(1) ; endif
              g31 = g(3)-g(1)
              g32 = g(3)-g(2)
              g41 = g(4)-g(1)
              g42 = g(4)-g(2)
              g43 = g(4)-g(3)

c              if(sum(abs(site))/=0) write(*,*) abs(g41),abs(g21),abs(g43)

              ! Case: g1=g2=g3=g4
              if(abs(g41)<1d-8) then
                do i = 1,4
                  kk             = kcorn(tetra(i,itetra))
                  wintgr_out(kk) = wintgr_out(kk) + e(1) / 6 / 4
                enddo
c                if(sum(abs(site))/=0) write(*,*) '->1',e(1)/6/4
                cycle
              ! Case g1=g2<g3=g4
              else if(abs(g21)<1d-8.and.abs(g43)<1d-8) then
                cdum = ( e(1) * (-g31**2+4*img*g31+6) + e(3) * ( 2*img*g31-6) ) / (2*g31**4)
                kk   = kcorn(tetra(pnt(1),itetra)) ; wintgr_out(kk) = wintgr_out(kk) + cdum
c                if(sum(abs(site))/=0) write(*,*) '->2',kk,'!!',abs(wintgr_out(kk))
                kk   = kcorn(tetra(pnt(2),itetra)) ; wintgr_out(kk) = wintgr_out(kk) + cdum
c                if(sum(abs(site))/=0) write(*,*) '->2',kk,'!!',abs(wintgr_out(kk))
                cdum = ( e(3) * (-g31**2-4*img*g31+6) + e(1) * (-2*img*g31-6) ) / (2*g31**4)
                kk   = kcorn(tetra(pnt(3),itetra)) ; wintgr_out(kk) = wintgr_out(kk) + cdum
c                if(sum(abs(site))/=0) write(*,*) '->2',kk,'!!',abs(wintgr_out(kk))
                kk   = kcorn(tetra(pnt(4),itetra)) ; wintgr_out(kk) = wintgr_out(kk) + cdum
c                if(sum(abs(site))/=0) write(*,*) '->2',kk,'!!',abs(wintgr_out(kk))
                cycle
              endif

              ! Other cases
              i2 = 0 ; i3 = 0
              j2 = 0 ; j3 = 0

                                i1 = -2*img*e(3) * efunc1(-img*g31) / (g21*g31*g41)
              if(abs(g32)>1d-8) i2 = -2*img*e(3) * efunc1(-img*g32) / (g21*g32*g42)
              if(abs(g43)>1d-8) i3 =  2*img*e(3) * efunc1( img*g43) / (g41*g42*g43)
                                j1 =  e(3) * ( (6-2*img*g(1)) * efunc2(-img*g31) + g(1) * g31**3/3 ) / (g21*g31*g41)
              if(abs(g32)>1d-8) j2 =  e(3) * ( (6-2*img*g(2)) * efunc2(-img*g32) + g(2) * g32**3/3 ) / (g21*g32*g42)
              if(abs(g43)>1d-8) j3 = -e(3) * ( (6-2*img*g(4)) * efunc2( img*g43) - g(4) * g43**3/3 ) / (g41*g42*g43)

                                 cdum =          3*i1 - (j1-i1*g(1)) * (1/g21+1/g31+1/g41)
              if(i2/=0.or.j2/=0) cdum = cdum - (   i2 - (j2-i2*g(1)) * (1/g21) )
              if(i3/=0.or.j3/=0) cdum = cdum + (      - (j3-i3*g(4)) * (1/g41) )

              kk = kcorn(tetra(pnt(1),itetra)) ; wintgr_out(kk) = wintgr_out(kk) + cdum / 6
c              if(sum(abs(site))/=0) write(*,*) '->3',kk,'!!',abs(wintgr_out(kk))

                                 cdum =                 (j1-i1*g(1)) * (1/g21)
              if(i2/=0.or.j2/=0) cdum = cdum - ( 2*i2 + (j2-i2*g(1)) * (1/g21) - (j2-i2*g(2)) * (1/g32+1/g42) )
              if(i3/=0.or.j3/=0) cdum = cdum + (      - (j3-i3*g(4)) * (1/g42) )

              kk = kcorn(tetra(pnt(2),itetra)) ; wintgr_out(kk) = wintgr_out(kk) + cdum / 6
c              if(sum(abs(site))/=0) write(*,*) '->3',kk,'!!',abs(wintgr_out(kk))

                                 cdum =                 (j1-i1*g(1)) * (1/g31)
              if(i2/=0.or.j2/=0) cdum = cdum - (        (j2-i2*g(2)) * (1/g32) )
              if(i3/=0.or.j3/=0) cdum = cdum + (      - (j3-i3*g(4)) * (1/g43) )

              kk = kcorn(tetra(pnt(3),itetra)) ; wintgr_out(kk) = wintgr_out(kk) + cdum / 6
c              if(sum(abs(site))/=0) write(*,*) '->3',kk,'!!',abs(wintgr_out(kk))

                                 cdum =                 (j1-i1*g(1)) * (1/g41)
              if(i2/=0.or.j2/=0) cdum = cdum - (        (j2-i2*g(2)) * (1/g42) )
              if(i3/=0.or.j3/=0) cdum = cdum + ( 3*i3 + (j3-i3*g(4)) * (1/g41+1/g42+1/g43) )

              kk = kcorn(tetra(pnt(4),itetra)) ; wintgr_out(kk) = wintgr_out(kk) + cdum / 6
c              if(sum(abs(site))/=0) write(*,*) '->3',kk,'!!',abs(wintgr_out(kk))

            enddo

          enddo
        enddo
      enddo

      wintgr_out = wintgr_out / nkpt ! take average over nkpt cubes

      deallocate ( kcorn )

      return

c Simple averaging over symmetry-equivalent k points does not make sense here because of R that would break the symmetry.

      ! Average over symmetry-equivalent k points
      allocate( kcorn(nsym) )
      do k1 = 1,nkpt
        do isym = 1,nsym ; if(any(matmul(sym(isym)%rot,site)-site/=0)) cycle
          k2 = kptsym(k1,isym)
          write(*,*) k1,k2,wintgr_out(k2)
        enddo
        read(*,*)
      enddo


# if 0
      n1 = 1
      n2 = nkpti
      do k1 = n1,n2
        n     = 0
        kcorn = 0
        do k2 = n1,n1+nkpt-1
          if(kptp(k2)==k1) then
            n        = n + 1
            kcorn(n) = mod(k2-1,nkpt)+1
          endif
        enddo
        wintgr1 = 0
        do i = 1,n
          wintgr1 = wintgr1 + wintgr_out(kcorn(i))
        enddo
        wintgr1 = wintgr1 / n
        do i = 1,n
          wintgr_out(kcorn(i)) = wintgr1
        enddo
      enddo
      deallocate ( kcorn )
# endif

      contains

      function efunc1(z)
      implicit none
      complex_dp             :: efunc1
      complex_dp, intent(in) :: z
      if(abs(z)>1d-4) then
        efunc1 = exp(z) - (1+z+z**2/2)
      else
        efunc1 = z**3/6 + z**4/24 + z**5/120
      endif
      end function efunc1

      function efunc2(z)
      implicit none
      complex_dp             :: efunc2
      complex_dp, intent(in) :: z
      if(abs(z)>1d-2) then
        efunc2 = exp(z) - (1+z+z**2/2+z**3/6)
      else
        efunc2 = z**4/24 + z**5/120 + z**6/720 + z**7/5040
      endif
      end function efunc2

      end

c -----------------

c     Returns tetrahedron weight function for integration of nonlinear function.
c
c     The weight function is returned for the state (iband,ikpt,ispin)
c     in the form of coefficients pfrq(0:3,nfrq-1), which define cubic polynomials
c     in the nfrq-1 intervals of the mesh frq(:nfrq). The mesh is defined relative
c     to efermi+ene0 (usually, ene(iband,ikpt,ispin)).
c
      subroutine tetrahedron_nonlin(nfrq,frq,pfrq,ene0,iband,ikpt,ispin,lchk)

      use global
      use util, only: chr
      use, intrinsic :: iso_fortran_env
      implicit none
      logical, intent(in)  :: lchk
      integer, intent(in)  :: iband,ikpt,ispin
      integer, intent(out) :: nfrq
      real_dp, intent(in)  :: ene0
      real_dp, intent(out) :: frq(27),pfrq(0:3,26)
      real_dp              :: ecube(8),etetra(4)
      real_dp              :: e21,e31,e41,e32,e42,e43,denom,rdum
      integer, parameter   :: mode    = 1                   ! get_tetra mode (1,2,3)
      integer, parameter   :: dintval = 72+((mode+2)/4)*216 ! 72 (mode=1), 288 (mode=2,3)
      real_dp              :: eintval(2,dintval),pintval(0:3,dintval),a(0:2),b(0:1)
      integer, allocatable :: pnt(:)
      integer              :: tetra(4,6*4),ntetra,kcorn(8)
      integer, pointer     :: pkpt0(:,:,:)
      integer              :: k0(3),k1,k2,k3
      integer              :: nintval
      integer              :: icube,itetra
      integer              :: i,j,j1,j2

      if(ikpt>nkpt) then ; pkpt0 => pkpt1 ; k0 = nint( (kpt(:,ikpt)-kptadd) * nkpt3 )
      else               ; pkpt0 => pkpt  ; k0 = nint(  kpt(:,ikpt)         * nkpt3 )
      endif

      call get_tetra(tetra,ntetra,mode)

      allocate(pnt(4))

      nintval = 0

      if(lchk) then
        write(*,*) 'ntetra',ntetra
        do i  = -1,1
        do j  = -1,1
        do j1 = -1,1
          k1 = k0(1) + i  ; if(k1<0) k1 = k1 + nkpt3(1)
          k2 = k0(2) + j  ; if(k2<0) k2 = k2 + nkpt3(2)
          k3 = k0(3) + j1 ; if(k3<0) k3 = k3 + nkpt3(3)
          kcorn(1) = pkpt0(k1,k2,k3)
          write(*,'(3I3,2X,3F10.5,F12.5,2I5)')
     &      i,j,j1,matmul(rlat,kpt(:,kcorn(1))),ene(iband,kcorn(1),ispin)-efermi,kptp(kcorn(1)),kcorn(1)
        enddo
        enddo
        enddo
        write(*,*)
      endif

      do icube = 0,7
        k1       = k0(1) - mod(icube,2)   ; if(k1<0) k1 = k1 + nkpt3(1)
        k2       = k0(2) - mod(icube/2,2) ; if(k2<0) k2 = k2 + nkpt3(2)
        k3       = k0(3) -    (icube/4)   ; if(k3<0) k3 = k3 + nkpt3(3)
        kcorn(1) = pkpt0(k1  ,k2  ,k3  ) !
        kcorn(2) = pkpt0(k1+1,k2  ,k3  ) !
        kcorn(3) = pkpt0(k1  ,k2+1,k3  ) !
        kcorn(4) = pkpt0(k1+1,k2+1,k3  ) ! Definition of the cube corners
        kcorn(5) = pkpt0(k1  ,k2  ,k3+1) !
        kcorn(6) = pkpt0(k1+1,k2  ,k3+1) !
        kcorn(7) = pkpt0(k1  ,k2+1,k3+1) !
        kcorn(8) = pkpt0(k1+1,k2+1,k3+1) !_____
        ecube    = ene(iband,kcorn,ispin)-(ene0+efermi) ! Get energies at cube corners (relative to efermi+ene0)

        if(kcorn(icube+1)/=ikpt) Bug('icube+1 is not ikpt.')

        if(lchk) then
          write(*,*) 'icube',icube+1
          write(*,*) 'kcent',-mod(icube,2),-mod(icube/2,2),icube/4
          do i = 1,8
            write(*,'(3I5,3F10.5,2X,3F10.5,F15.5)') i,kcorn(i),kptp(kcorn(i)),kpt(:,kcorn(i)),matmul(rlat,kpt(:,kcorn(i))),ecube(i)
          enddo
        endif

        do itetra = 1,ntetra

          if(all(tetra(:,itetra)/=icube+1)) cycle ! Nothing to be done if tetrahedron does not contain ikpt (input k point)

          etetra = ecube(tetra(:,itetra)) ; call rorderp(pnt,etetra,4) ! Get energies at tetrahedron corners and size-order them
c          write(*,'(4F10.5,4I3,2I5)') etetra,tetra(:,itetra),icube,itetra
          etetra = etetra(pnt)

          if(etetra(2)-etetra(1)>1d-8) e21 = (etetra(2)-etetra(1))**(-1)
          if(etetra(3)-etetra(1)>1d-8) e31 = (etetra(3)-etetra(1))**(-1)
          if(etetra(4)-etetra(1)>1d-8) e41 = (etetra(4)-etetra(1))**(-1)
          if(etetra(3)-etetra(2)>1d-8) e32 = (etetra(3)-etetra(2))**(-1)
          if(etetra(4)-etetra(2)>1d-8) e42 = (etetra(4)-etetra(2))**(-1)
          if(etetra(4)-etetra(3)>1d-8) e43 = (etetra(4)-etetra(3))**(-1)

          j = nintval+1

          do i = 1,3
            if(tetra(pnt(i),itetra)==icube+1) goto 1
          enddo
          i = 4 ; if(tetra(pnt(i),itetra)/=icube+1) Bug('ikpt not element of tetrahedron.')

 1        if(etetra(2)-etetra(1)>1d-8) then
            a = [ etetra(1)**2 , -etetra(1)*2 , 1d0 ]                                        ; denom = e21*e31*e41
            select case(i)
              case(1) ; b = [ etetra(2)*e21 + etetra(3)*e31 + etetra(4)*e41 , -e21-e31-e41 ] ; call intval(etetra(1:2),.true.)
              case(2) ; b = [-etetra(1)*e21 , e21 ]                                          ; call intval(etetra(1:2),.true.)
              case(3) ; b = [-etetra(1)*e31 , e31 ]                                          ; call intval(etetra(1:2),.true.)
              case(4) ; b = [-etetra(1)*e41 , e41 ]                                          ; call intval(etetra(1:2),.true.)
            end select
          endif
          if(etetra(3)-etetra(2)>1d-8) then
            a = [-etetra(1)*etetra(4) , etetra(1)+etetra(4) , -1d0 ]                         ; denom = e31*e41*e42
            select case(i)
              case(1) ; b = [ etetra(3)*e31 + etetra(4)*e41 , -e31-e41 ]                     ; call intval(etetra(2:3),.true.)
              case(2) ; b = [ etetra(4)*e42 , -e42 ]                                         ; call intval(etetra(2:3),.true.)
              case(3) ; b = [-etetra(1)*e31 ,  e31 ]                                         ; call intval(etetra(2:3),.true.)
              case(4) ; b = [-etetra(1)*e41 - etetra(2)*e42 ,  e41+e42 ]                     ; call intval(etetra(2:3),.true.)
            end select
            a = [-etetra(2)*etetra(3) , etetra(2)+etetra(3) , -1d0 ]                         ; denom = e31*e32*e42
            select case(i)
              case(1) ; b = [ etetra(3)*e31 , -e31 ]                                         ; call intval(etetra(2:3),.false.)
              case(2) ; b = [ etetra(3)*e32 + etetra(4)*e42 , -e32-e42 ]                     ; call intval(etetra(2:3),.false.)
              case(3) ; b = [-etetra(1)*e31 - etetra(2)*e32 ,  e31+e32 ]                     ; call intval(etetra(2:3),.false.)
              case(4) ; b = [-etetra(2)*e42 , e42 ]                                          ; call intval(etetra(2:3),.false.)
            end select
          endif
          if(etetra(4)-etetra(3)>1d-8) then
            a = [ etetra(4)**2 , -etetra(4)*2 , 1d0 ]                                        ; denom = e41*e42*e43
            select case(i)
              case(1) ; b = [ etetra(4)*e41 , -e41 ]                                         ; call intval(etetra(3:4),.true.)
              case(2) ; b = [ etetra(4)*e42 , -e42 ]                                         ; call intval(etetra(3:4),.true.)
              case(3) ; b = [ etetra(4)*e43 , -e43 ]                                         ; call intval(etetra(3:4),.true.)
              case(4) ; b = [-etetra(1)*e41 - etetra(2)*e42 - etetra(3)*e43 , e41+e42+e43 ]  ; call intval(etetra(3:4),.true.)
            end select
          endif

          if(lchk) then
            open(102,file='/usr/users/iff_th1/cfried/delete/fort.102',status='unknown')
            write(102,'(''#'',4F13.8,I3)') etetra,i
            write(*,'(2I2,2X,4I2,F15.5'NoA) itetra,i,tetra(pnt,itetra),etetra(1)
            do i = 2,4
              if(etetra(i)-etetra(i-1)<1d-8) then ; write(*,'(''  =='',F10.5'NoA) etetra(i)
              else                                ; write(*,'(''    '',F10.5'NoA) etetra(i)
              endif
            enddo
            a(0)  = 0
            do i = j,nintval
              do j = 0,100
                rdum = (eintval(1,i)*(100-j)+eintval(2,i)*j) / 100
                write(102,*) rdum,pintval(0,i)+rdum*(pintval(1,i)+rdum*(pintval(2,i)+rdum*pintval(3,i)))
              enddo
              rdum = eintval(2,i)
              a(0) = a(0) + rdum*(pintval(0,i)+rdum*(pintval(1,i)/2+rdum*(pintval(2,i)/3+rdum*pintval(3,i)/4)))
              rdum = eintval(1,i)
              a(0) = a(0) - rdum*(pintval(0,i)+rdum*(pintval(1,i)/2+rdum*(pintval(2,i)/3+rdum*pintval(3,i)/4)))
            enddo
            close(102)
            write(*,'(F15.6'NoA) a(0)
            read(*,*)
          endif

        enddo ! itetra
      enddo   ! icube

      nullify(pkpt0)

      if(nintval>dintval) Error('Increase dintval dimension to at least '//chr(nintval))

      deallocate(pnt)
      allocate(pnt(2*nintval))

      call rorderpf(pnt,eintval,2*nintval)
      nfrq = 0
      rdum = -huge(0d0)/2
      do i = 1,2*nintval
        j  = pnt(i) - 1
        j1 = mod(j,2)+1
        j2 = j/2+1
        if(eintval(j1,j2)-rdum>1d-8) then
          nfrq      = nfrq + 1
          frq(nfrq) = eintval(j1,j2)
          rdum      = eintval(j1,j2)
        endif
      enddo
      deallocate(pnt)

      pfrq(:,:nfrq-1) = 0
      do i = 1,nintval
        do j = 1,nfrq-1
          if(frq(j)+1d-8>=eintval(1,i).and.frq(j+1)-1d-8<=eintval(2,i)) pfrq(:,j) = pfrq(:,j) + pintval(:,i)
        enddo
      enddo
      pfrq(:,:nfrq-1) = pfrq(:,:nfrq-1) / (ntetra*nkpt)

      contains

      subroutine intval(etet,lnew)
      implicit none
      logical, intent(in) :: lnew
      real_dp, intent(in) :: etet(2)
      real_dp             :: p(0:3)
      p = [ a(0)*b(0) , a(1)*b(0)+a(0)*b(1) , a(1)*b(1)+a(2)*b(0) , a(2)*b(1) ] * denom
      if(any(abs(p)>1d9)) then ! Avoid large coefficients to improve numerical stability
        p(0)  = p(0) + p(1)*(etet(2)+etet(1))/2 + p(2)*(etet(2)**2+etet(2)*etet(1)+etet(1)**2)/3 +
     &                 p(3)*(etet(2)**3+etet(2)**2*etet(1)+etet(2)*etet(1)**2+etet(1)**3)/4
        p(1:) = 0
      endif
      if(lnew) then
        nintval = nintval + 1
        if(nintval<=dintval) then
          eintval(:,nintval) = etet
          pintval(:,nintval) = p
        endif
      else if(nintval<=dintval) then
        pintval(:,nintval)   = pintval(:,nintval) + p
      endif
      end subroutine intval

      end

c -----------------

c     Returns tetrahedron weight of a single state (iband,ikpt,ispin) with respect to an energy cut at ecut (relative to efermi).
c     (Same as tetrahedron_init for a single state if ecut=0)
c
      function tetrahedron_weight(iband,ikpt,ispin,ecut)
      use global
      implicit none
# include "interface/cuttetra.inc"
      real_dp             :: tetrahedron_weight
      real_dp, intent(in) :: ecut
      integer, intent(in) :: iband,ikpt,ispin
      integer, pointer    :: pkpt0(:,:,:)
      real_dp             :: ecube(8),etetra(4),weight
      real_dp             :: wtetra1(4,3),vtetra1(3)
      integer             :: ikpt1,k0(3),k1,k2,k3
      integer             :: tetra(4,6*4),ntetra,ntetra1
      integer             :: icube,itetra,kcorn(8)
      integer             :: deg1,deg2,iband1
      integer             :: i,j,k
      logical             :: below

      if(iband>maxeband) then
        if(ene(iband,ikpt,ispin)-efermi<ecut) then ; tetrahedron_weight = 1d0/nkpt
        else                                       ; tetrahedron_weight = 0
        endif
        return
      endif

      if(ikpt>nkpt) then ; pkpt0 => pkpt1
      else               ; pkpt0 => pkpt
      endif

      call getdeg(deg1,deg2,iband,ikpt,ispin)
      deg2 = min(deg2,maxeband)

      ! Catch trivial case: all adjacent cubes below or above ecut
      below = ene(iband,ikpt,ispin)-efermi < ecut
      if(ikpt>nkpt) then ; k0 = nint( (kpt(:,ikpt)-kptadd) * nkpt3 )
      else               ; k0 = nint(  kpt(:,ikpt)         * nkpt3 )
      endif
      do k = -1,1 ; k3 = k0(3) + k ; if(k3<0) k3 = k3 + nkpt3(3)
      do j = -1,1 ; k2 = k0(2) + j ; if(k2<0) k2 = k2 + nkpt3(2)
      do i = -1,1 ; k1 = k0(1) + i ; if(k1<0) k1 = k1 + nkpt3(1) ; if(all([i,j,k]==0)) cycle
        ikpt1 = pkpt0(k1,k2,k3)
        if( any( (ene(deg1:deg2,ikpt1,ispin)-efermi<ecut) .neqv. below ) ) goto 10
      enddo
      enddo
      enddo
      if(below) then ; tetrahedron_weight = 1d0/nkpt
      else           ; tetrahedron_weight = 0
      endif
      return

      ! Non-trivial case
 10   call get_tetra(tetra,ntetra,1)

      weight = 0

      do ikpt1 = 1,size(kpt,2) ; if(kptp(ikpt1)/=kptp(ikpt)) cycle
        if(ikpt1>nkpt) then ; k0 = nint( (kpt(:,ikpt1)-kptadd) * nkpt3 )
        else                ; k0 = nint(  kpt(:,ikpt1)         * nkpt3 )
        endif
        do iband1 = deg1,deg2
          do icube = 0,7
            k1       = k0(1) - mod(icube,2)   ; if(k1<0) k1 = k1 + nkpt3(1)
            k2       = k0(2) - mod(icube/2,2) ; if(k2<0) k2 = k2 + nkpt3(2)
            k3       = k0(3) -    (icube/4)   ; if(k3<0) k3 = k3 + nkpt3(3)
            kcorn(1) = pkpt0(k1  ,k2  ,k3  ) !
            kcorn(2) = pkpt0(k1+1,k2  ,k3  ) !
            kcorn(3) = pkpt0(k1  ,k2+1,k3  ) !
            kcorn(4) = pkpt0(k1+1,k2+1,k3  ) ! Definition of the cube corners
            kcorn(5) = pkpt0(k1  ,k2  ,k3+1) !
            kcorn(6) = pkpt0(k1+1,k2  ,k3+1) !
            kcorn(7) = pkpt0(k1  ,k2+1,k3+1) !
            kcorn(8) = pkpt0(k1+1,k2+1,k3+1) !_____
            ecube    = ene(iband1,kcorn,ispin) - efermi ! Get energies at cube corners (relative to efermi)
            if(kcorn(icube+1)/=ikpt1) Bug('icube+1 is not ikpt1.')

            do itetra = 1,ntetra

              if(all(tetra(:,itetra)/=icube+1)) cycle ! Nothing to be done if tetrahedron does not contain ikpt1

              etetra = ecube(tetra(:,itetra))
              if(all(etetra>=ecut-1d-10)) cycle
              if(all(etetra<=ecut+1d-10)) then
                weight = weight + 1
                cycle
              endif

              call cuttetra(wtetra1,vtetra1,ntetra1,etetra,ecut)

              do i = 1,4
                if(tetra(i,itetra)==icube+1) weight = weight + sum ( vtetra1(:ntetra1) * wtetra1(i,:ntetra1) ) ! condition true only once
              enddo

            enddo ! itetra
          enddo   ! icube
        enddo     ! iband1
      enddo       ! ikpt1

      nullify(pkpt0)

      tetrahedron_weight = weight / 4 / ntetra / nkpt / (deg2-deg1+1) / count(kptp==kptp(ikpt))

      end

