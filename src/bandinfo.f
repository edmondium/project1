c--------------------------------------------------------------------------------
c Copyright (c) 2021 Peter Grünberg Institut, Forschungszentrum Jülich, Germany
c This file is part of SPEX and available as free software under the conditions
c of the MIT license as expressed in the LICENSE file in more detail.
c--------------------------------------------------------------------------------

# include "cppmacro.h"

c -----------

c
c     Write projection table (and, if possible, PDOS) according to orbitals defined in charr(:)      
c
c     Optional arguments:
c       spin                         : spin index (1 or 2, or 0 for both spins)
c                                      magnetic systems     : projection onto both spin channels, otherwise projection only on "spin"
c                                      non-magnetic systems : output is identical for spin = 0 and 1 (l_soc: projection is spin-summed),
c                                                             spin = 2 gives an error.
c       cmtin(:maxlmindx,:ncent,:,:) : MT coefficients
c                                      third  dimension : bands
c                                      fourth dimension : spin index, must be
c                                        1      for spin = 1 and 2
c                                        nspin2 for spin = 0
c       band(:)       (with cmtin)   : alternative band indices for output (instead of 1,2,3,...)
c       result(:,:,:) (with cmtin)   : results are returned in "result" instead of written to "spex.binfo"
c
c     cmtin provided     : projection is calculated for all states in cmtin and appended to "spex.binfo"
c     cmtin not provided : projection is calculated for all states in global cmt.
c                          This allows PDOS to be calculated, which is written to "spex.dos".
c                          Projections are appended to "spex.binfo" only for selected k points, namely,
c                          for k=0, the additional "+" kpoint if available, and
c                          for all kpoints in the line JOB (e.g., "JOB GW 2:...").
c
c     title : title for output to "spex.binfo"
c
c begin interface
      subroutine band_info(charr,cmtin,spin,band,result,title)
      use global
      use util
      use file
      use, intrinsic :: iso_fortran_env !inc
      implicit none
      character(*), intent(in)                         :: charr(:)
      character(*), intent(in),               optional :: title      
      integer,      intent(in),               optional :: spin,band(:)
      real_dp,      intent(out), allocatable, optional :: result(:,:,:)
      complex_dp,   intent(in),               optional :: cmtin(:,:,:,:)
c end interface
      real_dp,      allocatable :: eulerorb(:,:)
      real_dp,      allocatable :: project(:,:),project_dos(:,:,:,:),weight(:,:),dos(:,:),dos0(:,:)
      real_dp,      parameter   :: incr = 0.0025d0 ! increment
      integer,      parameter   :: ndigit = 10     ! space available per value
      character(len=ndigit)     :: str
      character(20)             :: str1
      complex_dp,   allocatable :: cmt1(:,:,:)
      integer,      allocatable :: lorb(:),lmorb(:),centorb(:)
      integer                   :: ncharr,norb,norb1
      integer                   :: iunit
      integer                   :: spn,nspn,nbnd
      integer                   :: icent,itype,icent0,itype0
      integer                   :: ikpt,iband,ispin,iorb
      integer                   :: i,nn,s
      real_dp                   :: emax,emin,efermi1
      logical                   :: loutput
      real                      :: cputime
# include "interface/orbitals.inc"

      Load( if(.not.present(cmtin)) Error('Not implemented for LOAD.') )
      
      ncharr = size(charr) ; if(ncharr==0) Error('Orbital definition missing.')
      spn    = 0           ; if(present(spin)) spn  = spin
      nspn   = nspin       ; if(spn>0)         nspn = 1 ! nspn = number of projection values written

      if(present(cmtin)) then
        if(size(cmtin,1)/=maxlmindx)        Bug('Dimension 1 of cmtin is not maxlmindx.')
        if(size(cmtin,2)/=ncent)            Bug('Dimension 2 of cmtin is not cent.')
        if(size(cmtin,4)/=max(nspn,nspin3)) Bug('Dimension 4 of cmtin is not max(nspn,nspin3).')
        if(present(band)) then
          if(size(cmtin,3)/=size(band)) Bug('Band-index array has wrong size.')
        endif
      else if(present(result)) then ; Bug('result argument must not be present without cmtin.')
      endif

      ! Open file and write title
      if(.not.present(result)) then
        if(title==' ') Bug('Missing title for PROJECT output.')
        iunit = fopen('spex.binfo',status='unknown',position='append')
        nn    = 37+len_trim(title)
        write(iunit,'('//chr(nn)//'A)') ('#',i=1,nn)
        write(iunit,'(A)')  '# Orbital projections (PROJECT) - '//trim(title)//' #'
        write(iunit,'('//chr(nn)//'A)') ('#',i=1,nn)
        write(6,'(/A'NoA) 'Calculate band projections... '
      endif

      call cpu_time(cputime)

      ! Write projections
      ! - Define orbitals
      call orbitals(lorb,lmorb,centorb,eulerorb,norb,charr,ncharr)
      ! - Calculate and write values 
      allocate(project(norb,nspn))
      if(.not.present(cmtin)) allocate(cmt1(maxlmindx,ncent,max(nspn,nspin3)))
      do ikpt = 1,min(nkpt+1,size(kpt,2))
        loutput = present(cmtin)
        if(.not.loutput) then
          loutput = ikpt==1 .or. ikpt==nkpt+1
          do i = 1,njob
            if(allocated(job(i)%kpt)) then
              loutput = loutput .or. any(job(i)%kpt==ikpt)
            endif
          enddo
        endif
        loutput = loutput .and. .not.present(result)
        if(loutput) then
          ! - Write header
          if(present(cmtin)) then
            write(iunit,'(A'NoA) '#     '
          else
            write(iunit,'(/A/A/A'NoA) '# k point '//Chr(ikpt)//': ('//Chn(kpt(:,ikpt),',')//')','#','#               '
            if(spn==0) then ; do ispin = 2,nspin1 ; write(iunit,'(A'NoA) '          ' ; enddo ; endif
          endif
          call write_orbitals_header(iunit,charr,ncharr,nspn,ndigit,1) ; write(iunit,*)
          write(iunit,'(A'NoA) '# Band'
          if(.not.present(cmtin)) then
            write(iunit,'(A'NoA) '    Energy'
            if(spn==0) then ; do ispin = 2,nspin1 ; write(iunit,'(A'NoA) '          ' ; enddo ; endif
          endif
          call write_orbitals_header(iunit,charr,ncharr,nspn,ndigit,2) ; write(iunit,*)
        endif
        if(present(cmtin)) then ; nbnd = size(cmtin,3)
        else                    ; nbnd = minval(nband(ikpt,:))
        endif
        do iband = 1,nbnd
          if(present(cmtin)) then
            call get_projection(project,lmorb,centorb,eulerorb,norb,cmtin(:,:,iband,:),spn)
          else
            call wavefunction_mt(cmt1,maxlmindx,0,iband,ikpt,spn)
            call get_projection(project,lmorb,centorb,eulerorb,norb,cmt1,spn)
          endif
          do ispin = 1,nspn
            call merge_projection(project(:,ispin),norb1,norb,lorb,centorb)
          enddo
          if(.not.present(cmtin)) then
            if(.not.allocated(project_dos)) then
              allocate(project_dos(nkpt,maxband,norb1,nspn))
              project_dos = 0
            endif
            if(ikpt<=nkpt) project_dos(ikpt,iband,:,:) = project(:norb1,:)
          endif
          if(present(result)) then
            if(.not.allocated(result)) allocate(result(norb1,nspn,nbnd))
            result(:,:,iband) = project(:norb1,:)
          endif
          if(loutput) then
            if(present(band)) then ; write(iunit,'(I6'NoA) band(iband)
            else                   ; write(iunit,'(I6'NoA) iband
            endif              
            if(.not.present(cmtin)) then
              if(spn==0) then ; write(iunit,'('//Chr(nspin1)//'F10.4'NoA) ene(iband,ikpt,:)  *escale
              else            ; write(iunit,'(F10.4'NoA)                  ene(iband,ikpt,spn)*escale
              endif
            endif
            do iorb = 1,norb1
              write(iunit,'('//Chr(nspn)//'F'//Chr(ndigit)//'.4'NoA) project(iorb,:)
            enddo
            write(iunit,*)
          endif
        enddo
        if(present(cmtin)) exit
      enddo      
      if(.not.present(cmtin)) deallocate(cmt1)
      deallocate(lorb,lmorb,centorb,eulerorb)
      deallocate(project)

      if(.not.present(result)) then
        write(iunit,*)
        call fclose(iunit)
        call cpu_done(cputime)
        write(6,'(A)') 'Band projections written to "spex.binfo".'
      endif

      if(present(cmtin)) return

      ! Density of states
      write(6,'(A'NoA) 'Calculate (P)DOS... '
      call cpu_time(cputime)
      if(maxeband>maxband) Bug('maxeband>maxband')
      allocate ( weight(nkpt,maxeband),dos(nspn,0:norb1),dos0(nspn,0:norb1) )
      efermi1 = efermi
      iunit   = fopen('spex.dos',status='unknown')
      str     = ' '
      write(iunit,'(A,A)') '### (P)DOS ###','#'
      write(iunit,'(A'NoA) '#         '//str ; if(nspn==2) write(iunit,'(A'NoA) str
      call write_orbitals_header(iunit,charr,ncharr,nspn,ndigit,1) ; write(iunit,*)
      str = 'Total'
      write(iunit,'(A'NoA) '#   Energy'//adjustr(str) ; str = ' ' ; if(nspn==2) write(iunit,'(A'NoA) str
      call write_orbitals_header(iunit,charr,ncharr,nspn,ndigit,2) ; write(iunit,*)

      if(spn==0) then
        emin = minval(ene(:maxeband,:nkpt,:))
        emax = maxval(ene(:maxeband,:nkpt,:))
      else
        emin = minval(ene(:maxeband,:nkpt,spn))
        emax = maxval(ene(:maxeband,:nkpt,spn))
      endif
      nn   = (emax-emin)/incr
      dos0 = 0
      do i = 1,nn
        efermi = ( i*emax + (nn-i)*emin ) / nn
        write(iunit,'(F10.4'NoA) efermi * escale
        do s = 1,nspn
          if(spn==0) then ; ispin = s
          else            ; ispin = min(spn,nspin1)
          endif
          if(gauss(1)==0) then ; call tetrahedron_init(weight,nkpt,1,maxeband,efermi,ispin,.false.)
          else                 ; call gauss_init      (weight,nkpt,1,maxeband,efermi,ispin,.false.)
          endif
          dos(s,0) = sum(weight) / incr * 2/nspin2
          do iorb = 1,norb1
            dos(s,iorb) = sum(weight*project_dos(:,:maxeband,iorb,s)) / incr * 2/nspin2
          enddo
        enddo
        do iorb = 0,norb1
          do s = 1,nspn
            str1 = chr( dos(s,iorb)-dos0(s,iorb), 'F'//trim(chr(ndigit-1))//'.5')
            write(iunit,'('' '',A'NoA) trim(str1)
          enddo
        enddo
        write(iunit,*) 
        dos0 = dos
      enddo
      call fclose(iunit)
      call cpu_done(cputime)
      write(6,'(A)') '(P)DOS plot written to "spex.dos".'
      efermi = efermi1
      deallocate(weight,dos,dos0,project_dos)

      end

c -----------

c     Write header for table of orbital projections (two lines: iline=1,2)
c     No line break at end of line.
      subroutine write_orbitals_header(iunit,charr,ncharr,nspn,ndigit,iline)
      use global, only: atype,neq
      use util,   only: chr
      implicit none
      integer,      intent(in)  :: iunit,ncharr,ndigit,iline,nspn
      character(*), intent(in)  :: charr(ncharr)
      character(:), allocatable :: label
      character(len=ndigit)     :: str
      integer                   :: icent,icent0,itype,itype0
      integer                   :: i      
      i      = 0
      icent0 = 0
      itype0 = 0
      allocate( character(len=len(charr(1))) :: label )
      do
        i = i + 1
        call get_orbital_label(label,icent,i,charr,ncharr) ; if(label==' ') exit
        itype = atype(abs(icent))
        if(icent>0.or.-icent==sum(neq(:itype-1))+1) then
          if(iline==1) then
            str = ' '
            if(icent<0) then
              icent0 = 0
              if(itype/=itype0) then ; str = 'Type '//Chr(itype) ; itype0 = itype ; endif              
            else
              itype0 = 0
              if(icent/=icent0) then ; str = 'Atom '//Chr(icent) ; icent0 = icent ; endif
            endif
          else
            str = ' '//label ; if(len_trim(label)>=ndigit) str(ndigit:) = '*'
          endif
          write(iunit,'(A'NoA) adjustr(str)         
          if(nspn==2) then
            str = ' '
            write(iunit,'(A'NoA) str
          endif
        endif
      enddo
      deallocate(label)
      end

c -----------

c
c     Merge projection array project(:norb) to the ones for which an explicit label is given -> norb1, project(:norb1)
c     Examples:
c       - [s,p,d] gives, for two atoms in the atom type, 18 orbitals (s, px, py, pz, dxy, dyz, dxz, dx2y2, dz2
c         for each atom), but there are only three labels. So, norb = 18 and norb1 = 3. 
c       - (px,py,dxy) is not merged and norb = norb1 = 3.
c     The orbitals are defined in lorb(:norb) and centorb(:norb).
c     (This is an internal routine with no error checking. The orbital definitions in lorb/centorb must be correct.)
      subroutine merge_projection(project,norb1,norb,lorb,centorb)
      use global, only: atype,neq
      use, intrinsic :: iso_fortran_env
      implicit none
      integer, intent(in)    :: norb,lorb(norb),centorb(norb)      
      integer, intent(out)   :: norb1
      real_dp, intent(inout) :: project(norb)
      real_dp                :: p
      integer                :: itype,icent
      integer                :: i,j,k,k1,l,n
      norb1 = 0
      i     = 0
      j     = 0
      k     = 0
      k1    = 0
      do while(i<norb)
        l     = lorb(i+1)
        icent = centorb(i+1)
        itype = atype(abs(icent))
        if(l<100) then ; n = 2*l+1
        else           ; n = l/100
        endif
        p = sum(project(i+1:i+n))        
        if(icent>0.or.-icent==sum(neq(:itype-1))+1) then ! if projection on atom or first of atom type
          j          = j + 1
          project(j) = p
          if(neq(itype)>1.and.icent<0) k = k + 1 ! k -> count to number of entries for current atom type (if atom type has more than 1 atom)
        else if(icent<0) then
          k1              = k1 + 1 ! for each atom in atom type, k1 = index for current entry (counts from 1 to k)
          project(j-k+k1) = project(j-k+k1) + p ! rewind to proper index for current entry
          if(k1==k) then                       ! if last entry
            k1 = 0                             !   reset entry index for next atom
            if(-icent==sum(neq(:itype))) k = 0 !   if last atom, reset entry count
          endif
        endif
        i = i + n
      enddo
      norb1 = j

      end

c -----------

c
c     Returns in project(:norb,*) the projection of the state cmt1 onto the orbitals (lmorb,centorb,eulerorb).
c     spin = 0 ; return all spins in project(:norb,:nspin)
c     spin = 1 ; only return spin up in project(:norb,1)
c     spin = 2 ; only return spin down in project(:norb,1)
      subroutine get_projection(project,lmorb,centorb,eulerorb,norb,cmt1,spin)
      use global
      implicit none
      integer,    intent(in)  :: norb,lmorb(norb),centorb(norb),spin
      real_dp,    intent(in)  :: eulerorb(3,norb)
      real_dp,    intent(out) :: project(norb,*)
      complex_dp, intent(in)  :: cmt1(maxlmindx,ncent,*)
      integer,    parameter   :: MAXLORB = 4
      complex_dp              :: proj((MAXLORB+1)**2,norb)
      real_dp,    allocatable :: olap(:,:)
      real_dp                 :: rproj
      integer                 :: itype,icent,iorb
      integer                 :: i,j,k,l,m,n,lm,s,s1
      real_dp                 :: intgrf

      if(spin<0.or.spin>nspin) Error('ispin out of range.')
      if(spin==0) then ; project(:,:nspin) = 0
      else             ; project(:,:1)     = 0
      endif
      
      proj = 0
      do iorb = 1,norb
        lm           = abs(lmorb(iorb)) ; if(lm<100.and.int(sqrt(dble(lm))+1d-12)>MAXLORB) Bug('l exceeds MAXLORB.')
        proj(:,iorb) = 0
        do l = 0,MAXLORB
          call add_proj(proj(l**2+1,iorb),l,lmorb(iorb),eulerorb(:,iorb))
        enddo
      enddo

      do itype = 1,ntype
        if(all(atype(abs(centorb))/=itype)) cycle
        i = 0
        do l = 0,min(MAXLORB,lcut(itype))
          n = nindx(l,itype)
          if(any(proj(l**2+1:(l+1)**2,:)/=0)) then
            allocate ( olap(n,n) )
            s1 = 0
            do s = 1,nspin ; if(spin>0.and.spin/=s) cycle
              s1 = s1 + 1
              do j = 1,n
                do k = 1,n
                  olap(j,k) = intgrf( bas1(:,j,l,itype,s) * bas1(:,k,l,itype,s) +
     &                                bas2(:,j,l,itype,s) * bas2(:,k,l,itype,s) , itype )
                enddo
              enddo
              lm = l**2
              do m = -l,l
                lm = lm + 1
                j  = i + (m+l) * n
                do iorb = 1,norb
                  icent = abs(centorb(iorb))
                  if(itype==atype(icent)) then
                    if(proj(lm,iorb)/=0) then
                      rproj            = real(proj(lm,iorb))**2 + imag(proj(lm,iorb))**2
                      project(iorb,s1) = project(iorb,s1) + rproj *
     &                                   dot_product (       cmt1(j+1:j+n,icent,s1) ,
     &                                   matmul (     olap , cmt1(j+1:j+n,icent,s1) ))
                      if(l_soc.and.nspin==1) then
                        project(iorb,s1) = project(iorb,s1) + rproj *
     &                                     dot_product (       cmt1(j+1:j+n,icent,2) ,
     &                                     matmul (     olap , cmt1(j+1:j+n,icent,2) ))
                      endif
                    endif
                  endif
                enddo
              enddo
            enddo
            deallocate(olap)
          endif
          i = i + n*(2*l+1)
        enddo ! l
      enddo ! itype
      
      end

c -----------

c     Old band_info version.      
c      
      subroutine band_info_old

      use global
      use wrapper
      use util
      use file
      use, intrinsic :: iso_fortran_env

      implicit none

      real_dp, allocatable :: rarr(:),proj(:,:,:,:,:),olap(:,:),weight(:,:),dos(:,:,:,:)
      real_dp, parameter   :: incr = 0.0025d0 ! increment
      real_dp              :: rdum,rdum1,rdum2,acpw,emin,emax
      integer              :: iunit
      integer              :: ispin,icent,itype,ieq,iband,ikpt,ikindx
      integer              :: i,j,k,l,ll,m,n,nn
      real_dp              :: intgrf
      
      ! Orbital projections
      write(6,'(//A)') 'Band info:'
      write(6,'(A)')   '  decomposition into s,p,d,f,g character for each atom type and interstitial contribution'      
      allocate ( rarr(0:maxlcut) )
      if(lkptadd) then ; allocate ( proj(nkpt+1,maxband,nspin,0:maxlcut,ntype) )
      else             ; allocate ( proj(nkpt,  maxband,nspin,0:maxlcut,ntype) )
      endif
      proj = 0
      do ispin  = 1,nspin1
        icent = 0
        do itype = 1,ntype
          do ieq = 1,neq(itype)
            icent = icent + 1
            i     = 0
            do l = 0,lcut(itype)
              n = nindx(l,itype)
              allocate ( olap(n,n) )
              do j = 1,n
                do k = 1,n
                  olap(j,k) = intgrf( bas1(:,j,l,itype,ispin) * bas1(:,k,l,itype,ispin) +
     &                                bas2(:,j,l,itype,ispin) * bas2(:,k,l,itype,ispin) , itype )
                enddo
              enddo
              do m = -l,l
                do ikpt = 1,nkpt+1 ; if(.not.lkptadd.and.ikpt==nkpt+1) exit
                  ikindx = kindx(ikpt)
                  do iband = 1,nband(ikpt,ispin)
                    proj(ikpt,iband,ispin,l,itype) = proj(ikpt,iband,ispin,l,itype) +
     &                                                   dotprod ( cmt(i+1:i+n,icent,iband,ikindx,ispin),
     &                                               matmul ( olap,cmt(i+1:i+n,icent,iband,ikindx,ispin)))
                    if(l_soc) then
                      proj(ikpt,iband,ispin,l,itype) = proj(ikpt,iband,ispin,l,itype) +
     &                                                     dotprod ( cmt(i+1:i+n,icent,iband,ikindx,2),
     &                                                 matmul ( olap,cmt(i+1:i+n,icent,iband,ikindx,2)))
                    endif
                  enddo
                enddo
                i = i + n
              enddo
              deallocate ( olap )
            enddo
          enddo
        enddo
      enddo
      do ispin = 1,nspin1
        if(nspin1==2) then
          if(ispin==1) then ; write(6,'(/A)') '  Spin up'
          else              ; write(6,'(/A)') '  Spin down'
          endif
        endif
        do ikpt = 1,nkpt+1
          if(ikpt<=nkpti.or.ikpt==nkpt+1.and.lkptadd) then
            ikindx = kindx(ikpt)
            if(ikpt==nkpt+1) then ; write(6,'(/A)')    '  K point    0  (additional k point)'
            else                  ; write(6,'(/A,I5)') '  K point',ikpt
            endif
            write(6,'(''           '''NoA)
            do itype = 1,ntype
              write(6,'(A'NoA) '      s      p      d      f      g    rest'
            enddo
            write(6,'(A)') '      PW    av|G|'
            do iband = 1,nband(ikpt,ispin)
              rdum = 0
              write(6,'(I4,F9.5,'' '''NoA) iband,ene(iband,ikpt,ispin)
              do itype = 1,ntype
                ll   = min(4,lcut(itype))
                rarr = proj(ikpt,iband,ispin,:,itype)
                write(6,'('//chr(ll+2)//'F7.4,'' '''NoA) (rarr(l),l=0,ll),sum(rarr)-sum(rarr(:ll))
                rdum = rdum + sum(rarr)
              enddo
              rdum1 = 0
              rdum2 = 0
              do i = 1,ngpt(ikpt)
                acpw  = abs(cpw(i,iband,ikindx,ispin))**2
                rdum1 = rdum1 + acpw * sum(matmul(rlat,gpt(:,pgpt(i,ikpt)))**2)
                rdum2 = rdum2 + acpw
              enddo
              write(6,'(2F8.4)') 1-rdum,sqrt(rdum1/rdum2)
            enddo
          endif
        enddo
      enddo
      deallocate ( rarr )
      
      ! Density of states
      rdum1 = efermi
      iunit = fopen('spex.dos',status='new',numbered=.true.)
      allocate ( rarr(0:2+nspin*ntype*5),weight(nkpt,maxeband),dos(-2:3,ntype,nspin,2) )
      write(iunit,'(A)')   '### (P)DOS ###'
      write(iunit,'(A'NoA) '#                                                                   spin up'
      if(nspin==2) write(iunit,'('//chr(15*(5*ntype-1))//'X,A'NoA) '      spin down'
      write(iunit,'(/A'NoA) '#        energy      total DOS    spin-up DOS  spin-down DOS'
      do ispin = 1,nspin
        do itype = 1,ntype
          write(iunit,'(A,I3'NoA) '        Atom',itype
          do l = 0,3
            write(iunit,'(14X,A'NoA) lchar(l)
          enddo
        enddo
      enddo
      write(iunit,*)
      emin = minval(ene(:maxeband,:nkpt,:))
      emax = maxval(ene(:maxeband,:nkpt,:))
      nn   = (emax-emin)/incr
      dos  = 0
      do i = 1,nn
        rarr         = 0
        efermi       = ( i*emax + (nn-i)*emin ) / nn
        dos(:,:,:,2) = dos(:,:,:,1)
        do ispin = 1,nspin1
          if(gauss(1)==0) then ; call tetrahedron_init(weight,nkpt,1,maxeband,efermi,ispin,.false.)
          else                 ; call gauss_init      (weight,nkpt,1,maxeband,efermi,ispin,.false.)
          endif
          dos(-2,1,ispin,1)    = sum(weight) / incr
          rarr(0)              = rarr(0) + ( dos(-2,1,ispin,1) - dos(-2,1,ispin,2) ) * 2/nspin ! total
          rarr(ispin)          =             dos(-2,1,ispin,1) - dos(-2,1,ispin,2)
          if(nspin==1) rarr(2) = rarr(1)
          do itype = 1,ntype
            ll = min(3,lcut(itype))
            dos(-1,itype,ispin,1) = 0
            do l = 0,ll
              dos( l,itype,ispin,1) = sum(proj(:nkpt,:maxeband,ispin,l,itype)*weight) * 2/nspin / incr
              dos(-1,itype,ispin,1) = dos(-1,itype,ispin,1) + dos(l,itype,ispin,1)
            enddo
            do l = -1,ll
              rarr(4+(ispin-1)*ntype*5+(itype-1)*5+l) = dos(l,itype,ispin,1) - dos(l,itype,ispin,2)
            enddo
          enddo
        enddo
        write(iunit,'('//chr(4+nspin*ntype*5)//'F15.6)') efermi,rarr
      enddo
      call fclose(iunit)
      deallocate ( rarr,weight,dos,proj )
      write(6,'(A)') 'DOS plot written to spex.dos.XXX.'
      efermi = rdum1

      end
