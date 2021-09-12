c--------------------------------------------------------------------------------
c Copyright (c) 2021 Peter Grünberg Institut, Forschungszentrum Jülich, Germany
c This file is part of SPEX and available as free software under the conditions
c of the MIT license as expressed in the LICENSE file in more detail.
c--------------------------------------------------------------------------------

# include "cppmacro.h"
      program spex_extr
      use global, only: escale0 => escale
      use, intrinsic :: iso_fortran_env
      implicit none
      integer,        parameter   :: ninter=10 ! for spline-interpolated band structure
      integer                     :: i,k,i1,i2
      integer                     :: ifile,nfile,idat,ndat,idat0
      integer                     :: spin,block,ord,inter(2),job
      integer                     :: iband,ikpt,ispin,iblock
      integer                     :: ind,onull,offk,offg
      integer,        allocatable :: pnt(:),band(:),bdat(:),kdat(:),sdat(:),ldat(:),kpt(:),kpt1(:)
      integer                     :: d1,d2,d3,d4,s1,s2
      integer                     :: nkpth,ikpth,nkpoint,icol
      integer,        allocatable :: spth(:),kindx(:),spnt(:,:)
      real_dp,        allocatable :: kpth(:),frq(:)
      complex_dp,     allocatable :: dpth(:,:),dat(:),c(:),spec(:,:,:)
      complex_dp                  :: cdum
      real_dp,        parameter   :: hartree = 27.21138386d0 , bohr = 0.529189379d0
      real_dp                     :: latpar,lat(3,3),rlat(3,3)
      real_dp                     :: escale,kscale
      real_dp                     :: kpath,rdum,vec(3),null,efermi,eks
      real_dp,        allocatable :: kpoint(:,:),kvec(:,:)
      real_dp                     :: f,s
      character(2)                :: ty,ty1
      character(256)              :: line,title
      character(256), allocatable :: file(:)
      character(5)                :: chr
      logical                     :: comp,ldum,diag,vbm_tst(2),lskip,require,def,lboth,bandst,switch

      ! Default parameters
      allocate ( band(1),kpt(1) )      
      band   = 0
      kpt    = 0
      spin   = 0
      block  = 0
      ord    = 0
      null   = 0
      onull  = 0
      idat0  = 1
      inter  = 0
      icol   = 1
      comp   = .false.
      escale = 1
      kscale = 1
      bandst = .false.
      job    = -1
      switch = .false.

      ! Allocate file
      i     = 0
      line  = 'A'
      do while(line/=' ')
        i = i + 1
        call get_command_argument(i,line) ; if(line=='--help') call help
      enddo
      allocate(file(i-1))

      ! Command line
      call get_command_argument(1,ty)
      if(ty==' ') call help()
      if(all(ty/=['k ','h ','H ','g ','gl','gd','G ','s '])) call usage()
      if(ty=='g ') ty = 'gl'
      nfile = 0
      i     = 2
 1    call get_command_argument(i,line)
      select case(line(:2))
        case('-b')   ; bandst = .true. ! band structure
        case('-o')   ; if(ord /=0)     call usage() ; ord = 1
        case('-O')   ; if(ord /=0)     call usage() ; ord = 2
        case('-c')   ; comp   = .true.
        case('-k')   ; switch = .true.
        case('-i')   ; if(.not.bandst) call usage() ; if(line(4:)==' ') then ; inter = 2
                                                      else                   ; read(line(4:),*,err=1000) inter
                                                      endif
        case('-h')   ; escale = 1d0/hartree
        case('-a')   ; kscale = 1d0/bohr
        case('-B')   ; block  = 1
        case('n=')   ; if(band(1)/=0) call usage() ; deallocate(band) ; allocate ( band(valcount(line(3:))) )
                                                     read(line(3:),*,err=1000) band                                                       
        case('k=')   ; if(kpt(1)/=0)  call usage() ; deallocate(kpt)  ; allocate ( kpt(valcount(line(3:))) )
                                                     read(line(3:),*,err=1000) kpt
        case('s=')   ; if(spin/=0)    call usage() ; read(line(3:),*,err=1000) spin        
        case('0=')   ; if     (line=='0=ef')  then ; onull = 1
                       else if(line=='0=vbm') then ; onull = 2
                       else                        ; read(line(3:),*,err=1000) null
                       endif
        case('c=')   ; read(line(3:),*,err=1000) icol
        case('j=')   ; read(line(3:),*,err=1000) job
        case(' ')    ; goto 2
        case default ; nfile       = nfile + 1
                       file(nfile) = line
      end select
      i = i + 1
      goto 1
 2    continue

      if(nfile==0) call usage()
      
      write(6,'(A'NoA) '# Command:'
      k = 0
      do i1 = 0,i-1
        call get_command_argument(i1,line)
        if(any(line==file)) then
          k = k + 1
          if(k==2) write(6,'('' ...'''NoA)
          if(k>=2) cycle
        else
          k = 0
        endif
        write(6,'(A'NoA) ' '//trim(line)
      enddo
      write(6,*)
      write(6,*)
      
c
c     Spectral function

      if(ty=='s') then
        ! determine dimensions for kpoint, spec, frq
        allocate(kvec(3,2))
        kvec    = 0
        nkpoint = 0
        ndat    = 0
        s1      =  huge(0)
        s2      = -huge(0)
        k       = 0
        do ifile = 1,nfile
          inquire(file=file(ifile),exist=ldum)
          if(.not.ldum) then
            write(0,'(A)') 'File '//trim(file(ifile))//' does not exist.'
            stop
          endif
          open(1,file=file(ifile),status='old')
          call read_spec_file(.false.)
          close(1)
        enddo
        if(ndat==0) stop 'No data found.'
        if(spin==0) then ; allocate(spec(ndat,nkpoint,s1:s2))
        else             ; allocate(spec(ndat,nkpoint,spin:spin))
        endif
        allocate(kpoint(3,nkpoint))
        allocate(kindx(nkpoint))
        allocate(frq(ndat))
        spec = 0
        frq  = huge(0d0)
        ! read data -> kpoint, spec, frq
        kvec  = 0
        ispin = 0
        k     = 0
        do ifile = 1,nfile
          open(1,file=file(ifile),status='old')
          call read_spec_file(.true.)
          close(1)          
        enddo
        call get_rlat(rlat,lat)
        if(any(frq==huge(0d0))) stop 'Frequencies not fully defined. (bug?)'
        ! write to output
        escale = escale * hartree/escale0 ! frq might be in hartree, scale to eV
        write(6,'(A)') trim(title) ; line = title(2:)
        write(6,'(A,8X,A,11X,A,A)') '#','k path','freq',adjustr(line(:30))
        kpath = 0
        do i = 0,ndat*nkpoint-1
          if(switch) then ; k    = i / ndat    + 1 ; idat = i -    (k-1) * ndat    + 1 ; ldum = k>1.and.idat==1
          else            ; idat = i / nkpoint + 1 ; k    = i - (idat-1) * nkpoint + 1 ; ldum = k>1 ; if(k==1) kpath = 0
          endif
          if(ldum) then
            call nearest_k(kpoint(:,k),kpoint(:,k-1))
            kpath = kpath + sqrt(sum(matmul(rlat,kpoint(:,k)-kpoint(:,k-1))**2))
          endif
          if(comp) then ; write(6,'(2F15.7,'//chr(size(spec,3)*2)//'F30.10'NoA) kpath*kscale,frq(idat)*escale,spec(idat,k,:)
          else          ; write(6,'(2F15.7,'//chr(size(spec,3))  //'F30.10'NoA) kpath*kscale,frq(idat)*escale,real(spec(idat,k,:))
          endif
          if(kindx(k)/=0) write(6,'(A,3F10.5'NoA) '  # k point:',kpoint(:,k)
          write(6,*)
          if(switch) then ; if(idat==ndat) write(6,'(/A)')
          else            ; if(k==nkpoint) write(6,'(/A)')
          endif
        enddo
        deallocate(kpoint,frq,spec,kvec,kindx)
        stop
      endif

c
c     Read out files

      def  = .false.
 4    idat = 0      
      
      Lfile: do ifile = 1,nfile

        ! open file
        allocate(kpt1(size(kpt)))
        line = file(ifile)
        kpt1 = kpt
        vec  = huge(rdum)
        ind  = index(line,'(')
        if(ind/=0) then
          if(size(kpt)/=1) stop 'Multiple k points selected but ''(k=..'' extension allowed only for single k point.'
          file(ifile) = line(:ind-1)
          if(line(ind+1:ind+2)=='k=') then
            i1  = index(line(ind+1:),',') ; if(i1==0) i1 = len(line)
            i2  = index(line(ind+1:),')') ; if(i2==0) call usage()
            i   = min(i1,i2)              ; if(i<=3)  call usage()
            read(line(ind+3:ind+i-1),*) kpt1
            ind = ind+i
            if(line(ind:ind+2)==',p=') then
              i   = index(line(ind+1:),')') ; if(i<=3) call usage()
              read(line(ind+3:ind+i-1),*,iostat=i) vec
              if(i/=0) call usage()
            else if(line(ind+1:)/=' ') then
              call usage()
            endif
          else
            call usage()
          endif
        endif
        
        inquire(file=file(ifile),exist=ldum)
        if(.not.ldum) then
          write(0,'(A)') 'File '//trim(file(ifile))//' does not exist.'
          stop
        endif
        open(1,file=file(ifile),status='old')

        if(onull==2) then
          null    = -huge(null)
          vbm_tst = .false.
        endif

        ! read lattice parameter
        line = ' '
        do while(line(3:19)/='Lattice parameter')
          read(1,'(A)',iostat=i) line ; if(i/=0) stop 'Reached end of file while looking for lattice parameter'
        enddo
        read(line(23:),*) latpar   ; read(1,'(A)') line
        read(line(23:),*) lat(:,1) ; read(1,'(A)') line
        read(line(23:),*) lat(:,2) ; read(1,'(A)') line
        read(line(23:),*) lat(:,3)      
        call get_rlat(rlat,lat)

        ! Read Fermi energy
        if(onull/=0) then
          do 
            read(1,'(A)',iostat=i) line
            if(i/=0) stop 'Reached end of file while looking for Fermi energy.'
            if(line(:12)=='Fermi energy') then
              read(line(36:47),*) efermi
              exit
            endif
          enddo
        endif

        ! Skip to proper job
        if(job>=0) then
          do
            read(1,'(A)',iostat=i) line
            if(i/=0) stop 'Reached end of file while looking for job number.'
            if(line(:16)=='##### Job index:') then
              read(line(17:21),*) i
              if(i==job) exit
            endif
          enddo
        endif

        ! Read list of k points
        line = ' '
        do while(line/='List of k points')
          read(1,'(A)',iostat=i) line
          if(i/=0) then
            if(def) write(0,'(A)') '# Reached end of file '//trim(file(ifile))//' while looking for list of k points.'
            close(1)
            deallocate(kpt1)
            cycle Lfile
          endif
        enddo
        nkpoint = -1
        do while(line/=' ')
          read(1,'(A)',iostat=i) line ; if(i/=0) stop 'Reached end of file while defining k-point list.'
          nkpoint = nkpoint + 1    
        enddo
        if(nkpoint==0) stop 'No k points found.'
        allocate(kpoint(3,0:nkpoint),kindx(nkpoint))
        do k = 0,nkpoint ; backspace(1) ; enddo
        do k = 1,nkpoint
          read(1,*,iostat=i) kindx(k),kpoint(:,k) ; if(i/=0) stop 'Read error while reading k-point list'
        enddo

        ! read data
 10     read(1,'(A)',end=13) line
        if(line(:9)=='### BLOCK') then
          read(line(11:15),*,err=1001) iblock
          read(line(23:27),*,err=1001) ikpt
          if(line(38:46)=='spin down') then ; ispin = 2
          else                              ; ispin = 1
          endif
        else if(line(:11)=='### K POINT') then
          iblock = 0
          if(ty=='H'.or.ty=='G') then
            write(0,'(A)') 'No full solution in file '//trim(file(ifile))
            stop
          endif
          read(line(13:17),*,err=1001) ikpt
          if(line(20:28)=='spin down') then ; ispin = 2
          else                              ; ispin = 1
          endif
        else if(job>=0.and.line(:16)=='##### Job index:') then
          goto 13
        else
          goto 10
        endif
        
        require = (any(ikpt==kpt1).or.kpt1(1)==0).and.(ispin==spin.or.spin==0)
        do while(line(:26)/='--- DIAGONAL ELEMENTS [eV]')
          read(1,'(A)',iostat=i) line ; if(i/=0) stop 'Reached end of file while looking for diagonal elements.'
        enddo
        diag = any(ty==['k ','h ','gd','gl'])
        ty1  = ty
        if(.not.diag.and.line=='--- DIAGONAL ELEMENTS [eV] (identical to full solution) ---') then
          diag = .true.
          if(ty=='G') ty1 = 'gd'
          if(ty=='H') ty1 = 'h '
        endif

        if(diag) then

          do while(line(:7)/=' Bd    ')
            read(1,'(A)',iostat=i) line ; if(i/=0) stop 'Reached end of file while looking for "Bd ..."'
          enddo
          lskip = .true.
          offk  = 44
          offg  = 64
          if(index(line,'GW')==0.and.index(line,'GT')==0) then ! no GW and GT
            lskip = .false.
            if(index(line,'HF')==0.and.index(line,'COSX')==0.and.index(line,'SX')==0) then ; offk = 14
            else                                                                           ; offk = 24
            endif
            if(ty1(:1)=='g') then ; write(0,'(A)') 'No GW data in file '//trim(file(ifile)) ; stop ; endif
          endif
          if(index(line,'GT')/=0) then ! GT
            offk  = 24
            offg  = 34
            if(ty1=='h') then ; write(0,'(A)') 'No HF data in file '//trim(file(ifile)) ; stop ; endif
          endif
 11       read(1,'(A)',iostat=i) line ; if(i/=0) stop 'Reached end of file while reading data.'
          if(line/=' ') then
            if(index(line,'Warning!')/=0) goto 11
            read(line(:3),*,err=1001) iband
            read(line(offk:offk+9),*,err=1001) eks
            if     (ty1=='k ') then ; read(line(offk   :offk+9) ,*,err=1001) rdum ; cdum = rdum ; if(lskip) read(1,*)
            else if(ty1=='h ') then ; read(line(offk+10:offk+19),*,err=1001) rdum ; cdum = rdum ; if(lskip) read(1,*)
            else if(ty1=='gl') then ; read(line(offg   :offg+9) ,*,err=1001) rdum ; cdum = rdum ; read(1,*)
                                      read(line(offg+10:offg+19),*,err=1001) rdum ; cdum = cdum + (0d0,1d0)*rdum
            else if(ty1=='gd') then ; read(1,'(A)') line ; if(line(offg:offg+9)==' ') stop 'No gd values found.'
                                      read(line(offg   :offg+9) ,*,err=1001) rdum ; cdum = rdum
                                      read(line(offg+10:offg+19),*,err=1001) rdum ; cdum = cdum + (0d0,1d0)*rdum
            endif
            if(require.and.(any(iband==band).or.band(1)==0)) then
              idat = idat + 1
              if(def) then
                call define_dat() ; dat(idat) = cdum
              endif
            endif  
            if(onull==2) then
              if(eks<=efermi) then ; vbm_tst(1) = .true. ; null = max(null,real(cdum)) 
              else                 ; vbm_tst(2) = .true.
              endif
            endif
            goto 11
          endif
          
        else

          do while(line(:13)/=' Bd        KS')
            read(1,'(A)',iostat=i) line
            if(i/=0) then
              write(0,'(A)') 'No full solution in file '//trim(file(ifile))
              stop
            endif
          enddo
          if(ty=='G') then
            if(line(52:53)/='GW') then
              write(0,'(A)') 'No GW data in file ',trim(file(ifile))
              stop
            endif
          endif
 12       read(1,'(A)',iostat=i) line ; if(i/=0) stop 'Reached end of file while looking for data.'
          if(line/=' ') then
            if(index(line,'Warning!')/=0) goto 12
            read(line(:3),*,err=1001) iband
            read(line(4:13),*,err=1001) eks
            if     (ty=='H') then ; read(line(24:33),*,err=1001) rdum ; cdum = rdum
            else if(ty=='G') then ; read(line(44:53),*,err=1001) rdum ; cdum = rdum
                                    read(line(54:63),*,err=1001) rdum ; cdum = cdum + (0d0,1d0)*rdum
            endif
            if(require.and.(any(iband==band).or.band(1)==0)) then
              idat = idat + 1
              if(def) then
                call define_dat() ; dat(idat) = cdum
              endif
            endif
            if(onull==2) then
              if(eks<=efermi) then ; vbm_tst(1) = .true. ; null = max(null,real(cdum)) 
              else                 ; vbm_tst(2) = .true.
              endif
            endif
            goto 12              
          endif
            
        endif
          
        goto 10
 13     close(1)
        
        if(def) then
          if(onull==1) then
            null = efermi
          else if(onull==2) then
            if(null==huge(null).or.any(.not.vbm_tst)) then
              null = 0 ; write(0,'(A)') '# Did not find VBM; 0=vbm ignored.'
            endif
          endif
          dat(idat0:idat) = dat(idat0:idat) - null
          idat0           = idat + 1
        endif

        deallocate(kpoint,kindx,kpt1)
        
      enddo Lfile

      if(.not.def) then
        ndat = idat
        def  = .true.
        allocate(dat(ndat),bdat(ndat),kdat(ndat),sdat(ndat),ldat(ndat))
        if(bandst) allocate(kvec(3,ndat))
        goto 4
      endif

      if(ndat>0) then
        write(0,'(A,I4,A)') '# Read successful. (',ndat,' )'
      else
        write(0,'(A,I4,A)') '# Read successful. ( no data )'
        stop
      endif

      ! Band structure: Define unique kdat -> kindx
      if(bandst) then
        allocate(kpoint(3,ndat),kindx(ndat))
        nkpoint     = 1
        kpoint(:,1) = kvec(:,1)
        do i = 2,ndat
          if(any(kvec(:,i)/=kvec(:,i-1))) then
            nkpoint           = nkpoint + 1
            kpoint(:,nkpoint) = kvec(:,i)
          endif
        enddo
        do i = 1,ndat
          do k = 1,nkpoint
            if(all(kvec(:,i)==kpoint(:,k))) kindx(i) = k
          enddo
        enddo
        deallocate(kpoint)
      endif

      ! Order
      if(ord/=0) then
        i1 = 1
 20     i2 = i1
        do while(i2/=ndat)
          if(allocated(kindx)) then
            if(kindx(i2)/=kindx(i2+1).or.sdat(i2)/=sdat(i2+1)) exit
          else
            if(kdat(i2)/=kdat(i2+1).or.sdat(i2)/=sdat(i2+1)) exit
          endif
          i2 = i2 + 1
        enddo
        allocate ( pnt(i1:i2) )
        if(ord==1) then ; call rorderp(pnt,1d0*bdat(i1:i2),i2-i1+1)
        else            ; call rorderp(pnt,real(dat(i1:i2)),i2-i1+1)
        endif
        pnt         = pnt + i1 - 1
        ldat(i1:i2) = ldat(pnt)
        dat(i1:i2)  = dat(pnt)
        kdat(i1:i2) = kdat(pnt)
        sdat(i1:i2) = sdat(pnt)
        if(bandst.and.ord==2) then
          call rorderp(pnt,1d0*bdat(i1:i2),i2-i1+1)
          pnt = pnt + i1 - 1
        endif
        bdat(i1:i2) = bdat(pnt)
        deallocate ( pnt )
        i1 = i2 + 1
        if(i2/=ndat) goto 20
      endif

c
c     Output
      ! band structure
      if(bandst) then
        if(ord==2) write(6,'(A/)') '# Band indices changed due to -O'
        lboth = .false.
        ! if both spins are present try to write both spin data in one line (lboth=.true.)
        if(any(sdat(2:ndat)/=sdat(1))) then
          allocate(spnt(minval(bdat):maxval(bdat),minval(kindx):maxval(kindx)))
          spnt  = 0
          lboth = .true.          
          do i = 1,ndat
            spnt(bdat(i),kindx(i)) = spnt(bdat(i),kindx(i)) + (3-2*sdat(i))
          enddo
          if(all(spnt==0)) then ! data spin symmetric, write in one line
            do i = 1,ndat
              if(sdat(i)==2) spnt(bdat(i),kindx(i)) = i
            enddo
          else                    ! data not spin symmetriy, write spin index explicitly
            lboth = .false.
            write(0,'(A)') '# Data not spin symmetric, spin column added.'
            deallocate(spnt)
          endif
        endif
        if(all(sdat(2:ndat)==sdat(1)).or.lboth) then
          spin = sdat(1) ; if(lboth) spin = 1
          write(6,'(A/)') '#   k path    Energy'
        else
          write(6,'(A/)') '#   k path  s    Energy'
        endif
        if(lboth) then ; allocate(c(2))
        else           ; allocate(c(1))
        endif
        do ispin = 1,2
          if(spin/=0.and.spin/=ispin) cycle
          i1 = minval(bdat(:ndat))
          i2 = maxval(bdat(:ndat))
          do iband = i1,i2
            ! Get kpath values for iband
            nkpth = count(bdat(:ndat)==iband.and.sdat(:ndat)==ispin)
            allocate(kpth(nkpth),dpth(nkpth,size(c)),spth(nkpth),kpoint(3,nkpth))
            kpoint = huge(0d0)
            ikpth  = 0
            kpath  = 0
            vec    = kvec(:,1)
            do i = 1,ndat
              if(bdat(i)==iband.and.sdat(i)==ispin) then
                if(any(vec/=kvec(:,i))) then
                  call nearest_k(kvec(:,i),vec)
                  kpath = kpath + sqrt ( sum(matmul(rlat,kvec(:,i)-vec)**2) )
                  vec   = kvec(:,i)
                endif
                ikpth         = ikpth + 1
                kpth(ikpth)   = kpath
                dpth(ikpth,1) = dat(i) ; if(lboth) dpth(ikpth,2) = dat(spnt(bdat(i),kindx(i)))
                spth(ikpth)   = sdat(i)
                if(kdat(i)/=0) kpoint(:,ikpth) = kvec(:,i)
              endif
            enddo
            ! Write out
            write(6,'(A,I4)') '#',iband
            do ikpth = 1,nkpth
              if(inter(1)/=0.and.ikpth>1) then ! interpolate
                if(kpth(ikpth)>kpth(ikpth-1)) then                  
                  do i = 1,ninter-1
                    rdum = ( kpth(ikpth-1) * (ninter-i) + kpth(ikpth) * i ) / ninter
                    call           interpolate(c(1),cdum,rdum,dpth(:,1),kpth,nkpth,inter(1),inter(2))
                    if(lboth) call interpolate(c(2),cdum,rdum,dpth(:,2),kpth,nkpth,inter(1),inter(2))
                    write(6,'(F10.5'NoA) kscale * rdum
                    if(spin==0)     write(6,'(I3'NoA)   spth(ikpth)
                    if(comp) then ; write(6,'(4F10.5)') escale * c
                    else          ; write(6,'(2F10.5)') escale * real(c)
                    endif 
                  enddo
                endif
              endif
              write(6,'(F10.5'NoA) kscale * kpth(ikpth) ! explicit
              if(spin==0)     write(6,'(I3'NoA)                              spth(ikpth)
              if(comp) then ; write(6,'('//chr(size(dpth,2))//'(2F10.5)'NoA) escale * dpth(ikpth,:)
              else          ; write(6,'('//chr(size(dpth,2))//'  F10.5' NoA) escale * real(dpth(ikpth,:))
              endif
              if(all(kpoint(:,ikpth)/=huge(0d0))) write(6,'(A,3F10.5'NoA) '  # k point:',kpoint(:,ikpth)
              write(6,*)
            enddo
            deallocate(kpth,dpth,spth,kpoint)
            if(nkpth>0) write(6,'(/)')
          enddo
        enddo
        deallocate(c)
        if(lboth) deallocate(spnt)
      else
        if(all(kdat(2:ndat)==kdat(1))) kpt   = 1
        if(all(sdat(2:ndat)==sdat(1))) spin  = 1
        if(all(bdat(2:ndat)==bdat(1))) band  = 1
        if(all(ldat(2:ndat)==ldat(1))) block = 1
        line = ' '
        write(6,'(A'NoA) '#'
        d1 = len_trim(chr(maxval(bdat(:ndat))))
        if(block==0)
     &  d2 = len_trim(chr(maxval(ldat(:ndat))))
        d3 = len_trim(chr(maxval(kdat(:ndat))))
        d4 = len_trim(chr(maxval(sdat(:ndat))))
        write(6,'(A'NoA) line(:d1)//'n'
        if(block==0)
     &  write(6,'(A'NoA) line(:d2)//'b'
        write(6,'(A'NoA) line(:d3)//'k'
        write(6,'(A'NoA) line(:d4)//'s'
        write(6,'(A)') '    Energy'
        do i = 1,ndat
          write(6,'(A'NoA) ' '                     ; write(line,*) d1+1
          write(6,'(I'//trim(line)NoA) bdat(i)
          if(block==0) then                        ; write(line,*) d2+1
            write(6,'(I'//trim(line)NoA) ldat(i) 
          endif                                    ; write(line,*) d3+1
          write(6,'(I'//trim(line)NoA) kdat(i) ; write(line,*) d4+1
          write(6,'(I'//trim(line)NoA) sdat(i)
          if(comp) then   ; write(6,'(2F10.5)') dat(i)
          else            ; write(6,'( F10.5)') escale * real(dat(i))
          endif
        enddo
      endif

      deallocate(dat,bdat,kdat,sdat,ldat,file)
      deallocate(band,kpt)
      if(allocated(kvec)) deallocate(kvec)

      stop
 1000 stop 'Syntax error in command line. Not a value.'
 1001 write(0,'(A)') 'Read error in file '//trim(file(ifile))
      write(0,'(A)') 'Offending line: '//trim(line)
      stop
      
      contains

      subroutine define_dat()
      implicit none
      integer :: i
      bdat(idat) = iband
      kdat(idat) = ikpt
      sdat(idat) = ispin
      ldat(idat) = iblock
      if(bandst) then
        if(all(kindx/=ikpt)) then
          write(0,'(A)') 'Error in file '//trim(file(ifile))
          stop 'Did not find k point in k-point list. (bug?)'
        endif
c        if(all(kindx/=ikpt)) ikpt = kindx(1) !stop 'Did not find k point in k-point list. (bug?)'        
c# warning change back!!   
        do i = 1,nkpoint
          if(kindx(i)==ikpt) kvec(:,idat) = kpoint(:,i)
        enddo        
      endif
      end subroutine define_dat

c     Counts the values (comma-separated) in string
      function valcount(string)
      implicit none
      integer                   :: valcount
      character*(*), intent(in) :: string
      integer                   :: i
      if(string==' ') stop 'Empty string (indices missing).'
      valcount = 1
      do i = 1,len_trim(string)
        if(string(i:i)==',') then
          valcount = valcount + 1
        else
          if(string(i:i)<'0'.or.string(i:i)>'9') call usage()
        endif
      enddo
      end function valcount

c     Returns in kvec1 the k vector nearest to kvec2.
      subroutine nearest_k(kvec1,kvec2)
      implicit none
      real_dp, intent(out) :: kvec1(3)
      real_dp, intent(in)  :: kvec2(3)
      real_dp              :: vec(3),vec1(3),rdum,rdum1
      integer              :: i,j,k
      vec1 = kvec1
      rdum = sum(matmul(rlat,kvec1-kvec2)**2)
      do i = -2,2 ; do j = -2,2 ; do k = -2,2
        vec   = kvec1 - [i,j,k]
        rdum1 = sum(matmul(rlat,vec-kvec2)**2)
        if(rdum1+1d-12<rdum) then
          vec1 = vec
          rdum = rdum1 
        endif
      enddo ; enddo ; enddo
      kvec1 = vec1
      end subroutine nearest_k

c     Read next line
      subroutine nextline(i,mode)
      implicit none
      integer, intent(out) :: i
      integer, intent(in)  :: mode
      read(1,'(A)',iostat=i) line
      if(mode/=0.and.i/=0) stop 'Unexpected end of file.'
      end subroutine nextline

c     Read file from unit 1 containing spectral function 
c     ldat = .false. : update dimensions (update: kvec,nkpoint,ndat,s1,s2)
c     ldat = .true.  : read data         (update: kpoint,frq,dat,kindx)
      subroutine read_spec_file(ldat)
      implicit none
      logical, intent(in) :: ldat
      real_dp             :: lat1(3,3)
      integer             :: ispin,iband,ikindx,idat,j
      real_dp             :: s,si
      ! header
      call nextline(i,1)
      if(ifile==1)         then ; title = line ; if(line(1:1)/='#') stop 'Not a spectral file: title is not commented.'
      else if(line/=title) then ; write(0,'(A)') 'Title changed:'//trim(title(2:))//' ->'//trim(line(2:)) ; stop
      endif
      call nextline(i,1) ; if(line(:10)/='# lattvec:') stop 'Not a spectral file: lattvec not found.'
      read(line(11:),*,iostat=i) lat1 ; if(i/=0)       stop 'Error while reading lattice vectors.'
      if(ifile==1)            then ; lat = lat1
      else if(any(lat/=lat1)) then ; stop 'Lattice vectors disagree.'
      endif
      call nextline(i,1) ; if(line(:10)/='# k point:') stop 'Not a spectral file: k point not found.'
      rewind(1)
      ! body
      idat = 0
      do
        call nextline(i,0) ; if(i/=0) exit
        if(line(:10)=='# k point:') then
          ikindx = 0
          ispin  = 1
          iband  = 0
          read(line(11:),*,iostat=i) kvec(:,2) ; if(i/=0) stop 'Error while reading k point'
          call nextline(i,1)
          if(line(:10)=='# k index:') then ; read(line(11:),*,iostat=i) ikindx ; if(i/=0) stop 'Error while reading k index.'
          else                             ; backspace(1)
          endif
          call nextline(i,1)
          if(line(:7)=='# spin:') then ; read(line(8:),*,iostat=i) ispin ; if(i/=0) stop 'Error while reading spin.'
          else                         ; backspace(1)
          endif
          call nextline(i,1)
          if(line(:7)=='# band:') then ; read(line(8:),*,iostat=i) iband ; if(i/=0) stop 'Error while reading band index.'
          else                         ; backspace(1)
          endif
          if(.not.ldat) then
            s1   = min(s1,ispin)
            s2   = max(s2,ispin)
            ndat = max(ndat,idat)
          endif
          idat = 0
          if((spin==0.or.ispin==spin).and.(any(kvec(:,1)/=kvec(:,2)).or.k==0).and.(iband==0.or.band(1)==0.or.any(iband==band))) then
            k         = k + 1
            kvec(:,1) = kvec(:,2) ; if(ldat) then ; kpoint(:,k) = kvec(:,2) ; kindx(k) = ikindx ; endif
          endif
        else if(line/=' '.and.line(:1)/='#') then
          if((spin==0.or.ispin==spin).and.(iband==0.or.band(1)==0.or.any(iband==band))) then ; if(k==0) stop 'k index = 0. (bug?)'
            idat = idat + 1
            if(ldat) then
              if(comp) then
                read(line,*,iostat=i) f,(s,j=1,2*icol-1),si ; if(i/=0) stop 'Error while reading data.'
                spec(idat,k,ispin) = s + (0d0,1d0)*si
              else
                read(line,*,iostat=i) f,(s,j=1,icol) ; if(i/=0) stop 'Error while reading data.'
                spec(idat,k,ispin) = s
              endif
              if(frq(idat)==huge(0d0)) then ; frq(idat) = f
              else                          ; if(f/=frq(idat)) stop 'Frequencies disagree.'
              endif
            endif
          endif
        endif
      enddo
      if(.not.ldat) then
        ndat    = max(ndat,idat)
        nkpoint = k
      endif
      end subroutine read_spec_file

      end

c     Syntax error
      subroutine usage()
      write(0,'(A)') 'Usage: spex.extr [k|h|H|g|gl|gd|G] [options] out-file1 ...'
      write(0,'(A)') '       spex.extr s [options] spectral-file1 ...'
      write(0,'(A)') '       for help type "spex.extr --help"'
      stop
      end

c     Help      
      subroutine help()
      write(6,'(A)') 'Usage: spex.extr [k|h|H|g|gl|gd|G] <k=kpoint> <s=spin> <n=bands|-b> <-i<=m,n>> <-o|-O> <-c> <-h> <-a> '//
     &                       'out-files'
      write(0,'(A)') '       spex.extr s <s=spin> <n=bands> <c=column> <-h> <-a> spectral-files'
      write(6,'(A)') '       k|h|H|gl|gd|G|s : Kohn-Sham | Hartree-Fock | Hartree-Fock(full) |',
     &               '                         GW(linearized) | GW(direct) | GW(full) | spectral function',
     &               '       g        : identical to gl'
      write(6,'(A)') '       k=kpoint : Specifies k-point indices (kpoint=1,2,...)'
      write(6,'(A)') '       s=spin   : Specifies spin index (spin=1 or 2)'
      write(6,'(A)') '       n=bands  : Specifies band indices (band=1,2,...)'
      write(6,'(A)') '       0=null   : Energies relative to null (ef=Fermi energy, vbm=valence band max)'
      write(6,'(A)') '       c=column : Specifies column'
      write(6,'(A)') '       -b       : Prepares data for band-structure plot (only out-files)'
      write(6,'(A)') '       -k       : Switches kpoint-frequency axes (only spectral-files)'
      write(6,'(A)') '       -i<=m,n> : Uses spline interpolation for band-structure plot with mth and nth derivative set'
      write(6,'(A)') '                  to zero at the left and right plot edges, respectively (defaults to m=n=1)'
      write(6,'(A)') '       -o       : Orders energies according to band index'
      write(6,'(A)') '       -O       : Orders energies according to size'
      write(6,'(A)') '       -c       : Write energies as complex numbers.'
      write(6,'(A)') '       -h       : Write energies in Hartree (eV is default).'
      write(6,'(A)') '       -a       : Write k path in 1/Angstrom (1/Bohr is default).'      
      write(6,'(A)') '       To specify a k-point index for a particular file: file"(k=5)".'
      write(6,'(A)') '       To redefine the corresponding k vector: file"(k=5,p=0.1 0 0)".'
      stop
      end

c     Orders rarr(1:n) according to size and returns a correspondingly defined pointer in pnt.
      subroutine rorderp(pnt,rarr,n)
      use, intrinsic :: iso_fortran_env
      implicit none
      integer, intent(in)  :: n
      integer, intent(out) :: pnt(n)
      real_dp, intent(in)  :: rarr(n)
      integer              :: i,j,k
      do i=1,n
        pnt(i) = i
        do j=1,i-1
          if(rarr(pnt(j))>rarr(i)) then
            do k=i,j+1,-1
              pnt(k) = pnt(k-1)
            enddo
            pnt(j) = i
            exit
          endif
        enddo
      enddo
      end      

      subroutine vprod(a,b,c)
      use, intrinsic :: iso_fortran_env
      implicit none
      real_dp, intent(out) :: a(3)
      real_dp, intent(in)  :: b(3),c(3)
      a(1) = b(2)*c(3) - b(3)*c(2)
      a(2) = b(3)*c(1) - b(1)*c(3)
      a(3) = b(1)*c(2) - b(2)*c(1)
      end

      subroutine get_rlat(rlat,lat)
      use, intrinsic :: iso_fortran_env
      implicit none
      real_dp, intent(out) :: rlat(3,3)
      real_dp, intent(in)  :: lat(3,3)
      real_dp, parameter   :: pi = 3.14159265359d0
      call vprod(rlat(:,1),lat(:,2),lat(:,3))
      call vprod(rlat(:,2),lat(:,3),lat(:,1))
      call vprod(rlat(:,3),lat(:,1),lat(:,2))
      rlat = 2*pi * rlat / dot_product(lat(:,1),rlat(:,1))
      end

c     Returns an interpolated function value and derivative (in y and dy) at x
c     of a function given by n value pairs (xx,yy) using spline interpolation
c     with vanishing gradients at the end points.      
c     x, xx : real valued ; y, yy : complex valued.
      recursive subroutine interpolate(y,dy,x,yy,xx,n,m1,m2)
      use, intrinsic :: iso_fortran_env
      implicit none
      integer,    intent(in)  :: n,m1,m2
      real_dp,    intent(in)  :: xx(n),x
      complex_dp, intent(in)  :: yy(n)
      complex_dp, intent(out) :: y,dy
      complex_dp, parameter   :: img = (0d0,1d0)
      complex_dp              :: a(n-1),b(n-1)
      real_dp                 :: u,v
      integer                 :: i
      if(xx(1)>xx(n)) then
        call interpolate(y,dy,-x,yy,-xx,n,m1,m2)
        dy = -dy
        return
      endif
      call cubicspline_c(a,b,yy, xx,n,(0d0,0d0),m1,(0d0,0d0),m2)
      if(x<=xx(1)) then
        dy = (b(1)-3*yy(1)) / (xx(2)-xx(1))
        y  = yy(1) + dy*(x-xx(1))
      else if(x>=xx(n)) then
        dy = (3*yy(n)-a(n-1)) / (xx(n)-xx(n-1))
        y  = yy(n) + dy*(x-xx(n))
      else
        i = 1
        do while(x>xx(i+1))
          i = i + 1
        enddo
        u  = (x-xx(i)) / (xx(i+1)-xx(i))
        v  = 1 - u
        y  = yy(i)*v**3 + yy(i+1)*u**3 + a(i)*u**2*v + b(i)*u*v**2
        dy = (b(i)-3*yy(i)) * v**2 + (3*yy(i+1)-a(i)) * u**2 + 2*(a(i)-b(i)) * u*v
        dy = dy / (xx(i+1)-xx(i))
      endif
      end


c     ------------

c     Returns the coefficients a(i) and b(i), i=1..n-1, of the interpolation s(x) with cubic splines of the set
c     {x(i),y(i); i=1..n}. Within the interval [x(i),x(i+1)] the interpolation is given by
c
c     s(x) = y(i)*v**3 + y(i+1)*u**3 + a(i)*u**2*v + b(i)*u*v**2
c
c     with u=(x-x(i))/(x(i+1)-x(i)) and v=1-u.
c
c     There are two free parameters, which are fixed with the last four arguments:
c
c     (d/dx)^m1 s(x(1)) = y1 and
c     (d/dx)^mn s(x(n)) = yn for m1 = 1,2,3,4 and mn = 1,2,3,
c     and (d/dw)^3 s(x(w(1))) = y1 with x(w) = w/(1+w) for m1 = 4.
c
      subroutine cubicspline(a,b,y,x,n,y1,m1,yn,mn)
      use, intrinsic :: iso_fortran_env
      implicit none
      integer, intent(in)  :: n,m1,mn
      real_dp, intent(in)  :: x(n)
      real_dp, intent(in)  :: y(n),y1,yn
      real_dp, intent(out) :: a(n-1),b(n-1)
      real_dp              :: h(n-1)
      real_dp              :: d(n),dl(n-1),du(n-1),p(n)
      integer              :: i
      if(n<=1) stop 'cubicspline: at least two points are needed.'
      do i = 1,n-1
        if(x(i+1)<=x(i)) stop 'cubicspline: x array not ordered according to size. (bug?)'
        h(i) = x(i+1) - x(i)
      enddo      
      ! n-2 tridiagonal system of linear equations
      do i = 2,n-1
        d(i)    = 4/h(i) + 4/h(i-1)
        du(i)   = 2/h(i)
        dl(i-1) =          2/h(i-1)
        p(i)    = 6 * ( (y(i)-y(i-1)) / h(i-1)**2 + (y(i+1)-y(i)) / h(i)**2 )
      enddo
      ! imposing the additional boundary conditions
      if(m1==1) then
        d(1)  = 1
        du(1) = 0
        p(1)  = y1
      else if(m1==2) then
        d(1)  = 4/h(1)
        du(1) = 2/h(1)
        p(1)  = 6 * (y(2)-y(1)) / h(1)**2 - y1
      else if(m1==3) then
        d(1)  = 6/h(1)**2
        du(1) = 6/h(1)**2
        p(1)  = 12 * (y(2)-y(1)) / h(1)**3 + y1
      else if(m1==4) then
        d(1)  = 6/h(1)**2 + 24/h(1) + 6
        du(1) = 6/h(1)**2 + 12/h(1)
        p(1)  = 12*(3*h(1)+1) * (y(2)-y(1)) / h(1)**3 + y1
      else
        stop 'cubicspline: wrong m1 argument.'
      endif
      if(mn==1) then
        d(n)    = 1
        dl(n-1) = 0
        p(n)    = yn
      else if(mn==2) then
        d(n)    = 4/h(n-1)
        dl(n-1) = 2/h(n-1)
        p(n)    = 6 * (y(n)-y(n-1)) / h(n-1)**2 + yn
      else if(mn==3) then
        d(n)    = 6/h(n-1)**2
        dl(n-1) = 6/h(n-1)**2
        p(n)    = 12 * (y(n)-y(n-1)) / h(n-1)**3 + yn
      else
        stop 'cubicspline: wrong mn argument.'
      endif
      ! solve matrix equation
      call dgtsv(n,1,dl,d,du,p,n,i)
      if(i/=0) stop 'cubicspline: dgtsv failed.'
      do i = 1,n-1
        a(i) = 3*y(i+1) - h(i)*p(i+1)
        b(i) = 3*y(i)   + h(i)*p(i)
      enddo
      end

c     ------------

c     Complex version of cubicspline.
      subroutine cubicspline_c(a,b,y,x,n,y1,m1,yn,mn)
      use, intrinsic :: iso_fortran_env
      implicit none
      integer,    intent(in)  :: n,m1,mn
      real_dp,    intent(in)  :: x(n)
      complex_dp, intent(in)  :: y(n),y1,yn
      complex_dp, intent(out) :: a(n-1),b(n-1)
      complex_dp, parameter   :: img = (0d0,1d0)
      real_dp                 :: ar(n-1),br(n-1)
      real_dp                 :: ai(n-1),bi(n-1)
      call cubicspline(ar,br,real(y),x,n,real(y1),m1,real(yn),mn)
      call cubicspline(ai,bi,imag(y),x,n,imag(y1),m1,imag(yn),mn)
      a = ar + img*ai
      b = br + img*bi
      end
      
c     ------------

      function chr(i)
      implicit none
      character(5)        :: chr
      integer, intent(in) :: i
      write(chr,'(I5)') i
      chr = adjustl(chr)
      end
      
