c--------------------------------------------------------------------------------
c Copyright (c) 2021 Peter Grünberg Institut, Forschungszentrum Jülich, Germany
c This file is part of SPEX and available as free software under the conditions
c of the MIT license as expressed in the LICENSE file in more detail.
c--------------------------------------------------------------------------------

c     Routines for reading and writing dump files (files with intermediate results during k loop).
c

# include "cppmacro.h"
# include "jobtype.h"

c     ------------

c
c     Read dump files fname (and fname.* files from MPIKPT calculations)
c

c     lkpt = k-point map
c     typ >= 0 : read dump file of type typ and merge data
c     typ  = 0 : define frequency mesh freq (from typ=2 file), oselfc, and job1%full. (If block and freq are allocated, it is checked if they are identical to the definition in the file.)
c     typ  = 1 : read exchange matrix      (-> selfx)
c     typ  = 2 : read correlation matrix   (-> selfc)
c     typ  = 3 : read Coulomb matrices     (-> barew,screenw)
c     typ  = 4 : read contracted W matrix  (-> screenk_mt,screenk_pw)
c     typ  < 0 : remove dump files (not used)
c
# define FN ' File: '//trim(fname)
      subroutine read_dump(typ,lkpt,dkpt,fname,job1)
      use global
      use util, only: chr
      use arrays
      use file
      implicit none
      type(jobtype), intent(in)    :: job1
      integer,       intent(in)    :: typ,dkpt
      logical,       intent(inout) :: lkpt(dkpt)
      character(*),  intent(in)    :: fname
      logical                      :: ldum
      integer                      :: i,iunit,ios,nfile
      character(6)                 :: sfx
      inquire(file=trim(fname),exist=ldum)
      if(.not.ldum) return
      iunit = fopen(trim(fname),form='unformatted',status='old')
      read(iunit,iostat=ios) i,nfile ; if(ios/=0.or.nfile==0) nfile = 1
      call fclose(iunit)
      if(nfile<0) Bug('Number of sig files < 0.')
      if(typ>=0) then
        ! Read dump files and sum
        if(typ>0.and.nfile>1) write(6,'(/A)') 'Try to merge '//trim(fname)//'(.*) files ...'
        do i = 1,nfile
          sfx = ' ' ; if(i>1) sfx = '.'//trim(chr(i))
          inquire(file=trim(fname)//trim(sfx),exist=ldum)
          if(ldum) call read_dump_data(typ,lkpt,dkpt,trim(fname)//trim(sfx),job1,i>1)
        enddo
      else if(typ<0) then
        ! Remove old files
        write(6,'(/A)') 'Old '//trim(fname)//'(.*) files are removed.'
        do i = 1,nfile
          sfx = ' ' ; if(i>1) sfx = '.'//trim(chr(i))
          inquire(file=trim(fname)//trim(sfx),exist=ldum)
          if(ldum) then
            iunit = fopen(trim(fname)//trim(sfx))
            call fclose(iunit,status='delete')
          endif
        enddo
      endif
      end

c     ------------

c
c     Reads data from dump file fname
c
c     Current dump file version: 4  (versions 2 and 3 can be read)
c
c     ladd = .true. : Arrays are added (otherwise overwritten)
c
      subroutine read_dump_data(typ,lkpt,dkpt,fname,job1,ladd)
      use global
      use util, only: chr
      use arrays
      use file
      use, intrinsic :: iso_fortran_env
      implicit none
      type(jobtype), intent(in)    :: job1
      integer,       intent(in)    :: typ,dkpt
      logical,       intent(in)    :: ladd
      logical,       intent(inout) :: lkpt(dkpt)
      character(*),  intent(in)    :: fname
      logical                      :: ldum
      logical,       allocatable   :: undef(:)
      integer,       allocatable   :: iarr(:)
      real_dp,       allocatable   :: frq(:)
      complex_dp,    allocatable   :: carr2(:,:),carr3(:,:,:)
      MCOMPLEX_dp,   allocatable   :: marr(:),marr3(:,:,:)
      real_dp                      :: kvec(3)
      logical                      :: lhlp(dkpt)
      integer                      :: iunit,ios,iblock,version InvC(ifreq)
      integer                      :: i,j,k,ii,n,m,s,nb,ib
      integer                      :: nfrq
      character(160)               :: error
      complex_dp,    allocatable   :: carr4(:,:,:,:)
      complex_dp                   :: cfrq
      logical,       allocatable   :: undef2(:,:)
      integer                      :: jj,nsit,rsit(3)
      inquire(file=trim(fname),exist=ldum)
      if(.not.ldum) return
      if(typ/=0) write(6,'(A)') 'File "'//trim(fname)//'" exists and will be read ...'
      iunit = fopen(fname,form='unformatted',status='old')
      ! Read preamble and check parameters
      call param_dump(typ,version,error,iunit,job1)
      if(error/=' ') Error(trim(error)//FN)
      if(typ==0) then
        call fclose(iunit)
        return
      endif
      ! Read k-point map
      lhlp = .false.
      do
        read(iunit,iostat=ios) i ; if(ios/=0) Error('Read error (k map).'//FN)
        if(i==0) exit
        if(i<0)     Error('k-point index negative. (bug?)'//FN)
        if(i>dkpt)  Error('k-point index exceeds maximal index. Wrong BZ parameters?'//FN)
        if(lkpt(i)) Error('k point '//trim(chr(i))//' already exists. Wrong merge?'//FN)
        lkpt(i) = .true.
        lhlp(i) = .true. ; if(version<3.and.typ==1) then ; lkpt(dkpt) = .true. ; lhlp(dkpt) = .true. ; endif
      enddo
      if(any(lhlp)) then
        write(6,'(A'NoA) trim(fname)//':'
        do i = 1,nkpti
          if(lhlp(i)) then
            j = 0
            if(i>1)       then ; if(lhlp(i-1)) j = 1     ; endif
            if(i<nkpti)   then ; if(lhlp(i+1)) j = j + 2 ; endif
            if     (j==1) then ; write(6,'(A'NoA) '-'//trim(chr(i))
            else if(j/=3) then ; write(6,'(A'NoA) ' '//trim(chr(i))
            endif
          endif
        enddo
        if(dkpt>nkpti)   then ; if(lhlp(nkpti+1)) write(6,'('' C'''NoA) ; endif
        if(dkpt>nkpti+1) then ; if(lhlp(nkpti+2)) write(6,'('' c'''NoA) ; endif
        write(6,*)
      endif
      ! Read data
      if(typ<=2) then ! self-energy
        if(version>=3) then ; read(iunit,iostat=ios) nb ; if(ios/=0) Error('Read error (nb).'//FN)
        else                ; nb = nblock
        endif
        if(nb<=0) Error('Number of blocks not positive. (bug?)'//FN)
        allocate(undef(nblock),iarr(size(block,1)))
        undef = .true.
        do ib = 1,nb
          iarr = 0
          read(iunit,iostat=ios) n,k,s,iarr(:min(n,size(iarr))),kvec
          if(ios/=0) then
            backspace(iunit)
            read(iunit,iostat=ios) n,k,s,iarr(:min(n,size(iarr)))
            kvec = -1 ! flag for not defined
          endif
          if(ios/=0) Error('Error reading parameters for block '//trim(chr(ib))//FN)
          if(version>=3) then
            iblock = 0
            do i = 1,nblock
              m = sizeblock(i)
              if(all([n,k,s,iarr(:min(n,m))]==(/m,job1%kpt(block(1,i)),job1%spin(block(1,i)),
     &            job1%band(block(:min(n,m),i))/))) then
                if(all(kvec/=-1).and.any(abs(kvec-kpt(:,k))>1d-12))
     &            Error('JOB k point ('//Chn(kpt(:,k),',')//') differs from stored one ('//Chn(kvec,',')//').'//FN)
                iblock        = i
                undef(iblock) = .false.                
                exit
              endif
            enddo
          else
            if(n/=sizeblock(ib).or.k/=job1%kpt(block(1,ib)).or.s/=job1%spin(block(1,ib))
     &        .or.any(iarr/=job1%band(block(:,ib)))) Error('Parameters changed for block '//trim(chr(ib))//FN)
            iblock = ib
          endif
          k = 0
          do i = 1,iblock-1
            m = sizeblock(i)
            if(typ==1) then ; k = k + m*(m+1)/2
            else            ; k = k + m**2
            endif
          enddo
          if(typ==1) then
            if(iblock==0) then
              read(iunit,iostat=ios)
            else
              if(ladd) allocate(marr(n*(n+1)/2))
              if(ladd) marr = selfx(k+1:k+n*(n+1)/2)
              read(iunit,iostat=ios) selfx(k+1:k+n*(n+1)/2)
              if(ladd) selfx(k+1:k+n*(n+1)/2) = selfx(k+1:k+n*(n+1)/2) + marr
              if(ladd) deallocate(marr)
            endif
          else
            if(iblock==0) then
              read(iunit,iostat=ios)
            else
              if(ladd) allocate(carr2(n**2,size(selfc,2)))
              if(ladd) carr2 = selfc(k+1:k+n**2,:)
# ifdef INV
              read(iunit,iostat=ios) (((selfc(k+i+(j-1)*n,ifreq),i=1,j),j=1,n),ifreq=1,size(selfc,2))
              do j = 1,n
                do i = 1,j-1
                  selfc(k+j+(i-1)*n,:) = selfc(k+i+(j-1)*n,:)
                enddo
              enddo
# else
              read(iunit,iostat=ios) selfc(k+1:k+n**2,:)
# endif
              if(ladd) selfc(k+1:k+n**2,:) = selfc(k+1:k+n**2,:) + carr2
              if(ladd) deallocate(carr2)
            endif
          endif
          if(ios/=0) then
            if(iblock==0) then ; Error('Error skipping self-energy block'//FN)
            else               ; Error('Error reading self-energy for block '//trim(chr(iblock))//FN)
            endif
          endif
        enddo
        if(any(undef).and.version>=3) then
          write(0,'(A'NoA) 'read_dump_data: Missing or incomplete blocks:'
          do i = 1,nblock
            if(undef(i)) write(0,'(I6'NoA) i
          enddo
          write(0,*)
          Error('File does not contain the required data.'//FN)
        endif
        deallocate(undef,iarr)
      else if(typ==3) then ! Coulomb matrix
        allocate(undef2(nfreqc,nsite),undef(nsite))
        if(ladd) allocate(carr4(nwan,nwan,nwan,nwan))
        undef2 = .true.
        undef  = .true.
        read(iunit,iostat=ios) nsit,nfrq ; if(ios/=0) Error('Could not read number of RSITE vectors and frequencies.'//FN)
        do ii = 1,nsit
          read(iunit,iostat=ios) rsit    ; if(ios/=0) Error('Could not read RSITE vector.'//FN)
          i = 0 ; do j = 1,nsite ; if(all(rsit==rsite(:,j)).and.undef(j)) i = j ; enddo
          if(i==0) then
            call skip(2*nfrq+1)
          else
            if(ladd) carr4 = barew(:,:,:,:,i)
            read(iunit,iostat=ios) barew(:,:,:,:,i) ; if(ios/=0) Error('Read error (barew).'//FN)
            if(ladd) barew(:,:,:,:,i) = barew(:,:,:,:,i) + carr4
            undef(i) = .false.
            do jj = 1,nfrq
              read(iunit,iostat=ios) cfrq ; if(ios/=0) Error('Could not read frequency.'//FN)
              j = 0 ; do k = 1,nfreqc ; if(abs(cfrq-freqc(k))<1d-10.and.undef2(k,i)) j = k ; enddo
              if(j==0) then
                call skip(1)
              else
                if(ladd) carr4 = screenw(:,:,:,:,j,i)
                read(iunit,iostat=ios) screenw(:,:,:,:,j,i) ; if(ios/=0) Error('Read error (screenw).'//FN)
                if(ladd) screenw(:,:,:,:,j,i) = screenw(:,:,:,:,j,i) + carr4
                undef2(j,i) = .false.
              endif
            enddo
          endif
        enddo
        if(ladd) deallocate(carr4)
        if(any(undef2)) then
          write(0,'(A)') 'Data not found in dump file (rsite,freq):'
          do i = 1,nsite
            if(any(undef2(:,i))) then
              write(0,'(3I3)')     rsite(:,i)
              write(0,'(2F15.10)') pack(freqc,undef2(:,i))
            endif
          enddo
          Error('Requested data not found in dump file.'//FN)
        endif
        deallocate(undef2,undef)
      else if(typ==4) then ! W contraction
        allocate(undef(nfreq),frq(1))
        undef = .true.
        if(ladd) allocate(carr3(maxlmindxm,maxlmindxm,ncent))
        read(iunit,iostat=ios) nfrq     ; if(ios/=0) Error('Could not read number of frequencies.'//FN)
        do ii = 1,nfrq
          read(iunit,iostat=ios) frq(1) ; if(ios/=0) Error('Could not read frequency.'//FN)
          i = 0 ; do j = 1,nfreq ; if(abs(frq(1)-freq(j))<1d-10.and.undef(j)) i = j ; enddo
          if(i==0) then
            call skip(1)
          else
            if(ladd) carr3 = screenk_mt(:,:,:,i)
            read(iunit,iostat=ios) screenk_mt(:,:,:,i) ; if(ios/=0) Error('Read error (screenk_mt).'//FN)
            if(ladd) screenk_mt(:,:,:,i) = screenk_mt(:,:,:,i) + carr3
            undef(i) = .false.
          endif
        enddo
        if(ladd) deallocate(carr3)
        if(any(undef)) then
          write(0,'(A'NoA) 'Data not found in dump file (freq):'
          write(0,'(F15.10)') pack(freq,undef)
          Error('Requested data not found in file.'//FN)
        endif
        if(allocated(screenk_pw)) then
          if(ladd) allocate(marr3(size(screenk_pw,1),size(screenk_pw,2),size(screenk_pw,3)))
          if(ladd) marr3 = screenk_pw
          read(iunit,iostat=ios) screenk_pw ; if(ios/=0) Error('Read error (screenk_pw).'//FN)
          if(ladd) screenk_pw = screenk_pw + marr3
          if(ladd) deallocate(marr3)
        endif
        deallocate(undef,frq)
      endif
      call fclose(iunit)

      contains

      subroutine skip(n)
      implicit none
      integer, intent(in) :: n
      integer             :: i
      do i = 1,n
        read(iunit,iostat=ios)
        if(ios/=0) Error('Error while skipping record '//trim(chr(i))//' of '//trim(chr(n))//'.'//FN)
      enddo
      end subroutine skip

      end

c     ------------

c
c     Writes  dump file for kpoints kpt1(:nkpt1)
c     Deletes dump file if nkpt1 = 0
c
c     See above for typ.
c
# undef FN
# define FN ' File: '//trim(fname)//trim(sfx)
      subroutine write_dump(typ,kpt1,nkpt1,fname,job1 MpiC(color) MpiC(maxcolor) )
      use global
      use arrays
      use file
      use, intrinsic :: iso_fortran_env
      Mpi2( use util, only: chr )
      implicit none
      type(jobtype), intent(in) :: job1
      integer,       intent(in) :: typ,nkpt1,kpt1(nkpt1) MpiC(color) MpiC(maxcolor)
      character(*),  intent(in) :: fname
      character(160)            :: error
      logical                   :: lexist
      integer,       parameter  :: version = 4
      integer                   :: version1,iunit
      integer                   :: iblock,i,j,k,n,ios InvC(ifreq)
      character(6)              :: sfx
      sfx = ' ' ; Mpi( if(color>1) sfx = '.'//trim(chr(color)) )
      inquire(file=trim(fname)//sfx,exist=lexist)
      if(nkpt1==0) then
        if(lexist) then
          iunit = fopen(trim(fname)//sfx,form='unformatted',status='unknown')
          call fclose(iunit,status='delete')
        endif
        return
      endif
      iunit = fopen(trim(fname)//sfx,form='unformatted',status='unknown')
      if(lexist) then
        ! Attempt to read if file exists
        call param_dump(typ,version1,error,iunit,job1)
        lexist = error==' '.and.version1==version
      endif
      if(.not.lexist) then
        rewind(iunit)
        ! Write basic parameters if file does not exist or read attempt failed
        write(iunit) version MpiC(maxcolor)
        if(typ<=2) then
          i = 0 ; if(job1%full) i = 1
          write(iunit) typ,nkpti,i,job1%label,checksum
          if(typ==2) then
            if(oselfc==1) then
              write(iunit) nfreq,oselfc
              write(iunit) freq(:nfreq)
            else
              write(iunit) nfreqr,oselfc
              write(iunit) freqr(:nfreqr)
            endif
          endif
        else if(typ==3) then
          write(iunit) typ,nkpti
          write(iunit) nwan,spinw
        else if(typ==4) then ; iblock = 0 ; if(allocated(screenk_pw)) iblock = size(screenk_pw)
          write(iunit) typ,nkpti
          write(iunit) size(screenk_mt)/size(screenk_mt,4),iblock
        endif
      else
        ! Skip k points already present
        do
          read(iunit,iostat=ios) i ; if(ios/=0) Error('Read error.'//FN)
          if(i==0) exit
          if(i<0)          Error('k-point index negative. (bug?)'//FN)
          if(any(i==kpt1)) Error('Cannot overwrite existing k point.'//FN)
        enddo
        backspace(iunit)
      endif
      ! Write k points
      do k = 1,nkpt1
        write(iunit) kpt1(k)
      enddo
      write(iunit) 0
      if(typ<=2) then ! Write selfx(c)
        write(iunit) nblock,size(block,1),block
        k = 0
        do iblock = 1,nblock
          n = sizeblock(iblock)
          j = job1%kpt(block(1,iblock))
          write(iunit) n,j,job1%spin(block(1,iblock)),job1%band(block(:n,iblock)),kpt(:,j)
          if(typ==1) then
            write(iunit) selfx(k+1:k+n*(n+1)/2)
            k = k + n*(n+1)/2
          else
# ifdef INV
            write(iunit) (((selfc(k+i+(j-1)*n,ifreq),i=1,j),j=1,n),ifreq=1,size(selfc,2))
# else
            write(iunit) selfc(k+1:k+n**2,:)
# endif
            k = k + n**2
          endif
        enddo
      else if(typ==3) then ! Write Coulomb matrices
        write(iunit) nsite,nfreqc
        do i = 1,nsite
          write(iunit) rsite(:,i)
          write(iunit) barew(:,:,:,:,i)
          do j = 1,nfreqc
            write(iunit) freqc(j)
            write(iunit) screenw(:,:,:,:,j,i)
          enddo
        enddo
      else if(typ==4) then ! Write W contraction
        write(iunit) nfreq
        do j = 1,nfreq
          write(iunit) freq(j)
          write(iunit) screenk_mt(:,:,:,j)
        enddo
        if(allocated(screenk_pw)) write(iunit) screenk_pw
      endif
      call fclose(iunit)

      end

c     ------------

c
c     Read and test parameters in preamble of dump file (iunit,fname).
c     (see read_dump for the case typ=0.)
c
# define ERROR(arg) then; error = arg; return; endif
      subroutine param_dump(typ,version,error,iunit,job1)
      use global
      use arrays
      implicit none
      type(jobtype)               :: job1 ! job1%full modified if typ=0
      integer,        intent(out) :: version
      character(160), intent(out) :: error
      integer,        intent(in)  :: typ,iunit
      integer                     :: ios,typ1,nkpti1,i,siz(2),nfrq,oselfc1,nwan1,spin1
      logical                     :: lfull
      real_dp                     :: checksum1
      real_dp,        allocatable :: frq(:)
      character(6)                :: jlab
      error = ' '
      jlab  = ' '
      lfull = .false.
      read(iunit,iostat=ios) version
      if(ios/=0) ERROR('Read error (version).')
      if(typ<=2) then
        if(all(version/=[2,3,4])) ERROR('Incompatible file format.')
        if(version>=3) then ;                  read(iunit,iostat=ios) typ1,nkpti1,i,jlab,checksum1
          if(ios/=0) then ; backspace(iunit) ; read(iunit,iostat=ios) typ1,nkpti1,i,jlab ; checksum1 = 0
          if(ios/=0) then ; backspace(iunit) ; read(iunit,iostat=ios) typ1,nkpti1,i      ; endif ; endif
        else                                 ; read(iunit,iostat=ios) nkpti1,i,siz,typ1  ; write(6,'(A)') 'Old dump file version 2'
        endif
        lfull = i/=0
        if(ios/=0) ERROR('Read error (parameters)')
        if(typ==0) job1%full = lfull
        if(typ/=0.and.(lfull.neqv.job1%full)) ERROR('FULL/diagonal mismatch.')
      else
        if(version/=4) ERROR('Incompatible file format.')
        read(iunit,iostat=ios) typ1,nkpti1
      endif
      if(lfull) then
        if(checksum1==0) then
          Warn('Did not find checksum. Cannot check consistency of wave functions.')
        else if(abs(checksum1-checksum)>acc_checksum*sqrt(dble(size(cmt)+size(cpw)))) then
          error = 'Checksum failure. Current wave functions inconsistent with stored self-energy matrix. Remove dump file.' ; return
        endif
      endif
      if(ios/=0)                 ERROR('Read error (param).')
      if(nkpti1/=nkpti)          ERROR('Number of irr. k points does not agree.')
      if(typ/=0.and.typ1/=typ)   ERROR('Data type does not agree.')
      if(typ==0.and.typ1/=2)     ERROR('Cannot read frequency mesh from exchange file.')
      if(typ<=2) then
        if(version==2) then
          if(typ==1) then
            if(siz(1)/=size(selfx))   ERROR('Self-energy matrix size does not agree.')
            if(siz(2)/=1)             ERROR('Frequencies in exchange self-energy?')
          else
            if(siz(1)/=size(selfc))   ERROR('Self-energy matrix size does not agree.')
            if(siz(2)/=size(selfc,2)) ERROR('Number of frequencies does not agree.')
          endif
        else if(typ==2.or.typ==0) then
          read(iunit,iostat=ios) nfrq,oselfc1
          if(ios/=0)                  ERROR('Read error (param2).')
          allocate(frq(nfrq))
          read(iunit,iostat=ios) frq
          if(ios/=0)                  ERROR('Read error (frq).')
          if(typ==0) then
            if(oselfc/=oselfc1) then
              if(oselfc>1) then
                Warn('Self-energy type (CONTOUR) replaced by the one from sigc/t file (CONTINUE).')
                nfreqr = 0 ; deallocate(freqr)
                nfreqc = 0 ; if(allocated(freqc)) deallocate(freqc)
              else
                Warn('Self-energy type (CONTINUE) replaced by the one from sigc/t file (CONTOUR).')
              endif
            endif
            oselfc = oselfc1
            if(oselfc==1) then      
              if(allocated(freq)) then
                if(nfreq/=nfrq.or.any(freq/=frq)) then
                  Warn('Imag. frequency mesh replaced by the one from sigc/t file.')
                  deallocate(freq)
                endif
              endif
              if(.not.allocated(freq)) then
                nfreq = nfrq ; allocate ( freq(nfreq) )
                freq  = frq  ; deallocate ( frq )
              endif
            else
              if(allocated(freqr)) then
                if(nfreqr/=nfrq.or.any(freqr/=frq)) then
                  Warn('Real frequency mesh replaced by the one from the sigc/t file.')
                  deallocate(freqr)
                endif
              endif
              if(.not.allocated(freqr)) then
                nfreqr = nfrq ; allocate ( freqr(nfreqr) )
                freqr  = frq  ; deallocate ( frq )
              endif
            endif
            return
          else
            if(oselfc1/=oselfc)                ERROR('Self-energy type does not agree.')
            if(oselfc==1) then
              if(nfrq/=nfreq)                  ERROR('Imag. frequency mesh has different size.')
              if(maxval(abs(frq-freq))>1d-10)  ERROR('Imag. frequency mesh has different values.')
            else if(allocated(freqr)) then
              if(nfrq/=nfreqr)                 ERROR('Real frequency mesh has different size.')
              if(maxval(abs(frq-freqr))>1d-10) ERROR('Real frequency mesh has different values.')
            endif
            deallocate(frq)
          endif
        else if(typ==1.and.jlab/=' ') then
          if(jlab/=job1%label) then
            if(all(jlab/=['HF    ','GW    ','RPAENE','HFENE ','GWT   ']).or.all(job1%type/=[J_HF,J_GW,J_RPA,J_HFE]))
     &                                         ERROR('Cannot use '//trim(jlab)//' data for current job.')
          endif
        endif
      else if(typ==3) then
        read(iunit,iostat=ios) nwan1,spin1
        if(ios/=0)       ERROR('Read error.')
        if(nwan1/=nwan)  ERROR('Number of Wannier functions does not agree.')
        if(spin1/=spinw) ERROR('Spin index does not agree.')
      else if(typ==4) then
        read(iunit,iostat=ios) siz
        if(ios/=0)                                      ERROR('Read error (param2).')
        if(siz(1)/=size(screenk_mt)/size(screenk_mt,4)) ERROR('MT W contraction size incorrect.')
        if(allocated(screenk_pw)) then
          if(siz(2)==0)                                 ERROR('PW W contraction array not present.')
          if(siz(2)/=size(screenk_pw))                  ERROR('PW W contraction size incorrect.')
        endif
      endif
      end

# undef ERROR
