c--------------------------------------------------------------------------------
c Copyright (c) 2021 Peter Grünberg Institut, Forschungszentrum Jülich, Germany
c This file is part of SPEX and available as free software under the conditions
c of the MIT license as expressed in the LICENSE file in more detail.
c--------------------------------------------------------------------------------

c     Included from key.F .
c
# if   ARGTYPE >= 5
      subroutine getkey ARGTYPE (iunit,key,arg,section,default,status,writeout,require)
# elif ARGTYPE >= 3
      subroutine getkey ARGTYPE (iunit,key,arg,section,default,status,writeout,mini,maxi,mine,maxe,allow_eV,require)
# else
      subroutine getkey ARGTYPE (iunit,key,arg,section,default,status,writeout,mini,maxi,mine,maxe,require)
# endif

      use, intrinsic :: iso_fortran_env
      implicit none
      integer,      intent(in)              :: iunit
      character(*), intent(in)              :: key
      character(*), intent(in), optional    :: section
      logical,      intent(in), optional    :: writeout,require
      integer,      intent(out),optional    :: status
# if   ARGTYPE == 1
# define CTYPE i
      integer,      intent(inout)           :: arg
      integer,      intent(in), optional    :: default
      integer,      intent(in), optional    :: mini,mine
      integer,      intent(in), optional    :: maxi,maxe
# elif ARGTYPE == 2
# define CTYPE i
      integer,                  allocatable :: arg(:)
      integer,      intent(in), optional    :: default(:)
      integer,      intent(in), optional    :: mini,mine
      integer,      intent(in), optional    :: maxi,maxe
      integer                               :: i
# elif ARGTYPE == 3
# define CTYPE r
      real_dp,      intent(inout)           :: arg
      real_dp,      intent(in), optional    :: default
      real_dp,      intent(in), optional    :: mini,mine
      real_dp,      intent(in), optional    :: maxi,maxe
      logical,      intent(in), optional    :: allow_eV
# elif ARGTYPE == 4
# define CTYPE r
      real_dp,                  allocatable :: arg(:)
      real_dp,      intent(in), optional    :: default(:)
      real_dp,      intent(in), optional    :: mini,mine
      real_dp,      intent(in), optional    :: maxi,maxe
      logical,      intent(in), optional    :: allow_eV
      integer                               :: i
# elif ARGTYPE == 5
      logical,      intent(inout)           :: arg
      logical,      intent(in), optional    :: default
# elif ARGTYPE == 6
      character(*), intent(inout)           :: arg
      character(*), intent(in), optional    :: default
# elif ARGTYPE == 7
      character(:),             allocatable :: arg(:)
      character(*), intent(in), optional    :: default(:)
      integer                               :: i
# endif
      integer                               :: ios,ind,ind1,narg
      integer,                  parameter   :: maxlen = 65536
      character(40)                         :: keystr
      character(maxlen)                     :: line
      logical                               :: searchkey,writeout1
# if ARGTYPE != 5
      logical                               :: eV = .false.
# endif
      
# if ARGTYPE == 3 || ARGTYPE == 4
      if(present(allow_eV)) eV = allow_eV
# endif

      writeout1 = .false.
      if(present(writeout)) writeout1 = writeout

      if(.not.present(section)) then
        if(write_section/=' ') then
          write(6,'(A)') 'END'
        endif
      endif

      if(present(section)) then
        searchkey = .false. ; keystr = trim(key)//' (section '//trim(section)//')'
      else
        searchkey = .true.  ; keystr = trim(key)
      endif

# if ARGTYPE == 2 || ARGTYPE == 4 || ARGTYPE == 7
      narg = 0
      if(allocated(arg)) then
        narg = size(arg)
        if(narg==0) deallocate ( arg )
      endif
      if(present(default).and.allocated(arg)) then
        if(size(default)/=narg) Error('Number of arguments and default values do not agree.')
      endif
# endif

# if ARGTYPE == 5
      if(present(default)) then
        arg = default
      else
        arg = .false.
      endif
# endif

      rewind(iunit)

      do
        call getline(iunit,line,ind1,ios,.false.) ! ind1 unused
        if(ios/=0.or.line=='EXIT') exit
        line = adjustl(line) ! Remove leading spaces
        if(line(:8)=='SECTION '.and.key(:min(8,len(key)))/='SECTION ') then
          if(searchkey) searchkey = .false.
          if(present(section)) then
            if(adjustl(line(8:))==section) searchkey = .true.
            if(writeout1.and.write_section/=section) then
              if(write_section/=' ') write(6,'(A)') 'END'
              write(6,'(A)') 'SECTION '//trim(section)
              write_section = trim(section)//'      '
            endif
          endif
          cycle
        endif
        if(line=='END') then
          if(searchkey) searchkey = .false.
          if(.not.present(section)) searchkey = .true.
          cycle
        endif
        if(.not.searchkey) cycle
        ind = len_trim(key) + 1
        if(line(:ind)==key) then
          ind1 = index(line(ind:),'#') + ind - 1 ! truncate comments
          if(ind1==ind-1) ind1 = len(line)       !
# if ARGTYPE != 5
          if(present(status).and.line(ind:ind1)==' ') then
            status = 1
            if(writeout1) then
              if(present(section)) write(6,'(A'NoA) '  '
              write(6,'(A)') key
            endif
            goto 2
          endif
# endif          
# if ARGTYPE == 5
          arg = .not.arg
# else
          call getval ARGTYPE (line(ind:ind1),arg,narg,eV,keystr)
# endif
# if ARGTYPE < 5
          if(present(mini)) call check   ARGTYPE (arg,mini,-1,keystr)
          if(present(maxi)) call check   ARGTYPE (arg,maxi, 1,keystr)
          if(present(mine)) call check_e ARGTYPE (arg,mine,-1,keystr)
          if(present(maxe)) call check_e ARGTYPE (arg,maxe, 1,keystr)
# endif
          goto 3 ! write out and leave
        endif
      enddo

      if(present(status)) status = 0

 2    continue

      if(.not.present(default)) then
        if(present(status)) return ! arg remains undefined on return (status=0)
        if(.not.present(require)) then ; return ! only give error if keyword required
        else if(.not.require)     then ; return
        endif
# if ARGTYPE == 5
        Bug('No default given for logical key '//trim(keystr)//'.')
# endif
        Error('Keyword '//trim(keystr)//' not found.')
      endif

# if ARGTYPE == 2 || ARGTYPE == 4 || ARGTYPE == 7
      if(.not.allocated(arg)) allocate(arg,mold=default)
      if(size(arg)/=size(default)) Bug('Sizes of arg and default differ.')
# endif

      arg = default
      return

 3    if(present(status)) status = 2
      if(writeout1) then
# if   ARGTYPE == 5
        if(present(default)) then
          if(arg.neqv.default) then
            if(present(section)) write(6,'(A'NoA) '  '
            write(6,'(A)') key
          endif
        else
          if(arg) then
            if(present(section)) write(6,'(A'NoA) '  '
            write(6,'(A)') key
          endif
        endif
# else
        if(present(section)) write(6,'(A'NoA) '  '
        write(6,'(A'NoA) key//'       '(:7-len_trim(key))
# endif
# if   ARGTYPE == 1 || ARGTYPE == 3
        write(6,'(''  '''NoA) ; call write_ CTYPE (arg)
# elif ARGTYPE == 2 || ARGTYPE == 4
        do i = 1,narg
          write(6,'(''  '''NoA) ; call write_ CTYPE (arg(i))
        enddo
# elif ARGTYPE == 6
        write(6,'(3X,A'NoA) trim(arg)
# elif ARGTYPE == 7
        write(6,'('' '''NoA)
        do i = 1,narg
          write(6,'(2X,A'NoA) trim(arg(i))
        enddo
# endif
# if    ARGTYPE != 5
        write(6,*)
# endif
      endif

      end subroutine getkey ARGTYPE

c -------------------------
      
# if   ARGTYPE != 5
      subroutine getval ARGTYPE (line,arg,narg,eV,keystr)
      use util, only: chr,strlist
      use, intrinsic :: iso_fortran_env
      implicit none
      character(*), intent(in)               :: line,keystr
      logical,      intent(in)               :: eV
      integer,      allocatable              :: pos(:)
      integer                                :: narg,lenarg
      integer                                :: i,j,ios
      character(80)                          :: typ
      character(:), allocatable              :: word
# if   ARGTYPE == 1
      integer,      intent(out)              :: arg
      logical                                :: isinteger
      typ = ' integer'
# elif ARGTYPE == 2
      integer,                   allocatable :: arg(:)
      logical                                :: isinteger
      typ = ' integer'
# elif ARGTYPE == 3
      real_dp,      intent(out)              :: arg
      typ = ' real'
# elif ARGTYPE == 4
      real_dp,                   allocatable :: arg(:)
      typ = ' real'
# elif ARGTYPE == 6
      character(*), intent(out)              :: arg
      typ = ' character'
# elif ARGTYPE == 7
      character(:),              allocatable :: arg(:)      
      typ = ' character'
# endif
      
      call strlist(pos,line,' ')
      narg   = ubound(pos,1)
      lenarg = 1 ! avoid lenarg being undefined if narg==0
      if(narg>0) then
        lenarg = pos(1)
        do i = 2,narg
          lenarg = max(lenarg,pos(i)-pos(i-1)-1)
        enddo
        allocate(character(lenarg) :: word)
        word = adjustl(line(:pos(1)))
      endif

# if   ARGTYPE == 2 || ARGTYPE == 4 || ARGTYPE == 7
      
      if(allocated(arg)) then
        if(narg>size(arg)) then
          Error('Read '//Chr(narg)//' arguments after keyword '//trim(keystr)//' but needed only '//Chr(size(arg)))
        else if(narg<size(arg)) then
          Error('Read '//Chr(narg)//' arguments after keyword '//trim(keystr)//' but needed '//Chr(size(arg)))
        endif
      else
#   if ARGTYPE == 7
        allocate(character(lenarg) :: arg(narg))
#   else
        allocate(arg(narg))
#   endif
      endif

      do i = 1,narg
        if(i>1) word = adjustl(line(pos(i-1)+2:pos(i)))
#   if ARGTYPE == 2
        if(.not.isinteger(word)) Error('Argument '//trim(word)//' is not an integer.')
#   endif
#   if ARGTYPE < 5
        j = max(1, len_trim(word) - 1 )
        if(eV.and.word(j:)=='eV') then
          word(j:) = ' '
          read(word,*,iostat=ios) arg(i)
          arg(i) = arg(i) / hartree
        else
          read(word,*,iostat=ios) arg(i)
        endif
        if(ios/=0) Error('Could not read'//trim(typ)//' argument from '//trim(word))
#   else
        arg(i) = word
#   endif
      enddo
      
# else
      
      if(narg>1) Error('Read '//Chr(narg)//' arguments after keyword '//trim(keystr)//' but needed only one.')
      if(narg<1) Error('Could not read'//trim(typ)//' argument after keyword '//trim(keystr))
#   if ARGTYPE == 1
      if(.not.isinteger(word)) Error('Argument '//trim(word)//' after keyword '//trim(keystr)//' is not an integer.')
#   endif
#   if ARGTYPE < 5
      j = max(1, len_trim(word) - 1 )
      if(eV.and.word(j:)=='eV') then
        word(j:) = ' '
        read(word,*,iostat=ios) arg
        arg = arg / hartree
      else
        read(word,*,iostat=ios) arg
      endif
      if(ios/=0) Error('Could not read'//trim(typ)//' argument from '//trim(word)//' after keyword '//trim(keystr))
#   else
      arg = word
#   endif
# endif

      deallocate(pos)
      if(allocated(word)) deallocate(word)

      end subroutine getval ARGTYPE
# endif

# if   ARGTYPE < 5
      subroutine check ARGTYPE (arg,minmax,mode,keystr)
      use, intrinsic :: iso_fortran_env
      implicit none
      integer, intent(in) :: mode
      logical             :: ldum
      character(*)        :: keystr
      integer             :: i
# if   ARGTYPE == 1
      integer             :: arg,minmax,arg1
      ldum = arg*mode>minmax*mode      ; arg1 = arg
# elif ARGTYPE == 2
      integer             :: arg(:),minmax,arg1
      ldum = any(arg*mode>minmax*mode) ; arg1 = mode*maxval(mode*arg)
# elif ARGTYPE == 3
      real_dp             :: arg,minmax,arg1
      ldum = arg*mode>minmax*mode      ; arg1 = arg
# elif ARGTYPE == 4
      real_dp             :: arg(:),minmax,arg1
      ldum = any(arg*mode>minmax*mode) ; arg1 = mode*maxval(mode*arg)
# endif
      i = (mode + 1) / 2 + 1
      if(ldum) then
        write(0,'(A'NoA) 'Value out of range after keyword '//trim(keystr)
        write(0,'(A'NoA) ': ' ; call write_ CTYPE (arg1) ; write(0,'(1X,A'NoA) '<>'(i:i)//' ' ; call write_ CTYPE (minmax)
        write(0,'(A)')   ' .'
        Error('Value out of range after keyword '//trim(keystr))
      endif
      end subroutine check ARGTYPE
# endif

# if   ARGTYPE < 5
      subroutine check_e ARGTYPE (arg,minmax,mode,keystr)
      use, intrinsic :: iso_fortran_env
      implicit none
      integer, intent(in) :: mode
      logical             :: ldum
      character(*)        :: keystr
      integer             :: i
# if   ARGTYPE == 1
      integer             :: arg,minmax,arg1
      ldum = arg*mode>=minmax*mode      ; arg1 = arg
# elif ARGTYPE == 2
      integer             :: arg(:),minmax,arg1
      ldum = any(arg*mode>=minmax*mode) ; arg1 = mode*maxval(mode*arg)
# elif ARGTYPE == 3
      real_dp             :: arg,minmax,arg1
      ldum = arg*mode>=minmax*mode      ; arg1 = arg
# elif ARGTYPE == 4
      real_dp             :: arg(:),minmax,arg1
      ldum = any(arg*mode>=minmax*mode) ; arg1 = mode*maxval(mode*arg)
# endif
      i = mode + 2
      if(ldum) then
        write(0,'(A'NoA) 'Value out of range after keyword '//trim(keystr)
        write(0,'(A'NoA) ': ' ; call write_ CTYPE (arg1) ; write(0,'(1X,A'NoA) '<=>='(i:i+1)//' ' ; call write_ CTYPE (minmax)
        write(0,'(A)')   ' .'
        Error('Value out of range after keyword '//trim(keystr))
      endif
      end subroutine check_e ARGTYPE
# undef CTYPE
# endif
