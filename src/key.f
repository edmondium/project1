c--------------------------------------------------------------------------------
c Copyright (c) 2021 Peter Grünberg Institut, Forschungszentrum Jülich, Germany
c This file is part of SPEX and available as free software under the conditions
c of the MIT license as expressed in the LICENSE file in more detail.
c--------------------------------------------------------------------------------

c     The routine getkey reads values from a formatted input file.
c
c     The usage is
c     call getkey(<unitnumber>,<keyword>,<dest>[,section=<section>][,default=<default>][,status=<status>[,mini=<mini>][,maxi=<maxi>][,mine=<mine>][,maxe=<maxe>]
c
c     Looks for <keyword> in unit <unitnumber> and reads argument(s) <dest> behind it.
c     <dest> can be a single argument or an array of integer values, real_dp values or character strings.
c     NOTE: The array must be an allocatable array!
c
c     If <keyword> is not found:
c     <require> = true : the program stops with an error message
c     <default> given  : <dest>:=<default>.
c     <status>  given  : <status> takes one of the values
c                    0 : <keyword> not found,
c                    1 : <keyword> is found but no arguments are provided,
c                    2 : <keyword> and arguments are found. (If <default>=logical, status=0 if not found, status=2 if found.)
c                    In the cases 0 and 1, <dest>=<default> if <default> is present, otherwise <dest>=undefined on return.
c
c     <dest> can also be logical. Then, the occurrence of <keyword> switches its value
c     to .true., or if <default> is given to .not.<default>.
c
c     If <section> is given, <keyword> is searched for only between "SECTION <section>"
c     and "END".
c
c     With mini, maxi, mine and maxe a range for values can be specified:
c     <dest> >= <mini>, <dest> <= <maxi>, <dest> > <mine>, <dest> < <maxe>.
c
c     Everything in the line after "#" is treated as a comment.
c     Note that unknown keywords are not reported! Therefore, in order to prevent misspelled
c     keywords another routine should look through the input file and check the keywords.

# include "cppmacro.h"

      module key

      use global, only: hartree
      use util, only: chr
      use, intrinsic :: iso_fortran_env

      character(8), save :: write_section = ' '

      interface getkey
      module procedure  getkey1, getkey2, getkey3, getkey4, getkey5, getkey6, getkey7
      end interface

c     --------------------------

      contains
# define ARGTYPE 1
# include "getkey.inc"
# undef  ARGTYPE
# define ARGTYPE 2
# include "getkey.inc"
# undef  ARGTYPE
# define ARGTYPE 3
# include "getkey.inc"
# undef  ARGTYPE
# define ARGTYPE 4
# include "getkey.inc"
# undef  ARGTYPE
# define ARGTYPE 5
# include "getkey.inc"
# undef  ARGTYPE
# define ARGTYPE 6
# include "getkey.inc"
# undef  ARGTYPE
# define ARGTYPE 7
# include "getkey.inc"
# undef  ARGTYPE

c     ------------------

      subroutine write_i(arg)
      implicit none
      character(20)       :: form,line
      integer, intent(in) :: arg
      integer             :: i
      write(line,*) arg
      if(arg>=0) i = 1
      if(arg<0) i = 0
      write(form,'(A,I3,A)') '(I',len_trim(adjustl(line))+i,')'
      write(6,form NoA) arg
      end subroutine write_i

      subroutine write_r(arg)
      implicit none
      character(20)       :: form
      character(80)       :: line
      real_dp, intent(in) :: arg
      real_dp             :: lg
c      write(*,'(A,F10.5'NoA) '!!!',arg ; return
      if(nint(arg)==arg) then ; call write_i(nint(arg)) ; return ; endif
      if(arg==0) then ; lg = 0
      else            ; lg = log(abs(arg))/log(10d0) ; endif
      if(lg>=0) then
        write(form,'(A,I3.3,A,I3.3,A)') '(F',int(lg)+3+15,'.',15,')'
      else if(lg>=-3) then
        write(form,'(A,I3.3,A,I3.3,A)') '(F',int(-lg)+3+15,'.',int(-lg)+15,')'
      else
        form = '(ES9.2)'
      endif
      line = ' '
      write(line,form) arg
      if(lg>=-3) then
        do while (line(len_trim(line):len_trim(line))=='0')
          line = line(:len_trim(line)-1)
        enddo
      endif
      write(6,'(A'NoA) line(:len_trim(line))
      end subroutine write_r

c     ------------------

c     Checks the syntax of the GW input file ('spex.inp'). Does not check whether enough arguments are provided. (This is done in getkey.)
c     Sections within sections are not possible.
c     A leading '#' denotes a comment.

      subroutine checkkeys(iunit)
      use util, only: strlist
      use, intrinsic :: iso_fortran_env
      implicit none
      integer,       intent(in) :: iunit
      character(*),  parameter  :: key = 
     &                             'BZ NBAND CHKOLAP CHKMISM WRTKPT KPT GAUSS JOB WRITE MEM RESTART STOREIBZ ALIGNBD '//
     &                             'CORES BANDINFO DIPOLE ITERATE BANDOMIT KPTPATH MTACCUR CORESOC IBC PLUSSOC ENERGY '//
     &                             'BLOECHL CUTZERO MPIKPT MPISPLIT DELTAEX NOSYM TRSOFF STOREBZ TIMING PROJECT FIXPHASE '//
     &                             'SECTION MBASIS TOL SELECT LCUT GCUT CHKPROD WFADJUST OPTIMIZE NOAPW ADDBAS '//
     &                             'SECTION WFPROD MINCPW APPROXPW LCUT FFT MPIMT MPIPW '//
     &                             'SECTION COULOMB LEXP TSTCOUL NOSTORE STEPRAD CHKCOUL MULTIPOLE '//
     &                             'SECTION SUSCEP HILBERT TETRAF MULTDIFF DISORDER HUBBARD WGHTTHR PLASMA FSPEC FPADE '//
     &                             'SECTION SENERGY MESH FREQINT ZERO CONTINUE CONTOUR ALIGNVXC VXC MPIBLK MPISYM SPECTRAL '//
     &                                     'ORDER SMOOTH LOGDIV '//
     &                             'SECTION WANNIER UREAD WSCALE RSITE ORBITALS MAXIMIZE FROZEN DISENTGL BACKFOLD IRREP '//
     &                                     'CUTGOLD INTERPOL SUBSET WBLOCH PLOT '//
     &                             'SECTION LAPW GCUT LCUT EPAR LO '
      character(65536)          :: line
      character(256)            :: keyword,section
      integer, allocatable      :: pos(:)
      integer                   :: ios,ind,i,keyfirst,nkey
      logical, allocatable      :: exist(:)
      logical                   :: insection,found

      call strlist(pos,key,' ')

      nkey = ubound(pos,1)

      allocate(exist(nkey))

      exist     = .false.
      insection = .false.
      keyfirst  = 1
      rewind(iunit)

      do

        call getline(iunit,line,i,ios,.true.) ! i unused
        if(ios/=0.or.line=='EXIT') exit
        if(line==' ') cycle

        line    = adjustl(line) ! Remove leading spaces
        ind     = index(line,' ') ! Get keyword
        keyword = line(:ind-1)    !
        line    = adjustl(line(ind+1:))

        if(keyword=='END') then
          if(.not.insection) Error('Missing SECTION statement before END statement.')
          insection = .false.
          keyfirst  = 1
          write(6,'(A)') 'END'
          cycle
        endif
        if(keyword=='SECTION'.and.insection) Error('Missing END statement after SECTION statement.')

        ! Look for keyword in key-list (between two "SECTION" keys if insection=.true.)
# define KEY(i) key(pos(i-1)+2:pos(i))
        found = .false.
        i     = keyfirst
        do while(i<=nkey)
          if(KEY(i)==keyword) then
            if(keyword=='SECTION') then
              i         = i + 1
              section   = line
              if(section==' ') Error('Missing section keyword after SECTION.')
              if(section==KEY(i)) then
                insection = .true.
                keyfirst  = i+1
                found     = .true.
                write(6,'(A)') 'SECTION '//trim(KEY(i))
                exit
              endif
            else
              if(exist(i)) Error('Keyword '//trim(keyword)//' given more than once.')
              found    = .true.
              exist(i) = .true.
              if(insection) write(6,'(''  '''NoA)
              write(6,'(A)') keyword(:max(8,len_trim(keyword)))//' '//trim(line)
              exit
            endif
          endif
          if(KEY(i)=='SECTION') exit
          i = i + 1
        enddo

        if(.not.found) then
          if(keyword=='SECTION') then
            Error('Section keyword '//trim(section)//' unknown.')
          else
            if(insection) then
              Error('Keyword '//trim(keyword)//' unknown (in section '//trim(section)//').')
            else
              Error('Keyword '//trim(keyword)//' unknown (outside sections).')
            endif
          endif
        endif

      enddo

      deallocate(pos,exist)

      if(insection) Error('Missing END statement after SECTION statement.')
      write(6,*)

      end subroutine checkkeys

c     ------------------

c     Reads line from iunit. The line length is limited only by len(line). Comments #... are skipped.
c     Special comments ##... and ###... are written as comment to stdout and as a warning to stderr, respectively, if lcom is true.
c     Line continuation with a backslash (...\).
c     Input  : iunit - unit number
c              lcom  - if true, comments starting with "##" and "###" are written to stdout and stderr (see above).
c     Output : line  - line read from unit
c              llen  - length of line
c              ios   - I/O status (ios/=0 if eof has been reached)
      subroutine getline(iunit,line,llen,ios,lcom)
      implicit none
      integer,      intent(in)  :: iunit
      logical,      intent(in)  :: lcom
      character(*), intent(out) :: line
      integer,      intent(out) :: ios,llen
      integer,      parameter   :: leninp = 80 ! only used internally; does not limit readable line length
      character(leninp)         :: line1
      character(240)            :: msg
      integer                   :: ind,siz
      line = ' '
      llen = 0
      do
        read(iunit,'(A)',iostat=ios,iomsg=msg,advance='no',size=siz) line1
        line = line(:llen)//line1(:siz)
        llen = llen + siz ; if(llen>len(line)) Error('Line exceeded length of line string: '//chr(len(line)))
        if(ios/=0) then
          if(is_iostat_eor(ios)) then
            ios = 0
            ind = index(line,'#') ! strip comment
            if(ind/=0.and.ind<len(line)-1) then
              if(lcom.and.line(ind+1:ind+1)=='#') then
                if(line(ind+2:ind+2)=='#') then ; Warn('User-defined warning: '//trim(adjustl(line(ind+3:))))
                else                            ; write(6,'(A)') trim(adjustl(line(ind+1:)))
                endif
              endif                  
              line(ind:) = ' '
              llen       = ind - 1
            endif 
            ind = len_trim(line) ! line continuation
            if(ind>0) then
              if(line(ind:)=='\') then
                line(ind:) = ' '
                llen       = ind - 1
                cycle
              endif
            endif
            exit
          else if(is_iostat_end(ios)) then ; exit
          else                             ; Error('Unexpected error occurred while reading from input file: '//trim(msg))
          endif
        endif
      enddo
      end subroutine getline

c     ------------------

      end module key

