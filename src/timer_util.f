c--------------------------------------------------------------------------------
c Copyright (c) 2021 Peter Grünberg Institut, Forschungszentrum Jülich, Germany
c This file is part of SPEX and available as free software under the conditions
c of the MIT license as expressed in the LICENSE file in more detail.
c--------------------------------------------------------------------------------

# include "cppmacro.h"

      module timer_util

c     Timing routines

      Mpi( use Mwrapper )
      use, intrinsic :: iso_fortran_env
      integer, parameter :: maxtimer = 100, lenname = 30
      type timer_type
        character(len=lenname) :: name
        real                   :: cputime,cputime0,walltime
        integer_dp             :: walltime0
        integer                :: ncall
        logical                :: run
      end type timer_type
      type(timer_type)   :: timer(maxtimer) = [ ( timer_type(' ',0.,0.,0.,0,0_int64,.false.) , ntimer=1,maxtimer ) ]
      integer            :: ntimer
      
      contains

c
c     Returns index of timer with the given name.
c     Optional:
c       err  = false : program continues if timer name does not exist, see option "next" (default)
c              true  : program stops with an error if timer name does not exist
c       next = false : if timer name does not exist, the index is set to zero
c              true  : - " - , the number of timers (ntimer) is increased by one and the index is set to ntimer
      function timer_index(name,next,err)
      implicit none
      integer                       :: timer_index
      character(*),      intent(in) :: name
      logical, optional, intent(in) :: next,err
      if(len(name)>lenname) Bug('Timer name too long: '//trim(name))
      do timer_index = 1,ntimer
        if(timer(timer_index)%name==name) return
      enddo
      if(present(err)) then
        if(err) Bug('Timer name unknown: '//trim(name))
      endif
      if(present(next)) then
        if(next) then
          ntimer      = ntimer + 1
          timer_index = ntimer
          if(timer_index>maxtimer) then
            write(6,'(A)') 'List of timer names:',timer%name
            Bug('Timer limit exceeded. List written to stdout.')
          endif
          return
        endif
      endif
      timer_index = 0
      end function timer_index

c
c     Returns walltime (kind='wall') or cputime (kind else) of timer with the given name.
c     Optional err = true : program stops with an error if timer name does not exist (otherwise timer_get=0 is returned).
      function timer_get(name,kind,err)
      implicit none
      real                               :: timer_get
      character(*),           intent(in) :: name
      character(*), optional, intent(in) :: kind      
      logical,      optional, intent(in) :: err
      integer                            :: i
      i = timer_index(name)
      if(i==0) then
        if(present(err)) then
          if(err) Bug('Timer name does not exist: '//trim(name))
        endif
        timer_get = 0
      else
        if(present(kind)) then
          if(kind=='wall') then
            timer_get = timer(i)%walltime
            return
          endif
        endif
        timer_get = timer(i)%cputime
      endif
      end function timer_get
      
c
c     Starts a new timer with the given name(string) if the timer name does not exist yet.
c     Otherwise the previously stopped timer is restarted.
      subroutine timer_start(name)
      implicit none
      character(len=*), intent(in) :: name
      integer_dp                   :: wallrate,wallmax
      integer                      :: i
      i = timer_index(name,next=.true.)
      if(timer(i)%run) Bug('Timer still running: '//trim(name))
      timer(i)%name = name
      timer(i)%run  = .true.
      call cpu_time(timer(i)%cputime0)
      call system_clock(timer(i)%walltime0,wallrate,wallmax)
      end subroutine timer_start
    
c
c     Stops the timer with the given name (string). 
c     Elapsed time (cputime and walltime) in the global timer array.
c     Optional:
c       time  : walltime (kind='wall') or cputime (kind else) since timer_start (last time interval)
c       times : the same but accumulates the times
      subroutine timer_stop(name,time,times,kind)
      implicit none
      character(*),           intent(in)    :: name
      character(*), optional, intent(in)    :: kind
      real,         optional, intent(out)   :: time
      real,         optional, intent(inout) :: times      
      real                                  :: cputime,time1
      integer_dp                            :: walltime,wallrate,wallmax
      integer                               :: i      
      call cpu_time(cputime)
      call system_clock(walltime,wallrate,wallmax)
      i = timer_index(name,err=.true.)
      if(.not.timer(i)%run) Bug('Timer not running: '//trim(name))
      walltime = walltime - timer(i)%walltime0
      if(walltime<0) walltime = walltime + wallmax
      timer(i)%cputime  = timer(i)%cputime  + cputime - timer(i)%cputime0
      timer(i)%walltime = timer(i)%walltime + real(walltime)/wallrate
      timer(i)%ncall    = timer(i)%ncall    + 1
      timer(i)%run      = .false.
      if(present(time).or.present(times)) then
        if(kind=='wall') then ; time1 = real(walltime)/wallrate
        else                  ; time1 = cputime - timer(i)%cputime0
        endif
        if(present(time))  time  = time1
        if(present(times)) times = times + time1
      endif
      end subroutine timer_stop

c
c     Prints details of timers with the names stored in names(:) 
c     Optional:
c       title          : title to give to timing information
c       keep   = false : remove timers from list (default)
c                true  : keep timers
c       err    = 0     : skip printing timing if timer name does not exist (default)
c                1     : stop with error if timer name does not exist
c                2     : print zero timing if timer name does not exist
c       since  = true  : all timers since names(1) (size(names) must be 1)
c       select = true  : all timers that contain names(1) in their names (size(names) must be 1),
c                        the key names(1) is removed from the names
c       prnt   = true  : print information (default true)
      subroutine timer_print(names, title, keep, err, since, select, prnt, unit)
      use util, only: chr
      use file
      implicit none
      character(*),           intent(in) :: names(:)
      character(*), optional, intent(in) :: title
      integer,      optional, intent(in) :: err,unit
      logical,      optional, intent(in) :: keep,since,select,prnt
      logical                            :: since1,select1,prnt1
      character(lenname)                 :: name
      real                               :: cputime,walltime
      integer                            :: delete(size(names)),iunit
      integer                            :: ncall,i,j,k,l,j1,j2
      logical                            :: ldum
      delete  = 0
      since1  = .false.
      select1 = .false.      
      prnt1   = .true.
      iunit   = 6
      if(present(since))  since1  = since
      if(present(select)) select1 = select
      if(present(prnt))   prnt1   = prnt
      if(present(unit).and.prnt1)  iunit = fopen('spex.time.'//Chr(unit),status='unknown')
      if(present(title).and.prnt1) write(iunit,'(/'//chr(len_trim(title))//'A/A)') ('-',i=1,len_trim(title)),trim(title)
      l = 0
      if(since1) then
        if(size(names)>1) Bug('Expected only one timer name with since=.true.')
        j1 = timer_index(names(1),err=.true.)
        j2 = ntimer
      else if(select1) then
        if(size(names)>1) Bug('Expected only one timer name with since=.true.')
        j1 = 1
        j2 = ntimer
      else
        j1 = 1
        j2 = size(names)
      endif
      do j = j1,j2
        if(since1.or.select1) then
          i    = j
          name = timer(i)%name
          if(select1) then
            if(index(name,names(1))==0) cycle
            name = remove_string(name,names(1))
          endif
        else
          name = names(j)
          i    = timer_index(name,next=.false.)
        endif
        cputime = -1
        if(i==0) then
          if(present(err)) then
            if(err==1) then
              Bug('Timer name does not exist: '//trim(name))
            else if(err==2) then
              cputime  = 0
              walltime = 0
              ncall    = 0
            endif
          endif
        else
          if(timer(i)%run) Bug('Timer still running: '//trim(name))
          cputime  = timer(i)%cputime
          walltime = timer(i)%walltime
          ncall    = timer(i)%ncall
        endif
        if(cputime>=0) then
          l = l + 1
          if(prnt1) write(iunit,'(A)') name//'CPU: '//trim(chr(cputime,'F7.2'))//          
     &                                '      Wall: '//trim(chr(walltime,'F7.2'))//'      #calls: '//trim(chr(ncall,'I3'))
        endif
      enddo
      if(prnt1.and.l==0) Warn('Did not find a valid timer for timing information.')
      if(present(title).and.prnt1) write(iunit,'('//chr(len_trim(title))//'A/)') ('-',i=1,len_trim(title))
      if(present(unit).and.prnt1) close(iunit)
      if(present(keep)) then
        if(keep) return
      endif
      if(since1) then
        i = timer_index(names(1),err=.true.)
        do j = i,ntimer
          timer(j) = timer_type(' ',0.,0.,0.,0,0_int64,.false.)
        enddo
        ntimer = i - 1
      else if(select1) then
        j = 0
        do i = 1,ntimer
          name = timer(i)%name
          if(index(name,names(1))==0) then
            j        = j + 1
            timer(j) = timer(i)
          endif
        enddo
        do i = j+1,ntimer
          timer(i) = timer_type(' ',0.,0.,0.,0,0_int64,.false.)
        enddo
        ntimer = j
      else
        do j = 1,size(names)
          i = timer_index(names(i),next=.false.)
          if(i/=0) then
            do k = i,ntimer-1
              timer(k) = timer(k+1)
            enddo
            timer(ntimer) = timer_type(' ',0.,0.,0.,0,0_int64,.false.)            
          endif
        enddo
        ntimer = ntimer - size(names)
      endif

      contains

      function remove_string(string,key)
      implicit none
      character(*), intent(in) :: string
      character(*), intent(in) :: key
      character(len(string))   :: remove_string
      integer                  :: ind,lenkey
      lenkey = len(key)
      ind    = index(string,key)
      if(ind/=0) then
        remove_string = name(:ind-1)//name(ind+lenkey:)
      else
        remove_string = name
      endif
      end function remove_string
      
      end subroutine timer_print

c
c     Delete timer information
c     Optional:
c       names        : list of timer names to be deleted (otherwise all timers are deleted)
c       since = true : delete all timer names since names(1) (size(names) must be 1)
      subroutine timer_del(names,since)
      implicit none
      character(*), optional, intent(in) :: names(:)
      logical,      optional, intent(in) :: since
      integer                            :: i,j
      if(present(names)) then
        if(present(since)) then
          if(since) then
            if(size(names)>1) Bug('Expected only one timer name with since=.true.')
            i = timer_index(names(1),err=.true.)
            do j = i,ntimer
              timer(j) = timer_type(' ',0.,0.,0.,0,0_int64,.false.)
            enddo
            ntimer = i - 1
            return
          endif
        endif
        do j = 1,size(names)
          i        = timer_index(names(j),err=.true.)          
          timer(i) = timer_type(' ',0.,0.,0.,0,0_int64,.false.)
        enddo
        ntimer = ntimer - size(names)
      else
        do j = 1,ntimer
          timer(j) = timer_type(' ',0.,0.,0.,0,0_int64,.false.)
        enddo
        ntimer = 0
      endif
      end subroutine timer_del

c
c     Replaces cpu_done. Calls timer_stop.
      subroutine timer_done(name)
      use util, only: chr
      implicit none
      character(*), intent(in) :: name
      real                     :: cputime
      call timer_stop(name,time=cputime)
      write(6,'(A)') 'done   ( Timing: '//Chf(cputime,'F16.2')//' )'
      end subroutine timer_done

      end module timer_util
