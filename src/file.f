c--------------------------------------------------------------------------------
c Copyright (c) 2021 Peter Grünberg Institut, Forschungszentrum Jülich, Germany
c This file is part of SPEX and available as free software under the conditions
c of the MIT license as expressed in the LICENSE file in more detail.
c--------------------------------------------------------------------------------

c     The function fopen and the routine fclose take care that no unit number is duplicated.
c
c     Example:
c     iunit = fopen('filename',status='old',action='read',form='formatted')
c     ...
c     call fclose(iunit)      

# include "cppmacro.h"

      module file

      integer, parameter :: MINUNIT=10,MAXUNIT=99
      logical            :: unitopen(MINUNIT:MAXUNIT)=.false.

c     ------------

      contains
      function fopen(name,status,action,form,access,position,recl,numbered)

      use, intrinsic :: iso_fortran_env
      implicit none
      integer                   :: fopen
      character(len=*),optional :: name
      character(256)            :: name1
      character(240)            :: msg
      character(len=*),optional :: status, action, form, access, position
      character(11)             :: status1,action1,form1,access1,position1,number
      integer_dp,      optional :: recl
      integer_dp                :: recl1
      logical,         optional :: numbered
      logical                   :: numbered1,ldum
      integer                   :: i

      status1='old'; action1='readwrite'; form1='formatted'; access1='sequential'; position1='rewind' ! defaults (identical to Fortran standard except status and position)
      numbered1=.false.
      recl1  =0
      if(present(status))   status1   = status
      if(present(action))   action1   = action
      if(present(form))     form1     = form
      if(present(access))   access1   = access
      if(present(position)) position1 = position
      if(present(numbered)) numbered1 = numbered
      if(present(recl)) then
        recl1=recl
      else if (access1=='direct') then
        Error('No record length given for direct-access file (error in code).')
      endif
      if(present(name)) then
        if(status1=='scratch') Error('Filename given for scratch file.')
        name1 = name
      else
        if(status1/='scratch') Error('No filename given.')
        name1 = 'scratch file'
      endif
      if(all(unitopen)) Error('No unit available.')
      fopen = MINUNIT
      do while(unitopen(fopen))
        fopen = fopen + 1 ; if(fopen>MAXUNIT) Error('Reached unit count of 99.')
      enddo
 1    if(numbered1) then
        if(action1/='write'.and.action1/='readwrite') Bug('Cannot open "numbered" file for reading only.')
        i = 0
        do
          if(i<1000) then ; write(number,'(I3.3)') i
          else            ; write(number,'(I11)')  i ; number = adjustl(number)
          endif
          write(name1,'(A,I3.3)') name//'.'//trim(number)
          inquire(file=name1,exist=ldum)
          if(.not.ldum) exit
          i = i + 1
        enddo
      endif
      ! Some compilers don't like recl if access is sequential.
      if(access1=='direct') then
        if(present(name)) then
          open(fopen,file=name1,
     &         status=status1,action=action1,form=form1,access=access1,position=position1,recl=recl1,iostat=i,iomsg=msg)
        else
          open(fopen,
     &         status=status1,action=action1,form=form1,access=access1,position=position1,recl=recl1,iostat=i,iomsg=msg)
        endif
      else
        if(present(name)) then
          open(fopen,file=name1,
     &         status=status1,action=action1,form=form1,access=access1,position=position1,iostat=i,iomsg=msg)
        else
          open(fopen,
     &         status=status1,action=action1,form=form1,access=access1,position=position1,iostat=i,iomsg=msg)
        endif
      endif
      if(i/=0) then
        if(numbered1) goto 1 ! we have probably hit a race condition -> retry
        write(0,'(/2A)')  'fopen: Could not open ',trim(name1)
        write(0,'(A,I5)') '       Error code is',i
        write(0,'(A)')    '       status='//trim(status1)//', action='//trim(action1)//', form='//trim(form1)//
     &                         ', access='//trim(access1)//', position='//trim(position1)
        Error('Error while trying to open a file: '//trim(msg))
      endif
      unitopen(fopen)=.true.

      end function fopen

c     ------------

      subroutine fclose(iunit,status)

      use, intrinsic :: iso_fortran_env
      implicit none
      integer, intent(in)       :: iunit
      integer                   :: i
      character(240)            :: msg
      character(len=*),optional :: status

      if(present(status)) then ; close(iunit,iostat=i,iomsg=msg,status=status)
      else                     ; close(iunit,iostat=i,iomsg=msg)
      endif
      if(i/=0) then
        write(0,'(/A,I3)') 'fclose: Could not close unit',iunit
        write(0,'(A,I5)')  '        Error code is',i
        Error('Error while trying to close a file: '//trim(msg))
      endif
      if(iunit<MINUNIT.or.iunit>MAXUNIT) Bug('Unit out of range.')
      unitopen(iunit)=.false.

      end subroutine fclose

c     ------------

      end module file

