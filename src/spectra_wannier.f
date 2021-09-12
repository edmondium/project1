c--------------------------------------------------------------------------------
c Copyright (c) 2021 Peter Grünberg Institut, Forschungszentrum Jülich, Germany
c This file is part of SPEX and available as free software under the conditions
c of the MIT license as expressed in the LICENSE file in more detail.
c--------------------------------------------------------------------------------

# include "cppmacro.h"
# include "jobtype.h"
# include "restype.h"      

      subroutine spectra_wannier(job1)

      use global
      use wrapper
      use util
      use file
      Mpi( use Mwrapper )

      use, intrinsic :: iso_fortran_env
      implicit none
      type(jobtype), intent(in) :: job1
      complex_dp,  allocatable  :: frq(:)
      complex_dp,  allocatable  :: screenw(:,:,:,:)
      complex_dp,  allocatable  :: suscepw(:,:,:,:,:)
      real_dp                   :: kp(3)
      integer                   :: spin,nfrq,ikpt,ikpt_job
      logical                   :: ldum
      character(14)             :: fname
      character(7)              :: cjob

      if(job1%type/=J_SUSR.and.job1%type/=J_GOLD.or.job1%kern/='BSE') return

      Rwrite(6,'(//A)') '### subroutine: spectra_wannier ###'

      select case(job1%type)
        case(J_SUSR) ; Rwrite(6,'(/A)') '< Renormalized magnetic susceptibility '//trim(job1%kern)//' >'
        case(J_GOLD) ; Rwrite(6,'(/A)') '< Renormalized magnetic susceptibility '//trim(job1%kern)//' (Goldstone mode) >'
        case default ; Error('job number unknown.')
      end select

c
c     Define frequencies
      Rbegin
      nfrq = 0 ; if(allocated(frq)) deallocate ( frq )
      call getfreq(kp, nfrq,job1) ; allocate ( frq(nfrq) ) ! get nfrq and allocate
      call getfreq(frq,nfrq,job1)                          ! get frq
      Rend
      Mpi( call Mcast(nfrq); call Mcastl(frq) )

      spin = job1%spin(1)

c
c     Read Wannier Coulomb matrix
      call checkmem('suscepw',16d0*nwan**4*nfrq)
      call checkmem('screenw',16d0*nwan**4)
      allocate ( suscepw(nwan,nwan,nwan,nwan,nfrq) )
      allocate ( screenw(nwan,nwan,nwan,nwan)      )

      Rcall readcou(ikpt,suscepw,screenw,.false.,.true.,spin,[0,0,0],1,(0d0,0d0),1)
      Rcall reorder4(screenw)
      Rif(ikpt==0)     Error('spex.cou not found.')
      Rif(ikpt/=nkpti) Warn('W incomplete or calculated with different k-point set.')
      Mpi( call Mcast(screenw) )

      ! Loop over k points
      do ikpt_job = 1,size(job1%kpt)

        ikpt = job1%kpt(ikpt_job)
        if(ikpt==0) cycle

        kp = kpt(:,ikpt)

        if(maxval(job%indx)>1) then ; cjob = chr(job1%indx,'I3.3')
        else                        ; cjob = ' '
        endif

        Rwrite(6,'(/A)')             '========='
        Rwrite(6,'(/A,I4,A,3F10.5)') 'K point (',ikpt,' ):',kp
        if(size(job1%kpt)>1) cjob = trim(cjob)//'.'//chr(ikpt_job,'I3.3')
        if(iand(restart,R_spec)/=0) then
          ldum = .true.
          select case(job1%type)
            case(J_SUS)  ; fname = 'suscep'//cjob
            case(J_SUSR) ; fname = 'suscepR'//cjob
            case(J_DIEL) ; fname = 'dielec'//cjob
            case(J_SCR)  ; fname = 'screen'//cjob
            case(J_SCRW) ; ldum  = .false.
            case(J_GOLD) ; ldum  = .false.
          end select
          Rif(ldum) inquire(file=fname,exist=ldum)
          Mpi( call Mcast(ldum) )
          if(ldum) then
            Rwrite(6,'(/A)') 'Output file '//trim(fname)//' exists. Skipping...'
            cycle
          endif
        endif

c
c       Calculate susceptibility in Wannier basis for spin-wave calculations
        call susceptibility_wannier(ikpt,spin,allocated(irrep_wan),nfrq,frq,suscepw)

# if 0
        allocate ( sus(nwan**2,nwan**2) )
        do i = 1,nfrq
          sus = reshape ( suscepw(:,:,:,:,i) , [ nwan**2,nwan**2 ] )
          sus = (sus + transpose(conjg(sus))) / 2
c          write(*,*) 'sumimag',sum(abs(imag(sus)))
          write(*,*) real ( sum ( [ (sus(j,j),j=1,nwan**2) ] ) )
          call diagonalize(eig,sus)
          write(500,'(26F20.10)') imag(frq(i)),eig
        enddo
        Error(' ')
# endif

c        write(400,'(3F20.15)') (imag(frq(i)),sum( [ ((suscepw(n2,n1,n2,n1,i),n1=1,nwan),n2=1,nwan) ] ),i=1,nfrq)
c        call susceptibility_wannier(76,7-spin,.false.,nfrq,frq,suscepw)
c        write(401,'(3F20.15)') (imag(frq(i)),sum( [ ((suscepw(n1,n2,n1,n2,i),n1=1,nwan),n2=1,nwan) ] ),i=1,nfrq)
c        Error(' ')

c
c       Solve Bethe-Salpeter equation
        call bethesalpeter(suscepw,screenw,ikpt,frq,nfrq,spin,cjob)

      enddo ! end loop over k points

      call checkmem('suscepw',-16d0*nwan**4*nfrq)
      call checkmem('screenw',-16d0*nwan**4)
      deallocate ( suscepw,screenw )

      end

c --------------



