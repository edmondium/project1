c--------------------------------------------------------------------------------
c Copyright (c) 2021 Peter Grünberg Institut, Forschungszentrum Jülich, Germany
c This file is part of SPEX and available as free software under the conditions
c of the MIT license as expressed in the LICENSE file in more detail.
c--------------------------------------------------------------------------------

# ifdef LOAD

c     --------

c
c (0) Read bands band(:) at kpt1_in(:) (by calling read_wavef) and transform if necessary.
      subroutine read_wavef0(band,kpt1_in,ispin,cmtout,cpwout,lstr)
      use, intrinsic :: iso_fortran_env
      implicit none
      integer,     intent(in)            :: band(:),kpt1_in(:),ispin
      logical,     intent(in),  optional :: lstr
      complex_dp,  intent(out), optional :: cmtout(:,:,:,:,:)
      MCOMPLEX_dp, intent(out), optional :: cpwout(:,:,:,:)
      complex_dp,  allocatable           :: cmthlp(:,:,:)
      MCOMPLEX_dp, allocatable           :: cpwhlp(:,:)
      logical                            :: lmt,lpw,lbcnt
      integer,     allocatable           :: kpt1(:)
      integer                            :: str
      integer                            :: i,k,ib,b1,b2,k1,k2,nb,mb,ikpt1,ikptp,nk
      str  = 1             ; Mpi( if(present(lstr)) str = Msize )
      nk   = size(kpt1_in) ; allocate(kpt1(nk))
      kpt1 = kptp(kpt1_in) ; where(kpt1>nkpt) kpt1 = kpt1 - nkpt + nkpti
      k1   = 0
      do k = 1,nk
        if(all(kpt1(k)/=kpt1(:k1))) then
          k1       = k1 + 1
          kpt1(k1) = kpt1(k)
        endif
      enddo
      nk    = k1
      lmt   = present(cmtout) ; if(lmt) cmtout = 0
      lpw   = present(cpwout) ; if(lpw) cpwout = 0
      lbcnt = all ( [( band(i+1)-band(i)==1 , i=1,size(band)-1 )] )
      if(lbcnt) then
        b1 = band(1)
        b2 = band(size(band))
      else if(str/=1) then
        Error('not implemented: str/=1 and band(:).')
      endif
      if(.not.(lmt.or.lpw)) Bug('neither cmtout nor cpwout given.')
      if     (.not.lpw) then ; call read_wavef(band,kpt1(:nk),ispin,cmtout,       str=str)
      else if(.not.lmt) then ; call read_wavef(band,kpt1(:nk),ispin,cpwout=cpwout,str=str)
      else                   ; call read_wavef(band,kpt1(:nk),ispin,cmtout,cpwout,str=str)
      endif
      ! get maximum band dimension      
      if(lbcnt) then
        mb = 0
        do k = 1,nk
          ikpt1 = kpt1_in(k)
          nb    = min( nband(ikpt1,ispin) , b2 ) - b1 + 1 ; if(str/=1) nb = (nb-1) / str + 1          
          mb    = max(mb,nb)
        enddo
      else
        mb = size(band)
      endif
      ! predefine cmtout/cpwout with parents of each kpt1
      if(lmt) then
        allocate(cmthlp(maxlmindx,nk,nspin3))
        do ib = 1,mb
          do i = 1,ncent
            cmthlp = cmtout(:,i,ib,:nk,:)
            do k = 1,size(kpt1_in)
              ikptp = kptp(kpt1_in(k))
              if(ikptp>nkpt) ikptp = ikptp - nkpt + nkpti
              if(count(kpt1(:nk)==ikptp)/=1) Bug('More than one or no parent.')
              do k1 = 1,nk
                if(kpt1(k1)==ikptp) cmtout(:,i,ib,k,:) = cmthlp(:,k1,:)
              enddo
            enddo
          enddo
        enddo
        deallocate(cmthlp)
      endif
      if(lpw) then
        allocate(cpwhlp(nk,nspin3))
        do ib = 1,mb
          do i = 1,maxgpt
            cpwhlp = cpwout(i,ib,:nk,:)
            do k = 1,size(kpt1_in)
              ikptp = kptp(kpt1_in(k))
              if(ikptp>nkpt) ikptp = ikptp - nkpt + nkpti
              if(count(kpt1(:nk)==ikptp)/=1) Bug('More than one or no parent.')
              do k1 = 1,nk
                if(kpt1(k1)==ikptp) cpwout(i,ib,k,:) = cpwhlp(k1,:)
              enddo
            enddo
          enddo
        enddo
        deallocate(cpwhlp)
      endif
      ! transform all equivalents
# ifdef old_trafo
      if(.not.lmt) allocate(cmthlp(maxlmindx,ncent,1)) ! dummy array
      if(.not.lpw) allocate(cpwhlp(maxgpt,1))          ! dummy array
# else
      if(lmt) allocate(cmthlp(1,1,1)) ! dummy array
      if(lpw) allocate(cpwhlp(1,1))   ! dummy array
# endif
      do k = 1,size(kpt1_in)
        ikpt1 = kpt1_in(k)
        if(ikpt1/=kptp(ikpt1)) then
          do ib = 1,mb
# ifndef old_trafo
            if(lmt) then
              call waveftrafo_mt_io(cmtout(:,:,ib,k,1),cmthlp,maxlmindx,0,0,symkpt(ikpt1),-kptp(ikpt1))
              if(l_soc) then
                call waveftrafo_mt_io(cmtout(:,:,ib,k,2),cmthlp,maxlmindx,0,0,symkpt(ikpt1),-kptp(ikpt1))
                call waveftrafo_soc(cmtout(:,:,ib,k,:),maxlmindx*ncent,symkpt(ikpt1))
              endif
            endif
            if(lpw) then
              call waveftrafo_pw_io(cpwout(:,ib,k,1),cpwhlp,symkpt(ikpt1),-kptp(ikpt1))
              if(l_soc) then
                call waveftrafo_pw_io(cpwout(:,ib,k,2),cpwhlp,symkpt(ikpt1),-kptp(ikpt1))
                call waveftrafo_soc(cpwout(:,ib,k,:),maxgpt,symkpt(ikpt1))
              endif
            endif
# else
            if     (.not.lpw) then ; call waveftrafo_io(cmtout(:,:,ib,k,:),cpwhlp,          ikpt1,1)
            else if(.not.lmt) then ; call waveftrafo_io(cmthlp,            cpwout(:,ib,k,:),ikpt1,2)
            else                   ; call waveftrafo_io(cmtout(:,:,ib,k,:),cpwout(:,ib,k,:),ikpt1,3)
            endif
# endif
          enddo
        endif
      enddo
# if old_trafo
      if(.not.lmt) deallocate(cmthlp)
      if(.not.lpw) deallocate(cpwhlp)
# else
      if(lmt) deallocate(cmthlp)
      if(lpw) deallocate(cpwhlp)
# endif
      end subroutine read_wavef0

c     --------

c (1) Read bands bo1..bo2 at kpt1 and bu1..bu2 at kpt1+ikpt
c
c     lmt = .true. : read cmt coefficients
c     lpw = .true. : read cpw coefficients
c
c     storeibz = .false :      
c       -> cmt(:maxlmindx,:ncent,      :nbo     ,:nkpt1,:nspin3)      ! nbo = bo2 - bo1 + 1 ; k-points according to kpt1(:)
c          cmt(:maxlmindx,:ncent, nbo+1:nbo+nbu ,:nkpt1,:nspin3)      ! nbu = bu2 - bu1 + 1 ; k-points according to kpt1(:)+ikpt
c       -> cpw(:maxgpt,                :nbo     ,:nkpt1,:nspin3)
c          cpw(:maxgpt,           nbo+1:nbo+nbu ,:nkpt1,:nspin3)
c
c     storeibz = .true. :
c       -> cmt(:maxlmindx,:ncent,      :nbo     ,:nk   ,:nspin3)      ! nk = nkpti if ikpt<=nkpt, nk = nkpti2 otherwise
c          cmt(:maxlmindx,:ncent, nbo+1:nbo+nbu ,:nk   ,:nspin3)
c       -> cpw(:maxgpt,                :nbo     ,:nk   ,:nspin3)
c          cpw(:maxgpt,           nbo+1:nbo+nbu ,:nk   ,:nspin3)
c
      subroutine read_wavef1(bo1,bo2,bu1,bu2,kpt1,nkpt1,ikpt,ispin1,ispin2,lmt,lpw)
      use, intrinsic :: iso_fortran_env
      implicit none
      integer, intent(in)      :: bo1,bo2,bu1,bu2,nkpt1,kpt1(nkpt1),ikpt,ispin1,ispin2
      logical, intent(in)      :: lmt,lpw
      integer, allocatable     :: kpt1_(:)
      integer                  :: b1,b2
      integer                  :: i,k,ikptp,ikpt1,ikpt2,ib,k1,nk
      complex_dp,  allocatable :: cmthlp(:,:,:),cmth(:,:,:,:,:)
      MCOMPLEX_dp, allocatable :: cpwhlp(:,:),cpwh(:,:,:,:)
      integer                  :: kptsum
      if(.not.(lmt.or.lpw)) Bug('lmt and lpw both false.')
      if(nkpt1<nkpti)       Bug('nkpt1 < nkpti.')
      if(any(kpt1>nkpt))    Bug('kpt1 out of range.')
      if(any([bo1,bu1]<1).or.any([bo2,bu2]>maxband).or.any([bo1,bu1]>[bo2,bu2]))
     &                      Bug('Band indices (bo1,bo2,bu1,bu2) out of range.')
      if(storeibz) then
        if(ikpt<=nkpt) then ; nk = nkpti
        else                ; nk = nkpti2
        endif
        allocate(kpt1_(nkpti))
        kpt1_ = [(k,k=1,nkpti)]
        if(lmt) allocate ( cmt(maxlmindx,ncent,bu2-bu1+bo2-bo1+2,nk,nspin3) )
        if(lpw) allocate ( cpw(maxgpt,         bu2-bu1+bo2-bo1+2,nk,nspin3) )
        kindx(:nkpt) = kptp(:nkpt)
        if(nkpt2>nkpt) kindx(nkpt+1:) = kptp(nkpt+1:) - nkpt
      else
        allocate(kpt1_(nkpt1))
        kpt1_ = kpt1
        do k = 1,nkpti
          if(all(kpt1/=k)) Bug('irreducible k point not included in kpt1 list.')
        enddo
        if(lmt) allocate ( cmt(maxlmindx,ncent,bu2-bu1+bo2-bo1+2,nkpt1,nspin3) )
        if(lpw) allocate ( cpw(maxgpt,         bu2-bu1+bo2-bo1+2,nkpt1,nspin3) )
        if(lmt) allocate ( cmth(maxlmindx,ncent,bu2-bu1+bo2-bo1+2,nkpt1,nspin3) )
        if(lpw) allocate ( cpwh(maxgpt,         bu2-bu1+bo2-bo1+2,nkpt1,nspin3) )        
        kindx(kpt1) = [(k,k=1,nkpt1)]
      endif
      if(lmt) cmt = 0
      if(lpw) cpw = 0
      if(lmt) cmth = 0
      if(lpw) cpwh = 0
c     Occupied states (bo1,bo2)
      b2 = bo2 - bo1 + 1
      if     (.not.lpw) then ; call read_wavef0([(i,i=bo1,bo2)],kpt1_,ispin1,cmt(:,:,:b2,:,:))
      else if(.not.lmt) then ; call read_wavef0([(i,i=bo1,bo2)],kpt1_,ispin1,        cpwout = cpw(:,:b2,:,:))
      else                   ; call read_wavef0([(i,i=bo1,bo2)],kpt1_,ispin1,cmt(:,:,:b2,:,:),cpw(:,:b2,:,:))
      endif
c     Unoccupied states (bu1,bu2)
      if(bu2-bu1>=0) then
        b1 = bo2 - bo1 + 2
        b2 = bu2 - bu1 + b1
        if(ikpt<=nkpt) then ; k1 = 1       ; nk = nkpti
        else                ; k1 = nkpti+1 ; nk = nkpti2
        endif
        if(storeibz.and.ikpt>nkpt) then
          deallocate(kpt1_)
          allocate(kpt1_(nkpti2))
          kpt1_ = [(k,k=nkpt+1,nkpt+nkpti2)]
        else if(.not.storeibz) then
          do k = 1,size(kpt1_)
            kpt1_(k) = kptsum(kpt1_(k),ikpt)
          enddo
        endif
        if     (.not.lpw) then ; call read_wavef0([(i,i=bu1,bu2)],kpt1_,ispin2,cmt(:,:,b1:b2,:,:))
        else if(.not.lmt) then ; call read_wavef0([(i,i=bu1,bu2)],kpt1_,ispin2,          cpwout = cpw(:,b1:b2,:,:))
        else                   ; call read_wavef0([(i,i=bu1,bu2)],kpt1_,ispin2,cmt(:,:,b1:b2,:,:),cpw(:,b1:b2,:,:))
        endif
      endif
      deallocate(kpt1_)
      end subroutine read_wavef1

c     --------

c (2) Read bands band(:) at ikpt,ispin -> cmtout,cpwout (included to allow for array structure with single k point)
      subroutine read_wavef2(band,n,ikpt,ispin,cmtout,cpwout)
      use, intrinsic :: iso_fortran_env
      implicit none
      integer,     intent(in)            :: n,band(n),ikpt,ispin
      complex_dp,  intent(out), optional :: cmtout(maxlmindx,ncent,n,1,nspin3)
      MCOMPLEX_dp, intent(out), optional :: cpwout(maxgpt,n,1,nspin3)      
      integer                            :: ikptp,i
      complex_dp                         :: cdum
      if     (.not.present(cpwout)) then ; call read_wavef0(band,[ikpt],ispin,cmtout)
      else if(.not.present(cmtout)) then ; call read_wavef0(band,[ikpt],ispin,cpwout=cpwout)
      else                               ; call read_wavef0(band,[ikpt],ispin,cmtout,cpwout)
      endif
      end subroutine read_wavef2

c     --------

c (3) Read bands b1..b2 at kpt1 -> cmt,cpw  [as (1) without "unocc." states; preserves band indices]
      subroutine read_wavef3(b1,b2,kpt1,nkpt1,ispin)
      Mpi( use Mwrapper )
      use, intrinsic :: iso_fortran_env
      implicit none
      integer, intent(in)           :: b1,b2,nkpt1,kpt1(nkpt1),ispin
      integer                       :: b2_,k,ikindx,ikindx1,ikpt1,ib,k1,i
      integer(HID_T)                :: Hfile
      integer(HSIZE_T), allocatable :: Hdim(:)
      character(8)                  :: ckpt
      complex_dp                    :: cmthlp(maxlmindx,ncent,nspin3)
      MCOMPLEX_dp                   :: cpwhlp(maxgpt,nspin3)
      Hpos = 0 ; if(l_soc) Hpos = 1
      if(any(kpt1>nkpt)) Error('k-point index of shifted set.')
      if(storeibz) then
        allocate ( cmt(maxlmindx,ncent,b1:b2,nkpti,nspin3) ) ; cmt = 0
        allocate ( cpw(maxgpt,         b1:b2,nkpti,nspin3) ) ; cpw = 0
        kindx(:nkpt) = kptp(:nkpt) ; if(nkpt2>nkpt) kindx(nkpt+1:) = kindx(nkpt+1:) - nkpt + nkpti
        call read_wavef0([(i,i=b1,b2)],[(k,k=1,nkpti)],ispin,cmt,cpw)
      else
        allocate ( cmt(maxlmindx,ncent,b1:b2,nkpt1,nspin3) ) ; cmt = 0
        allocate ( cpw(maxgpt,         b1:b2,nkpt1,nspin3) ) ; cpw = 0
        kindx(kpt1) = [(k,k=1,nkpt1)]
        call read_wavef0([(i,i=b1,b2)],kpt1,ispin,cmt,cpw)
      endif
      end subroutine read_wavef3

c     --------

# endif

c
c     Routines to reorder the k-point list (kpt) and related arrays to simplify the writing of the read-in and write-out routines.
c     (Do not change.)      
c
c     kpt_reorder1 : reorders kpt(:,:) (and other arrays) to have all irreducible kpoints at the beginning [kpt(:,:nkpti+nkpti2)]
c     kpt_reorder2 : undoes reordering (->normal order)
c       mode == 0 : only reorder kpt(:,:)
c       mode == 1 : also reorder ngpt(:)
c       mode == 2 : also reorder nband(:,:)
c       mode == 3 : also reorder ngpt(:) and nband(:,:)
c
      subroutine kpt_reorder1(mode)
      use, intrinsic :: iso_fortran_env
      implicit none
      integer, intent(in) :: mode
      real_dp             :: vec(3)
      integer             :: ihlp(nspin2)
      integer             :: i,idum
      if(all(mode/=[0,1,2,3])) Bug('unknown mode')
      if(kpt_order/=0)         Bug('kpt is already reordered.')
      do i = 1,nkpti2
        vec            = kpt(:,nkpti+i)
        kpt(:,nkpti+i) = kpt(:,nkpt+i)
        kpt(:,nkpt+i)  = vec
        if(iand(mode,1)/=0) then
          if(.not.allocated(ngpt))  Bug('not allocated: ngpt')
          idum             = ngpt(nkpti+i)
          ngpt(nkpti+i)    = ngpt(nkpt+i)
          ngpt(nkpt+i)     = idum
        endif
        if(iand(mode,2)/=0) then
          if(.not.allocated(nband)) Bug('not allocated: nband')
          ihlp             = nband(nkpti+i,:)
          nband(nkpti+i,:) = nband(nkpt+i,:)
          nband(nkpt+i,:)  = ihlp
        endif       
      enddo
      kpt_order = mode + 1
      end subroutine kpt_reorder1

      subroutine kpt_reorder2(mode)
      use, intrinsic :: iso_fortran_env
      implicit none
      integer, intent(in) :: mode
      real_dp             :: vec(3)
      integer             :: ihlp(nspin2)
      integer             :: i,idum
      if(all(mode/=[0,1,2,3]))          Bug('unknown mode')
      if     (kpt_order==0)      then ; Bug('kpt is not reordered.')
      else if(kpt_order/=mode+1) then ; Bug('wrong mode or kpt wrongly ordered')
      endif
      do i = nkpti2,1,-1
        vec            = kpt(:,nkpt+i)
        kpt(:,nkpt+i)  = kpt(:,nkpti+i)
        kpt(:,nkpti+i) = vec
        if(iand(mode,1)/=0) then
          if(.not.allocated(ngpt))  Bug('not allocated: ngpt')
          idum             = ngpt(nkpt+i)
          ngpt(nkpt+i)     = ngpt(nkpti+i)
          ngpt(nkpti+i)    = idum
        endif
        if(iand(mode,2)/=0) then
          if(.not.allocated(nband)) Bug('not allocated: nband')
          ihlp             = nband(nkpt+i,:)
          nband(nkpt+i,:)  = nband(nkpti+i,:)
          nband(nkpti+i,:) = ihlp
        endif
      enddo
      kpt_order = 0
      end subroutine kpt_reorder2

c     --------
      
