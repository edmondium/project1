c--------------------------------------------------------------------------------
c Copyright (c) 2021 Peter Grünberg Institut, Forschungszentrum Jülich, Germany
c This file is part of SPEX and available as free software under the conditions
c of the MIT license as expressed in the LICENSE file in more detail.
c--------------------------------------------------------------------------------

c Define symmetry arrays.
c
c sym(i)%rrot  : rotation matrix wrt reciprocal lattice vectors [ = transpose ( sym(sym(i)%inv)%rot ) ]
c sym(i)%inv   : inverse symmetry operation
c sym(i)%trs   : true if symmetry operation involves time-reversal
c pcent(c,i)   : / rotation of atom center c through operation i leads to
c tcent(:,c,i) : \ cent(:,pcent(c,i)) + tcent(:,c,i)  [ tcent is a lattice vector ]
c gensym(i)    : array of generators
c symgen(:,i)  : generators of symmetry operation i: sym(i) = gensym(symgen(1,i)) * gensym(symgen(2,i)) * ...
c
c If the system has no inversion symmetry, combinations with time reversal symmetry are appended.
c
# include "cppmacro.h"
      subroutine def_symmetry
      use readwrite, only: invs
      use util,      only: reallocate,chr
      use global,    only: pi,l_soc,trsoff,invsym,sym,nsym,nsymt,pcent,tcent,symtab,gensym,ngensym,symgen,ntype,neq,cent,ncent,
     &                     symcent,symtype,lat,rlat,modulo1r
      use, intrinsic :: iso_fortran_env
      implicit none
      type(symtype), allocatable :: sym1(:)
      real_dp                    :: kvec(3),axis(3),rvec(3)
      real_dp                    :: det,angle
      integer,       allocatable :: iarr(:)
      integer                    :: nsym1
      integer                    :: itype,ineq,ineq1,icent,icent1,isym
      integer                    :: i,j,k,l
      invsym = 0
      do i = 1,nsym
        if(all(sym(i)%rot==reshape([-1,0,0,0,-1,0,0,0,-1],[3,3])).and.all(sym(i)%transl==0)) invsym = i
      enddo
      nsym1   = nsym
      nsymt   = nsym
      sym%trs = .false.
      call determine_symmetry(0)
      if(invs.and.invsym==0) Error('Did not find index of inversion symmetry operation.')
      if((invsym==0.or.l_soc).and..not.trsoff) then ! ".not.trsoff" only needed if we change the default above
        allocate   ( sym1(nsym) ) ; sym1         = sym
        deallocate ( sym )        ; nsym         = 2*nsym
        allocate   ( sym(nsym) )  ; sym(:nsym/2) = sym1
      else if(invsym/=0) then
        if(.not.invs) then
          Warn('Found inversion symmetry operation. You seem to use an executable compiled without -DINV')
          write(0,'(A)') '          on a system with inversion symmetry. You should better use an executable compiled with -DINV.'
          write(0,'(A)') '          Anyway, time reversal symmetry will not be used.'
        endif
        trsoff = .true.
      endif
      allocate ( pcent(ncent,nsym),tcent(3,ncent,nsym),symtab(nsym,nsym) )
      if(any(sym(1)%rot-reshape([1,0,0,0,1,0,0,0,1],[3,3])/=0))
     &  Error('First symmetry operation is not the identity.')
      ! Write out properties of symmetry operations
      write(6,'(/A)')   'Symmetry information:'
      write(6,'(A,I3)') '  Number of symmetry operations =',nsym
      write(6,'(A)')    '   #  Det           Axis           Angle   Syp  Inv'
      symtab = 0
      do i = 1,nsym1
        write(6,'(2X,I2,''  '''NoA) i
        call rotprop(det,angle,axis,i)
        if(angle==0) then
          write(6,'(I3,1X,3A7,A9'NoA) nint(det),'    -- ','    -- ','    -- ','     --  '
        else
          axis = axis / maxval(abs(axis)) ! originally: rarr = rarr / abs(rarr(minloc(abs(rarr),1,abs(rarr)>1d-10)))
          write(6,'(I3,1X,3F7.2,F9.2'NoA) nint(det),axis,angle/pi*180
        endif
        if(any(sym(i)%transl/=0)) then
          write(6,'(3X,''No '''NoA)
        else
          write(6,'(3X,''Yes'''NoA)
        endif
        sym(i)%inv = 0
        do j = 1,nsym1
          do k = 1,nsym1
            kvec = matmul(sym(i)%rot,sym(j)%transl) + sym(i)%transl - sym(k)%transl
            if(all(matmul(sym(i)%rot,sym(j)%rot)==sym(k)%rot).and.sum(abs(kvec-nint(kvec)))<1d-10) then
              if(symtab(i,j)/=0) Error('Error in symmetry multiplication table.')
              symtab(i,j) = k
              if(k==1) then
                write(6,'(I4)') j
                if(sym(i)%inv/=0) Error('inverse operation already defined.')
                sym(i)%inv  = j
                sym(i)%rrot = transpose ( sym(j)%rot ) ! temporary fix for ifc
              endif
            endif
          enddo

          if(symtab(i,j)==0) Error('Could not define symmetry multiplication table.')
          if(all(matmul(sym(i)%rot,sym(j)%rot)==reshape([1,0,0,0,1,0,0,0,1],[3,3]))
     &        .and.all(modulo1r(matmul(sym(i)%rot,sym(j)%transl)+sym(i)%transl)<1d-10)) then
            if(j/=sym(i)%inv) Error('error.')
          endif
        enddo

        if(sym(i)%inv==0) Error('inverse operation not found.')
        ! Check symmetry operations
        icent = 0
        do itype = 1,ntype
          do ineq = 1,neq(itype)
            icent          = icent + 1
            rvec           = matmul(sym(i)%rot,cent(:,icent)) + sym(i)%transl
            icent1         = sum(neq(:itype-1))
            pcent(icent,i) = 0
            do ineq1 = 1,neq(itype)
              icent1 = icent1 + 1
              if( all( abs( rvec-cent(:,icent1)-nint(rvec-cent(:,icent1)) ) <1d-10 ) ) then
                if(pcent(icent,i)/=0) Error('Double atom mapping.')
                pcent(icent,i)   = icent1
                tcent(:,icent,i) = nint ( rvec-cent(:,icent1) )
                if(maxval(abs(cent(:,icent1)+tcent(:,icent,i)-rvec))>1d-6)
     &            Error('Could not define atom mapping.')
              endif
            enddo
            if(pcent(icent,i)==0) then
              write(6,'(A,I3,A)') 'getinput: Symmetry operation',i,' is inconsistent with the crystal structure.'
              write(6,'(A,I3,A)') '          Mapping of atom',icent,' failed.'
              Error('Symmetry operation inconsistent with crystal structure.')
            endif
          enddo
        enddo
      enddo

      ! Define symcent (symcent(ic): symmetry that rotates to ic from parent center)
      allocate(symcent(ncent))
      symcent = 0
      icent   = 1
      do itype = 1,ntype
        do isym = 1,nsym1
          icent1 = pcent(icent,isym)
          if(icent1<1.or.icent1>ncent) Error('pcent pointer incorrect.')
          if(symcent(icent1)==0) symcent(icent1) = isym
        enddo
        icent = icent + neq(itype)
      enddo
      if(any(symcent==0)) Error('Definition of array symcent failed.')

      ! Check sym()%rrot
      do j = 1,nsym1
        if( sum ( abs ( matmul( matmul(transpose( lat), lat) / (2*pi) , matmul ( sym(j)%rot,
     &                          matmul(transpose(rlat),rlat) / (2*pi) ) ) - sym(j)%rrot ) ) > 1d-12 ) then
          Error('Non-integer matrix element in sym()%rrot. Symmetry broken. Check geometry.')
        endif
      enddo

      ! Define combinations with time reversal symmetry
      if(.not.trsoff) then
        write(6,'(2X,I2,A,I2,A)') nsym/2+1,' - ',nsym,'  Time reversal symmetries'
        do isym = 1,nsym/2
          sym(nsym/2+isym)%transl     =  sym(isym)%transl
          sym(nsym/2+isym)%rot        =  sym(isym)%rot
          sym(nsym/2+isym)%rrot       = -sym(isym)%rrot
          sym(nsym/2+isym)%inv        =  sym(isym)%inv + nsym/2
          sym(nsym/2+isym)%trs        =  .true.
          pcent(  :,nsym/2+isym)      =  pcent(  :,isym)
          tcent(:,:,nsym/2+isym)      =  tcent(:,:,isym)
          sym(nsym/2+isym)%esoc       =  matmul ( reshape ( [ 0,1,-1,0 ] , [2,2] ) , sym(isym)%esoc )
          symtab(nsym/2+1:,nsym/2+1:) =  symtab(:nsym/2,:nsym/2)
          symtab(nsym/2+1:,:nsym/2)   =  symtab(:nsym/2,:nsym/2) + nsym/2
          symtab(:nsym/2,nsym/2+1:)   =  symtab(:nsym/2,:nsym/2) + nsym/2
        enddo
      endif

      ! Determine generators
      allocate(iarr(nsym),gensym(nsym),symgen(nsym,nsym))
      symgen   = 0
      j        = 0
      iarr(1)  = 1 ! identity operator is not a generator
      iarr(2:) = 0
      do while(any(iarr==0))
        isym      = minloc(iarr,1)
        j         = j + 1
        gensym(j) = isym                    
        i         = isym
        k         = 1
        do
          if(iarr(i)==0) then
            iarr(i)      = 1
            symgen(:k,i) = [ (j,l=1,k) ]
          endif
          k = k + 1
          i = symtab(i,isym)
          if(i==isym) exit
        enddo
        do i = 1,nsym
          if(iarr(i)==1) then
            k = symtab(isym,i)
            if(iarr(k)==0) then
              iarr(k)      = 1
              l            = count( [ symgen(:,i) , symgen(:,isym) ] /= 0 ) ; if(l>size(symgen,1)) Bug('Out of bounds of symgen.')
              symgen(:l,k) = pack ( [ symgen(:,i) , symgen(:,isym) ] , [ symgen(:,i) , symgen(:,isym) ] /= 0 )
            endif
            k = symtab(i,isym)
            if(iarr(k)==0) then
              iarr(k)      = 1
              l            = count( [ symgen(:,isym) , symgen(:,i) ] /= 0 ) ; if(l>size(symgen,1)) Bug('Out of bounds of symgen.')
              symgen(:l,k) = pack ( [ symgen(:,isym) , symgen(:,i) ] , [ symgen(:,isym) , symgen(:,i) ] /= 0 )
            endif
          endif
        enddo
      enddo
      deallocate(iarr)
      ngensym = j
      l       = maxval( [ (count(symgen(:,isym)/=0),isym=1,nsym) ] )
      call reallocate(gensym,ngensym)
      call reallocate(symgen,l,nsym)

      ! Check generators
      do isym = 1,nsym
        k = 1
        do i = 1,size(symgen,1)
          if(symgen(i,isym)==0) exit
          k = symtab( gensym(symgen(i,isym)) , k )
        enddo
        if(k/=isym) Bug('Error in determining symmetry generators.')
      enddo
      if(any( gensym(pack(symgen,symgen/=0)) == 1 )) Bug('Found identity operator among generators.')
      write(6,'(A)') '  Symmetry generators: '//Chn(gensym,', ')
      
      end

c     ------------------

c     Calculates determinant, rotation axis, and rotation angle of symmetry operation isym.
c     Also defines SU2 rotation matrix (only used if l_soc=.true.).
      subroutine rotprop(det,angle,axis,isym)
      use global,  only: pi,lat,rlat,sym,l_soc,sqa,nspin
      use util,    only: chr
      use wrapper, only: identity
      use, intrinsic :: iso_fortran_env
      implicit none
      integer, intent(in)  :: isym
      real_dp, intent(out) :: det,angle,axis(3)
      real_dp              :: rot(3,3),rot1(3,3)
      real_dp              :: determinant
      real_dp              :: cosa,sina,sqaxis(3)
      integer              :: i,j,k
      rot = matmul(lat,matmul(sym(isym)%rot,transpose(rlat)))/(2*pi)
      ! Calculate determinant
      det = determinant(rot)
      if(abs(1-abs(det))>1d-10) Error('determinant neither 1 nor -1')
      if( any ( abs( matmul(transpose(rot),rot)-identity(3) )>1d-10 ) ) then
        write(0,'(/A)')'rotprop: Defective symmetry transformation matrix.'
        write(0,'(A)') '         Not a proper rotation matrix:'
        write(0,'(9X,3F20.15)') rot
        Error('Defective symmetry transformation matrix.')
      endif
      rot1 = rot / det ! rot1 is proper rotation

      cosa = ( rot1(1,1)+rot1(2,2)+rot1(3,3) - 1 ) / 2
      if(abs(1-cosa)<1d-10) then
        angle   = 0
      else
        if(abs(1+cosa)<1d-10) cosa = -1
        i       = maxloc ( [ rot1(1,1),rot1(2,2),rot1(3,3) ] ,1)
        j       = mod(i  ,3) + 1
        k       = mod(i+1,3) + 1
        axis(i) = sqrt ( (rot1(i,i)-cosa) / (1-cosa) )
        sina    = (rot1(k,j)-rot1(j,k))/(2*axis(i))
        if(sina<0) then ! rotation angle defined between 0 and 180 degrees
          sina    = -sina
          axis(i) = -axis(i)
        endif
        axis(j) = (rot1(i,j)+rot1(j,i)) / (2*(1-cosa)*axis(i))
        axis(k) = (rot1(i,k)+rot1(k,i)) / (2*(1-cosa)*axis(i))
        angle   = acos(cosa)
        if(abs(1-sum(axis**2))             >1d-12) Error('Axis vector not normalized.')
        if(sum(abs(matmul(rot1,axis)-axis))>1d-12) Error('Determination of rotation axis failed.')
      endif

      if(nspin==2) then
        sqaxis = [ sin(sqa(1)) * cos(sqa(2)) , sin(sqa(1)) * sin(sqa(2)) , cos(sqa(1)) ]
        if(l_soc.and.any(abs(matmul(rot1,sqaxis)-sqaxis)>1d-12))
     &    Error('Symmetry operation '//trim(chr(isym))//' rotates spin quantization axis.')
      endif
      call rot_su2(sym(isym)%esoc,rot)

      end

c     ------------------

c     Returns SU(2) rotation matrix.
      subroutine rot_su2(su2,rot)
      use global, only: img,sqa
      use wrapper
      use, intrinsic :: iso_fortran_env
      implicit none
      complex_dp, intent(out) :: su2(2,2)
      real_dp,    intent(in)  :: rot(3,3)
      real_dp                 :: rot1(3,3),axis(3)
      real_dp                 :: det,angle,cosa,sina
      real_dp                 :: determinant
      integer                 :: i,j,k
      det = determinant(rot)
      if(abs(1-abs(det))>1d-10) Error('determinant neither 1 nor -1')
      if( any ( abs( matmul(transpose(rot),rot)-identity(3) )>1d-10 ) ) then
        write(0,'(/A)')'rot_su2: Defective symmetry transformation matrix.'
        write(0,'(A)') '         Not a proper rotation matrix:'
        write(0,'(9X,3F20.15)') rot
        Error('Defective symmetry transformation matrix.')
      endif
      if(any(sqa/=0)) then        
        call rotation_sph(rot1,sqa(1),sqa(2))
        rot1 = matmul ( transpose(rot1) , matmul( rot / det , rot1) ) ! rot1 is proper rotation (+trafo to align with sqaxis)
      else
        rot1 = rot / det                                              ! rot1 is proper rotation
      endif

      cosa = ( rot1(1,1)+rot1(2,2)+rot1(3,3) - 1 ) / 2
      if(abs(1-cosa)<1d-10) then
        su2     = reshape ( [ 1,0,0,1 ] , [ 2,2 ] )
      else
        if(abs(1+cosa)<1d-10) cosa = -1
        i       = maxloc ( [ rot1(1,1),rot1(2,2),rot1(3,3) ] ,1)
        j       = mod(i  ,3) + 1
        k       = mod(i+1,3) + 1
        axis(i) = sqrt ( (rot1(i,i)-cosa) / (1-cosa) )
        sina    = (rot1(k,j)-rot1(j,k))/(2*axis(i))
        if(sina<0) then ! rotation angle defined between 0 and 180°
          sina    = -sina
          axis(i) = -axis(i)
        endif
        axis(j) = (rot1(i,j)+rot1(j,i)) / (2*(1-cosa)*axis(i))
        axis(k) = (rot1(i,k)+rot1(k,i)) / (2*(1-cosa)*axis(i))
        angle   = acos(cosa)
        if(abs(1-sum(axis**2))             >1d-12) Error('Axis vector not normalized.')
        if(sum(abs(matmul(rot1,axis)-axis))>1d-12) Error('Determination of rotation axis failed.')
        ! rotated spin eigenvectors in terms of original spin frame (as rows)
        cosa = cos(angle/2)
        sina = sin(angle/2)
        su2  = reshape ( [ cosa - img*axis(3)*sina , (-img*axis(1)-axis(2))*sina ,
     &                 (-img*axis(1)+axis(2))*sina , cosa + img*axis(3)*sina ] , [ 2,2 ] )
      endif
      end

c     ------------------

c
c     Determines system symmetry
c     mode=0 : Checks sym(:). Only spatial symmetries are checked (unless SOC & nspin==2).
c     mode=1 : Defines sym(:). Includes time-reversal symmetries (trs) (order: first w/o trs, then with trs).
      subroutine determine_symmetry(mode)
      use global,  only: pi,sym,nsym,symtype,lat,rlat,cent,ntype,neq,vol,l_soc,nspin,invsym,trsoff
      use util,    only: chr
      use wrapper, only: identity,diagonal
      use, intrinsic :: iso_fortran_env     
      implicit none
      integer, intent(in)        :: mode
      type(symtype), allocatable :: sym1(:)
      integer,       allocatable :: iarr(:)
      integer, parameter         :: nn  = 3                        ! -nn..nn values used for lattice vector search
      integer, parameter         :: div = 2                        ! translational vector may components n/div with n=integer
      real_dp, parameter         :: accuracy(2) = [ 1d-12 , 1d-4 ] ! accuracy for normal and sloppy run
      real_dp                    :: lat1(3,3),rlat1(3,3),rot(3,3),dist(3),trans(3),trans1(3),transv(3,maxval(neq))
      integer                    :: mat(3,3),mat1(3,3),mrot(3,3),latvec(3,3,(2*nn+1)**3),nlatvec(3)
      integer                    :: nsym1,iturn,itrans,ntrans
      integer                    :: ix,iy,iz,itype,ieq,ieq1,ic
      integer                    :: i,j
      complex_dp                 :: su2(2,2)
      character(12)              :: str
      character(5)               :: ctrs(3)      
      logical                    :: found,ltrs
      if(all(mode/=[0,1])) Bug('Unknown mode.')
      ! Backfold lattice vectors (if necessary)
      call backfold(lat1,lat)
      call vprod(rlat1(:,1),lat1(:,2),lat1(:,3))
      call vprod(rlat1(:,2),lat1(:,3),lat1(:,1))
      call vprod(rlat1(:,3),lat1(:,1),lat1(:,2))
      if( abs( dot_product(lat1(:,1),rlat1(:,1)) - vol ) > 1d-12 ) Bug('Routine backfold breaks unit cell.')
      rlat1 = 2*pi * rlat1 / vol
      mat   = nint( matmul ( transpose(rlat)  / (2*pi) , lat1 ) ) ! lat1 expressed in lat
      mat1  = nint( matmul ( transpose(rlat1) / (2*pi) , lat  ) ) ! lat  expressed in lat1
      if( any(abs( matmul(lat,mat)-lat1  ) > 1d-12) ) Bug('Could not determine integer backfolding matrix.')
      if( any(abs( matmul(lat1,mat1)-lat ) > 1d-12) ) Bug('Could not determine integer backfolding matrix.')      
      ! iturn=1 : normal run with accuracy(1), iturn=2 : sloppy run with accuracy(2)
      do iturn = 1,2
        ! determine all Bravais vectors with the same length as lat1(:,i), i=1,2,3
        nlatvec = 0
        do ix = -nn,nn
          do iy = -nn,nn
            do iz = -nn,nn
              if(all([ix,iy,iz]==0)) cycle
              do i = 1,3
                if( abs( sqrt(sum(matmul(lat1,[ix,iy,iz])**2)) - sqrt(sum(lat1(:,i)**2)) ) < accuracy(iturn) ) then
                  nlatvec(i)             = nlatvec(i) + 1
                  latvec(:,i,nlatvec(i)) = [ix,iy,iz]
                endif
              enddo
            enddo
          enddo
        enddo
        ! Run over all combinations of three Bravais vectors to determine possible symmetry operations
        ltrs  = .false.
        nsym1 = 0
        do
          do ix = 1,nlatvec(1)
            do iy = 1,nlatvec(2)
      izloop: do iz = 1,nlatvec(3)
                mrot  = reshape ( [ latvec(:,1,ix) , latvec(:,2,iy) , latvec(:,3,iz) ] , [ 3,3 ] )
                mrot  = matmul(mat, matmul(mrot,mat1) )
                rot   = matmul ( lat , matmul ( mrot , transpose(rlat) / (2*pi) ) ) ! rot = rotation matrix
                ! Check validity of symmetry operation
                if( any(abs( matmul(rot,transpose(rot)) - identity(3) ) > accuracy(iturn)) ) cycle ! skip if rot is not an orthogonal matrix
                ! Check if symmetry operation (mrot,*) is consistent with basis (atom centers)
                ! (1) determine possible translation vectors by rotation of first atom (->transv)
                itrans = 0
                do ieq = 1,neq(1)
                  dist = cent(:,ieq) - matmul( mrot , cent(:,1) )
                  if( any(abs( div*dist - nint(div*dist) ) <= accuracy(iturn) )) then ! this is a valid translation vector
                    itrans           = itrans + 1
                    trans            = dble(nint(div*dist)) / div
                    transv(:,itrans) = trans - nint(trans-1d-6) ! -1d-6 is to avoid 0.5 to turn into -0.5
                  endif
                enddo
                ntrans = itrans
                ! (2) check if any of the transv forms a valid symmorphic operation in combination with mrot                
                trans = huge(0d0)
         tloop: do itrans = 1,ntrans
                  ic = 0
                  do itype = 1,ntype
                    do ieq = 1,neq(itype)
                      found = .false.
                      do ieq1 = 1,neq(itype)
                        dist = cent(:,ic+ieq1) - ( matmul( mrot , cent(:,ic+ieq) ) + transv(:,itrans) )
                        if( all(abs( dist - nint(dist) ) <= accuracy(iturn) ) ) found = .true.
                      enddo
                      if(.not.found) cycle tloop
                    enddo
                    ic = ic + neq(itype)
                  enddo
                  trans = transv(:,itrans)
                  exit
                enddo tloop
                if(trans(1)==huge(0d0)) cycle izloop ! next mrot if no valid translation vector has been found
                ! (3) check if SQA is invariant
                if(l_soc.and.nspin==2) then
                  call rot_su2(su2,rot)
                  if(ltrs) su2 = matmul ( reshape ( [ 0,1,-1,0 ] , [2,2] ) , su2 )
                  if( any( abs(su2-diagonal([su2(1,1),su2(2,2)])) > accuracy(iturn) ) ) cycle izloop
                endif
                ! (4) check for double definitions (should not happen)
                if(allocated(sym1)) then
                  do i = 1,nsym1
                    dist = sym1(i)%transl - trans
                    if( all(sym1(i)%rot==mrot) .and. all(abs(dist-nint(dist))<1d-12) .and. (sym1(i)%trs.eqv.ltrs) )
     &                Error('Double symmetry operation!')
                  enddo
                endif
                ! (5) define symmetry operation
                nsym1 = nsym1 + 1
                if(allocated(sym1)) then
                  sym1(nsym1)%rot    = mrot
                  sym1(nsym1)%transl = trans
                  sym1(nsym1)%trs    = ltrs
                endif
              enddo izloop
            enddo
          enddo
          ! loop for trs round in the case (SOC & nspin==2)
          if(l_soc.and.nspin==2.and..not.trsoff.and..not.ltrs) then
            ltrs = .true.
            cycle
          endif
          ! allocate and define sym1(:) in the next round
          if(allocated(sym1)) exit
          if(mode==1.and.(invsym==0.or.l_soc.and.nspin==1).and..not.trsoff) then ; allocate(sym1(nsym1*2)) ! allow space for time-reversal symmetries (see below)
          else                                                                   ; allocate(sym1(nsym1))
          endif
          nsym1 = 0
          ltrs  = .false.
        enddo

        if(mode==0.or.iturn==2) then ! checking mode
          if(.not.allocated(sym)) Bug('sym(:) not allocated.')
          if     (nsym1<nsym) then ; str = 'less'
          else if(nsym1>nsym) then ; str = 'more'
          endif
          if(.not.(l_soc.and.nspin==2)) str(len_trim(str)+2:) = 'spatial'
          if(nsym1/=nsym) then
            if(iturn==1) then ; Warn('Found '//trim(str)//' symmetries: '//Chr(nsym1)//' vs '//Chr(nsym)//' in normal run.')
            else              ; Warn('Found '//trim(str)//' symmetries: '//Chr(nsym1)//' vs '//Chr(nsym)//' in sloppy run.')
            endif
          endif
          allocate(iarr(nsym))
          iarr = 0
          do i = 1,nsym
            do j = 1,nsym1
              dist = sym(i)%transl - sym1(j)%transl
              if( all(sym1(j)%rot==sym(i)%rot) .and. all(abs(dist-nint(dist))<1d-12) .and. (sym1(j)%trs.eqv.sym(i)%trs) ) iarr(i)=j
            enddo
          enddo
          if(any(iarr==0)) then
            write(6,'(/A'NoA) 'Input symmetry operations not found'
            if(iturn==1) then ; write(6,'(A)') ' in normal run:'
            else              ; write(6,'(A)') ' in sloopy run:'
            endif
            do i = 1,nsym
              ctrs = ' '
              if(sym(i)%trs) ctrs(1) = '  TRS'
              if(iarr(i)==0) write(6,'(3(/2X,3I3,F10.5,A))') (sym(i)%rot(j,:),sym(i)%transl(j),ctrs(j),j=1,3)
            enddo
          endif
          found = .false.
          do i = 1,nsym1
            found = found .or. all(iarr/=i)
          enddo
          if(found) then
            write(6,'(/A'NoA) 'Extra symmetry operations found'
            if(iturn==1) then ; write(6,'(A)') ' in normal run:'
            else              ; write(6,'(A)') ' in sloopy run:'
            endif
            do i = 1,nsym1
              ctrs = ' '
              if(sym1(i)%trs) ctrs(1) = '  TRS'              
              if(all(iarr/=i)) write(6,'(3(/2X,3I3,F10.5,A))') (sym1(i)%rot(j,:),sym1(i)%transl(j),ctrs(j),j=1,3)
            enddo
          endif
          deallocate(iarr)
        else if(iturn==1) then ! defining mode
          if(allocated(sym)) Bug('sym(:) already allocated.')
          ! include time-reversal here  [ cases (inversion symmetry & noSOC) and (SOC & nspin==1) ]
          if((invsym==0.or.l_soc.and.nspin==1).and..not.trsoff) then
            do i = 1,nsym1
              sym1(nsym+i)     = sym1(i)
              sym1(nsym+i)%trs = .true.
            enddo
            nsym1 = 2*nsym1
          endif
          ! put identity on the first index
          do i = 1,nsym1
            if(all(sym1(i)%rot==identity(3)).and..not.sym1(i)%trs) then
              mrot           = sym1(1)%rot
              trans          = sym1(1)%transl
              ltrs           = sym1(1)%trs
              sym1(1)        = sym1(i)
              sym1(i)%rot    = mrot
              sym1(i)%transl = trans
              sym1(i)%trs    = ltrs
              exit
            endif
          enddo
          if(any(sym1(1)%rot/=identity(3))) Bug('First symmetry operation is not the identity.')
          nsym = nsym1 ; allocate(sym(nsym))
          sym  = sym1
        endif
        deallocate(sym1)
        
      enddo
      
      end
      
