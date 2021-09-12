/**** Preprocessor macros for Spex. ****/


/* Undefine macros from a previous include (avoids "redefine" warnings). */
  
# undef real_dp
# undef complex_dp
# undef integer_dp
# undef int_dp
# undef MCOMPLEX_dp
# undef MINV
# undef MCONJG
# undef MUNIT
# undef MBYTES
# undef MMCOMPLEX
# undef Inv
# undef NoInv
# undef InvC
# undef NoInvC
# undef ifInv
# undef Nacc
# undef Nacc1
# undef pointer_cnt
# undef target_cnt


/* Define macros. */
  
# ifdef unbug82065 /* avoids bug82065 of GFortran */
#   define real_dp     real(8)
#   define complex_dp  complex(8)
#   define integer_dp  integer(8)
#   define int_dp      8
# else
#   define real_dp     real(real64)
#   define complex_dp  complex(real64)
#   define integer_dp  integer(int64)
#   define int_dp      int64
# endif
# define pointer_cnt  pointer, contiguous
# define target_cnt   target,  contiguous

# ifdef noDOLLAR
#   define NoA //')',advance="no"
# else
#   define NoA //',$)'
# endif

# ifdef old_trafo
#   warning old_trafo
# endif

# ifdef INV
#   define MINV        .true.
#   define MCOMPLEX_dp real_dp
#   define MCONJG
#   define MUNIT       1d0
#   define MBYTES      8d0
#   define Inv(arg)    arg
#   define NoInv(arg)
#   define InvC(arg)   ,arg
#   define NoInvC(arg)
#   define ifInv(arg1,arg2) arg1
# else
#   define MINV        .false.
#   define MCOMPLEX_dp complex_dp
#   define MCONJG      conjg
#   define MUNIT       (1d0,0d0)
#   define MBYTES      16d0
#   define Inv(arg)
#   define NoInv(arg)  arg
#   define InvC(arg)
#   define NoInvC(arg) ,arg
#   define ifInv(arg1,arg2) arg2
# endif 

# ifdef CHECK
#   define Check(arg) arg
# else
#   define Check(arg)
# endif

# ifdef TIMING
#   define Timing(arg) arg
# else
#   define Timing(arg)
# endif

# ifdef MPI
#   ifdef INV
#     define MMCOMPLEX mpi_double_precision
#   else
#     define MMCOMPLEX mpi_double_complex
#   endif
# endif  

# if ( defined(INV) && defined(SOC) )
#   error "Spin-orbit coupling; compile without -DINV!"
# endif

# ifdef CPPTRAD
#   define STRINGIFY(arg) "arg" /* for cpp -traditional */
# else  
#   define STRINGIFY(arg) #arg  /* for fpp/cpp */
# endif

# ifdef M_SIZEOF
#   define sizeof(arg) (size(arg)*(storage_size(arg)/8))
# endif

# ifdef noMSYNC
#   define mpi_f_sync_reg(arr) mpi_get_address(arr,Maddress_dummy,Merr)
# elif defined(noMSYNC2)
#   define mpi_f_sync_reg(arr) donothing(arr)  
# endif

# define Allocate_(array,dim)      allocate(array dim);call checkmem(STRINGIFY(array),1d0*sizeof(array))
# define tAllocate(array,dim) then;allocate(array dim);call checkmem(STRINGIFY(array),1d0*sizeof(array));endif
# define Deallocate_(array)      call checkmem(STRINGIFY(array),-1d0*sizeof(array));deallocate(array)
# define tDeallocate(array) then;call checkmem(STRINGIFY(array),-1d0*sizeof(array));deallocate(array);endif

# define CHKMEM(array)  call checkmem(STRINGIFY(array), 1d0*sizeof(array))
# define CHKMEM_(array) call checkmem(STRINGIFY(array),-1d0*sizeof(array))
# ifdef WRTMEM
#   define Wrtmem(arg) arg
#   define CHKMEM0(info) call checkmem(STRINGIFY(info),0d0)
# else
#   define Wrtmem(arg)
#   define CHKMEM0(info)
# endif

# define Warn(arg)  write(0,'(A,I4.4,A)') 'SPEX-WARNING ('//__FILE__//':',__LINE__,') '//arg
# define Info(arg)  write(0,'(A,I4.4,A)') 'SPEX-INFO ('//__FILE__//':',__LINE__,') '//arg
# define Error(arg) call spex_error('ERROR',__FILE__,__LINE__,arg)
# define Bug(arg)   call spex_error('BUG',__FILE__,__LINE__,arg)
# define Chr(arg)           trim(chr(arg))
# define Chn(arg,sep)       trim(adjustl(chr(arg,sep)))
# define Chf(arg,form)      trim(adjustl(chr(arg,form)))
# define Chfn(arg,form,sep) trim(adjustl(chr(arg,form,sep)))


/* This is an MPI limit for message sizes, which affects parallel HDF5. Spex divides datasets if they exceed the limit. */
# define MaxChunk huge(0)
# define MaxChunkBytes


/* MPI macros */

# ifdef MPI
#   define Mpi(arg) arg
#   define Mpi2(arg1,arg2) arg1,arg2
#   define NoMpi(arg)
#   define MpiC(arg) ,arg
#   define MpiR(arg) if(Mrank.eq.0) then; arg ;endif
#   define MpiO(arg) if(Nrank.eq.0) then; arg ;endif
#   define MnoR(arg) if(Mrank.ne.0) then; arg ;endif 
#   define MnoO(arg) if(Nrank.ne.0) then; arg ;endif 
#   define ifR      if(Mrank.eq.0)
#   define ifO      if(Nrank.eq.0)
#   define Rbegin   if(Mrank.eq.0) then
#   define Obegin   if(Nrank.eq.0) then
#   define Rend     endif
#   define Oend     endif
#   define andR     .and.Mrank.eq.0
#   define Rwrite   if(Mrank.eq.0) write
#   define Rread    if(Mrank.eq.0) read
#   define Rcall    if(Mrank.eq.0) call
#   define Ocall    if(Nrank.eq.0) call
#   define Rif(arg) if(Mrank.eq.0.and.(arg))
#   define Oif(arg) if(Nrank.eq.0.and.(arg))  
#   define beginSingle       if(Mrank.ne.0) then ; call begin_split(mpi_undefined) ; else ; call begin_split(0)
#   define endSingle         endif ; call end_split
#   define beginSplit(color) call begin_split(color)
#   define changeSplit(diff) call change_split(diff)
#   define endSplit          call end_split
#   define beginSplitNodes(color) call begin_split_nodes(color)
#   define endSplitNodes          call end_split_nodes
#   define ifMpi(arg1,arg2) arg1
#   define lMOD(arg) .and.mod(arg,Msize)==Mrank
#   define ifMOD(arg) if(mod(arg,Msize).eq.Mrank) then
#   define ifNMOD(arg) if(mod(arg,Nsize).eq.Nrank) then
#   define ifMODP(arg) arg=arg+1;if(mod(arg,Msize).eq.Mrank) then
#   define elseMOD(arg) else ; arg ; endif
#   define endMOD endif
#   define Mcycle(arg) if(mod(arg,Msize).ne.Mrank) cycle
#   define McycleP(arg)   arg=arg+1;if(mod(arg,Msize).ne.Mrank) cycle
#   define Mstride(m1,m2) m1+Mrank,m2,Msize
#   define MrangeRank1(i,n) (((i)*Msize-1)/(n))
#   define MrangeDef1(m1,m2,n) m1=1+Mrank*(n)/Msize;m2=(1+Mrank)*(n)/Msize
#   define MrangeDef(m1,m2,n1,n2) m1=n1+Mrank*(n2-n1+1)/Msize;m2=n1-1+(1+Mrank)*(n2-n1+1)/Msize
#   define NrangeDef(m1,m2,n1,n2) m1=n1+Nrank*(n2-n1+1)/Nsize;m2=n1-1+(1+Nrank)*(n2-n1+1)/Nsize
#   define Mrange1(n)   (Mrank*(n)/Msize)+1,((1+Mrank)*(n)/Msize)
#   define Mcol1(n)     (Mrank*(n)/Msize)+1:((1+Mrank)*(n)/Msize)
#   define Mnum1(n)    (-Mrank*(n)/Msize   +((1+Mrank)*(n)/Msize))
#   define Mrange(n1,n2) (Mrank*(n2-n1+1)/Msize)+n1,n1-1+((1+Mrank)*(n2-n1+1)/Msize)
#   define Nrange(n1,n2) (Nrank*(n2-n1+1)/Nsize)+n1,n1-1+((1+Nrank)*(n2-n1+1)/Nsize)
#   define Mcol(n1,n2)   (Mrank*(n2-n1+1)/Msize)+n1:n1-1+((1+Mrank)*(n2-n1+1)/Msize)
#   define Ncol(n1,n2)   (Nrank*(n2-n1+1)/Nsize)+n1:n1-1+((1+Nrank)*(n2-n1+1)/Nsize)
#   define Mdistr(array,indx,n1,n2) do indx=n1,n2;call Mcast(array,mod(indx,Msize));enddo
#   define MrangeDistr(array,irnk)          do irnk=0,Msize-1;call Mcast(array,irnk);enddo
#   define OrangeDistr(array,irnk) Obegin ; do irnk=0,Msize-1;call Mcast(array,Opnt(irnk),comm=Ocomm);enddo ; Oend
#   define McolD1(n,irnk)     (irnk*(n)/Msize)+1:((1+irnk)*(n)/Msize)
#   define McolD(n1,n2,irnk)  (irnk*(n2-n1+1)/Msize)+n1:n1-1+((1+irnk)*(n2-n1+1)/Msize)
#   define Ncycle(arg) if(mod(arg,Nsize).ne.Nrank) cycle
#   define NcycleP(arg) arg=arg+1;if(mod(arg,Nsize).ne.Nrank) cycle
#   define Ocycle(arg) if(mod(arg,Osize).ne.Orank) cycle
#   define OcycleP(arg) arg=arg+1;if(mod(arg,Osize).ne.Orank) cycle
#   define Pcycle(indx,n,l)     if(pcycle_skip(indx,n,0d0,l)) cycle
#   define PcycleM(indx,n,xm,l) if(pcycle_skip(indx,n,xm,l))  cycle
#   define Pcycle_not(indx,n,l) if(.not.pcycle_skip(indx,n,0d0,l)) cycle
#   define PDistr(indx,indx1,arr,step) do irank_ = 0,Msize-1,step ; indx1 = indx - Mrank/step + irank_/step ; color_ = 0 ; if(Mrank-irank_>0.and.Mrank-irank_<step) color_ = mpi_undefined ; beginSplit(color_) ; if(color_==0) call Mcast(arr,rank=irank_) ; endSplit ; enddo
#   define RError(arg) ifR Error(arg)
#   define RWarn(arg)  ifR Warn(arg)
#   define RInfo(arg)  ifR Info(arg)
#   define RBug(arg)   ifR Bug(arg)
#   define OError(arg) ifO Error(arg)
#   define OWarn(arg)  ifO Warn(arg)
#   define OInfo(arg)  ifO Info(arg)
#   define OBug(arg)   ifO Bug(arg)
#   define Win(array) array,win_ array
# else
#   define Mpi(arg)
#   define Mpi2(arg1,arg2)
#   define NoMpi(arg) arg
#   define MpiC(arg)
#   define MpiR(arg)
#   define MpiO(arg) continue
#   define MnoR(arg)
#   define MnoO(arg)
/* #   define NoOMpiR(arg) */
#   define ifR
#   define ifO
#   define Rbegin continue
#   define Obegin continue
#   define Rend
#   define Oend
#   define andR
#   define Rwrite   write
#   define Rread    read
#   define Rcall    call
#   define Ocall    call
#   define Rif(arg) if(arg)
#   define Oif(arg) if(arg)
#   define beginSingle continue
#   define endSingle
#   define beginSplit(color)
#   define changeSplit(diff)					  
#   define endSplit continue
#   define beginSplitNodes(color)
#   define endSplitNodes continue
#   define ifMpi(arg1,arg2) arg2
#   define lMOD(arg)
#   define ifMOD(arg)
#   define ifNMOD(arg)					  
#   define ifMODP(arg)
#   define elseMOD(arg)
#   define endMOD
#   define Mcycle(arg)
#   define McycleP(arg)
#   define Ncycle(arg)
#   define NcycleP(arg)
#   define Ocycle(arg)
#   define OcycleP(arg)
#   define Pcycle(indx,n,l)
#   define PcycleM(indx,n,xm,l)
#   define Pcycle_not(indx,n,l)
#   define PDistr(indx,indx1,arr,step)
#   define Mstride(m1,m2) m1,m2
#   define MrangeDef1(m1,m2,n) m1=1;m2=n
#   define MrangeDef(m1,m2,n1,n2) m1=n1;m2=n2
#   define NrangeDef(m1,m2,n1,n2) m1=n1;m2=n2
#   define Mrange1(n) 0+1,n
#   define Mcol1(n)   0+1:n
#   define Mnum1(n)   n
#   define Mrange(n1,n2) 0+n1,n2
#   define Nrange(n1,n2) 0+n1,n2
#   define Mcol(n1,n2)   0+n1:n2
#   define Ncol(n1,n2)   0+n1:n2
#   define Mdistr(array,indx,n1,n2)
#   define MrangeDistr(array,irnk)
#   define OrangeDistr(array,irnk)  
#   define Mdiagonalize diagonalize
#   define Minverse inverse
#   define RError(arg) Error(arg)
#   define RWarn(arg)  Warn(arg)
#   define RInfo(arg)  Info(arg)
#   define RBug(arg)   Bug(arg)
#   define OError(arg) Error(arg)
#   define OWarn(arg)  Warn(arg)
#   define OInfo(arg)  Info(arg)
#   define OBug(arg)   Bug(arg)
#   define Win(array) array
# endif


/* If defined, calls to MATMAT (wrapping *GEMM) are replaced by the internal MATMUL if the call is in capitals. */
/* # define MATMAT matmul */


/* MPI shared-memory macros */

# ifdef MPI
#   define Ncommon(var) call node_allocate(win_ var,ptr,int(storage_size(var)/8),[1]);call c_f_pointer(ptr,var);Nfence_(var)
#   define Nfree(var)   call mpi_win_free(win_ var,Merr)
#   define Nfence_(arr)  call mpi_f_sync_reg(arr) ; call mpi_win_fence(0,win_ arr,Merr) ; call mpi_f_sync_reg(arr)
#   define Nacc1_r_(arr,pos,buf) call mpi_accumulate(buf,1,MPI_DOUBLE_PRECISION,0, transfer(c_loc(arr pos),1_MPI_ADDRESS_KIND)-transfer(c_loc(arr),1_MPI_ADDRESS_KIND),        1,MPI_DOUBLE_PRECISION,MPI_SUM,win_ arr,Merr)
# else
#   define Ncommon(var)
#   define Nfree(var)
#   define Nfence_(arr) continue
#   define Nacc1_r_(arr,pos,buf) arr pos = arr pos + (buf)
# endif

/*  def_Ninit: Initialization of shared arrays is done by all processes to optimize (balance) memory access. */
/*  Note: Using def_Ninit for initializing an array (e.g., array=0) was observed to lead to memory doubling. */
/* # define def_Ninit */

# if defined(MPI) && !defined(noSHARED)
#   define Nfence(arr)  Nfence_(arr)
#   define Nallocate0(arr,dim) call node_allocate(win_ arr,ptr,int(storage_size(arr)/8),dim);call c_f_pointer(ptr,arr,dim);Ocall checkmem('#'//adjustl(STRINGIFY(arr)),1d0*sizeof(arr));Nfence(arr);Ninit(arr,0)
#   define Nallocate(arr,dim)  call node_allocate(win_ arr,ptr,int(storage_size(arr)/8),dim);call c_f_pointer(ptr,arr,dim);Ocall checkmem('#'//adjustl(STRINGIFY(arr)),1d0*sizeof(arr));Nfence(arr)
#   define Ndeallocate(arr) Ocall checkmem('#'//adjustl(STRINGIFY(arr)),-1d0*sizeof(arr));call mpi_win_free(win_ arr,Merr);nullify(arr)
#   define tNdeallocate(arr) then;Ndeallocate(arr);endif
#   define Nallocated(arr) associated(arr)
#   define S_ /
#   define Nacc_r(arr,pos,buf)  call mpi_accumulate(buf,size(buf),MPI_DOUBLE_PRECISION,0, transfer(c_loc(arr pos),1_MPI_ADDRESS_KIND)-transfer(c_loc(arr),1_MPI_ADDRESS_KIND),size(buf),MPI_DOUBLE_PRECISION,MPI_SUM,win_ arr,Merr)
#   define Nacc_c(arr,pos,buf)  call mpi_accumulate(buf,size(buf),MPI_DOUBLE_COMPLEX,  0, transfer(c_loc(arr pos),1_MPI_ADDRESS_KIND)-transfer(c_loc(arr),1_MPI_ADDRESS_KIND),size(buf),MPI_DOUBLE_COMPLEX,  MPI_SUM,win_ arr,Merr)
#   define Nacc1_r(arr,pos,buf) Nacc1_r_(arr,pos,buf)
#   define Nacc1_c(arr,pos,buf) call mpi_accumulate(buf,        1,MPI_DOUBLE_COMPLEX,  0, transfer(c_loc(arr pos),1_MPI_ADDRESS_KIND)-transfer(c_loc(arr),1_MPI_ADDRESS_KIND),        1,MPI_DOUBLE_COMPLEX,  MPI_SUM,win_ arr,Merr)
#   define Macc1_r(arr,pos,rk,buf) call mpi_accumulate(buf,     1,MPI_DOUBLE_PRECISION,rk,transfer(c_loc(arr pos),1_MPI_ADDRESS_KIND)-transfer(c_loc(arr),1_MPI_ADDRESS_KIND),        1,MPI_DOUBLE_PRECISION,MPI_SUM,win_ arr,Merr)
#   ifdef def_Ninit
#     define Ninit(array,val) Ninit_ array(1:size(array,kind=c_size_t)) => array ; Ninit_ array(Ncol(1,size(array,kind=c_size_t))) = val ; nullify(Ninit_ array) ; Nfence(array)
#   else
#     define Ninit(array,val) ifO array = val ; Nfence(array)
#   endif
# else
#   define Nfence(arr) continue
#   define Nallocate0(arr,dim) allocate(arr dim);call checkmem(STRINGIFY(arr),1d0*sizeof(arr)) ; arr=0
#   define Nallocate(arr,dim)  allocate(arr dim);call checkmem(STRINGIFY(arr),1d0*sizeof(arr)) 
#   define Ndeallocate(arr) call checkmem(STRINGIFY(arr),-1d0*sizeof(arr));deallocate(arr)
#   define tNdeallocate(arr) then;Ndeallocate(arr);endif
#   define Nallocated(arr) allocated(arr)
#   define S_
#   define Nacc_r(arr,pos,buf)  arr pos = arr pos + (buf)
#   define Nacc_c(arr,pos,buf)  arr pos = arr pos + (buf)
#   define Nacc1_r(arr,pos,buf) arr pos = arr pos + (buf)
#   define Nacc1_c(arr,pos,buf) arr pos = arr pos + (buf)
#   define Macc1_r(arr,pos,rk,buf) arr pos = arr pos + (buf)
#   define Ninit(array,val) array = val
# endif

# ifdef INV
#   define Nacc(arr,pos,buf)  Nacc_r(arr,pos,buf)
#   define Nacc1(arr,pos,buf) Nacc1_r(arr,pos,buf)
# else
#   define Nacc(arr,pos,buf)  Nacc_c(arr,pos,buf)
#   define Nacc1(arr,pos,buf) Nacc1_c(arr,pos,buf)
# endif

# ifdef LOAD
#   define Load(arg) arg
#   define LoadC(arg) ,arg
#   define NoLoad(arg)
#   define ifLoad(arg1,arg2) arg1
# else
#   define Load(arg)
#   define LoadC(arg)
#   define NoLoad(arg) arg
#   define ifLoad(arg1,arg2) arg2
# endif

# ifdef myLAPACK
#   define zheevx my_zheevx
#   define ZHEEVX my_zheevx
#   define ZLANHE my_zlanhe
# endif

# if defined(MPI) && defined(mySCALAPACK)
#   define pzhegvx   my_pzhegvx
#   define PZHEGVX   my_pzhegvx
#   define pzheevx   my_pzheevx
#   define PZHEEVX   my_pzheevx
#   define PZUNMTR   my_pzunmtr
#   define PCHK1MAT  my_pchk1mat
#   define PCHK2MAT  my_pchk2mat
#   define GLOBCHK   my_globchk
#   define PZSTEIN   my_pzstein
#   define PZLAEVSWP my_pzlaevswp
#   define PZUNMQL   my_pzunmql
#   define PZLARFB   my_pzlarfb
# endif

/* Macro for                  */
/* Intel:    __INTEL_COMPILER */
/* Gfortran: __GFORTRAN__     */
/* NAG:      NAGFOR           */

