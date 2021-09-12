c--------------------------------------------------------------------------------
c Copyright (c) 2021 Peter Grünberg Institut, Forschungszentrum Jülich, Germany
c This file is part of SPEX and available as free software under the conditions
c of the MIT license as expressed in the LICENSE file in more detail.
c--------------------------------------------------------------------------------

c Template "readwrite" module.
c      
c ----------------------------------
c
c This module contains file input and output routines for the transfer of data between Spex and
c the FLAPW-DFT program.
c
c The variables and arrays to be read or written are described at the beginning of each routine.
c Arrays that are not allocated in getinput must be allocated here. See below and "global.f" for
c the datatypes of the variables and arrays.
c
c Name this file "readwrite_PROG.f" where PROG is the name of the DFT program. Configuration with
c the option "--with-dft=PROG" will prepare this file for compilation instead of the default
c "readwrite.f".
c
c
c Notes:
c
c - Formatted and unformatted files should be opened with the function "fopen" (which returns a file unit
c   number) and closed with the routine "fclose" (see "file.f").
c - For HDF5 files, it is recommended to use the wrapper routines of "Hwrapper.f".
c - Parameters that are needed for the data transfer with the DFT code can be defined globally for
c   all routines below (*).
c - MPI: All routines are called by the root process. (In order for parallel HDF5 to work properly,
c   the current communicator Mcomm contains root as a single process.) Data distribution (broadcasting)
c   is done in "getinput". The routine "read_wavef" (there are two variants) is an exception:
c   It is called collectively by a number of processes. See respective notes!
c - The call to "kpt_reorder1" at the beginning of a routine reorders the array kpt(:,:) [and, depending
c   on "mode", also ngpt(:) and nband(:,:)] in such a way that all irreducible kpoints are at the
c   beginning of the array (arrays) (indices 1:nkpt1; nkpt1=nkpti+nkpti2), which (hopefully) helps to 
c   simplify the writing of the routines. The call to "kpt_reorder2" at the end of a routine undoes the 
c   reordering and restores the "normal" array structure.
c - Routines that have a call to "Bug(...)" are essential. There are several routines that are optional.
c   These have a call to "Error(...)" and are needed only in particular cases (e.g., for certain keywords
c   or special macros).      
c

# include "cppmacro.h"

      module readwrite

      use global
      use file
      use hdf5
      use Hwrapper
      use, intrinsic :: iso_fortran_env

c Name of DFT program (PROG)
      character(*), parameter :: dftlong    = 'PROG full name and version number'
      character(*), parameter :: dftshort   = 'PROG' ! short name to be queried (e.g., in shell scripts)
      character(*), parameter :: file_wavef = 'FILE with wavefunction data' ! name of file holding wavefunction data
c Parameters that must also be present in routine getinput
      logical                 :: invs
      integer, allocatable    :: gpt1(:,:,:)
      integer, allocatable    :: nlo(:),llo(:,:)
      integer                 :: nkpt1,nband0,kpt_order
      real_dp                 :: fac_kptadd,maxene
c Other global parameters that are needed for input and output can go here (*)

      contains

c ----------------------------------

c
c BASIC LAPW PARAMETERS
c
c Define:      
c   nspin         = number of spins (1 or 2)
c   ncent         = number of atoms (centers) in the unit cell
c   ntype         = number of atom types (groups of symmetry-equivalent atoms)
c
c Allocate:
c   neq(ntype), lcut(ntype), ztype(ntype), cent(3,ncent), nlo(ntype)
c
c Define:      
c   latpar        = lattice parameter in Bohr (only needed for interpretation of cartesian Bloch vectors "[...]"; preset to 1.0)
c   lat(3,3)      = lattice basis vectors (as columns: lat(:,1),lat(:,2),lat(:,3)) in Bohr
c   invs          = inversion symmetry (.true. if plane-wave coefficients are real-valued)
c   l_soc         = spin-orbit coupling (.true. or .false.)
c   sqa(2)        = spin quantization axis (two Euler angles in radian)
c   gcut          = G cutoff (planewave cutoff)
c   neq(ntype)    = number of equivalent atoms (per atom type)
c   lcut(ntype)   = l cutoff (per atom type)
c   ztype(ntype)  = atomic number (per atom type)
c   cent(3,ncent) = atomic position vector (per atom type) in internal coordinates [i.e., relative to lattice basis vectors (lat)]
c   nlo(ntype)    = number of local orbitals (per atom type)
c
c Allocate:
c   llo(maxval(nlo),ntype)
c
c Define:
c   llo(maxval(nlo),ntype) = l quantum numbers for the local orbitals (per atom type)      
c      
      subroutine read_param
      implicit none
      ! ...
      Bug('Routine is a stub: read_param.')
      ! ...
      end subroutine read_param

c ----------------------------------

c
c SYMMETRY INFORMATION
c
c Define: 
c   nsym                = number of spatial symmetry operations
c
c Allocate:
c   sym(nsym)           = definition of symmetry operations (sym(1) must be the identity operation)
c
c Define:
c   sym(nsym)%rot(3,3)  = (integer) rotation matrix in internal coordinates (i.e., with respect to lattice vectors)
c   sym(nsym)%transl(3) = translational vector of non-symmorphic symmetry operation in internal coordinates
c      
      subroutine read_symmetry
      implicit none
      ! ...
      Bug('Routine is a stub: read_symmetry.')
      ! ...
      end subroutine read_symmetry

c ----------------------------------

c
c RADIAL GRID AND BASIS FUNCTIONS (u, udot, ulo...)
c
c Define:
c   grid(ntype)%number      = number of radial grid points (per atom type)
c   grid(ntype)%radius      = MT radii in Bohr
c   maxgrid                 = maximal number of radial grid points [maxval(grid%number)]
c   lexpgrid                = .true. if grids {r_n} are exponential: r_n = r_0*exp(h*n), n=0,1,2,3,...,N (N=%number, r_0=%first, h=%increment),
c                           otherwise grids are interpolated to an "optimal" exponential grid.
c ---
c
c Depending on the definition of the radial grid, do (1), (2), or (3):
c
c (1) Radial grid {r_n} is exponential: r_n = r_0*exp(h*n), n=0,1,2,3,...,N (N=grid(:)%number)
c
c   Define:
c     grid(ntype)%first     = first grid point r_0 (nonzero) in Bohr
c     grid(ntype)%increment = natural logarithm h of multiplicative increment exp(h)
c
c (2) Radial grid {r_n} is not exponential. DFT code can provide integration weights. !!NOTE: This option is currently disabled!     
c
c   Allocate:
c     rgrid(maxgrid,ntype)
c     gridf(maxgrid,ntype)      
c
c   Define:
c     rgrid(maxgrid,ntype)  = radial grid points in Bohr
c     gridf(maxgrid,ntype)  = radial integration weights
c           
c (3) Radial grid {r_n} is not exponential. DFT code does not provide integration weights.
c
c   Allocate:
c     rgrid(maxgrid,ntype)
c
c   Define:
c     rgrid(maxgrid,ntype)  = radial grid points in Bohr
c
c ---
c
c Define:
c   bas1(maxgrid,maxindx,0:maxlcut,ntype,nspin) = large component of radial functions (multiplied with r)
c   bas2(maxgrid,maxindx,0:maxlcut,ntype,nspin) = small component of radial functions (multiplied with r)
c   ubas(maxindx,0:maxlcut,ntype,nspin)         = value of radial function at MT radius (without r)
c   dubas(maxindx,0:maxlcut,ntype,nspin)        = value of gradient of radial function at MT radius (without r) (in Bohr^-1)
c
c Notes:
c   - The option (1) corresponds to the "standard" definition (of Fleur). Option (2) allows external integration weights
c     to be employed (thus avoiding a loss of precision by interpolation). Option (3) will interpolate the present
c     radial grid [rgrid(:,:)] to the "standard" definition (causing a loss of numerical precision).
c   - The routine getinput "knows" what option has been chosen [(1), (2), or (3)] from the allocation state of the arrays
c     rgrid(:,:) and gridf(:,:).
c   - In case of interpolation (3), the array rgrid(:,:) will be replaced by an exponential grid. The original grid is
c     copied to rgrid0(:,:), which might be helpful for the routine read_pot.      
c   - maxlcut is the maximal l-quantum number per atom: maxlcut = maxval(lcut) with lcut(1:ntype) defined in read_param.
c   - maxindx is the maximal number of radial functions per atom. If there are no local orbitals, then maxindx=2.
c     If there is one local orbital (or a shell of local orbitals), then maxindx=3, etc.
c   - The radial-function index n (second argument of bas1 and bas2, first argument of ubas and dubas), is ordered
c     according to
c     n = 1 : solution of radial Schroedinger (or Dirac) equation [usually called u_l(r) or \phi_l(r)],
c     n = 2 : energy derivative (solution of inhomogeneous equation) [usually called udot_l(r) or \phidot_l(r)],
c     n = 3 : first local orbital,
c     n = 4 : second local orbital, etc.
c     n <= nindx(1:ntype) [maxindx=maxval(nindx)].    
c   - The radial functions with indices n=1,3,4... are assumed to be normalized. The energy derivative (n=2) is not
c     normalized. It is uniquely defined by the inhomogeneous radial equation and the condition of orthogonality to
c     to the first radial function (n=1).      
c      
      subroutine read_radial
      implicit none
      ! ...
      Bug('Routine is a stub: read_radial.')
      ! ...
      end subroutine read_radial

c ----------------------------------

c
c CORE STATES
c
c Define:
c   lcutc(ntype)                                   = l cutoff for core electrons (per atom type)
c   maxlcutc                                       = maximal l cutoff [maxval(lcutc)]
c   maxindxc                                       = maximal number of core states (per l quantum number and atom)
c                                                    (e.g., 2p, 3p... for lcore_soc=.false.,
c                                                    2p3/2, 2p1/2, 3p3/2,... for lcore_soc=.true., see below)
c Allocate:
c   nindxc(0:maxlcutc,ntype)
c   ecore(maxindxc,0:maxlcutc,ntype,nspin)
c   core1(maxgrid,maxindxc,0:maxlcutc,ntype,nspin)
c   core2(maxgrid,maxindxc,0:maxlcutc,ntype,nspin)
c
c Define:
c   nindxc(0:maxlcutc,ntype)                       = number of core electrons (per l quantum number and atom type)
c   ecore(maxindxc,0:maxlcutc,ntype,nspin)         = core-state energy (per core state) in Hartree
c   core1(maxgrid,maxindxc,0:maxlcutc,ntype,nspin) = large component of core wave functions
c   core2(maxgrid,maxindxc,0:maxlcutc,ntype,nspin) = small component of core wave functions
c   lcore_soc                                      = .true. if core states are SOC splitted (e.g., p1/2, p3/2)
c 
c Notes:
c   - If lcore_soc=.false., the core-state index (first dimension of ecore) corresponds to the principle number,
c                           e.g., 2p = 1, 3p = 2, ...
c   - If lcore_soc=.true.,  the core-state index labels the SOC splitted states:
c                           1s1/2 = 1, 2s1/2 = 2, 3s1/2 = 3, ... for l=0
c                           2p3/2 = 1, 2p1/2 = 2, 3p3/2 = 3, ... for l=1
c                           3d5/2 = 1, 3d3/2 = 2, 4d5/2 = 3, ... for l=2 etc.
c   - In the case of spin polarization (nspin=2), two sets of core states are required, one for spin up, the other for spin down,
c     even if lcore_soc is defined as true. ("Two sets of core states" means that one set is calculated with the spin-up potential,
c     as if it were the potential of a non-spin-polarized system, and the other set with the spin-down potential.) In the latter case
c     and if keyword CORESOC is given, the core states are refined in the routine getinput as general four-component spinors including
c     SOC and spin polarization. 
c   - The variable "lcore_soc" is later redefined by the keyword CORESOC. If lcore_soc is defined as false here, CORESOC cannot be used.
c      
      subroutine read_core
      implicit none
      ! ...
      Bug('Routine is a stub: read_core.')
      ! ...      
      end subroutine read_core

c ----------------------------------

c
c RECIPROCAL LATTICE VECTORS (G POINTS)
c
c Define:
c   ngpt(nkpt1)          = number of G points (per k point)      
c   maxgpt               = maximal number of G points [maxval(ngpt)]
c
c Allocate:
c   gpt1(3,maxgpt,nkpt1)
c
c Define:
c   gpt1(3,maxgpt,nkpt1) = G vectors
c      
c Notes:
c   - nkpt1  is the number of irreducible k points = nkpti + nkpti2
c   - nkpti                -- " --                 in the unshifted k-point set
c   - nkpti2               -- " --                 in the shifted k-point set (only for KPT +=...)
c   - kpt(3,nkpt1)  irreducible k points
c   - The array kpt is allocated as kpt(3,nkpt2) with nkpt2>nkpt1, because it also includes all
c     symmetry-equivalent k points.      
c
      subroutine read_gpt
      implicit none
      ! ...
      call kpt_reorder1(0)
      ! ...
      Bug('Routine is a stub: read_gpt.')
      ! ...
      call kpt_reorder2(0)
      end subroutine read_gpt

c ----------------------------------
      
c
c BAND ENERGIES
c
c Define:
c   nband(nkpt1,nspin1)       = number of bands (per k point and spin)
c   maxband                   = maximal number of bands [maxval(nband)]
c
c Allocate:
c   ene(maxband,nkpt1,nspin1)
c
c Define:
c   ene(maxband,nkpt1,nspin1) = (size-ordered) energies of eigenstates (per k point and spin) in Hartree
c
c Notes:
c   - nspin1=nspin if spin is a good quantum number (i.e., without SOC: l_soc=.false.), otherwise nspin1=1.
c   - Cutting of degenerate subspaces should be avoided.
c
      subroutine read_ene
      implicit none
      ! ...
      call kpt_reorder1(0)
      ! ...
      Bug('Routine is a stub: read_ene.')
      ! ...
      call kpt_reorder2(0)
      end subroutine read_ene

c ----------------------------------

c
c POTENTIALS
c      
c Argument:
c   mode = 1 : Read full potential and exchange-correlation (xc) potential
c   mode = 2 : Read full potential and exchange (x) potential (for PBE0)
c
c Define:
c   lfirst (local)    = .true. if this is the first call and arrays are to be allocated [allocated(llh)]
c   nlh(ntype)        = number of lattice harmonics (per atom type)
c   maxlh             = maximum number of lattice harmonics [maxval(nlh)]
c
c Allocate (*):
c   llh(maxlh,ntype)
c   nmlh(maxlh,ntype)
c
c Define:      
c   llh(maxlh,ntype)  = l quantum number of lattice harmonics
c   nmlh(maxlh,ntype) = number of "members" in lattice harmonics (i.e., different values of m)
c   maxmlh            = maximum number of "members" [maxval(nmlh)]
c
c Allocate (*):
c   mlh(maxmlh,maxlh,ntype)
c   clh(maxmlh,maxlh,ntype)
c   vmt(maxgrid,maxlh,ntype,nspin)
c   vmt_xc(maxgrid,maxlh,ntype,nspin)
c
c Read:      
c   mlh(maxmlh,maxlh,ntype)           = m quantum number of "members"
c   clh(maxmlh,maxlh,ntype)           = lattice-harmonic expansion coefficient for each "member"
c   vmt(maxgrid,maxlh,ntype,nspin)    = total MT potential in Hartree
c   vmt_xc(maxgrid,maxlh,ntype,nspin) = xc MT potential (or x MT potential) in Hartree
c
c Allocate (*) (**):      
c   vpw(:,:,:,nspin)
c   vpw_xc(:,:,:,nspin)
c
c Define:
c   vpw(:,:,:,nspin)                  = total interstitial potential in Hartree
c   vpw_xc(:,:,:,nspin)               = xc interstitial potential (or x interstitial potential) in Hartree
c
c Notes:
c   (*) Only if lfirst=.true.
c   - The routine may be called more than once (only PBE0). Array allocation only in the first call (lfirst=.true.).
c     When the routine is called again later (for example with mode=2), only vmt_xc and vpw_xc need to be redefined.
c   - The arrays vmt and vmt_xc are defined in terms of lattice harmonics, which are linear combinations of
c     spherical harmonics Y_lm(r) according to
c     L_n(r) = SUM(m) C_nlm * Y_lm(r),
c     where n = 1,...,nlh(itype) (itype=atom type),
c           l = llh(n,itype),
c           m = mlh(i,n,itype) with i=1,...,nmlh(n,itype)
c           C_nlm = clh(i,n,itype)
c   - Important exception: The spherical (l=0) term of the full potential [vmt(:,1,:,:)], corresponding to the first
c     (n=1) lattice harmonic, must be defined as if L_1(r)=1. This means that the factor Y_00=1/sqrt(4pi) must be
c     contained in vmt(:,1,:,:).      
c   - The arrays vpw and vpw_xc have the (integer) G vector (m,n,k) in the first three dimensions: G = m*g1 + n*g2 + k*g3,
c     where g1, g2, and g3 are the basis vectors of the reciprocal lattice.
c   (**) The arrays have to be allocated large enough (see previous point).
c   - The interstitial potential is V(r)=SUM(G) vpw(G)*exp(iGr).
c   - The coefficients vpw(G) are assumed to have been convoluted with the step function, i.e., V(r)=0 if r points into
c     a MT sphere. (Here, V(r)=0 is approximate because the Fourier series is finite.)      
c   - If the radial grid is interpolated to an exponential grid [option (3) of read_radial], then the array rgrid0(:,:)
c     contains the original radial grid.
c      
      subroutine read_pot(mode)
      implicit none
      integer, intent(in) :: mode
      logical             :: lfirst
      ! ...
      if(all(mode/=[1,2])) Bug('Mode for read_pot unknown.')
      lfirst = .not.allocated(llh)
      ! ...
      Bug('Routine is a stub: read_pot.')
      ! ...
      end subroutine read_pot

c ----------------------------------
      
c
c VXC MATRIX
c
c This routine is optional. By default, Spex calculates the exchange-correlation matrix elements 
c <phi_n|vxc|phi_m> using the xc potential read by read_pot. Alternatively, they can be read from
c an external file provided by the DFT program.
c      
c Arguments:
c   ikpt (local)  = k-point index
c   ispin         = spin index
c   nb            = Dimension of vxc matrix
c   band(nb)      = List of band indices      
c   l_qsgw        = .true. if <phi_n|SIGMA(QSGW)|phi_m> should be added to vxc; SIGMA(QSGW) is the
c                   hermitianized self-energy of QSGW      
c   vxcmat(nb,nb) = (Output) vxc matrix <phi_n|vxc|phi_m> in Hartree with n,m = band(1)..band(nb)
c
c Note:
c   - The argument ikpt_in is the k-point index according to the "normal" indexing used internally in
c     Spex. It is transformed to ikpt, which is a k-point index between 1 and nkpt1 (nkpt1=nkpti+nkpti2),
c     the same order as the one used in the other routines.      
c
      subroutine read_vxc(vxcmat,band,nb,ikpt_in,ispin,l_qsgw)
      use, intrinsic :: iso_fortran_env
      implicit none
      integer,     intent(in)  :: nb,band(nb),ikpt_in,ispin
      logical,     intent(in)  :: l_qsgw
      MCOMPLEX_dp, intent(out) :: vxcmat(nb,nb)
      integer                  :: ikpt
      ! ...
      call kpt_reorder1(3)
      if     (ikpt_in<=nkpti)  then ; ikpt = ikpt_in
      else if(ikpt_in==nkpt+1) then ; ikpt = nkpti + 1
      else                          ; Error('k-point index out of range.')
      endif
      ! ...
      Error('Routine is not available: read_vxc.')
      ! ...
      call kpt_reorder2(3)
      if(l_qsgw) call read_qsgw(vxcmat,band,nb,ikpt_in,ispin)
      end subroutine read_vxc

c ----------------------------------
      
c
c SIGMA(QSGW) MATRIX      
c
c Arguments:
c   ikpt (local)  = k-point index (see note in "read_vxc")
c   ispin         = spin index
c   nb            = Dimension of vxc matrix
c   band(nb)      = List of band indices      
c   vxcmat(nb,nb) = (In/output) The matrix <phi_n|SIGMA(QSGW)|phi_m> with n,m = band(1)..band(nb)
c                   is added to the array vxcmat. SIGMA(QSGW) is the hermitianized QSGW self-energy in Hartree.
c                   The matrix must be provided by the DFT program.      
c
c Note:
c   - In the case of self-consistent HF and PBE0 calculations, SIGMA(QSGW) stands for SIGMA(HF) or
c     SIGMA(PBE0).
c      
      subroutine read_qsgw(vxcmat,band,nb,ikpt_in,ispin)
      use, intrinsic :: iso_fortran_env
      implicit none
      integer,     intent(in)    :: nb,band(nb),ikpt_in,ispin
      MCOMPLEX_dp, intent(inout) :: vxcmat(nb,nb)
      integer                    :: ikpt
      ! ...
      call kpt_reorder1(3)
      if     (ikpt_in<=nkpti)  then ; ikpt = ikpt_in
      else if(ikpt_in==nkpt+1) then ; ikpt = nkpti + 1
      else                          ; Error('k-point index out of range.')
      endif
      ! ...
      Error('Routine is not available: read_qsgw.')
      ! ...
      call kpt_reorder2(3)
      end subroutine read_qsgw      

c ----------------------------------

c
c WANNIER U MATRICES
c
c This routine is optional. Spex usually constructs the Wannier U matrices itself or calls the Wannier90
c library for the maximal localization procedure. The U matrices can alternatively be read from an external
c file provided by another program. See the keyword UREAD.
c
c Define
c   nwan     = number of Wannier functions
c   nwanband = number of bands from which Wannier functions are constructed (nwanband=wanbandf-wanbandi+1)
c   wanbandi = lowest band index for Wannier construction
c   wanbandf = highest band index for Wannier construction
c
c Allocate:
c   uwan(wanbandi:wanbandf,nwan,nkpt,nspin2)      
c      
c Define:
c   uwan(wanbandi:wanbandf,nwan,nkpt,nspin2) = Wannier U matrices
c
c Notes:
c   - nkpt is the number of k points in the unshifted set (e.g., 64 for a 4x4x4 sampling).
c   - nspin2=nspin if l_soc=.false. (no SOC), otherwise nspin2=2.
c         
      subroutine read_wannier
      implicit none
      ! ...
      Error('Routine is not available: read_wannier.')
      ! ...
      end subroutine read_wannier

c ----------------------------------

c
c K POINT SET      
c
c Spex requires the eigensolutions (wavefunction and energies) on a special k-point set.
c This routine writes the k-point set (only the irreducible part) to a file. The DFT
c program should read the k-points from this file, and solve the eigenvalue problem.      
c
c Write:
c   kpt(3,nkpt1) = irreducible k points
c
c Note:
c   - nkpt1 = nkpti + nkpti2 (nkpti|nkpti2 = number of irreducible k points in unshifted|shifted set)
c   - The array kpt is allocated as a larger array: kpt(3,nkpt2) with nkpt2>nkpt1.
c   - It might be helpful to use the variable "fac_kptadd" (defined globally in this module). This variable is the denominator in the
c     definition of "KPT +". Example: KPT +=1/3*(0,0,1) - in this case: fac_kptadd=3.
c     The variable can be used as a scaling factor for the shifted k points [kpt(:,nkpti+1:nkpt1)].
c     In the example, kpt(:,nkpti+1)=[0,0,0.33333333..]. To avoid rounding errors (in an ASCII file), one can write instead
c     fac_kptadd*kpt(:,nkpti+1)=[0,0,1] and, of course, the variable "fac_kptadd" into the k-point list.
c
      subroutine write_kpt
      implicit none
      ! ...
      call kpt_reorder1(0)
      ! ...
      Bug('Routine is a stub: write_kpt.')
      ! ...
      call kpt_reorder2(0)
      end subroutine write_kpt

c ----------------------------------

c
c Q POINT PATH
c
c Optional write routine. Needed if the q-point path generated by Spex should be used in the DFT code.
c Should be empty if unused.
c
c Write
c   qpt(3,nqpt) = q-point path for band-structure calculations
c   label(nqpt) = q-point labels (examples: "G", "X", "L")
c
      subroutine write_qpt(qpt,nqpt,label)
      implicit none
      integer,      intent(in) :: nqpt
      real_dp,      intent(in) :: qpt(3,nqpt)
      character(*), intent(in) :: label(nqpt)
      ! ...
      end subroutine write_qpt

c ----------------------------------
      
c
c BASIC LAPW PARAMETERS (write)
c
c Only needed for ITERATE and PLUSSOC.
c
c Argument:
c   lsoc = overrides the SOC flag (l_soc)
c
c Write:
c   all parameters and arrays described for read_param above.
c      
      subroutine write_param(lsoc)
      implicit none
      logical, intent(in) :: lsoc
      ! ...
      Error('Routine is not available: write_param.')
      ! ...
      end subroutine write_param

c ----------------------------------
      
c
c WAVE FUNCTIONS (write)
c
c Only needed for ITERATE and PLUSSOC.
c
c Arguments:
c   ikpt                               = k-point index [k point: kpt(:,ikpt)]
c   ispin                              = spin index (1 or 2)
c   ngpt1                              = number of G vectors
c   gpt1(3,ngpt1)                      = G vectors
c   eig(:)                             = eigenvalues in Hartree
c   cmtin(maxlmindx,ncent,nbnd,nspin3) = MT wave-function coefficients
c   cpwin(maxgpt,nbnd,nspin3)          = PW wave-function coefficients (according to gpt1(:,:))
c
c Other:
c   nbnd                               = number of bands (eigenvalues)
c      
c Write:
c   ngpt1      
c   gpt1(3,ngpt1)
c   nbnd
c   eig(nbnd)
c   cmtin(maxlmindx,ncent,nbnd,nspin3)
c   cpwin(maxgpt,nbnd,nspin3)      
c
c Notes:
c   - k-point index 1 (ikpt=1 and ispin=1) starts a new output file. Overwrite old files if necessary!
c   - nspin3 = number of spin components per state; nspin3=2 if l_soc=.true., otherwise nspin3=1.
c       
      subroutine write_wavef(ikpt,ispin,ngpt1,gpt1,eig,cmtin,cpwin,cpwin_c)
      use global
      use Hwrapper
      use hdf5
      implicit none
      integer,     intent(in)           :: ikpt,ispin,ngpt1
      integer,     intent(in)           :: gpt1(3,ngpt1)
      real_dp,     intent(in)           :: eig(:)
      complex_dp,  intent(in)           :: cmtin(maxlmindx,ncent,size(eig),nspin3)
      MCOMPLEX_dp, intent(in), optional :: cpwin(maxgpt,size(eig),nspin3)
      complex_dp,  intent(in), optional :: cpwin_c(maxgpt,size(eig),nspin3)
      integer                           :: nbnd
      if(present(cpwin).eqv.present(cpwin_c)) Bug('cpwin and cpwin_c both present or missing.')
      if(ngpt1>maxgpt)                        Bug('ngpt1 > maxgpt.')
      nbnd = size(eig)
      ! ...
      Error('Routine not available: write_wavef.')
      ! ...
      end subroutine write_wavef

c ----------------------------------
      
# ifndef LOAD

c ----------------------------------
      
c
c WAVE FUNCTIONS (default)
c
c This is the default routine for reading the wave functions (macro LOAD is not set).
c
c Define:
c   cpw(maxgpt,maxband,nkpt1,nspin2)          = PW wavefunction coefficients
c   cmt(maxlmindx,ncent,maxband,nkpt1,nspin2) = MT wavefunction coefficients
c
c Notes:      
c   - nkpt1 = nkpti + nkpti2 (see above)
c   - If the keyword STOREBZ is given, cpw and cmt are allocated as larger arrays, namely
c     cpw(maxgpt,maxband,nkpt2,nspin2) and cmt(maxlmindx,ncent,maxband,nkpt2,nspin2) with nkpt2>nkpt1 (see above),
c     but only the subarrays cpw(:,:,:nkpt1,:) and cmt(:,:,:,:nkpt1,:) need to be defined.
c   - The first dimension of cmt is an index i counting l- and m-quantum numbers and radial functions (index n).
c     The order is: i=1 : l=0, m=0,  n=1 (u)
c                   i=2 : l=0, m=0,  n=2 (udot)
c                   i=3 : l=1, m=-1, n=1 (u)
c                   i=4 : l=1, m=-1, n=2 (udot)
c                   i=5 : l=1, m=-1, n=3 (ulo)
c                   i=6 : l=1, m=0,  n=1 (u) etc.
c     corresponding to the nested loops (for a given atom-type index "itype"):
c     i = 0
c     do l = 0,lcut(itype)
c       do m = -l,l
c         do n = 1,nindx(l,itype)
c           i = i + 1
c           ...
c         enddo
c       enddo
c     enddo
c   - l- and m-quantum numbers refer to spherical harmonics, which are defined with respect to the global
c     coordinate system.   
c   - nspin2=nspin if l_soc=.false. (no SOC), otherwise nspin2=2.
c   - MPI: The routine is called collectively by all processes. The arrays cpw and cmt are allocated as
c     shared memory (MPI-3) in each node. Every CPU of a node has access to the shared memory (it is the same
c     physical memory). If data is read by a CPU in node 1, it is seen by all CPUs of that node (after the
c     Nfence call), but it is not present in node 2. The data has to be broadcasted in this case. You can use the
c     wrapper routine Mcast of "Mwrapper.f", e.g., Ocall Mcast(cpw,comm=Ocomm). "Ocall" means that Mcast is called
c     only by the processes in Ocomm. (The MPI communicators are: Mcomm=all processes, Ncomm=processes of a given
c     node, Ocomm=the root processes of all nodes).
c     If HDF5 is used, the data can be read collectively, and no broadcasting is necessary because the HDF5 library
c     takes care that all processes (and therefore all nodes) get the data.
c   - MPI: "Nfence" is a data synchronization call that separates memory writing and reading extents (to avoid race
c     conditions).      
c   - Hpos (see below) is the array dimension that offset and stride refer to (rightmost is zeroth, second from right
c     is first, etc.), see Hwrapper routines.      
c      
      subroutine read_wavef
      Mpi( use Mwrapper )
      implicit none
      Mpi ( integer :: Merr )
      ! ...
      call kpt_reorder1(1)
      ! Hpos = ?
      ! ...
      Bug('Routine is a stub: read_wavef.')
      ! ...
      Nfence(cmt)
      Nfence(cpw)
      call kpt_reorder2(1)
      end subroutine read_wavef

c ----------------------------------
      
# else

c ----------------------------------

c
c WAVE FUNCTIONS (-DLOAD)
c
c This routine is only needed if the code is compiled with the macro LOAD. With this macro, the wave functions
c are only read from harddisc when needed. So, they occupy less main memory.      
c 
c Arguments:
c   band(:)                            = list of band indices
c   kpt1(:)                            = list of k-point indices
c   ispin                              = spin-index
c   cmtout(maxlmindx,ncent,:,:,nspin3) = MT wavefunction coefficients (optional)
c   cpwout(maxgpt,:,:,nspin3)          = PW wavefunction coefficients (optional)
c   str                                = stride for band indices: str=1 or str=Msize (Msize = number of processes in current communicator Mcomm)
c                                        If str=Msize, then each process Mrank in Mcomm gets a different set of band indices
c                                        Mrank = 0 : 1 , Msize   , 2*Msize   , ...
c                                        Mrank = 1 : 2 , Msize+1 , 2*Msize+1 , ... etc.
c
c Other:
c   lmt   = .true. if MT coefficients should be read (cmtout is present)
c   lpw   = .true. if PW coefficients should be read (cpwout is present)
c   lbcnt = .true. if band indices are continuous (from b1 to b2) (might be helpful for the read-in process)
c
c Define:
c   cmtout(maxlmindx,ncent,:,:,nspin3)
c   cpwout(maxgpt,:,:,nspin3)      
c      
c Notes:
c   - The arrays cmtout and cpwout are NOT in shared memory. Each process gets its "share" of wave functions.
c   - The third (second) and fourth (third) dimensions of cmtout (cpwout) are larger than size(band) and size(kpt1), respectively.
c
      subroutine read_wavef(band,kpt1,ispin,cmtout,cpwout,str)
      use, intrinsic :: iso_fortran_env
      implicit none
      integer,     intent(in)            :: band(:),kpt1(:),ispin,str
      complex_dp,  intent(out), optional :: cmtout(:,:,:,:,:)
      MCOMPLEX_dp, intent(out), optional :: cpwout(:,:,:,:)
      logical                            :: lmt,lpw,lbcnt
      integer                            :: b1,b2
      b1    = 0
      b2    = 0
      lmt   = present(cmtout) ; if(lmt) then ; cmtout = 0 ; b1 =     size(cmtout,3)     ; b2 =     size(cmtout,4)     ;  endif
      lpw   = present(cpwout) ; if(lpw) then ; cpwout = 0 ; b1 = min(size(cpwout,2),b1) ; b2 = min(size(cpwout,3),b2) ; endif
      lbcnt = all ( [( band(i+1)-band(i)==1 , i=1,size(band)-1 )] )
      if(str/=1 Mpi(.and.str/=Msize) ) Bug('Wrong stride (str) argument.')
      if(b1<size(band))                Bug('Band dimension too small.')
      if(b2<size(kpt1))                Bug('kpoint dimension too small.')
      call kpt_reorder1(3)
      if(lbcnt) then
        b1 = band(1)
        b2 = band(size(band))
      else if(str/=1) then
        Error('not implemented: str/=1 and band(:).')
      endif
      ! Hpos = ?
      ! ...
      Error('Routine is not available: read_wavef (LOAD version).')
      ! ...
      call kpt_reorder2(3)
      end subroutine read_wavef
      
c ----------------------------------
      
# endif

c ----------------------------------

# include "readwrite.inc"

c ----------------------------------      

      end module readwrite
