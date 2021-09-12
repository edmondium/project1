c--------------------------------------------------------------------------------
c Copyright (c) 2021 Peter Grünberg Institut, Forschungszentrum Jülich, Germany
c This file is part of SPEX and available as free software under the conditions
c of the MIT license as expressed in the LICENSE file in more detail.
c--------------------------------------------------------------------------------

# include "cppmacro.h"

      module arrays

      use, intrinsic :: iso_fortran_env

c Frequency meshes
      integer                  :: nfreq,nfreqr,nfreqc
      real_dp,     allocatable :: freq(:),freqr(:)
      complex_dp,  allocatable :: freqc(:)

c Self-energy
      MCOMPLEX_dp, allocatable :: selfx(:)
      complex_dp,  pointer_cnt :: selfc(:,:)
      Mpi( integer             :: win_selfc )

c Wannier Coulomb matrices
      complex_dp, allocatable  :: screenw(:,:,:,:,:,:),barew(:,:,:,:,:)
      integer,    allocatable  :: rsite(:,:)
      integer                  :: nsite

      end module arrays

