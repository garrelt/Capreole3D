module scaling

  ! Module for Capreole
  ! Author: Garrelt Mellema
  ! Date: 2003-06-01

  ! This module defines the scaling parameters.

  ! Version: No scaling

  use precision, only: dp
  private

  real(kind=dp),parameter,public :: scleng=1d0
  real(kind=dp),parameter,public :: scmass=1d0
  real(kind=dp),parameter,public :: sctime=1d0
  real(kind=dp),parameter,public :: scvelo=scleng/sctime
  real(kind=dp),parameter,public :: scdens=(scmass/scleng)/(scleng*scleng)
  real(kind=dp),parameter,public :: scmome=scdens*scvelo
  real(kind=dp),parameter,public :: scener=scdens*scvelo*scvelo
  real(kind=dp),parameter,public :: sccool=scener/sctime

end module scaling
