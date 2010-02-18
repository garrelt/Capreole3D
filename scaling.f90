!>
!! \brief This module defines the scaling parameters.
!!
!! Module for Capreole
!!
!! \b Author: Garrelt Mellema
!!
!! \b Date: 2010-02-04 (2003-06-01)
!!
!! Version: cgs units, circumstellar scales
!!
!! This module is also accepted by the F compiler (Dec 9, 2003)

module scaling

  use precision, only: dp
  private

  real(kind=dp),parameter,public :: scleng=1.0e16_dp !< 10^16 cm
  real(kind=dp),parameter,public :: scmass=1.0e27_dp !< 10^27 g
  real(kind=dp),parameter,public :: sctime=1.0e8_dp  !< 10^8 s
  real(kind=dp),parameter,public :: scvelo=scleng/sctime
  real(kind=dp),parameter,public :: scdens=(scmass/scleng)/(scleng*scleng)
  real(kind=dp),parameter,public :: scmome=scdens*scvelo
  real(kind=dp),parameter,public :: scener=scdens*scvelo*scvelo
  real(kind=dp),parameter,public :: sccool=scener/sctime

end module scaling
